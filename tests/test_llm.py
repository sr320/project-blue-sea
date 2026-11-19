import json

import pytest

from bluesea.literature import BLINDING_NOTE, build_prompt, redact
from bluesea.llm import BudgetExceeded, CachedClient

from .conftest import make_record


class FakeClient:
    """Stands in for the Anthropic client; counts calls so the cache can be
    shown to actually prevent them."""

    def __init__(self, response: str = '{"score": 0.7}'):
        self.response = response
        self.calls = 0

    class _Message:
        def __init__(self, text, in_tokens=1000, out_tokens=200):
            block = type("Block", (), {"type": "text", "text": text})()
            self.content = [block]
            self.usage = type(
                "Usage", (), {"input_tokens": in_tokens, "output_tokens": out_tokens}
            )()

    class _Messages:
        def __init__(self, outer):
            self.outer = outer

        def create(self, **kwargs):
            self.outer.calls += 1
            return FakeClient._Message(self.outer.response)

    @property
    def messages(self):
        return FakeClient._Messages(self)


def make_cached_client(tmp_path, **kwargs):
    defaults = dict(
        cache_dir=tmp_path / "cache",
        model="test-model",
        budgets={"literature": 10.0},
        dry_run=False,
        client=FakeClient(),
    )
    defaults.update(kwargs)
    return CachedClient(**defaults)


def test_dry_run_makes_no_calls(tmp_path):
    fake = FakeClient()
    client = make_cached_client(tmp_path, dry_run=True, client=fake)
    assert client.complete("literature", "sys", "prompt") == ""
    assert fake.calls == 0
    report = client.report()
    assert report["dry_run"] is True
    assert report["stages"]["literature"]["calls"] == 1   # projected, not spent
    assert report["total_cost_usd"] > 0


def test_cache_prevents_a_second_call(tmp_path):
    fake = FakeClient()
    client = make_cached_client(tmp_path, client=fake)
    first = client.complete("literature", "sys", "prompt")
    second = client.complete("literature", "sys", "prompt")
    assert first == second
    assert fake.calls == 1


def test_cache_survives_a_new_client(tmp_path):
    """A rerun of the workflow must re-spend nothing it has already paid for."""
    make_cached_client(tmp_path).complete("literature", "sys", "prompt")
    fake = FakeClient()
    fresh = make_cached_client(tmp_path, client=fake)
    assert fresh.complete("literature", "sys", "prompt") == '{"score": 0.7}'
    assert fake.calls == 0


def test_changing_the_prompt_is_a_cache_miss(tmp_path):
    fake = FakeClient()
    client = make_cached_client(tmp_path, client=fake)
    client.complete("literature", "sys", "prompt one")
    client.complete("literature", "sys", "prompt two")
    assert fake.calls == 2


def test_budget_ceiling_halts_the_stage(tmp_path):
    client = make_cached_client(tmp_path, budgets={"literature": 0.001})
    with pytest.raises(BudgetExceeded):
        for index in range(100):
            client.complete("literature", "sys", f"prompt {index}")


def test_complete_json_extracts_from_prose(tmp_path):
    fake = FakeClient('Here is my answer:\n```json\n{"score": 0.4}\n```\nDone.')
    client = make_cached_client(tmp_path, client=fake)
    assert client.complete_json("literature", "sys", "p") == {"score": 0.4}


def test_complete_json_returns_empty_on_garbage(tmp_path):
    client = make_cached_client(tmp_path, client=FakeClient("no json at all"))
    assert client.complete_json("literature", "sys", "p") == {}


def test_provenance_is_written_with_the_response(tmp_path):
    client = make_cached_client(tmp_path)
    client.complete("literature", "system prompt", "user prompt")
    files = list((tmp_path / "cache").rglob("*.json"))
    assert len(files) == 1
    record = json.loads(files[0].read_text())
    assert record["model"] == "test-model"
    assert record["prompt"] == "user prompt"
    assert record["system"] == "system prompt"
    assert record["response"]


def test_blinding_removes_vaccine_terms_from_the_prompt():
    """A positive control must not be scored highly just for being annotated
    as a vaccine antigen."""
    record = make_record("LSAL_P0")
    record.description = "ribosomal protein P0, vaccine antigen, immunization target"
    blinded = build_prompt(
        record, blinded=True, blinded_terms=["vaccine", "immunization"]
    )
    assert "[redacted]" in blinded
    assert BLINDING_NOTE in blinded
    body = blinded.replace(BLINDING_NOTE, "")
    assert "vaccine antigen" not in body.lower()


def test_redact_catches_inflections():
    assert "vaccin" not in redact("vaccinated vaccine vaccination", ["vaccin"]).lower()


def test_unblinded_prompt_keeps_the_annotation():
    record = make_record("LSAL_P0")
    record.description = "ribosomal protein P0, vaccine antigen"
    prompt = build_prompt(record, blinded=False, blinded_terms=["vaccine"])
    assert "vaccine antigen" in prompt
    assert BLINDING_NOTE not in prompt

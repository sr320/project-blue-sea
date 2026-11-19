"""The API stage: cached, budgeted, and separable from the compute nodes.

Hyak compute nodes cannot be assumed to have outbound network access, so this
stage runs on its own: it reads tool output from disk, makes its calls, and
writes JSON back. That split is good practice regardless -- it makes the
expensive stage independently restartable -- but here it is likely mandatory.

Three controls keep spending bounded:

* a disk cache keyed on the exact prompt, so a rerun re-spends nothing;
* a per-stage ceiling that halts the stage rather than overrunning it;
* a dry-run mode that reports projected calls and tokens without calling,
  which is the default for any config change.

For every call the model id, the full prompt and the raw response are kept.
That record is what makes an LLM-judged score defensible: a reviewer who
doubts a literature score can read exactly what was asked and what came back.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# Per-million-token rates used only for budget projection; they do not bill
# anything and should be updated when pricing changes.
DEFAULT_RATES_USD = {"input": 5.0, "output": 25.0}


class BudgetExceeded(RuntimeError):
    """Raised when a stage would cross its ceiling. Halts the stage."""


@dataclass
class Usage:
    calls: int = 0
    input_tokens: int = 0
    output_tokens: int = 0
    cached_calls: int = 0
    cost_usd: float = 0.0

    def add(self, input_tokens: int, output_tokens: int, rates: dict[str, float]) -> None:
        self.calls += 1
        self.input_tokens += input_tokens
        self.output_tokens += output_tokens
        self.cost_usd += (
            input_tokens * rates["input"] + output_tokens * rates["output"]
        ) / 1_000_000


@dataclass
class CallRecord:
    """Full provenance for one call, written next to the response."""

    key: str
    stage: str
    model: str
    prompt: str
    system: str
    response: str
    input_tokens: int
    output_tokens: int
    timestamp: float = field(default_factory=time.time)

    def to_dict(self) -> dict[str, Any]:
        return {
            "key": self.key,
            "stage": self.stage,
            "model": self.model,
            "system": self.system,
            "prompt": self.prompt,
            "response": self.response,
            "input_tokens": self.input_tokens,
            "output_tokens": self.output_tokens,
            "timestamp": self.timestamp,
        }


def estimate_tokens(text: str) -> int:
    """Rough token count for projection only: ~4 characters per token."""
    return max(1, len(text) // 4)


class CachedClient:
    """Anthropic client wrapper with an on-disk cache and spend ceilings.

    The cache key is a hash of model, system prompt, user prompt and the
    sampling parameters, so changing any of them is a cache miss and changing
    none of them never costs twice.
    """

    def __init__(
        self,
        cache_dir: str | Path,
        model: str,
        budgets: dict[str, float],
        dry_run: bool = True,
        max_tokens: int = 4096,
        temperature: float = 0.0,
        rates: dict[str, float] | None = None,
        client: Any | None = None,
    ) -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.model = model
        self.budgets = dict(budgets)
        self.dry_run = dry_run
        self.max_tokens = max_tokens
        self.temperature = temperature
        self.rates = rates or dict(DEFAULT_RATES_USD)
        self.usage: dict[str, Usage] = {}
        self._client = client
        self.projection: dict[str, Usage] = {}

    # -- cache -------------------------------------------------------------

    def _key(self, system: str, prompt: str) -> str:
        payload = json.dumps(
            {
                "model": self.model,
                "system": system,
                "prompt": prompt,
                "max_tokens": self.max_tokens,
                "temperature": self.temperature,
            },
            sort_keys=True,
        )
        return hashlib.sha256(payload.encode()).hexdigest()

    def _cache_path(self, key: str) -> Path:
        # Shard by prefix: a flat directory of tens of thousands of files is
        # slow on the cluster filesystem.
        return self.cache_dir / key[:2] / f"{key}.json"

    def _read_cache(self, key: str) -> CallRecord | None:
        path = self._cache_path(key)
        if not path.exists():
            return None
        with open(path) as handle:
            return CallRecord(**json.load(handle))

    def _write_cache(self, record: CallRecord) -> None:
        path = self._cache_path(record.key)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        with open(tmp, "w") as handle:
            json.dump(record.to_dict(), handle, indent=2)
        tmp.replace(path)   # atomic, so a killed job cannot leave a half file

    # -- budget ------------------------------------------------------------

    def _usage(self, stage: str) -> Usage:
        return self.usage.setdefault(stage, Usage())

    def _check_budget(self, stage: str, projected_cost: float) -> None:
        ceiling = self.budgets.get(stage)
        if ceiling is None:
            return
        spent = self._usage(stage).cost_usd
        if spent + projected_cost > ceiling:
            raise BudgetExceeded(
                f"stage '{stage}' would reach ${spent + projected_cost:.2f} "
                f"against a ceiling of ${ceiling:.2f}"
            )

    # -- calling -----------------------------------------------------------

    def complete(self, stage: str, system: str, prompt: str) -> str:
        """Return the model's text, from cache when possible.

        In dry-run mode nothing is called: the projected call is tallied and
        an empty string comes back, so a caller can walk the whole stage and
        report what a real run would cost.
        """
        key = self._key(system, prompt)
        cached = self._read_cache(key)
        if cached is not None:
            self._usage(stage).cached_calls += 1
            return cached.response

        input_tokens = estimate_tokens(system + prompt)
        projected_cost = (
            input_tokens * self.rates["input"]
            + self.max_tokens * self.rates["output"]
        ) / 1_000_000

        if self.dry_run:
            projected = self.projection.setdefault(stage, Usage())
            projected.add(input_tokens, self.max_tokens, self.rates)
            return ""

        self._check_budget(stage, projected_cost)
        response, in_tokens, out_tokens = self._call_api(system, prompt)
        self._usage(stage).add(in_tokens, out_tokens, self.rates)
        self._write_cache(
            CallRecord(
                key=key,
                stage=stage,
                model=self.model,
                prompt=prompt,
                system=system,
                response=response,
                input_tokens=in_tokens,
                output_tokens=out_tokens,
            )
        )
        return response

    def _call_api(self, system: str, prompt: str) -> tuple[str, int, int]:
        client = self._client
        if client is None:
            try:
                import anthropic
            except ImportError as error:   # pragma: no cover - env dependent
                raise RuntimeError(
                    "anthropic package is required for a live run; "
                    "install it or keep llm.dry_run enabled"
                ) from error
            if not os.environ.get("ANTHROPIC_API_KEY"):
                raise RuntimeError("ANTHROPIC_API_KEY is not set")
            client = self._client = anthropic.Anthropic()

        message = client.messages.create(
            model=self.model,
            max_tokens=self.max_tokens,
            temperature=self.temperature,
            system=system,
            messages=[{"role": "user", "content": prompt}],
        )
        text = "".join(
            block.text
            for block in message.content
            if getattr(block, "type", "") == "text"
        )
        return text, message.usage.input_tokens, message.usage.output_tokens

    def complete_json(self, stage: str, system: str, prompt: str) -> dict[str, Any]:
        """Parse a JSON object out of the response.

        A model that wraps its JSON in prose or a fenced block is common
        enough that failing the whole stage over it would be wasteful, so the
        outermost brace-delimited span is extracted. A genuinely unparseable
        response returns {} and the caller flags the record rather than
        inventing a score.
        """
        text = self.complete(stage, system, prompt)
        if not text:
            return {}
        start, end = text.find("{"), text.rfind("}")
        if start == -1 or end <= start:
            return {}
        try:
            return json.loads(text[start : end + 1])
        except json.JSONDecodeError:
            return {}

    # -- reporting ---------------------------------------------------------

    def report(self) -> dict[str, Any]:
        source = self.projection if self.dry_run else self.usage
        return {
            "dry_run": self.dry_run,
            "model": self.model,
            "stages": {
                stage: {
                    "calls": usage.calls,
                    "cached_calls": self._usage(stage).cached_calls,
                    "input_tokens": usage.input_tokens,
                    "output_tokens": usage.output_tokens,
                    "cost_usd": round(usage.cost_usd, 2),
                    "ceiling_usd": self.budgets.get(stage),
                }
                for stage, usage in source.items()
            },
            "total_cost_usd": round(sum(u.cost_usd for u in source.values()), 2),
        }


def client_from_config(config: dict[str, Any], paths: dict[str, Any], **overrides: Any):
    llm_cfg = dict(config["llm"])
    llm_cfg.update(overrides)
    return CachedClient(
        cache_dir=Path(paths["cache"]),
        model=llm_cfg["model"],
        budgets=llm_cfg["budget_usd"],
        dry_run=llm_cfg["dry_run"],
        max_tokens=llm_cfg["max_tokens"],
        temperature=llm_cfg["temperature"],
    )

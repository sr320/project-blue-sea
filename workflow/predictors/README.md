This directory is copied to `/opt/predictors` during the container build.
Stage licensed installations and model files here, plus an `install.sh` that
installs their dependencies into the appropriate local environments. Those
files are ignored by Git. Record tool/source versions and required local
files in `config/config.yaml`; see `docs/tools.md`.

The repository does not distribute licensed predictors or model weights.
The baseline container is insufficient for a predictor run. The build script
returns a nonzero status if its post-build predictor check fails.

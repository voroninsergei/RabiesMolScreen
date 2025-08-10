# Architecture (10/10)

Goals:
- Clear boundary with external binaries (vina/smina/obabel) via `ports/` + `adapters/`.
- Standardized results contract with versioned schema (`config/results.schema.json` + `rabiesmol/io/contract.py`).
- Single CLI with explicit steps `prepare -> dock -> rescore -> report` (`rabiesmol/cli/main.py`).
- Everything else (HTML template, configs, examples) is in well-known places.

## Package layout
- `rabiesmol/ports` – abstract interfaces (`DockingEngine`).
- `rabiesmol/adapters` – concrete adapters for external tools (`VinaEngine`, `SminaEngine`, `obabel` helper).
- `rabiesmol/services` – orchestration of steps (pure python, testable).
- `rabiesmol/io` – IO helpers and contract validators.
- `rabiesmol/report` – Jinja2 templates + renderer.
- `rabiesmol/utils` – logging + subprocess runner with crisp errors.
- `config/` – defaults + JSON schema for outputs.
- `examples/` – tiny demo data.
- `reports/` – generated reports (ignored in VCS).
- `scripts/` – optional convenience scripts.
- `tests/` – contract smoke-tests.

## Result contract
Parquet/CSV must have the following columns (types in parentheses):

- `schema_version (string)` – value `1.0.0` today.
- `ligand_id (string)`
- `protein_id (string)`
- `engine (string)`
- `pose_id (int)`
- `score (float)`
- `unit (string)`
- `ligand_smiles (string)`
- `metadata (string JSON)`

Any extra columns are allowed for internal use but will be flagged by the validator.

## Plugins
Adapters can be extended by shipping packages that register an entry-point
group `rabiesmol.plugins` and expose a `DockingEngine` subclass.

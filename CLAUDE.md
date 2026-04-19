# CLAUDE.md

## Package structure

`src/` layout — the importable package lives at `src/gpaw_weaver/`:

| File | Contents |
|---|---|
| `params.py` | `make_pw_params`, `make_fd_params`, `get_mode_filestr` — pure functions, no I/O |
| `log.py` | `extract_scf_convergence` — parses GPAW SCF log text |
| `calculations.py` | `run_and_store_gpaw_calculation`, `load_gpaw_calculation` — DB-backed run/load |
| `__init__.py` | Re-exports the full public API |

All public symbols are importable flat: `from gpaw_weaver import make_pw_params`.

## Important conventions

- **`db` is an optional keyword argument** to both `run_and_store_gpaw_calculation` and `load_gpaw_calculation`. It accepts an ASE `Database` object, a `str`/`Path` file path (auto-connected), or `None` to use `calculations.db` in the working directory.
- `gpw_dir` and `gpw_logs` are optional `Path` parameters defaulting to `gpw_files/` and `gpw_logs/` relative to the working directory.
- `legacy_gpaw` (bool) selects the old vs. new GPAW implementation. It is stored in the database so `load_gpaw_calculation` can reload with the same implementation.
- Log files and GPW files are named after the **converged** DB row ID (not the initial one). The initial-ID log is renamed once the converged entry is written.
- `_serialize_calc_params` converts non-primitive values to JSON strings with `sort_keys=True` so DB key-value pairs are always consistent and queryable.

## Testing

Tests live in `tests/`. Run with `pytest` (configured in `pyproject.toml`).

| File | What it tests |
|---|---|
| `helpers.py` | `make_log()` and `make_fake_gpaw_class()` — shared test utilities |
| `conftest.py` | Shared fixtures: `fe_atom`, `pw_params`, `db`, `work_dirs` |
| `test_params.py` | Pure unit tests for `params.py` |
| `test_log.py` | Log parser for all three magnetism formats |
| `test_calculations.py` | Integration tests with mocked GPAW |

**All calculation tests are parametrized over `legacy_gpaw=True/False`** — add the `@LEGACY` decorator to any new calculation test.

`FakeGPAW` (from `helpers.make_fake_gpaw_class`) is an ASE `Calculator` subclass that:
- Writes fake SCF log content to the `txt` path on `__init__` (matching real GPAW behaviour)
- Returns canned energy/forces/magmoms from `calculate()`
- Implements `get_number_of_spins()` and `write()` as needed

## Development

```bash
pip install -e ".[dev]"
pytest
```

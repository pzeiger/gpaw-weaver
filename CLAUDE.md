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

### Function signatures

```python
run_and_store_gpaw_calculation(atoms_initial, calc_params, system,
                               db=None, save_gpw=False, save_gpw_mode='calculation',
                               legacy_gpaw=True, gpw_dir=..., gpw_logs=...)

load_gpaw_calculation(atoms_initial, system,
                      db=None, calc_params=None, legacy_gpaw=None, gpw_logs=...)
```

- `system` is a **required positional argument** to both functions.
- `atoms_initial` is also required in `load_gpaw_calculation` — its hash is used to identify the stored entry (see below).

### Database (`db`)

- `db` is an **optional keyword argument** accepting an ASE `Database` object, a `str`/`Path` file path (auto-connected), or `None` to use `calculations.db` in the working directory.
- Paths without a file extension are automatically given `.db`.

### Structure identification (`atoms_hash`)

Every stored entry carries an `atoms_hash` key-value pair computed by `_atoms_hash(atoms, calc_params)`:
- SHA-256 over atomic numbers, positions (rounded to 6 d.p.), cell, pbc, and magnetic moments.
- If `calc_params['magmoms']` is present it takes precedence over `atoms.get_initial_magnetic_moments()` for the magmom contribution — pass the same `calc_params` to `load_gpaw_calculation` to reproduce the hash.
- This ensures different phases, magnetic configurations (FM vs AFM), and non-collinear states stored under the same `system` label are always distinguished.

### GPAW implementation dispatch (`legacy_gpaw`)

- `legacy_gpaw=True` (default) uses `gpaw.calculator.GPAW` (old implementation).
- `legacy_gpaw=False` uses `gpaw.new.ase_interface.GPAW` (new implementation).
- Non-collinear calculations (`magmoms` as 3-vectors, `soc`) are **new GPAW only**.
- For old GPAW, `magmoms` from `calc_params` is applied via `atoms.set_initial_magnetic_moments()` (not the constructor). For new GPAW it is passed directly.
- `legacy_gpaw` is **never** passed to the GPAW constructor — it is a gpaw_weaver bookkeeping value stored in the DB so `load_gpaw_calculation` can reload with the correct implementation.

### ASE reserved keys

`_safe_db_key` prefixes any `calc_params` key that clashes with ASE reserved database column names (e.g. `magmoms`) with `gpaw_` before storage and querying. Both run and load go through `_serialize_calc_params`, so the renamed key is consistent on both sides.

### Other conventions

- `gpw_dir` and `gpw_logs` are optional `Path` parameters defaulting to `gpw_files/` and `gpw_logs/` relative to the working directory.
- Log files and GPW files are named after the **converged** DB row ID (not the initial one). The initial-ID log is renamed once the converged entry is written.
- `_serialize_calc_params` converts non-primitive values to JSON strings with `sort_keys=True` so DB key-value pairs are always consistent and queryable.
- `load_gpaw_calculation` raises `LookupError` when no match is found and `ValueError` when more than one converged row matches — narrow with `legacy_gpaw` or a more specific `system` label.

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
- Tests patch **both** `gpaw_weaver.calculations.GPAW` and `gpaw_weaver.calculations._NewGPAW` with `FakeGPAW`

## Development

```bash
pip install -e ".[dev]"
pytest
```

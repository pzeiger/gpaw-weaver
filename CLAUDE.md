# CLAUDE.md

## Backward compatibility (first-class requirement)

**Every change must keep existing databases readable and reproducible — a
newer gpaw-weaver must read and write an older DB with no migration.** This is
a primary design constraint for all future work here, not an afterthought:

- Treat `_calculation_hash` and `_serialize_calc_params` as a stable contract.
  Never change how existing keys are hashed/serialised; a stored `atoms_hash`
  must stay reproducible so old entries remain findable.
- Add capabilities *additively*: new parameters default to a no-op, and new
  behaviour activates only when a new key/argument is explicitly used, so
  every existing call site and stored row behaves exactly as before.
- If a change genuinely cannot preserve identity, it needs a deliberate,
  documented, versioned migration — never a silent hash shift.
- Guard the invariant with tests (the existing suite exercising the unchanged
  path is part of the proof).

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
run_and_store_gpaw_calculation(atoms_initial, calc_params,
                               db=None, label=None,
                               save_gpw=False, save_gpw_mode='calculation',
                               legacy_gpaw=True, parallel=None, vdw_factory=None,
                               gpw_dir=..., gpw_logs=...)

load_gpaw_calculation(atoms_initial,
                      db=None, calc_params=None, legacy_gpaw=None, gpw_logs=...)

query_gpaw_calculations(atoms_initial, db=None, calc_params=None)
```

- Both `atoms_initial` and `calc_params` are required positional arguments for `run_and_store_gpaw_calculation`.
- `atoms_initial` is required in `load_gpaw_calculation` — it is hashed to identify the stored entry.
- `label` is optional in both functions. It is stored in the DB for human readability but is **not** used as a query filter.
- `query_gpaw_calculations` returns a list of all `AtomsRow` objects whose stored structure (numbers, positions, cell, pbc, magmoms) and calc_params are a superset of the supplied values. Keys absent from the supplied `calc_params` are ignored (partial match). Returns both initial and converged rows.

### Database (`db`)

- `db` is an **optional keyword argument** accepting an ASE `Database` object, a `str`/`Path` file path (auto-connected), or `None` to use `calculations.db` in the working directory.
- Paths without a file extension are automatically given `.db`.

### Calculation identification (`atoms_hash`, `_calculation_hash`)

Every stored entry carries an `atoms_hash` key-value pair computed by `_calculation_hash(atoms, calc_params)`:
- SHA-256 over: atomic numbers, positions (rounded to 6 d.p.), cell, pbc, magnetic moments, and the full serialised `calc_params`.
- Magnetic moments: if `calc_params['magmoms']` is present it takes precedence over `atoms.get_initial_magnetic_moments()`.
- Because `calc_params` is included in the hash, different structures, phases, magnetic configurations (FM vs AFM, collinear vs non-collinear), and calculation settings (XC functional, k-points, …) all produce distinct hashes automatically.
- `calc_params` must be passed to `load_gpaw_calculation` to reproduce the correct hash. Omitting it hashes only the atomic structure and atoms magmoms.

`load_gpaw_calculation` raises `LookupError` when no match is found and `ValueError` when more than one converged row matches — narrow with `legacy_gpaw` or store distinct entries under different `label` values.

### van der Waals wrapper corrections (`vdw` / `vdw_factory`)

Two kinds of vdW correction:
- **Non-local `xc` functionals** (vdW-DF, vdW-DF-cx, …) need nothing special — they are just an `xc` value in `calc_params`, hashed and applied like any GPAW parameter.
- **Wrapper-style corrections** (Tkatchenko-Scheffler, DFT-D3, DFT-D4) are ASE calculators that *wrap* GPAW. weaver builds only a bare GPAW, so these use the `vdw` / `vdw_factory` pair:
  - `calc_params['vdw']` is a serialisable **descriptor** (e.g. `{'name': 'ts09', 'xc': 'PBE'}`). It carries the correction's **identity** — it is serialised, hashed, and stored like any other key, so a vdW run is distinct from the same calculation without it, and distinct per scheme.
  - `vdw_factory(dft_calc, atoms, descriptor)` is the **runtime mechanism** that turns the descriptor into a wrapper calculator around the freshly built GPAW. It is *not* serialised — weaver stays independent of any specific dispersion backend; the caller (project) owns the factory.
  - Required whenever `calc_params` has a `vdw` key (else `ValueError` — never a silent bare-GPAW run). The inner GPAW (`dft_calc`), not the wrapper, is used for the spin check and the `.gpw` restart.

**Backward compatibility (must be preserved by any future change here).** The `vdw` feature is strictly additive and **fully backward compatible in both directions** — a new gpaw-weaver reads *and* writes an old database with **no migration**:
- **No schema change.** ASE `db` rows hold arbitrary key-value pairs; old rows simply lack a `vdw` key. Nothing about the table changes.
- **Hash invariance.** `_calculation_hash` / `_serialize_calc_params` are **untouched**. `vdw` only enters the hash *when the key is present*, so any `calc_params` without a `vdw` key serialises byte-for-byte as before and produces the **identical** `atoms_hash` — old entries are still found by `load_gpaw_calculation` / `query_gpaw_calculations`, and a new run without vdW lands on the same hash an old weaver would have computed.
- **Inert unless used.** `vdw_factory` defaults to `None` and the wrapping branch runs only when `calc_params.get('vdw')` is set; every pre-existing call site and the entire pre-existing test suite hit exactly the old code path (the internal rename `calc`→`dft_calc` has no external effect).
- **Old weaver + new DB** also reads fine: an old weaver ignores the extra `vdw` kv, and non-vdW lookups never collide with vdW rows (distinct hashes). It merely cannot *re-run* a vdW entry (no hook) — never data corruption.

Do not alter `_calculation_hash`, `_serialize_calc_params`, or the "no `vdw` key ⇒ old behaviour" invariant without a deliberate, documented migration — existing databases depend on it.

### GPAW implementation dispatch (`legacy_gpaw`)

- `legacy_gpaw=True` (default) uses `gpaw.calculator.GPAW` (old implementation).
- `legacy_gpaw=False` uses `gpaw.new.ase_interface.GPAW` (new implementation).
- Non-collinear calculations (`magmoms` as 3-vectors, `soc`) are **new GPAW only**.
- For old GPAW, `magmoms` from `calc_params` is applied via `atoms.set_initial_magnetic_moments()` (not the constructor). For new GPAW it is passed directly.
- `legacy_gpaw` is **never** passed to the GPAW constructor — it is a gpaw_weaver bookkeeping value stored in the DB so `load_gpaw_calculation` can reload with the correct implementation.

### ASE reserved keys

`_safe_db_key` prefixes any `calc_params` key that clashes with ASE reserved database column names (e.g. `magmoms`) with `gpaw_` before storage. Both run and load go through `_serialize_calc_params`, so the renamed key is consistent on both sides.

### Other conventions

- `gpw_dir` and `gpw_logs` are optional `Path` parameters defaulting to `gpw_files/` and `gpw_logs/` relative to the working directory.
- Log files and GPW files are named after the **calculation hash** (same value stored as `atoms_hash` in the DB): `<gpw_logs>/<hash>.txt` and `<gpw_dir>/<hash>.gpw`. This makes them uniquely identifiable without a DB row ID and avoids any rename step.
- **DB entries are written only after the calculation converges successfully.** A failed or interrupted run leaves no partial rows in the database.
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
- Tests patch **both** `gpaw_weaver.calculations.GPAW` and `gpaw_weaver.calculations._NewGPAW` with `FakeGPAW`

## Development

```bash
pip install -e ".[dev]"
pytest
```

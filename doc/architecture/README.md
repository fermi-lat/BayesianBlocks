---
package: BayesianBlocks
summary: Standalone C++ Bayesian Blocks (Jackson, Scargle et al.) for unbinned events, binned data and point measurements, exposed to Python by SWIG, plus a pure-Python twin and old helper scripts.
layer: analysis
depends_on: []
used_by: []
last_reviewed_commit: 3f1404a
---

# BayesianBlocks Architecture

**BayesianBlocks** is one C++ class (about 200 lines in
`src/BayesianBlocks.cxx`) that finds the optimal piecewise-constant
segmentation of 1-D data with the dynamic-programming algorithm of
Jackson, Scargle et al. (2005) and the false-positive `ncp_prior`
calibration of Scargle et al. (2013). It is meant to be used from Python:
the SWIG module `BayesianBlocks` (built as `_BayesianBlocks`) is the
package's real interface. The package also installs a pure-Python
reimplementation and some 2011–2013 analysis scripts.

There are **no `gt*` tools** and no tests that run in CI.

`src/mainpage.h` (a short description of the algorithm with dead
`trotsky.arc.nasa.gov` links) and `doc/release.notes` (last entry
2013) are accurate as far as they go, but they do not describe the Python
side.

## Purpose and scope

**Does:**
- Partition data into blocks in three modes (chosen by constructor):
  - **Unbinned events.** Each event gets a Voronoi cell between the
    midpoints to its neighbours, with the first and last cells closed at
    the given `tstart`/`tstop`.
  - **Binned data.** Bin contents and bin widths.
  - **Point measurements.** Values *y ± dy* at positions *x*.
- Optionally replace the cell "sizes" used in the cost function with
  other weights (`setCellSizes`), for example per-cell exposure.
- Compute `ncp_prior` from the number of events and a false-positive
  fraction.

**Does not:**
- Read or write files. The C++ code uses only the standard library.
- Use any IRF, exposure or likelihood code, even though CMake links them
  (see below). The exposure weighting in `bb_analysis.py` comes from a
  `gtexposure` output *file*.
- Share code with `evtbin::BayesianBinner` (used by `gtburstfit`; see
  `src/evtbin/doc/architecture/` and `src/burstFit/doc/architecture/`) or
  with `pyBurstAnalysisGUI`'s pure-Python `GtBurst.BayesianBlocks`
  (used by `gtburst`). These are **three independent implementations**
  of the same algorithm. They differ in their priors: evtbin adds a
  hidden 6.7 to `ncp_prior`, and GtBurst uses its own calibration.
  Results are not interchangeable.

## Position in the system

| Direction | Package | Why |
|-----------|---------|-----|
| Link only | `Likelihood`, `irfs` (`irfInterface`), `astro`, `st_facilities` | `target_link_libraries(BayesianBlocks PUBLIC astro st_facilities Likelihood irfInterface)`, but `BayesianBlocks.h`/`.cxx` include only `<deque>`, `<vector>`, `<cmath>`, `<numeric>`, `<stdexcept>`, `<iostream>` (and `<fenv.h>`). No Fermi symbol is used. These links are probably left over from the SCons build. Because they are `PUBLIC`, anything that linked `BayesianBlocks` would also pull in all of `Likelihood`. |
| Depends on (external) | Python 3, SWIG | `py_BayesianBlocks` links `Python3::Python`. |
| Runtime (scripts only) | `evtbin`, `pyExposure`, `sane` | `create_lc_file.py` runs `gtbindef` and `gtbin` (evtbin) and `gtexposure` (pyExposure) through `sane`'s `GtApp`. These are subprocess calls, not code dependencies. |
| Used by | — | No package includes `BayesianBlocks/` or imports the `BayesianBlocks` module. |

So the package's real dependencies are none. Its place in the
`analysis` tier, and its depth 7 in `docs/architecture/dependency-layers.md`,
come entirely from the unused links. If those links were removed, it
could build right after the externals.

## Key abstractions

| Type | Role |
|------|------|
| `BayesianBlocks` (global namespace) | Holds the cells: `m_cellSizes`, `m_cellContent`, `m_cellErrors` (point mode), and prefix sums (`std::deque`, with a leading 0) of sizes and contents, so that `blockSize(i,j)`/`blockContent(i,j)` are O(1). `globalOpt(ncp_prior, xx, yy)` runs the O(N²) dynamic program, backtracks the change points and fills the step-function light curve. |
| `BlockCost` (private, abstract) | `operator()(imin, imax)` → log-fitness of one block of cells `imin..imax` (inclusive). |
| `BlockCostEvent` | For events and binned data: `N·(ln(N/T) − 1)`, where N = block content and T = block size (scaled). An empty block costs 0. The `−N` term sums to a constant over any partition, so it does not change the optimum. |
| `BlockCostPoint` | For point measurements: `−½·Σw · Var_w(y)` with w = 1/dy². It loops over the cells in the block, so point mode is O(N³) overall. |

Each block costs `fitness − ncp_prior`. `ncp_prior(N, p)` =
`4 − ln(p / (0.0136·N^0.478))`, which is eq. 21 of Scargle et al. 2013
(0.0136 = 1/73.53).

Two size sequences are kept on purpose:
- `m_cellSizePartialSums` holds the **scaled** sizes (the cost uses them,
  and `setCellSizes` replaces them).
- `m_unscaledCellSizePartialSums` holds the **original** widths. The
  output x-coordinates and y = content / unscaled width use them.

So after `setCellSizes(exposure)` the block *boundaries* are
exposure-weighted, but the returned *y* is still counts per unit time,
not flux.

## Data flow and control flow

```python
import BayesianBlocks                          # SWIG module (fermitools site-packages dir)
bb = BayesianBlocks.BayesianBlocks(times, tstart, tstop)   # unbinned; times must be sorted
# bb = BayesianBlocks.BayesianBlocks(tstart, counts, widths)   # binned
# bb = BayesianBlocks.BayesianBlocks(x, y, dy)                 # point measurements
bb.setCellSizes(exposures)                     # optional; len == number of cells
p = BayesianBlocks.BayesianBlocks.ncp_prior(len(times), 1e-3)
xx, yy = bb.lightCurve(p)                      # %extend wrapper around globalOpt
# xx = (b0, b1, b1, b2, ...), yy = (r0, r0, r1, r1, ...): two points per block
```

The script workflow (`python/`), from 2013:
1. `create_lc_file.py evfile scfile lcfile irfs`: writes one time bin per
   event's Voronoi cell (`gtbindef`), bins the FT1 file into a light curve
   with those bins (`gtbin`), then adds an `EXPOSURE` column (`gtexposure`).
2. `bb_analysis.py evfile lcfile`: runs unbinned BB with
   `setCellSizes(lc.EXPOSURE)` and `fp_frac=1e-3`. It writes
   `change_points.txt` (x, y, counts/exposure per block).

**Neither script runs in a current Fermitools environment**: both import
`pyfits` (through `FitsNTuple.py`), which is not installed (confirmed at
runtime with Fermitools 2.5.1: `ModuleNotFoundError: No module named
'pyfits'`).

## Public interfaces

### Python (installed to `${FERMI_INSTALL_PYTHON}` = `site-packages/fermitools`)

| Module | Status |
|--------|--------|
| `BayesianBlocks` (+ `_BayesianBlocks` extension) | Works (confirmed at runtime for all three modes). Exposes the class, `DoubleVector`-style STL templates, `lightCurve(ncp_prior)` → `(xx, yy)` tuples, static `ncp_prior` (also flattened as `BayesianBlocks_ncp_prior`), `enableFPE` (`BayesianBlocks_enableFPE`). C++ exceptions become `RuntimeError`. |
| `BayesianBlocks_python` | Pure-Python version of the same algorithm, used as a reference. Its **binned constructor argument order differs**: `(content, sizes, tstart)` vs. C++ `(tstart, content, sizes)`. Its unbinned constructor takes only the events: it sorts them itself, and closes the outer cells by reflecting the neighbouring midpoint instead of using `tstart`/`tstop`. |
| `bb_analysis`, `create_lc_file`, `FitsNTuple` | Broken: need `pyfits`. `FitsNTuple.makeNTuple` also needs HippoDraw and Python 2 (`xrange`, `sys.maxint`). |
| `test_BB`, `test_BB_speed` | Interactive comparison and timing scripts. They need `pylab_plotter`, `distributions` and `read_data`, none of which ship, and `test_BB_speed` hard-codes `/u/gl/jchiang/python`. They are installed into site-packages anyway. |

### C++ library

`BayesianBlocks/BayesianBlocks.h`, a SHARED library `libBayesianBlocks`
(alias `Fermitools::BayesianBlocks`). The class is in the global
namespace, not in `namespace BayesianBlocks`. Nothing in ScienceTools
links it except the SWIG module.

### Environment variables

None.

### File formats

The library reads and writes no files. The scripts read FT1 (`EVENTS`
`TIME`, `TSTART`/`TSTOP`, DSS `ENERGY` range) and a `gtbin` LC file with
`EXPOSURE`, and write a text `change_points.txt`.

## Build and test

- **CMake targets:** `BayesianBlocks` (SHARED) and `py_BayesianBlocks`
  (SWIG, module name `_BayesianBlocks`, with warnings 362/389/401/402/503/509
  suppressed).
- **Install:**
  - headers go to `include/BayesianBlocks`
  - `python/*` goes to `site-packages/fermitools`
  - the generated `BayesianBlocks.py` goes to `site-packages/fermitools`
    (the same `install(FILES …)` appears **twice** in `CMakeLists.txt`;
    harmless)
  - there is no `data/` and no `pfiles/`
- **Tests:** none run automatically. There is no C++ test program, no
  CTest registration, and no entry in `recipe/tests/ST-unit-test`.
  `python/test_BB.py` compares the C++ and Python implementations
  visually, but cannot run (see above). Nothing checks this package's
  results in CI. A quick manual check:
  `python -c "import BayesianBlocks as B; print(B.BayesianBlocks([0.1,0.2,0.5,0.9],0.,1.).lightCurve(1.))"`.
- **Legacy:** `SConscript`, `BayesianBlocksLib.py`, `doc/Doxyfile`
  (points at a nonexistent `../cmt`) and `azure-pipelines.yml` are
  unused by the CMake build.
- **Branches:** `origin/f26-cpp17-compliance` (trailing return types
  only), `origin/py-stable-abi` (SWIG link flags) and `origin/Issue17`
  (optionally static libraries) are unmerged. None of them changes
  behaviour.

## Invariants and gotchas

**Input requirements (no checks)**
- **Unbinned input must be sorted.** The C++ constructor does not sort
  (the Python twin does). Unsorted times give negative cell sizes and a
  silently wrong segmentation (confirmed at runtime: 400 events in 3 true
  blocks gave 9 blocks when shuffled).
- **Unbinned and point modes need at least 2 points.**
  `generateCells` reads `arrival_times[1]` and `[npts-2]` without
  checking. With one event this is undefined behaviour. At runtime it
  returned a block `(0, 1.5)` for `tstart=0, tstop=2` without an error
  (confirmed at runtime).
- `tstart`/`tstop` must enclose the events. Otherwise the first or last
  cell has a negative size.
- `setCellSizes` throws if its length differs from the number of cells,
  which is the number of events in unbinned mode.
- Point mode needs `dy > 0`. A zero error gives an infinite weight.

**Memory and copying**
- `m_blockCost` is allocated with `new` and never deleted (there is no
  destructor), so every `BayesianBlocks` object leaks one small object.
  The implicit copy constructor copies the pointer. A copy's cost
  functor therefore still refers to the **original** object's data, and
  dangles once the original is destroyed. SWIG does not copy, so Python
  users are not affected.
- In binned mode, `m_tstop` is never initialized. It is not read in that
  mode.

**Semantics**
- `lightCurve`/`globalOpt` return **two points per block** (step-plot
  form), not edges and values.
- After `setCellSizes`, *y* is still content / unscaled width, as
  described above. `bb_analysis.py` computes the flux separately from the
  LC file's `EXPOSURE`.
- `enableFPE()` turns on `FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW` for the
  **whole process** (Python included). On macOS it throws, because
  `TRAP_FPE` is defined only for non-Apple builds (top-level
  `CMakeLists.txt`).
- The algorithm is O(N²) time and memory per call (O(N³) in point mode).
  For tens of thousands of events, expect minutes.

## Extension points

- **New data mode / fitness function:** add a `BlockCost` subclass and a
  constructor that sets `m_blockCost` to it and fills
  `m_cellSizes`/`m_cellContent`. Then call `cellPartialSums()`.
  `lightCurve` has a point-mode branch that you may need to extend.
- **Python additions** go in the `%extend` block of
  `src/BayesianBlocks.i`.
- If you fix the scripts, replace `pyfits` with `astropy.io.fits` and
  `FitsNTuple` with direct table reads. `GtApp` is still provided by
  `sane`.
- Before adding features, consider whether this package or
  `evtbin::BayesianBinner` should be the single shared implementation.
  The CMake links to `Likelihood`/`irfs`/`astro`/`st_facilities` can be
  dropped without code changes (from reading the code; not built).

## Further reading

- B. Jackson, J. D. Scargle et al. 2005, *An algorithm for optimal
  partitioning of data on an interval*, IEEE SPL 12, 105
  (arXiv:math/0309285). This is the reference cited in
  `BayesianBlocks_python.py`.
- J. D. Scargle et al. 2013, ApJ 764, 167, *Studies in Astronomical Time
  Series Analysis VI* (the `ncp_prior` calibration).
- `src/evtbin/doc/architecture/` (`BayesianBinner`) and
  `src/burstFit/doc/architecture/`.
- `src/pyBurstAnalysisGUI/python/GtBurst/BayesianBlocks.py`, the
  implementation `gtburst` actually uses.

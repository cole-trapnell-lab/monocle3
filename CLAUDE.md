# monocle3

R package v1.4.29: clustering, DE, and trajectory analysis — the `cell_data_set` (CDS) core the
whole stack is built on. Upstream of hooke, platt, and zscapetools.

## Commands
**There is no Makefile here.** Root `make fast` prints `(no Makefile; skipping)` and skips this repo.
```bash
Rscript -e 'testthat::test_local(".", reporter="summary")'   # ~24 test files
Rscript -e 'testthat::test_file("tests/testthat/test-fit_models.R")'   # single file
R CMD check --no-manual .                                    # the real gate; slow
```

## Gotchas
- **Tests branch on a `TRAVIS` env var** to select expected target values and skip interactive
  code. CI sets `TRAVIS: true`. Results can differ with and without it.
- `tests/testthat.R` uses `test_check("monocle3")`, which needs the package **installed**;
  `test_local(".")` is the loop that works against the working tree.
- `LinkingTo: Rcpp` + `src/` C++ — edits there need a recompile before tests mean anything.
- ~60 `Imports:` including Bioconductor (`batchelor`, `HDF5Array`, `limma`, `S4Vectors`) and
  heavy system-dependent packages (`sf`, `spdep`, `leidenbase`, `BPCells`). Adding a dependency
  here is expensive for every downstream repo.
- `00travis.yml` is retired notes, not live config.

## Release notes
`NEWS.md` is the stack's reference format: `# Monocle3 <version>` / `### Changes` / one bullet per
user-visible change, newest section first. **It currently lags `DESCRIPTION`** — DESCRIPTION is at
1.4.29, NEWS.md's newest entry is 1.4.25. Write the NEWS entry in the same commit as the bump.

## CI
`.github/workflows/check_on_push.yml`, `on: [push]`, `R CMD check` in
`ghcr.io/cole-trapnell-lab/monocle3_depend`. `.github/CODEOWNERS` gates workflow edits.

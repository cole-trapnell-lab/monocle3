# Codex instructions (Hooke)

## What this repo is
Hooke is an R package used in the Trapnell Lab stack for statistical analysis of single-cell perturbation data.

# Codex instructions (Hooke)

Hooke is an R package used in the Trapnell Lab stack for statistical analysis of single-cell perturbation data.

## Cost & safety guardrails
- Safe to run R package checks and unit tests (minutes).
- Do NOT run long benchmarks or full-dataset analyses unless explicitly requested.
- Never commit secrets, credentials, or private datasets.

## Change expectations
- Prefer small PRs (<300 LOC changed when possible).
- If changing exported functions, update:
  - roxygen docs (if used)
  - tests
  - NEWS.md (if present)
- If changing object structures or outputs, add a backward-compat test or a migration note.
- When adding dependencies, keep them minimal and declare them in `DESCRIPTION`; avoid heavy, non-CRAN dependencies unless justified.

## Validation tiers
- **FAST (default; safe):** `make fast` (runs `testthat::test_local` with summary reporter; should finish in minutes and not build the package).
- **SMOKE (opt-in):** Small toy example / vignette using tiny fixtures only.
- **FULL (manual only):** Real dataset runs / benchmarks.

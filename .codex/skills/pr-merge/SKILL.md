---
name: pr-merge
description: Review branch changes against main before merge, ensure docs and changelog updates are present, and report merge-readiness warnings.
---

# PR Merge Readiness

Use this skill when the user wants a final pre-merge pass for a feature branch.

## Goals

- Compare branch diff against `main`.
- Confirm user-visible changes have matching docs updates.
- Add or refine `CHANGELOG.md` entries under `Unreleased`.
- Surface warnings/risks before merge.

## Workflow

1. Gather diff scope:
   - `git diff --name-status main...HEAD`
   - `git diff --stat main...HEAD`
2. Classify changed files:
   - Core solver (`source/`), bindings, tests, docs, CI, packaging.
3. Docs consistency check:
   - If behavior/API changed, verify `README.md` and/or `docs/source/` updates exist.
   - If docs are missing, flag as warning and propose exact files to update.
4. Changelog update:
   - Add concise user-visible bullets under `## [Unreleased]` in `CHANGELOG.md`.
   - Use sections already present (`Changed`, `Fixed`, `Added`).
5. Sanity checks:
   - Confirm no accidental edits to generated docs/build outputs.
   - Highlight risky areas: missing tests, interface changes, scientific/numerical behavior impacts.
6. Produce merge recommendation:
   - `Ready`, `Ready with warnings`, or `Not ready`.

## Output Format

- Findings first, ordered by severity, with file references.
- Proposed `CHANGELOG.md` additions.
- Final merge-readiness verdict and blockers.

## Guardrails

- Do not merge branches automatically.
- Do not rewrite changelog history outside `Unreleased` unless user asks.
- Keep warnings concrete and evidence-based.

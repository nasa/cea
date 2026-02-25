---
name: docs-cleanup
description: Audit repository documentation for stale, conflicting, or incomplete install/run guidance and outdated examples; report concrete file-level gaps and recommended fixes.
---

# Docs Cleanup

Use this skill when the user asks for a documentation audit, cleanup pass, or stale-doc review.

## Goals

- Find conflicting or outdated install/run instructions.
- Check whether examples still match current CLI/API behavior.
- Flag missing guidance that blocks onboarding, testing, or release tasks.
- Keep output focused on actionable findings with file references.

## Workflow

1. Collect candidate docs:
   - `README.md`
   - `docs/source/**/*.rst`
   - `CONTRIBUTING.md`
   - `.github/` templates/workflows when relevant to user instructions.
2. Validate install/run instructions against current repo behavior:
   - Build and test commands in docs vs commands in `CMakeLists.txt`, `pyproject.toml`, CI workflows, and `scripts/`.
   - Confirm paths, filenames, and prerequisites are still valid.
3. Validate examples:
   - Compare documented examples with `samples/`, `test/`, and current interfaces.
   - Check that command arguments and expected outputs are still plausible.
4. Identify gaps:
   - Missing prerequisites, platform caveats, troubleshooting notes, or release-process notes.
5. Report findings first, ordered by severity:
   - `High`: wrong instruction likely to fail.
   - `Medium`: confusing/outdated but recoverable.
   - `Low`: clarity/consistency improvements.

## Output Format

- Findings list with `severity`, `file`, and exact issue.
- Proposed fix for each finding.
- Open questions or assumptions.
- If no findings: explicitly state that and note residual risk areas not fully validated.

## Guardrails

- Prefer small documentation edits over broad rewrites.
- Avoid speculative claims; cite the source file/command that contradicts the docs.
- Do not edit generated docs under `docs/_build/` or `docs/doxygen/`.

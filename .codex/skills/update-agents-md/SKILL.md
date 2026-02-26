---
name: update-agents-md
description: Update AGENTS.md to reflect recent repository workflow, architecture, and policy changes while preserving existing project constraints.
---

# Update AGENTS.md

Use this skill when AGENTS guidance has drifted from the codebase or team workflow.

## Goals

- Keep `AGENTS.md` current with architecture boundaries, testing expectations, and repo workflows.
- Preserve critical existing constraints unless the user explicitly asks to change them.

## Workflow

1. Read current `AGENTS.md` fully.
2. Scan recent project changes for process-impacting updates:
   - `README.md`, `CONTRIBUTING.md`, `CHANGELOG.md`
   - CI workflows in `.github/workflows/`
   - Build/package entry points (`CMakeLists.txt`, `pyproject.toml`)
   - Interface boundaries under `source/`.
3. Propose targeted edits only where guidance is stale or missing.
4. Apply minimal, focused wording changes:
   - Keep rules concrete and enforceable.
   - Avoid aspirational or ambiguous language.
5. Validate consistency:
   - No contradictions with docs or CI.
   - Keep scientific correctness and compatibility priorities explicit.

## Output Format

- Summary of what changed in `AGENTS.md`.
- Why each change was needed.
- Any deferred items requiring user decision.

## Guardrails

- Do not remove existing hard constraints unless explicitly instructed.
- Avoid adding style-only or noisy policy text.
- Keep edits short and maintainable.

---
name: version-increment
description: Prepare the repository for a new release tag by updating hard-coded version references, promoting CHANGELOG Unreleased entries, and running release-readiness checks.
---

# Version Increment

Use this skill when the user provides a new version number and wants the repo prepared for tagging.

## Inputs

- `new_version` (required), semver-like string (example: `3.0.4`).
- `release_date` (optional). If omitted, use local current date (`YYYY-MM-DD`).

## Known Hard-Coded Version Locations (current repo)

- `CMakeLists.txt` (`project(... VERSION x.y.z)`).
- `source/bind/python/cea/__init__.py` (`__version__ = "x.y.z"`).
- `docs/source/conf.py` (`release = '...'`).

Update these locations and verify no additional hard-coded version literals remain outside generated docs.

## Workflow

1. Validate `new_version` format and ensure it is greater than the current version.
2. Update the three known hard-coded version locations.
3. Search for additional stale literals:
   - Prefer `rg -n "<old_version>|release\s*=|__version__\s*=|VERSION\s+[0-9]+\.[0-9]+\.[0-9]+"`
   - Exclude generated docs: `docs/_build/`, `docs/doxygen/`.
4. Update `CHANGELOG.md`:
   - Move current `## [Unreleased]` entries into a new release section `## [<new_version>] - <release_date>`.
   - Recreate an empty `## [Unreleased]` section at the top with standard subsections used by this repo.
5. Run sanity checks:
   - Ensure changed files are limited to release metadata/docs unless user asked for more.
   - Optionally run targeted tests/docs build if requested.
6. Report exact files changed and any warnings before tagging.

## Output Format

- New version and date used.
- List of updated files with one-line reason each.
- Warnings/blockers before tag creation.
- Suggested tag command: `git tag v<new_version>` (only if user asks to tag).

## Guardrails

- Do not create/push tags unless explicitly requested.
- Do not edit generated docs under `docs/_build/` or `docs/doxygen/`.
- Keep changelog wording user-visible and avoid implementation-only noise.

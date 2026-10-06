# SAGE publication projects

The SAGE GUI discovers every immediate subfolder that contains a
`project.json` file. Click **Refresh** on the Projects tab after adding or
editing a project.

Each project may contain:

- `project.json` — title, publication reference, notes, status, and links.
- `setup.json` — configuration values to merge over the user's current GUI
  configuration. Machine-specific data and result roots are preserved.
- a canonical `run_SAGE_export.m`, referenced by `source_export`, when one
  is available.

Use `"kind": "sage"` for a normal SAGE setup and `"kind": "standalone"`
for a paper with a separate analysis pipeline. A SAGE project becomes
loadable when its status is `ready` and its setup file exists.

This data-only format means browsing a project never executes MATLAB code.

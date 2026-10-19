# Tooling used to produce the AI-edited manuscript drafts (2026-10-03)

**Note (2026-10 revision):** these scripts predate the decision to remove Figure 1c, the "super-additive" / interaction-information calibration and the within-class negative control from the paper; they still contain that wording and do not reproduce the current drafts. They are kept only for traceability of the earlier drafts.

Not part of the analysis. Kept for traceability of how `NPJ_review/MS_final_reviewed_AI_revision_2026-10-03.docx`
and `Reviewer_comments_AI_responses_2026-10-03.docx` were produced; paths inside are specific to the machine and
session where they ran (the docx skill scripts and a scratch directory).

- `docx_edit.py` – tracked-change engine: purple (`7030A0`) `<w:ins>` insertions, `<w:del>` deletions, Word comments.
- `ms_numbers.py` – reads every number from the repository's output files (no number is typed by hand).
- `apply_manuscript_edits.py` / `apply_reviewer_responses.py` – the edits, each with the comment explaining it.
- `write_revision_changes.py` – generates `manuscript/REVISION_CHANGES.md`.
- `finalize.sh` – numbers index → supplementary files → both docx drafts (validated against the originals).
- `overnight_rerun.sh` – the driver of the 2026-10-02/03 rerun (phases A–D).

Note (2026-10-05): the supplementary tables were renumbered by order of first citation after these scripts
were run (old -> new: 5->1, 1->2, 8->3, 2->4, 9->5, 3->6, 4->7, 6->8, 7->9). The scripts here and
`manuscript/TEXT_EDITS.md` still use the old numbers; the repository files, notebooks, METHODS.md and the
supplementary workbook use the new ones.

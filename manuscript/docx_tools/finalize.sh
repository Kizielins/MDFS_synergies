#!/usr/bin/env bash
# After the overnight driver: numbers index, manuscript + reviewer docx drafts, supplementary files, previews.
set -uo pipefail
W="$(cd "$(dirname "$0")" && pwd)"
REPO="/Users/kizie/Library/Mobile Documents/com~apple~CloudDocs/OneDrive-UJ/Sonata2/mdfs_combine_taxonomy_and_functions/CRC_segata_pooling/framework_repository"
NPJ="/Users/kizie/Library/Mobile Documents/com~apple~CloudDocs/OneDrive-UJ/Sonata2/NPJ_review"
SK="/Users/kizie/Library/Application Support/Claude/local-agent-mode-sessions/skills-plugin/8879dd4d-1049-447a-89bb-74aba2ceb3dd/36e6cf99-bd73-4f25-b6e1-acaca24f5810/skills/docx/scripts"
STAMP=$(date +%Y-%m-%d)
cd "$REPO" && python3 manuscript/collect_numbers.py > /dev/null && echo "numbers index ok"
python3 "$W/write_revision_changes.py" > /dev/null && echo "REVISION_CHANGES.md ok"
python3 manuscript/build_supplementary_files.py && echo "supplementary files ok"
cd "$W"
# manuscript
rm -rf ms_final && cp -r un ms_final && python3 apply_manuscript_edits.py ms_final | head -1 && (cd ms_final && rm -f ../MS_AI.docx && zip -q -Xr ../MS_AI.docx .) \
 && python3 "$SK/office/validate.py" MS_AI.docx --original orig.docx --author "Claude (AI assistant)" | tail -1 \
 && cp MS_AI.docx "$NPJ/MS_final_reviewed_AI_revision_$STAMP.docx" && echo "manuscript draft -> NPJ_review/MS_final_reviewed_AI_revision_$STAMP.docx"
# reviewer responses
rm -rf rc_final && cp -r rc_un rc_final && python3 apply_reviewer_responses.py rc_final | tail -1 && (cd rc_final && rm -f ../RC_AI.docx && zip -q -Xr ../RC_AI.docx .) \
 && python3 "$SK/office/validate.py" RC_AI.docx --original rc_orig.docx --author "Claude (AI assistant)" | tail -1 \
 && cp RC_AI.docx "$NPJ/Reviewer_comments_AI_responses_$STAMP.docx" && echo "responses draft -> NPJ_review/Reviewer_comments_AI_responses_$STAMP.docx"
cp "$REPO/manuscript/supplementary/Supplementary_tables.xlsx" "$NPJ/Supplementary_tables_AI_revision_$STAMP.xlsx"
cp "$REPO/manuscript/supplementary/Supplementary_Figures.docx" "$NPJ/Supplementary_Figures_AI_revision_$STAMP.docx"
pandoc --track-changes=accept -t plain --wrap=none MS_AI.docx > MS_AI_accepted.txt; pandoc --track-changes=accept -t plain --wrap=none RC_AI.docx > RC_AI_accepted.txt
echo "finalize done"

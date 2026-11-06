#!/usr/bin/env bash
# Runs every analysis for the revised manuscript, in order. Logs: outputs/logs/<script>.log.
# Stops at the first failure. Needs data_raw/CF_registry_master.xlsx (not included; see README.md).
cd "$(dirname "$0")"
if [ ! -f data_raw/CF_registry_master.xlsx ]; then echo "Missing data_raw/CF_registry_master.xlsx (see README.md)"; exit 1; fi
L=outputs/logs; mkdir -p "$L"
# gt writes .docx tables through pandoc; use RStudio's copy on macOS if pandoc is not on the PATH
if ! command -v pandoc >/dev/null 2>&1 && [ -z "$RSTUDIO_PANDOC" ]; then
  for d in /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64 \
           /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64 \
           /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools; do
    if [ -x "$d/pandoc" ]; then export RSTUDIO_PANDOC="$d"; break; fi
  done
fi
scripts=("A01.R" "A02.R" "Supplementary Scripts/Build_Exac_Anchored_Tidies.R" "A03.R" "A03.1.R" "A04.R" "A05.2.R"
         "A08_revision_fixes.R" "A09_revision_outputs.R" "A10_fev1_reference_check.R" "A12_revision_sensitivity.R"
         "A13_revision_figures.R" "A14_sampling_sensitivity.R" "A15_exacerbations.R" "A16_figure1.R"
         "A18_chronic_status_unknown.R" "A21_revision_consistency.R" "A22_figureS2_infections.R" "A23_table1_extra.R"
         "A19_column_audit.R" "A24_age_check.R")
for s in "${scripts[@]}"; do
  b=$(basename "$s" .R); t0=$(date +%s)
  echo "$(date '+%H:%M:%S') start $s"
  if ! Rscript "$s" > "$L/$b.log" 2>&1; then echo "$(date '+%H:%M:%S') FAILED $s (see $L/$b.log)"; exit 1; fi
  echo "$(date '+%H:%M:%S') done  $s ($(( $(date +%s) - t0 )) s)"
done
echo "$(date '+%H:%M:%S') ALL DONE"

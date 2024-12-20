#!/bin/bash

# Define the base path
BASE_PATH="reconstruction/alert/src/main/java/org/jlab/rec/atof"

# Directories to process
DIRECTORIES=(
  "ML_INPUT"
  "BAR_WEDGE_CLUSTERING"
  "ATOF_RECON"
  "ATOF_RECON_CLUSTER_BANK"
  "ATOF_RECON_CLUSTERING"
  "ATOF_RECON_ZPHITIME"
  "BarWedgeClusters"
  "Bar_Only_Clusters"
  "RESIDUALS_CALCULATIONS"
)

# Process each directory
for dir in "${DIRECTORIES[@]}"; do
  FULL_PATH="$BASE_PATH/$dir"

  # Check if .git exists and rename it
  if [ -d "$FULL_PATH/.git" ]; then
    echo "Renaming .git to .git_backup in $FULL_PATH"
    mv "$FULL_PATH/.git" "$FULL_PATH/.git_backup"
  else
    echo "No .git found in $FULL_PATH"
  fi
done

# Stage the changes in Git
git add "$BASE_PATH"

# Commit the changes
git commit -m "Renamed .git directories to .git_backup for all subdirectories"

# Push the changes
git push origin development

echo "All directories processed and changes pushed!"

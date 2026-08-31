#########
# Master Script to automate submitting jobs in BAF,
# Automate File Handling 
# with a repo_name based standard folder structure
#########



#!/usr/bin/env bash
set -euo pipefail


# global automated parameters:

repo="mbl_dtc"   





################################ 
# --- sanity checks ---
################################

if [[ ! -d "$repo" ]]; then
  echo "ERROR: repo directory not found: $repo" >&2
  exit 1
fi

if [[ ! -f "./$repo/cat_env.sh" ]]; then
  echo "ERROR: missing or non-executable: $repo/cat_env.sh" >&2
  exit 1
fi

if [[ ! -f "./$repo/cat_run.sh" ]]; then
  echo "ERROR: missing or non-executable: $repo/cat_run.sh" >&2
  exit 1
fi


if [[ ! -f "./$repo/cat_jl.sh" ]]; then
  echo "ERROR: missing or non-executable: $repo/cat_jl.sh" >&2
  exit 1
fi



################################ 
# --- Git Sync ---
################################

git pull

cd $BUDDY/BAF/
git pull

cd ~/local/BAF/

################################ 
# --- Script Execution ---
################################


bash ./$repo/cat_run.sh
bash ./$repo/cat_env.sh
bash ./$repo/cat_jl.sh

bash ./$repo/condor_submit.sh


echo "Done."

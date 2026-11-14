#!/bin/bash
# removes from the condor queue the jobs of one job list (e.g. jobs_v2.1_L14.txt), and only those:
#   ./rm_jobs.sh jobs_v2.1_L14.txt            show what would be removed, then ask
#   ./rm_jobs.sh jobs_v2.1_L14.txt --idle     only the idle (not yet started) jobs
#   ./rm_jobs.sh jobs_v2.1_L14.txt --yes      no question (scripts)
# A job belongs to the list when it was submitted from a batch directory ~/local/batches/mbldtc_<set>_L<L>_theta<th>_<stamp>
# (cat_batch.sh), with <set>, <L>, <th> from the lines of the list (all stamps of that set/L/theta: earlier submissions too).
# So jobs_test_v21.txt (L=14, theta=0, v21) also matches the production L=14 theta=0 v21 batches: read the per batch
# counts it prints before typing yes. A list line without a set name is set polfed (common.sh).
# The JobBatchName ("mbldtc_polfed: L=.., theta=..") is the same for every set, so it is NOT used.
# condor_q/condor_rm talk to the scheduler of this machine; run it where you submitted.

LIST=$1; shift
[[ -f "$LIST" ]] || { echo "usage: $0 <job list> [--idle] [--yes]" >&2; exit 1; }
IDLE=0; YES=0
for a in "$@"; do case $a in --idle) IDLE=1;; --yes) YES=1;; *) echo "unknown option $a" >&2; exit 1;; esac; done

# one pattern per (set, L, theta) of the list: L = column 1, theta = column 2, set = the first non-numeric column
mapfile -t PATS < <(grep -v '^\s*#' "$LIST" | awk 'NF>=3 {s=""; for(i=3;i<=NF;i++) if($i !~ /^[0-9.]+$/){s=$i; break}; if(s=="") s="polfed"; print "mbldtc_" s "_L" $1 "_theta" $2 "_"}' | sort -u)
(( ${#PATS[@]} )) || { echo "no job lines in $LIST" >&2; exit 1; }
C=""
for p in "${PATS[@]}"; do C+="${C:+ || }regexp(\"/${p}\", Iwd)"; done
C="($C) && Owner == \"$USER\""
(( IDLE )) && C="$C && JobStatus == 1"

echo "job list $LIST -> batch directories ${PATS[*]}*"
echo "constraint: $C"
# what is in the queue for it, by batch and state (1 idle, 2 running, 5 held)
condor_q -constraint "$C" -af Iwd JobStatus | awk '{n=split($1,a,"/"); s[a[n]" "($2==1?"idle":$2==2?"running":$2==5?"held":"state"$2)]++} END {for (k in s) print "  " k ": " s[k]}' | sort
N=$(condor_q -constraint "$C" -af ClusterId | wc -l)
(( N )) || { echo "nothing in the queue for it"; exit 0; }
R=$(condor_q -constraint "$C && JobStatus == 2" -af ClusterId | wc -l)
(( R )) && echo "WARNING: $R of them are RUNNING (their compute so far is lost); --idle removes the idle ones only"
if (( ! YES )); then
  read -r -p "remove these $N jobs? type yes: " ans
  [[ "$ans" == "yes" ]] || { echo "nothing removed"; exit 0; }
fi
condor_rm -constraint "$C"
echo "left in the queue for it: $(condor_q -constraint "$C" -af ClusterId | wc -l)"

#!/bin/bash
# Produces plain run–subrun lists that reflect the beam-quality cut (“_wcut”) bookkeeping.
# No env setup, no arguments, no stdout chatter. Outputs only text files.

set -euo pipefail
: "${LC_ALL:=C}"; export LC_ALL

DBROOT="/exp/uboone/data/uboonebeam/beamdb"
SLIP_DIR="/exp/uboone/app/users/guzowski/slip_stacking"

RUN_DB="$DBROOT/run.db"

# Prefer v2, fallback to v1
BNB_DB="$DBROOT/bnb_v2.db";   [[ -f "$BNB_DB"  ]] || BNB_DB="$DBROOT/bnb_v1.db"
NUMI_DB="$DBROOT/numi_v2.db"; [[ -f "$NUMI_DB" ]] || NUMI_DB="$DBROOT/numi_v1.db"
NUMI_V4_DB="$SLIP_DIR/numi_v4.db"

# Run periods used in your previous script
declare -A RUN_START RUN_END
RUN_START[1]="2015-10-01T00:00:00"; RUN_END[1]="2016-07-31T23:59:59"
RUN_START[2]="2016-10-01T00:00:00"; RUN_END[2]="2017-07-31T23:59:59"
RUN_START[3]="2017-10-01T00:00:00"; RUN_END[3]="2018-07-31T23:59:59"
RUN_START[4]="2018-10-01T00:00:00"; RUN_END[4]="2019-07-31T23:59:59"
RUN_START[5]="2019-10-01T00:00:00"; RUN_END[5]="2020-03-31T23:59:59"

# Output files (overwrite cleanly)
out_bnb="bnb_good_runs.txt"
out_numi_all="numi_good_runs.txt"
out_numi_fhc="numi_FHC_good_runs.txt"
out_numi_rhc="numi_RHC_good_runs.txt"
: > "$out_bnb"
: > "$out_numi_all"
: > "$out_numi_fhc"
: > "$out_numi_rhc"

# Helper: unique+sorted append
append_unique_sorted() {
  local infile="$1" outfile="$2"
  if [[ -s "$infile" ]]; then
    # merge, sort unique
    { cat "$outfile" "$infile" 2>/dev/null || true; } \
      | sort -n | uniq > "${outfile}.tmp"
    mv "${outfile}.tmp" "$outfile"
  fi
}

for r in 1 2 3 4 5; do
  S="${RUN_START[$r]}"; E="${RUN_END[$r]}"

  # BNB good subruns: join to BNB table and require any quality-cut metric > 0
  tmp="tmp_bnb_run${r}.txt"
  sqlite3 -noheader -batch "$RUN_DB" >"$tmp" <<SQL
ATTACH '$BNB_DB' AS bnb;
.mode list
.separator " "
SELECT r.run, r.subrun
FROM runinfo AS r
JOIN bnb.bnb AS b
  ON r.run=b.run AND r.subrun=b.subrun
WHERE r.begin_time >= '$S' AND r.end_time <= '$E'
  AND (IFNULL(b.E1DCNT,0)>0 OR IFNULL(b.tor860,0)>0 OR IFNULL(b.tor875,0)>0)
ORDER BY r.run, r.subrun;
SQL
  append_unique_sorted "$tmp" "$out_bnb"
  rm -f "$tmp"

  # NuMI good subruns (TOTAL): require any *_wcut accounting present (>0)
  tmp="tmp_numi_all_run${r}.txt"
  sqlite3 -noheader -batch "$RUN_DB" >"$tmp" <<SQL
ATTACH '$NUMI_DB' AS numi;
.mode list
.separator " "
SELECT r.run, r.subrun
FROM runinfo AS r
JOIN numi.numi AS n
  ON r.run=n.run AND r.subrun=n.subrun
WHERE r.begin_time >= '$S' AND r.end_time <= '$E'
  AND (IFNULL(n.EA9CNT,0)>0 OR IFNULL(n.tor101,0)>0 OR IFNULL(n.tortgt,0)>0)
ORDER BY r.run, r.subrun;
SQL
  append_unique_sorted "$tmp" "$out_numi_all"
  rm -f "$tmp"

  # If horn-split DB exists, also write FHC and RHC lists
  if [[ -f "$NUMI_V4_DB" ]]; then
    tmp="tmp_numi_fhc_run${r}.txt"
    sqlite3 -noheader -batch "$RUN_DB" >"$tmp" <<SQL
ATTACH '$NUMI_V4_DB' AS n4;
.mode list
.separator " "
SELECT r.run, r.subrun
FROM runinfo AS r
JOIN n4.numi AS n
  ON r.run=n.run AND r.subrun=n.subrun
WHERE r.begin_time >= '$S' AND r.end_time <= '$E'
  AND IFNULL(n.EA9CNT_fhc,0)>0
ORDER BY r.run, r.subrun;
SQL
    append_unique_sorted "$tmp" "$out_numi_fhc"
    rm -f "$tmp"

    tmp="tmp_numi_rhc_run${r}.txt"
    sqlite3 -noheader -batch "$RUN_DB" >"$tmp" <<SQL
ATTACH '$NUMI_V4_DB' AS n4;
.mode list
.separator " "
SELECT r.run, r.subrun
FROM runinfo AS r
JOIN n4.numi AS n
  ON r.run=n.run AND r.subrun=n.subrun
WHERE r.begin_time >= '$S' AND r.end_time <= '$E'
  AND IFNULL(n.EA9CNT_rhc,0)>0
ORDER BY r.run, r.subrun;
SQL
    append_unique_sorted "$tmp" "$out_numi_rhc"
    rm -f "$tmp"
  fi
done

# Remove empty horn-split files if numi_v4.db wasn't present
[[ -f "$NUMI_V4_DB" ]] || { rm -f "$out_numi_fhc" "$out_numi_rhc"; : > /dev/null; }

# The script prints nothing, it only leaves these files:
#   - bnb_good_runs.txt
#   - numi_good_runs.txt
#   - numi_FHC_good_runs.txt (if numi_v4.db exists)
#   - numi_RHC_good_runs.txt (if numi_v4.db exists)
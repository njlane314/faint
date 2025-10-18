#!/bin/bash
: "${LC_ALL:=C}"; export LC_ALL

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup sam_web_client

command -v sqlite3 >/dev/null || { echo "sqlite3 not found"; exit 1; }

DBROOT="/exp/uboone/data/uboonebeam/beamdb"
SLIP_DIR="/exp/uboone/app/users/guzowski/slip_stacking"

BNB_DB="$DBROOT/bnb_v2.db";   [[ -f "$BNB_DB"  ]] || BNB_DB="$DBROOT/bnb_v1.db"
NUMI_DB="$DBROOT/numi_v2.db"; [[ -f "$NUMI_DB" ]] || NUMI_DB="$DBROOT/numi_v1.db"
RUN_DB="$DBROOT/run.db"
NUMI_V4_DB="$SLIP_DIR/numi_v4.db"

have_numi_v4=false; [[ -f "$NUMI_V4_DB" ]] && have_numi_v4=true

declare -A RUN_START RUN_END
RUN_START[1]="2015-10-01T00:00:00"; RUN_END[1]="2016-07-31T23:59:59"
RUN_START[2]="2016-10-01T00:00:00"; RUN_END[2]="2017-07-31T23:59:59"
RUN_START[3]="2017-10-01T00:00:00"; RUN_END[3]="2018-07-31T23:59:59"
RUN_START[4]="2018-10-01T00:00:00"; RUN_END[4]="2019-07-31T23:59:59"
RUN_START[5]="2019-10-01T00:00:00"; RUN_END[5]="2020-03-31T23:59:59"

bar() { printf '%*s\n' 100 '' | tr ' ' '='; }
hdr() { printf "%-7s | %-6s | %14s %14s %14s %14s %14s %14s %14s %14s\n" \
              "Run" "$1" "EXT" "Gate" "Cnt" "TorA" "TorB/Target" "Cnt_wcut" "TorA_wcut" "TorB/Target_wcut"; }
row() { printf "%-7s | %-6s | %14.1f %14.1f %14.1f %14.4g %14.4g %14.1f %14.4g %14.4g\n" "$@"; }

bnb_totals_sql() {
  local S="$1" E="$2"
  sqlite3 -noheader -list "$RUN_DB" <<SQL
ATTACH '$BNB_DB' AS bnb;

WITH rset AS (
  SELECT run, subrun
  FROM runinfo
  WHERE begin_time >= '$S' AND end_time <= '$E'
),
rinfo_sums AS (
  SELECT
    IFNULL(SUM(EXTTrig), 0.0) AS ext_sum,
    IFNULL(SUM(Gate2Trig), 0.0) AS gate2_sum,
    IFNULL(SUM(E1DCNT), 0.0) AS e1dcnt_sum,
    IFNULL(SUM(tor860)*1e12, 0.0) AS tor860_sum,
    IFNULL(SUM(tor875)*1e12, 0.0) AS tor875_sum
  FROM runinfo
  WHERE begin_time >= '$S' AND end_time <= '$E'
),
bnb_dedup AS (
  SELECT run, subrun,
         MAX(E1DCNT) AS E1DCNT,
         MAX(tor860) AS tor860,
         MAX(tor875) AS tor875
  FROM bnb.bnb
  GROUP BY run, subrun
),
bnb_sums AS (
  SELECT
    IFNULL(SUM(d.E1DCNT), 0.0) AS e1dcnt_w,
    IFNULL(SUM(d.tor860)*1e12, 0.0) AS tor860_w,
    IFNULL(SUM(d.tor875)*1e12, 0.0) AS tor875_w
  FROM bnb_dedup d
  JOIN rset r USING(run, subrun)
)
SELECT
  rinfo_sums.ext_sum,
  rinfo_sums.gate2_sum,
  rinfo_sums.e1dcnt_sum,
  rinfo_sums.tor860_sum,
  rinfo_sums.tor875_sum,
  bnb_sums.e1dcnt_w,
  bnb_sums.tor860_w,
  bnb_sums.tor875_w
FROM rinfo_sums, bnb_sums;
SQL
}

numi_totals_sql() {
  local S="$1" E="$2"
  sqlite3 -noheader -list "$RUN_DB" <<SQL
ATTACH '$NUMI_DB' AS numi;

WITH rset AS (
  SELECT run, subrun
  FROM runinfo
  WHERE begin_time >= '$S' AND end_time <= '$E'
),
rinfo_sums AS (
  SELECT
    IFNULL(SUM(EXTTrig), 0.0) AS ext_sum,
    IFNULL(SUM(Gate1Trig), 0.0) AS gate1_sum,
    IFNULL(SUM(EA9CNT), 0.0) AS ea9_sum,
    IFNULL(SUM(tor101)*1e12, 0.0) AS tor101_sum,
    IFNULL(SUM(tortgt)*1e12, 0.0) AS tortgt_sum
  FROM runinfo
  WHERE begin_time >= '$S' AND end_time <= '$E'
),
numi_dedup AS (
  SELECT run, subrun,
         MAX(EA9CNT) AS EA9CNT,
         MAX(tor101) AS tor101,
         MAX(tortgt) AS tortgt
  FROM numi.numi
  GROUP BY run, subrun
),
numi_sums AS (
  SELECT
    IFNULL(SUM(d.EA9CNT), 0.0) AS ea9_w,
    IFNULL(SUM(d.tor101)*1e12, 0.0) AS tor101_w,
    IFNULL(SUM(d.tortgt)*1e12, 0.0) AS tortgt_w
  FROM numi_dedup d
  JOIN rset r USING(run, subrun)
)
SELECT
  rinfo_sums.ext_sum,
  rinfo_sums.gate1_sum,
  rinfo_sums.ea9_sum,
  rinfo_sums.tor101_sum,
  rinfo_sums.tortgt_sum,
  numi_sums.ea9_w,
  numi_sums.tor101_w,
  numi_sums.tortgt_w
FROM rinfo_sums, numi_sums;
SQL
}

numi_horns_sql() {
  local S="$1" E="$2" which="$3"
  local col_sfx; [[ "$which" == "FHC" ]] && col_sfx="fhc" || col_sfx="rhc"
  sqlite3 -noheader -list "$RUN_DB" <<SQL
ATTACH '$NUMI_V4_DB' AS n4;

WITH rset AS (
  SELECT run, subrun
  FROM runinfo
  WHERE begin_time >= '$S' AND end_time <= '$E'
),
n4_dedup AS (
  SELECT run, subrun,
         MAX(EA9CNT_${col_sfx}) AS EA9CNT_${col_sfx},
         MAX(tor101_${col_sfx}) AS tor101_${col_sfx},
         MAX(tortgt_${col_sfx}) AS tortgt_${col_sfx}
  FROM n4.numi
  GROUP BY run, subrun
)
SELECT
  IFNULL(SUM(d.EA9CNT_${col_sfx}), 0.0),
  IFNULL(SUM(d.tor101_${col_sfx})*1e12, 0.0),
  IFNULL(SUM(d.tortgt_${col_sfx})*1e12, 0.0)
FROM n4_dedup d
JOIN rset r USING(run, subrun);
SQL
}

numi_ext_split_sql() {
  local S="$1" E="$2" which="$3"
  local col_sfx; [[ "$which" == "FHC" ]] && col_sfx="fhc" || col_sfx="rhc"
  sqlite3 -noheader -list "$RUN_DB" <<SQL
ATTACH '$NUMI_V4_DB' AS n4;

WITH rset AS (
  SELECT run, subrun, EXTTrig
  FROM runinfo
  WHERE begin_time >= '$S' AND end_time <= '$E'
),
n4_tag AS (
  SELECT run, subrun
  FROM n4.numi
  WHERE IFNULL(EA9CNT_${col_sfx},0) > 0
  GROUP BY run, subrun
)
SELECT IFNULL(SUM(r.EXTTrig), 0.0)
FROM rset r
JOIN n4_tag t USING(run, subrun);
SQL
}

numi_runinfo_split_sql() {
  local S="$1" E="$2" which="$3"
  local col_sfx; [[ "$which" == "FHC" ]] && col_sfx="fhc" || col_sfx="rhc"
  sqlite3 -noheader -list "$RUN_DB" <<SQL
ATTACH '$NUMI_V4_DB' AS n4;

WITH rset AS (
  SELECT run, subrun, Gate1Trig, EA9CNT, tor101, tortgt
  FROM runinfo
  WHERE begin_time >= '$S' AND end_time <= '$E'
),
n4_tag AS (
  SELECT run, subrun
  FROM n4.numi
  WHERE IFNULL(EA9CNT_${col_sfx},0) > 0
  GROUP BY run, subrun
)
SELECT
  IFNULL(SUM(r.Gate1Trig), 0.0),
  IFNULL(SUM(r.EA9CNT), 0.0),
  IFNULL(SUM(r.tor101)*1e12, 0.0),
  IFNULL(SUM(r.tortgt)*1e12, 0.0)
FROM rset r
JOIN n4_tag t USING(run, subrun);
SQL
}

bar
echo "BNB & NuMI POT / Toroid by run period"
echo "FAST SQL path (sqlite3, duplicate-safe); DBROOT: $DBROOT"
$have_numi_v4 && echo "NuMI FHC/RHC splits from $NUMI_V4_DB" || echo "NuMI FHC/RHC splits disabled (numi_v4.db not found)"
bar

for r in 1 2 3 4 5; do
  S="${RUN_START[$r]}"; E="${RUN_END[$r]}"

  echo "Run $r  (${S} → ${E})"
  echo

  hdr "BNB"
  IFS='|' read -r EXT Gate2 Cnt TorA TorB CntW TorAW TorBW < <(bnb_totals_sql "$S" "$E")
  row "$r" "TOTAL" "$EXT" "$Gate2" "$Cnt" "$TorA" "$TorB" "$CntW" "$TorAW" "$TorBW"

  echo
  hdr "NuMI"
  IFS='|' read -r EXT1 Gate1 CntN Tor101 Tortgt CntNW Tor101W TortgtW < <(numi_totals_sql "$S" "$E")
  row "$r" "TOTAL" "$EXT1" "$Gate1" "$CntN" "$Tor101" "$Tortgt" "$CntNW" "$Tor101W" "$TortgtW"

  if $have_numi_v4; then
    IFS='|' read -r cntF_w tor101F_w tortgtF_w < <(numi_horns_sql "$S" "$E" "FHC")
    IFS='|' read -r cntR_w tor101R_w tortgtR_w < <(numi_horns_sql "$S" "$E" "RHC")

    IFS='|' read -r gateF cntF tor101F tortgtF < <(numi_runinfo_split_sql "$S" "$E" "FHC")
    IFS='|' read -r gateR cntR tor101R tortgtR < <(numi_runinfo_split_sql "$S" "$E" "RHC")

    extF=$(numi_ext_split_sql "$S" "$E" "FHC")
    extR=$(numi_ext_split_sql "$S" "$E" "RHC")

    row "$r" "FHC" "$extF" "$gateF" "$cntF" "$tor101F" "$tortgtF" "$cntF_w" "$tor101F_w" "$tortgtF_w"
    row "$r" "RHC" "$extR" "$gateR" "$cntR" "$tor101R" "$tortgtR" "$cntR_w" "$tor101R_w" "$tortgtR_w"
  fi

  echo
  bar
done

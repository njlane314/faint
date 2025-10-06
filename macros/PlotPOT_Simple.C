#include "TCanvas.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TSystem.h"
#include <sqlite3.h>
#include <ctime>
#include <cstring>
#include <string>
#include <iostream>
#include <algorithm>

static std::string DBROOT() {
  if (const char* e = gSystem->Getenv("DBROOT")) return e;
  return "/exp/uboone/data/uboonebeam/beamdb";
}
static std::string SLIP() {
  if (const char* e = gSystem->Getenv("SLIP_DIR")) return e;
  return "/exp/uboone/app/users/guzowski/slip_stacking";
}

static bool parse_time(const char* s, std::tm& out)
{
  if (!s || !*s) return false;
  const char* formats[] = {
    "%Y-%m-%dT%H:%M:%S",
    "%Y-%m-%d %H:%M:%S"
  };
  for (const char* fmt : formats) {
    std::tm tm{};
    if (strptime(s, fmt, &tm)) { out = tm; return true; }
  }
  return false;
}

static bool iso2utc(const char* s, time_t& out) {
  std::tm tm{};
  if (!parse_time(s, tm)) return false;
  out = timegm(&tm);
  return true;
}
static time_t sunday_after_or_on(time_t t) {
  std::tm gm = *gmtime(&t);
  int add = (7 - gm.tm_wday) % 7;
  gm.tm_hour = gm.tm_min = gm.tm_sec = 0;
  time_t day0 = timegm(&gm);
  return day0 + add*86400;
}

void PlotPOT_Simple(const char* outstem = "pot_timeline")
{
  if (gSystem->Load("libsqlite3.so") < 0) {
    std::cerr << "Failed to load libsqlite3.so\n";
    return;
  }
  const std::string run_db  = DBROOT() + "/run.db";
  const std::string bnb_db  = gSystem->AccessPathName((DBROOT()+"/bnb_v2.db").c_str()) ?
                              DBROOT()+"/bnb_v1.db" : DBROOT()+"/bnb_v2.db";
  const std::string numi_db = gSystem->AccessPathName((DBROOT()+"/numi_v2.db").c_str()) ?
                              DBROOT()+"/numi_v1.db" : DBROOT()+"/numi_v2.db";
  const std::string n4_db   = SLIP()   + "/numi_v4.db";

  sqlite3* db = nullptr;
  if (sqlite3_open(run_db.c_str(), &db) != SQLITE_OK) {
    std::cerr << "Could not open " << run_db << "\n"; return;
  }
  auto exec = [&](const std::string& sql){
    char* err=nullptr; int rc = sqlite3_exec(db, sql.c_str(), nullptr, nullptr, &err);
    if (rc!=SQLITE_OK && err) { std::cerr << "SQL error: " << err << "\n"; sqlite3_free(err); }
  };
  exec("ATTACH DATABASE '"+bnb_db +"' AS bnb;");
  exec("ATTACH DATABASE '"+numi_db+"' AS numi;");
  exec("ATTACH DATABASE '"+n4_db  +"' AS n4;");

  auto fetch_bound = [&](const char* sql){
    sqlite3_stmt* st=nullptr;
    time_t result = 0;
    if (sqlite3_prepare_v2(db, sql, -1, &st, nullptr)==SQLITE_OK) {
      if (sqlite3_step(st)==SQLITE_ROW) {
        const char* s = (const char*)sqlite3_column_text(st,0);
        time_t tmp;
        if (s && iso2utc(s, tmp) && tmp>0) result = tmp;
      }
    }
    sqlite3_finalize(st);
    return result;
  };

  const char* sql_min =
    "SELECT begin_time FROM runinfo "
    "WHERE begin_time IS NOT NULL "
    "  AND begin_time NOT IN ('1970-01-01T00:00:00','1970-01-01 00:00:00') "
    "ORDER BY begin_time ASC LIMIT 1;";
  const char* sql_max =
    "SELECT begin_time FROM runinfo "
    "WHERE begin_time IS NOT NULL "
    "  AND begin_time NOT IN ('1970-01-01T00:00:00','1970-01-01 00:00:00') "
    "ORDER BY begin_time DESC LIMIT 1;";
  time_t tmin = fetch_bound(sql_min);
  time_t tmax = fetch_bound(sql_max);
  if (!tmin || !tmax || tmin>=tmax) { std::cerr<<"No time range\n"; sqlite3_close(db); return; }

  const double W = 7.0*86400.0;
  double xlo = (double)(sunday_after_or_on(tmin) - 7*86400);
  double xhi = (double)(sunday_after_or_on(tmax) + 7*86400);
  int nbins = std::max(1, int((xhi - xlo)/W + 0.5));

  TH1D hBNB("hBNB","",nbins,xlo,xhi), hFHC("hFHC","",nbins,xlo,xhi), hRHC("hRHC","",nbins,xlo,xhi);
  hBNB.SetDirectory(nullptr); hFHC.SetDirectory(nullptr); hRHC.SetDirectory(nullptr);
  hBNB.SetFillColorAlpha(kGreen+2,0.95); hBNB.SetLineColor(kGreen+2);
  hFHC.SetFillColorAlpha(kOrange+7,0.95); hFHC.SetLineColor(kOrange+7);
  hRHC.SetFillColorAlpha(kRed+1,0.95);   hRHC.SetLineColor(kRed+1);

  auto fill_from = [&](const char* sql, TH1D& h){
    sqlite3_stmt* s=nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &s, nullptr)!=SQLITE_OK) return;
    while (sqlite3_step(s)==SQLITE_ROW) {
      const char* bt = (const char*)sqlite3_column_text(s,0);
      double pot = sqlite3_column_double(s,1);
      if (!bt || pot<=0) continue;
      time_t tt;
      if (!iso2utc(bt, tt)) continue;
      double t = (double)tt;
      if (t<=0) continue;
      h.Fill(t, pot/1e18);
    }
    sqlite3_finalize(s);
  };

  const char* q_bnb =
    "SELECT r.begin_time, 1e12 * ("
    " CASE WHEN IFNULL(b.tor875,0)>0 THEN IFNULL(b.tor875,r.tor875) "
    "      WHEN IFNULL(r.tor875,0)>0 THEN r.tor875 "
    "      WHEN IFNULL(b.tor860,0)>0 THEN IFNULL(b.tor860,r.tor860) "
    "      ELSE IFNULL(r.tor860,0) END ) AS pot "
    "FROM runinfo r LEFT JOIN bnb.bnb b ON r.run=b.run AND r.subrun=b.subrun "
    "WHERE (IFNULL(b.tor875,0)+IFNULL(r.tor875,0)+IFNULL(b.tor860,0)+IFNULL(r.tor860,0))>0;";

  const char* q_fhc =
    "SELECT r.begin_time, 1e12 * (CASE WHEN IFNULL(n.tortgt_fhc,0)>0 THEN n.tortgt_fhc "
    "                                  ELSE IFNULL(n.tor101_fhc,0) END) AS pot "
    "FROM runinfo r JOIN n4.numi n ON r.run=n.run AND r.subrun=n.subrun "
    "WHERE IFNULL(n.EA9CNT_fhc,0)>0;";
  const char* q_rhc =
    "SELECT r.begin_time, 1e12 * (CASE WHEN IFNULL(n.tortgt_rhc,0)>0 THEN n.tortgt_rhc "
    "                                  ELSE IFNULL(n.tor101_rhc,0) END) AS pot "
    "FROM runinfo r JOIN n4.numi n ON r.run=n.run AND r.subrun=n.subrun "
    "WHERE IFNULL(n.EA9CNT_rhc,0)>0;";

  fill_from(q_bnb, hBNB);
  fill_from(q_fhc, hFHC);
  fill_from(q_rhc, hRHC);

  sqlite3_close(db);

  std::vector<double> x(nbins), cum(nbins), scaled(nbins);
  double maxStack = 0, sum=0, maxCum=0;
  for (int i=1;i<=nbins;++i) {
    double s = hBNB.GetBinContent(i)+hFHC.GetBinContent(i)+hRHC.GetBinContent(i);
    maxStack = std::max(maxStack, s);
    sum += s*1e18;
    double c = sum/1e20;
    cum[i-1]=c; maxCum=std::max(maxCum,c);
    x[i-1]=hBNB.GetXaxis()->GetBinCenter(i);
  }
  double yMax = (maxStack>0)? maxStack*1.15 : 1.0;
  double scale = (maxCum>0)? yMax/maxCum : 1.0;
  for (int i=0;i<nbins;++i) scaled[i]=cum[i]*scale;

  TCanvas c("c","POT timeline",1600,500);
  c.SetMargin(0.08,0.10,0.12,0.06);

  THStack hs("hs","");
  hs.Add(&hBNB); hs.Add(&hFHC); hs.Add(&hRHC);
  hs.Draw("hist");
  hs.GetXaxis()->SetTimeDisplay(1);
  hs.GetXaxis()->SetTimeFormat("%Y");
  hs.GetYaxis()->SetTitle("Protons per week  (#times 10^{18})");
  hs.SetMaximum(yMax);
  hs.SetMinimum(0);

  TGraph g(nbins, x.data(), scaled.data());
  g.SetLineColor(kBlue+1); g.SetLineWidth(3);
  g.Draw("L SAME");

  TGaxis right(hs.GetXaxis()->GetXmax(), 0, hs.GetXaxis()->GetXmax(), yMax,
               0, maxCum, 510, "+L");
  right.SetLabelColor(kBlue+1);
  right.SetTitleColor(kBlue+1);
  right.SetTitle("Total Protons  (#times 10^{20})");
  right.Draw();

  TLegend leg(0.12,0.75,0.40,0.90);
  leg.SetBorderSize(0); leg.SetFillStyle(0);
  leg.AddEntry(&hBNB,"BNB (\\nu)","f");
  leg.AddEntry(&hFHC,"NuMI FHC (\\nu)","f");
  leg.AddEntry(&hRHC,"NuMI RHC (\\bar{\\nu})","f");
  leg.AddEntry(&g,   "Total POT (cumulative)","l");
  leg.Draw();

  c.SaveAs(Form("%s.png", outstem));
  c.SaveAs(Form("%s.pdf", outstem));
}


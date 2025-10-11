#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TH2D.h>
#include <TMatrixDSym.h>
#include <TMatrixD.h>
#include <TLatex.h>

#include "rarexsec/Hub.hh"
#include "rarexsec/syst/Systematics.hh"

using namespace rarexsec;

// simple correlation TH2 builder
static TH2D* makeCorrTH2(const TMatrixDSym& C, const char* name) {
    const int n = C.GetNrows();
    TH2D* h = new TH2D(name, name, n, 0.5, n+0.5, n, 0.5, n+0.5);
    for (int i=0;i<n;++i) {
        const double sii = C(i,i) > 0 ? std::sqrt(C(i,i)) : 0.0;
        for (int j=0;j<n;++j) {
            const double sjj = C(j,j) > 0 ? std::sqrt(C(j,j)) : 0.0;
            const double rho = (sii>0 && sjj>0) ? (C(i,j)/(sii*sjj)) : 0.0;
            h->SetBinContent(i+1, j+1, rho);
        }
    }
    return h;
}

static std::unique_ptr<TGraphAsymmErrors> bandFromCov(const TH1D& H, const TMatrixDSym& C) {
    const int n = H.GetNbinsX();
    auto g = std::make_unique<TGraphAsymmErrors>(n);
    for (int i=1;i<=n;++i) {
        const double x1 = H.GetXaxis()->GetBinLowEdge(i);
        const double x2 = H.GetXaxis()->GetBinUpEdge(i);
        const double xc = 0.5*(x1+x2);
        const double halfw = 0.5*(x2-x1);
        const double y  = H.GetBinContent(i);
        const double sy = C(i-1,i-1) > 0 ? std::sqrt(C(i-1,i-1)) : 0.0;
        g->SetPoint(i-1, xc, y);
        g->SetPointError(i-1, halfw, halfw, sy, sy);
    }
    return g;
}

void make_phase_space_combined_plots() {
    ROOT::EnableImplicitMT();

    if (gSystem->Load("librarexsec") < 0) {
        throw std::runtime_error("Failed to load librexsec");
    }

    const std::string config_path = "data/samples.json";
    const std::string beamline    = "numi-fhc";
    const std::vector<std::string> periods = {"run1"};

    Hub hub(config_path);
    const auto sim  = hub.simulation_entries(beamline, periods);
    const auto data = hub.data_entries(beamline, periods);

    std::vector<const Entry*> A_beam, B_strange;
    for (auto* e : sim) {
        if (!e) continue;
        if (e->kind == sample::origin::beam)         A_beam.push_back(e);
        else if (e->kind == sample::origin::strangeness) B_strange.push_back(e);
    }

    plot::Histogram1DSpec spec;                 // fill to your liking
    spec.id     = "enu";
    spec.title  = "Reconstructed E_{#nu};E_{#nu}^{reco} [GeV];Events";
    spec.nbins  = 20;
    spec.xmin   = 0.0;
    spec.xmax   = 2.0;
    spec.expr   = "reco_nu_energy";    // <-- set your variable
    spec.weight = "evt_weight";        // <-- your nominal weight
    spec.sel    = "";                  // optional selection string

    auto H_A = syst::make_total_mc_hist(spec, A_beam, "_A");
    auto H_B = syst::make_total_mc_hist(spec, B_strange, "_B");

    auto H_data = syst::make_total_mc_hist(spec, data, "_data"); // same helper works for data

    // --- Build block covariances (A⊕B) by source ---
    // MC statistics
    TMatrixDSym C_stat_block = syst::block_diag_stat(*H_A, *H_B);

    // PPFX (NuMI): 600 universes, central ppfx_cv
    const int Nppfx = 600;
    TMatrixDSym C_ppfx_block = syst::block_cov_from_weight_vector_ushort_scaled(
        spec, A_beam, spec, B_strange, "weightsPPFX", Nppfx, "ppfx_cv", 1.0/1000.0);

    // GENIE: All_UBGenie stored in map<string, vector<double>> under key "weights"
    const int Ngenie = 500;
    TMatrixDSym C_genie_block = syst::block_cov_from_map_weight_vector(
        spec, A_beam, spec, B_strange, "weights", "All_UBGenie", Ngenie,
        "weightSplineTimesTune");

    // Reinteractions (Geant4): from map key "reint_all" (if present)
    // Set Nreint to the vector length you produce (often 100)
    const int Nreint = 100;
    TMatrixDSym C_reint_block = syst::block_cov_from_map_weight_vector(
        spec, A_beam, spec, B_strange, "weights", "reint_all", Nreint, "");

    // POT normalization (e.g., 2% fully correlated across A and B)
    TMatrixDSym C_pot_block = syst::pot_cov_block(*H_A, *H_B, 0.02);

    // Sum block covariance
    std::vector<const TMatrixDSym*> src_blocks = {
        &C_stat_block, &C_ppfx_block, &C_genie_block, &C_reint_block, &C_pot_block
    };
    TMatrixDSym C_block_total = syst::sum(src_blocks);

    // --- Propagate to "summed categories" spectrum (per-bin A+B) ---
    auto H_sum = syst::sum_same_binning(*H_A, *H_B, "h_sum");
    TMatrixDSym C_sum_total = syst::sum_covariance_block_same_binning(C_block_total,
                                                                      H_A->GetNbinsX(),
                                                                      H_B->GetNbinsX());

    // Also keep per-source sums for breakdown
    TMatrixDSym C_sum_stat   = syst::sum_covariance_block_same_binning(C_stat_block,   H_A->GetNbinsX(), H_B->GetNbinsX());
    TMatrixDSym C_sum_ppfx   = syst::sum_covariance_block_same_binning(C_ppfx_block,   H_A->GetNbinsX(), H_B->GetNbinsX());
    TMatrixDSym C_sum_genie  = syst::sum_covariance_block_same_binning(C_genie_block,  H_A->GetNbinsX(), H_B->GetNbinsX());
    TMatrixDSym C_sum_reint  = syst::sum_covariance_block_same_binning(C_reint_block,  H_A->GetNbinsX(), H_B->GetNbinsX());
    TMatrixDSym C_sum_pot    = syst::sum_covariance_block_same_binning(C_pot_block,    H_A->GetNbinsX(), H_B->GetNbinsX());

    // Optional: shape-only band (remove single global norm mode)
    // auto C_sum_shape = rarexsec::syst::shape_only(C_sum_total, *H_sum);

    // ----------------- Plot 1: Data vs MC with 1σ band -----------------
    TCanvas c1("c1","Prediction with 1σ band", 900, 700);
    gStyle->SetOptStat(0);

    H_sum->SetLineWidth(2);
    H_sum->SetLineColor(kBlack);
    H_sum->SetMarkerStyle(20);
    H_sum->SetMarkerColor(kBlack);
    H_sum->Draw("hist");

    auto band = bandFromCov(*H_sum, C_sum_total);
    band->SetFillStyle(3004);
    band->SetFillColor(kGray+1);
    band->SetLineColor(kGray+1);
    band->Draw("e2 same");

    if (H_data) {
        H_data->SetMarkerStyle(20);
        H_data->SetMarkerColor(kBlue+2);
        H_data->SetLineColor(kBlue+2);
        H_data->Draw("e1 same");
    }

    TLegend leg1(0.62,0.70,0.88,0.88);
    leg1.AddEntry(H_sum.get(), "MC (A#oplusB)", "l");
    leg1.AddEntry(band.get(),  "Total uncert. (1#sigma)", "f");
    if (H_data) leg1.AddEntry(H_data.get(), "Data", "lep");
    leg1.Draw();

    c1.Update();
    c1.SaveAs("plot_prediction_band.pdf");

    // ------------- Plot 2: fractional per-source uncertainties ----------
    TCanvas c2("c2","Fractional by source", 900, 700);
    auto frac_from = [&](const TMatrixDSym& C) {
        auto H = std::make_unique<TH1D>(*H_sum);
        H->SetName(("frac_"+std::string(C.GetName())).c_str());
        for (int i=1;i<=H->GetNbinsX();++i) {
            const double y = H_sum->GetBinContent(i);
            const double s = C(i-1,i-1)>0 ? std::sqrt(C(i-1,i-1)) : 0.0;
            H->SetBinContent(i, (y>0) ? s/y : 0.0);
            H->SetBinError(i, 0.0);
        }
        return H;
    };
    auto Hf_stat  = frac_from(C_sum_stat);
    auto Hf_ppfx  = frac_from(C_sum_ppfx);
    auto Hf_genie = frac_from(C_sum_genie);
    auto Hf_reint = frac_from(C_sum_reint);
    auto Hf_pot   = frac_from(C_sum_pot);

    Hf_stat->SetLineColor(kGreen+2);
    Hf_ppfx->SetLineColor(kOrange-3);
    Hf_genie->SetLineColor(kRed+1);
    Hf_reint->SetLineColor(kMagenta+2);
    Hf_pot->SetLineColor(kBlue+2);

    Hf_stat->SetMaximum(0.5); // adjust to your scale
    Hf_stat->GetYaxis()->SetTitle("Fractional uncertainty");
    Hf_stat->Draw("hist");
    Hf_ppfx->Draw("hist same");
    Hf_genie->Draw("hist same");
    Hf_reint->Draw("hist same");
    Hf_pot->Draw("hist same");

    TLegend leg2(0.62,0.66,0.88,0.88);
    leg2.AddEntry(Hf_stat.get(),  "MC stat", "l");
    leg2.AddEntry(Hf_ppfx.get(),  "Flux (PPFX)", "l");
    leg2.AddEntry(Hf_genie.get(), "GENIE", "l");
    leg2.AddEntry(Hf_reint.get(), "Reint (Geant4)", "l");
    leg2.AddEntry(Hf_pot.get(),   "POT norm", "l");
    leg2.Draw();

    c2.Update();
    c2.SaveAs("plot_fractional_by_source.pdf");

    // ---------------- Plot 3: correlation matrix heatmap ----------------
    TCanvas c3("c3","Correlation matrix", 800, 700);
    TH2D* Hcorr = makeCorrTH2(C_sum_total, "corr");
    Hcorr->GetZaxis()->SetRangeUser(-1.0, 1.0);
    Hcorr->SetTitle("Correlation (A#oplusB summed bins)");
    Hcorr->Draw("colz");
    c3.Update();
    c3.SaveAs("plot_correlation.pdf");

    // --------------- Optional: integrated impact ranking ----------------
    // (simple ranking by sqrt(sum_ij C_ij)/sum_y)
    double mu = 0.0; for (int i=1;i<=H_sum->GetNbinsX();++i) mu += H_sum->GetBinContent(i);
    auto int_sigma = [&](const TMatrixDSym& C){ double s2=0; const int n=C.GetNrows(); for(int i=0;i<n;++i)for(int j=0;j<n;++j)s2+=C(i,j); return (s2>0)?std::sqrt(s2):0.0; };
    const double s_stat  = int_sigma(C_sum_stat)/std::max(1.0,mu);
    const double s_ppfx  = int_sigma(C_sum_ppfx)/std::max(1.0,mu);
    const double s_genie = int_sigma(C_sum_genie)/std::max(1.0,mu);
    const double s_reint = int_sigma(C_sum_reint)/std::max(1.0,mu);
    const double s_pot   = int_sigma(C_sum_pot)/std::max(1.0,mu);

    TCanvas c4("c4","Integrated impact", 700, 500);
    TH1D rank("rank","Integrated fractional impact;Source;Fraction",5,0.5,5.5);
    rank.GetXaxis()->SetBinLabel(1,"MC stat");
    rank.GetXaxis()->SetBinLabel(2,"PPFX");
    rank.GetXaxis()->SetBinLabel(3,"GENIE");
    rank.GetXaxis()->SetBinLabel(4,"Reint");
    rank.GetXaxis()->SetBinLabel(5,"POT");
    rank.SetBinContent(1, s_stat);
    rank.SetBinContent(2, s_ppfx);
    rank.SetBinContent(3, s_genie);
    rank.SetBinContent(4, s_reint);
    rank.SetBinContent(5, s_pot);
    rank.SetFillColor(kAzure-9);
    rank.Draw("bar2");
    c4.Update();
    c4.SaveAs("plot_impact_ranking.pdf");
}

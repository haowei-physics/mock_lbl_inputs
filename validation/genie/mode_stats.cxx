// =============================================================================
//  mode_stats.cxx  -  interaction-mode composition of the 12 stock GENIE
//                     G18_10b_02_11b HK samples (no reweighting).
//
//    mode_stats <working_dir> <data_dir> <outdir> [nevents]
//        (<working_dir>, <data_dir>: see ds::kRepo, ds::kData in dial_scan.cxx)
//        One pass over each sample. Writes <outdir>/mode_stats.root,
//        <outdir>/mode_stats.pdf and the tables in <outdir>/mode_stats.md.
//
//  dial_scan.cxx is #included (inside namespace ds) so the sample list, the
//  Category() classification, the Enu/Q2 binning and the constants are the
//  same code. The per-event normalisation, NuFast oscillation, FakeRecon drops
//  and selection curves below are line-for-line the nominal pass of
//  ds::Scan(), which is not callable on its own (it also reweights).
//
//  Modes: dial_scan's 12 categories with DIS split at the tune's Wcut
//  (1.809 GeV) into the non-resonant background (nonRES) and DIS. Charm and
//  strange events stay in "other", as in Category().
// =============================================================================

// Everything dial_scan.cxx and NuFast_LBL.cpp include, included here first so
// that their include guards keep these headers out of namespace ds.
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/XSecSplineList.h"

#include "RwFramework/GReWeight.h"
#include "RwFramework/GSyst.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"
#include "RwCalculators/GReWeightAGKY.h"
#include "RwCalculators/GReWeightDISNuclMod.h"
#include "RwCalculators/GReWeightDeltaradAngle.h"
#include "RwCalculators/GReWeightFGM.h"
#include "RwCalculators/GReWeightNonResonanceBkg.h"
#include "RwCalculators/GReWeightNuXSecCCQE.h"
#include "RwCalculators/GReWeightNuXSecCCQEELFF.h"
#include "RwCalculators/GReWeightNuXSecCCQEaxial.h"
#include "RwCalculators/GReWeightNuXSecCCQEvec.h"
#include "RwCalculators/GReWeightNuXSecCCRES.h"
#include "RwCalculators/GReWeightNuXSecCOH.h"
#include "RwCalculators/GReWeightNuXSecDIS.h"
#include "RwCalculators/GReWeightNuXSecNC.h"
#include "RwCalculators/GReWeightNuXSecNCEL.h"
#include "RwCalculators/GReWeightNuXSecNCRES.h"
#include "RwCalculators/GReWeightResonanceDecay.h"
#include "RwCalculators/GReWeightXSecEmpiricalMEC.h"
#include "RwCalculators/GReWeightXSecMEC.h"

namespace ds {
#include "dial_scan.cxx"   // its main() becomes ds::main and is never called
}

using namespace genie;

// ------------------------------------------------------------------ modes
const double kWcut = 1.809;   // CommonParam Wcut of G18_10b_02_11b
const double kWbkg = 2.0;     // GReWeightNonResonanceBkg acts below this W
const int kNMode = 14;
const char *kMode[kNMode] = {"CCQE", "CCMEC", "CCRES", "CCnonRES", "CCDIS", "CCCOH", "CCother",
                             "NCEL", "NCMEC", "NCRES", "NCnonRES", "NCDIS", "NCCOH", "NCother"};
const int kCol[kNMode] = {kBlue + 1, kAzure + 6, kGreen + 2, kOrange + 1, kRed + 1, kMagenta + 1, kGray + 1,
                          kBlue - 9, kCyan - 7, kGreen - 8, kOrange - 9, kRed - 9, kMagenta - 9, kGray};
// ds::Category index (12 categories) -> mode, with DIS split by W
int Mode(int cat, double W)
{
  static const int m[ds::kNCat] = {0, 1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13};
  if ((cat == 3 || cat == 9) && W < kWcut) return m[cat] - 1;
  return m[cat];
}

const int kNWBins = 64;   // true W, 0-4 GeV in 62.5 MeV bins

// Per-sample tallies. Index conventions: tgt 0=H 1=O; nuc 0=p 1=n 2=other
// (MEC cluster, COH); npi 0=1pi 1=2pi 2=other (NonResonanceBkg definition).
struct Tally {
  double n[kNMode] = {}, w[kNMode] = {};                 // unweighted, weighted
  double tgt[kNMode][2] = {}, nuc[kNMode][3] = {};       // weighted
  double tgtN[kNMode][2] = {}, nucN[kNMode][3] = {};     // unweighted
  double nonres[2][2][3] = {}, nonresN[2][2][3] = {};    // [CC/NC][p/n][npi], W < Wcut
  double dis = 0, disW2 = 0, disDial = 0;               // all DIS; W < 2; dial-affected (W < 2)
  double disN = 0, disW2N = 0, disDialN = 0;
  double res = 0, res1g = 0, res1eta = 0, resN = 0, res1gN = 0, res1etaN = 0;
  double sel[ds::kNSel][kNMode] = {};                    // weighted x selection expectation
  long nev = 0, nkeep = 0, nmis = 0;
  // GENIE process content of NCEL / CCother / NCother: key -> {weighted, unweighted}
  std::map<std::string, std::array<double, 2>> proc[3];
};
const int kProcMode[3] = {7, 6, 13};   // NCEL, CCother, NCother

int main(int argc, char **argv)
{
  if (argc < 4) { std::cerr << "usage: mode_stats <working_dir> <data_dir> <outdir> [nevents]" << std::endl; return 1; }
  ds::kRepo = argv[1];
  ds::kData = argv[2];
  if (!ds::LoadSamples()) return 1;
  const std::string outdir = argv[3];
  const Long64_t nmax = argc > 4 ? atoll(argv[4]) : -1;

  ds::InitGenie(false);
  gStyle->SetOptStat(0);

  // selection curves, as in ds::Scan
  const char *resp_cat[9] = {"numuCCQE", "numubarCCQE", "numuCCnonQE", "numubarCCnonQE",
                             "nueCCQE", "nuebarCCQE", "nueCCnonQE", "nuebarCCnonQE", "NC"};
  const char *mix[3] = {"1Rmu", "1Re0de", "1Re1de"};
  TFile fr((ds::kRepo + "/lbl_config/hk/response.root").c_str());
  TH1D *heff[9], *hmix[9][3];
  for (int c = 0; c < 9; ++c) {
    heff[c] = (TH1D*)fr.Get(Form("h_eff_%s", resp_cat[c]));
    for (int s = 0; s < 3; ++s) hmix[c][s] = (TH1D*)fr.Get(Form("h_mix_%s_%s", resp_cat[c], mix[s]));
  }
  auto lookup = [](TH1D *h, double x) {
    return h->GetBinContent(std::min(std::max(h->FindBin(x), 1), h->GetNbinsX())); };

  TFile fout((outdir + "/mode_stats.root").c_str(), "RECREATE");
  std::vector<Tally> T(ds::kSamples.size());
  std::vector<std::array<std::array<TH1D*, kNMode>, 3>> H(ds::kSamples.size());   // [var][mode]

  for (size_t is = 0; is < ds::kSamples.size(); ++is) {
    const ds::Sample S = ds::kSamples[is];
    const std::string name = S.name;
    const bool isRHC = name.substr(name.size() - 3) == "rhc";
    Tally &t = T[is];
    TStopwatch clock;

    TFile ff((ds::kRepo + "/" + S.flux).c_str());
    const double flux_per_pot = ((TH1*)ff.Get(S.hist.c_str()))->Integral("width") / 0.05 / 1.0e21;

    TFile fg(ds::GhepFile(name).c_str());
    TTree *gtree = (TTree*)fg.Get("gtree");
    NtpMCEventRecord *mcrec = nullptr;
    gtree->SetBranchAddress("gmcrec", &mcrec);
    TFile fl((ds::kData + "/" + S.stem + "_flattree.root").c_str());
    TTree *flat = (TTree*)fl.Get("FlatTree_VARS");
    if (flat->GetEntries() != gtree->GetEntries()) {
      std::cerr << "[ERROR] " << name << ": flat/GHEP entry mismatch" << std::endl; return 1; }
    Long64_t nev = gtree->GetEntries();
    if (nmax > 0 && nmax < nev) nev = nmax;

    double fScaleFactor; float Enu_true, ELep, Q2; int PDGnu, PDGLep; Char_t cc; Bool_t flagCCQE;
    flat->SetBranchAddress("fScaleFactor", &fScaleFactor);
    flat->SetBranchAddress("Enu_true", &Enu_true);
    flat->SetBranchAddress("ELep", &ELep);
    flat->SetBranchAddress("Q2", &Q2);
    flat->SetBranchAddress("PDGnu", &PDGnu);
    flat->SetBranchAddress("PDGLep", &PDGLep);
    flat->SetBranchAddress("cc", &cc);
    flat->SetBranchAddress("flagCCQE", &flagCCQE);

    for (int m = 0; m < kNMode; ++m) {
      H[is][0][m] = new TH1D(Form("enu_%s", kMode[m]), Form("%s %s;E_{#nu}^{true} [GeV];events", name.c_str(), kMode[m]),
                             ds::kNBins, ds::kBins);
      H[is][1][m] = new TH1D(Form("q2_%s", kMode[m]), Form("%s %s;Q^{2}_{true} [GeV^{2}];events", name.c_str(), kMode[m]),
                             ds::kNQ2, ds::kQ2Bins);
      H[is][2][m] = new TH1D(Form("w_%s", kMode[m]), Form("%s %s;W_{true} [GeV];events", name.c_str(), kMode[m]),
                             kNWBins, 0, 4);
      for (int v = 0; v < 3; ++v) H[is][v][m]->SetDirectory(nullptr);
    }

    std::vector<double> nfE(1);
    std::vector<std::array<std::array<double, 3>, 3>> nfP(1);
    for (Long64_t i = 0; i < nev; ++i) {
      gtree->GetEntry(i);
      flat->GetEntry(i);
      EventRecord &ev = *(mcrec->event);
      Interaction *in = ev.Summary();
      ++t.nev;
      // alignment check, as in ds::Scan
      const double enu = ev.Probe()->E();
      if (std::fabs(enu - Enu_true) > 1e-4 * std::max(1.f, Enu_true) || ev.Probe()->Pdg() != PDGnu ||
          (bool)cc != in->ProcInfo().IsWeakCC()) ++t.nmis;

      const bool keep = !(Q2 < -999.0) && !(cc && (PDGLep == 0 || ELep < 0.0));
      if (!keep) { mcrec->Clear(); continue; }
      ++t.nkeep;

      // normalisation and oscillation, as in ds::Scan
      const double pot_w = fScaleFactor * flux_per_pot * S.pot * ds::kNucleons;
      const int a = (std::abs(S.flux_pdg) == 12) ? 0 : 1;
      const int b = (std::abs(PDGnu) == 12) ? 0 : 1;
      double osc;
      if (cc) {
        nfE[0] = (S.flux_pdg > 0) ? Enu_true : -Enu_true;
        ds::NuFast::Probability_Matter_LBL(0.307, 0.0218, 0.528, -1.601, 7.53e-5, 2.509e-3 + 7.53e-5,
                                           295.0, nfE, 2.6 * 0.5, 0, nfP);
        osc = nfP[0][a][b];
      } else {
        osc = (a != b) ? 0.0 : 1.0;
      }
      const double base = pot_w * osc;

      // selection, as in ds::Scan
      int rc;
      if (!cc)                   rc = 8;
      else if (abs(PDGnu) == 12) rc = (flagCCQE ? 4 : 6) + (PDGnu > 0 ? 0 : 1);
      else                       rc = (flagCCQE ? 0 : 2) + (PDGnu > 0 ? 0 : 1);
      const double eff = lookup(heff[rc], Enu_true);
      double sel[ds::kNSel] = {0};
      for (int s = 0; s < 3; ++s) {
        if (isRHC && s == 2) continue;
        sel[(isRHC ? 3 : 0) + s] = eff * lookup(hmix[rc][s], Enu_true);
      }

      // true W: the generated W where the generator set it (RES, DIS),
      // otherwise (q + p_hit)^2 with the struck nucleon / cluster; -1 if none (COH)
      double W = -1;
      if (in->Kine().KVSet(kKVSelW)) W = in->Kine().W(true);
      else if (ev.HitNucleon() && ev.FinalStatePrimaryLepton())
        W = (*ev.HitNucleon()->P4() + *ev.Probe()->P4() - *ev.FinalStatePrimaryLepton()->P4()).M();

      const int m = Mode(ds::Category(in), W);
      t.n[m] += 1; t.w[m] += base;
      const int tg = in->InitState().Tgt().A() == 1 ? 0 : 1;
      const int hn = in->InitState().Tgt().HitNucPdg();
      const int nc = hn == kPdgProton ? 0 : hn == kPdgNeutron ? 1 : 2;
      t.tgt[m][tg] += base; t.tgtN[m][tg] += 1;
      for (int k = 0; k < 3; ++k) if (m == kProcMode[k]) {
        const ProcessInfo &pi = in->ProcInfo();
        const XclsTag &xt = in->ExclTag();
        std::string key = std::string(ScatteringType::AsString(pi.ScatteringTypeId())) + " " +
                          InteractionType::AsString(pi.InteractionTypeId()) + "; " +
                          (tg == 0 ? "H" : "O") + ", hit " + (hn == 0 ? "none" : std::to_string(hn));
        if (xt.IsCharmEvent())   key += "; charm " + std::to_string(xt.CharmHadronPdg());
        if (xt.IsStrangeEvent()) key += "; strange " + std::to_string(xt.StrangeHadronPdg());
        t.proc[k][key][0] += base; t.proc[k][key][1] += 1;
      }
      t.nuc[m][nc] += base; t.nucN[m][nc] += 1;
      for (int s = 0; s < ds::kNSel; ++s) t.sel[s][m] += base * sel[s];
      H[is][0][m]->Fill(enu, base);
      H[is][1][m]->Fill(Q2, base);
      H[is][2][m]->Fill(W, base);

      // NonResonanceBkg: hadron multiplicity as GReWeightNonResonanceBkg::CalcWeight
      if (in->ProcInfo().IsDeepInelastic()) {
        int nhad = 0, nnuc = 0, npi = 0;
        for (int k = 0; k < ev.GetEntries(); ++k) {
          GHepParticle *p = ev.Particle(k);
          if (p->FirstMother() < 0 || ev.Particle(p->FirstMother())->Pdg() != kPdgHadronicSyst) continue;
          ++nhad;
          if (pdg::IsNucleon(p->Pdg())) ++nnuc;
          if (pdg::IsPion(p->Pdg()))    ++npi;
        }
        const int ip = (nhad >= 2 && nhad <= 3 && nnuc == 1 && (npi == 1 || npi == 2)) ? npi - 1 : 2;
        const double Wd = in->Kine().W(true);
        t.dis += base; t.disN += 1;
        if (Wd < kWbkg) {
          t.disW2 += base; t.disW2N += 1;
          if (ip < 2 && nc < 2) { t.disDial += base; t.disDialN += 1; }
        }
        if ((m == 3 || m == 10) && nc < 2) {
          t.nonres[m == 3 ? 0 : 1][nc][ip] += base;
          t.nonresN[m == 3 ? 0 : 1][nc][ip] += 1;
        }
      }

      // resonance decays as GReWeightResonanceDecay::RewBR
      if (in->ProcInfo().IsResonant()) {
        bool g1 = false, e1 = false;
        for (int k = 0; k < ev.GetEntries(); ++k) {
          GHepParticle *p = ev.Particle(k);
          if (!utils::res::IsBaryonResonance(p->Pdg())) continue;
          int ng = 0, ne = 0;
          for (int d = p->FirstDaughter(); d >= 0 && d <= p->LastDaughter(); ++d) {
            ng += ev.Particle(d)->Pdg() == kPdgGamma;
            ne += ev.Particle(d)->Pdg() == kPdgEta;
          }
          g1 |= ng == 1; e1 |= ne == 1;
        }
        t.res += base; t.resN += 1;
        if (g1) { t.res1g += base;   t.res1gN += 1; }
        if (e1) { t.res1eta += base; t.res1etaN += 1; }
      }
      mcrec->Clear();
    }
    std::cout << "[INFO] " << name << ": " << t.nev << " events, " << t.nev - t.nkeep
              << " dropped as in FakeRecon, " << t.nmis << " GHEP/flat mismatches, "
              << clock.RealTime() << " s" << std::endl;
    if (t.nmis > 0) { std::cerr << "[ERROR] GHEP and flat tree are not aligned for " << name << std::endl; return 1; }

    // ---- root output: <sample>/<mode>/{enu,q2,w} + <sample>/summary histograms
    TDirectory *ds_ = fout.mkdir(name.c_str());
    for (int m = 0; m < kNMode; ++m) {
      ds_->mkdir(kMode[m])->cd();
      for (int v = 0; v < 3; ++v) H[is][v][m]->Write();
    }
    ds_->cd();
    auto modeHist = [&](const char *nm, const char *ti, const double *x) {
      TH1D h(nm, ti, kNMode, 0, kNMode);
      for (int m = 0; m < kNMode; ++m) { h.GetXaxis()->SetBinLabel(m + 1, kMode[m]); h.SetBinContent(m + 1, x[m]); }
      h.Write();
    };
    modeHist("mode_weighted", "oscillated expected events by mode", t.w);
    modeHist("mode_unweighted", "generated events (after FakeRecon drops) by mode", t.n);
    double x[kNMode];
    const char *tn[2] = {"H", "O"}, *nn[3] = {"p", "n", "other"};
    for (int k = 0; k < 2; ++k) {
      for (int m = 0; m < kNMode; ++m) x[m] = t.tgt[m][k];
      modeHist(Form("mode_weighted_tgt%s", tn[k]), Form("expected events on %s", tn[k]), x);
    }
    for (int k = 0; k < 3; ++k) {
      for (int m = 0; m < kNMode; ++m) x[m] = t.nuc[m][k];
      modeHist(Form("mode_weighted_hit%s", nn[k]), Form("expected events, hit nucleon %s", nn[k]), x);
    }
    for (int s = 0; s < ds::kNSel; ++s)
      modeHist(Form("mode_sel_%s", ds::kSel[s]), Form("expected selected events, %s", ds::kSel[s]), t.sel[s]);
    TH1D hnr("nonres_weighted", "non-resonant bkg (W<1.809), expected events", 12, 0, 12);
    const char *ccn[2] = {"CC", "NC"}, *pin[3] = {"1pi", "2pi", "other"};
    for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int c = 0; c < 3; ++c) {
      hnr.GetXaxis()->SetBinLabel(a * 6 + b * 3 + c + 1, Form("%s_%s_%s", ccn[a], nn[b], pin[c]));
      hnr.SetBinContent(a * 6 + b * 3 + c + 1, t.nonres[a][b][c]);
    }
    hnr.Write();
  }

  // ---- PDF: one page per sample, stacked by mode
  gStyle->SetPaperSize(28, 20);
  TCanvas cv("cv", "", 1400, 1000);
  const std::string pdf = outdir + "/mode_stats.pdf";
  cv.Print((pdf + "[").c_str());
  const char *vt[3] = {"E_{#nu}^{true} [GeV]", "Q^{2}_{true} [GeV^{2}]", "W_{true} [GeV]"};
  for (size_t is = 0; is < ds::kSamples.size(); ++is) {
    cv.Clear(); cv.Divide(2, 2);
    std::vector<THStack*> st;
    TLegend *leg = new TLegend(0.05, 0.05, 0.95, 0.95);
    leg->SetNColumns(2);
    leg->SetHeader(Form("%s (oscillated, %.1f expected events)", ds::kSamples[is].name.c_str(),
                        std::accumulate(T[is].w, T[is].w + kNMode, 0.0)));
    for (int v = 0; v < 3; ++v) {
      cv.cd(v + 1);
      THStack *s = new THStack(Form("st%zu_%d", is, v), Form("%s;%s;expected events / %s", ds::kSamples[is].name.c_str(), vt[v], v == 1 ? "GeV^{2}" : "GeV"));
      for (int m = 0; m < kNMode; ++m) {
        TH1D *h = H[is][v][m];   // drawn per unit x; the ROOT file keeps events / bin
        if (h->Integral() <= 0) continue;
        h->SetFillColor(kCol[m]); h->SetLineColor(kCol[m]);
        TH1D *hd = (TH1D*)h->Clone(); hd->SetDirectory(nullptr); hd->Scale(1.0, "width");
        s->Add(hd);
        if (v == 0) leg->AddEntry(h, Form("%s (%.1f%%)", kMode[m],
                                          100 * T[is].w[m] / std::accumulate(T[is].w, T[is].w + kNMode, 0.0)), "f");
      }
      s->Draw("hist");
      if (v == 0) s->GetXaxis()->SetRangeUser(0, 5);
      st.push_back(s);
    }
    cv.cd(4); leg->Draw();
    cv.Print(pdf.c_str());
    for (auto s : st) delete s;
    delete leg;
  }
  cv.Print((pdf + "]").c_str());
  fout.Close();

  // ---- markdown tables
  std::ofstream md(outdir + "/mode_stats.md");
  auto pc = [](double a, double b) { return b > 0 ? 100.0 * a / b : 0.0; };
  auto tot = [](const double *x) { return std::accumulate(x, x + kNMode, 0.0); };
  auto header = [&](const char *first, const char *last) {
    md << "| " << first << " |";
    for (int m = 0; m < kNMode; ++m) md << " " << kMode[m] << " |";
    md << " " << last << " |\n|---|";
    for (int m = 0; m <= kNMode; ++m) md << "---:|";
    md << "\n";
  };
  auto modeTable = [&](bool weighted) {
    header("sample", weighted ? "total (events)" : "total (generated)");
    for (size_t is = 0; is < T.size(); ++is) {
      const double *x = weighted ? T[is].w : T[is].n;
      md << "| " << ds::kSamples[is].name << " |";
      for (int m = 0; m < kNMode; ++m) md << Form(" %.1f |", pc(x[m], tot(x)));
      md << (weighted ? Form(" %.1f |\n", tot(x)) : Form(" %.0f |\n", tot(x)));
    }
  };

  md << "# Interaction-mode composition of the 12 stock GENIE samples\n\n"
        "Tune G18_10b_02_11b. nonRES = DIS with W < 1.809 GeV (the tune's Wcut), DIS = W >= 1.809 GeV; "
        "charm/strange events are in \"other\" (as dial_scan's Category()). "
        "Weighted = oscillated expected events with the dial_scan.cxx normalisation. "
        "All numbers are after the FakeRecon drops. Fractions in %.\n\n";
  md << "## Mode fractions, weighted (%)\n\n"; modeTable(true);
  md << "\n## Mode fractions, unweighted (%)\n\n"; modeTable(false);

  md << "\n## Non-resonant background (W < 1.809 GeV) by hit nucleon and pion multiplicity\n\n"
        "Weighted, % of the whole sample. 1pi / 2pi as GReWeightNonResonanceBkg counts them "
        "(2-3 primary hadrons, exactly one nucleon, 1 or 2 pions); \"oth\" = everything else (dial gives w = 1).\n\n"
        "| sample | CC p1pi | CC p2pi | CC n1pi | CC n2pi | CC oth | NC p1pi | NC p2pi | NC n1pi | NC n2pi | NC oth |\n"
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n";
  for (size_t is = 0; is < T.size(); ++is) {
    const Tally &t = T[is]; const double all = tot(t.w);
    md << "| " << ds::kSamples[is].name << " |";
    for (int a = 0; a < 2; ++a) {
      for (int b = 0; b < 2; ++b) for (int c = 0; c < 2; ++c) md << Form(" %.2f |", pc(t.nonres[a][b][c], all));
      md << Form(" %.2f |", pc(t.nonres[a][0][2] + t.nonres[a][1][2], all));
    }
    md << "\n";
  }
  md << "\n## DIS below 2.0 GeV (where GReWeightNonResonanceBkg acts) and RES decays\n\n"
        "| sample | DIS (all W) % of sample | DIS with W < 2.0, % of DIS | NonRESBG-dial-affected, % of sample "
        "| RES % of sample | RES with resonance -> 1 gamma, % of RES | RES with resonance -> 1 eta, % of RES |\n"
        "|---|---:|---:|---:|---:|---:|---:|\n";
  for (size_t is = 0; is < T.size(); ++is) {
    const Tally &t = T[is]; const double all = tot(t.w);
    md << Form("| %s | %.2f | %.1f | %.2f | %.1f | %.3f | %.2f |\n", ds::kSamples[is].name.c_str(), pc(t.dis, all),
               pc(t.disW2, t.dis), pc(t.disDial, all), pc(t.res, all), pc(t.res1g, t.res), pc(t.res1eta, t.res));
  }
  md << "\nUnweighted counts (generated events): ";
  for (size_t is = 0; is < T.size(); ++is)
    md << Form("%s: %.0f DIS-W<2 dial-affected, %.0f RES->1gamma, %.0f RES->1eta; ", ds::kSamples[is].name.c_str(),
               T[is].disDialN, T[is].res1gN, T[is].res1etaN);
  md << "\n";

  auto split = [&](const char *title, int which) {
    md << "\n## " << title << "\n\n";
    header("sample", "all");
    for (size_t is = 0; is < T.size(); ++is) {
      const Tally &t = T[is]; double num = 0, den = 0;
      md << "| " << ds::kSamples[is].name << " |";
      for (int m = 0; m < kNMode; ++m) {
        const double a = which == 0 ? t.tgt[m][0] : t.nuc[m][0];
        const double d = which == 0 ? t.tgt[m][0] + t.tgt[m][1] : t.nuc[m][0] + t.nuc[m][1];
        num += a; den += d;
        md << (d > 0 ? Form(" %.1f |", pc(a, d)) : " - |");
      }
      md << Form(" %.1f |\n", pc(num, den));
    }
  };
  split("Fraction on hydrogen, per mode (weighted, %; rest is oxygen)", 0);
  split("Fraction with a proton as hit nucleon, per mode (weighted, %; of events with a single hit p or n)", 1);
  md << "\nEvents without a single hit nucleon (MEC 2N clusters, COH), weighted % of each mode: ";
  {
    double o[kNMode] = {}, a[kNMode] = {};
    for (auto &t : T) for (int m = 0; m < kNMode; ++m) { o[m] += t.nuc[m][2]; a[m] += t.w[m]; }
    for (int m = 0; m < kNMode; ++m) if (o[m] > 0) md << Form("%s %.1f; ", kMode[m], pc(o[m], a[m]));
  }
  md << "\n";

  md << "\n## Analysis-sample composition (expectation values, all 12 samples summed, weighted %)\n\n";
  header("selection", "total (events)");
  for (int s = 0; s < ds::kNSel; ++s) {
    double x[kNMode] = {};
    for (auto &t : T) for (int m = 0; m < kNMode; ++m) x[m] += t.sel[s][m];
    md << "| " << ds::kSel[s] << " |";
    for (int m = 0; m < kNMode; ++m) md << Form(" %.1f |", pc(x[m], tot(x)));
    md << Form(" %.1f |\n", tot(x));
  }
  md << "\nPer (sample, selection) numbers are in mode_stats.root, <sample>/mode_sel_<selection>.\n";

  // ---- GENIE process content of NCEL, CCother, NCother
  md << "\n## GENIE process content of NCEL, CCother and NCother\n\n"
        "Key = ScatteringType InteractionType; hit target, hit nucleon PDG; charm/strange tag with the "
        "tagged hadron PDG. Weighted % of the mode summed over all 12 samples, % of the mode in "
        "numubar_x_nuebar_rhc, and generated events (all samples).\n";
  for (int k = 0; k < 3; ++k) {
    std::map<std::string, std::array<double, 3>> sum;   // weighted all, weighted numubar_x_nuebar_rhc, unweighted
    double wall = 0, wsel = 0;
    for (size_t is = 0; is < T.size(); ++is)
      for (auto &e : T[is].proc[k]) {
        const bool sel = std::string(ds::kSamples[is].name) == "numubar_x_nuebar_rhc";
        sum[e.first][0] += e.second[0]; sum[e.first][1] += sel ? e.second[0] : 0; sum[e.first][2] += e.second[1];
        wall += e.second[0]; wsel += sel ? e.second[0] : 0;
      }
    std::vector<std::pair<std::string, std::array<double, 3>>> v(sum.begin(), sum.end());
    std::sort(v.begin(), v.end(), [](auto &a, auto &b) { return a.second[2] > b.second[2]; });
    md << "\n### " << kMode[kProcMode[k]] << "\n\n| process | % (all samples, weighted) | % (numubar_x_nuebar_rhc) | generated |\n|---|---:|---:|---:|\n";
    for (auto &e : v)
      md << Form("| %s | %.2f | %.2f | %.0f |\n", e.first.c_str(), pc(e.second[0], wall), pc(e.second[1], wsel), e.second[2]);
  }
  md.close();
  std::cout << "[INFO] wrote " << outdir << "/mode_stats.{root,pdf,md}" << std::endl;
  return 0;
}

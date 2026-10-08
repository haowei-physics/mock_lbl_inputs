// =============================================================================
//  dial_scan.cxx  -  one-at-a-time GENIE Reweight dial scan for the stock
//                    G18_10b_02_11b HK samples of mock_lbl_inputs.
//
//  Two modes (<working_dir>, <data_dir>: see kRepo, kData below):
//    dial_scan scan <working_dir> <data_dir> <sample> <outdir> [nevents] [nshape]
//        Reweights one sample for every dial in kDials at -2,-1,+1,+2 sigma
//        (0 sigma is the nominal, weight = 1 by construction) and writes
//        <outdir>/part_<sample>.root. The 5 dials in kShapeDials run on the
//        first nshape events only (default 5000): in their "shape" mode the
//        calculator does two numerical cross-section integrals per event and
//        point (0.3-0.8 s per event for RES, 60 ms for CCQE; the BY shape
//        dials are fast enough to run in full). Those subsample rates are scaled up to the full sample
//        and compared with a nominal built from the same subsample.
//    dial_scan summarize <working_dir> <data_dir> <outdir>
//        Adds the 12 part files, writes <outdir>/genie_dial_sensitivity.root
//        and <outdir>/ratio_plots.pdf, and prints markdown tables to stdout.
//
//  Reweighting follows grwght1p (GENIE-MC/Reweight R-1_04_02,
//  src/Apps/gRwght1Param.cxx): same tune, same spline file, same calculator
//  classes, same calculator mode switches, same MEC Tl/ctl fix. The
//  differences: only the calculators that handle the dial are adopted (the
//  others return 1 at nominal anyway), and the four tweak points are held in
//  four GReWeight objects so the GHEP file is read once per dial, not once per
//  point.
//
//  Rates are normalised exactly as src/FakeRecon.C does
//  (POT x fScaleFactor x flux x nucleons, NuFast oscillation at T2K Asimov A,
//  the same dropped events), and the analysis-sample yields use the
//  response.root efficiency and mixing curves as FakeSelect.C does, but as
//  expectation values instead of a random draw.
// =============================================================================

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
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
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

// NuFast_LBL.cpp carries its own test main(); rename it out of the way.
#define main nufast_unused_main
#include "NuFast_LBL.cpp"
#undef main

using namespace genie;
using namespace genie::rew;

// ----------------------------------------------------------------- inputs
// Set in main() from the command line: the repository ($working_dir: sample
// manifest, Messenger config, flux and response files) and the data directory
// of the GENIE chain ($data_dir; layout in docs/genie/README.md).
std::string kRepo, kData;
const std::string kTune    = "G18_10b_02_11b";
std::string SplinesXml() { return kData + "/genie/xsec/hk_v1_stock.xml"; }
std::string MsgXml()     { return kRepo + "/gen/genie/Messenger_production.xml"; }
std::string GhepFile(const std::string &name) { return kData + "/genie/" + name + "/" + name + ".ghep.root"; }

struct Sample { std::string name, stem, flux, hist; double pot; int flux_pdg; };
// The rows of gen/genie/samples.txt, in its order (filled by LoadSamples()).
std::vector<Sample> kSamples;
bool LoadSamples()
{
  const std::string path = kRepo + "/gen/genie/samples.txt";
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    line = line.substr(0, line.find('#'));
    line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
    if (line.empty()) continue;
    std::vector<std::string> f;
    std::stringstream ss(line);
    for (std::string c; std::getline(ss, c, '|');) f.push_back(c);
    if (f.size() != 10) { std::cerr << "[ERROR] bad row in " << path << ": " << line << std::endl; return false; }
    kSamples.push_back({f[0], f[1], f[3], f[4], std::stod(f[8]), std::stoi(f[9])});
  }
  if (kSamples.empty()) std::cerr << "[ERROR] no samples in " << path << std::endl;
  return !kSamples.empty();
}
const double kNucleons = 1.132e35;

// The 40 applicable + 13 unclear dials (see docs/genie/report.md for the classification).
const std::vector<std::string> kDials = {
  "MaNCEL", "EtaNCEL",
  "MaCCQE", "NormCCQE", "MaCCQEshape", "CoulombCCQE", "CCQEPauliSupViaKF",
  "NormCCRES", "MaCCRESshape", "MvCCRESshape", "MaCCRES", "MvCCRES",
  "NormNCRES", "MaNCRESshape", "MvNCRESshape", "MaNCRES", "MvNCRES",
  "MaCOHpi", "R0COHpi", "NormCCCOH", "NormNCCOH",
  "NonRESBGvpCC1pi", "NonRESBGvpCC2pi", "NonRESBGvpNC1pi", "NonRESBGvpNC2pi",
  "NonRESBGvnCC1pi", "NonRESBGvnCC2pi", "NonRESBGvnNC1pi", "NonRESBGvnNC2pi",
  "NonRESBGvbarpCC1pi", "NonRESBGvbarpCC2pi", "NonRESBGvbarpNC1pi", "NonRESBGvbarpNC2pi",
  "NonRESBGvbarnCC1pi", "NonRESBGvbarnCC2pi", "NonRESBGvbarnNC1pi", "NonRESBGvbarnNC2pi",
  "AhtBY", "BhtBY", "CV1uBY", "CV2uBY",
  "AhtBYshape", "BhtBYshape", "CV1uBYshape", "CV2uBYshape",
  "AGKYxF1pi", "AGKYpT1pi", "RDecBR1gamma", "RDecBR1eta",
  "NormCCMEC", "NormNCMEC", "FracPN_CCMEC", "FracDelta_CCMEC",
};
const std::vector<std::string> kShapeDials = {
  "MaCCQEshape", "MaCCRESshape", "MvCCRESshape", "MaNCRESshape", "MvNCRESshape",
};
const std::vector<double> kPoints = {-2, -1, 1, 2};   // 0 is the nominal
const char *kPtName[] = {"m2", "m1", "p1", "p2"};

// ---------------------------------------------------------- categories
const int kNFlav = 4;
const char *kFlav[kNFlav] = {"numu", "numubar", "nue", "nuebar"};
const int kNCat = 12;
const char *kCat[kNCat] = {"CCQE", "CCMEC", "CCRES", "CCDIS", "CCCOH", "CCother",
                           "NCEL", "NCMEC", "NCRES", "NCDIS", "NCCOH", "NCother"};
const int kNSel = 5;   // RHC has no 1Re1de, as in FakeSelect.C
const char *kSel[kNSel] = {"FHC_1Rmu", "FHC_1Re0de", "FHC_1Re1de", "RHC_1Rmu", "RHC_1Re0de"};

const int kNBins = 23;
const double kBins[kNBins + 1] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                                  1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0,
                                  5.0, 7.0, 10.0};

const int kNQ2 = 14;
const double kQ2Bins[kNQ2 + 1] = {0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,
                                  1.0, 1.25, 1.5, 2.0, 3.0};

int FlavIndex(int pdg) { return pdg == 14 ? 0 : pdg == -14 ? 1 : pdg == 12 ? 2 : 3; }

int Category(const Interaction *in)
{
  const ProcessInfo &pi = in->ProcInfo();
  const bool exotic = in->ExclTag().IsCharmEvent() || in->ExclTag().IsStrangeEvent();
  if (pi.IsWeakCC()) {
    if (exotic)                    return 5;
    if (pi.IsQuasiElastic())       return 0;
    if (pi.IsMEC())                return 1;
    if (pi.IsResonant())           return 2;
    if (pi.IsDeepInelastic())      return 3;
    if (pi.IsCoherentProduction()) return 4;
    return 5;
  }
  if (exotic)                    return 11;
  if (pi.IsQuasiElastic())       return 6;
  if (pi.IsMEC())                return 7;
  if (pi.IsResonant())           return 8;
  if (pi.IsDeepInelastic())      return 9;
  if (pi.IsCoherentProduction()) return 10;
  return 11;
}

// ------------------------------------------------------------ GENIE setup
void InitGenie(bool splines)
{
  utils::app_init::MesgThresholds(MsgXml());
  utils::app_init::RandGen(12345);
  GHepRecord::SetPrintLevel(0);
  RunOpt::Instance()->SetTuneName(kTune);
  RunOpt::Instance()->BuildTune();
  if (splines) XSecSplineList::Instance()->LoadFromXml(SplinesXml());
  // GReWeight::CalcWeight logs every calculator at NOTICE for every event.
  Messenger::Instance()->SetPriorityLevel("ReW", pWARN);
}

// The calculators grwght1p adopts (FZone and INuke excluded: they exit() for hN).
std::vector<std::pair<std::string, std::function<GReWeightI*()>>> Calculators()
{
  return {
    {"xsec_ncel",       [] { return new GReWeightNuXSecNCEL; }},
    {"xsec_ccqe",       [] { return new GReWeightNuXSecCCQE; }},
    {"xsec_ccqe_elff",  [] { return new GReWeightNuXSecCCQEELFF; }},
    {"xsec_ccqe_axial", [] { return new GReWeightNuXSecCCQEaxial; }},
    {"xsec_ccqe_vec",   [] { return new GReWeightNuXSecCCQEvec; }},
    {"xsec_ccres",      [] { return new GReWeightNuXSecCCRES; }},
    {"xsec_ncres",      [] { return new GReWeightNuXSecNCRES; }},
    {"xsec_nonresbkg",  [] { return new GReWeightNonResonanceBkg; }},
    {"xsec_coh",        [] { return new GReWeightNuXSecCOH; }},
    {"xsec_dis",        [] { return new GReWeightNuXSecDIS; }},
    {"nuclear_qe",      [] { return new GReWeightFGM; }},
    {"hadro_res_decay", [] { return new GReWeightResonanceDecay; }},
    {"hadro_agky",      [] { return new GReWeightAGKY; }},
    {"nuclear_dis",     [] { return new GReWeightDISNuclMod; }},
    {"xsec_nc",         [] { return new GReWeightNuXSecNC; }},
    {"xsec_empmec",     [] { return new GReWeightXSecEmpiricalMEC; }},
    {"xsec_mec",        [] { return new GReWeightXSecMEC; }},
    {"delta_rad",       [] { return new GReWeightDeltaradAngle; }},
  };
}

bool NeedsMode(GSyst_t s)
{
  return s == kXSecTwkDial_MaCCQE || s == kXSecTwkDial_MaCCRES || s == kXSecTwkDial_MvCCRES ||
         s == kXSecTwkDial_MaNCRES || s == kXSecTwkDial_MvNCRES ||
         s == kXSecTwkDial_AhtBYshape || s == kXSecTwkDial_BhtBYshape ||
         s == kXSecTwkDial_CV1uBYshape || s == kXSecTwkDial_CV2uBYshape;
}

// The mode switches of gRwght1Param.cxx:268-307.
void ApplyMode(GReWeightI *calc, GSyst_t s)
{
  if (s == kXSecTwkDial_MaCCQE)
    if (auto c = dynamic_cast<GReWeightNuXSecCCQE*>(calc)) c->SetMode(GReWeightNuXSecCCQE::kModeMa);
  if (s == kXSecTwkDial_MaCCRES || s == kXSecTwkDial_MvCCRES)
    if (auto c = dynamic_cast<GReWeightNuXSecCCRES*>(calc)) c->SetMode(GReWeightNuXSecCCRES::kModeMaMv);
  if (s == kXSecTwkDial_MaNCRES || s == kXSecTwkDial_MvNCRES)
    if (auto c = dynamic_cast<GReWeightNuXSecNCRES*>(calc)) c->SetMode(GReWeightNuXSecNCRES::kModeMaMv);
  if (s == kXSecTwkDial_AhtBYshape || s == kXSecTwkDial_BhtBYshape ||
      s == kXSecTwkDial_CV1uBYshape || s == kXSecTwkDial_CV2uBYshape)
    if (auto c = dynamic_cast<GReWeightNuXSecDIS*>(calc)) c->SetMode(GReWeightNuXSecDIS::kModeABCV12uShape);
}

// ================================================================== scan
struct Evt { bool keep; int flav, cat; float enu, q2; double base; double sel[kNSel]; };

int Scan(const std::string &name, const std::string &outdir, Long64_t nmax, Long64_t nshape)
{
  auto it = std::find_if(kSamples.begin(), kSamples.end(),
                         [&](const Sample &s) { return name == s.name; });
  if (it == kSamples.end()) { std::cerr << "unknown sample " << name << std::endl; return 1; }
  const Sample S = *it;
  const bool isRHC = name.substr(name.size() - 3) == "rhc";

  InitGenie(true);

  // ---- flux normalisation, as FakeRecon.C
  TFile ff((kRepo + "/" + S.flux).c_str());
  TH1 *hflux = (TH1*)ff.Get(S.hist.c_str());
  const double flux_per_pot = hflux->Integral("width") / 0.05 / 1.0e21;

  // ---- selection curves, as FakeRecon.C / FakeSelect.C
  const char *resp_cat[9] = {"numuCCQE", "numubarCCQE", "numuCCnonQE", "numubarCCnonQE",
                             "nueCCQE", "nuebarCCQE", "nueCCnonQE", "nuebarCCnonQE", "NC"};
  const char *mix[3] = {"1Rmu", "1Re0de", "1Re1de"};
  TFile fr((kRepo + "/lbl_config/hk/response.root").c_str());
  TH1D *heff[9], *hmix[9][3];
  for (int c = 0; c < 9; ++c) {
    heff[c] = (TH1D*)fr.Get(Form("h_eff_%s", resp_cat[c]));
    for (int s = 0; s < 3; ++s) hmix[c][s] = (TH1D*)fr.Get(Form("h_mix_%s_%s", resp_cat[c], mix[s]));
  }
  auto lookup = [](TH1D *h, double x) {
    return h->GetBinContent(std::min(std::max(h->FindBin(x), 1), h->GetNbinsX())); };

  // ---- inputs
  TFile fg(GhepFile(name).c_str());
  TTree *gtree = (TTree*)fg.Get("gtree");
  NtpMCEventRecord *mcrec = nullptr;
  gtree->SetBranchAddress("gmcrec", &mcrec);
  TFile fl((kData + "/" + S.stem + "_flattree.root").c_str());
  TTree *flat = (TTree*)fl.Get("FlatTree_VARS");
  if (flat->GetEntries() != gtree->GetEntries()) {
    std::cerr << "[ERROR] flat/GHEP entry mismatch" << std::endl; return 1; }
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

  TFile fout((outdir + "/part_" + name + ".root").c_str(), "RECREATE");

  auto book = [&](const std::string &tag) {
    std::vector<TH1D*> h;
    for (int f = 0; f < kNFlav; ++f) for (int c = 0; c < kNCat; ++c)
      h.push_back(new TH1D(Form("h_%s_%s_%s", tag.c_str(), kFlav[f], kCat[c]), "", kNBins, kBins));
    for (int s = 0; s < kNSel; ++s)
      h.push_back(new TH1D(Form("s_%s_%s", tag.c_str(), kSel[s]), "", kNBins, kBins));
    h.push_back(new TH1D(Form("q_%s", tag.c_str()), "", kNQ2, kQ2Bins));
    for (auto x : h) { x->Sumw2(false); x->SetDirectory(nullptr); }
    return h;
  };
  auto fill = [&](std::vector<TH1D*> &h, const Evt &e, double w) {
    h[e.flav * kNCat + e.cat]->Fill(e.enu, e.base * w);
    for (int s = 0; s < kNSel; ++s) if (e.sel[s] > 0) h[kNFlav * kNCat + s]->Fill(e.enu, e.base * e.sel[s] * w);
    h[kNFlav * kNCat + kNSel]->Fill(e.q2, e.base * w);
  };
  auto save = [&](std::vector<TH1D*> &h, const std::string &dir) {
    TDirectory *d = fout.GetDirectory(dir.c_str());
    if (!d) d = fout.mkdir(dir.c_str());
    d->cd();
    for (auto x : h) { x->Write(); delete x; }
    h.clear();
  };

  // ---- nominal pass: event bookkeeping, GHEP/flat alignment check
  std::vector<Evt> evts(nev);
  std::vector<TH1D*> hnom = book("nom");
  std::vector<double> nfE(1);
  std::vector<std::array<std::array<double, 3>, 3>> nfP(1);
  Long64_t n_mismatch = 0, n_dropped = 0;
  for (Long64_t i = 0; i < nev; ++i) {
    gtree->GetEntry(i);
    flat->GetEntry(i);
    EventRecord &ev = *(mcrec->event);
    Evt &e = evts[i];
    e.enu  = ev.Probe()->E();
    e.q2   = Q2;
    e.flav = FlavIndex(ev.Probe()->Pdg());
    e.cat  = Category(ev.Summary());
    if (std::fabs(e.enu - Enu_true) > 1e-4 * std::max(1.f, Enu_true) || ev.Probe()->Pdg() != PDGnu ||
        (bool)cc != ev.Summary()->ProcInfo().IsWeakCC()) ++n_mismatch;
    mcrec->Clear();

    // FakeRecon.C drops: Pauli blocked (Q2 = -999) and CC without a lepton.
    e.keep = !(Q2 < -999.0) && !(cc && (PDGLep == 0 || ELep < 0.0));
    if (!e.keep) { ++n_dropped; continue; }

    const double pot_w = fScaleFactor * flux_per_pot * S.pot * kNucleons;
    const int a = (std::abs(S.flux_pdg) == 12) ? 0 : 1;
    const int b = (std::abs(PDGnu) == 12) ? 0 : 1;
    double osc;
    if (cc) {
      nfE[0] = (S.flux_pdg > 0) ? Enu_true : -Enu_true;
      NuFast::Probability_Matter_LBL(0.307, 0.0218, 0.528, -1.601, 7.53e-5, 2.509e-3 + 7.53e-5,
                                     295.0, nfE, 2.6 * 0.5, 0, nfP);
      osc = nfP[0][a][b];
    } else {
      osc = (a != b) ? 0.0 : 1.0;
    }
    e.base = pot_w * osc;

    int rc;
    if (!cc)                   rc = 8;
    else if (abs(PDGnu) == 12) rc = (flagCCQE ? 4 : 6) + (PDGnu > 0 ? 0 : 1);
    else                       rc = (flagCCQE ? 0 : 2) + (PDGnu > 0 ? 0 : 1);
    const double eff = lookup(heff[rc], Enu_true);
    for (int s = 0; s < kNSel; ++s) e.sel[s] = 0;
    for (int s = 0; s < 3; ++s) {
      if (isRHC && s == 2) continue;          // 1Re1de is FHC only
      e.sel[(isRHC ? 3 : 0) + s] = eff * lookup(hmix[rc][s], Enu_true);
    }
    fill(hnom, e, 1.0);
  }
  save(hnom, "nominal");
  std::cout << "[INFO] " << name << ": " << nev << " events, " << n_dropped
            << " dropped as in FakeRecon, " << n_mismatch << " GHEP/flat mismatches" << std::endl;
  if (n_mismatch > 0) { std::cerr << "[ERROR] GHEP and flat tree are not aligned" << std::endl; return 1; }

  // ---- which calculators handle which dial (decided once; a fresh calculator
  //      is only built when the dial needs a non-default mode)
  auto calcs = Calculators();
  std::map<std::string, std::vector<int>> handlers;
  {
    std::vector<GReWeightI*> shared;
    for (auto &c : calcs) shared.push_back(c.second());
    for (auto &d : kDials) {
      GSyst_t s = GSyst::FromString(d);
      for (size_t k = 0; k < shared.size(); ++k) {
        GReWeightI *c = NeedsMode(s) ? calcs[k].second() : shared[k];
        ApplyMode(c, s);
        if (c->IsHandled(s)) handlers[d].push_back(k);
        if (c != shared[k]) delete c;
      }
    }
    for (auto c : shared) delete c;
  }

  // ---- per-dial weight statistics (owned by fout)
  fout.cd();
  TTree &wstats = *new TTree("wstats", "per dial and tweak point weight statistics");
  char dname[64]; double twk, nkeep, nneq1, wmin, wmax, nneg, nnan, nbig, sumw_all; int ncalc;
  Long64_t nproc;
  wstats.Branch("dial", dname, "dial/C");
  wstats.Branch("nproc", &nproc);
  wstats.Branch("twk", &twk);
  wstats.Branch("ncalc", &ncalc);
  wstats.Branch("nkeep", &nkeep);
  wstats.Branch("nneq1", &nneq1);
  wstats.Branch("wmin", &wmin);
  wstats.Branch("wmax", &wmax);
  wstats.Branch("nneg", &nneg);
  wstats.Branch("nnan", &nnan);
  wstats.Branch("nbig", &nbig);
  wstats.Branch("sumw_all", &sumw_all);   // unweighted sum over the processed events, as rw_summary.C

  for (auto &d : kDials) {
    TStopwatch clock;
    GSyst_t s = GSyst::FromString(d);
    const int np = kPoints.size();
    const bool capped = std::find(kShapeDials.begin(), kShapeDials.end(), d) != kShapeDials.end();
    const Long64_t nrun = capped ? std::min(nev, nshape) : nev;
    const double scale = double(nev) / nrun;   // the subsample stands for the full sample
    std::vector<GReWeight*> rw(np);
    for (int p = 0; p < np; ++p) {
      rw[p] = new GReWeight;
      for (int k : handlers[d]) {
        GReWeightI *c = calcs[k].second();
        ApplyMode(c, s);
        rw[p]->AdoptWghtCalc(calcs[k].first, c);
      }
      rw[p]->Systematics().Init(s);
      rw[p]->Systematics().Set(s, kPoints[p]);
      rw[p]->Reconfigure();
    }
    std::vector<std::vector<TH1D*>> h(np);
    for (int p = 0; p < np; ++p) h[p] = book(kPtName[p]);
    std::vector<TH1D*> hsub;                    // nominal of the same subsample
    if (capped) hsub = book("nom");
    std::vector<double> vneq1(np, 0), vmin(np, 1e30), vmax(np, -1e30), vneg(np, 0),
                        vnan(np, 0), vbig(np, 0), vsum(np, 0);

    for (Long64_t i = 0; i < nrun; ++i) {
      gtree->GetEntry(i);
      if (capped && evts[i].keep) fill(hsub, evts[i], scale);
      EventRecord &ev = *(mcrec->event);
      // As grwght1p: set Tl and ctl for MEC reweighting.
      Interaction *in = ev.Summary();
      const TLorentzVector &p4l = in->KinePtr()->FSLeptonP4();
      in->KinePtr()->SetKV(kKVTl, p4l.E() - in->FSPrimLepton()->Mass());
      in->KinePtr()->SetKV(kKVctl, p4l.CosTheta());
      for (int p = 0; p < np; ++p) {
        const double w = handlers[d].empty() ? 1.0 : rw[p]->CalcWeight(ev);
        if (std::isnan(w)) { ++vnan[p]; continue; }
        vsum[p] += w;
        if (!evts[i].keep) continue;
        if (std::fabs(w - 1.0) > 1e-9) ++vneq1[p];
        vmin[p] = std::min(vmin[p], w);
        vmax[p] = std::max(vmax[p], w);
        if (w < 0)  ++vneg[p];
        if (w > 10) ++vbig[p];
        fill(h[p], evts[i], w * scale);
      }
      mcrec->Clear();
    }
    if (capped) save(hsub, d);
    Long64_t nkeep_run = 0;
    for (Long64_t i = 0; i < nrun; ++i) nkeep_run += evts[i].keep;
    for (int p = 0; p < np; ++p) {
      save(h[p], d);
      snprintf(dname, sizeof(dname), "%s", d.c_str());
      twk = kPoints[p]; ncalc = handlers[d].size(); nproc = nrun;
      nkeep = nkeep_run; nneq1 = vneq1[p]; wmin = vmin[p]; wmax = vmax[p];
      nneg = vneg[p]; nnan = vnan[p]; nbig = vbig[p]; sumw_all = vsum[p];
      wstats.Fill();
      delete rw[p];
    }
    std::cout << "[INFO] " << name << " " << d << ": " << handlers[d].size() << " calculator(s), w!=1 in "
              << vneq1[2] << " events at +1 sigma, w in [" << vmin[2] << ", " << vmax[2]
              << "], <w>(-1,+1) = " << vsum[1] / nrun << ", " << vsum[2] / nrun
              << " over " << nrun << " events, " << clock.RealTime() << " s" << std::endl;
  }
  fout.cd();
  wstats.Write();
  fout.Close();
  return 0;
}

// ============================================================= summarize
// Sum a histogram name over all part files.
struct Parts {
  std::vector<TFile*> f;
  TH1D *Get(const std::string &dir, const std::string &h) {
    TH1D *sum = nullptr;
    for (auto x : f) {
      TH1D *a = (TH1D*)x->Get((dir + "/" + h).c_str());
      if (!a) continue;
      if (!sum) { sum = (TH1D*)a->Clone(); sum->SetDirectory(nullptr); } else sum->Add(a);
    }
    return sum;
  }
};

// Spread of an Enu ratio around its mean, relative to the size of the change:
// 0 = pure normalisation, 1 = pure shape (rate-preserving).
double ShapeIndex(TH1D *nom, TH1D *var)
{
  const int nb = nom->GetNbinsX();
  const double N = nom->Integral(0, nb + 1), V = var->Integral(0, nb + 1);
  if (N <= 0) return -1;
  const double R = V / N;
  double s2 = 0, sn = 0;
  for (int b = 1; b <= nb; ++b) {
    const double n = nom->GetBinContent(b);
    if (n < 0.005 * N) continue;               // ignore near-empty bins
    const double r = var->GetBinContent(b) / n;
    s2 += n * (r - R) * (r - R); sn += n;
  }
  if (sn <= 0) return -1;
  const double rms = std::sqrt(s2 / sn);
  if (rms + std::fabs(R - 1) < 1e-6) return -1;   // no effect
  return rms / (rms + std::fabs(R - 1));
}

int Summarize(const std::string &outdir)
{
  InitGenie(false);
  GSystUncertainty *unc = GSystUncertainty::Instance();
  gStyle->SetOptStat(0);

  Parts P;
  for (auto &s : kSamples) {
    TFile *f = TFile::Open((outdir + "/part_" + s.name + ".root").c_str());
    if (!f || f->IsZombie()) { std::cerr << "[ERROR] missing part file for " << s.name << std::endl; return 1; }
    P.f.push_back(f);
  }

  TFile fout((outdir + "/genie_dial_sensitivity.root").c_str(), "RECREATE");

  // Sum the per-(flav,cat) histograms of one tag into total / per-flav / per-cat / per-sel.
  struct Sums { TH1D *tot; TH1D *flav[kNFlav]; TH1D *cat[kNCat]; TH1D *sel[kNSel]; TH1D *q; };
  auto sums = [&](const std::string &dir, const std::string &tag) {
    Sums S{};
    for (int f = 0; f < kNFlav; ++f) for (int c = 0; c < kNCat; ++c) {
      TH1D *h = P.Get(dir, Form("h_%s_%s_%s", tag.c_str(), kFlav[f], kCat[c]));
      if (!S.tot)     S.tot     = (TH1D*)h->Clone(Form("%s_total", tag.c_str()));
      else            S.tot->Add(h);
      if (!S.flav[f]) S.flav[f] = (TH1D*)h->Clone(Form("%s_%s", tag.c_str(), kFlav[f]));
      else            S.flav[f]->Add(h);
      if (!S.cat[c])  S.cat[c]  = (TH1D*)h->Clone(Form("%s_%s", tag.c_str(), kCat[c]));
      else            S.cat[c]->Add(h);
      delete h;
    }
    for (int s = 0; s < kNSel; ++s) {
      S.sel[s] = P.Get(dir, Form("s_%s_%s", tag.c_str(), kSel[s]));
      S.sel[s]->SetName(Form("%s_%s", tag.c_str(), kSel[s]));
    }
    S.q = P.Get(dir, Form("q_%s", tag.c_str()));
    S.q->SetName(Form("%s_Q2", tag.c_str()));
    S.tot->SetDirectory(nullptr);
    for (auto h : S.flav) h->SetDirectory(nullptr);
    for (auto h : S.cat)  h->SetDirectory(nullptr);
    return S;
  };
  auto I = [](TH1D *h) { return h->Integral(0, kNBins + 1); };
  auto pct = [](double v, double n) { return n > 0 ? 100.0 * (v / n - 1.0) : 0.0; };

  Sums nom = sums("nominal", "nom");
  TDirectory *dn = fout.mkdir("nominal");
  dn->cd();
  nom.tot->Write();
  for (auto h : nom.flav) h->Write();
  for (auto h : nom.cat)  h->Write();
  for (auto h : nom.sel)  h->Write();
  nom.q->Write();

  // ---- summary tree (owned by fout)
  fout.cd();
  TTree &tree = *new TTree("summary", "rate change in percent relative to nominal");
  char dname[64]; double twk, sig, dtot, dflav[kNFlav], dcat[kNCat], dsel[kNSel];
  double wmin, wmax, nneq1, nneg, nnan, nbig, shape_tot, shape_main, shape_q2; int ncalc;
  Long64_t nproc;
  tree.Branch("dial", dname, "dial/C");
  tree.Branch("twk", &twk);
  tree.Branch("sigma_frac", &sig);
  tree.Branch("dtot", &dtot);
  tree.Branch("dflav", dflav, Form("dflav[%d]/D", kNFlav));
  tree.Branch("dcat", dcat, Form("dcat[%d]/D", kNCat));
  tree.Branch("dsel", dsel, Form("dsel[%d]/D", kNSel));
  tree.Branch("ncalc", &ncalc);
  tree.Branch("nneq1", &nneq1);
  tree.Branch("wmin", &wmin);
  tree.Branch("wmax", &wmax);
  tree.Branch("nneg", &nneg);
  tree.Branch("nnan", &nnan);
  tree.Branch("nbig", &nbig);
  tree.Branch("shape_tot", &shape_tot);
  tree.Branch("shape_main", &shape_main);
  tree.Branch("shape_q2", &shape_q2);
  tree.Branch("nproc", &nproc);   // events processed, summed over samples

  struct Row { std::string dial; double sp, sm, d[4], flav[2][kNFlav], cat[2][kNCat], sel[2][kNSel];
               std::string main, type, q2type, flags; double nneq1; double wmin, wmax; };
  std::vector<Row> rows;

  gStyle->SetPaperSize(28, 20);   // same aspect as the canvas, so the page is filled
  TCanvas cv("cv", "", 1400, 1000);
  const std::string pdf = outdir + "/ratio_plots.pdf";
  cv.Print((pdf + "[").c_str());
  const int col[4] = {kBlue + 1, kAzure + 7, kOrange + 7, kRed + 1};

  for (auto &d : kDials) {
    GSyst_t s = GSyst::FromString(d);
    Row r{}; r.dial = d;
    r.sp = unc->OneSigmaErr(s, +1); r.sm = unc->OneSigmaErr(s, -1);
    TDirectory *dd = fout.mkdir(d.c_str());

    // weight statistics, merged over samples
    double st_neq1[4] = {0}, st_min[4], st_max[4], st_neg[4] = {0}, st_nan[4] = {0}, st_big[4] = {0};
    int st_ncalc = 0; Long64_t st_nproc = 0;
    for (int p = 0; p < 4; ++p) { st_min[p] = 1e30; st_max[p] = -1e30; }
    for (auto f : P.f) {
      TTree *w = (TTree*)f->Get("wstats");
      char n[64]; double t, a, b, c, e, g, h; int nc; Long64_t np;
      w->SetBranchAddress("nproc", &np);
      w->SetBranchAddress("dial", n); w->SetBranchAddress("twk", &t); w->SetBranchAddress("ncalc", &nc);
      w->SetBranchAddress("nneq1", &a); w->SetBranchAddress("wmin", &b); w->SetBranchAddress("wmax", &c);
      w->SetBranchAddress("nneg", &e); w->SetBranchAddress("nnan", &g); w->SetBranchAddress("nbig", &h);
      for (Long64_t j = 0; j < w->GetEntries(); ++j) {
        w->GetEntry(j);
        if (d != n) continue;
        int p = std::find(kPoints.begin(), kPoints.end(), t) - kPoints.begin();
        st_neq1[p] += a; st_min[p] = std::min(st_min[p], b); st_max[p] = std::max(st_max[p], c);
        st_neg[p] += e; st_nan[p] += g; st_big[p] += h; st_ncalc = nc;
        if (p == 0) st_nproc += np;
      }
      w->ResetBranchAddresses();
    }

    // Subsampled (shape) dials carry their own nominal over the same events.
    Sums loc{};
    TH1D *probe = P.Get(d, "q_nom");
    if (probe) { delete probe; loc = sums(d, "nom"); }
    const Sums &N = probe ? loc : nom;

    std::vector<Sums> var;
    for (int p = 0; p < 4; ++p) var.push_back(sums(d, kPtName[p]));

    // main categories: share of the absolute rate change at +1 sigma
    std::vector<std::pair<double, int>> ch;
    double chsum = 0;
    for (int c = 0; c < kNCat; ++c) {
      const double x = std::fabs(I(var[2].cat[c]) - I(N.cat[c]));
      ch.push_back({x, c}); chsum += x;
    }
    std::sort(ch.rbegin(), ch.rend());
    for (auto &x : ch) if (chsum > 0 && x.first >= 0.2 * chsum)
      r.main += std::string(r.main.empty() ? "" : ", ") + kCat[x.second];
    const int cmain = ch[0].second;

    for (int p = 0; p < 4; ++p) {
      snprintf(dname, sizeof(dname), "%s", d.c_str());
      twk = kPoints[p];
      sig = kPoints[p] > 0 ? r.sp : r.sm;
      dtot = pct(I(var[p].tot), I(N.tot));
      for (int f = 0; f < kNFlav; ++f) dflav[f] = pct(I(var[p].flav[f]), I(N.flav[f]));
      for (int c = 0; c < kNCat; ++c)  dcat[c]  = pct(I(var[p].cat[c]),  I(N.cat[c]));
      for (int q = 0; q < kNSel; ++q)  dsel[q]  = pct(I(var[p].sel[q]),  I(N.sel[q]));
      ncalc = st_ncalc; nneq1 = st_neq1[p]; wmin = st_min[p]; wmax = st_max[p];
      nneg = st_neg[p]; nnan = st_nan[p]; nbig = st_big[p];
      shape_tot  = ShapeIndex(N.tot, var[p].tot);
      shape_main = ShapeIndex(N.cat[cmain], var[p].cat[cmain]);
      shape_q2   = ShapeIndex(N.q, var[p].q);
      nproc      = st_nproc;
      tree.Fill();

      r.d[p] = dtot;
      if (p == 1 || p == 2) {
        const int k = p == 1 ? 0 : 1;
        for (int f = 0; f < kNFlav; ++f) r.flav[k][f] = dflav[f];
        for (int c = 0; c < kNCat; ++c)  r.cat[k][c]  = dcat[c];
        for (int q = 0; q < kNSel; ++q)  r.sel[k][q]  = dsel[q];
      }
      if (p == 2) {
        r.nneq1 = nneq1;
        if (shape_main < 0) r.type = "none";
        else r.type = shape_main < 0.15 ? "norm" : shape_main > 0.6 ? "shape" : "norm+shape";
        r.q2type = shape_q2 < 0 ? "-" : Form("%.2f", shape_q2);
      }
    }
    r.wmin = *std::min_element(st_min, st_min + 4);
    r.wmax = *std::max_element(st_max, st_max + 4);
    if (probe) r.flags += Form("subsample: %lld events over 12 samples; ", st_nproc);
    if (st_ncalc == 0) r.flags += "no calculator handles it; ";
    if (st_neq1[1] + st_neq1[2] == 0) r.flags += "all weights = 1; ";
    for (int p = 0; p < 4; ++p) {
      if (st_neg[p] > 0) r.flags += Form("w<0 at %+gs (%.0f ev); ", kPoints[p], st_neg[p]);
      if (st_nan[p] > 0) r.flags += Form("NaN at %+gs (%.0f ev); ", kPoints[p], st_nan[p]);
      if (st_big[p] > 0) r.flags += Form("w>10 at %+gs (%.0f ev); ", kPoints[p], st_big[p]);
    }
    rows.push_back(r);

    // ---- ratio histograms
    dd->cd();
    auto ratio = [&](TH1D *v, TH1D *n, const std::string &nm) {
      TH1D *x = (TH1D*)v->Clone(nm.c_str());
      x->Divide(n);
      x->SetTitle(Form("%s;E_{#nu}^{true} [GeV];ratio to nominal", d.c_str()));
      x->Write();
      return x;
    };
    std::vector<TH1D*> rtot(4), rcat(kNCat), rsel(kNSel), rq2(4);
    for (int p = 0; p < 4; ++p) {
      var[p].tot->SetName(Form("rate_total_%s", kPtName[p])); var[p].tot->Write();
      rtot[p] = ratio(var[p].tot, N.tot, Form("ratio_total_%s", kPtName[p]));
      for (int f = 0; f < kNFlav; ++f) ratio(var[p].flav[f], N.flav[f], Form("ratio_%s_%s", kFlav[f], kPtName[p]));
      for (int c = 0; c < kNCat; ++c) {
        TH1D *x = ratio(var[p].cat[c], N.cat[c], Form("ratio_%s_%s", kCat[c], kPtName[p]));
        if (p == 2) rcat[c] = x;
      }
      for (int q = 0; q < kNSel; ++q) {
        TH1D *x = ratio(var[p].sel[q], N.sel[q], Form("ratio_%s_%s", kSel[q], kPtName[p]));
        if (p == 2) rsel[q] = x;
      }
      rq2[p] = ratio(var[p].q, N.q, Form("ratio_Q2_%s", kPtName[p]));
      rq2[p]->GetXaxis()->SetTitle("Q^{2}_{true} [GeV^{2}]");
      rq2[p]->Write("", TObject::kOverwrite);
    }

    // ---- one PDF page per dial
    cv.Clear(); cv.Divide(2, 2);
    auto frame = [&](int pad, const char *title, std::vector<TH1D*> hs) {
      cv.cd(pad); gPad->SetLogx(false); gPad->SetGridy();
      double lo = 1, hi = 1;
      for (auto h : hs) for (int b = 1; b <= h->GetNbinsX(); ++b) if (h->GetBinContent(b) != 0) {
        lo = std::min(lo, h->GetBinContent(b)); hi = std::max(hi, h->GetBinContent(b)); }
      const double m = std::max(0.002, 0.15 * (hi - lo));
      hs[0]->SetTitle(title);
      hs[0]->GetYaxis()->SetRangeUser(lo - m, hi + m);
      if (hs[0]->GetXaxis()->GetXmax() > 5) hs[0]->GetXaxis()->SetRangeUser(0, 5);
    };
    frame(1, Form("%s: all events;E_{#nu}^{true} [GeV];ratio to nominal", d.c_str()), rtot);
    TLegend *l1 = new TLegend(0.65, 0.7, 0.88, 0.88);
    for (int p = 0; p < 4; ++p) {
      rtot[p]->SetLineColor(col[p]); rtot[p]->SetLineWidth(2);
      rtot[p]->Draw(p ? "hist same" : "hist");
      l1->AddEntry(rtot[p], Form("%+g#sigma (%+.2f%%)", kPoints[p], r.d[p]), "l");
    }
    l1->Draw();
    std::vector<TH1D*> shown;
    for (int c = 0; c < kNCat; ++c)
      if (I(N.cat[c]) > 1e-3 * I(N.tot) && std::fabs(I(var[2].cat[c]) / I(N.cat[c]) - 1) > 1e-5)
        shown.push_back(rcat[c]);
    TLegend *l2 = new TLegend(0.65, 0.6, 0.88, 0.88);
    if (!shown.empty()) {
      frame(2, Form("+1#sigma by interaction type;E_{#nu}^{true} [GeV];ratio to nominal"), shown);
      for (size_t k = 0; k < shown.size(); ++k) {
        shown[k]->SetLineColor(k + 1 < 5 ? k + 1 : k + 2); shown[k]->SetLineWidth(2);
        shown[k]->Draw(k ? "hist same" : "hist");
        std::string nm = shown[k]->GetName();
        l2->AddEntry(shown[k], nm.substr(6, nm.size() - 9).c_str(), "l");
      }
      l2->Draw();
    }
    frame(3, Form("+1#sigma by analysis sample;E_{#nu}^{true} [GeV];ratio to nominal"), rsel);
    TLegend *l3 = new TLegend(0.6, 0.65, 0.88, 0.88);
    for (int q = 0; q < kNSel; ++q) {
      rsel[q]->SetLineColor(col[q % 4] + (q >= 4 ? 2 : 0)); rsel[q]->SetLineStyle(q < 3 ? 1 : 2);
      rsel[q]->SetLineWidth(2);
      rsel[q]->Draw(q ? "hist same" : "hist");
      l3->AddEntry(rsel[q], kSel[q], "l");
    }
    l3->Draw();
    frame(4, "all events vs true Q^{2};Q^{2}_{true} [GeV^{2}];ratio to nominal", rq2);
    for (int p = 0; p < 4; ++p) {
      rq2[p]->SetLineColor(col[p]); rq2[p]->SetLineWidth(2);
      rq2[p]->Draw(p ? "hist same" : "hist");
    }
    cv.Print(pdf.c_str());
    delete l1; delete l2; delete l3;
  }
  cv.Print((pdf + "]").c_str());
  fout.cd();
  tree.Write();
  fout.Close();

  // ---- markdown tables
  std::sort(rows.begin(), rows.end(), [](const Row &a, const Row &b) {
    return std::max(std::fabs(a.d[1]), std::fabs(a.d[2])) > std::max(std::fabs(b.d[1]), std::fabs(b.d[2])); });
  printf("\n### Ranking by total event rate\n\n");
  printf("| # | dial | 1sigma (+/-) | -2s | -1s | +1s | +2s | main interaction types | Enu type | Q2 shape idx | weight range | flags |\n");
  printf("|---|---|---|---|---|---|---|---|---|---|---|---|\n");
  int k = 0;
  for (auto &r : rows)
    printf("| %d | %s | %.3g/%.3g | %+.3f%% | %+.3f%% | %+.3f%% | %+.3f%% | %s | %s | %s | [%.3g, %.3g] | %s |\n",
           ++k, r.dial.c_str(), r.sp, r.sm, r.d[0], r.d[1], r.d[2], r.d[3],
           r.main.empty() ? "-" : r.main.c_str(), r.type.c_str(), r.q2type.c_str(), r.wmin, r.wmax, r.flags.c_str());

  printf("\n### By neutrino flavour, -1s / +1s (%%)\n\n| dial |");
  for (int f = 0; f < kNFlav; ++f) printf(" %s |", kFlav[f]);
  printf("\n|---|"); for (int f = 0; f < kNFlav; ++f) printf("---|"); printf("\n");
  for (auto &r : rows) {
    printf("| %s |", r.dial.c_str());
    for (int f = 0; f < kNFlav; ++f) printf(" %+.2f / %+.2f |", r.flav[0][f], r.flav[1][f]);
    printf("\n");
  }
  printf("\n### By interaction type, -1s / +1s (%%)\n\n| dial |");
  for (int c = 0; c < kNCat; ++c) printf(" %s |", kCat[c]);
  printf("\n|---|"); for (int c = 0; c < kNCat; ++c) printf("---|"); printf("\n");
  for (auto &r : rows) {
    printf("| %s |", r.dial.c_str());
    for (int c = 0; c < kNCat; ++c) {
      if (r.cat[0][c] == 0 && r.cat[1][c] == 0) printf(" 0 |");
      else printf(" %+.2f / %+.2f |", r.cat[0][c], r.cat[1][c]);
    }
    printf("\n");
  }
  printf("\n### By analysis sample, -1s / +1s (%%)\n\n| dial |");
  for (int q = 0; q < kNSel; ++q) printf(" %s |", kSel[q]);
  printf("\n|---|"); for (int q = 0; q < kNSel; ++q) printf("---|"); printf("\n");
  for (auto &r : rows) {
    printf("| %s |", r.dial.c_str());
    for (int q = 0; q < kNSel; ++q) printf(" %+.2f / %+.2f |", r.sel[0][q], r.sel[1][q]);
    printf("\n");
  }

  printf("\n### Nominal rates (events, all 12 samples, oscillated, no selection)\n\n");
  printf("| total | %.1f |\n|---|---|\n", I(nom.tot));
  for (int f = 0; f < kNFlav; ++f) printf("| %s | %.1f |\n", kFlav[f], I(nom.flav[f]));
  for (int c = 0; c < kNCat; ++c)  printf("| %s | %.1f |\n", kCat[c], I(nom.cat[c]));
  for (int q = 0; q < kNSel; ++q)  printf("| %s (selected) | %.1f |\n", kSel[q], I(nom.sel[q]));
  return 0;
}

int main(int argc, char **argv)
{
  const std::string mode = argc > 1 ? argv[1] : "";
  if ((mode == "scan" && argc >= 6) || (mode == "summarize" && argc >= 5)) {
    kRepo = argv[2];
    kData = argv[3];
    if (!LoadSamples()) return 1;
  }
  if (mode == "scan" && argc >= 6)
    return Scan(argv[4], argv[5], argc > 6 ? atoll(argv[6]) : -1, argc > 7 ? atoll(argv[7]) : 5000);
  if (mode == "summarize" && argc >= 5)
    return Summarize(argv[4]);
  std::cerr << "usage: dial_scan scan <working_dir> <data_dir> <sample> <outdir> [nevents] [nshape]\n"
               "       dial_scan summarize <working_dir> <data_dir> <outdir>" << std::endl;
  return 1;
}

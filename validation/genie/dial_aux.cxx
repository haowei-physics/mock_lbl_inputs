// =============================================================================
//  dial_aux.cxx  -  auxiliary checks on the dial_scan.cxx output, for docs/genie/report.md.
//
//    dial_aux <working_dir> <data_dir> <sample> <outfile> [nshape] [nmax]
//        (nmax: test on the first events only; <working_dir>, <data_dir>: see
//        kRepo, kData below)
//
//  Recomputes, for one sample, the weights of every dial that dial_scan.cxx runs
//  on the full sample (all but the 5 subsampled shape dials), with the same tune,
//  calculators, mode switches, MEC Tl/ctl fix and FakeRecon normalisation, and
//  writes plain-text lines:
//
//    Z <sample> <dial> <twk> nkeep nnan nneg nzero nbig nsmall
//        weight counts (nzero: w == 0 exactly; nsmall: 0 < w < 0.1)
//    R <sample> <dial> <twk> <scope> S0 S1 S2  [cat: S0 S1 S2]x12  [sel: S0 S1 S2]x5
//        S0 = sum b, S1 = sum b (w-1), S2 = sum b^2 (w-1)^2 over the kept events,
//        b = FakeRecon event weight (x selection expectation for sel). scope = all
//        or sub (first nshape events, scaled to the sample as dial_scan does).
//    N <sample> <dial> <twk> <entry> w Enu Q2 Z hitnuc phasespace dxsec
//        one line per negative weight
//
//  Extra tweak points beyond -2,-1,+1,+2: +-0.01 sigma for the dials whose
//  calculator compares with a dxsec it does not recompute (offset test), and the
//  analogue points of the 5 shape dials (same physical shift as 1 sigma of the
//  shape dial: MaCCQE +-0.833, Ma/MvCCRES and Ma/MvNCRES +-0.5).
//  Code copied from dial_scan.cxx; that file is not modified.
// =============================================================================

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
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

#define main nufast_unused_main
#include "NuFast_LBL.cpp"
#undef main

using namespace genie;
using namespace genie::rew;

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

// dial -> tweak points (the full-sample dials of dial_scan.cxx)
std::vector<std::pair<std::string, std::vector<double>>> Plan()
{
  const std::vector<double> std4 = {-2, -1, 1, 2};
  const std::vector<double> tiny = {-2, -1, -0.01, 0.01, 1, 2};
  const std::vector<std::string> offset = {
    "AhtBY", "BhtBY", "CV1uBY", "CV2uBY", "AhtBYshape", "BhtBYshape", "CV1uBYshape", "CV2uBYshape",
    "MaCOHpi", "R0COHpi", "MaNCEL", "EtaNCEL", "CoulombCCQE", "CCQEPauliSupViaKF",
    "FracPN_CCMEC", "FracDelta_CCMEC"};
  const std::vector<std::string> all = {
    "MaNCEL", "EtaNCEL", "MaCCQE", "NormCCQE", "CoulombCCQE", "CCQEPauliSupViaKF",
    "NormCCRES", "MaCCRES", "MvCCRES", "NormNCRES", "MaNCRES", "MvNCRES",
    "MaCOHpi", "R0COHpi", "NormCCCOH", "NormNCCOH",
    "NonRESBGvpCC1pi", "NonRESBGvpCC2pi", "NonRESBGvpNC1pi", "NonRESBGvpNC2pi",
    "NonRESBGvnCC1pi", "NonRESBGvnCC2pi", "NonRESBGvnNC1pi", "NonRESBGvnNC2pi",
    "NonRESBGvbarpCC1pi", "NonRESBGvbarpCC2pi", "NonRESBGvbarpNC1pi", "NonRESBGvbarpNC2pi",
    "NonRESBGvbarnCC1pi", "NonRESBGvbarnCC2pi", "NonRESBGvbarnNC1pi", "NonRESBGvbarnNC2pi",
    "AhtBY", "BhtBY", "CV1uBY", "CV2uBY",
    "AhtBYshape", "BhtBYshape", "CV1uBYshape", "CV2uBYshape",
    "RDecBR1gamma", "RDecBR1eta",
    "NormCCMEC", "NormNCMEC", "FracPN_CCMEC", "FracDelta_CCMEC",
    "AGKYxF1pi", "AGKYpT1pi"};   // slowest last
  std::vector<std::pair<std::string, std::vector<double>>> plan;
  for (auto &d : all) {
    std::vector<double> p = std::find(offset.begin(), offset.end(), d) != offset.end() ? tiny : std4;
    if (d == "MaCCQE") p = {-2, -1, -0.833, -0.01, 0.01, 0.833, 1, 2};
    if (d == "MaCCRES" || d == "MvCCRES" || d == "MaNCRES" || d == "MvNCRES")
      p = {-2, -1, -0.5, -0.01, 0.01, 0.5, 1, 2};
    plan.push_back({d, p});
  }
  return plan;
}

const int kNCat = 12;
const int kNSel = 5;

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

struct Evt { bool keep; int cat; double base; double sel[kNSel]; };

// sum b, sum b(w-1), sum b^2(w-1)^2
struct Acc {
  double s[3] = {0, 0, 0};
  void add(double b, double w) { s[0] += b; s[1] += b * (w - 1); s[2] += b * b * (w - 1) * (w - 1); }
};
struct Sums { Acc tot, cat[kNCat], sel[kNSel]; };

int main(int argc, char **argv)
{
  if (argc < 5) { std::cerr << "usage: dial_aux <working_dir> <data_dir> <sample> <outfile> [nshape] [nmax]" << std::endl; return 1; }
  kRepo = argv[1];
  kData = argv[2];
  if (!LoadSamples()) return 1;
  const std::string name = argv[3];
  const Long64_t nshape = argc > 5 ? atoll(argv[5]) : 5000;
  auto it = std::find_if(kSamples.begin(), kSamples.end(),
                         [&](const Sample &s) { return name == s.name; });
  if (it == kSamples.end()) { std::cerr << "unknown sample " << name << std::endl; return 1; }
  const Sample S = *it;
  const bool isRHC = name.substr(name.size() - 3) == "rhc";

  utils::app_init::MesgThresholds(MsgXml());
  utils::app_init::RandGen(12345);
  GHepRecord::SetPrintLevel(0);
  RunOpt::Instance()->SetTuneName(kTune);
  RunOpt::Instance()->BuildTune();
  XSecSplineList::Instance()->LoadFromXml(SplinesXml());
  Messenger::Instance()->SetPriorityLevel("ReW", pWARN);

  TFile ff((kRepo + "/" + S.flux).c_str());
  TH1 *hflux = (TH1*)ff.Get(S.hist.c_str());
  const double flux_per_pot = hflux->Integral("width") / 0.05 / 1.0e21;

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

  TFile fg(GhepFile(name).c_str());
  TTree *gtree = (TTree*)fg.Get("gtree");
  NtpMCEventRecord *mcrec = nullptr;
  gtree->SetBranchAddress("gmcrec", &mcrec);
  TFile fl((kData + "/" + S.stem + "_flattree.root").c_str());
  TTree *flat = (TTree*)fl.Get("FlatTree_VARS");
  if (flat->GetEntries() != gtree->GetEntries()) { std::cerr << "[ERROR] flat/GHEP entry mismatch" << std::endl; return 1; }
  Long64_t nev = gtree->GetEntries();
  if (argc > 6 && atoll(argv[6]) > 0) nev = std::min(nev, atoll(argv[6]));

  double fScaleFactor; float Enu_true, ELep, Q2; int PDGnu, PDGLep; Char_t cc; Bool_t flagCCQE;
  flat->SetBranchAddress("fScaleFactor", &fScaleFactor);
  flat->SetBranchAddress("Enu_true", &Enu_true);
  flat->SetBranchAddress("ELep", &ELep);
  flat->SetBranchAddress("Q2", &Q2);
  flat->SetBranchAddress("PDGnu", &PDGnu);
  flat->SetBranchAddress("PDGLep", &PDGLep);
  flat->SetBranchAddress("cc", &cc);
  flat->SetBranchAddress("flagCCQE", &flagCCQE);

  std::vector<Evt> evts(nev);
  std::vector<float> enu(nev), q2(nev);
  std::vector<double> nfE(1);
  std::vector<std::array<std::array<double, 3>, 3>> nfP(1);
  for (Long64_t i = 0; i < nev; ++i) {
    gtree->GetEntry(i);
    flat->GetEntry(i);
    Evt &e = evts[i];
    e.cat = Category(mcrec->event->Summary());
    enu[i] = Enu_true; q2[i] = Q2;
    mcrec->Clear();
    e.keep = !(Q2 < -999.0) && !(cc && (PDGLep == 0 || ELep < 0.0));
    if (!e.keep) continue;
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
      if (isRHC && s == 2) continue;
      e.sel[(isRHC ? 3 : 0) + s] = eff * lookup(hmix[rc][s], Enu_true);
    }
  }

  FILE *out = fopen(argv[4], "w");
  auto calcs = Calculators();
  const double scale = double(nev) / std::min(nev, nshape);

  for (auto &dp : Plan()) {
    const std::string &d = dp.first;
    const std::vector<double> &pts = dp.second;
    GSyst_t s = GSyst::FromString(d);
    std::vector<int> handlers;
    for (size_t k = 0; k < calcs.size(); ++k) {
      GReWeightI *c = calcs[k].second();
      ApplyMode(c, s);
      if (c->IsHandled(s)) handlers.push_back(k);
      delete c;
    }
    const int np = pts.size();
    std::vector<GReWeight*> rw(np);
    for (int p = 0; p < np; ++p) {
      rw[p] = new GReWeight;
      for (int k : handlers) {
        GReWeightI *c = calcs[k].second();
        ApplyMode(c, s);
        rw[p]->AdoptWghtCalc(calcs[k].first, c);
      }
      rw[p]->Systematics().Init(s);
      rw[p]->Systematics().Set(s, pts[p]);
      rw[p]->Reconfigure();
    }
    std::vector<Sums> all(np), sub(np);
    std::vector<std::array<long, 6>> cnt(np, {0, 0, 0, 0, 0, 0});   // keep nan neg zero big small

    for (Long64_t i = 0; i < nev; ++i) {
      gtree->GetEntry(i);
      EventRecord &ev = *(mcrec->event);
      Interaction *in = ev.Summary();
      const TLorentzVector &p4l = in->KinePtr()->FSLeptonP4();
      in->KinePtr()->SetKV(kKVTl, p4l.E() - in->FSPrimLepton()->Mass());
      in->KinePtr()->SetKV(kKVctl, p4l.CosTheta());
      const Evt &e = evts[i];
      for (int p = 0; p < np; ++p) {
        const double w = handlers.empty() ? 1.0 : rw[p]->CalcWeight(ev);
        if (!e.keep) continue;
        auto &c = cnt[p];
        ++c[0];
        if (std::isnan(w)) { ++c[1]; continue; }
        if (w < 0)  ++c[2];
        if (w == 0) ++c[3];
        if (w > 10) ++c[4];
        if (w > 0 && w < 0.1) ++c[5];
        if (w < 0)
          fprintf(out, "N %s %s %+g %lld %.6g %.4f %.4f %d %d %d %.6g\n", name.c_str(), d.c_str(), pts[p],
                  i, w, enu[i], q2[i], in->InitState().Tgt().Z(), in->InitState().Tgt().HitNucPdg(),
                  (int)ev.DiffXSecVars(), ev.DiffXSec());
        Sums *dst[2] = {&all[p], i < nshape ? &sub[p] : nullptr};
        for (int k = 0; k < 2; ++k) {
          if (!dst[k]) continue;
          const double f = k ? scale : 1.0;
          dst[k]->tot.add(e.base * f, w);
          dst[k]->cat[e.cat].add(e.base * f, w);
          for (int q = 0; q < kNSel; ++q) if (e.sel[q] > 0) dst[k]->sel[q].add(e.base * e.sel[q] * f, w);
        }
      }
      mcrec->Clear();
    }
    for (int p = 0; p < np; ++p) {
      auto &c = cnt[p];
      fprintf(out, "Z %s %s %+g %ld %ld %ld %ld %ld %ld\n", name.c_str(), d.c_str(), pts[p],
              c[0], c[1], c[2], c[3], c[4], c[5]);
      for (int k = 0; k < 2; ++k) {
        const Sums &x = k ? sub[p] : all[p];
        fprintf(out, "R %s %s %+g %s %.8g %.8g %.8g", name.c_str(), d.c_str(), pts[p], k ? "sub" : "all",
                x.tot.s[0], x.tot.s[1], x.tot.s[2]);
        for (int q = 0; q < kNCat; ++q) fprintf(out, " %.8g %.8g %.8g", x.cat[q].s[0], x.cat[q].s[1], x.cat[q].s[2]);
        for (int q = 0; q < kNSel; ++q) fprintf(out, " %.8g %.8g %.8g", x.sel[q].s[0], x.sel[q].s[1], x.sel[q].s[2]);
        fprintf(out, "\n");
      }
      delete rw[p];
    }
    fflush(out);
    std::cout << "[INFO] " << name << " " << d << " done" << std::endl;
  }
  fclose(out);
  std::cout << "[INFO] " << name << " finished" << std::endl;
  return 0;
}

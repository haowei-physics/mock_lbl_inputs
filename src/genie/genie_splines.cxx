// =============================================================================
//  genie_splines.cxx  -  per-event GENIE Reweight response graphs for GUNDAM,
//                        for the stock G18_10b_02_11b HK samples.
//
//  Usage:
//      genie_splines [--interp] <ghep.root> <flattree.root> <fakerecon.root>
//                    <xsec.xml> <messenger.xml> <out.root> [nmax]
//  (--interp: check mode, see kInterpPts)
//  Reads
//      <ghep.root>       GHEP events (scripts/genie/gen.sh), for reweighting
//      <flattree.root>   their flat tree (scripts/genie/flat.sh), for alignment
//      <fakerecon.root>  FakeRecon output (scripts/recon.sh)
//      <xsec.xml>        GENIE cross-section splines (scripts/genie/xsec.sh)
//      <messenger.xml>   message thresholds (gen/genie/Messenger_production.xml)
//  Writes <out.root>, in production <data_dir>/<stem>_spline.root:
//          FlatTree_VARS : clone of the fakerecon tree + one TGraph branch per dial
//          GenieDialMeta : one entry per dial (nominal, sigma, knots, bounds, counts)
//  <stem> is the NEUT-style name (numu/numu_x_numu_200K_FHC, ...; see
//  gen/genie/samples.txt), so that src/FakeSelect.C reads the output unchanged
//  with input_suffix "_spline.root". scripts/genie/splines.sh runs this for
//  every sample.
//
//  ---------------------------------------------------------------------------
//  Decision log (agreed 2026-10-06; see docs/genie/report.md)
//
//  * Follow GENIE's own settings, not the T2KReWeight/NEUT parameter set.
//    Central values: the stock G18_10b_02_11b tune, read at run time from the
//    model registry exactly as the calculators read them.
//  * Knots at -2, -1, 0, +1, +2 sigma with the GENIE Reweight default sigmas
//    (GSystUncertainty), converted to physical values: x = nominal*(1 + n*sigma),
//    which is how every calculator used here forms its tweaked value.
//    The GENIE Reweight README says these defaults are not the systematic error
//    of any tune; they are a starting point and may change later.
//  * 9 spline dials (5 evenly spaced knots, catmull-rom in GUNDAM).
//  * 21 normalisation dials. FakeSelect carries no GENIE interaction mode, so a
//    Norm dial with an applyCondition is not possible without forking
//    FakeSelect.C. Instead each affected event gets a graph through
//    (1-sigma_-, w), (1, 1), (1+sigma_+, w) with w = x exactly (GUNDAM "Graph",
//    linear, extrapolates linearly). GENIE's own calculator decides which
//    events are affected, so no event classification is re-implemented here.
//    Lower bound 0, upper bound 1 + 3 sigma.
//  * NormCCQE has an asymmetric sigma (+0.20/-0.15). The graph uses both sides
//    (it is exact either way); the Gaussian prior uses sigma = 0.20.
//  * Dropped dials: Bodek-Yang (nominal offset from plain QPM vs stored
//    KNOTuned dsigma), MaCOHpi/R0COHpi (t factor in stored COH dsigma), AGKY,
//    RDecBR, FracPN/FracDelta (final state only; FakeSelect cannot see them),
//    MaCCQEshape and the RES shape dials (alternative parametrisations),
//    NormCCRES/NormNCRES (would double count MaCCRES/MaNCRES).
//  * Per-event rules copied from src/MakeSplines.C:
//      - nominal knot written with weight exactly 1;
//      - all weights exactly 1 -> empty graph (no response, no GUNDAM dial);
//      - any weight NaN/inf/negative, or all weights 0 -> empty graph, counted
//        as rejected. This includes the negative MaCCQE weights of some nubar
//        CCQE events: they lose the MaCCQE response, as in the NEUT chain.
//      - TGraph branches with split level 0.
//  * Branch and GUNDAM parameter names are GENIE_<GSyst name>. NEUT uses the
//    bare names (MaCCQE, nominal 1.21 GeV); the prefix makes a NEUT config read
//    against GENIE data (or the reverse) fail on a missing leaf instead of
//    silently binding the wrong prior to a graph, and lets both coexist in one
//    tree, since FakeSelect carries TGraph branches by class, not by name.
//  * Events are matched by entry. GHEP and flat tree have equal entry counts;
//    FakeRecon's drop rules (Q2 < -999; CC without a lepton) are applied here
//    too, and every kept event is checked against its fakerecon row (Enu_true,
//    PDGnu, cc, bitwise). Any mismatch is fatal.
//
//  * (Stage A note) "root -q" exits with the macro's return value, and
//    FakeRecon.C (and FakeSelect.C) return 1 on success, so a "set -e" script
//    treats a successful run as a failure. Not fixed in the shared macros or
//    scripts; scripts/genie/run_all.sh ignores the exit status of
//    scripts/recon.sh and scripts/select.sh and checks their outputs instead.
//  * (Stage A note) The AhrensNCELPXSec.xml of the "neutlike" config variant
//    (not part of this repository) sets EL-Axial-Eta = 0, which would make
//    EtaNCEL a null dial (x = 0*(1 + n*sigma)). This chain sets no GXMLPATH or
//    --xml-path, so the stock GENIE config (0.12) is used, as in the dial scan. Guarded: the run stops unless the registry
//    gives EtaNCEL = 0.12 (expect value in kDials).
//
//  * (Stage C decision) CCQEPauliSupViaKF stays plain catmull-rom, like the
//    other spline dials. Its per-event response is a step: for
//    QELEventGenerator events GReWeightFGM returns exactly 0 or 1 (0 if the
//    final nucleon momentum is below the tweaked k_F), so a blocked event has
//    knots (1, 1, 1, 0, 0) and no 5-knot interpolation is exact between knots
//    (all are exact at the knots). Three interpolations were compared with
//    GENIE at -1.5, -0.5, +0.5, +1.5 sigma (numu_x_numu_fhc; largest rate
//    difference over the 3 FHC samples, % of the rate):
//        twk     catmull-rom   linear (Graph)   catmull-rom, monotonic
//        -1.5    2e-5          5e-5             5e-5
//        -0.5    0.177         1.2e-4           9e-6
//        +0.5    0.018         0.266            0.266
//        +1.5    0.068         0.211            0.346
//        max     0.177         0.266            0.346
//    Rule agreed beforehand: use monotonic only if it beats both on the
//    largest gap, otherwise keep plain catmull-rom. Catmull-rom has the
//    smallest largest gap; its error is an overshoot below nominal, where the
//    dial does almost nothing, while linear and monotonic underestimate the
//    curvature above nominal, where it acts.
//    Physical relevance uncertain, flagged for mentor review: Nieves + LFG
//    applies its own Pauli blocking in the cross section, while this dial
//    takes k_F from GENIE's generic PauliBlocker (GReWeightFGM calls this
//    "a bit of a hack").
//
//  Reweighting setup copied from validation/genie/dial_scan.cxx (itself
//  following grwght1p of GENIE Reweight R-1_04_02): same tune, spline file,
//  calculator list, mode switches and MEC Tl/ctl fix.
// =============================================================================

#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <TFile.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTree.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/XSecSplineList.h"

#include "RwFramework/GReWeight.h"
#include "RwFramework/GSyst.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"
#include "RwCalculators/GReWeightFGM.h"
#include "RwCalculators/GReWeightNonResonanceBkg.h"
#include "RwCalculators/GReWeightNuXSecCCQE.h"
#include "RwCalculators/GReWeightNuXSecCCRES.h"
#include "RwCalculators/GReWeightNuXSecCOH.h"
#include "RwCalculators/GReWeightNuXSecNCEL.h"
#include "RwCalculators/GReWeightNuXSecNCRES.h"
#include "RwCalculators/GReWeightXSecEmpiricalMEC.h"
#include "RwCalculators/GReWeightXSecMEC.h"

using namespace genie;
using namespace genie::rew;

// ----------------------------------------------------------------- inputs
const std::string kTune    = "G18_10b_02_11b";

// ------------------------------------------------------------------ dials
// model/key: where the calculator reads its nominal value (empty = the dial
// scales something by (1 + n*sigma), so the nominal is 1).
// expect: tune value from docs/genie/report.md; the run stops if the registry disagrees
// (0 = not checked, only logged).
struct DialDef {
  const char *name;
  bool spline;
  const char *model, *key;
  double expect;
};
const std::vector<DialDef> kDials = {
  // ---- spline dials: 5 knots at -2..+2 sigma, x in physical units
  {"MaCCQE",            true, "QEL-CC", "FormFactorsAlg/AxialFormFactorModel/QEL-Ma", 0.994989}, // CCQE dipole axial mass [GeV], sigma 3%
  {"CCQEPauliSupViaKF", true, "", "", 0},          // Fermi-momentum scale for CCQE Pauli suppression (one-sided response);
                                                   // step-like per event, flagged for mentor review (decision log)
  {"CoulombCCQE",       true, "", "", 0},          // scale of the Coulomb potential correction in CCQE
  {"MaCCRES",           true, "RES-CC", "RES-Ma", 1.088962}, // CC RES axial mass [GeV], sigma 20%; tune value, not model default 1.12
  {"MvCCRES",           true, "RES-CC", "RES-Mv", 0.840},    // CC RES vector mass [GeV], sigma 10%
  {"MaNCRES",           true, "RES-NC", "RES-Ma", 1.088962}, // NC RES axial mass [GeV], sigma 20%
  {"MvNCRES",           true, "RES-NC", "RES-Mv", 0.840},    // NC RES vector mass [GeV], sigma 10%
  {"MaNCEL",            true, "QEL-NC", "QEL-Ma", 0},        // NC elastic axial mass [GeV]
  {"EtaNCEL",           true, "QEL-NC", "EL-Axial-Eta", 0.12}, // NC elastic strange axial form factor eta; stock 0.12
                                                             // (the neutlike config sets 0; not used here)
  // ---- normalisation dials: graph y = x through -1, 0, +1 sigma; nominal 1
  {"NormCCQE",  false, "", "", 0},   // CCQE (not charm/strange), sigma +0.20/-0.15; prior sigma 0.20
  {"NormCCMEC", false, "", "", 0},   // CC 2p2h (Valencia), sigma 0.5
  {"NormNCMEC", false, "", "", 0},   // NC 2p2h, sigma 0.5
  {"NormCCCOH", false, "", "", 0},   // CC coherent pi (Berger-Sehgal), sigma 1.0; the -1 sigma point is 0
  {"NormNCCOH", false, "", "", 0},   // NC coherent pi, sigma 1.0
  // non-resonant background: DIS with W < 2 GeV, one nucleon + 1 or 2 pions,
  // weight R_curr/R_def = 1 + n*0.5 (independent of the tune's R_def)
  {"NonRESBGvpCC1pi",    false, "", "", 0},  // nu p CC 1pi
  {"NonRESBGvpCC2pi",    false, "", "", 0},  // nu p CC 2pi
  {"NonRESBGvpNC1pi",    false, "", "", 0},  // nu p NC 1pi
  {"NonRESBGvpNC2pi",    false, "", "", 0},  // nu p NC 2pi
  {"NonRESBGvnCC1pi",    false, "", "", 0},  // nu n CC 1pi
  {"NonRESBGvnCC2pi",    false, "", "", 0},  // nu n CC 2pi (largest in FHC 1Rmu)
  {"NonRESBGvnNC1pi",    false, "", "", 0},  // nu n NC 1pi
  {"NonRESBGvnNC2pi",    false, "", "", 0},  // nu n NC 2pi
  {"NonRESBGvbarpCC1pi", false, "", "", 0},  // nubar p CC 1pi
  {"NonRESBGvbarpCC2pi", false, "", "", 0},  // nubar p CC 2pi
  {"NonRESBGvbarpNC1pi", false, "", "", 0},  // nubar p NC 1pi
  {"NonRESBGvbarpNC2pi", false, "", "", 0},  // nubar p NC 2pi
  {"NonRESBGvbarnCC1pi", false, "", "", 0},  // nubar n CC 1pi
  {"NonRESBGvbarnCC2pi", false, "", "", 0},  // nubar n CC 2pi
  {"NonRESBGvbarnNC1pi", false, "", "", 0},  // nubar n NC 1pi
  {"NonRESBGvbarnNC2pi", false, "", "", 0},  // nubar n NC 2pi
};
const std::vector<double> kSplinePts = {-2, -1, 0, 1, 2};
const std::vector<double> kNormPts   = {-1, 0, 1};
// interp mode (a check, not production): GENIE's own weights between and
// beyond the inner knots, to test the spline interpolation in GUNDAM
const std::vector<double> kInterpPts = {-1.5, -0.5, 0.5, 1.5};
const double kNormCCQEPriorSigma = 0.20;   // decision: Gaussian prior uses the + side
const double kNormUpperSigmas    = 3.0;    // decision: norm upper bound 1 + 3 sigma

// ------------------------------------------------------------ GENIE setup
void InitGenie(const std::string &msg_xml, const std::string &xsec_xml)
{
  utils::app_init::MesgThresholds(msg_xml);
  utils::app_init::RandGen(12345);
  GHepRecord::SetPrintLevel(0);
  RunOpt::Instance()->SetTuneName(kTune);
  RunOpt::Instance()->BuildTune();
  XSecSplineList::Instance()->LoadFromXml(xsec_xml);
  Messenger::Instance()->SetPriorityLevel("ReW", pWARN);
}

// The calculators that can handle the dials above (dial_scan.cxx's list,
// restricted to these dials). Which of them handles a dial is asked at run time.
std::vector<std::pair<std::string, std::function<GReWeightI*()>>> Calculators()
{
  return {
    {"xsec_ncel",      [] { return new GReWeightNuXSecNCEL; }},
    {"xsec_ccqe",      [] { return new GReWeightNuXSecCCQE; }},
    {"xsec_ccres",     [] { return new GReWeightNuXSecCCRES; }},
    {"xsec_ncres",     [] { return new GReWeightNuXSecNCRES; }},
    {"xsec_nonresbkg", [] { return new GReWeightNonResonanceBkg; }},
    {"xsec_coh",       [] { return new GReWeightNuXSecCOH; }},
    {"nuclear_qe",     [] { return new GReWeightFGM; }},
    {"xsec_empmec",    [] { return new GReWeightXSecEmpiricalMEC; }},
    {"xsec_mec",       [] { return new GReWeightXSecMEC; }},
  };
}

// The mode switches of gRwght1Param.cxx, as in dial_scan.cxx.
void ApplyMode(GReWeightI *calc, GSyst_t s)
{
  if (s == kXSecTwkDial_MaCCQE)
    if (auto c = dynamic_cast<GReWeightNuXSecCCQE*>(calc)) c->SetMode(GReWeightNuXSecCCQE::kModeMa);
  if (s == kXSecTwkDial_MaCCRES || s == kXSecTwkDial_MvCCRES)
    if (auto c = dynamic_cast<GReWeightNuXSecCCRES*>(calc)) c->SetMode(GReWeightNuXSecCCRES::kModeMaMv);
  if (s == kXSecTwkDial_MaNCRES || s == kXSecTwkDial_MvNCRES)
    if (auto c = dynamic_cast<GReWeightNuXSecNCRES*>(calc)) c->SetMode(GReWeightNuXSecNCRES::kModeMaMv);
}

// The nominal a calculator starts from, read the way its Init() reads it.
double ModelParam(const std::string &model, const std::string &key)
{
  Registry *gpl = AlgConfigPool::Instance()->GlobalParameterList();
  AlgId id(gpl->GetAlg("XSecModel@genie::EventGenerator/" + model));
  Algorithm *alg = AlgFactory::Instance()->AdoptAlgorithm(id);
  alg->AdoptSubstructure();
  const double v = alg->GetConfig().GetDouble(key);
  delete alg;
  return v;
}

// ================================================================== main
struct Dial {
  DialDef def;
  GSyst_t syst;
  double nominal, sig_p, sig_m;          // sigmas are fractional, as in the table
  std::vector<double> twk, x;            // twk-dial value and physical value per knot
  int inom;                              // index of the nominal knot
  std::vector<GReWeight*> rw;            // one per knot (nullptr at the nominal)
  std::vector<GReWeight*> rw_interp;     // interp mode only: one per kInterpPts
  std::vector<double> w_interp;
  double prior_sigma, lower, upper;
  Long64_t n_resp = 0, n_flat = 0, n_bad = 0, n_xy = 0;
  TGraph *graph = nullptr;
};

int main(int argc, char **argv)
{
  const bool interp = argc > 1 && std::string(argv[1]) == "--interp";
  const int a = interp ? 2 : 1;
  if (argc - a < 6) {
    std::cerr << "usage: genie_splines [--interp] <ghep.root> <flattree.root> <fakerecon.root> "
                 "<xsec.xml> <messenger.xml> <out.root> [nmax]" << std::endl;
    return 1;
  }
  const std::string ghep_file = argv[a], flat_file = argv[a + 1], rec_file = argv[a + 2];
  const std::string xsec_xml = argv[a + 3], msg_xml = argv[a + 4], out_file = argv[a + 5];
  const Long64_t nmax = argc > a + 6 ? atoll(argv[a + 6]) : -1;

  InitGenie(msg_xml, xsec_xml);
  GSystUncertainty *unc = GSystUncertainty::Instance();

  // ---- dial setup: nominal, sigma, knots, one GReWeight per non-nominal knot
  auto calcs = Calculators();
  std::vector<Dial> dial;
  for (auto &def : kDials) {
    Dial d;
    d.def = def;
    d.syst = GSyst::FromString(def.name);
    if (d.syst == kNullSystematic) { std::cerr << "[ERROR] unknown dial " << def.name << std::endl; return 1; }
    d.nominal = (*def.model) ? ModelParam(def.model, def.key) : 1.0;
    if (def.expect > 0 && std::fabs(d.nominal / def.expect - 1) > 1e-5) {
      std::cerr << "[ERROR] " << def.name << " nominal " << d.nominal << " from the registry, expected "
                << def.expect << std::endl;
      return 1;
    }
    d.sig_p = unc->OneSigmaErr(d.syst, +1);
    d.sig_m = unc->OneSigmaErr(d.syst, -1);
    // catmull-rom needs evenly spaced knots, so spline dials must be symmetric
    if (def.spline && std::fabs(d.sig_p - d.sig_m) > 1e-12) {
      std::cerr << "[ERROR] " << def.name << " has asymmetric sigma " << d.sig_p << "/" << d.sig_m << std::endl;
      return 1;
    }
    d.twk = def.spline ? kSplinePts : kNormPts;
    for (size_t k = 0; k < d.twk.size(); ++k) {
      const double n = d.twk[k];
      d.x.push_back(d.nominal * (1 + n * (n < 0 ? d.sig_m : d.sig_p)));
      if (n == 0) d.inom = k;
    }
    if (def.spline) {
      d.prior_sigma = d.nominal * d.sig_p;
      d.lower = d.x.front();
      d.upper = d.x.back();
    } else {
      d.prior_sigma = (d.syst == kXSecTwkDial_NormCCQE) ? kNormCCQEPriorSigma : d.sig_p;
      d.lower = 0.0;
      d.upper = 1.0 + kNormUpperSigmas * d.prior_sigma;
    }

    // calculators that handle this dial (asked of fresh instances with the mode set)
    std::vector<int> handlers;
    for (size_t c = 0; c < calcs.size(); ++c) {
      GReWeightI *calc = calcs[c].second();
      ApplyMode(calc, d.syst);
      if (calc->IsHandled(d.syst)) handlers.push_back(c);
      delete calc;
    }
    if (handlers.empty()) { std::cerr << "[ERROR] no calculator handles " << def.name << std::endl; return 1; }

    auto make_rw = [&](double twk) {
      GReWeight *rw = new GReWeight;
      for (int c : handlers) {
        GReWeightI *calc = calcs[c].second();
        ApplyMode(calc, d.syst);
        rw->AdoptWghtCalc(calcs[c].first, calc);
      }
      rw->Systematics().Init(d.syst);
      rw->Systematics().Set(d.syst, twk);
      rw->Reconfigure();
      return rw;
    };
    for (size_t k = 0; k < d.twk.size(); ++k)
      d.rw.push_back((int)k == d.inom ? nullptr : make_rw(d.twk[k]));
    if (interp && def.spline)
      for (double t : kInterpPts) d.rw_interp.push_back(make_rw(t));

    printf("[INFO] GENIE_%-20s %-6s nominal %.6g  sigma +%.4g/-%.4g  %zu calculator(s)  knots:",
           def.name, def.spline ? "spline" : "norm", d.nominal, d.sig_p, d.sig_m, handlers.size());
    for (double x : d.x) printf(" %.6g", x);
    printf("  prior %.6g +- %.6g  bounds [%.6g, %.6g]\n", d.nominal, d.prior_sigma, d.lower, d.upper);
    dial.push_back(d);
  }

  // ---- inputs
  TFile fg(ghep_file.c_str());
  TTree *gtree = (TTree*)fg.Get("gtree");
  NtpMCEventRecord *mcrec = nullptr;
  gtree->SetBranchAddress("gmcrec", &mcrec);

  TFile fl(flat_file.c_str());
  TTree *flat = (TTree*)fl.Get("FlatTree_VARS");
  float Enu_true, ELep, Q2; int PDGnu, PDGLep; Char_t cc;
  flat->SetBranchAddress("Enu_true", &Enu_true);
  flat->SetBranchAddress("ELep", &ELep);
  flat->SetBranchAddress("Q2", &Q2);
  flat->SetBranchAddress("PDGnu", &PDGnu);
  flat->SetBranchAddress("PDGLep", &PDGLep);
  flat->SetBranchAddress("cc", &cc);

  TFile fr(rec_file.c_str());
  TTree *rec = (TTree*)fr.Get("FlatTree_VARS");
  if (!gtree || !flat || !rec) { std::cerr << "[ERROR] missing input tree" << std::endl; return 1; }
  if (flat->GetEntries() != gtree->GetEntries()) {
    std::cerr << "[ERROR] flat/GHEP entry mismatch" << std::endl; return 1; }
  float r_Enu; int r_PDGnu; Char_t r_cc;
  rec->SetBranchAddress("Enu_true", &r_Enu);
  rec->SetBranchAddress("PDGnu", &r_PDGnu);
  rec->SetBranchAddress("cc", &r_cc);
  for (auto &d : dial)
    if (rec->GetBranch(("GENIE_" + std::string(d.def.name)).c_str())) {
      std::cerr << "[ERROR] GENIE_" << d.def.name << " already exists in the fakerecon tree" << std::endl; return 1; }

  Long64_t nev = gtree->GetEntries();
  if (nmax > 0 && nmax < nev) nev = nmax;

  // ---- output: clone of the fakerecon tree + one TGraph per dial (split level 0).
  // In interp mode the output (written to a file of its own by the check
  // scripts, so production files are not touched) also has a GenieInterp tree:
  // per spline dial, GENIE's weights at kInterpPts, row by row with FlatTree_VARS.
  TFile fout(out_file.c_str(), "RECREATE");
  fout.cd();
  TTree *out = rec->CloneTree(0);
  for (auto &d : dial) {
    d.graph = new TGraph();
    out->Branch(("GENIE_" + std::string(d.def.name)).c_str(), "TGraph", &d.graph, 32000, 0);
  }
  TTree *tinterp = nullptr;
  if (interp) {
    tinterp = new TTree("GenieInterp", "GENIE weights at -1.5, -0.5, +0.5, +1.5 sigma");
    for (auto &d : dial)
      if (d.def.spline) tinterp->Branch(("GENIE_" + std::string(d.def.name)).c_str(), &d.w_interp);
  }

  // ---- event loop
  TStopwatch clock;
  Long64_t j = 0, n_dropped = 0;
  std::vector<double> w;
  for (Long64_t i = 0; i < nev; ++i) {
    gtree->GetEntry(i);
    flat->GetEntry(i);
    EventRecord &ev = *(mcrec->event);
    if (std::fabs(ev.Probe()->E() - Enu_true) > 1e-4 * std::max(1.f, Enu_true) ||
        ev.Probe()->Pdg() != PDGnu || (bool)cc != ev.Summary()->ProcInfo().IsWeakCC()) {
      std::cerr << "[ERROR] GHEP/flat mismatch at entry " << i << std::endl; return 1; }

    // FakeRecon.C drops: Pauli blocked (Q2 = -999) and CC without a lepton.
    if (Q2 < -999.0 || (cc && (PDGLep == 0 || ELep < 0.0))) { ++n_dropped; mcrec->Clear(); continue; }

    if (j >= rec->GetEntries()) { std::cerr << "[ERROR] more kept events than fakerecon entries" << std::endl; return 1; }
    rec->GetEntry(j);
    if (r_Enu != Enu_true || r_PDGnu != PDGnu || r_cc != cc) {
      std::cerr << "[ERROR] fakerecon row " << j << " does not match GHEP entry " << i << std::endl; return 1; }

    // As grwght1p: set Tl and ctl for MEC reweighting.
    Interaction *in = ev.Summary();
    const TLorentzVector &p4l = in->KinePtr()->FSLeptonP4();
    in->KinePtr()->SetKV(kKVTl, p4l.E() - in->FSPrimLepton()->Mass());
    in->KinePtr()->SetKV(kKVctl, p4l.CosTheta());

    for (auto &d : dial) {
      const size_t np = d.x.size();
      w.assign(np, 1.0);
      bool all_one = true, all_zero = true, bad = false;
      for (size_t k = 0; k < np; ++k) {
        if ((int)k == d.inom) continue;              // weight 1 by construction
        w[k] = d.rw[k]->CalcWeight(ev);
        if (w[k] != 1.0) all_one = false;
        if (w[k] != 0.0) all_zero = false;
        if (!std::isfinite(w[k]) || w[k] < 0.0) bad = true;
      }
      if (all_zero) bad = true;
      if (bad)          { d.graph->Set(0); ++d.n_bad;  continue; }
      if (all_one)      { d.graph->Set(0); ++d.n_flat; continue; }
      d.graph->Set(np);
      for (size_t k = 0; k < np; ++k) d.graph->SetPoint(k, d.x[k], w[k]);
      ++d.n_resp;
      // a norm graph must be the line y = x
      if (!d.def.spline)
        for (size_t k = 0; k < np; ++k) if (std::fabs(w[k] - d.x[k]) > 1e-9) { ++d.n_xy; break; }
    }
    out->Fill();
    if (tinterp) {
      for (auto &d : dial) {
        d.w_interp.clear();
        for (auto *rw : d.rw_interp) d.w_interp.push_back(rw->CalcWeight(ev));
      }
      tinterp->Fill();
    }
    ++j;
    mcrec->Clear();
    if (i % 20000 == 0) {
      std::cout << "[INFO] entry " << i << " / " << nev << ", " << clock.RealTime() << " s" << std::endl;
      clock.Continue();
    }
  }
  if (nmax < 0 && j != rec->GetEntries()) {
    std::cerr << "[ERROR] " << j << " kept events but " << rec->GetEntries() << " fakerecon entries" << std::endl; return 1; }

  // ---- metadata (GENIE's counterpart of TweakMetaData), read by MakeGenieXsecCov.C
  TTree meta("GenieDialMeta", "GENIE dial definitions");
  std::string m_name, m_kind;
  double m_nom, m_sp, m_sm, m_prior_sigma, m_lo, m_hi;
  std::vector<double> m_x, m_twk;
  Long64_t m_resp, m_flat, m_bad, m_xy;
  meta.Branch("name", &m_name);
  meta.Branch("kind", &m_kind);
  meta.Branch("nominal", &m_nom);
  meta.Branch("sigma_plus", &m_sp);
  meta.Branch("sigma_minus", &m_sm);
  meta.Branch("prior_sigma", &m_prior_sigma);
  meta.Branch("lower", &m_lo);
  meta.Branch("upper", &m_hi);
  meta.Branch("twk", &m_twk);
  meta.Branch("knots", &m_x);
  meta.Branch("n_resp", &m_resp);
  meta.Branch("n_flat", &m_flat);
  meta.Branch("n_bad", &m_bad);
  meta.Branch("n_xy", &m_xy);
  for (auto &d : dial) {
    m_name = "GENIE_" + std::string(d.def.name); m_kind = d.def.spline ? "spline" : "norm";
    m_nom = d.nominal; m_sp = d.sig_p; m_sm = d.sig_m; m_prior_sigma = d.prior_sigma;
    m_lo = d.lower; m_hi = d.upper; m_twk = d.twk; m_x = d.x;
    m_resp = d.n_resp; m_flat = d.n_flat; m_bad = d.n_bad; m_xy = d.n_xy;
    meta.Fill();
  }
  fout.cd();
  out->Write("", TObject::kOverwrite);
  meta.Write("", TObject::kOverwrite);
  if (tinterp) tinterp->Write("", TObject::kOverwrite);
  fout.Close();

  // ---- summary
  std::cout << "[INFO] " << ghep_file << ": " << nev << " GHEP entries, " << n_dropped << " dropped as in FakeRecon, "
            << j << " written to " << out_file << std::endl;
  printf("[INFO] %-26s %10s %10s %10s %8s\n", "dial", "response", "no resp.", "rejected", "y!=x");
  int fail = 0;
  for (auto &d : dial) {
    printf("[INFO] GENIE_%-20s %10lld %10lld %10lld %8lld%s\n", d.def.name, d.n_resp, d.n_flat, d.n_bad, d.n_xy,
           d.n_bad ? "  <-- rejected" : "");
    if (d.n_xy) ++fail;
  }
  if (fail) { std::cerr << "[ERROR] " << fail << " norm dial(s) with y != x" << std::endl; return 1; }
  return 0;
}

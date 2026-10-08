// =============================================================================
//  check_genie_splines.cxx  -  file-level checks of one <stem>_spline.root.
//
//  Usage:
//      check_genie_splines <sample> <spline.root> <fakerecon.root> <scan part_<sample>.root> <response.root>
//      check_genie_splines select <fakeselect.root> <fakeselect_nodial.root> <spline files in FakeSelect order>
//  (the second mode is described above CheckSelect below)
//
//  Checks of the first mode (decisions in genie_splines.cxx):
//   1. The spline tree has the fakerecon entries, row by row (Enu_true, PDGnu,
//      cc, POTWeight, PLep_rec equal bitwise).
//   2. Every GENIE_ branch is a TGraph with split level 0.
//   3. Every non-empty graph has the knots of GenieDialMeta, evenly spaced for
//      spline dials, y = 1 exactly at the nominal knot, no NaN or negative y,
//      and y = x for norm dials.
//   4. Rates: the per-sample rates obtained by weighting each event with its
//      graph value at a knot reproduce the dial scan (part_<sample>.root) at
//      the same tweak. Rates use the scan's expectation weighting:
//      POTWeight * OscWeight * SelEffWeight * mix(RespCat, Enu_true), i.e. the
//      FakeSelect probabilities instead of its random draw. Events whose graph
//      was rejected count with weight 1 here but with their raw weight in the
//      scan, so rejected events show up as a difference.
//   5. Events with a response vs the scan's count of w != 1 at +1 sigma.
// =============================================================================

#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <TChain.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TKey.h>
#include <TTree.h>

// ---------------------------------------------------------------------------
//  Mode "select": check_genie_splines select <fakeselect.root> <fakeselect_nodial.root> <spline1> ... <spline6>
//   a. Same entries and identical sample assignment and kinematics as a
//      FakeSelect of the same events without dials (same seed): the graphs
//      must not change the draw.
//   b. Every GENIE_ graph arrives as a TGraph and equals the one in its
//      _spline.root input (FakeSelect chains its inputs in the order given,
//      which must be FakeSelect.C's file order).
int CheckSelect(int argc, char **argv)
{
  int nfail = 0;
  TFile fa(argv[2]), fb(argv[3]);
  TTree *a = (TTree*)fa.Get("FakeSelect"), *b = (TTree*)fb.Get("FakeSelect");
  if (!a || !b) { std::cerr << "[ERROR] missing FakeSelect tree" << std::endl; return 1; }
  TChain in("FlatTree_VARS");
  for (int i = 4; i < argc; ++i) in.Add(argv[i]);
  const Long64_t n = a->GetEntries();
  std::cout << "[INFO] entries: with dials " << n << ", without " << b->GetEntries()
            << ", spline inputs " << in.GetEntries() << std::endl;
  if (b->GetEntries() != n || in.GetEntries() != n) { std::cout << "[FAIL] entry counts differ" << std::endl; return 1; }

  int ta, tb, sa, sb; float ea, eb; double wa, wb;
  a->SetBranchAddress("EventType", &ta);   b->SetBranchAddress("EventType", &tb);
  a->SetBranchAddress("SampleIndex", &sa); b->SetBranchAddress("SampleIndex", &sb);
  a->SetBranchAddress("Enu_rec", &ea);     b->SetBranchAddress("Enu_rec", &eb);
  a->SetBranchAddress("SelEffWeight", &wa); b->SetBranchAddress("SelEffWeight", &wb);

  std::vector<std::string> names;
  for (auto *o : *a->GetListOfBranches()) {
    TBranch *br = (TBranch*)o;
    if (std::string(br->GetName()).rfind("GENIE_", 0) != 0) continue;
    TClass *cl = TClass::GetClass(br->GetClassName());
    if (!cl || !cl->InheritsFrom(TGraph::Class()) || br->GetSplitLevel() != 0) {
      std::cout << "[FAIL] " << br->GetName() << " is not a split-level-0 TGraph" << std::endl; ++nfail; }
    names.push_back(br->GetName());
  }
  std::cout << "[INFO] " << names.size() << " GENIE_ graph branches in the FakeSelect tree" << std::endl;
  std::vector<TGraph*> ga(names.size(), nullptr), gi(names.size(), nullptr);
  for (size_t d = 0; d < names.size(); ++d) {
    a->SetBranchAddress(names[d].c_str(), &ga[d]);
    in.SetBranchAddress(names[d].c_str(), &gi[d]);
  }

  Long64_t bad_sel = 0, bad_graph = 0;
  std::vector<Long64_t> nresp(names.size(), 0);
  std::map<int, Long64_t> per_type;
  for (Long64_t i = 0; i < n; ++i) {
    a->GetEntry(i); b->GetEntry(i); in.GetEntry(i);
    if (ta != tb || sa != sb || ea != eb || wa != wb) ++bad_sel;
    ++per_type[ta];
    for (size_t d = 0; d < names.size(); ++d) {
      const int np = ga[d]->GetN();
      if (np > 0 && ta >= 0) ++nresp[d];
      bool same = np == gi[d]->GetN();
      for (int k = 0; same && k < np; ++k)
        same = ga[d]->GetX()[k] == gi[d]->GetX()[k] && ga[d]->GetY()[k] == gi[d]->GetY()[k];
      if (!same) ++bad_graph;
    }
  }
  printf("[INFO] entries whose EventType/SampleIndex/Enu_rec/SelEffWeight differ from the no-dial FakeSelect: %lld\n", bad_sel);
  printf("[INFO] graphs differing from their _spline.root input: %lld\n", bad_graph);
  if (bad_sel) ++nfail;
  if (bad_graph) ++nfail;
  printf("[INFO] entries per EventType:");
  for (auto &p : per_type) printf("  %d: %lld", p.first, p.second);
  printf("\n[INFO] selected events (EventType >= 0) with a response, per dial:\n");
  for (size_t d = 0; d < names.size(); ++d) printf("  %-26s %8lld\n", names[d].c_str(), nresp[d]);
  std::cout << "\n" << (nfail ? "[RESULT] FAILED " + std::to_string(nfail) + " check(s)" : "[RESULT] all checks passed")
            << std::endl;
  return nfail ? 1 : 0;
}

// ---------------------------------------------------------------------------
//  Mode "interp": check_genie_splines interp <sample> <interp file> <response.root> [linear dials] [monotonic dials]
//  ([linear dials]: comma-separated GENIE_ names that GUNDAM reads as dialType
//  Graph, evaluated linearly here; optional 6th argument, likewise, for dials
//  read as "catmull-rom, monotonic"; all others with plain catmull-rom. Pass ""
//  for an empty list.)
//  The interp file is written by "genie_splines ... <nmax> <interp_file>": the
//  usual FlatTree_VARS with graphs, plus GenieInterp holding GENIE's own
//  weights at -1.5, -0.5, +0.5, +1.5 sigma for the 9 spline dials. Here each
//  graph is evaluated there with GUNDAM's catmull-rom (CompactSpline, copied
//  from gundam Utils/include/CalculateCompactSpline.h, lower bound 0) and
//  compared with GENIE, per event and as per-sample rate changes (expectation
//  weighting, as in the first mode). Tests the interpolation between knots and
//  half a sigma beyond the inner knots, which the knot-level checks cannot.
double CompactSpline(double x, const TGraph *g, bool monotonic = false)
{
  const int dim = g->GetN();
  const double *y = g->GetY();
  const double low = g->GetX()[0], step = g->GetX()[1] - g->GetX()[0];
  const double xx = (x - low) / step;
  const int ix = (xx < 0) ? xx - 1 : xx;
  auto clampi = [&](int i) { return std::min(std::max(i, 0), dim - 2); };
  const int d21_0 = clampi(ix - 1), d32_0 = clampi(ix), d43_0 = clampi(ix + 1);
  const double p2 = y[d32_0], p3 = y[d32_0 + 1], fx = xx - d32_0;
  const double d21 = y[d21_0 + 1] - y[d21_0], d32 = p3 - p2, d43 = y[d43_0 + 1] - y[d43_0];
  double m2 = 0.5 * (d21 + d32), m3 = 0.5 * (d32 + d43);
  if (monotonic) {   // "catmull-rom, monotonic": Fritsch-Carlson limits, as CalculateMonotonicSpline.h
    if (d32 * d21 <= 0.0) m2 = 0.0;
    if (d43 * d32 <= 0.0) m3 = 0.0;
    const double delta2 = 3.0 * std::min(std::fabs(d21), std::fabs(d32));
    const double delta3 = 3.0 * std::min(std::fabs(d32), std::fabs(d43));
    m2 = std::min(std::max(m2, -delta2), delta2);
    m3 = std::min(std::max(m3, -delta3), delta3);
  }
  double v = (((2.0 * p2 - 2.0 * p3 + m3 + m2) * fx + 3.0 * p3 - 3.0 * p2 - m3 - 2.0 * m2) * fx + m2) * fx + p2;
  return std::max(v, 0.0);
}

// GUNDAM's Graph dial (Utils/include/CalculateGraph.h): linear between knots,
// linear extrapolation from the end segments; response bounded below by 0.
double LinearGraph(double x, const TGraph *g)
{
  const int n = g->GetN();
  int k = 0;
  while (k < n - 2 && x > g->GetX()[k + 1]) ++k;
  const double x1 = g->GetX()[k], x2 = g->GetX()[k + 1], y1 = g->GetY()[k], y2 = g->GetY()[k + 1];
  return std::max(y1 + (y2 - y1) * (x - x1) / (x2 - x1), 0.0);
}

int CheckInterp(int argc, char **argv)
{
  const std::string sample = argv[2];
  const bool isRHC = sample.substr(sample.size() - 3) == "rhc";
  // optional 5th argument: comma-separated dials read as dialType Graph in GUNDAM
  const std::string linear_list = argc > 5 ? std::string(",") + argv[5] + "," : "";
  // optional 6th argument: comma-separated dials read as "catmull-rom, monotonic"
  const std::string mono_list = argc > 6 ? std::string(",") + argv[6] + "," : "";
  TFile f(argv[3]), fresp(argv[4]);
  TTree *t = (TTree*)f.Get("FlatTree_VARS"), *ti = (TTree*)f.Get("GenieInterp"), *meta = (TTree*)f.Get("GenieDialMeta");
  if (!t || !ti || !meta) { std::cerr << "[ERROR] missing tree" << std::endl; return 1; }
  const double pts[4] = {-1.5, -0.5, 0.5, 1.5};

  // spline dials and their knot spacing
  struct D { std::string name; double nominal, sigma; TGraph *g = nullptr; std::vector<double> *wi = nullptr; };
  std::vector<D> dial;
  {
    std::string *name = nullptr, *kind = nullptr; double nom, sp;
    meta->SetBranchAddress("name", &name); meta->SetBranchAddress("kind", &kind);
    meta->SetBranchAddress("nominal", &nom); meta->SetBranchAddress("sigma_plus", &sp);
    for (Long64_t i = 0; i < meta->GetEntries(); ++i) {
      meta->GetEntry(i);
      if (*kind == "spline") dial.push_back({*name, nom, sp});
    }
  }
  for (auto &d : dial) { t->SetBranchAddress(d.name.c_str(), &d.g); ti->SetBranchAddress(d.name.c_str(), &d.wi); }
  double potw, osc, eff; int rc; float enu;
  t->SetBranchAddress("POTWeight", &potw); t->SetBranchAddress("OscWeight", &osc);
  t->SetBranchAddress("SelEffWeight", &eff); t->SetBranchAddress("RespCat", &rc); t->SetBranchAddress("Enu_true", &enu);

  const char *cat[9] = {"numuCCQE", "numubarCCQE", "numuCCnonQE", "numubarCCnonQE",
                        "nueCCQE", "nuebarCCQE", "nueCCnonQE", "nuebarCCnonQE", "NC"};
  const char *mixn[3] = {"1Rmu", "1Re0de", "1Re1de"};
  TH1D *hmix[9][3];
  for (int c = 0; c < 9; ++c) for (int s = 0; s < 3; ++s)
    hmix[c][s] = (TH1D*)fresp.Get(Form("h_mix_%s_%s", cat[c], mixn[s]));
  auto lookup = [](TH1D *h, double x) { return h->GetBinContent(std::min(std::max(h->FindBin(x), 1), h->GetNbinsX())); };
  const int nsel = isRHC ? 2 : 3;
  const char *selname[5] = {"FHC_1Rmu", "FHC_1Re0de", "FHC_1Re1de", "RHC_1Rmu", "RHC_1Re0de"};
  const int sel0 = isRHC ? 3 : 0;

  const size_t nd = dial.size();
  std::vector<double> r0(nsel, 0);
  std::vector<std::vector<std::vector<double>>> rs(nd, std::vector<std::vector<double>>(4, std::vector<double>(nsel, 0))), rg = rs;
  std::vector<std::vector<double>> maxdev(nd, std::vector<double>(4, 0)), sumdev = maxdev;
  std::vector<Long64_t> naff(nd, 0), nempty_resp(nd, 0);
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i); ti->GetEntry(i);
    double sw[3];
    for (int s = 0; s < nsel; ++s) { sw[s] = potw * osc * eff * lookup(hmix[rc][s], enu); r0[s] += sw[s]; }
    for (size_t d = 0; d < nd; ++d) {
      const bool has = dial[d].g->GetN() > 0;
      if (has) ++naff[d];
      bool empty_resp = false;
      for (int k = 0; k < 4; ++k) {
        const double x = dial[d].nominal * (1 + pts[k] * dial[d].sigma);
        const bool lin = linear_list.find("," + dial[d].name + ",") != std::string::npos;
        const bool mono = mono_list.find("," + dial[d].name + ",") != std::string::npos;
        const double ys = has ? (lin ? LinearGraph(x, dial[d].g) : CompactSpline(x, dial[d].g, mono)) : 1.0;
        const double yg = (*dial[d].wi)[k];
        if (!has && yg != 1.0) empty_resp = true;
        if (has) { maxdev[d][k] = std::max(maxdev[d][k], std::fabs(ys - yg)); sumdev[d][k] += std::fabs(ys - yg); }
        for (int s = 0; s < nsel; ++s) { rs[d][k][s] += sw[s] * ys; rg[d][k][s] += sw[s] * yg; }
      }
      if (empty_resp) ++nempty_resp[d];
    }
  }

  printf("[INFO] %s: %lld events. Rate change (%%) per sample, GUNDAM interpolation (catmull-rom; linear or monotonic for the listed dials) / GENIE direct\n",
         sample.c_str(), t->GetEntries());
  printf("%-24s %5s", "dial", "twk");
  for (int s = 0; s < nsel; ++s) printf(" %22s", selname[sel0 + s]);
  printf("  %10s %12s %12s\n", "max|drate|", "max|dw| evt", "mean|dw| evt");
  for (size_t d = 0; d < nd; ++d) {
    for (int k = 0; k < 4; ++k) {
      printf("%-24s %+5.1f", dial[d].name.c_str(), pts[k]);
      double worst = 0;
      for (int s = 0; s < nsel; ++s) {
        const double a = 100 * (rs[d][k][s] / r0[s] - 1), b = 100 * (rg[d][k][s] / r0[s] - 1);
        printf("   %9.4f / %9.4f", a, b);
        worst = std::max(worst, std::fabs(a - b));
      }
      printf("  %10.2e %12.2e %12.2e\n", worst, maxdev[d][k], naff[d] ? sumdev[d][k] / naff[d] : 0.0);
    }
  }
  printf("[INFO] events with an empty graph but a GENIE weight != 1 at an interp point (rejected events):");
  for (size_t d = 0; d < nd; ++d) if (nempty_resp[d]) printf("  %s %lld", dial[d].name.c_str(), nempty_resp[d]);
  printf("\n");
  return 0;
}

int main(int argc, char **argv)
{
  if (argc >= 5 && std::string(argv[1]) == "select") return CheckSelect(argc, argv);
  if (argc >= 5 && std::string(argv[1]) == "interp") return CheckInterp(argc, argv);
  if (argc < 6) { std::cerr << "usage: check_genie_splines <sample> <spline> <fakerecon> <scan part> <response>" << std::endl; return 1; }
  const std::string sample = argv[1];
  const bool isRHC = sample.substr(sample.size() - 3) == "rhc";
  int nfail = 0;
  auto fail = [&](const std::string &msg) { std::cout << "[FAIL] " << msg << std::endl; ++nfail; };

  TFile fs(argv[2]), frec(argv[3]), fscan(argv[4]), fresp(argv[5]);
  TTree *t = (TTree*)fs.Get("FlatTree_VARS");
  TTree *meta = (TTree*)fs.Get("GenieDialMeta");
  TTree *rec = (TTree*)frec.Get("FlatTree_VARS");
  if (!t || !meta || !rec) { std::cerr << "[ERROR] missing tree" << std::endl; return 1; }

  // ---- dial definitions
  struct D { std::string name, kind; double nominal; std::vector<double> twk, x; long long bad;
             TGraph *g = nullptr; long long n_resp = 0; };
  std::vector<D> dial;
  {
    std::string *name = nullptr, *kind = nullptr; double nom; std::vector<double> *twk = nullptr, *x = nullptr;
    Long64_t bad;
    meta->SetBranchAddress("name", &name); meta->SetBranchAddress("kind", &kind);
    meta->SetBranchAddress("nominal", &nom); meta->SetBranchAddress("twk", &twk);
    meta->SetBranchAddress("knots", &x); meta->SetBranchAddress("n_bad", &bad);
    for (Long64_t i = 0; i < meta->GetEntries(); ++i) {
      meta->GetEntry(i);
      dial.push_back({*name, *kind, nom, *twk, *x, bad});
    }
  }
  std::cout << "[INFO] " << sample << ": " << dial.size() << " dials, " << t->GetEntries() << " entries" << std::endl;

  // ---- 1. row-by-row match with fakerecon
  if (t->GetEntries() != rec->GetEntries()) fail("entry count differs from fakerecon");
  float e1, e2, p1, p2; int pdg1, pdg2; Char_t cc1, cc2; double w1, w2;
  t->SetBranchAddress("Enu_true", &e1);   rec->SetBranchAddress("Enu_true", &e2);
  t->SetBranchAddress("PLep_rec", &p1);   rec->SetBranchAddress("PLep_rec", &p2);
  t->SetBranchAddress("PDGnu", &pdg1);    rec->SetBranchAddress("PDGnu", &pdg2);
  t->SetBranchAddress("cc", &cc1);        rec->SetBranchAddress("cc", &cc2);
  t->SetBranchAddress("POTWeight", &w1);  rec->SetBranchAddress("POTWeight", &w2);
  double osc, eff; int rc;
  t->SetBranchAddress("OscWeight", &osc);
  t->SetBranchAddress("SelEffWeight", &eff);
  t->SetBranchAddress("RespCat", &rc);

  // ---- 2. branch type and split level
  for (auto &d : dial) {
    TBranch *b = t->GetBranch(d.name.c_str());
    if (!b) { fail(d.name + ": no branch"); continue; }
    TClass *cl = TClass::GetClass(b->GetClassName());
    if (!cl || !cl->InheritsFrom(TGraph::Class())) fail(d.name + ": not a TGraph");
    if (b->GetSplitLevel() != 0) fail(d.name + ": split level " + std::to_string(b->GetSplitLevel()));
    if (d.kind == "spline")
      for (size_t k = 2; k < d.x.size(); ++k)
        if (std::fabs((d.x[k] - d.x[k-1]) - (d.x[1] - d.x[0])) > 1e-9 * std::fabs(d.x[1] - d.x[0]))
          fail(d.name + ": knots not evenly spaced");
    t->SetBranchAddress(d.name.c_str(), &d.g);
  }

  // ---- selection curves, as FakeSelect.C / dial_scan.cxx
  const char *cat[9] = {"numuCCQE", "numubarCCQE", "numuCCnonQE", "numubarCCnonQE",
                        "nueCCQE", "nuebarCCQE", "nueCCnonQE", "nuebarCCnonQE", "NC"};
  const char *mixn[3] = {"1Rmu", "1Re0de", "1Re1de"};
  TH1D *hmix[9][3];
  for (int c = 0; c < 9; ++c) for (int s = 0; s < 3; ++s)
    hmix[c][s] = (TH1D*)fresp.Get(Form("h_mix_%s_%s", cat[c], mixn[s]));
  auto lookup = [](TH1D *h, double x) {
    return h->GetBinContent(std::min(std::max(h->FindBin(x), 1), h->GetNbinsX())); };
  const int nsel = isRHC ? 2 : 3;
  const char *selname[5] = {"FHC_1Rmu", "FHC_1Re0de", "FHC_1Re1de", "RHC_1Rmu", "RHC_1Re0de"};
  const int sel0 = isRHC ? 3 : 0;

  // rate[d][k][s]
  std::vector<std::vector<std::vector<double>>> rate(dial.size());
  for (size_t d = 0; d < dial.size(); ++d)
    rate[d].assign(dial[d].x.size(), std::vector<double>(nsel, 0.0));
  std::vector<double> rnom(nsel, 0.0);

  // ---- 3. per-event graph checks, and the rate sums
  long long bad_x = 0, bad_nom = 0, bad_y = 0, bad_xy = 0, bad_n = 0;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    rec->GetEntry(i);
    if (e1 != e2 || p1 != p2 || pdg1 != pdg2 || cc1 != cc2 || w1 != w2) {
      fail("row " + std::to_string(i) + " differs from fakerecon"); break; }
    double sw[3];
    for (int s = 0; s < nsel; ++s) {
      sw[s] = w1 * osc * eff * lookup(hmix[rc][s], e1);
      rnom[s] += sw[s];
    }
    for (size_t d = 0; d < dial.size(); ++d) {
      D &dd = dial[d];
      const int n = dd.g->GetN();
      if (n > 0) {
        ++dd.n_resp;
        if (n != (int)dd.x.size()) { ++bad_n; continue; }
        for (int k = 0; k < n; ++k) {
          const double x = dd.g->GetX()[k], y = dd.g->GetY()[k];
          if (x != dd.x[k]) ++bad_x;
          if (dd.twk[k] == 0 && y != 1.0) ++bad_nom;
          if (!std::isfinite(y) || y < 0) ++bad_y;
          if (dd.kind == "norm" && std::fabs(y - x) > 1e-9) ++bad_xy;
        }
      }
      for (size_t k = 0; k < dd.x.size(); ++k) {
        const double y = n > 0 ? dd.g->GetY()[k] : 1.0;
        for (int s = 0; s < nsel; ++s) rate[d][k][s] += sw[s] * y;
      }
    }
  }
  printf("[INFO] graph checks: wrong #points %lld, x != knot %lld, y != 1 at nominal %lld, "
         "NaN/negative y %lld, norm y != x %lld\n", bad_n, bad_x, bad_nom, bad_y, bad_xy);
  if (bad_n || bad_x || bad_nom || bad_y || bad_xy) fail("graph checks");

  // ---- 4./5. comparison with the dial scan
  auto integral = [&](const std::string &path) {
    TH1D *h = (TH1D*)fscan.Get(path.c_str());
    return h ? h->Integral(0, h->GetNbinsX() + 1) : std::nan("");
  };
  std::vector<double> snom(nsel);
  for (int s = 0; s < nsel; ++s) snom[s] = integral(std::string("nominal/s_nom_") + selname[sel0 + s]);
  printf("\n[INFO] nominal rates (this file vs scan):");
  for (int s = 0; s < nsel; ++s) printf("  %s %.3f / %.3f", selname[sel0 + s], rnom[s], snom[s]);
  printf("\n");
  for (int s = 0; s < nsel; ++s)
    if (std::fabs(rnom[s] / snom[s] - 1) > 1e-6) fail(std::string("nominal rate ") + selname[sel0 + s]);

  // scan w != 1 counts at +1 sigma
  std::map<std::string, double> neq1;
  if (TTree *ws = (TTree*)fscan.Get("wstats")) {
    char dn[64]; double twk, n;
    ws->SetBranchAddress("dial", dn); ws->SetBranchAddress("twk", &twk); ws->SetBranchAddress("nneq1", &n);
    for (Long64_t i = 0; i < ws->GetEntries(); ++i) { ws->GetEntry(i); if (twk == 1) neq1[dn] = n; }
  }

  const char *ptname[5] = {"m2", "m1", "", "p1", "p2"};
  printf("\n[INFO] rate change (%%) per sample: this file / scan, and the largest |difference| (in %% of the rate)\n");
  printf("%-26s %5s", "dial", "twk");
  for (int s = 0; s < nsel; ++s) printf(" %22s", selname[sel0 + s]);
  printf(" %10s\n", "max|diff|");
  double worst_all = 0;
  for (size_t d = 0; d < dial.size(); ++d) {
    D &dd = dial[d];
    const std::string bare = dd.name.substr(6);
    for (size_t k = 0; k < dd.x.size(); ++k) {
      if (dd.twk[k] == 0) continue;
      const int ip = (int)std::lround(dd.twk[k]) + 2;
      printf("%-26s %+5.0f", dd.name.c_str(), dd.twk[k]);
      double worst = 0;
      for (int s = 0; s < nsel; ++s) {
        const double mine = 100 * (rate[d][k][s] / rnom[s] - 1);
        const double scan = 100 * (integral(bare + "/s_" + ptname[ip] + "_" + selname[sel0 + s]) / snom[s] - 1);
        printf("   %9.4f / %9.4f", mine, scan);
        worst = std::max(worst, std::fabs(mine - scan));
      }
      printf(" %10.2e%s\n", worst, (worst > 1e-4 && dd.bad == 0) ? "  <--" : "");
      if (dd.bad == 0) worst_all = std::max(worst_all, worst);
      if (worst > 1e-4 && dd.bad == 0) fail(dd.name + " rate differs from the scan");
    }
  }
  printf("[INFO] largest |difference| over dials with no rejected events: %.2e %% of the rate\n", worst_all);

  printf("\n[INFO] events with a response / rejected / scan's w!=1 at +1 sigma\n");
  for (auto &dd : dial) {
    const std::string bare = dd.name.substr(6);
    printf("  %-26s %8lld %6lld %8.0f\n", dd.name.c_str(), dd.n_resp, dd.bad, neq1.count(bare) ? neq1[bare] : -1.0);
  }

  std::cout << "\n" << (nfail ? "[RESULT] FAILED " + std::to_string(nfail) + " check(s)" : "[RESULT] all checks passed")
            << std::endl;
  return nfail ? 1 : 0;
}

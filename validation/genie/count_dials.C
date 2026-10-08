// =============================================================================
//  count_dials.C  -  independent recount of the GUNDAM dial counts.
//
//  Usage (validation/genie/run_checks.sh gundamcheck runs it):
//      root -l -b -q 'count_dials.C+("<fakeselect files, space separated>", "<GUNDAM config>/inputs/samples/binning")'
//
//  For each of the 5 samples (EventType 0..4) it counts the events GUNDAM keeps
//  (selection cut EventType == k, nominal weight
//  POTWeight*SelEffWeight*((cc != 0 || SampleIndex < 4) ? 1 : 0) non-zero, and
//  inside one bin of the sample's binning, lower edge included) and, per
//  GENIE_ branch, how many of them carry a graph GUNDAM turns into a dial
//  (non-empty, and not flat at 1 within 2*FLT_EPSILON as in DialUtils::isFlat).
//  Compare with the "Nb of event affected by parameters" table of a gundamFitter
//  dry run. The binning variables are those of the T2K binning files:
//  Enu_rec/Theta_rec for mu-like, PLep_rec/Theta_rec for e-like.
//
//  Written for the Stage B check (FHC); generalised to all 5 samples for Stage C.
// =============================================================================
#include <cfloat>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TChain.h"
#include "TGraph.h"

struct B { double e0, e1, t0, t1; };
const char *kBinFile[5] = {"binning_T2K_FHC_MuLike.txt", "binning_T2K_FHC_ELike.txt", "binning_T2K_FHC_EPiLike.txt",
                           "binning_T2K_RHC_MuLike.txt", "binning_T2K_RHC_ELike.txt"};
const bool kUsesP[5] = {false, true, true, false, true};   // e-like in PLep_rec, mu-like in Enu_rec

std::vector<B> ReadBins(const std::string &f) {
  std::vector<B> v; std::ifstream in(f); std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line.rfind("variables", 0) == 0) continue;
    std::istringstream ss(line); B b; if (ss >> b.e0 >> b.e1 >> b.t0 >> b.t1) v.push_back(b);
  }
  return v;
}

// Event selection shared with bin_check.C. Returns the sample index (0..4) or -1.
struct Sel {
  std::vector<B> bins[5];
  int et, si; float erec, th, plep, etrue; Char_t cc; double potw, selw;
  void Init(const std::string &bdir, TChain &t) {
    for (int k = 0; k < 5; ++k) bins[k] = ReadBins(bdir + "/" + kBinFile[k]);
    t.SetBranchAddress("EventType", &et); t.SetBranchAddress("Enu_rec", &erec); t.SetBranchAddress("Theta_rec", &th);
    t.SetBranchAddress("PLep_rec", &plep); t.SetBranchAddress("SampleIndex", &si); t.SetBranchAddress("cc", &cc);
    t.SetBranchAddress("POTWeight", &potw); t.SetBranchAddress("SelEffWeight", &selw); t.SetBranchAddress("Enu_true", &etrue);
  }
  double Weight() const { return potw * selw * ((cc != 0 || si < 4) ? 1 : 0); }
  int Sample() const {
    if (et < 0 || et > 4 || Weight() == 0) return -1;
    const double x = kUsesP[et] ? plep : erec;
    for (auto &b : bins[et]) if (x >= b.e0 && x < b.e1 && th >= b.t0 && th < b.t1) return et;
    return -1;
  }
};

void count_dials(const char *files, const char *bdir) {
  TChain t("FakeSelect");
  std::istringstream fs(files); std::string f; while (fs >> f) t.Add(f.c_str());
  Sel sel; sel.Init(bdir, t);
  std::vector<std::string> names; std::vector<TGraph*> g;
  for (auto *o : *t.GetListOfBranches()) if (std::string(o->GetName()).rfind("GENIE_", 0) == 0) names.push_back(o->GetName());
  g.assign(names.size(), nullptr);
  for (size_t d = 0; d < names.size(); ++d) t.SetBranchAddress(names[d].c_str(), &g[d]);
  long long nin[5] = {0, 0, 0, 0, 0};
  std::vector<std::vector<long long>> nd(names.size(), std::vector<long long>(5, 0));
  for (Long64_t i = 0; i < t.GetEntries(); ++i) {
    t.GetEntry(i);
    const int s = sel.Sample();
    if (s < 0) continue;
    ++nin[s];
    for (size_t d = 0; d < names.size(); ++d) {
      const int n = g[d]->GetN(); if (n == 0) continue;
      bool flat = true;
      for (int k = 0; k < n; ++k) if (std::abs(g[d]->GetY()[k] - g[d]->GetY()[0]) > 2 * FLT_EPSILON) flat = false;
      if (flat && std::abs(g[d]->GetY()[0] - 1) < 2 * FLT_EPSILON) continue;
      ++nd[d][s];
    }
  }
  printf("events in binning: %lld %lld %lld %lld %lld\n", nin[0], nin[1], nin[2], nin[3], nin[4]);
  for (size_t d = 0; d < names.size(); ++d) {
    long long tot = 0; for (int s = 0; s < 5; ++s) tot += nd[d][s];
    printf("%-26s %8lld %8lld %8lld %8lld %8lld %8lld\n", names[d].c_str(), tot,
           nd[d][0], nd[d][1], nd[d][2], nd[d][3], nd[d][4]);
  }
}

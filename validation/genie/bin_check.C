// =============================================================================
//  bin_check.C  -  GUNDAM's prefit one-sigma variations vs the graphs.
//
//  Usage (validation/genie/run_checks.sh gundamcheck runs it):
//      root -l -b -q 'bin_check.C+("<fakeselect files>", "<binning dir>", "<xsecCovGENIE.root>", "<gundam one-sigma output.root>")'
//
//  For every GENIE parameter: the relative deviation (%) of each sample in bins
//  of true Enu (inputs/output/binning_Enu_true.txt) when the parameter is moved
//  to prior + 1 sigma, computed directly from the graphs in the FakeSelect
//  files (TGraph::Eval, which is exact here: +1 sigma is a spline knot and the
//  norm graphs are straight lines), with GUNDAM's selection, weights and
//  binning (count_dials.C) and no oscillation. Compared bin by bin with the
//  "Relative Deviation (%)" histograms GUNDAM writes under
//  FitterEngine/preFit/oneSigma/.../canvas/Enu_true/samples_n0_TCanvas, from a
//  run with override_noOsc.yaml (GENIE GUNDAM config). GUNDAM's scan moves
//  parameters by +1 sigma only.
//
//  Written for the Stage B check (FHC); generalised to all 5 samples for Stage C.
// =============================================================================
#include "count_dials.C"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TMatrixD.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TPad.h"
#include "TVectorD.h"

void bin_check(const char *files, const char *bdir, const char *covf, const char *gfile) {
  const double eb[12] = {0, 0.4, 0.5, 0.6, 0.7, 1.0, 1.5, 2.5, 3.5, 5.0, 7.0, 30.0};
  TFile fc(covf);
  auto *nm = (TObjArray*)fc.Get("xsec_param_names"); auto *pr = (TVectorD*)fc.Get("xsec_param_priors");
  auto *cv = (TMatrixD*)fc.Get("xsec_param_cov");
  TChain t("FakeSelect");
  std::istringstream fs(files); std::string f; while (fs >> f) t.Add(f.c_str());
  Sel sel; sel.Init(bdir, t);
  const int n = nm->GetEntries();
  std::vector<TGraph*> g(n, nullptr); std::vector<double> hi(n);
  for (int d = 0; d < n; ++d) {
    t.SetBranchAddress(((TObjString*)nm->At(d))->GetName(), &g[d]);
    hi[d] = (*pr)[d] + std::sqrt((*cv)(d, d));
  }
  double r0[5][11] = {};
  std::vector<std::vector<std::vector<double>>> rp(n, std::vector<std::vector<double>>(5, std::vector<double>(11, 0)));
  for (Long64_t i = 0; i < t.GetEntries(); ++i) {
    t.GetEntry(i);
    const int s = sel.Sample();
    if (s < 0) continue;
    int ib = -1; for (int k = 0; k < 11; ++k) if (sel.etrue >= eb[k] && sel.etrue < eb[k+1]) ib = k;
    if (ib < 0) continue;
    const double w = sel.Weight();
    r0[s][ib] += w;
    for (int d = 0; d < n; ++d) rp[d][s][ib] += w * (g[d]->GetN() > 0 ? std::max(0.0, g[d]->Eval(hi[d])) : 1.0);
  }
  TFile fg(gfile);
  double worst = 0; std::string wname;
  int ncompared = 0;
  for (int d = 0; d < n; ++d) {
    TString pn = Form("#%d_%s", d, ((TObjString*)nm->At(d))->GetName());
    auto c = (TCanvas*)fg.Get("FitterEngine/preFit/oneSigma/Cross-Section Parameters/" + pn + "/canvas/Enu_true/samples_n0_TCanvas");
    if (!c) { printf("[FAIL] no one-sigma canvas for %s\n", pn.Data()); continue; }
    double wd = 0; int ip = 0;
    for (auto p : *c->GetListOfPrimitives()) {     // one pad per sample, in sample order
      auto pad = dynamic_cast<TPad*>(p); if (!pad || ip > 4) continue;
      auto h = (TH1*)pad->GetListOfPrimitives()->At(0);
      if (h) for (int k = 0; k < 11; ++k) {
        const double mine = r0[ip][k] > 0 ? 100 * (rp[d][ip][k] / r0[ip][k] - 1) : 0;
        wd = std::max(wd, std::abs(mine - h->GetBinContent(k + 1)));
        ++ncompared;
      }
      ++ip;
    }
    printf("%-28s max|diff| %.2e %%\n", pn.Data(), wd);
    if (wd > worst) { worst = wd; wname = pn.Data(); }
  }
  printf("largest per-bin |difference| over %d (parameter, sample, bin) values: %.2e %% (%s)\n",
         ncompared, worst, wname.c_str());
}

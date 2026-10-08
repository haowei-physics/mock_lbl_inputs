// =============================================================================
//  MakeGenieXsecCov.C  -  GENIE cross-section parameter file for GUNDAM, in the
//                         layout of mock_lbl_inputs' xsecCovHK.root.
//
//  Usage:
//      root -l -b -q 'MakeGenieXsecCov.C("<spline1.root> <spline2.root> ...", "<out.root>")'
//
//  Reads the GenieDialMeta tree that genie_splines.cxx writes into every
//  <stem>_spline.root (GENIE's counterpart of T2KReWeight's TweakMetaData) and
//  writes
//      xsec_param_names   TObjArray of TObjString (single key)
//      xsec_param_priors  TVectorD   nominal, physical units
//      xsec_param_lower   TVectorD
//      xsec_param_upper   TVectorD
//      xsec_param_cov     TMatrixD   diagonal, prior_sigma^2
//  with the same key names as xsecCovHK.root, so the GUNDAM parameter-set yaml
//  only changes the file path.
//
//  Decisions (see the decision log in genie_splines.cxx):
//   * priors = the stock G18_10b_02_11b values; sigma = GENIE Reweight default
//     (not a tune's systematic error; a starting point);
//   * diagonal: GENIE ships no correlations;
//   * spline dials: bounds = outer knots (-2, +2 sigma), so GUNDAM never
//     extrapolates a spline;
//   * norm dials: prior 1, bounds [0, 1 + 3 sigma]; NormCCQE's Gaussian sigma
//     is 0.20 (the + side of +0.20/-0.15).
//  All of these numbers are decided in genie_splines.cxx and only copied here.
//  Every input file must carry identical definitions; any difference is fatal,
//  since it would mean the samples were produced with different settings.
// =============================================================================

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TMatrixD.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TTree.h"
#include "TVectorD.h"

int MakeGenieXsecCov(std::string spline_files, std::string cov_out)
{
  struct Par { std::string name; double prior, sigma, lo, hi; std::vector<double> knots; };
  std::vector<Par> ref;

  std::istringstream ss(spline_files);
  std::string path;
  int nfile = 0;
  while (ss >> path) {
    TFile f(path.c_str());
    TTree *meta = f.IsZombie() ? nullptr : (TTree*)f.Get("GenieDialMeta");
    if (!meta) { std::cout << "[ERROR] no GenieDialMeta in " << path << std::endl; return 0; }
    std::string *name = nullptr; std::vector<double> *knots = nullptr;
    double nom, sig, lo, hi;
    meta->SetBranchAddress("name", &name);
    meta->SetBranchAddress("nominal", &nom);
    meta->SetBranchAddress("prior_sigma", &sig);
    meta->SetBranchAddress("lower", &lo);
    meta->SetBranchAddress("upper", &hi);
    meta->SetBranchAddress("knots", &knots);
    std::vector<Par> cur;
    for (Long64_t i = 0; i < meta->GetEntries(); ++i) {
      meta->GetEntry(i);
      cur.push_back({*name, nom, sig, lo, hi, *knots});
    }
    if (nfile == 0) ref = cur;
    else {
      bool same = cur.size() == ref.size();
      for (size_t i = 0; same && i < cur.size(); ++i)
        same = cur[i].name == ref[i].name && cur[i].prior == ref[i].prior && cur[i].sigma == ref[i].sigma &&
               cur[i].lo == ref[i].lo && cur[i].hi == ref[i].hi && cur[i].knots == ref[i].knots;
      if (!same) { std::cout << "[ERROR] dial definitions in " << path << " differ from the first file" << std::endl; return 0; }
    }
    ++nfile;
  }
  if (nfile == 0) { std::cout << "[ERROR] no input files" << std::endl; return 0; }

  const int n = ref.size();
  TObjArray names;
  names.SetOwner(true);
  TVectorD priors(n), lower(n), upper(n);
  TMatrixD cov(n, n);
  for (int i = 0; i < n; ++i) {
    names.Add(new TObjString(ref[i].name.c_str()));
    priors[i] = ref[i].prior;
    lower[i]  = ref[i].lo;
    upper[i]  = ref[i].hi;
    cov(i, i) = ref[i].sigma * ref[i].sigma;
  }

  TFile out(cov_out.c_str(), "RECREATE");
  if (out.IsZombie()) { std::cout << "[ERROR] cannot write " << cov_out << std::endl; return 0; }
  names.Write("xsec_param_names", TObject::kSingleKey);
  priors.Write("xsec_param_priors");
  lower.Write("xsec_param_lower");
  upper.Write("xsec_param_upper");
  cov.Write("xsec_param_cov");
  out.Close();

  std::cout << "[INFO] " << n << " GENIE parameters, identical in " << nfile << " file(s)" << std::endl;
  printf("%-26s %10s %10s %8s   %s\n", "parameter", "prior", "sigma", "rel", "bounds");
  for (auto &p : ref)
    printf("%-26s %10.6g %10.6g %7.1f%%   [%.6g, %.6g]\n", p.name.c_str(), p.prior, p.sigma,
           100 * p.sigma / p.prior, p.lo, p.hi);
  std::cout << "[INFO] wrote " << cov_out << std::endl;
  return 1;
}

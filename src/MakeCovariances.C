// Writes the parameter covariances GUNDAM reads through parameterDefinitionFilePath,
// into <out_dir>/{oscProb,det,xsec,flux}/. Run from the repo root.
//
//   root -l -b -q 'src/MakeCovariances.C("<a _Weight.flat.root>","<out_dir>")'   all four
//   root -l -b -q -e '.L src/MakeCovariances.C' -e 'MakeDetCov(0.98,"out.root","out.txt")'
//
// Every matrix is written as TMatrixD, never TMatrixDSym: GenericToolbox picks the
// reader by prefix-matching the class name against "TMatrixT", which also matches
// "TMatrixTSym<double>" and would C-cast a sibling type.
//
// Where a binning file is written next to a covariance, GUNDAM pairs dial i with
// parameter i by position and sorts neither, so both come out of the same loop.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TAxis.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixT.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TTree.h"
#include "TVectorD.h"

// ---------------------------------------------------------------------------
// Oscillation parameters: diagonal external constraints, priors at the values
// the events were generated at (arXiv:2505.15019 Tab. 1).
int MakeOscCov (std::string output_file) {

    struct Par { const char *name; double prior; double sigma; const char *note; };
    const std::vector<Par> par = {
        {"PMNS_SIN_SQUARED_12",       0.307,     0.013,   "solar, PDG"},
        {"PMNS_SIN_SQUARED_13",       0.0218,    0.0007,  "reactor, PDG"},
        {"PMNS_SIN_SQUARED_23",       0.528,     0.021,   "free - HK measures it"},
        {"PMNS_DELTA_MASS_SQUARED_21",7.53e-5,   1.8e-6,  "solar, PDG"},
        {"PMNS_DELTA_MASS_SQUARED_32",2.509e-3,  3.4e-5,  "free - HK measures it"},
        {"PMNS_DELTA_CP",            -1.601,     1.17,    "free - HK measures it"},
        {"PMNS_SIGN_MASS_SQUARED_32", 1.0,       0.1,     "fixed - normal ordering"}
    };
    const int n = (int)par.size();

    TObjArray *names = new TObjArray(n);
    names->SetName("osc_param_names");
    TVectorD priors(n);
    TMatrixD cov(n, n);
    for (int i = 0; i < n; ++i) {
        names->Add(new TObjString(par[i].name));
        priors(i) = par[i].prior;
        cov(i, i) = par[i].sigma * par[i].sigma;
    }

    TFile *FileOut = TFile::Open(output_file.c_str(), "recreate");
    names->Write("osc_param_names", TObject::kSingleKey);
    priors.Write("osc_param_priors");
    cov.Write("osc_param_cov");
    FileOut->Close();

    printf("\n %-30s %12s %12s %10s   %s\n", "parameter", "prior", "sigma", "relative", "role");
    for (int i = 0; i < n; ++i)
        printf(" %-30s %12.5g %12.5g %9.2f%%   %s\n",
               par[i].name, par[i].prior, par[i].sigma,
               par[i].prior != 0.0 ? 100.0 * par[i].sigma / std::fabs(par[i].prior) : 0.0,
               par[i].note);
    printf("\n[INFO] wrote %s\n", output_file.c_str());
    return 1;
}

// ---------------------------------------------------------------------------
// Detector: one normalisation per sample, sigmas from arXiv:2505.15019 Tab. 3
// (detector row, T2K-systematics column). rho correlates the FHC/RHC pair of the
// same topology; 0.98 is what the Tab. 3 FHC/RHC 1Re ratio column implies.
int MakeDetCov (
    double rho,
    std::string cov_out,
    std::string binning_out) {

    if (rho < -1.0 || rho > 1.0) {
        std::cout << "[ERROR] rho must be in [-1, 1], got " << rho << std::endl;
        return 0;
    }

    struct Par { const char *name; int event_type; double sigma; };
    const std::vector<Par> par = {
        {"DET_NORM_FHC_1RMU",    0, 0.021},   // improved: 0.008
        {"DET_NORM_FHC_1RE_0DE", 1, 0.031},   // improved: 0.011
        {"DET_NORM_FHC_1RE_1DE", 2, 0.132},   // improved: 0.049
        {"DET_NORM_RHC_1RMU",    3, 0.019},   // improved: 0.007
        {"DET_NORM_RHC_1RE_0DE", 4, 0.039}    // improved: 0.015
    };
    const int n = (int)par.size();

    // FHC <-> RHC counterparts of the same topology
    const int pair_a[2] = {0, 1};
    const int pair_b[2] = {3, 4};

    TMatrixD cov(n, n);
    for (int i = 0; i < n; ++i) cov(i, i) = par[i].sigma * par[i].sigma;
    for (int p = 0; p < 2; ++p) {
        const int a = pair_a[p], b = pair_b[p];
        const double c = rho * par[a].sigma * par[b].sigma;
        cov(a, b) = c;
        cov(b, a) = c;
    }

    // positive-definite check, since GUNDAM's failure message would not name this file
    TMatrixDSym sym(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) sym(i, j) = 0.5 * (cov(i, j) + cov(j, i));
    TVectorD eig;
    sym.EigenVectors(eig);
    double lo = eig[0], hi = eig[0];
    for (int i = 1; i < n; ++i) { if (eig[i] < lo) lo = eig[i]; if (eig[i] > hi) hi = eig[i]; }
    std::cout << "[INFO] rho = " << rho << "   eigenvalues " << lo << " .. " << hi
              << "   condition " << hi / lo << std::endl;
    if (lo <= 0.0) {
        std::cout << "[ERROR] not positive definite, refusing to write" << std::endl;
        return 0;
    }

    TObjArray names;
    names.SetOwner(true);
    TVectorD priors(n);
    for (int i = 0; i < n; ++i) {
        names.Add(new TObjString(par[i].name));
        priors[i] = 1.0;
    }

    TFile *FileOut = TFile::Open(cov_out.c_str(), "recreate");
    if (!FileOut || FileOut->IsZombie()) {
        std::cout << "[ERROR] cannot write " << cov_out << std::endl;
        return 0;
    }
    names.Write("det_param_names", TObject::kSingleKey);
    priors.Write("det_param_priors");
    cov.Write("det_param_cov");
    FileOut->Close();

    // one bin per parameter; EventType listed once is a point match, not a range
    FILE *fb = fopen(binning_out.c_str(), "w");
    if (!fb) {
        std::cout << "[ERROR] cannot write " << binning_out << std::endl;
        return 0;
    }
    fprintf(fb, "# Detector normalisation parameters, one per sample, in the order of\n");
    fprintf(fb, "# det_param_cov. GUNDAM pairs them by position. Do not reorder.\n");
    fprintf(fb, "\n");
    fprintf(fb, "variables: EventType\n");
    for (int i = 0; i < n; ++i)
        fprintf(fb, "%d    # %-22s sigma %5.1f%%\n",
                par[i].event_type, par[i].name, 100.0 * par[i].sigma);
    fclose(fb);

    std::cout << "\n[INFO] " << n << " detector normalisation parameters" << std::endl;
    printf("%-24s %10s %8s\n", "parameter", "EventType", "sigma");
    for (int i = 0; i < n; ++i)
        printf("%-24s %10d %7.1f%%\n", par[i].name, par[i].event_type, 100.0 * par[i].sigma);
    std::cout << "\n[INFO] wrote " << cov_out << std::endl;
    std::cout << "[INFO] wrote " << binning_out << std::endl;
    return 1;
}

// ---------------------------------------------------------------------------
// Cross-section: diagonal, one sigma per dial recovered from the T2KReWeight tweak
// points in TweakMetaData. The multiples below are copied from neut_dial_tweaks.toml;
// the macro checks every knot implies the same sigma. Bounds are the outer knots,
// beyond which the spline has no response.
int MakeXsecCov (
    std::string weight_file,
    std::string cov_out) {

    struct Par { const char *name; std::vector<double> mult; };
    const std::vector<Par> par = {
        {"MaCCQE",                 {-2, -1, 1, 2}},
        {"MaRES",                  {-2, -1, 1, 2}},
        {"CA5RES",                 {-2, -1, 1, 2}},
        {"BgSclRES",               {-2, -1, 1, 2}},
        {"NucleonFSI_Double",      {-1, -0.5, 0.5, 1}},
        {"NucleonFSI_Elastic",     {-1, -0.5, 0.5, 1}},
        {"NucleonFSI_Single",      {-1, -0.5, 0.5, 1}},
        {"NucleonFSI_Total",       {-1, -0.5, 0.5, 1}},
        {"PionFSI_AbsProb",        {-1, -0.5, 0.5, 1}},
        {"PionFSI_QELowMomProb",   {-1, -0.5, 0.5, 1}},
        {"PionFSI_QEHighMomProb",  {-1, -0.5, 0.5, 1}},
        {"PionFSI_CExLowMomProb",  {-1, -0.5, 0.5, 1}},
        {"PionFSI_CExHighMomProb", {-1, -0.5, 0.5, 1}},
        {"PionFSI_InelProb",       {-0.9, -0.5, 0.5, 1}},
        {"SRCFrac_O",              {0.25, 0.5}}
    };
    const int n = (int)par.size();

    TFile *FileW = TFile::Open(weight_file.c_str(), "read");
    if (!FileW || FileW->IsZombie()) {
        std::cout << "[ERROR] could not open " << weight_file << std::endl;
        return 0;
    }
    TTree *meta = (TTree*)FileW->Get("TweakMetaData");
    if (!meta) {
        std::cout << "[ERROR] no TweakMetaData in " << weight_file << std::endl;
        return 0;
    }

    std::vector<double> nominal(n, 0.0), sigma(n, 0.0), lo(n, 0.0), hi(n, 0.0);
    for (int i = 0; i < n; ++i) {
        std::vector<double> *tweaks = nullptr;
        double from = 0.0;
        if (!meta->GetBranch(par[i].name)) {
            std::cout << "[ERROR] " << par[i].name << " is not in TweakMetaData" << std::endl;
            return 0;
        }
        meta->SetBranchAddress(par[i].name, &tweaks);
        meta->SetBranchAddress((std::string(par[i].name) + "_from").c_str(), &from);
        meta->GetEntry(0);

        if (!tweaks || tweaks->size() != par[i].mult.size()) {
            std::cout << "[ERROR] " << par[i].name << " has " << (tweaks ? tweaks->size() : 0)
                      << " tweak points but the table expects " << par[i].mult.size() << std::endl;
            return 0;
        }

        double s_sum = 0.0, s_min = 1e300, s_max = -1e300;
        for (size_t k = 0; k < tweaks->size(); ++k) {
            const double s = ((*tweaks)[k] - from) / par[i].mult[k];
            s_sum += s;
            if (s < s_min) s_min = s;
            if (s > s_max) s_max = s;
        }
        const double s_mean = s_sum / tweaks->size();
        const double spread = (s_mean != 0.0) ? (s_max - s_min) / std::abs(s_mean) : 0.0;
        if (spread > 1e-3) {
            std::cout << "[ERROR] " << par[i].name << ": knots imply sigma " << s_min
                      << " .. " << s_max << ", the multiples do not match the file" << std::endl;
            return 0;
        }

        nominal[i] = from;
        sigma[i]   = s_mean;
        lo[i] = *std::min_element(tweaks->begin(), tweaks->end());
        hi[i] = *std::max_element(tweaks->begin(), tweaks->end());
        if (from < lo[i]) lo[i] = from;
        if (from > hi[i]) hi[i] = from;

        meta->ResetBranchAddresses();
    }
    FileW->Close();

    TMatrixD cov(n, n);
    for (int i = 0; i < n; ++i) cov(i, i) = sigma[i] * sigma[i];

    TObjArray names;
    names.SetOwner(true);
    TVectorD priors(n), lbound(n), ubound(n);
    for (int i = 0; i < n; ++i) {
        names.Add(new TObjString(par[i].name));
        priors[i] = nominal[i];
        lbound[i] = lo[i];
        ubound[i] = hi[i];
    }

    TFile *FileOut = TFile::Open(cov_out.c_str(), "recreate");
    if (!FileOut || FileOut->IsZombie()) {
        std::cout << "[ERROR] cannot write " << cov_out << std::endl;
        return 0;
    }
    names.Write("xsec_param_names", TObject::kSingleKey);
    priors.Write("xsec_param_priors");
    lbound.Write("xsec_param_lower");
    ubound.Write("xsec_param_upper");
    cov.Write("xsec_param_cov");
    FileOut->Close();

    std::cout << "\n[INFO] " << n << " cross-section parameters" << std::endl;
    printf("%-24s %10s %10s %8s   %s\n", "parameter", "nominal", "sigma", "rel", "bounds");
    for (int i = 0; i < n; ++i)
        printf("%-24s %10.4f %10.4f %7.1f%%   [%.4f, %.4f]\n",
               par[i].name, nominal[i], sigma[i],
               nominal[i] != 0.0 ? 100.0 * sigma[i] / nominal[i] : 0.0, lo[i], hi[i]);
    std::cout << "\n[INFO] wrote " << cov_out << std::endl;
    return 1;
}

// ---------------------------------------------------------------------------
// Flux: the SK half of the T2K 2020 public flux covariance (fractional, so priors
// are 1.0). Two files: C, the SK block on its own (no near detector), and
// C - s * B^T A^-1 B, the SK block after an ND280 measurement of strength s in
// [0,1] (A = ND280 block, B = cross block). Both share one binning file.
int MakeFluxCov (
    double nd_strength,
    std::string flux_cov_file,
    std::string cov_out,
    std::string cov_nd_out,
    std::string binning_out) {

    if (nd_strength < 0.0 || nd_strength > 1.0) {
        printf("[ERROR] nd_strength must be in [0,1], got %g\n", nd_strength);
        return 0;
    }

    // the eight SK blocks in released order; isRHC and flux_pdg are the FakeSelect
    // branches that identify each block
    struct Block { const char *axis; const char *tag; int isRHC; int flux_pdg; };
    const std::vector<Block> sk_block = {
        {"sk_numode_numu_bins",   "FHC_NUMU",   0,  14},
        {"sk_numode_numub_bins",  "FHC_NUMUB",  0, -14},
        {"sk_numode_nue_bins",    "FHC_NUE",    0,  12},
        {"sk_numode_nueb_bins",   "FHC_NUEB",   0, -12},
        {"sk_anumode_numu_bins",  "RHC_NUMU",   1,  14},
        {"sk_anumode_numub_bins", "RHC_NUMUB",  1, -14},
        {"sk_anumode_nue_bins",   "RHC_NUE",    1,  12},
        {"sk_anumode_nueb_bins",  "RHC_NUEB",   1, -12}
    };
    // the ND280 blocks give the row offset of the SK half
    const std::vector<std::string> nd280_axis = {
        "nd280_numode_numu_bins",  "nd280_numode_numub_bins",
        "nd280_numode_nue_bins",   "nd280_numode_nueb_bins",
        "nd280_anumode_numu_bins", "nd280_anumode_numub_bins",
        "nd280_anumode_nue_bins",  "nd280_anumode_nueb_bins"
    };

    TFile *FileIn = TFile::Open(flux_cov_file.c_str(), "read");
    if (!FileIn || FileIn->IsZombie()) {
        printf("[ERROR] could not open %s\n", flux_cov_file.c_str());
        return 0;
    }
    TMatrixT<double> *full = (TMatrixT<double> *)FileIn->Get("total_flux_cov");
    if (!full) {
        printf("[ERROR] no total_flux_cov in %s\n", flux_cov_file.c_str());
        return 0;
    }

    int offset = 0;
    for (size_t i = 0; i < nd280_axis.size(); ++i) {
        TAxis *a = (TAxis *)FileIn->Get(nd280_axis[i].c_str());
        if (!a) {
            printf("[ERROR] no %s in %s\n", nd280_axis[i].c_str(), flux_cov_file.c_str());
            return 0;
        }
        offset += a->GetNbins();
    }

    struct Par { std::string name; int row; int isRHC; int flux_pdg; double lo; double hi; };
    std::vector<Par> par;
    int row = offset;
    for (size_t i = 0; i < sk_block.size(); ++i) {
        TAxis *a = (TAxis *)FileIn->Get(sk_block[i].axis);
        if (!a) {
            printf("[ERROR] no %s in %s\n", sk_block[i].axis, flux_cov_file.c_str());
            return 0;
        }
        for (int b = 1; b <= a->GetNbins(); ++b) {
            Par p;
            p.name     = Form("FLUX_%s_%02d", sk_block[i].tag, b - 1);
            p.row      = row;
            p.isRHC    = sk_block[i].isRHC;
            p.flux_pdg = sk_block[i].flux_pdg;
            p.lo       = a->GetBinLowEdge(b);
            p.hi       = a->GetBinUpEdge(b);
            par.push_back(p);
            ++row;
        }
    }
    const int n = (int)par.size();

    if (row != full->GetNrows()) {
        printf("[ERROR] SK block ends at row %d but the matrix is %d x %d\n",
               row, full->GetNrows(), full->GetNcols());
        return 0;
    }
    if (offset != n) {
        printf("[ERROR] ND280 half is %d bins and the SK half is %d\n", offset, n);
        return 0;
    }

    TMatrixD A(n, n), B(n, n), C(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            A(i, j) = (*full)(i, j);
            B(i, j) = (*full)(i, par[j].row);
            C(i, j) = (*full)(par[i].row, par[j].row);
        }

    TMatrixD Ainv(A);
    double detA = 0.0;
    Ainv.Invert(&detA);
    if (detA == 0.0) {
        printf("[ERROR] the ND280 block is singular\n");
        return 0;
    }
    TMatrixD Bt(TMatrixD::kTransposed, B);
    TMatrixD G = Bt * Ainv * B;
    TMatrixD Cnd(C);
    Cnd -= nd_strength * G;

    // symmetrise and check positive definite; A^-1 is numerical
    TMatrixDSym CndSym(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) CndSym(i, j) = 0.5 * (Cnd(i, j) + Cnd(j, i));
    TMatrixDSymEigen eig(CndSym);
    const double lambda_min = eig.GetEigenValues()(n - 1);
    const double lambda_max = eig.GetEigenValues()(0);
    if (lambda_min <= 0.0) {
        printf("[ERROR] constrained matrix is not positive definite: lambda_min = %.4e\n", lambda_min);
        return 0;
    }

    TObjArray *names = new TObjArray(n);
    names->SetName("flux_param_names");
    TVectorD priors(n);
    for (int i = 0; i < n; ++i) {
        names->Add(new TObjString(par[i].name.c_str()));
        priors(i) = 1.0;
    }
    TMatrixD CndOut(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) CndOut(i, j) = CndSym(i, j);

    // same object names in both files, so configParSet.yaml switches by path only
    TFile *FileOut = TFile::Open(cov_out.c_str(), "recreate");
    names->Write("flux_param_names", TObject::kSingleKey);
    priors.Write("flux_param_priors");
    C.Write("flux_param_cov");
    FileOut->Close();

    TFile *FileOutND = TFile::Open(cov_nd_out.c_str(), "recreate");
    names->Write("flux_param_names", TObject::kSingleKey);
    priors.Write("flux_param_priors");
    CndOut.Write("flux_param_cov");
    FileOutND->Close();

    // isRHC and flux_pdg listed once are point matches; Enu_true twice is a range
    FILE *txt = fopen(binning_out.c_str(), "w");
    if (!txt) {
        printf("[ERROR] could not write %s\n", binning_out.c_str());
        return 0;
    }
    fprintf(txt, "# SK flux normalisation bins, T2K 2020 public release, in the order of\n");
    fprintf(txt, "# flux_param_cov. GUNDAM pairs them by position. Do not reorder.\n");
    fprintf(txt, "variables: isRHC flux_pdg Enu_true Enu_true\n");
    std::string current;
    for (int i = 0; i < n; ++i) {
        std::string tag = par[i].name.substr(5, par[i].name.size() - 8);
        if (tag != current) {
            current = tag;
            fprintf(txt, "# %s\n", current.c_str());
        }
        fprintf(txt, "%d %d %9.4f %9.4f\n", par[i].isRHC, par[i].flux_pdg, par[i].lo, par[i].hi);
    }
    fclose(txt);

    printf("\n SK flux normalisations, rows %d-%d of the %d x %d released matrix\n",
           par.front().row, par.back().row, full->GetNrows(), full->GetNcols());
    printf(" ND constraint strength s = %.2f\n\n", nd_strength);
    printf(" %-20s %8s %10s %12s %10s %10s %8s\n",
           "parameter", "isRHC", "flux_pdg", "Enu_true", "no ND", "with ND", "ratio");
    double sum_marg = 0.0, sum_nd = 0.0, lo_nd = 1e9, hi_nd = 0.0;
    for (int i = 0; i < n; ++i) {
        const double s_marg = std::sqrt(C(i, i));
        const double s_nd   = std::sqrt(CndSym(i, i));
        sum_marg += s_marg;
        sum_nd   += s_nd;
        if (s_nd < lo_nd) lo_nd = s_nd;
        if (s_nd > hi_nd) hi_nd = s_nd;
        printf(" %-20s %8d %10d %5.2f-%5.2f %9.2f%% %9.2f%% %8.2f\n",
               par[i].name.c_str(), par[i].isRHC, par[i].flux_pdg,
               par[i].lo, par[i].hi, 100.0 * s_marg, 100.0 * s_nd, s_nd / s_marg);
    }
    printf("\n mean sigma       %.2f%%  ->  %.2f%%\n", 100.0 * sum_marg / n, 100.0 * sum_nd / n);
    printf(" range with ND    %.2f%% to %.2f%%\n", 100.0 * lo_nd, 100.0 * hi_nd);
    printf(" eigenvalues      lambda_min = %.4e  lambda_max = %.4e  cond = %.3e\n",
           lambda_min, lambda_max, lambda_max / lambda_min);
    printf("\n[INFO] %d parameters\n", n);
    printf("[INFO] wrote %s\n", cov_out.c_str());
    printf("[INFO] wrote %s\n", cov_nd_out.c_str());
    printf("[INFO] wrote %s\n", binning_out.c_str());
    return 1;
}

// ---------------------------------------------------------------------------
int MakeCovariances (std::string weight_file, std::string out_dir) {
    int ok = 1;
    ok &= MakeOscCov(out_dir + "/oscProb/oscCovHK.root");
    ok &= MakeDetCov(0.0, out_dir + "/det/detCovHK.root", out_dir + "/det/detBinning.txt");
    ok &= MakeXsecCov(weight_file, out_dir + "/xsec/xsecCovHK.root");
    ok &= MakeFluxCov(0.8, "lbl_config/hk/t2kflux_2020_covariance.root", out_dir + "/flux/fluxCovHK.root", out_dir + "/flux/fluxCovHK_ND.root", out_dir + "/flux/fluxBinning.txt");
    return ok;
}

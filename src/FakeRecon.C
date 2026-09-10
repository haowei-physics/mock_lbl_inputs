#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TKey.h"
#include "TRandom3.h"
#include "TTree.h"

#include "NuFast_LBL.cpp"

// The goal of this script is to compute the POT, oscillation and selection
// efficiency weights, and make fake recon variables, for the input flattree.
// The weights are written out one per branch and are never combined here.
int FakeRecon (std::string input_file,
               std::string output_file,
               std::string recon_smearing_matrix = "lbl_config/hk/smearing.root",
               std::string response_file = "lbl_config/hk/response.root",
               UInt_t rng_seed = 1234,

               std::string flux_file = "lbl_config/hk/t2kflux_2020_plus250kA_uniform.root",
               std::string flux_hist = "enu_sk_numu",
               // 10 years of HK running is 2.7e22 POT total, split 1:3 nu:nubar,
               // so a nu-mode sample only sees a quarter of it.
               double pot = 6.75e21,
               double n_nucleons = 1.132e35,  // 188 kt HK fiducial volume of water

               int flux_pdg = 14) {

    // T2K Asimov A oscillation parameters, normal ordering. NuFast wants
    // Dmsq31 = Dmsq32 + Dmsq21.
    const double s12sq  = 0.307;
    const double s13sq  = 0.0218;
    const double s23sq  = 0.528;
    const double delta  = -1.601;      // radians
    const double Dmsq21 = 7.53e-5;     // eV^2
    const double Dmsq31 = 2.509e-3 + Dmsq21;
    const double baseline_km = 295.0;
    const double rhoYe  = 2.6 * 0.5;   // Earth crust density x electron fraction



    // True masses, GeV.
    const double m_mu = 0.105658;
    const double m_e  = 0.000511;


    const double michel_eff_plus  = 0.96;   // mu+
    const double michel_eff_minus = 0.83;   // mu-, capture included
    const double michel_eff_pip   = 0.45;   // pi+ -> mu+ -> e+
    const double flux_pot_norm     = 1.0e21;
    const double flux_bin_norm_GeV = 0.05;

    // ---------------------------------------------------------------- input
    TFile *FileIn = TFile::Open(input_file.c_str(), "read");
    if (!FileIn || FileIn->IsZombie()) {
        std::cout << "[ERROR] could not open input file " << input_file << std::endl;
        return 0;
    }

    TTree *flattree = (TTree*)FileIn->Get("FlatTree_VARS");
    if (!flattree) {
        std::cout << "[ERROR] no FlatTree_VARS in " << input_file << std::endl;
        return 0;
    }

    const int kMaxFSP = 100;
    double fScaleFactor;
    float  Enu_true, ELep, CosLep, Q2;
    int    PDGnu;
    Char_t cc;
    int    PDGLep;
    Bool_t flagCCQE;
    int    nfsp;
    int    pdg[kMaxFSP];
    flattree->SetBranchAddress("nfsp",         &nfsp);
    flattree->SetBranchAddress("pdg",          pdg);
    flattree->SetBranchAddress("cc",           &cc);
    flattree->SetBranchAddress("PDGLep",       &PDGLep);
    flattree->SetBranchAddress("flagCCQE",     &flagCCQE);
    flattree->SetBranchAddress("fScaleFactor", &fScaleFactor);
    flattree->SetBranchAddress("Enu_true",     &Enu_true);
    flattree->SetBranchAddress("ELep",         &ELep);
    flattree->SetBranchAddress("CosLep",       &CosLep);
    flattree->SetBranchAddress("PDGnu",        &PDGnu);
    flattree->SetBranchAddress("Q2",           &Q2);



    // ------------------------------------------------------------- flux norm
    TFile *FileFlux = TFile::Open(flux_file.c_str(), "read");
    if (!FileFlux || FileFlux->IsZombie()) {
        std::cout << "[ERROR] could not open flux file " << flux_file << std::endl;
        return 0;
    }

    TH1 *hflux = 0;
    if (!flux_hist.empty()) {
        hflux = (TH1*)FileFlux->Get(flux_hist.c_str());
    } else {
        TIter nextkey(FileFlux->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)nextkey())) {
            TObject *obj = key->ReadObj();
            if (obj->InheritsFrom("TH1")) { hflux = (TH1*)obj; break; }
        }
    }
    if (!hflux) {
        std::cout << "[ERROR] no flux histogram found in " << flux_file << std::endl;
        return 0;
    }

    // nu / cm^2 / POT
    const double flux_per_pot =
        hflux->Integral("width") / flux_bin_norm_GeV / flux_pot_norm;
    std::cout << "[INFO] flux hist " << hflux->GetName()
              << " -> " << flux_per_pot << " nu/cm^2/POT" << std::endl;

    // ------------------------------------------------------------- smearing
    // Both histograms are binned in the true quantity, with the bin content
    // giving the absolute Gaussian width to apply there. See MakeSmearing.C.
    TFile *FileSmear = 0;
    TH1D  *h_dp = 0, *h_dcos = 0;
    if (!recon_smearing_matrix.empty()) {
        FileSmear = TFile::Open(recon_smearing_matrix.c_str(), "read");
        if (!FileSmear || FileSmear->IsZombie()) {
            std::cout << "[ERROR] could not open smearing file "
                      << recon_smearing_matrix << std::endl;
            return 0;
        }
        h_dp   = (TH1D*)FileSmear->Get("h_dp_vs_p");
        h_dcos = (TH1D*)FileSmear->Get("h_dcos_vs_cos");
        if (!h_dp || !h_dcos) {
            std::cout << "[ERROR] " << recon_smearing_matrix
                      << " is missing h_dp_vs_p or h_dcos_vs_cos" << std::endl;
            return 0;
        }
        std::cout << "[INFO] smearing widths from " << recon_smearing_matrix << std::endl;
    } else {
        std::cout << "[INFO] no smearing file given, recon variables = truth" << std::endl;
    }

    // ------------------------------------------------------------- response
    // The FCFV and single-ring efficiency against true neutrino energy, one
    // curve per interaction category. The category order here has to match
    // MakeResponse.C exactly, because RespCat is written as an index into it and
    // FakeSelect uses that index to pick the mixing curves.
    const int kNCat = 9;
    const char *cat_name[kNCat] = {"numuCCQE", "numubarCCQE",
                                   "numuCCnonQE", "numubarCCnonQE",
                                   "nueCCQE", "nuebarCCQE",
                                   "nueCCnonQE", "nuebarCCnonQE", "NC"};
    TFile *FileResp = TFile::Open(response_file.c_str(), "read");
    if (!FileResp || FileResp->IsZombie()) {
        std::cout << "[ERROR] could not open response file " << response_file
                  << " - run MakeResponse.C first" << std::endl;
        return 0;
    }
    TH1D *h_eff[kNCat];
    for (int c = 0; c < kNCat; ++c) {
        h_eff[c] = (TH1D*)FileResp->Get(Form("h_eff_%s", cat_name[c]));
        if (!h_eff[c]) {
            std::cout << "[ERROR] " << response_file << " has no h_eff_"
                      << cat_name[c] << std::endl;
            return 0;
        }
    }
    std::cout << "[INFO] response from " << response_file << ", "
              << kNCat << " categories" << std::endl;

    std::cout << "[INFO] flux PDG " << flux_pdg
              << ": CC weight is P(flux flavour -> PDGnu), and NC is dropped if"
              << " those differ" << std::endl;

    TRandom3 rng(rng_seed);

    // --------------------------------------------------------------- output
    TFile *FileOut = TFile::Open(output_file.c_str(), "recreate");
    FileOut->cd();

    // Empty clone: keeps every input branch, addresses are shared with the
    // input tree so a GetEntry()/Fill() pair copies the event across.
    TTree *outtree = flattree->CloneTree(0);

    // Kept as separate branches, never multiplied together here.
    double POTWeight, OscWeight, SelEffWeight;
    float  PLep_rec, CosLep_rec;
    int    true_ME, recon_ME, RespCat;
    // The interaction category, as an index into the list MakeResponse.C writes.
    // Written here rather than recomputed downstream so that the efficiency and
    // the mixing cannot end up disagreeing about what an event is.
    outtree->Branch("RespCat",       &RespCat,       "RespCat/I");
    outtree->Branch("true_ME",       &true_ME,       "true_ME/I");
    outtree->Branch("recon_ME",      &recon_ME,      "recon_ME/I");
    outtree->Branch("POTWeight",     &POTWeight,     "POTWeight/D");
    outtree->Branch("OscWeight",     &OscWeight,     "OscWeight/D");
    outtree->Branch("SelEffWeight",  &SelEffWeight,  "SelEffWeight/D");
    outtree->Branch("PLep_rec",      &PLep_rec,      "PLep_rec/F");
    outtree->Branch("CosLep_rec",    &CosLep_rec,    "CosLep_rec/F");

    // ----------------------------------------------------------- event loop
    const Long64_t nentries = flattree->GetEntries();
    double sum_weight = 0.0, sum_weight_cc = 0.0, sum_weight_osc = 0.0;
    Long64_t n_dropped_no_lepton = 0;

    // NuFast is vectorised over energy; one entry at a time is plenty here.
    std::vector<double> nufast_E(1);
    std::vector<std::array<std::array<double, 3>, 3>> nufast_P(1);
    for (Long64_t i = 0; i < nentries; ++i) {
        flattree->GetEntry(i);

        // Dropped: Pauli blocked events. The struck nucleon comes out into a
        // momentum state below the Fermi surface that is already occupied, so
        // there is no nucleon to follow out of the nucleus. That is a property
        // of the nuclear final state and is independent of the lepton vertex -
        // it happens to CC and NC alike.
        //
        // They carry Q2 = -999 because the flat tree writer gives up on them
        // before any of the derived kinematics are computed


        if (Q2 < -999.0) { ++n_dropped_no_lepton; continue; }
        if (cc && (PDGLep == 0 || ELep < 0.0)) { ++n_dropped_no_lepton; continue; }

        // fScaleFactor is <xsec>/N_generated in cm^2/nucleon, so this is
        // xsec x fluence x targets and comes out as a number of interactions
        // in the fiducial volume, with no selection applied.
        POTWeight = fScaleFactor * flux_per_pot * pot * n_nucleons;
        sum_weight += POTWeight;

        // Oscillation channel: initial flavour from the flux, final flavour
        // from the interaction. 0 = e and 1 = mu in NuFast's probability
        // matrix, so nufast_P[0][1][0] is P(numu->nue).
        const int  a_flav = (abs(flux_pdg) == 12) ? 0 : 1;
        const int  b_flav = (abs(PDGnu)    == 12) ? 0 : 1;
        const bool is_appearance = (a_flav != b_flav);

        // CC takes P(flux flavour -> interaction flavour) at the T2K Asimov A
        // point. For an intrinsic sample that is the survival probability; for
        // an appearance sample the two flavours differ and it is the
        // appearance probability instead.
        if (cc) {
            // NuFast takes a negative energy to mean an antineutrino, and the
            // beam species sets the sign.
            nufast_E[0] = (flux_pdg > 0) ? Enu_true : -Enu_true;
            NuFast::Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31,
                                           baseline_km, nufast_E, rhoYe, 0, nufast_P);
            OscWeight = nufast_P[0][a_flav][b_flav];
        } else if (is_appearance) {
            // NC in an appearance sample would be double counting. The NC cross
            // section is flavour blind, so by unitarity the NC events in the
            // intrinsic sample built on this same flux already carry the whole
            // NC rate from it, including the part that arrived as another
            // flavour. Counting them here as well would add them twice.
            OscWeight = 0.0;
        } else {
            // NC in an intrinsic sample carries the full, unoscillated rate,
            // for the reason just given.
            OscWeight = 1.0;
        }

        sum_weight_osc += POTWeight * OscWeight;
        if (cc && PDGnu == 14) sum_weight_cc += POTWeight;

        // The interaction category, and the FCFV + single-ring efficiency read
        // off that category's curve at this event's true energy. Order must
        // match cat_name above and the list in MakeResponse.C.
        if (!cc)                          RespCat = 8;                    // NC
        else if (abs(PDGnu) == 12)        RespCat = (flagCCQE ? 4 : 6) + (PDGnu > 0 ? 0 : 1);
        else                              RespCat = (flagCCQE ? 0 : 2) + (PDGnu > 0 ? 0 : 1);

        {
            TH1D *h = h_eff[RespCat];
            const int ie = std::min(std::max(h->FindBin(Enu_true), 1), h->GetNbinsX());
            SelEffWeight = h->GetBinContent(ie);
        }

        // ------------------------------------------------ fake reconstruction
        const int nfsp_use = std::min(nfsp, kMaxFSP);

        // Smear the lepton kinematics. Enu_rec gets built from these downstream,
        // under whichever hypothesis the mixing assigns, so a non-QE event picks
        // up the usual quasi-elastic reconstruction bias by construction rather
        // than having it put in by hand.
        const double m_lep = (abs(PDGLep) == 13) ? m_mu
                           : (abs(PDGLep) == 11) ? m_e : 0.0;
        const double p_true = (ELep > m_lep) ? std::sqrt(ELep * ELep - m_lep * m_lep) : 0.0;

        double p_smear = p_true, cos_smear = CosLep;
        if (h_dp) {
            const int ip = std::min(std::max(h_dp->FindBin(p_true), 1), h_dp->GetNbinsX());
            p_smear = p_true + rng.Gaus(0.0, h_dp->GetBinContent(ip));
            if (p_smear < 0.0) p_smear = 0.0;

            const int ic = std::min(std::max(h_dcos->FindBin(CosLep), 1), h_dcos->GetNbinsX());
            cos_smear = CosLep + rng.Gaus(0.0, h_dcos->GetBinContent(ic));
            if (cos_smear >  1.0) cos_smear =  1.0;
            if (cos_smear < -1.0) cos_smear = -1.0;
        }

        PLep_rec   = p_smear;
        CosLep_rec = cos_smear;

        // Michel electrons. true_ME counts the decay-capable parents; recon_ME
        // tags each one independently at the efficiency for its charge, which
        // is the binomial thinning a real sub-event search does. Capture is not
        // drawn separately because it is already inside michel_eff_minus.
        true_ME  = 0;
        recon_ME = 0;
        for (int ip = 0; ip < nfsp_use; ++ip) {
            double eff;
            if      (pdg[ip] == -13) eff = michel_eff_plus;
            else if (pdg[ip] ==  13) eff = michel_eff_minus;
            else if (pdg[ip] == 211) eff = michel_eff_pip;
            else continue;                       // pi- is captured, no Michel
            ++true_ME;
            if (rng.Uniform() < eff) ++recon_ME;
        }


        outtree->Fill();

        if (i % 10000 == 0) std::cout << "[INFO] entry " << i << " / " << nentries << std::endl;
    }

    // Carry the flux and event histograms over so the output stays a valid
    // flattree.
    TH1D *flux = (TH1D*)FileIn->Get("FlatTree_FLUX");
    TH1D *evt  = (TH1D*)FileIn->Get("FlatTree_EVT");
    FileOut->cd();
    if (flux) flux->Clone("FlatTree_FLUX")->Write();
    if (evt)  evt->Clone("FlatTree_EVT")->Write();

    outtree->Write("", TObject::kOverwrite);
    FileOut->Close();
    if (FileSmear) FileSmear->Close();
    FileResp->Close();
    FileFlux->Close();
    FileIn->Close();

    if (n_dropped_no_lepton > 0)
        std::cout << "[INFO] dropped " << n_dropped_no_lepton << " entries ("
                  << 100.0 * n_dropped_no_lepton / nentries
                  << "%) that are Pauli blocked or have no lepton - this count "
                     "must equal nentries minus the T2KReWeight tree entries"
                  << std::endl;

    // A zero here means a weight never got assigned in the loop. It is worth a
    // hard failure rather than a printed zero, because every rate downstream
    // silently becomes zero and the tables still render.
    if (sum_weight <= 0.0 || sum_weight_osc <= 0.0) {
        std::cout << "[ERROR] sum_weight " << sum_weight << ", sum_weight_osc "
                  << sum_weight_osc << " - a weight was never assigned" << std::endl;
        return 0;
    }


    std::cout << "[INFO] wrote " << nentries - n_dropped_no_lepton
              << " entries to " << output_file
              << " (from " << nentries << " generated)" << std::endl;
    std::cout << "[INFO] --- interaction rates for " << pot << " POT ---" << std::endl;
    std::cout << "[INFO] all interactions in FV : " << sum_weight << std::endl;
    std::cout << "[INFO] numu CC                : " << sum_weight_cc
              << "  (" << 100.0 * sum_weight_cc / sum_weight << "% of interactions)" << std::endl;
    std::cout << "[INFO] oscillated             : " << sum_weight_osc
              << "  (ratio to unoscillated " << sum_weight_osc / sum_weight << ")"
              << std::endl;
    std::cout << "[INFO] data events per generated event : "
              << sum_weight_osc / nentries << std::endl;

    return 1;
}

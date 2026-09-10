#include <string>
#include <vector>

#include "TChain.h"
#include "TClass.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"

// Chains every fakerecon file of one beam mode, assigns each event to a sample,
// and writes a slim tree with the surviving events and their reconstructed
// kinematics. The inputs have different POT and fluxes behind them, but
// POTWeight already carries each sample's own normalisation, so chaining them
// is safe and nothing below needs to know which file an entry came from -
// except the appearance/intrinsic split, which does, for the reason given
// there.
//
// The assignment is a draw, not a cut. response.root holds, for each
// interaction category and each true neutrino energy, the probability that a
// surviving event ends up 1-ring mu-like, 1-ring e-like with no decay electron,
// or 1-ring e-like with one; whatever is left over is rejected. FakeRecon has
// already applied the FCFV and single-ring efficiency from the same file as a
// weight, so between them the two stages are the entire selection.
//
// This replaces a chain of cuts on invented PID likelihoods. A mixing
// probability is what such a cut produces once it has been integrated over
// everything except energy, so modelling the likelihood in order to threshold
// it was a long way round to the same number, with the thresholds and the
// template shapes both free and only their product observable.
//
// isRHC is an argument rather than a local so both modes can be run without
// editing the macro.
// input_suffix picks which stage of the chain to read. "_spline.root" is the
// output of MakeSplines.C and carries the cross-section response graphs;
// "_fakerecon.root" is the stage before that and has none. Either works - the
// graph branches are discovered rather than assumed, so dropping back to
// fakerecon just produces an output tree with no cross-section dials in it.
int FakeSelect (bool isRHC,
                std::string input_dir,
                std::string output_file,
                std::string response_file = "lbl_config/hk/response.root",
                UInt_t rng_seed = 20250731,
                std::string input_suffix = "_spline.root") {
    // The three sample types of arXiv:2505.15019 Tab. 2. The beam mode is
    // already fixed by isRHC, so it is not part of the enum; 1Re1De exists in
    // FHC only, as in the paper.
    enum SampleType {k1Rmu, k1Re, k1Re1De, kNSample};
    const int kNotSelected = -2;

    // Interaction categories, in the order MakeResponse.C writes them and
    // FakeRecon indexes them. Nothing here recomputes the category from truth:
    // it is read from the RespCat branch, so the efficiency that was applied
    // upstream and the mixing applied here cannot disagree about what an event
    // is.
    const int kNCat = 9;
    const char *cat_name[kNCat] = {"numuCCQE", "numubarCCQE",
                                   "numuCCnonQE", "numubarCCnonQE",
                                   "nueCCQE", "nuebarCCQE",
                                   "nueCCnonQE", "nuebarCCnonQE", "NC"};
    const char *mix_name[kNSample] = {"1Rmu", "1Re0de", "1Re1de"};

    // Every sample of a given mode, intrinsic and appearance alike. The
    // <flux>_x_<interaction> name is what makes a composition table readable, so
    // the lists are explicit rather than globbed. Both lists are in the same
    // order, one file per Tab. 2 column, which is what col_of_file relies on.
    std::vector<std::string> FHC_list = {
        input_dir + "/numu/numu_x_numu_200K_FHC",
        input_dir + "/numubar/numubar_x_numubar_50K_FHC",
        input_dir + "/nue/nue_x_nue_200K_FHC",
        input_dir + "/nuebar/nuebar_x_nuebar_50K_FHC",
        input_dir + "/numu/numu_x_nue_50K_FHC",          // appearance numu->nue
        input_dir + "/numubar/numubar_x_nuebar_50K_FHC"  // appearance numub->nueb
    };
    std::vector<std::string> RHC_list = {
        input_dir + "/numu/numu_x_numu_50K_RHC",
        input_dir + "/numubar/numubar_x_numubar_200K_RHC",
        input_dir + "/nue/nue_x_nue_50K_RHC",
        input_dir + "/nuebar/nuebar_x_nuebar_200K_RHC",
        input_dir + "/numu/numu_x_nue_50K_RHC",          // appearance numu->nue
        input_dir + "/numubar/numubar_x_nuebar_50K_RHC"  // appearance numub->nueb
    };

    // CCQE energy reconstruction constants. The binding energy is the T2K value
    // for oxygen.
    const double m_n = 0.93957;   // GeV
    const double m_p = 0.93827;
    const double E_b = 0.027;
    const double m_mu = 0.105658;
    const double m_e  = 0.000511;

    const std::string mode = isRHC ? "RHC" : "FHC";
    std::vector<std::string> &file_list = isRHC ? RHC_list : FHC_list;


    const int kNCol = 6;
    const char *col_name[kNCol] = {"beam numu", "beam nue", "beam numub",
                                   "beam nueb", "numu->nue", "numub->nueb"};
    // file index -> column. Both lists are ordered the same way, so this is
    // the same map for either beam mode.
    const int col_of_file[6] = {0, 2, 1, 3, 4, 5};


    const int flux_pdg_of_file[6] = {14, -14, 12, -12, 14, -14};

    // ------------------------------------------------------------- response
    TFile *FileResp = TFile::Open(response_file.c_str(), "read");
    if (!FileResp || FileResp->IsZombie()) {
        std::cout << "[ERROR] could not open response file " << response_file
                  << " - run MakeResponse.C first" << std::endl;
        return 0;
    }
    TH1D *h_mix[kNCat][kNSample];
    for (int c = 0; c < kNCat; ++c) {
        for (int s = 0; s < kNSample; ++s) {
            h_mix[c][s] = (TH1D*)FileResp->Get(Form("h_mix_%s_%s", cat_name[c], mix_name[s]));
            if (!h_mix[c][s]) {
                std::cout << "[ERROR] " << response_file << " has no h_mix_"
                          << cat_name[c] << "_" << mix_name[s] << std::endl;
                return 0;
            }
        }
    }

    // ----------------------------------------------------------------- chain
    TChain events("FlatTree_VARS");
    for (size_t i = 0; i < file_list.size(); ++i) {
        file_list[i] += input_suffix;
        if (gSystem->AccessPathName(file_list[i].c_str())) {
            std::cout << "[ERROR] missing input " << file_list[i]
                      << " - run FakeRecon, then MakeSplines, on every " << mode
                      << " sample first" << std::endl;
            return 0;
        }
        events.Add(file_list[i].c_str());
    }

    const Long64_t nentries = events.GetEntries();
    std::cout << "[INFO] " << mode << ": " << file_list.size() << " files, "
              << nentries << " entries" << std::endl;

    // Every branch FakeRecon writes must be present. Without this check,
    // SetBranchAddress on a missing branch leaves the variable uninitialised
    // and every weight silently comes out zero - which is exactly what happens
    // if FakeRecon has not been rerun after its branch list changed.
    const char *needed[] = {"POTWeight", "OscWeight", "SelEffWeight", "PLep_rec",
                            "CosLep_rec", "RespCat", "recon_ME", "Enu_true",
                            "PDGnu", "cc", "flagCCQE"};
    for (const char *b : needed) {
        if (!events.GetBranch(b)) {
            std::cout << "[ERROR] input has no branch '" << b
                      << "' - rerun FakeRecon (bash recon.sh)" << std::endl;
            return 0;
        }
    }

    double POTWeight, OscWeight, SelEffWeight;
    float  PLep_rec, CosLep_rec, Enu_true;
    int    PDGnu, recon_ME, RespCat;
    Char_t cc;
    Bool_t flagCCQE;
    events.SetBranchAddress("POTWeight",        &POTWeight);
    events.SetBranchAddress("OscWeight",        &OscWeight);
    events.SetBranchAddress("SelEffWeight",     &SelEffWeight);
    events.SetBranchAddress("PLep_rec",         &PLep_rec);
    events.SetBranchAddress("CosLep_rec",       &CosLep_rec);
    events.SetBranchAddress("RespCat",          &RespCat);
    events.SetBranchAddress("recon_ME",         &recon_ME);
    events.SetBranchAddress("Enu_true",         &Enu_true);
    events.SetBranchAddress("PDGnu",            &PDGnu);
    events.SetBranchAddress("cc",               &cc);
    events.SetBranchAddress("flagCCQE",         &flagCCQE);

    TRandom3 rng(rng_seed);

    // ---------------------------------------------------------------- output
    // A slim tree rather than a clone of the chain: downstream needs the
    // EventType, the weights and the recon kinematics, plus truth for breakdowns.
    TFile *FileOut = TFile::Open(output_file.c_str(), "recreate");
    FileOut->cd();
    TTree *outtree = new TTree("FakeSelect", "sample-assigned fake data");

    int   EventTypeOut, SampleIndex, PDGLep_rec, flux_pdg;
    float Theta_rec;
    float Enu_rec, Q2_rec;
    // Constant for the whole file, but written per event anyway: the two beam
    // modes get chained together downstream and nothing in a chained tree
    // remembers which file it came from. GUNDAM needs it to split the samples
    // and to pick the right flux and detector parameters.
    Bool_t isRHC_out = isRHC;
    outtree->Branch("isRHC",            &isRHC_out,        "isRHC/O");
    outtree->Branch("EventType",        &EventTypeOut,     "EventType/I");
    outtree->Branch("SampleIndex",      &SampleIndex,      "SampleIndex/I");
    outtree->Branch("flux_pdg",         &flux_pdg,         "flux_pdg/I");
    outtree->Branch("RespCat",          &RespCat,          "RespCat/I");
    outtree->Branch("PDGLep_rec",       &PDGLep_rec,       "PDGLep_rec/I");
    outtree->Branch("recon_ME",         &recon_ME,         "recon_ME/I");
    outtree->Branch("Enu_rec",          &Enu_rec,          "Enu_rec/F");
    outtree->Branch("Q2_rec",           &Q2_rec,           "Q2_rec/F");
    outtree->Branch("POTWeight",        &POTWeight,        "POTWeight/D");
    outtree->Branch("OscWeight",        &OscWeight,        "OscWeight/D");
    outtree->Branch("SelEffWeight",     &SelEffWeight,     "SelEffWeight/D");
    outtree->Branch("PLep_rec",         &PLep_rec,         "PLep_rec/F");
    outtree->Branch("CosLep_rec",       &CosLep_rec,       "CosLep_rec/F");
    // The same angle in degrees. T2K bins its samples in the lepton angle to
    // the beam, not its cosine, and the binning files are copied from T2K
    // rather than rewritten - so the branch is provided in the units they use.
    outtree->Branch("Theta_rec",        &Theta_rec,        "Theta_rec/F");
    outtree->Branch("Enu_true",         &Enu_true,         "Enu_true/F");
    outtree->Branch("PDGnu",            &PDGnu,            "PDGnu/I");
    outtree->Branch("cc",               &cc,               "cc/B");

    // ------------------------------------------- cross-section response graphs
    // Written by MakeSplines.C, one TGraph branch per T2KReWeight dial, holding
    // this event's weight as a function of that dial. They have to be carried
    // through by hand: this output tree is built branch by branch rather than
    // cloned, so nothing propagates on its own.
    //
    // Discovered by class rather than from a list of dial names, so adding or
    // removing a dial in the T2KReWeight config needs no edit here. Anything
    // that is a TGraph is a response graph - nothing else in the chain is one.
    //
    // Split level 0 on both sides. GUNDAM reads the leaf as a TObject* and asks
    // it its class name (DataDispenser.cpp:1344), so the graph has to arrive as
    // one object rather than as scattered member branches.
    std::vector<std::string> graph_name;
    std::vector<TGraph*>     graph_buf;
    TObjArray *in_branches = events.GetListOfBranches();
    for (int ib = 0; ib < (in_branches ? in_branches->GetEntries() : 0); ++ib) {
        TBranch *br = (TBranch*)in_branches->At(ib);
        TClass  *cl = TClass::GetClass(br->GetClassName());
        if (!cl || !cl->InheritsFrom(TGraph::Class())) continue;
        graph_name.emplace_back(br->GetName());
    }
    graph_buf.assign(graph_name.size(), nullptr);
    for (size_t g = 0; g < graph_name.size(); ++g) {
        events.SetBranchAddress(graph_name[g].c_str(), &graph_buf[g]);
        outtree->Branch(graph_name[g].c_str(), "TGraph", &graph_buf[g], 32000, 0);
    }
    if (graph_name.empty()) {
        std::cout << "[INFO] no response graphs in the input - the output will "
                     "carry no cross-section dials" << std::endl;
    } else {
        std::cout << "[INFO] carrying " << graph_name.size()
                  << " cross-section response graphs:";
        for (auto &n : graph_name) std::cout << " " << n;
        std::cout << std::endl;
    }

    // ------------------------------------------------------- Tab. 2 accumulator
    // [sample][column], exactly the layout of arXiv:2505.15019 Tab. 2.
    const char *sample_name[kNSample] = {"1-ring mu-like",
                                         "1-ring e-like + 0 decay e",
                                         "1-ring e-like + 1 decay e"};
    std::vector<std::vector<double>> tab2(kNSample, std::vector<double>(kNCol, 0.0));
    std::vector<std::vector<double>> bycat(kNSample, std::vector<double>(kNCat, 0.0));
    std::vector<Long64_t> count(kNSample, 0);

    for (Long64_t i = 0; i < nentries; ++i) {
        events.GetEntry(i);

        // Which input file this entry came from. Needed because truth cannot
        // separate an appeared nue from an intrinsic beam nue - both are CC with
        // PDGnu 12 - so the appearance signal is only identifiable by its flux,
        // i.e. by which sample it was read out of.
        SampleIndex = events.GetTreeNumber();
        const int col = col_of_file[SampleIndex];
        flux_pdg = flux_pdg_of_file[SampleIndex];

        if (RespCat < 0 || RespCat >= kNCat) {
            std::cout << "[ERROR] entry " << i << " has RespCat " << RespCat
                      << ", outside 0.." << kNCat - 1 << std::endl;
            return 0;
        }

        // ------------------------------------------------- the sample draw
        // One uniform per event against the running sum of the three
        // probabilities, so the three samples are mutually exclusive by
        // construction and the leftover is the rejected fraction. The draw is
        // made whatever the weight, and the weight is not consulted here: an
        // event with a tiny SelEffWeight is still assigned a sample, it just
        // arrives carrying almost nothing.
        int sample = kNotSelected;
        {
            TH1D *h0 = h_mix[RespCat][0];
            const int ie = std::min(std::max(h0->FindBin(Enu_true), 1), h0->GetNbinsX());
            const double r = rng.Uniform();
            double acc = 0.0;
            for (int s = 0; s < kNSample; ++s) {
                acc += h_mix[RespCat][s]->GetBinContent(ie);
                if (r < acc) { sample = s; break; }
            }
            // 1Re1De is a nu-mode sample only, as in the paper. In nubar mode a
            // draw that lands there is a rejection instead.
            if (isRHC && sample == k1Re1De) sample = kNotSelected;
        }

        // ------------------------------------------------- recon kinematics
        // Built from the *assigned* sample, never from PDGLep. An NC event drawn
        // as mu-like gets the muon mass and the muon-neutrino CCQE formula,
        // which is exactly the mis-reconstruction the sample is supposed to
        // contain. Done for rejected events too, so we can still look at what
        // the selection threw away.
        PDGLep_rec = (sample == k1Rmu) ? (isRHC ? -13 : 13)
                                       : (isRHC ? -11 : 11);
        {
            const double m_lep = (abs(PDGLep_rec) == 13) ? m_mu : m_e;
            const double E_lep = std::sqrt(PLep_rec * PLep_rec + m_lep * m_lep);

            // CCQE quasi-elastic assumption. In neutrino mode the struck
            // nucleon is a neutron and the outgoing one a proton; in
            // antineutrino mode they swap. Which one to assume is a beam-mode
            // choice, not a truth-level one, so isRHC is legitimate here.
            const double m_tgt = isRHC ? m_p : m_n;   // struck
            const double m_out = isRHC ? m_n : m_p;   // outgoing
            const double m_eff = m_tgt - E_b;

            const double num = m_out * m_out - m_eff * m_eff - m_lep * m_lep
                             + 2.0 * m_eff * E_lep;
            const double den = 2.0 * (m_eff - E_lep + PLep_rec * CosLep_rec);

            // The denominator vanishes for very forward, high-momentum leptons,
            // where the QE assumption has no solution. Flag rather than let an
            // infinity into the tree.
            if (den != 0.0 && num / den > 0.0) {
                Enu_rec = num / den;
                Q2_rec  = -m_lep * m_lep
                        + 2.0 * Enu_rec * (E_lep - PLep_rec * CosLep_rec);
            } else {
                Enu_rec = -1.0;
                Q2_rec  = -1.0;
            }
        }

        // The only place the beam mode enters the sample index. FakeSelect still
        // works in the three mixing outcomes throughout - the mixing matrix does
        // not know about beam modes and should not - but the value written out is
        // offset so the five Tab. 2 rows are five distinct numbers: nu mode keeps
        // 0, 1, 2 and nubar mode becomes 3, 4. A GUNDAM sample is then one
        // equality on one branch instead of a mode cut and a type cut.
        Theta_rec = TMath::RadToDeg() * std::acos(std::min(1.0f, std::max(-1.0f, CosLep_rec)));

        EventTypeOut = (sample == kNotSelected) ? kNotSelected
                     : (isRHC ? sample + kNSample : sample);
        if (sample != kNotSelected) {
            const double w = POTWeight * OscWeight * SelEffWeight;
            tab2[sample][col]      += w;
            bycat[sample][RespCat] += w;
            ++count[sample];
        }

        outtree->Fill();

        if (i % 100000 == 0)
            std::cout << "[INFO] entry " << i << " / " << nentries << std::endl;
    }

    outtree->Write("", TObject::kOverwrite);
    FileOut->Close();
    FileResp->Close();

    // --------------------------------------------------------------- summary
    printf("\n%-38s", "");
    for (int c2 = 0; c2 < kNCol; ++c2) printf("%13s", col_name[c2]);
    printf("%13s%12s\n", "Total", "(entries)");

    for (int t = 0; t < kNSample; ++t) {
        if (isRHC && t == k1Re1De) continue;   // FHC only, as in the paper
        double tot = 0.0;
        char row[128];
        snprintf(row, sizeof(row), "%s-mode, %s", isRHC ? "nubar" : "nu", sample_name[t]);
        printf("%-38s", row);
        for (int c2 = 0; c2 < kNCol; ++c2) { printf("%13.1f", tab2[t][c2]); tot += tab2[t][c2]; }
        printf("%13.1f%12lld\n", tot, count[t]);
    }

    // The same rates split by interaction category rather than by parent
    // flavour. This is the view that says which response curve to move when a
    // cell is wrong, which the column view cannot: a column mixes CC and NC.
    printf("\n%-38s", "by interaction category");
    for (int c = 0; c < kNCat; ++c) printf("%16s", cat_name[c]);
    printf("\n");
    for (int t = 0; t < kNSample; ++t) {
        if (isRHC && t == k1Re1De) continue;
        char row[128];
        snprintf(row, sizeof(row), "%s-mode, %s", isRHC ? "nubar" : "nu", sample_name[t]);
        printf("%-38s", row);
        for (int c = 0; c < kNCat; ++c) printf("%16.1f", bycat[t][c]);
        printf("\n");
    }

    std::cout << "\n[INFO] wrote " << output_file << std::endl;
    return 1;
}

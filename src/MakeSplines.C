#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "TFile.h"
#include "TGraph.h"
#include "TObjArray.h"
#include "TTree.h"

// Cross-section response splines for GUNDAM, from T2KReWeight output.
//
// Reads a _fakerecon.root and its matching _Weight.flat.root and writes a copy
// of the fakerecon tree with one TGraph branch per cross-section dial. Each
// graph is that event's weight as a function of that dial's value.
//
// ---------------------------------------------------------------------------
// Why a TGraph per event and not a table of numbers
//
// GUNDAM reads event-by-event dials through dialLeafName (aliased from
// treeExpression, DialCollection.cpp:41). DataDispenser.cpp:1339 looks the leaf
// up in the dataset tree, takes it as a TObject* and hands it to
// DialCollection::makeDial, which ends in DialUtils::getPointList
// (DialUtils.cpp:14-19):
//
//   if( src_->InheritsFrom(TGraph::Class()) )   { ... }
//   if( src_->InheritsFrom(TSpline3::Class()) ) { ... }
//   return {};
//
// Anything else returns an empty point list, which is silently treated as "no
// dial" rather than as an error. So a pair of vector<double> branches would
// load without complaint and do nothing at all. It has to be a TGraph or a
// TSpline3, and it has to live in the tree GUNDAM reads events from - there is
// no AddFriend anywhere in the dataset loader, so a sidecar file is not an
// option.
//
// The knots do not need sorting here: checkDialPointList sorts by x
// (DialCollection.cpp:1233). They are sorted anyway, because a graph that reads
// in ascending order is easier to check by eye.
//
// ---------------------------------------------------------------------------
// The nominal has to be a knot, and T2KReWeight does not provide it
//
// TweakMetaData stores each dial's tweak points and <dial>_from, the value the
// reweighting started from. For most dials the nominal is not among the tweak
// points - MaCCQE has from = 1.21 and tweaks at 0.818, 1.014, 1.406, 1.602 -
// so the weight at nominal is never stored. It is 1.0 by construction, and this
// macro inserts it.
//
// Without that insertion the interpolated response at the prior would not be
// exactly 1, so the pre-fit MC would sit slightly off nominal before the fit
// even starts. Where from does coincide with a tweak point (RadCorr_rate has
// from = 0 and a tweak at 0) no point is inserted, and the stored weight there
// is checked against 1.0 instead - the printout reports the largest deviation
// found, which is a direct test of the assumption.
//
// ---------------------------------------------------------------------------
// Flat dials are written as empty graphs
//
// Most dials do not apply to most events: MaCCQE leaves about two thirds of
// them untouched, since it only moves CCQE. Those arrive as a vector of exact
// 1.0s. An empty TGraph gives an empty point list, which makeDialShortCircuit
// (DialCollection.cpp:1247) turns into no dial at all - the same outcome as a
// flat graph, at a fraction of the file size.
//
// The test is "every weight is exactly 1.0", not "the weights are all equal",
// and the difference matters. A response that is flat at some other value is a
// real constant rescaling, and GUNDAM represents it with a Shift dial
// (DialCollection.cpp:1259). Those graphs are written out in full.
//
// ---------------------------------------------------------------------------
// The join is by row position, and that is only safe because both files drop
// the same events
//
// Neither tree carries an event ID. T2KReWeight drops Pauli blocked events
// (NIWGFlatTree.h:470, EventWeightGen.cxx:396-399) and FakeRecon drops the same
// ones on Q2 = -999, so entry i of one is entry i of the other. The counts are
// compared below and a mismatch is fatal, because nothing downstream would
// notice events carrying each other's cross-section weights.

int MakeSplines (std::string recon_file,
                 std::string weight_file,
                 std::string output_file) {

    // ------------------------------------------------------------- open input
    TFile *FileRec = TFile::Open(recon_file.c_str(), "read");
    if (!FileRec || FileRec->IsZombie()) {
        std::cout << "[ERROR] could not open " << recon_file << std::endl;
        return 0;
    }
    TTree *rectree = (TTree*)FileRec->Get("FlatTree_VARS");
    if (!rectree) {
        std::cout << "[ERROR] no FlatTree_VARS in " << recon_file << std::endl;
        return 0;
    }

    TFile *FileWgt = TFile::Open(weight_file.c_str(), "read");
    if (!FileWgt || FileWgt->IsZombie()) {
        std::cout << "[ERROR] could not open " << weight_file << std::endl;
        return 0;
    }
    TTree *wgttree = (TTree*)FileWgt->Get("T2KRWFlatTree");
    TTree *metatree = (TTree*)FileWgt->Get("TweakMetaData");
    if (!wgttree || !metatree) {
        std::cout << "[ERROR] " << weight_file
                  << " needs both T2KRWFlatTree and TweakMetaData" << std::endl;
        return 0;
    }

    const Long64_t nentries = rectree->GetEntries();
    if (nentries != wgttree->GetEntries()) {
        std::cout << "[ERROR] " << nentries << " entries in " << recon_file
                  << " but " << wgttree->GetEntries() << " in " << weight_file
                  << " - these are joined by row position, so this must match "
                     "exactly. Rerun FakeRecon and T2KReWeight from the same "
                     "generated file." << std::endl;
        return 0;
    }

    // ------------------------------------------------------------- dial setup
    // Knot positions are the same for every event, so the ordering is worked out
    // once here. src[k] is where the weight for sorted knot k comes from: an
    // index into the event's weight vector, or -1 meaning the inserted nominal.
    struct Dial {
        std::string          name;
        std::vector<double>  x;
        std::vector<int>     src;
        std::vector<double> *wbuf = nullptr;
        double               nominal = 0.0;
        Long64_t             n_flat = 0;   // no response, legitimately
        Long64_t             n_bad  = 0;   // response rejected, see below
    };
    std::vector<Dial> dial;

    TObjArray *meta_branches = metatree->GetListOfBranches();
    for (int ib = 0; ib < meta_branches->GetEntries(); ++ib) {
        std::string bname = meta_branches->At(ib)->GetName();
        if (bname.size() > 5 && bname.compare(bname.size() - 5, 5, "_from") == 0) continue;
        if (!metatree->GetBranch((bname + "_from").c_str())) continue;

        std::vector<double> *tweaks = nullptr;
        double from = 0.0;
        metatree->SetBranchAddress(bname.c_str(), &tweaks);
        metatree->SetBranchAddress((bname + "_from").c_str(), &from);
        metatree->GetEntry(0);
        if (!tweaks || tweaks->empty()) {
            metatree->ResetBranchAddresses();
            continue;
        }

        // A dial name that collides with a branch already in the recon tree
        // would make the output tree ambiguous, and GUNDAM resolves leaves by
        // name. Refuse rather than write something that reads back wrong.
        if (rectree->GetBranch(bname.c_str())) {
            std::cout << "[ERROR] dial " << bname << " collides with an existing "
                      << "branch in " << recon_file << std::endl;
            return 0;
        }

        Dial d;
        d.name    = bname;
        d.nominal = from;

        // (position, source) pairs, nominal included, then sorted by position.
        std::vector<std::pair<double, int>> pt;
        pt.reserve(tweaks->size() + 1);
        for (size_t k = 0; k < tweaks->size(); ++k) pt.emplace_back((*tweaks)[k], (int)k);

        bool from_is_tweak = false;
        for (size_t k = 0; k < tweaks->size(); ++k)
            if (std::abs((*tweaks)[k] - from) < 1e-9 * (1.0 + std::abs(from)))
                from_is_tweak = true;
        if (!from_is_tweak) pt.emplace_back(from, -1);

        std::sort(pt.begin(), pt.end(),
                  [](const std::pair<double,int> &a, const std::pair<double,int> &b) {
                      return a.first < b.first;
                  });
        for (auto &p : pt) { d.x.push_back(p.first); d.src.push_back(p.second); }

        dial.push_back(d);
        metatree->ResetBranchAddresses();
    }

    if (dial.empty()) {
        std::cout << "[ERROR] no dials found in TweakMetaData" << std::endl;
        return 0;
    }

    for (auto &d : dial) wgttree->SetBranchAddress(d.name.c_str(), &d.wbuf);

    std::cout << "[INFO] " << dial.size() << " dials, " << nentries << " entries"
              << std::endl;
    for (auto &d : dial) {
        printf("[INFO]   %-24s nominal %8.4f   %zu knots:", d.name.c_str(),
               d.nominal, d.x.size());
        for (size_t k = 0; k < d.x.size(); ++k)
            printf(" %g%s", d.x[k], d.src[k] < 0 ? "*" : "");
        printf("\n");
    }
    std::cout << "[INFO] * marks the inserted nominal, where the weight is 1 "
                 "by construction" << std::endl;

    // ----------------------------------------------------------- output tree
    TFile *FileOut = TFile::Open(output_file.c_str(), "recreate");
    if (!FileOut || FileOut->IsZombie()) {
        std::cout << "[ERROR] could not write " << output_file << std::endl;
        return 0;
    }
    FileOut->cd();
    TTree *outtree = rectree->CloneTree(0);

    // Split level 0: the graph has to reach GUNDAM as one object, because
    // DataDispenser reads the leaf as a TObject* and asks it what class it is.
    std::vector<TGraph*> graph(dial.size(), nullptr);
    for (size_t d = 0; d < dial.size(); ++d) {
        graph[d] = new TGraph();
        outtree->Branch(dial[d].name.c_str(), "TGraph", &graph[d], 32000, 0);
    }

    // ------------------------------------------------------------ event loop
    double worst_nominal = 0.0;
    std::string worst_dial;
    for (Long64_t i = 0; i < nentries; ++i) {
        rectree->GetEntry(i);
        wgttree->GetEntry(i);

        for (size_t d = 0; d < dial.size(); ++d) {
            Dial &dd = dial[d];
            const size_t npt = dd.x.size();

            bool all_one  = true;   // no response at all
            bool all_zero = true;   // response is dead, see below
            bool bad      = false;  // response is unusable, see below

            for (size_t k = 0; k < npt; ++k) {
                if (dd.src[k] < 0) continue;          // the inserted nominal
                const double y = (*dd.wbuf)[dd.src[k]];

                if (y != 1.0) all_one = false;
                if (y != 0.0) all_zero = false;

                // A negative or non-finite weight is not a response, it is a
                // failed calculation.
                if (!std::isfinite(y) || y < 0.0) bad = true;

                // Where the nominal happens to coincide with a tweak point, the
                // stored weight there must be 1: reweighting a dial from its
                // own value to itself changes nothing, by construction. If it
                // is not 1 the engine did not produce a usable response for
                // this event, whatever the other knots say.
                if (std::abs(dd.x[k] - dd.nominal)
                        < 1e-9 * (1.0 + std::abs(dd.nominal))) {
                    const double dev = std::abs(y - 1.0);
                    if (dev > worst_nominal) { worst_nominal = dev; worst_dial = dd.name; }
                    if (dev > 1e-6) bad = true;
                }
            }

            // Every stored knot is zero. This happens for whole events at a
            // time - all dials at once, and cv_weight is zero for them too -
            // because T2KReWeight::CalcWeight multiplies over its weight
            // engines and one of them rejects the event at nominal. Which
            // engine was not tracked down; the NEUT engines all return 1 when
            // untweaked, so it is on the NIWG side.
            //
            // These events are given weight 1, which is what an empty graph
            // means: no dial, no response, event unchanged. The one thing that
            // must not happen is passing the zeros through, because a graph
            // that is flat at zero reaches makeDialShortCircuit
            // (DialCollection.cpp:1249-1263), is read as a genuine constant
            // rescaling, and becomes a Shift dial of value 0 - which deletes
            // the event from the fit without saying so.
            if (all_zero) bad = true;

            if (bad) {
                graph[d]->Set(0);
                ++dd.n_bad;
                continue;
            }

            if (all_one) {
                // No dial for this event. An empty graph is how that is said.
                graph[d]->Set(0);
                ++dd.n_flat;
                continue;
            }

            graph[d]->Set((Int_t)npt);
            for (size_t k = 0; k < npt; ++k) {
                const double y = (dd.src[k] < 0) ? 1.0 : (*dd.wbuf)[dd.src[k]];
                graph[d]->SetPoint((Int_t)k, dd.x[k], y);
            }
        }

        outtree->Fill();
        if (i % 50000 == 0)
            std::cout << "[INFO] entry " << i << " / " << nentries << std::endl;
    }

    FileOut->cd();
    outtree->Write("", TObject::kOverwrite);
    FileOut->Close();
    FileRec->Close();
    FileWgt->Close();

    // -------------------------------------------------------------- printout
    std::cout << "[INFO] wrote " << nentries << " entries to " << output_file
              << std::endl;
    // "no response" is the normal case and costs nothing. "rejected" is not
    // normal: those events keep no dial for that systematic, so it is a real
    // hole in the coverage and the number needs to be looked at rather than
    // buried.
    std::cout << "[INFO] per dial: events with no response / events rejected"
              << std::endl;
    Long64_t total_bad = 0;
    for (auto &d : dial) {
        total_bad += d.n_bad;
        printf("[INFO]   %-24s no response %8lld (%5.1f%%)   rejected %6lld (%.3f%%)%s\n",
               d.name.c_str(), d.n_flat, 100.0 * d.n_flat / nentries,
               d.n_bad, 100.0 * d.n_bad / nentries, d.n_bad ? "  <--" : "");
    }

    if (!worst_dial.empty())
        printf("[INFO] largest |weight - 1| at a nominal that is also a tweak "
               "point: %.3g (%s)\n", worst_nominal, worst_dial.c_str());
    else
        std::cout << "[INFO] no dial has its nominal among the tweak points"
                  << std::endl;

    if (total_bad > 0)
        std::cout << "[INFO] " << total_bad << " dial-event responses were "
                     "given weight 1 because T2KReWeight returned zero at every "
                     "knot, including the nominal. That happens for whole events "
                     "rather than for single dials." << std::endl;

    return 1;
}

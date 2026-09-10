#include <string>

#include "TFile.h"
#include "TH1D.h"

// Writes the true-energy grid that the tabulated oscillation dial fills.
//
// NuFast recomputes probabilities onto this grid every time a PMNS parameter
// moves, and each event then looks its own probability up by true energy. The
// grid is therefore the resolution of the oscillation calculation, and it has to
// be fine where the probability moves fastest - at Hyper-K that is the first
// oscillation maximum at 600 MeV, where P(numu->numu) swings from 1 to nearly 0
// across a few hundred MeV.
//
// The DUNE grid shipped with gudi-lbl-tdr25a runs 0.1 to 100 GeV, which suits a
// 1300 km baseline and leaves nothing below 100 MeV. The T2K flux starts at
// 10 MeV, so this one starts at zero instead and stops at 10 GeV, above which
// the T2K beam has no events worth the bins.
int MakeOscGrid (std::string output_file,
                 int    nbins = 2000,
                 double emin  = 0.0,
                 double emax  = 10.0) {

    TFile *FileOut = TFile::Open(output_file.c_str(), "recreate");
    TH1D *h = new TH1D("energyBinning", "oscillation grid;E_{#nu} [GeV];", nbins, emin, emax);
    h->Write();
    FileOut->Close();

    printf("[INFO] wrote %s: %d bins of %.1f MeV from %.2f to %.2f GeV\n",
           output_file.c_str(), nbins, 1000.0 * (emax - emin) / nbins, emin, emax);
    return 1;
}

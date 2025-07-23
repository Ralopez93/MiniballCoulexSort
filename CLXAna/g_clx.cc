
// Sorting for Miniball Coulex data and building histograms
// Cuts are done in doppler.hh using function, doppler::Cut()
// Liam Gaffney (liam.gaffney@cern.ch) - 26/10/2016

#define g_clx_cxx

// Select settings by uncommenting one of the following:
// #define TWOPART 		// define to plot every 2p combination
// #define GEANG			// define for plotting Ge angles per
// cluster #define SPEDEGEOMETRY	// define to overwrite Spede angles with
// SpedeGeometry routine

#ifndef __SpedeGeometry_HH__
#include "SpedeGeometry.hh"
#endif
#ifndef hist_hh
#include "hists.hh"
#endif
#ifndef g_clx_hh
#include "g_clx.hh"
#endif

void g_clx::Loop(string outputfilename) {
  if (fChain == 0)
    return;

  // Output file name
  TFile *out = new TFile(outputfilename.c_str(), "RECREATE");

  // Create doppler instance and set experimental parameters
  doppler dc;
  dc.ExpDefs(Zb, Ab, Zt, At, Eb, Ex, thick, depth, cddist, cdoffset, deadlayer,
             contaminant, spededist, Bcut, Tcut, srim, usekin, calfile);
  dc.mbAngles(); // re-define MB angles
  // Create stopping power curves from the srim output files
  // Comment out to use the default parameters in doppler.hh
  // stoppingpowers( BT, TT, BA, TA, BC, TC )
  if (!dc.stoppingpowers(true, true, true, true, false, false)) {
    cout << "Definition of stopping powers failed" << endl;
    return;
  }
  // Test if it's an electron or gamma. Note: currently unused?
  bool electron;

  // Include errors on histograms (required for correct bg subtraction)
  TH1::SetDefaultSumw2();

  // Declare the histograms here and initialise!
  hists h;
  cout << "initializing hists" << endl;
  h.Initialise(dc);
  cout << "Hists initialized" << endl;

  // Particle-particle time difference (from tppdiff)
  h.Set_ppwin(100.);

  // Maximum strip number where recoils and projectiles are separable
  // Applicable mostly for inverse kinematics. To turn off, set to 16
  h.Set_maxrecoil(16);

  // Minimum strip number to define high CoM recoils
  // For normal kinematics, this is likely to be the limit of the safe condition
  // i.e. strips that are >= minrecoil are unsafe. Only low CoM solution is safe
  h.Set_minrecoil(0);

  // Loop over events
  cout << "Looping over events...\n";
  Int_t nbytes = 0, nbs = 0;
  Int_t skipFactor = 1;
  for (Long64_t jentry = 0; jentry < fChain->GetEntries() / skipFactor; jentry++) {

    Long64_t ientry = LoadTree(jentry);

    if (ientry < 0)
      break;

    nbs = fChain->GetEntry(jentry);
    nbytes += nbs;

    if (jentry % 5000 == 0) {
      cout << " " << jentry << "/" << fChain->GetEntries() << "  ("
           << jentry * 100. / fChain->GetEntries() << "%)    \r";
      cout << flush;
    }
    if (clu_tune != -1 && cluid != clu_tune) {
      continue;
    }

    // Is it an electron or gamma? Note: unused?
    if (cluid < 8)
      electron = false;
    else if (cluid == 8)
      electron = true;
    else {
      cout << "Unexpected cluid: " << cluid << endl;
      break; // shouldn't be anything else
    }

    // fill CLX tree according to standard convention

    // adjust MB angles according to cid, sid
    tha = dc.GetGTh(cid, sid);
    pha = dc.GetGPh(cid, sid);
    for (int j = 0; j < gcor_gen.size(); j++) {
      gcor_tha[j] = dc.GetGTh(gcor_cid[j], gcor_sid[j]);
      gcor_pha[j] = dc.GetGPh(gcor_cid[j], gcor_sid[j]);
    }

    h.FillTree(gen, tha, pha, cluid, cid, sid, // single gamma
               gcor_gen, gcor_tha, gcor_pha, gcor_cluid, gcor_cid, gcor_sid, gcor_gtd, // correlated gamma
               laser, pen, nf, nb, sector, det, td, time); // particle info

  }

  out->Write();
}

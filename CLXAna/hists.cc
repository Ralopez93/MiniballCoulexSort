#ifndef hists_cxx
#define hists_cxx

#define GBINS 5000 // number of bins in gamma spectra
#define GMAX 5000  // maximum energy in gamma spectra
#define EBINS 1500 // number of bins in electron spectra
#define EMAX 3000  // maximum energy in electron spectra
#define PBINS 300  // number of bins in particle spectra
#define PMAX 1200  // maximum energy in particle spectra
#define TBINS 242  // number of bins in tdiff spectra
#define TMAX 1525  // maximum time in tdiff spectra

#define PTCUT_P 2.45e7 // Proton time cut limit for prompt
#define PTCUT_D 4.90e7 // Proton time cut limit for delayed
#define PD_FRAC -0.85  // scaling factor of delayed window for subtraction

#ifndef hist_hh
#include "hists.hh"
#endif

#define PID_BEAM 0
#define PID_TARG 1

void hists::Initialise(doppler dc_) {

  /// Initialise all of the histograms that you want to fill
  /// If you add a new one, be careful not to break your memory limit.
  /// This causes unknown behaviour and may not give any errors.
  /// In the best case, you have a segmentation violation and you cannot
  /// work out where it comes from... Try deleting some histograms.
  /// Whatever you do, make sure you declare them in hists.hh file.
  dc = dc_;

  cout << "Initialising tree and histograms...\n";

  string hname, htitle;

  // particle branches
  tree = new TTree("doppler", "doppler");
  tree->Branch("laser", &laser, "laser/I");
  tree->Branch("np", &np, "np/I"); //
  tree->Branch("tdpp", &tdpp, "tdpp/D"); // time difference for 2p events; 0 for 1p event
  tree->Branch("pid", pid, "pid[np]/I");    // particle ID: 0: target-like, 1: beam-like as defined by cuts
  tree->Branch("quad", quad, "quad[np]/I"); // quadrant ID
  tree->Branch("ring", ring, "ring[np]/I"); // ring ID
  tree->Branch("sect", sect, "sect[np]/I"); // sector ID

  tree->Branch("ep", ep, "ep[np]/D");    // calibrated particle energy in MeV
  tree->Branch("er", er, "er[np]/D");    // calibrated recoil particle energy in MeV (either reconstructed or 2p)
  tree->Branch("thp", thp, "thp[np]/D"); // theta of particle in lab frame
  tree->Branch("php", php, "php[np]/D"); // phi of particle in lab frame
  tree->Branch("thr", thr, "thr[np]/D"); // theta of recoil in lab frame
  tree->Branch("phr", phr, "phr[np]/D"); // phi of recoil in lab frame
  // gamma branches
  tree->Branch("ng", &ng, "ng/I");
  tree->Branch("td", td, "td[ng]/D"); // (smallest) time difference between particle and gamma
  tree->Branch("eg", eg, "eg[ng]/D"); // calibrated, not dc'ed gamma energy in keV
  tree->Branch("clu", clu, "clu[ng]/I"); // cluster ID from 0 to 7
  tree->Branch("cry", cry, "cry[ng]/I"); // crystal ID from 0 to 23
  tree->Branch("seg", seg, "seg[ng]/I"); // segment ID from 0 to 6 (0: core-only/ambiguous events)
  tree->Branch("thg", thg, "thg[ng]/D"); // theta of gamma in lab frame
  tree->Branch("phg", phg, "phg[ng]/D"); // phi of gamma in lab frame
  tree->Branch("tpg", tpg, "tpg[ng]/I"); // smallest particle-gamma time difference (if 2 particles)

  tree->Branch("abg", abg, "abg[ng]/D"); // angle between beam-like particle and gamma in degrees
  tree->Branch( "atg", atg, "atg[ng]/D"); // angle between target-like particle and gamma in degrees
  tree->Branch("ebg", ebg, "ebg[ng]/D"); // calibrated and dc-ed gamma energy to beam kinematics in keV
  tree->Branch("etg", etg, "etg[ng]/D"); // calibrated and dc-ed gamma energy to target kinematics in keV

  // Default values for g_clx definitions
  ppwin = 300.;
  maxrecoil = 16;
  minrecoil = 0;

  return;
}

void hists::Set_ppwin(float user_ppwin) { ppwin = user_ppwin; }

void hists::Set_maxrecoil(int user_maxrecoil) { maxrecoil = user_maxrecoil; }

void hists::Set_minrecoil(int user_minrecoil) { minrecoil = user_minrecoil; }

/**
 * @brief Fill the ROOT tree.
 * 
 * @param GEn Gamma energy [keV].
 * @param GTh Gamma lab theta angle [rad].
 * @param GPh Gamma lab phi angle [rad].
 * @param GCluid Gamma HPGE cluster ID [0-7]. 
 * @param GCid Gamma HPGE core ID [0-23].
 * @param GSid Gamma HPGE segment ID, where 0 is the core [0-6].
 * @param GCor_GEn Vector of correlated gamma-ray energies [keV].
 * @param GCor_GTh Vector of correlated theta angles {[rad]}.
 * @param GCor_GPh Vector of correlated phi angles {[rad]}.
 * @param GCor_GCluID Vector of correlated cluster IDs {[0-7]}.
 * @param GCor_GCid Vector of correlated core IDs {[0-23]}.
 * @param GCor_GSid Vector of correlated segment IDs {[0-6]}.
 * @param GCor_Gtd Vector of time-difference to original gamma-ray {[25ns]}
 * @param Laser Laser flag (on = 1, off = 0).
 * @param PEn Particle energies [keV].
 * @param Pnf Annular (front) strip ID of particle (0 = outer; 15 inner). A.K.A front CD Ring ID.
 * @param Pnb Secular (back) strip ID of particle (0 to 12; clockwise wrt beam).
 * @param Psec Sector of C-REX (0 = FCD; 1 = FBarrel; 2 = BBarrel; 3 = BCD).
 * @param Pquad Detector (quadrant) number of particle
 * @param Ptd Particle-gamma time difference in 25 ns timestamps.
 */
void hists::FillTree(float GEn, float GTh, float GPh, int GCluid, int GCid,
                     int GSid, vector<float> GCor_GEn, vector<float> GCor_GTh,
                     vector<float> GCor_GPh, vector<int> GCor_GCluID,
                     vector<int> GCor_GCid, vector<int> GCor_GSid,
                     vector<float> GCor_Gtd, vector<int> Laser,
                     vector<float> PEn, vector<int> Pnf, vector<int> Pnb,
                     vector<int> Psec, vector<int> Pquad, vector<float> Ptd) {
  laser_passed.resize(0);
  PEn_passed.resize(0);
  Pnf_passed.resize(0);
  Pnb_passed.resize(0);
  Pquad_passed.resize(0);
  Psec_passed.resize(0);
  Ptd_passed.resize(0);
  Ppid_passed.resize(0);

  for (int i = 0; i < PEn.size(); i++) {
    int pid = dc.Cut(PEn[i], Pnf[i], Psec[i]);

    if (pid >= 0) {
      laser_passed.push_back(Laser[i]);
      PEn_passed.push_back(PEn[i]);
      Pnf_passed.push_back(Pnf[i]);
      Pnb_passed.push_back(Pnb[i]);
      Pquad_passed.push_back(Pquad[i]);
      Psec_passed.push_back(Psec[i]);
      Ptd_passed.push_back(Ptd[i]);

      // JP convention: Target PID == 1, Beam PID == 0.
      if (pid == PID_BEAM)
        Ppid_passed.push_back(PID_BEAM);
      else if (PID_TARG)
        Ppid_passed.push_back(PID_TARG);
      else
        throw std::runtime_error("Invalid PID");
    }
  }

  int np_passed = PEn_passed.size();
  if (np_passed == 0)
    return;
  // passed vector is ordered in quadrants from 0 to 3

  if (np_passed == 1) { // 1-particle case, PID identified.
    // Here [0] signifies the detected particle. Can be beam or target.
    np = 1;
    tdpp = 0.;
    laser = laser_passed[0];
    ep[0] = PEn_passed[0];
    quad[0] = Pquad_passed[0];
    ring[0] = Pnf_passed[0];
    sect[0] = Pnb_passed[0];
    thp[0] = dc.GetPTh(Pnf_passed[0], Psec_passed[0]);
    php[0] = dc.GetPPhi(Pquad_passed[0], Pnb_passed[0], Psec_passed[0]);
    phr[0] = dc.GetQPhi(Pquad_passed[0], Pnb_passed[0], Psec_passed[0]);
    pid[0] = Ppid_passed[0];

    ng = 1 + GCor_GEn.size();
    td[0] = Ptd_passed[0];
    eg[0] = GEn;
    thg[0] = GTh;
    phg[0] = GPh;
    clu[0] = GCluid;
    cry[0] = GCid;
    seg[0] = GSid;
    for (int i = 0; i < ng - 1; i++) {
      eg[i + 1] = GCor_GEn[i];
      td[i + 1] = GCor_Gtd[i];
      thg[i + 1] = GCor_GTh[i];
      phg[i + 1] = GCor_GPh[i];
      clu[i + 1] = GCor_GCluID[i];
      cry[i + 1] = GCor_GCid[i];
      seg[i + 1] = GCor_GSid[i];
    }

    if (pid[0] == PID_TARG) { // Target detected, doppler corrections.
      // Target angles

      if (dc.UseKin()) { // Use the two-body kinematics. Not using this right?
        ep[0] = dc.GetTEnKinT(thp[0]);
        thr[0] = dc.GetBThLabT(thp[0]);
        er[0] = dc.GetBEnKinT(thp[0]);
      } else { // Use the particle energy and angle.
        ep[0] += dc.GetELoss(ep[0], dc.GetCDDeadLayer(), 1, "TS");
        er[0] = dc.GetBEn(PEn_passed[0], Pnf_passed[0], Psec_passed[0]);
        thr[0] = dc.GetBTh(Pnf_passed[0], Psec_passed[0]);
      }

      for (int i = 0; i < ng; i++) { // loop through gammas for angles and doppler correction.
        atg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
        abg[i] = dc.GammaAng(thr[0], phr[0], thg[i], phg[i]);
        etg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAt());
        ebg[i] = eg[i] * dc.DC(er[0], thr[0], phr[0], thg[i], phg[i], dc.GetAb());
      }
      // End of target detected.
    } else if (pid[0] == PID_BEAM) {
      
      if (dc.UseKin()) { // Use the two-body kinematics. Not using this right?
        ep[0] = dc.GetBEnKinB(thp[0]);
        thr[0] = dc.GetTThLabB(thp[0]);
        er[0] = dc.GetTEnKinB(thp[0]);
      } else { // Use the particle energy and angle.
        ep[0] += dc.GetELoss(ep[0], dc.GetCDDeadLayer(), 1, "BS");
        er[0] = dc.GetTEn(PEn_passed[0], Pnf_passed[0], Psec_passed[0]);
        thr[0] = dc.GetTTh(Pnf_passed[0], PEn_passed[0], Psec_passed[0]);
      }

      for (int i = 0; i < ng; i++) { // loop through gammas for angles and doppler correction 
        abg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
        atg[i] = dc.GammaAng(thr[0], phr[0], thg[i], phg[i]);
        ebg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAb());
        etg[i] = eg[i] * dc.DC(er[0], thr[0], phr[0], thg[i], phg[i], dc.GetAt());
      }
      // End of beam detected.
    } else
      throw std::runtime_error("Invalid PID during FillTree");

    tree->Fill();
    // 1p case done.
  } else if (np_passed == 2) {
    // Start of 2-particle case.
    // check quadrant correlation (diff = 2) and in.
    // Start checking if "good" 2p candidate.
    // Here [0] signifies the beam-like particle, [1] target-like particle.
    // time (ppwin, in ns) ; separate if 1n
    float time_diff = TMath::Abs(Ptd_passed[0] - Ptd_passed[1]); // 2p time difference in ns
    int quad_diff = TMath::Abs(Pquad_passed[0] - Pquad_passed[1]); // quadrant number difference

    // returns 0 for beam/target passed, 1 for target/beam passed, 
    // -1 for small 2p angles (ring > 10 (innermost = 16) for both)?
    int cut2 = dc.Cut_2p(PEn_passed[0], Pnf_passed[0], Psec_passed[0],
                         PEn_passed[1], Pnf_passed[1], Psec_passed[1]);  

    if ((quad_diff == 2) && (time_diff <= ppwin) && (cut2 >= 0)) { // we have good 2p candidate
      int ib, it;

      if (cut2 == PID_BEAM) {
        ib = 0;
        it = 1;
      } else if (cut2 == PID_TARG){
        ib = 1;
        it = 0;
      } else
        throw std::runtime_error("Invalid PID");

      laser = laser_passed[ib];
      np = 2;
      tdpp = Ptd_passed[it] - Ptd_passed[ib];
      int Bnf = Pnf_passed[ib];
      int Tnf = Pnf_passed[it];
      // ordering: 0 for beam, 1 for target (as for 110Sn)
      quad[1] = Pquad_passed[it];
      ring[1] = Pnf_passed[it];
      sect[1] = Pnb_passed[it];
      quad[0] = Pquad_passed[ib];
      ring[0] = Pnf_passed[ib];
      sect[0] = Pnb_passed[ib];

      ep[0] = PEn_passed[ib];
      ep[0] += dc.GetELoss(ep[ib], dc.GetCDDeadLayer(), 1, "BS");
      ep[1] = PEn_passed[it];
      ep[1] += dc.GetELoss(ep[it], dc.GetCDDeadLayer(), 1, "TS");
      er[0] = ep[1];
      er[1] = ep[0];

      thp[0] = dc.GetPTh(Bnf, Psec_passed[ib]);
      thp[1] = dc.GetPTh(Tnf, Psec_passed[it]);
      thr[0] = thp[1];
      thr[1] = thp[0];

      php[0] = dc.GetPPhi(Pquad_passed[ib], Pnb_passed[ib], Psec_passed[ib]);
      php[1] = dc.GetPPhi(Pquad_passed[it], Pnb_passed[it], Psec_passed[it]);
      phr[0] = php[1];
      phr[1] = php[0];

      // Correct PID scheme?
      pid[0] = Ppid_passed[ib];
      pid[1] = Ppid_passed[it];

      ng = 1 + GCor_GEn.size();
      td[0] = 0.5 * (Ptd_passed[0] + Ptd_passed[1]); // average
      eg[0] = GEn;
      thg[0] = GTh;
      phg[0] = GPh;
      clu[0] = GCluid;
      cry[0] = GCid;
      seg[0] = GSid;

      for (int i = 0; i < ng - 1; i++) {
        eg[i + 1] = GCor_GEn[i];
        td[i + 1] = GCor_Gtd[i];
        thg[i + 1] = GCor_GTh[i];
        phg[i + 1] = GCor_GPh[i];
        clu[i + 1] = GCor_GCluID[i];
        cry[i + 1] = GCor_GCid[i];
        seg[i + 1] = GCor_GSid[i];
      }

      // loop through gammas for angles and doppler correction
      for (int i = 0; i < ng; i++) { 
        abg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
        atg[i] = dc.GammaAng(thp[1], php[1], thg[i], phg[i]);
        ebg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAb());
        etg[i] = eg[i] * dc.DC(ep[1], thp[1], php[1], thg[i], phg[i], dc.GetAt());
      }

      tree->Fill();

    } else { // handle "broken" 2p event: either adjacent quads, time diff outside window, or identical pid
      // break into two 1p events
      for (int j = 0; j < np_passed; j++) {
        laser = laser_passed[j];
        np = 1;
        //	 b2p = 0;
        tdpp = 0.;
        ep[0] = PEn_passed[j];
        quad[0] = Pquad_passed[j];
        ring[0] = Pnf_passed[j];
        sect[0] = Pnb_passed[j];
        thp[0] = dc.GetPTh(Pnf_passed[j], Psec_passed[j]);
        php[0] = dc.GetPPhi(Pquad_passed[j], Pnb_passed[j], Psec_passed[j]);
        phr[0] = dc.GetQPhi(Pquad_passed[j], Pnb_passed[j], Psec_passed[j]);
        pid[0] = Ppid_passed[j];

        ng = 1 + GCor_GEn.size();
        td[0] = Ptd_passed[j]; // in particle
        eg[0] = GEn;
        thg[0] = GTh;
        phg[0] = GPh;
        clu[0] = GCluid;
        cry[0] = GCid;
        seg[0] = GSid;
        for (int i = 0; i < ng - 1; i++) {
          eg[i + 1] = GCor_GEn[i];
          td[i + 1] = GCor_Gtd[i];
          thg[i + 1] = GCor_GTh[i];
          phg[i + 1] = GCor_GPh[i];
          clu[i + 1] = GCor_GCluID[i];
          cry[i + 1] = GCor_GCid[i];
          seg[i + 1] = GCor_GSid[i];
        }

        if (pid[0] == PID_TARG) { // target detected, doppler corrections
          // Target angles

          // Use the two-body kinematics
          if (dc.UseKin()) {

            ep[0] = dc.GetTEnKinT(thp[0]);
            thr[0] = dc.GetBThLabT(thp[0]);
            er[0] = dc.GetBEnKinT(thp[0]);
          }

          // Or use the particle energy and angle
          else {
            ep[0] += dc.GetELoss(ep[0], dc.GetCDDeadLayer(), 1, "TS");
            er[0] = dc.GetBEn(PEn_passed[j], Pnf_passed[j], Psec_passed[j]);
            thr[0] = dc.GetBTh(Pnf_passed[j], Psec_passed[j]);
          }
          for (int i = 0; i < ng;
               i++) { // loop through gammas for angles and doppler correction
            atg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
            abg[i] = dc.GammaAng(thr[0], phr[0], thg[i], phg[i]);
            etg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAt());
            ebg[i] = eg[i] * dc.DC(er[0], thr[0], phr[0], thg[i], phg[i], dc.GetAb());
          }

        } else if (pid[0] == PID_BEAM) { // pid[0]==1, Beam detected

          //   // Use the two-body kinematics
          if (dc.UseKin()) {

            ep[0] = dc.GetBEnKinB(thp[0]);
            thr[0] = dc.GetTThLabB(thp[0]);
            er[0] = dc.GetTEnKinB(thp[0]);

          }
          //   // Or use the particle energy and angle
          else {
            ep[0] += dc.GetELoss(ep[0], dc.GetCDDeadLayer(), 1, "BS");
            er[0] = dc.GetTEn(PEn_passed[j], Pnf_passed[j], Psec_passed[j]);
            thr[0] = dc.GetTTh(Pnf_passed[j], PEn_passed[j], Psec_passed[j]);
          }
          for (int i = 0; i < ng;
               i++) { // loop through gammas for angles and doppler correction
            abg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
            atg[i] = dc.GammaAng(thr[0], phr[0], thg[i], phg[i]);
            ebg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAb());
            etg[i] = eg[i] * dc.DC(er[0], thr[0], phr[0], thg[i], phg[i], dc.GetAt());
          }
        } else {
          throw std::runtime_error("Invalid PID during sorting, broken 2p events.");
        }
        tree->Fill();
      }
    }

  // End of 2p events (both 'good' and 'broken').
  } else { // 3-4p events, loop through for 2p correlations and sort uncorrelated events as 1p.
    vector<pair<int, int>> v2p;
    vector<int> v2p_cut2;
    vector<int> v1p;
    v2p.clear();
    v2p_cut2.clear();
    v1p.clear();
    int qpattern = 0;
    vector<bool> unmatched;

    for (int j = 0; j < np_passed; j++) {
      unmatched.push_back(true);
    }

    for (int j = 0; j < np_passed; j++) {
      for (int k = j + 1; k < np_passed; k++) {
        if (j == k)
          continue;
        float time_diff = TMath::Abs(Ptd_passed[j] - Ptd_passed[k]); // 2p time difference in ns
        int quad_diff = TMath::Abs(Pquad_passed[j] - Pquad_passed[k]); // quadrant number difference

        // returns 0 for target-beam passed, 1 for beam/target passed, -1 for small 2p angles
        // (ring > 10 (innermost = 16) for both)
        int cut2 = dc.Cut_2p(PEn_passed[j], Pnf_passed[j], Psec_passed[j],
                             PEn_passed[k], Pnf_passed[k], Psec_passed[k]); 

        if (quad_diff == 2 && time_diff <= ppwin && cut2 >= 0) { // we have good 2p candidate
          v2p.push_back(make_pair(j, k));
          v2p_cut2.push_back(cut2);
          unmatched[j] = false;
          unmatched[k] = false;
        }
      }
    }

    for (int j = 0; j < np_passed; j++) {
      if (unmatched[j])
        v1p.push_back(j);
    } // checked all Np events

    // fill 1p events
    for (int j = 0; j < v1p.size(); j++) {
      laser = laser_passed[v1p[j]];
      np = 1;
      //	 b2p = 0;
      tdpp = 0.;
      ep[0] = PEn_passed[v1p[j]];
      quad[0] = Pquad_passed[v1p[j]];
      ring[0] = Pnf_passed[v1p[j]];
      sect[0] = Pnb_passed[v1p[j]];
      thp[0] = dc.GetPTh(Pnf_passed[v1p[j]], Psec_passed[v1p[j]]);
      php[0] = dc.GetPPhi(Pquad_passed[v1p[j]], Pnb_passed[v1p[j]], Psec_passed[v1p[j]]);
      phr[0] = dc.GetQPhi(Pquad_passed[v1p[j]], Pnb_passed[v1p[j]], Psec_passed[v1p[j]]);
      pid[0] = Ppid_passed[v1p[j]];

      ng = 1 + GCor_GEn.size();
      td[0] = Ptd_passed[v1p[j]]; // in particle
      eg[0] = GEn;
      thg[0] = GTh;
      phg[0] = GPh;
      clu[0] = GCluid;
      cry[0] = GCid;
      seg[0] = GSid;

      for (int i = 0; i < ng - 1; i++) {
        eg[i + 1] = GCor_GEn[i];
        td[i + 1] = GCor_Gtd[i];
        thg[i + 1] = GCor_GTh[i];
        phg[i + 1] = GCor_GPh[i];
        clu[i + 1] = GCor_GCluID[i];
        cry[i + 1] = GCor_GCid[i];
        seg[i + 1] = GCor_GSid[i];
      }

      if (pid[0] == PID_TARG) { // target detected, doppler corrections
        // Target angles

        if (dc.UseKin()) { // Use the two-body kinematics
          ep[0] = dc.GetTEnKinT(thp[0]);
          thr[0] = dc.GetBThLabT(thp[0]);
          er[0] = dc.GetBEnKinT(thp[0]);
        } else { // Or use the particle energy and angle
          ep[0] += dc.GetELoss(ep[0], dc.GetCDDeadLayer(), 1, "TS");
          er[0] = dc.GetBEn(PEn_passed[v1p[j]], Pnf_passed[v1p[j]], Psec_passed[v1p[j]]);
          thr[0] = dc.GetBTh(Pnf_passed[v1p[j]], Psec_passed[v1p[j]]);
        }

        for (int i = 0; i < ng; i++) { // loop through gammas for angles and doppler correction
          atg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
          abg[i] = dc.GammaAng(thr[0], phr[0], thg[i], phg[i]);
          etg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAt());
          ebg[i] = eg[i] * dc.DC(er[0], thr[0], phr[0], thg[i], phg[i], dc.GetAb());
        }

      } else if (pid[0] == PID_BEAM) { // pid[0]==1, Beam detected

        if (dc.UseKin()) { // Use the two-body kinematics
          ep[0] = dc.GetBEnKinB(thp[0]);
          thr[0] = dc.GetTThLabB(thp[0]);
          er[0] = dc.GetTEnKinB(thp[0]);
        } else {// Or use the particle energy and angle
          ep[0] += dc.GetELoss(ep[0], dc.GetCDDeadLayer(), 1, "BS");
          er[0] = dc.GetTEn(PEn_passed[v1p[j]], Pnf_passed[v1p[j]], Psec_passed[v1p[j]]);
          thr[0] = dc.GetTTh(Pnf_passed[v1p[j]], PEn_passed[v1p[j]], Psec_passed[v1p[j]]);
        }

        for (int i = 0; i < ng; i++) { // loop through gammas for angles and doppler correction
          abg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
          atg[i] = dc.GammaAng(thr[0], phr[0], thg[i], phg[i]);
          ebg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAb());
          etg[i] = eg[i] * dc.DC(er[0], thr[0], phr[0], thg[i], phg[i], dc.GetAt());
        }

      } else {
        throw std::runtime_error("Invalid PID during sorting, 1p events.");
      }
      tree->Fill();
    }

    // fill matched 2p events
    for (int j = 0; j < v2p.size(); j++) {
      int ib, it;
      if (v2p_cut2[j] == PID_BEAM) {
        ib = v2p[j].first;
        it = v2p[j].second;
      } else if (PID_TARG) {
        ib = v2p[j].second;
        it = v2p[j].first;
      } else
        throw std::runtime_error("Invalid PID");
      laser = laser_passed[ib]; // take ib as default
      np = 2;
      tdpp = Ptd_passed[it] - Ptd_passed[ib];
      int Bnf = Pnf_passed[ib];
      int Tnf = Pnf_passed[it];
      // ordering: 0 for beam, 1 for target (as for 110Sn)
      quad[1] = Pquad_passed[it];
      ring[1] = Pnf_passed[it];
      sect[1] = Pnb_passed[it];
      quad[0] = Pquad_passed[ib];
      ring[0] = Pnf_passed[ib];
      sect[0] = Pnb_passed[ib];

      ep[0] = PEn_passed[ib];
      ep[0] += dc.GetELoss(ep[ib], dc.GetCDDeadLayer(), 1, "BS");
      ep[1] = PEn_passed[it];
      ep[1] += dc.GetELoss(ep[it], dc.GetCDDeadLayer(), 1, "TS");
      er[0] = ep[1];
      er[1] = ep[0];
      thp[0] = dc.GetPTh(Bnf, Psec_passed[ib]);
      thp[1] = dc.GetPTh(Tnf, Psec_passed[it]);
      thr[0] = thp[1];
      thr[1] = thp[0];
      php[0] = dc.GetPPhi(Pquad_passed[ib], Pnb_passed[ib], Psec_passed[ib]);
      php[1] = dc.GetPPhi(Pquad_passed[it], Pnb_passed[it], Psec_passed[it]);
      phr[0] = php[1];
      phr[1] = php[0];
      pid[0] = Ppid_passed[v2p[j].first];
      pid[1] = Ppid_passed[v2p[j].second];

      ng = 1 + GCor_GEn.size();
      td[0] = 0.5 * (Ptd_passed[v2p[j].first] + Ptd_passed[v2p[j].second]); // average
      eg[0] = GEn;
      thg[0] = GTh;
      phg[0] = GPh;
      clu[0] = GCluid;
      cry[0] = GCid;
      seg[0] = GSid;

      for (int i = 0; i < ng - 1; i++) {
        eg[i + 1] = GCor_GEn[i];
        td[i + 1] = GCor_Gtd[i];
        thg[i + 1] = GCor_GTh[i];
        phg[i + 1] = GCor_GPh[i];
        clu[i + 1] = GCor_GCluID[i];
        cry[i + 1] = GCor_GCid[i];
        seg[i + 1] = GCor_GSid[i];
      }

      // Loop through gammas for angles and doppler correction.
      for (int i = 0; i < ng; i++) {
        abg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
        atg[i] = dc.GammaAng(thp[1], php[1], thg[i], phg[i]);
        ebg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAb());
        etg[i] = eg[i] * dc.DC(ep[1], thp[1], php[1], thg[i], phg[i], dc.GetAt());
      }

      tree->Fill();
    }

  // End of 3-4p events.
  }
  // End of FillTree
}

#endif

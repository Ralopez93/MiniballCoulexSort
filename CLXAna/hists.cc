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

#define PID_BEAM 0  // Beam, such as Sn.
#define PID_TARG 1  // Target, such as Pb

void hists::Initialise(doppler dc_) {

  /// Initialise all of the histograms that you want to fill
  /// If you add a new one, be careful not to break your memory limit.
  /// This causes unknown behaviour and may not give any errors.
  /// In the best case, you have a segmentation violation and you cannot
  /// work out where it comes from... Try deleting some histograms.
  /// Whatever you do, make sure you declare them in hists.hh file.
  dc = dc_;

  cout << "Initialising tree and histograms..." << endl;

  string hname, htitle;

  // particle branches
  tree = new TTree("doppler", "doppler");
  tree->Branch("run", &run_nbr, "run/I"); // The current run number or 0 if not provided by user.
  tree->Branch("laser", &laser, "laser/I");
  tree->Branch("np", &np, "np/I");          //
  tree->Branch("tdpp", &tdpp, "tdpp/D");    // time difference for 2p events; 0 for 1p event
  tree->Branch("time", time, "time[np]/D"); // Particle timestamp.
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
 * @brief Reset variables between each fill.
 * 
 */
void hists::resetVar() {
  laser_passed.resize(0);
  PEn_passed.resize(0);
  Pnf_passed.resize(0);
  Pnb_passed.resize(0);
  Pquad_passed.resize(0);
  Psec_passed.resize(0);
  Ptd_passed.resize(0);
  Ppid_passed.resize(0);
  PTheta_passed.resize(0);
  times_passed.resize(0);

  // Particle variables need to be reset,
  // since we have both 1P and 2P cases.
  for (int i = 0; i < 2; i++) {
    time[i] = -1;
    pid[i] = -1;
    quad[i] = -1;
    ring[i] = -1;
    sect[i] = -1;
    ep[i] = -1;
    er[i] = -1;
    thp[i] = -1;
    php[i] = -1;
    thr[i] = -1;
    phr[i] = -1;
  }
}

/**
 * @brief Check if pair of particles constitutes a "good" 2p event.
 * 
 * @param quad_diff 2p quadrant number difference.
 * @param time_diff 2p time difference in ticks of 25 ns.
 * @param ppwin Particle-particle time window?
 * @param cut2 Result from cut2().
 * 
 * @return True if criteria is fulfilled, else false.
 */
bool hists::isGood2p(int quad_diff, float time_diff, float ppwin, int cut2) {
  bool retval =  (
    (quad_diff == 2) &&
    (time_diff <= ppwin) &&
    (cut2 >= 0)
);

  return retval;
}

/**
 * @brief Go through 1p event.
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
 * @param GCor_Gtd Vector of time-difference to original gamma-ray {[25 ns]}
 */
void hists::doRoutine1P(float GEn, float GTh, float GPh, int GCluid, int GCid,
                        int GSid, vector<float> GCor_GEn, vector<float> GCor_GTh,
                        vector<float> GCor_GPh, vector<int> GCor_GCluID,
                        vector<int> GCor_GCid, vector<int> GCor_GSid,
                        vector<float> GCor_Gtd) {

  // Here [0] signifies the detected particle. Can be beam or target.
  laser = laser_passed[0];
  np = 1;
  tdpp = 0.;

  time[0] = times_passed[0];
  ep[0] = PEn_passed[0];
  er[0] = 0;

  quad[0] = Pquad_passed[0];
  ring[0] = Pnf_passed[0];
  sect[0] = Pnb_passed[0];

  thp[0] = PTheta_passed[0];
  thr[0] = 0;

  php[0] = dc.GetPPhi(Pquad_passed[0], Pnb_passed[0], Psec_passed[0]);
  phr[0] = dc.GetQPhi(Pquad_passed[0], Pnb_passed[0], Psec_passed[0]);
  pid[0] = Ppid_passed[0];

  // Gamma related variables.
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

  if (pid[0] == PID_BEAM) { // Beam detected, doppler corrections.
    if (dc.UseKin()) {
      // Use the two-body kinematics.
      thr[0] = dc.GetTThLabB(thp[0]);                             // Target angle.
      ep[0] = dc.GetBEnKinB(thp[0]);                              // Beam energy.
      er[0] = dc.GetTEnKinB(thr[0], thp[0]);                      // Target energy.
    } else {
      // Use the particle energy and angle.
      thr[0] = dc.GetTTh(PEn_passed[0], thp[0]);                  // Target angle.
      ep[0] += dc.GetELoss(ep[0], dc.GetCDDeadLayer(), 1, "BS");  // Beam energy.
      er[0] = dc.GetTEn(PEn_passed[0], thp[0]);                   // Target energy.
    }

    for (int i = 0; i < ng; i++) {
      // loop through gammas for angles and doppler correction 
      abg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
      atg[i] = dc.GammaAng(thr[0], phr[0], thg[i], phg[i]);
      ebg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAb());
      etg[i] = eg[i] * dc.DC(er[0], thr[0], phr[0], thg[i], phg[i], dc.GetAt());
    }
  } else if (pid[0] == PID_TARG) { // Target detected, doppler corrections.
    if (dc.UseKin()) {
      // Use the two-body kinematics.
      thr[0] = dc.GetBThLabT(thp[0]);                             // Beam angle
      ep[0] = dc.GetTEnKinT(thp[0]);                              // Target energy
      er[0] = dc.GetBEnKinT(thr[0], thp[0]);                      // Beam energy
    } else {
      // Use the particle energy and angle.
      thr[0] = dc.GetBTh(thp[0]);                                 // Beam angle.
      ep[0] += dc.GetELoss(ep[0], dc.GetCDDeadLayer(), 1, "TS");  // Target energy.
      er[0] = dc.GetBEn(PEn_passed[0], thp[0]);                   // Beam energy.
    }

    for (int i = 0; i < ng; i++) {
      // loop through gammas for angles and doppler correction.
      atg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
      abg[i] = dc.GammaAng(thr[0], phr[0], thg[i], phg[i]);
      etg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAt());
      ebg[i] = eg[i] * dc.DC(er[0], thr[0], phr[0], thg[i], phg[i], dc.GetAb());
    }
    // End of beam detected.
  } else
    throw std::runtime_error("Invalid PID during FillTree");

  // To remain consistent, populate second index of each particle list.
  ep[1] = er[0];
  er[1] = ep[0];
  
  thp[1] = thr[0];
  thr[1] = thp[0];

  php[1] = phr[0];
  phr[1] = php[0];

  tree->Fill();
}

/**
 * @brief Go through the 2p event. Some may be "bad" and need to be broken up into 1p events instead.
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
 * @param GCor_Gtd Vector of time-difference to original gamma-ray {[25 ns]}
 */
void hists::doRoutine2P(float GEn, float GTh, float GPh, int GCluid, int GCid,
                        int GSid, vector<float> GCor_GEn, vector<float> GCor_GTh,
                        vector<float> GCor_GPh, vector<int> GCor_GCluID,
                        vector<int> GCor_GCid, vector<int> GCor_GSid, 
                        vector<float> GCor_Gtd) {
  int np_passed = 2;
  // Start of 2-particle case. Check quadrant correlation (diff = 2) and in. Start checking if "good" 2p candidate.
  // Here [PID_BEAM] signifies the beam-like particle, [PID_TARG] target-like particle.
  float time_diff = TMath::Abs(Ptd_passed[PID_BEAM] - Ptd_passed[PID_TARG]); // 2p time difference in ticks of 25 ns
  int quad_diff = TMath::Abs(Pquad_passed[PID_BEAM] - Pquad_passed[PID_TARG]); // quadrant number difference

  // returns 0 for beam/target passed, 1 for target/beam passed, 
  // -1 for small 2p angles (ring > 10 (innermost = 16) for both)?
  int cut2 = dc.Cut_2p(PEn_passed[PID_BEAM], Pnf_passed[PID_BEAM], PTheta_passed[PID_BEAM],
                        PEn_passed[PID_TARG], Pnf_passed[PID_TARG], PTheta_passed[PID_TARG]);  

  if (isGood2p(quad_diff, time_diff, ppwin, cut2)) { 
    // we have good 2p candidate.
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

    time[PID_BEAM] = times_passed[ib];
    time[PID_TARG] = times_passed[it];

    // ordering: 0 for beam, 1 for target (as for 110Sn)
    quad[PID_TARG] = Pquad_passed[it];
    ring[PID_TARG] = Pnf_passed[it];
    sect[PID_TARG] = Pnb_passed[it];

    quad[PID_BEAM] = Pquad_passed[ib];
    ring[PID_BEAM] = Pnf_passed[ib];
    sect[PID_BEAM] = Pnb_passed[ib];

    ep[PID_BEAM] = PEn_passed[ib];
    ep[PID_TARG] = PEn_passed[it];

    thp[PID_BEAM] = PTheta_passed[ib];
    thp[PID_TARG] = PTheta_passed[it];

    if (dc.UseKin()) {
      // Use the two-body kinematics. Not using this right?
      ep[PID_BEAM] = dc.GetBEnKinB(thp[PID_BEAM]);
      ep[PID_TARG] = dc.GetTEnKinT(thp[PID_TARG]);
    } else {
      // Use the particle energy and angle.
      ep[PID_BEAM] += dc.GetELoss(ep[ib], dc.GetCDDeadLayer(), 1, "BS");
      ep[PID_TARG] += dc.GetELoss(ep[it], dc.GetCDDeadLayer(), 1, "TS");
    }

    // For recoil, just switch the indices.
    thr[PID_BEAM] = thp[PID_TARG];
    thr[PID_TARG] = thp[PID_BEAM];
    er[PID_BEAM] = ep[PID_TARG];
    er[PID_TARG] = ep[PID_BEAM];

    php[PID_BEAM] = dc.GetPPhi(Pquad_passed[ib], Pnb_passed[ib], Psec_passed[ib]);
    php[PID_TARG] = dc.GetPPhi(Pquad_passed[it], Pnb_passed[it], Psec_passed[it]);
    phr[PID_BEAM] = php[PID_TARG];
    phr[PID_TARG] = php[PID_BEAM];

    pid[PID_BEAM] = Ppid_passed[ib];
    pid[PID_TARG] = Ppid_passed[it];

    // Gamma related variables.
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
      abg[i] = dc.GammaAng(thp[PID_BEAM], php[PID_BEAM], thg[i], phg[i]);
      atg[i] = dc.GammaAng(thp[PID_TARG], php[PID_TARG], thg[i], phg[i]);
      ebg[i] = eg[i] * dc.DC(ep[PID_BEAM], thp[PID_BEAM], php[PID_BEAM], thg[i], phg[i], dc.GetAb());
      etg[i] = eg[i] * dc.DC(ep[PID_TARG], thp[PID_TARG], php[PID_TARG], thg[i], phg[i], dc.GetAt());
    }

    tree->Fill();

  } else {
    // handle "broken" 2p event: either adjacent quads, time diff outside window, or identical pid
    for (int j = 0; j < np_passed; j++) {
      // break into two 1p events
      laser = laser_passed[j];
      np = 1;

      tdpp = 0.;
      time[0] = times_passed[j];
      ep[0] = PEn_passed[j];
      quad[0] = Pquad_passed[j];
      ring[0] = Pnf_passed[j];
      sect[0] = Pnb_passed[j];
      thp[0] = PTheta_passed[j];
      php[0] = dc.GetPPhi(Pquad_passed[j], Pnb_passed[j], Psec_passed[j]);
      phr[0] = dc.GetQPhi(Pquad_passed[j], Pnb_passed[j], Psec_passed[j]);
      pid[0] = Ppid_passed[j];

      // Gamma related variables.
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

      if (pid[0] == PID_BEAM) {
        if (dc.UseKin()) {
          // Use the two-body kinematics
          thr[0] = dc.GetTThLabB(thp[0]);
          ep[0] = dc.GetBEnKinB(thp[0]);
          er[0] = dc.GetTEnKinB(thr[0], thp[0]);
        } else {
          // Or use the particle energy and angle
          ep[0] += dc.GetELoss(ep[0], dc.GetCDDeadLayer(), 1, "BS");
          er[0] = dc.GetTEn(PEn_passed[j], PTheta_passed[j]);
          thr[0] = dc.GetTTh(PEn_passed[j], PTheta_passed[j]);
        }

        for (int i = 0; i < ng; i++) { 
          // loop through gammas for angles and doppler correction
          abg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
          atg[i] = dc.GammaAng(thr[0], phr[0], thg[i], phg[i]);
          ebg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAb());
          etg[i] = eg[i] * dc.DC(er[0], thr[0], phr[0], thg[i], phg[i], dc.GetAt());
        }
      } else if (pid[0] == PID_TARG) {
        if (dc.UseKin()) {
          thr[0] = dc.GetBThLabT(thp[0]);

          // Use the two-body kinematics
          ep[0] = dc.GetTEnKinT(thp[0]);
          er[0] = dc.GetBEnKinT(thr[0], thp[0]);
        } else {
          // Or use the particle energy and angle
          ep[0] += dc.GetELoss(ep[0], dc.GetCDDeadLayer(), 1, "TS");
          er[0] = dc.GetBEn(PEn_passed[j], PTheta_passed[j]);
          thr[0] = dc.GetBTh(PTheta_passed[j]);
        }

        for (int i = 0; i < ng; i++) {
          // loop through gammas for angles and doppler correction
          atg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
          abg[i] = dc.GammaAng(thr[0], phr[0], thg[i], phg[i]);
          etg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAt());
          ebg[i] = eg[i] * dc.DC(er[0], thr[0], phr[0], thg[i], phg[i], dc.GetAb());
        }
      }

      tree->Fill();
    }
  }
}

/**
 * @brief Go through the 3p+ event. Most likely a pile-up issue.
 * 
 * Should be separated into 2p and 1p events.
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
 * @param GCor_Gtd Vector of time-difference to original gamma-ray {[25 ns]}
 * @param np_passed Number of particles that passed.
 */
void hists::doRoutineXP(float GEn, float GTh, float GPh, int GCluid, int GCid,
                        int GSid, vector<float> GCor_GEn, vector<float> GCor_GTh,
                        vector<float> GCor_GPh, vector<int> GCor_GCluID,
                        vector<int> GCor_GCid, vector<int> GCor_GSid,
                        vector<float> GCor_Gtd, int np_passed) {
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
      float time_diff = TMath::Abs(Ptd_passed[j] - Ptd_passed[k]); // 2p time difference in ticks of 25 ns
      int quad_diff = TMath::Abs(Pquad_passed[j] - Pquad_passed[k]); // quadrant number difference

      // returns 0 for target-beam passed, 1 for beam/target passed, -1 for small 2p angles
      // (ring > 10 (innermost = 16) for both)
      int cut2 = dc.Cut_2p(PEn_passed[j], Pnf_passed[j], Psec_passed[j],
                           PEn_passed[k], Pnf_passed[k], Psec_passed[k]); 

      if (isGood2p(quad_diff, time_diff, ppwin, cut2)) {
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
    time[0] = times_passed[v1p[j]];
    ep[0] = PEn_passed[v1p[j]];
    quad[0] = Pquad_passed[v1p[j]];
    ring[0] = Pnf_passed[v1p[j]];
    sect[0] = Pnb_passed[v1p[j]];
    thp[0] = PTheta_passed[v1p[j]];
    php[0] = dc.GetPPhi(Pquad_passed[v1p[j]], Pnb_passed[v1p[j]], Psec_passed[v1p[j]]);
    phr[0] = dc.GetQPhi(Pquad_passed[v1p[j]], Pnb_passed[v1p[j]], Psec_passed[v1p[j]]);
    pid[0] = Ppid_passed[v1p[j]];

    // Gamma related variables.
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

    if (pid[0] == PID_TARG) {
      if (dc.UseKin()) {
        // Use the two-body kinematics.
        thr[0] = dc.GetBThLabT(thp[0]);

        ep[0] = dc.GetTEnKinT(thp[0]);
        er[0] = dc.GetBEnKinT(thr[0], thp[0]);
      } else {
        // Or use the particle energy and angle.
        ep[0] += dc.GetELoss(ep[0], dc.GetCDDeadLayer(), 1, "TS");
        er[0] = dc.GetBEn(PEn_passed[v1p[j]], PTheta_passed[Pnf_passed[v1p[j]]]);
        thr[0] = dc.GetBTh(PTheta_passed[Pnf_passed[v1p[j]]]);
      }

      for (int i = 0; i < ng; i++) { // loop through gammas for angles and doppler correction
        atg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
        abg[i] = dc.GammaAng(thr[0], phr[0], thg[i], phg[i]);
        etg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAt());
        ebg[i] = eg[i] * dc.DC(er[0], thr[0], phr[0], thg[i], phg[i], dc.GetAb());
      }

    } else if (pid[0] == PID_BEAM) {
      if (dc.UseKin()) {
        // Use the two-body kinematics.
        thr[0] = dc.GetTThLabB(thp[0]);

        ep[0] = dc.GetBEnKinB(thp[0]);
        er[0] = dc.GetTEnKinB(thr[0], thp[0]);
      } else {
        // Or use the particle energy and angle.
        ep[0] += dc.GetELoss(ep[0], dc.GetCDDeadLayer(), 1, "BS");
        er[0] = dc.GetTEn(PEn_passed[v1p[j]], PTheta_passed[Pnf_passed[v1p[j]]]);
        thr[0] = dc.GetTTh(PEn_passed[v1p[j]], PTheta_passed[Pnf_passed[v1p[j]]]);
      }

      for (int i = 0; i < ng; i++) {
        // loop through gammas for angles and doppler correction
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
    } else if (v2p_cut2[j] == PID_TARG) {
      ib = v2p[j].second;
      it = v2p[j].first;
    } else
      throw std::runtime_error("Invalid PID");

    laser = laser_passed[ib]; // take ib as default
    np = 2;
    tdpp = Ptd_passed[it] - Ptd_passed[ib];
    time[PID_BEAM] = times_passed[ib];
    time[PID_TARG] = times_passed[it];

    // ordering: 0 for beam, 1 for target (as for 110Sn)
    quad[PID_TARG] = Pquad_passed[it];
    ring[PID_TARG] = Pnf_passed[it];
    sect[PID_TARG] = Pnb_passed[it];

    quad[PID_BEAM] = Pquad_passed[ib];
    ring[PID_BEAM] = Pnf_passed[ib];
    sect[PID_BEAM] = Pnb_passed[ib];

    ep[PID_BEAM] = PEn_passed[ib];
    ep[PID_TARG] = PEn_passed[it];

    if (dc.UseKin()) {
      // Use the two-body kinematics.
      ep[PID_BEAM] = dc.GetBEnKinB(thp[PID_BEAM]);
      ep[PID_TARG] = dc.GetTEnKinT(thp[PID_TARG]);
    } else {
      // Use the particle energy and angle.
      ep[PID_BEAM] += dc.GetELoss(ep[ib], dc.GetCDDeadLayer(), 1, "BS");
      ep[PID_TARG] += dc.GetELoss(ep[it], dc.GetCDDeadLayer(), 1, "TS");
    }

    thp[PID_BEAM] = PTheta_passed[ib];
    thp[PID_TARG] = PTheta_passed[it];

    php[PID_BEAM] = dc.GetPPhi(Pquad_passed[ib], Pnb_passed[ib], Psec_passed[ib]);
    php[PID_TARG] = dc.GetPPhi(Pquad_passed[it], Pnb_passed[it], Psec_passed[it]);

    pid[PID_BEAM] = Ppid_passed[v2p[j].first];
    pid[PID_TARG] = Ppid_passed[v2p[j].second];

    // For recoil, just switch the indices.
    er[PID_BEAM] = ep[PID_TARG];
    er[PID_TARG] = ep[PID_BEAM];

    thr[PID_BEAM] = thp[PID_TARG];
    thr[PID_TARG] = thp[PID_BEAM];

    phr[PID_BEAM] = php[PID_TARG];
    phr[PID_TARG] = php[PID_BEAM];

    // Gamma related variables.
    ng = 1 + GCor_GEn.size();
    td[PID_BEAM] = 0.5 * (Ptd_passed[v2p[j].first] + Ptd_passed[v2p[j].second]); // average
    eg[PID_BEAM] = GEn;
    thg[PID_BEAM] = GTh;
    phg[PID_BEAM] = GPh;
    clu[PID_BEAM] = GCluid;
    cry[PID_BEAM] = GCid;
    seg[PID_BEAM] = GSid;

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
      abg[i] = dc.GammaAng(thp[PID_BEAM], php[PID_BEAM], thg[i], phg[i]);
      atg[i] = dc.GammaAng(thp[PID_TARG], php[PID_TARG], thg[i], phg[i]);
      ebg[i] = eg[i] * dc.DC(ep[PID_BEAM], thp[PID_BEAM], php[PID_BEAM], thg[i], phg[i], dc.GetAb());
      etg[i] = eg[i] * dc.DC(ep[PID_TARG], thp[PID_TARG], php[PID_TARG], thg[i], phg[i], dc.GetAt());
    }

    tree->Fill();
  }
}

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
 * @param GCor_Gtd Vector of time-difference to original gamma-ray {[25 ns]}
 * @param Laser Laser flag (on = 1, off = 0).
 * @param PEn Particle energies [keV].
 * @param Pnf Annular (front) strip ID of particle (0 = outer; 15 inner). A.K.A front CD Ring ID.
 * @param Pnb Secular (back) strip ID of particle (0 to 12; clockwise wrt beam).
 * @param Psec Sector of C-REX (0 = FCD; 1 = FBarrel; 2 = BBarrel; 3 = BCD).
 * @param Pquad Detector (quadrant) number of particle
 * @param Ptd Particle-gamma time difference in 25 ns timestamps.
 * @param Ptimes Particle times in 25 ns timestamps.
 * @param cur_run_nbr The current run number or 0 if not provided by user.
 */
void hists::FillTree(float GEn, float GTh, float GPh, int GCluid, int GCid,
                     int GSid, vector<float> GCor_GEn, vector<float> GCor_GTh,
                     vector<float> GCor_GPh, vector<int> GCor_GCluID,
                     vector<int> GCor_GCid, vector<int> GCor_GSid,
                     vector<float> GCor_Gtd, vector<int> Laser,
                     vector<float> PEn, vector<int> Pnf, vector<int> Pnb,
                     vector<int> Psec, vector<int> Pquad, vector<float> Ptd,
                     vector<double> Ptimes, int cur_run_nbr) {

  resetVar();

  for (int i = 0; i < PEn.size(); i++) {
    // Saving theta early. Since it is randomly assigned (within some interval),
    // multiple calls to GetPTh() may prove problematic.

    float PTheta = dc.GetPTh(Pnf[i], Psec[i]);
    int pid = dc.Cut(PEn[i], Pnf[i], PTheta);

    if ((pid == PID_BEAM) || (pid == PID_TARG)) {
      laser_passed.push_back(Laser[i]);
      PEn_passed.push_back(PEn[i]);
      Pnf_passed.push_back(Pnf[i]);
      Pnb_passed.push_back(Pnb[i]);
      Pquad_passed.push_back(Pquad[i]);
      Psec_passed.push_back(Psec[i]);
      Ptd_passed.push_back(Ptd[i]);
      Ppid_passed.push_back(pid);
      PTheta_passed.push_back(PTheta);
      times_passed.push_back(Ptimes[i]);
    }
  }

  int np_passed = PEn_passed.size();
  if (np_passed == 0)
    return;
  // passed vector is ordered in quadrants from 0 to 3

  run_nbr = cur_run_nbr;

  // Introduce a flag that lets us only sort 1P events, skipping 2+?

  if (np_passed == 1) {
    // 1-particle case, PID identified.
    doRoutine1P(GEn, GTh, GPh, GCluid, GCid, GSid, GCor_GEn, GCor_GTh,
                GCor_GPh, GCor_GCluID, GCor_GCid, GCor_GSid, GCor_Gtd);
  } else if (np_passed == 2) {
    // 2-particle case.
    doRoutine2P(GEn, GTh, GPh, GCluid, GCid, GSid, GCor_GEn, GCor_GTh,
                GCor_GPh, GCor_GCluID, GCor_GCid, GCor_GSid, GCor_Gtd);
  } else {
    // 3-4p events, loop through for 2p correlations and sort uncorrelated events as 1p.
    doRoutineXP(GEn, GTh, GPh, GCluid, GCid, GSid, GCor_GEn, GCor_GTh,
                GCor_GPh, GCor_GCluID, GCor_GCid, GCor_GSid, GCor_Gtd,
                np_passed);
  }
} // End of FillTree

#endif

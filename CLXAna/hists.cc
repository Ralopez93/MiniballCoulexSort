#ifndef hists_cxx
#define hists_cxx

#ifndef hist_hh
#include "hists.hh"
#endif

#define PID_UNKNOWN -1
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

  tree->Branch("com", com, "com[np]/D"); // center-of-mass angle of particle

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

  ng = 0;
  laser = 0;
  run_nbr = 0;
  tdpp = 0;

  // Reset particle variables.
  for (int i = 0; i < 2; i++) {
    pid[i] = 0;
    quad[i] = 0;
    ring[i] = 0;
    sect[i] = 0;
    ep[i] = 0;
    er[i] = 0;
    thp[i] = 0;
    php[i] = 0;
    thr[i] = 0;
    phr[i] = 0;
    com[i] = 0;
  }

  // Reset gamma variables.
  for (int i = 0; i < 24; i++) {
    td[i] = 0;
    eg[i] = 0;
    ebg[i] = 0;
    etg[i] = 0;
    clu[i] = 0;
    cry[i] = 0;
    seg[i] = 0;
    thg[i] = 0;
    phg[i] = 0;
    tpg[i] = 0;
    abg[i] = 0;
    atg[i] = 0;
  }

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
bool hists::doRoutine1P(float GEn, float GTh, float GPh, int GCluid, int GCid,
                        int GSid, vector<float> GCor_GEn, vector<float> GCor_GTh,
                        vector<float> GCor_GPh, vector<int> GCor_GCluID,
                        vector<int> GCor_GCid, vector<int> GCor_GSid,
                        vector<float> GCor_Gtd) {

  // Here [0] signifies the detected particle. Can be beam or target.
  laser = laser_passed[0];
  np = 1;
  tdpp = 0.;

  time[0] = times_passed[0];
  ep[0] = 0.;
  er[0] = 0.;

  quad[0] = Pquad_passed[0];
  ring[0] = Pnf_passed[0];
  sect[0] = Pnb_passed[0];

  thp[0] = PTheta_passed[0];
  thr[0] = 0.;

  php[0] = dc.GetPPhi(Pquad_passed[0], Pnb_passed[0], Psec_passed[0]);
  phr[0] = dc.GetQPhi(Pquad_passed[0], Pnb_passed[0], Psec_passed[0]);
  pid[0] = Ppid_passed[0];

  com[0] = 0.;
  com[1] = 0.;

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

  if (pid[0] == PID_BEAM) { /***************************************/
    pid[1] = PID_TARG;
    if (dc.UseKin()) {
      // Use the two-body kinematics.
      ep[0] = dc.GetBEnKinB(thp[0]);                // Beam energy.
      er[0] = dc.GetTEnKinB(thp[0]);                // Target energy.
      thr[0] = dc.GetTThLabB(thp[0]);               // Target angle.
    } else {
      // Use the detected beam energy and angle.
      ep[0] = dc.GetBEnB(PEn_passed[0], thp[0]);    // Beam energy.
      er[0] = dc.GetTEnB(PEn_passed[0], thp[0]);    // Target energy.
      thr[0] = dc.GetTTh(thp[0]);                   // Target angle.
    }

    com[0] = dc.GetBThCoM(thp[0]);
    com[1] = dc.GetTThCoM(thr[0]);

    for (int i = 0; i < ng; i++) {
      // loop through gammas for angles and doppler correction 
      abg[i] = dc.GammaAng(thp[0], php[0], thg[i], phg[i]);
      atg[i] = dc.GammaAng(thr[0], phr[0], thg[i], phg[i]);
      ebg[i] = eg[i] * dc.DC(ep[0], thp[0], php[0], thg[i], phg[i], dc.GetAb());
      etg[i] = eg[i] * dc.DC(er[0], thr[0], phr[0], thg[i], phg[i], dc.GetAt());
    }
  } else if (pid[0] == PID_TARG) { /********************************/
    pid[1] = PID_BEAM;
    if (dc.UseKin()) {
      // Use the two-body kinematics.
      ep[0] = dc.GetTEnKinT(thp[0]);                // Target energy
      er[0] = dc.GetBEnKinT(thp[0]);                // Beam energy
      thr[0] = dc.GetBThLabT(thp[0]);               // Beam angle
    } else {
      // Use the detected target energy and angle.
      ep[0] = dc.GetTEnT(PEn_passed[0], thp[0]);    // Target energy.
      er[0] = dc.GetBEnT(PEn_passed[0], thp[0]);    // Beam energy.
      thr[0] = dc.GetBTh(thp[0]);                   // Beam angle.
    }

    com[0] = dc.GetTThCoM(thp[0]);
    com[1] = dc.GetBThCoM(thr[0]);

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

  return true;
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
bool hists::doRoutine2P(float GEn, float GTh, float GPh, int GCluid, int GCid,
                        int GSid, vector<float> GCor_GEn, vector<float> GCor_GTh,
                        vector<float> GCor_GPh, vector<int> GCor_GCluID,
                        vector<int> GCor_GCid, vector<int> GCor_GSid, 
                        vector<float> GCor_Gtd) {
  // Start of 2-particle case. Check quadrant correlation (diff = 2).
  // Here [PID_BEAM] signifies the beam-like particle, [PID_TARG] target-like particle.
  float time_diff = TMath::Abs(Ptd_passed[0] - Ptd_passed[1]); // 2p time difference in ticks of 25 ns
  int quad_diff = TMath::Abs(Pquad_passed[0] - Pquad_passed[1]); // quadrant number difference

  // Returns the identity of the first particle.
  // Or -1 if not a good 2p event.
  int cut2 = dc.Cut_2p(PEn_passed[0], Pnf_passed[0], PTheta_passed[0],
                       PEn_passed[1], Pnf_passed[1], PTheta_passed[1]);  

  // Small lambda function since the negative of the conditional statement got a bit convoluted.
  auto isGood2p = [] (int qd, float td, float ppwin, int cut2) {
    return ((qd == 2)     &&
            (td <= ppwin) &&
            ((cut2 == PID_BEAM) || (cut2 == PID_TARG)));
  };

  if (!isGood2p(quad_diff, time_diff, ppwin, cut2)) {
    std::cout << "Encountered broken 2p event, probably pileup. Rejecting event." << std::endl;
    return false;
  }

  // we have good 2p candidate.
  int ib = PID_UNKNOWN;
  int it = PID_UNKNOWN;

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

  ep[PID_BEAM] = 0.0;
  ep[PID_TARG] = 0.0;

  thp[PID_BEAM] = PTheta_passed[ib];
  thp[PID_TARG] = PTheta_passed[it];

  php[PID_BEAM] = dc.GetPPhi(Pquad_passed[ib], Pnb_passed[ib], Psec_passed[ib]);
  php[PID_TARG] = dc.GetPPhi(Pquad_passed[it], Pnb_passed[it], Psec_passed[it]);

  com[PID_BEAM] = dc.GetBThCoM(thp[PID_BEAM]);
  com[PID_TARG] = dc.GetTThCoM(thp[PID_TARG]);

  pid[PID_BEAM] = Ppid_passed[ib];
  pid[PID_TARG] = Ppid_passed[it];

  if (dc.UseKin()) {
    // Use the two-body kinematics. Not using this right?
    ep[PID_BEAM] = dc.GetBEnKinB(thp[PID_BEAM]);
    ep[PID_TARG] = dc.GetTEnKinT(thp[PID_TARG]);
  } else {
    // Use the particle energy and angle.
    ep[PID_BEAM] = dc.GetBEnB(PEn_passed[ib], thp[PID_BEAM]);
    ep[PID_TARG] = dc.GetTEnT(PEn_passed[it], thp[PID_TARG]);
  }

  // For recoil, just switch the indices.
  thr[0] = thp[PID_TARG];
  thr[1] = thp[PID_BEAM];

  er[0] = ep[PID_TARG];
  er[1] = ep[PID_BEAM];

  phr[0] = php[PID_TARG];
  phr[1] = php[PID_BEAM];

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

  return true;
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
 *              FCD = Forward CD? BCD = Backward CD?
 * @param Pquad Detector (quadrant) number of particle
 * @param Ptd Particle-gamma time difference in 25 ns timestamps.
 * @param Ptimes Particle times in 25 ns timestamps.
 * @param cur_run_nbr The current run number or 0 if not provided by user.
 */
bool hists::FillTree(float GEn, float GTh, float GPh, int GCluid, int GCid,
                     int GSid, vector<float> GCor_GEn, vector<float> GCor_GTh,
                     vector<float> GCor_GPh, vector<int> GCor_GCluID,
                     vector<int> GCor_GCid, vector<int> GCor_GSid,
                     vector<float> GCor_Gtd, vector<int> Laser,
                     vector<float> PEn, vector<int> Pnf, vector<int> Pnb,
                     vector<int> Psec, vector<int> Pquad, vector<float> Ptd,
                     vector<double> Ptimes, int cur_run_nbr) {

  resetVar();

  for (size_t i = 0; i < PEn.size(); i++) {
    // Saving theta early. Since it is randomly assigned (within some interval),
    // multiple calls to GetPTh() may prove problematic.

    float PTheta = dc.GetPTh(Pnf[i], Psec[i]);
    int cut = dc.Cut(PEn[i], Pnf[i], PTheta);

    if ((cut == PID_BEAM) || (cut == PID_TARG)) {
      laser_passed.push_back(Laser[i]);
      PEn_passed.push_back(PEn[i]);
      Pnf_passed.push_back(Pnf[i]);
      Pnb_passed.push_back(Pnb[i]);
      Pquad_passed.push_back(Pquad[i]);
      Psec_passed.push_back(Psec[i]);
      Ptd_passed.push_back(Ptd[i]);
      Ppid_passed.push_back(cut);
      PTheta_passed.push_back(PTheta);
      times_passed.push_back(Ptimes[i]);
    }
  }

  run_nbr = cur_run_nbr;
  if (GCor_GEn.size() > 23) {
    std::string err_msg = "Too many gammas (" + std::to_string(ng) + "): must adjust container size.";
    throw std::runtime_error(err_msg);
  }

  bool retval = false;
  int np_passed = PEn_passed.size();
  if (np_passed == 0) {
    std::cout << "Encountered 0p event. Rejected." << std::endl;
  } else if (np_passed == 1) {
    retval = doRoutine1P(GEn, GTh, GPh, GCluid, GCid, GSid, GCor_GEn, GCor_GTh,
                         GCor_GPh, GCor_GCluID, GCor_GCid, GCor_GSid, GCor_Gtd);
  } else if (np_passed == 2) {
    retval = doRoutine2P(GEn, GTh, GPh, GCluid, GCid, GSid, GCor_GEn, GCor_GTh,
                         GCor_GPh, GCor_GCluID, GCor_GCid, GCor_GSid, GCor_Gtd);
  } else {
    std::cout << "Encountered " << np_passed << "p event, probably pileup. Rejected." << std::endl;
  }

  return retval;
} // End of FillTree

#endif

#ifndef doppler_cxx
#define doppler_cxx

#ifndef doppler_hh
#include "doppler.hh"
#endif

#define PID_UNKNOWN -1
#define PID_BEAM 0
#define PID_TARG 1

TRandom3 doppler::rand = 1;

void doppler::ExpDefs(int Zb_, float Ab_, int Zt_, float At_, float Eb_,
                      float Ex_, float thick_, float depth_, float cddist_,
                      float cdoffset_, float deadlayer_, float contaminant_,
                      float spededist_, TCutG *Bcut_, TCutG *Tcut_,
                      string srimfile_, bool usekin_, bool usekinloss_, string calfile_) {

  /// Initialisation of experimental definitions from command line of config
  /// file
  Zb = Zb_;
  Ab = Ab_;
  Zt = Zt_;
  At = At_;
  Eb = Eb_;
  Ex = Ex_;
  thick = thick_;
  depth = depth_;
  cddist = cddist_;
  cdoffset = cdoffset_;
  deadlayer = deadlayer_;
  contaminant = contaminant_;
  spededist = spededist_;
  Bcut = Bcut_;
  Tcut = Tcut_;
  srimfile = srimfile_;
  usekin = usekin_;
  usekinloss = usekinloss_;
  calfile = calfile_;

  return;
}

void doppler::reactionEnergy() {

  /// Calculate the energy at interaction point
  Ereac = Eb * Ab;
  if (contaminant > 0)
    Ereac -= GetELoss(Ereac, contaminant, 0, "BC");
  Ereac -= GetELoss(Ereac, depth, 0, "BT");

  cout << "Reaction energy = " << Ereac / 1000. << " MeV";
  cout << " or " << Ereac / 1000. / Ab << " MeV/u" << endl;

  return;
}

void doppler::mbAngles(std::vector<Cluster> &clusters) {
  // call calibration
  Cal = new Calibration(calfile);
  for (int i = 0; i < 8; i++) { // loop over clusters
    Cluster cluster;

    mbg.SetupCluster(Cal->ClusterTheta(i), Cal->ClusterPhi(i),
                     Cal->ClusterAlpha(i), Cal->ClusterR(i), Cal->ZOffset());

    for (unsigned int j = 0; j < 3; j++) { // loop over cores
      Crystal crystal;

      gamma_theta[i][j][0] = mbg.GetCoreTheta(j) * TMath::DegToRad();
      gamma_phi[i][j][0] = mbg.GetCorePhi(j) * TMath::DegToRad();

      for (int k = 0; k < 6; k++) { // loop over segments

        gamma_theta[i][j][k + 1] = mbg.GetSegTheta(j, k) * TMath::DegToRad();
        gamma_phi[i][j][k + 1] = mbg.GetSegPhi(j, k) * TMath::DegToRad();
      }

      crystal.theta = gamma_theta[i][j][0] * TMath::RadToDeg();
      crystal.phi = gamma_phi[i][j][0] * TMath::RadToDeg();

      cluster.crystals.push_back(crystal);
    }
    clusters.push_back(cluster);
  }
}

bool doppler::stoppingpowers(bool BT, bool TT, bool BS, bool TS, bool BC, bool TC) {

  /// Initialisation of stopping powers
  bool success = true;
  gErrorIgnoreLevel = kWarning;
  for (int i = 0; i < 6; i++)
    gSP[i] = new TGraph();
  if (BT)
    success *= stoppingpowers("BT");
  if (TT)
    success *= stoppingpowers("TT");
  if (BS)
    success *= stoppingpowers("BS");
  if (TS)
    success *= stoppingpowers("TS");
  //  cout<<success<<" "<<endl;
  if (BC && contaminant > 0)
    success *= stoppingpowers("BC");
  if (TC && contaminant > 0)
    success *= stoppingpowers("TC");
  //  cout<<success<<" "<<contaminant<<endl;

  if (success)
    reactionEnergy();
  cout << "reactionEnergy done" << endl;

  return success;
}

bool doppler::stoppingpowers(string opt) {

  /// Open stopping power files and make TGraphs of data
  /// naming convention of files...

  string gElName[110] = {
      "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",
      "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti",
      "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
      "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",
      "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs",
      "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
      "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir",
      "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra",
      "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
      "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds"};

  unsigned int index = 0; // BT = 0, TT = 1, BS = 2, TS = 3
  string title = "Stopping powers for ";
  string srimfilename = srimfile;

  // Beam or target like..?
  if (opt.substr(0, 1) == "B") {

    srimfilename += "/" + convertInt(Ab + 0.5) + gElName[Zb - 1];
    title += convertInt(Ab + 0.5) + gElName[Zb - 1];

  }

  else if (opt.substr(0, 1) == "T") {

    srimfilename += "/" + convertInt(At + 0.5) + gElName[Zt - 1];
    title += convertInt(At + 0.5) + gElName[Zt - 1];
    index++;

  }

  else {

    cout << "opt must equal BT, TT, BS, TS, BC or TC \n";
    return false;
  }

  // Target, contaminant or silicon dead layer..?
  if (opt.substr(1, 1) == "T") {

    srimfilename += "_" + convertInt(At + 0.5) + gElName[Zt - 1] + ".txt";
    title += " in " + convertInt(At + 0.5) + gElName[Zt - 1];
    title += ";Ion energy [keV];Stopping power [MeV/(mg/cm^2)]";

  }

  else if (opt.substr(1, 1) == "S") {

    srimfilename += "_Si.txt";
    title += " in the Si dead layer";
    title += ";Ion energy [keV];Stopping power [MeV/(mg/cm^2)]";
    index += 2;

  }

  else if (opt.substr(1, 1) == "C") {

    srimfilename += "_contaminant.txt";
    title += " in the contaminant layer";
    title += ";Ion energy [keV];Stopping power [MeV/(mg/cm^2)]";
    index += 4;

  }

  else {

    cout << "opt must equal BT, TT, BS or TS \n";
    return false;
  }

  ifstream infile;
  infile.open(srimfilename.c_str(), ios::in);

  if (!infile.is_open()) {

    cout << "Cannot open " << srimfilename << endl;
    return false;
  }
  cout << "Using SRIM file: " << srimfilename << ". Title: " << title.c_str() << endl;
  gSP[index]->SetTitle(title.c_str());

  string line, units, tmp_str;
  stringstream line_ss;
  bool endflag = false;
  double BEn, nucl, elec, total, tmp_dbl;
  int p = 0;

  // Test file format
  getline(infile, line);
  if (line.substr(0, 5) == " ====") {

    while (line.substr(0, 5) != "  ---") // Err, what if you use SRIM-2006,
                                         // which doesn't have the spaces?
      getline(infile, line);

    getline(infile, line); // read first line of data
    //		cout<<line<<endl;
  }

  while (!infile.eof() && !endflag) {

    // Read in data
    line_ss.str("");
    line_ss << line;
    line_ss >> BEn >> units >> nucl >> elec >> tmp_dbl >> tmp_str >> tmp_dbl >>
        tmp_str;

    if (units == "eV")
      BEn *= 1E-3;
    else if (units == "keV")
      BEn *= 1E0;
    else if (units == "MeV")
      BEn *= 1E3;
    else if (units == "GeV")
      BEn *= 1E6;

    // MeV / ( mg / cm^2 )
    total = nucl + elec; 
    gSP[index]->SetPoint(p, BEn, total);

    // Get next line
    getline(infile, line);
    p++;

    // If we've reached the end, stop
    if (line.substr(0, 9) == "---------")
      endflag = true;
    if (line.substr(0, 9) == " Multiply")
      endflag = true;
  }
  infile.close();
  cout << "stopping powers defined" << endl;
  TCanvas *c = new TCanvas();
  gSP[index]->Draw("A*");
  string pdfname = srimfilename.substr(0, srimfilename.find_last_of(".")) + ".pdf";
  c->SetLogx();
  c->SaveAs(pdfname.c_str());

  delete c;

  return true;
}

/**
 * @brief Check if entry passes any particle gates. Graphical cuts are used
 * if they are given in the config file or on the command line with the -cut option.
 * 
 * @param PEn Particle energy.
 * @param ring Annular (front) strip ID for particle.
 * @param PTheta Particle scattering angle.
 *
 * @return PID of particle, -1 if outside gates.
 */
int doppler::Cut(float PEn, float ring, float PTheta) {

  int identity = PID_UNKNOWN;
  float ang = PTheta * TMath::RadToDeg();

  if (Bcut->GetN() > 0 && Tcut->GetN() > 0) {
    // Graphical cuts given at command line.
    if (Bcut->IsInside(ang, PEn / 1000.))
      identity = PID_BEAM;
    else if (Tcut->IsInside(ang, PEn / 1000.))
      identity = PID_TARG;

  } else if (Ab > At) {
    // inverse kinematics, include overlap region
    double a = 349.07, b = -4.997, c = -0.0145;

    if (PEn / 1000. <= (a + b * ang + c * ang * ang) && ring < 15)
      identity = PID_BEAM;
    else if (ring < 15)
      identity = PID_TARG;
  } else {
    // normal kinematics with Beam/Target separation
    double a = 497.602, b = -4.67677, c = -0.0274333, l = 0;
    double d = 435.186, e = -7.84811, f = 0.0199164, k = 0;
    double g = 0, h = 0, i = 0, n = 0;

    if (PEn / 1000. <= (a + b * ang + c * ang * ang + l * ang * ang * ang) &&
        PEn / 1000. >= (d + e * ang + f * ang * ang + k * ang * ang * ang))

      identity = PID_BEAM;

    else if (PEn / 1000. <= (d + e * ang + f * ang * ang + k * ang * ang * ang) &&
             PEn / 1000. >= (g + h * ang + i * ang * ang + n * ang * ang * ang))

      identity = PID_TARG;
  }

  return identity;
}

/**
 * @brief Check if entry passes the 2 particle condition. It calls Cut() twice
 * with each of the two particles passed to this function. If one of them is a
 * particle and one of them is a target, you get a good return.
 * 
 * @param PEn1 Energy for particle 1.
 * @param ring1 Annular (front) strip ID for particle 1.
 * @param PTheta1 Scattering angle for particle 1.
 * @param PEn2 Energy for particle 2.
 * @param ring2 Annular (front) strip ID for particle 2.
 * @param PTheta2 Scattering angle for particle 2.
 *
 * @return PID of first particle or -1 if conditions are not fulfilled.
 */
int doppler::Cut_2p(float PEn1, float ring1, float PTheta1,
                    float PEn2, float ring2, float PTheta2) {

  int identity = PID_UNKNOWN;

  if ((Cut(PEn1, ring1, PTheta1) == PID_BEAM) &&
      (Cut(PEn2, ring2, PTheta2) == PID_TARG)) {
    identity = PID_BEAM;
  } else if ((Cut(PEn1, ring1, PTheta1) == PID_TARG) &&
             (Cut(PEn2, ring2, PTheta2) == PID_BEAM)) {
    identity = PID_TARG;
  }

    // JP: test and check later
    // If the angle is small, it's unlikely to be a real 2h event
    // if( nf1 > 10 || nf2 > 10 ){
    //   cout<<"in doppler.cc: rings "<<nf1<<" vs "<<nf2<<endl;
    //   identity = -1;
    // }

  return identity;
}

float doppler::GetCDOffset() {

  /// Return offset of the CD in the phi rotation from vertical in degrees

  return cdoffset;
}

float doppler::GetSpedeDist() {

  /// Return distance of Spede detector in mm

  return spededist;
}

float doppler::GetCDDeadLayer() {

  /// Return dead layer of the Si in mm

  return deadlayer;
}

int doppler::GetZb() {

  /// Return Z of the projectile as an int

  return Zb;
}

float doppler::GetAb() {

  /// Return A of the projectile as a float

  return Ab;
}

int doppler::GetZt() {

  /// Return Z of the target as an int

  return Zt;
}

float doppler::GetAt() {

  /// Return A of the target as a float

  return At;
}

/**
 * @brief Get particle scattering angle based on ring number and sector.
 * Should only be called once per entry, since angle is normalized.
 * 
 * @param nf 
 * @param sector 
 *
 * @return The particle scattering angle in radians. 
 */
float doppler::GetPTh(float ring, int sector) {

  /// Returns theta angle from ann strip number in radians

  float angle = 0.0;

  // Forward CD - Standard CD
  if (sector == 4) {
    // r0 = 9.0 mm
    // nf is ring number, goes from 0 to 15.
    // So for example, r1 = r0 + 0.5 * 2 (basically adding half the pitch).
    double angle_lower = TMath::ATan((9.0 + (15.5 - ring) * 2.0 - 1.0) / cddist);
    double angle_upper = TMath::ATan((9.0 + (15.5 - ring) * 2.0 + 1.0) / cddist);

    // Avoid using global seed. Just want to keep the Theta calculations consistent.
    // Set to static to avoid resetting the RNG.
    // Variable lifetime is now same as program and not function.
    static TRandom3 rng(1);
    angle = rng.Uniform(angle_lower, angle_upper);
  }

  // Forward CD - CREX
  if (sector == 0)
    angle = TMath::ATan((9.0 + (ring + 0.5) * 2.0) / cddist);

  // Forwards Barrel
  if (sector == 1)
    angle = 0.5 * TMath::Pi() - TMath::ATan((8.0 + (ring + 0.5) * 3.125) / 29.0);

  // Backwards Barrel
  if (sector == 2)
    angle = 0.5 * TMath::Pi() + TMath::ATan((8.0 + (ring + 0.5) * 3.125) / 29.0);

  // Backwards CD
  if (sector == 3)
    angle = TMath::Pi() - TMath::ATan((9.0 + (ring + 0.5) * 2.0) / 64.0);

  if (sector != 4)
    std::cout << "Sector: " << sector << ", ring: " << ring << std::endl;

  return angle;
}

/**
 * @brief Get phi angle from quadrant and ohm strip number in radians.
 * 
 * @param quad CD quadrant.
 * @param seg CD Segment
 * @param sector CD sector (forward, backward, ...)?
 * 
 * @return Phi angle. 
 */
float doppler::GetPPhi(int quad, int seg, int sector) {

  float ph_det[4];
  if (sector == 4) { // standard CD
    ph_det[0] = 0.0 + cdoffset;   // top
    ph_det[1] = 90.0 + cdoffset;  // right
    ph_det[2] = 180.0 + cdoffset; // bottom
    ph_det[3] = 270.0 + cdoffset; // left
  }

  else { // CREX and TREX
    ph_det[0] = 0.0 + cdoffset;   // top
    ph_det[1] = 180.0 + cdoffset; // bottom
    ph_det[2] = 270.0 + cdoffset; // left
    ph_det[3] = 90.0 + cdoffset;  // right
  }

  float pphi = ph_det[quad];
  if (sector == 4)
    pphi += seg * 7.0;   // standard CD
  else if (sector < 4) { // CREX and TREX

    pphi += 1.75; // centre of first strip
    if (seg < 4)
      pphi += seg * 3.5; // first 4 strips singles (=4 segs)
    else if (seg < 12)
      pphi += 14. + (seg - 4) * 7.0; // middle 16 strips doubles (=8 segs)
    else
      pphi += 70. + (seg - 12) * 3.5; // last 4 strips singles (=4 segs)
  }

  if (pphi < 360.)
    return pphi * TMath::DegToRad();
  else
    return (pphi - 360.) * TMath::DegToRad();
}

/**
 * @brief Get phi angle of B/T using angle of T/B.
 * 
 * @param quad CD quadrant.
 * @param seg CD Segment
 * @param sector CD sector (forward, backward, ...)?
 * 
 * @return Phi angle. 
 */
float doppler::GetQPhi(int quad, int seg, int sector) {
  return GetPPhi(quad, seg, sector) + TMath::Pi();
}

/**
 * @brief Get target scattering angle from beam energy and 
 * beam scattering angle
 * 
 * @param BEn 
 * @param BTheta 
 * 
 * @return Target scattering angle.
 */
float doppler::GetTTh(float BEn, float BTheta) {
  float tau = Ab / At;
  float Eprime = Eb * Ab - Ex * (1 + tau);
  float epsilon = TMath::Sqrt(Eb * Ab / Eprime);
  float x, y, TTh;
  if (tau > 1) { // inverse kinematics: maximum scattering angle may be exceeded...
    y = TMath::ASin(1. / (tau * epsilon)); // maximum projectile angle in lab
    if (BTheta < y)
      y = BTheta;
    y = TMath::Tan(y);
  } else {
    y = TMath::Tan(BTheta); // y = tan(Theta_projlab)
  }
  if (tau > 1 && rand.Gaus(BEn, 30000.) < 50000.) {
    x = (y * y * epsilon * tau + TMath::Sqrt(-y * y * epsilon * epsilon * tau * tau + y * y + 1)) / (1 + y * y);
  } else {
    x = (y * y * epsilon * tau -
         TMath::Sqrt(-y * y * epsilon * epsilon * tau * tau + y * y + 1)) /
        (1 + y * y);
  }
  TTh = TMath::ATan(
      TMath::Sqrt(1 - x * x) /
      (epsilon + x)); // choose kinematic branch using energy cut... as I
                      // haven't a clue how to do it any other way?!
  if (TTh < 0)
    TTh += TMath::Pi();
  return TTh;
}

/**
 * @brief Get beam scattering angle from target angle.
 * 
 * @param TTheta Target scattering angle.
 * 
 * @return Beam scattering angle.
 */
float doppler::GetBTh(float TTheta) {

  float tau = Ab / At;
  float Eprime = Eb * Ab - Ex * (1 + tau);
  float epsilon = TMath::Sqrt(Eb * Ab / Eprime);
  float x, y, BTh;
  y = TMath::Tan(TTheta);
  x = (y * y * epsilon - TMath::Sqrt(-y * y * epsilon * epsilon + y * y + 1)) / (1 + y * y);
  BTh = TMath::ATan(TMath::Sqrt(1 - x * x) / (tau * epsilon + x));
  if (BTh < 0)
    BTh += TMath::Pi();
  return BTh;
}

/**
 * @brief Get target energy from beam energy and beam scattering angle.
 * 
 * @param BEn Beam energy.
 * @param BTheta Beam scattering angle. 
 *
 * @return Target energy.
 */
float doppler::GetTEn(float BEn, float BTheta) {
  // Correct for dead layer loss
  float dist = TMath::Abs(deadlayer / TMath::Cos(BTheta));
  float Eproj = BEn + GetELoss(BEn, dist, 1, "BS");

  // Trace energy loss back through target to get energy at interaction point
  dist = TMath::Abs((thick - depth) / TMath::Cos(BTheta));
  Eproj += GetELoss(Eproj, dist, 1, "BT");

  float Etarg = Ereac - Eproj;
  if (Etarg < 0.1)
    return 0.1; // recoil is stopped in target

  float angle = GetTTh(BEn, BTheta);
  if (angle < 0.501 * TMath::Pi() && angle > 0.499 * TMath::Pi())
    return 0.1; // stopped

  dist = TMath::Abs((thick - depth) / TMath::Cos(angle));
  Etarg -= GetELoss(Etarg, dist, 0, "TT");

  if (Etarg < 0.1)
    return 0.1;
  else
    return Etarg;
}

/**
 * @brief Get beam energy from target energy and target scattering angle.
 * 
 * @param TEn Target energy.
 * @param TTheta Target scattering angle.
 * 
 * @return Beam energy.
 */
float doppler::GetBEn(float TEn, float TTheta) {
  // Correct for dead layer loss
  float dist = TMath::Abs(deadlayer / TMath::Cos(TTheta));
  float Etarg = TEn + GetELoss(TEn, dist, 1, "TS");

  // Trace energy loss back through target to get energy at interaction point
  dist = TMath::Abs((thick - depth) / TMath::Cos(TTheta));
  Etarg += GetELoss(Etarg, dist, 1, "TT");

  float Eproj = Ereac - Etarg;
  if (Eproj < 0.1)
    return 0.1; // projectile is stopped in target

  float angle = GetBTh(TTheta);
  if (angle < 0.501 * TMath::Pi() && angle > 0.499 * TMath::Pi())
    return 0.1; // stopped

  dist = TMath::Abs((thick - depth) / TMath::Cos(angle));
  Eproj -= GetELoss(Eproj, dist, 0, "BT");

  if (Eproj < 0.1)
    return 0.1;
  else
    return Eproj;
}
  /** 
   */
/**
 * @brief Returns the energy loss at a given initial energy and distance travelled in the target,
 *  the contaminant layer or Si dead layer Ei is the initial energy in keV,
 *  return value is also in keV dist is the distance travelled in the target in mg/cm2.
 * 
 *  opt = 0 calculates normal energy loss as particle moves through target (default)
 *  opt = 1 calculates energy increase, i.e. tracing particle back
 *  to reaction point combo = "BT", "TT", "BC", "TC", "BS" or "TS" for the
 *  beam in target, target in target, beam in contaminant, target in
 *  contaminant, beam in Si or target in Si, respectively.
 * 
 *  Stopping power data  is taken from SRIM the output files must be placed in the
 *  './srim/' folder with the format 62Fe_109Ag.txt, 62Fe_Si.txt, 109Ag_109Ag.txt or
 *  109Ag_Si.txt, for combo = "BT", "TT", "BS" and "TS", repsectively. The
 *  srim file should be in units of MeV/(mg/cm^2)
 * 
 * @param Ei Initial energy.
 * @param dist distance traveled.
 * @param opt Option.
 * @param combo Reaction combo.
 * 
 * @return Energy loss.
 */
float doppler::GetELoss(float Ei, float dist, int opt, string combo) {

  double dedx = 0;
  int Nmeshpoints = 20; // number of steps to take in integration
  double dx = dist / (double)Nmeshpoints;
  double E = Ei;

  for (int i = 0; i < Nmeshpoints; i++) {

    if (E < 1000.)
      break; // when we fall below 1 MeV we assume maximum energy loss

    if (combo == "BT")
      dedx = gSP[0]->Eval(E);
    else if (combo == "TT")
      dedx = gSP[1]->Eval(E);
    else if (combo == "BS")
      dedx = gSP[2]->Eval(E);
    else if (combo == "TS")
      dedx = gSP[3]->Eval(E);
    else if (combo == "BC")
      dedx = gSP[4]->Eval(E);
    else if (combo == "TC")
      dedx = gSP[5]->Eval(E);

    if (opt == 1)
      E += 1000. * dedx * dx;
    else
      E -= 1000. * dedx * dx;
  }

  if (opt == 0)
    return Ei - E;
  else
    return E - Ei;
}

/**
 * @brief Calculate the energy loss for the beam in target.
 * Should only be used together with Catkin 2B procedure (usekin)..
 * 
 * @param BEn Beam energy.
 * @param BTh Beam scattering angle
 * 
 * @return Calcualted energy loss for beam in target.
 */
float doppler::GetBKinLoss(float BEn, float BTh) {
  float dist = TMath::Abs((thick - depth) / TMath::Cos(BTh));

 return GetELoss(BEn, dist, 0, "BT");
}

/**
 * @brief Calculate the energy loss for the target in target.
 * Should only be used together with Catkin 2B procedure (usekin).
 * 
 * @param TEn Target energy.
 * @param TTh Target scattering angle.
 * 
 * @return Calcualted energy loss for beam in target.
 */
float doppler::GetTKinLoss(float TEn, float TTh) {
  float dist = TMath::Abs((thick - depth) / TMath::Cos(TTh));

 return GetELoss(TEn, dist, 0, "TT");
}

/**
 * @brief Calculate the beam angle in the lab from the centre of mass angle in radians.
 * 
 * @param CoM theta angle of the beam/target in the centre of mass frame.
 * 
 * @return The beam angle. 
 */
float doppler::GetBThLab(float CoM) {
  float tau = Ab / At;
  float Eprime = Ereac - Ex * (1 + tau);
  float epsilon = TMath::Sqrt(Ereac / Eprime);

  float y = TMath::Sin(CoM) / (TMath::Cos(CoM) + tau * epsilon);

  float BTh = TMath::ATan(y);
  if (BTh < 0.)
    BTh += TMath::Pi();

  return BTh;
}

/**
 * @brief Calculate the target angle in the lab from the centre of mass angle in radians.
 * 
 * @param CoM theta angle of the beam/target in the centre of mass frame.
 * 
 * @return The target angle. 
 */
float doppler::GetTThLab(float CoM) {
  float tau = Ab / At;
  float Eprime = Ereac - Ex * (1 + tau);
  float epsilon = TMath::Sqrt(Ereac / Eprime);

  // y = tan(theta_lab)
  float y =
      TMath::Sin(TMath::Pi() - CoM) / (TMath::Cos(TMath::Pi() - CoM) + epsilon);

  float TTh = TMath::ATan(y);
  if (TTh < 0.)
    TTh += TMath::Pi();

  return TTh;
}

/**
 * @brief Calculate the beam angle in the lab from the target lab angle
 * 
 * @param TTh theta angle of the target in the laboratory frame.
 * @param kinflag kinematics flag such that true is the backwards solution
 *                (i.e. CoM > 90 deg)
 * @return Beam angle from target angle.
 */
float doppler::GetBThLabT(float TTh, bool kinflag) {
  return GetBThLab(GetTThCoM(TTh, kinflag));
}

/**
 * @brief Calculate the target angle in the lab from the beam lab angle
 * 
 * @param TTh theta angle of the target in the laboratory frame.
 * @param kinflag kinematics flag such that true is the backwards solution
 *                (i.e. CoM > 90 deg)
 * @return Beam angle from target angle.
 */
float doppler::GetTThLabB(float BTh, bool kinflag) {
  return GetTThLab(GetBThCoM(BTh, kinflag));
}

/**
 * @brief Calculates CoM scattering angle from the
 * beam laboratory angle in radians.
 * 
 * @param BTh theta angle of the beam in laboratory frame
 * @param kinflag kinematics flag such that true is the backwards solution
 *                (i.e. CoM > 90 deg)
 *
 * @return The Center-of-Mass scattering angle. 
 */
float doppler::GetBThCoM(float BTh, bool kinflag) {
  float tau = Ab / At;
  float Eprime = Ereac - Ex * (1 + tau);
  float epsilon = TMath::Sqrt(Ereac / Eprime);

  // maximum scattering angle may be exceeded...
  float maxang = TMath::ASin(1. / (tau * epsilon));
  if (tau * epsilon > 1 && BTh > maxang)
    BTh = maxang;

  if (kinflag && tau * epsilon < 1) {

    cerr << "Only one solution for the beam, kinflag = false" << endl;
    kinflag = false;
  }

  float y = epsilon * tau * TMath::Sin(BTh);
  if (kinflag)
    y = TMath::ASin(-y);
  else
    y = TMath::ASin(y);

  float CoM = BTh + y;

  if (CoM < 0.)
    CoM += TMath::Pi();
  if (CoM > TMath::Pi())
    CoM -= TMath::Pi();

  return CoM;
}

/**
 * @brief Calculates CoM scattering angle from the
 * target laboratory angle in radians.
 * 
 * @param TTh theta angle of the target in laboratory frame
 * @param kinflag kinematics flag such that true is the backwards solution
 *                (i.e. CoM > 90 deg)
 *
 * @return The Center-of-Mass scattering angle. 
 */
float doppler::GetTThCoM(float TTh, bool kinflag) {
  float tau = Ab / At;
  float Eprime = Ereac - Ex * (1 + tau);
  float epsilon = TMath::Sqrt(Ereac / Eprime);

  // maximum scattering angle may be exceeded...
  float maxang = TMath::ASin(1. / (epsilon));
  if (TTh > maxang)
    TTh = maxang;

  float y = epsilon * TMath::Sin(TTh);
  if (kinflag)
    y = TMath::ASin(-y);
  else
    y = TMath::ASin(y);

  float CoM = TTh + y;
  CoM = TMath::Pi() - CoM;
  if (CoM < 0.)
    CoM += TMath::Pi();
  if (CoM > TMath::Pi())
    CoM -= TMath::Pi();

  return CoM;
}

/**
 * @brief Calculate the beam energy for a given centre of mass angle
 * using two-body kinematics calculations only.
 * 
 * @param CoM theta angle of the beam/particle in the centre of mass frame.
 * 
 * @return The beam energy.
 */
float doppler::GetBEnKin(float CoM) {
  float tau = Ab / At;
  float Eprime = Ereac - Ex * (1 + tau);
  float epsilon = TMath::Sqrt(Ereac / Eprime);

  float Eproj = TMath::Power(At / (At + Ab), 2.0);
  Eproj *= 1. + tau * tau * epsilon * epsilon + 2. * tau * epsilon * TMath::Cos(CoM);
  Eproj *= Eprime;

  return Eproj;
}

/**
 * @brief Calculate the target energy for a given centre of mass angle
 * using two-body kinematics calculations only.
 * 
 * @param CoM theta angle of the beam/particle in the centre of mass frame.
 * 
 * @return The target energy.
 */
float doppler::GetTEnKin(float CoM) {
  float tau = Ab / At;
  float Eprime = Ereac - Ex * (1 + tau);
  float epsilon = TMath::Sqrt(Ereac / Eprime);

  float Etarg = (At * Ab) / TMath::Power((At + Ab), 2.0);
  Etarg *= 1. + epsilon * epsilon + 2. * epsilon * TMath::Cos(TMath::Pi() - CoM);
  Etarg *= Eprime;

  return Etarg;
}

/**
 * @brief Calculate beam energy for a given beam scattering angle
 * using two-body kinematics calculations only.
 * 
 * This is used when beam was detected.
 * 
 * @param BTh theta angle of the beam in laboratory frame (detected).
 * @param kinflag kinematics flag such that true is the backwards solution (i.e. CoM > 90 deg).
 * 
 * @return The calculated beam energy. 
 */
float doppler::GetBEnKinB(float BTh, bool kinflag) {
  float BEn = GetBEnKin(GetBThCoM(BTh, kinflag));

  if (usekinloss) {
    BEn -= GetBKinLoss(BEn, BTh);
  }

  return BEn;
}

/**
 * @brief Calculate beam energy for a given target scattering angle
 * using two-body kinematics calculations only.
 * 
 * This is used when target was detected, and not beam.
 * 
 * @param TTh theta angle of the target in laboratory frame (detected).
 * @param BTh theta angle of the beam in laboratory frame (calculated).
 * @param kinflag kinematics flag such that true is the backwards solution (i.e. CoM > 90 deg).
 * 
 * @return The calculated beam energy. 
 */
float doppler::GetBEnKinT(float TTh, float BTh, bool kinflag) {
  float BEn = GetTEnKin(GetTThCoM(TTh, kinflag));

  if (usekinloss) {
    BEn -= GetTKinLoss(BEn, BTh);
  }

  return BEn;
}

/**
 * @brief Calculate target energy for a given beam scattering angle
 * using two-body kinematics calculations only.
 * 
 * This is used when beam was detected, and not target.
 * 
 * @param BTh theta angle of the beam in laboratory frame (detected).
 * @param TTh theta angle of the target in laboratory frame (calculated).
 * @param kinflag kinematics flag such that true is the backwards solution (i.e. CoM > 90 deg).
 * 
 * @return The calculated beam energy using target angle.
 */
float doppler::GetTEnKinB(float BTh, float TTh, bool kinflag) {

  float TEn = GetTEnKin(GetBThCoM(BTh, kinflag));

  if (usekinloss) {
    TEn -= GetTKinLoss(TEn, TTh);
  }

  return TEn;
}

/**
 * @brief Calculate target energy for a given target scattering angle
 * using two-body kinematics calculations only.
 * 
 * This is used when target was detected.
 * 
 * @param TTh theta angle of the target in laboratory frame (detected).
 * @param kinflag kinematics flag such that true is the backwards solution (i.e. CoM > 90 deg).
 * 
 * @return The calculated target energy using target angle. 
 */
float doppler::GetTEnKinT(float TTh, bool kinflag) {
  float TEn = GetBEnKin(GetTThCoM(TTh, kinflag));

  if (usekinloss) {
    TEn -= GetBKinLoss(TEn, TTh);
  }

  return TEn;
}

float doppler::GammaAng(float PTh, float PPhi, float GTh, float GPhi) {

  /// Returns angle between particle and gamma	in radians

  double costheta =
      sin(PTh) * sin(GTh) * cos(PPhi - GPhi) + (cos(PTh) * cos(GTh));

  return TMath::ACos(costheta);
}

/**
 * @brief Calculates the beta factor after third order Taylor expansion.
 * 
 * @param Ek Kinetic energy.
 * @param m Mass.
 * 
 * @return Relativistic beta.
 */
float doppler::Beta(float Ek, float m) {
  double beta2 = -0.5 * m + TMath::Sqrt(m * (0.25 * m + 1.5 * Ek));
  beta2 /= 0.75 * m;

  return TMath::Sqrt(beta2);
}

/**
 * @brief Get gamma theta angle in radians.
 * 
 * @param cid Crystal ID.
 * @param sid Segment ID.
 * 
 * @return Gamma theta angle. 
 */
float doppler::GetGTh(int cid, int sid) {
  return gamma_theta[cid / 3][cid % 3][sid];
}

/**
 * @brief Get gamma phi angle in radians.
 * 
 * @param cid Crystal ID.
 * @param sid Segment ID.
 * 
 * @return Gamma phi angle. 
 */
float doppler::GetGPh(int cid, int sid) {
  return gamma_phi[cid / 3][cid % 3][sid];
}

/**
 * @brief Returns Doppler correction factor for given particle and gamma angular combination.
 * Velocity is calculated from detected particle energy, unless
 * the 'usekin' flag is set to true, in which case it uses the velocity
 * calculated from two-body kinematics
 * 
 * @param PEn Particle energy.
 * @param PTh Particle theta angle.
 * @param PPhi Particle phi angle.
 * @param GTh Gamma theta angle.
 * @param GPhi Gamma phi angle.
 * @param A Atomic mass.
 * 
 * @return Doppler correction factor. 
 */
float doppler::DC(float PEn, float PTh, float PPhi, float GTh, float GPhi, float A) {
  float beta = Beta(PEn, A * u_mass());
  float gamma = 1. / TMath::Sqrt(1. - beta * beta);
  float costheta = sin(PTh) * sin(GTh) * cos(PPhi - GPhi) + (cos(PTh) * cos(GTh));

  float corr = 1. - beta * costheta;
  corr *= gamma;

  return corr;
}

float doppler::DC_elec(float een, float PEn, float PTh, float PPhi, float GTh,
                       float GPhi, float A) {

  /// Returns Doppler correction factor for given particle and electron
  /// angular combination.  Factors in detected particle energy too

  float beta = TMath::Sqrt(2 * PEn / (A * u_mass()));
  float mass_e = 511.;
  float gamma = 1. / TMath::Sqrt(1. - beta * beta);
  float costheta =
      sin(PTh) * sin(GTh) * cos(PPhi - GPhi) + (cos(PTh) * cos(GTh));

  float energy = een + mass_e -
                 beta * costheta * TMath::Sqrt(een * een + 2. * mass_e * een);
  energy /= gamma;
  energy -= mass_e;

  return energy;
}

string doppler::convertInt(int number) {

  /// Convert an integer into a string

  stringstream ss;
  ss << number;
  return ss.str();
}

string doppler::convertFloat(float number) {

  /// Convert an float into a string

  stringstream ss;
  ss << number;
  return ss.str();
}
#endif

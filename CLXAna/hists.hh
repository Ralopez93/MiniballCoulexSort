#ifndef hists_hh
#define hists_hh

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TTree.h"
#include "TObject.h"

#include <iostream>
#include <string>
using namespace std;

// Headers for doppler
#ifndef doppler_hh
# include "doppler.hh"
#endif

using namespace std;

// Particle array length is 2 to accomodate beam and target.
#define P_ARR_LEN 2

// Gamma array length is currently set to 24. So far no entry has contained that many.
#define G_ARR_LEN 24

/// A class for making Coulex analysis histograms
/// The hists class is used to define all histograms for the analysis.
/// Crucially it defines a set of functions to be called, rather than
/// filling each histogram individually with every call in g_clx::Loop.
/// This has the advantage of not repeating fill commands and therefore
/// reducing the potential for copy/paste errors.
/// Please try to use these functions as much as is reasonably possible.

class hists {
public:

  TTree  *tree; // JP: tree output for quick analysis
  // tree branch variables
  int laser; // 0 for laser off, 1 for laser on
  int np;
  int run_nbr;    // For a lack of better solutions, going to tag each entry with the run number. Woo.
  double tdpp;    // Particle-Particle time difference in 25 ns timestamps.
  double time[P_ARR_LEN]; // Particle timestamp in 25 ns timestamps.
  int pid[P_ARR_LEN];
  int quad[P_ARR_LEN];
  int ring[P_ARR_LEN];
  int sect[P_ARR_LEN];
  double ep[P_ARR_LEN];
  double er[P_ARR_LEN];
  double thp[P_ARR_LEN];
  double php[P_ARR_LEN];
  double thr[P_ARR_LEN];
  double phr[P_ARR_LEN];
	double com[P_ARR_LEN];

  int ng;
  double td[G_ARR_LEN];  // Particle-gamma time difference in 25 ns timestamps
  double eg[G_ARR_LEN]; // calibrated, not dc'ed gamma energy in keV
  double ebg[G_ARR_LEN]; // calibrated and dc-ed gamma energy to beam kinematics in keV
  double etg[G_ARR_LEN]; // calibrated and dc-ed gamma energy to target kinematics in keV
  int clu[G_ARR_LEN];
  int cry[G_ARR_LEN]; // crystal ID from 0 to 23
  int seg[G_ARR_LEN]; // segment ID from 0 to 6 (0: core-only/ambiguous events)
  double thg[G_ARR_LEN]; // theta of gamma in lab frame
  double phg[G_ARR_LEN]; // phi of gamma in lab frame
  int tpg[G_ARR_LEN]; // smallest particle-gamma time difference (if 2 particles)
  double abg[G_ARR_LEN]; // angle difference between beam-like particle and gamma in degrees, 0-180
  double atg[G_ARR_LEN]; // angle difference between target-like particle and gamma in degrees, 0-180
  
	// Variables to be set in g_clx.C via Set_xxx functions
	float ppwin;
	int maxrecoil;
	int minrecoil;

	// Doppler instance
	doppler dc;

	// functions
	void Initialise( doppler dc_ );
	void Set_ppwin( float user_ppwin );
	void Set_maxrecoil( int user_maxrecoil );
	void Set_minrecoil( int user_minrecoil );

	// fill functions
  // fill Tree
  bool FillTree(float GEn, float GTh, float GPh, int GCluid, int GCid, int GSid, vector <float> GCor_GEn, vector <float> GCor_GTh,
		vector <float> GCor_GPh, vector <int> GCor_GCluID, vector <int> GCor_GCid, vector <int> GCor_GSid, vector <float> GCor_Gtd,
		vector <int> Laser, vector <float> PEn, vector<int> Pnf, vector<int> Pnb, vector<int> Psec,
		vector <int> Pquad, vector <float> Ptd, vector<double> Ptimes, int cur_run_nbr);
  
private:
  vector<int> laser_passed;
  vector<float> PEn_passed;
  vector<float> Pnf_passed;
  vector<float> Pnb_passed;
  vector<float> Pquad_passed;
  vector<float> Psec_passed;
  vector<float> Ptd_passed;
  vector<int> Ppid_passed;
  vector<float> PTheta_passed;
  vector<double> times_passed;

  void resetVar();
  bool doRoutine1P(float GEn, float GTh, float GPh, int GCluid, int GCid,
                   int GSid, vector<float> GCor_GEn, vector<float> GCor_GTh,
                   vector<float> GCor_GPh, vector<int> GCor_GCluID,
                   vector<int> GCor_GCid, vector<int> GCor_GSid,
                   vector<float> GCor_Gtd);
  bool doRoutine2P(float GEn, float GTh, float GPh, int GCluid, int GCid,
                   int GSid, vector<float> GCor_GEn, vector<float> GCor_GTh,
                   vector<float> GCor_GPh, vector<int> GCor_GCluID,
                   vector<int> GCor_GCid, vector<int> GCor_GSid,
                   vector<float> GCor_Gtd);
};
#endif

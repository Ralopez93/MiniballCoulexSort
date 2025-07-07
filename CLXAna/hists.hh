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

/// A class for making Coulex analysis histograms

/// The hists class is used to define all histograms for the analysis.
/// Crucially it defines a set of functions to be called, rather than
/// filling each histogram individually with every call in g_clx::Loop.
/// This has the advantage of not repeating fill commands and therefore
/// reducing the potential for copy/paste errors.
/// Please try to use these functions as much as is reasonably possible.

class hists {

	// Declare histos : 
	// B = Beam detection 
	// T = Target detection 
	// 2hB = both detected - Beam part
	// 2hT = both detected - Target part
	// p=prompt r=random

public:

  TTree  *tree; // JP: tree output for quick analysis
  // tree branch variables
  int laser; // 0 for laser off, 1 for laser on
  int np;
  double tdpp;
  int pid[2];
  int quad[2];
  int ring[2];
  int sect[2];
  double ep[2];
  double er[2];
  double thp[2];
  double php[2];
  double thr[2];
  double phr[2];
	
  int ng;
  double td[24];
  double eg[24]; // calibrated, not dc'ed gamma energy in keV
  double ebg[24]; // calibrated and dc-ed gamma energy to beam kinematics in keV
  double etg[24]; // calibrated and dc-ed gamma energy to target kinematics in keV
  int clu[24];
  int cry[24]; // crystal ID from 0 to 23
  int seg[24]; // segment ID from 0 to 6 (0: core-only/ambiguous events)
  double thg[24]; // theta of gamma in lab frame
  double phg[24]; // phi of gamma in lab frame
  int tpg[24]; // smallest particle-gamma time difference (if 2 particles)
  double abg[24]; // angle difference between beam-like particle and gamma in degrees, 0-180
  double atg[24]; // angle difference between target-like particle and gamma in degrees, 0-180
  
	// Variables to be set in g_clx.C via Set_xxx functions
	float ppwin;
	int maxrecoil;
	int minrecoil;

	// Array of cd angles for histogram bins
	double cd_angles[65];

	// Doppler instance
	doppler dc;

	// functions
	void Initialise( doppler dc_ );
	void Set_ppwin( float user_ppwin );
	void Set_maxrecoil( int user_maxrecoil );
	void Set_minrecoil( int user_minrecoil );

	// fill functions
  // fill Tree
  void FillTree(float GEn, float GTh, float GPh, int GCluid, int GCid, int GSid, vector <float> GCor_GEn, vector <float> GCor_GTh,
		vector <float> GCor_GPh, vector <int> GCor_GCluID, vector <int> GCor_GCid, vector <int> GCor_GSid,
		vector <float> GCor_Gtd,
		vector <int> Laser, vector <float> PEn, vector<int> Pnf, vector<int> Pnb, vector<int> Psec,
		vector <int> Pquad, vector <float> Ptd);
  
private:
  vector<int> laser_passed;
  vector<float> PEn_passed;
  vector<float> Pnf_passed;
  vector<float> Pnb_passed;
  vector<float> Pquad_passed;
  vector<float> Psec_passed;
  vector<float> Ptd_passed;
  vector<int> Ppid_passed;
	//ClassDef(hists,1);

  bool isGood2p(int quad_diff, float time_diff, float ppwin, int cut2);
};
#endif

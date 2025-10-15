#ifndef doppler_hh
#define doppler_hh

#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TCutG.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TObject.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <ctime>
using namespace std;

// Experimental definitions
//#ifndef ExpDefs_hh
//# include "ExpDefs.hh"
//#endif

// JP: update MB angles and distances for tuning
// Header for calibration
#ifndef __CALIBRATION_HH
# include "Calibration.hh"
#endif
// Header for MB geometry
#ifndef __MBGEOMETRY_HH__
# include "MBGeometry.hh"
#endif

typedef struct _Crystal {
  double theta;
  double phi;
} Crystal;

typedef struct _Cluster {
  std::vector<Crystal> crystals;
} Cluster;

/// A class for performing all aspects of the Doppler correction

/// The doppler class has a whole suite of functions for performing the
/// Doppler correction, including: calculating kinematics, energy loss,
/// particle gating, gamma-ray angles, particle angles, etc.

class doppler : public TObject {

public:

  void	ExpDefs( int Zb_, float Ab_, int Zt_, float At_, float Eb_, float Ex_, float thick_,
		 float depth_, float cddist_, float cdoffset_, float deadlayer_,
		 float contaminant_, float spededist_, TCutG *Bcut_, TCutG *Tcut_,
		 string srimfile_, bool usekin_, string calfile_ );
  int Cut(float PEn, float ring, float PTheta );
  int Cut_2p(float PEn1, float ring1, float PTheta1,
			       float PEn2, float ring2, float PTheta2 );
  bool	CutG_en2hit( float BEn, float TEn );
  int		GetZb();
  int		GetZt();
  float	GetAb();
  float	GetAt();
  float	GetCDOffset();
  float	GetCDDeadLayer();
  float	GetSpedeDist();
  float	GetPTh( float ring, int sector );
  float	GetPPhi( int quad, int seg, int sector );
  float	GetTTh( float BEn, float BTheta );
  float	GetBTh( float TTheta );
  float	GetQPhi( int quad, int seg, int sector );
  float	GetTEn( float BEn, float BTheta );
  float	GetBEn( float TEn, float TTheta );
  float	GetELoss( float Ei, float dist, int opt, string combo );
  float	GetBThLab( float CoM );
  float	GetTThLab( float CoM );
  float	GetBThCoM( float BTh, bool kinflag = false );
  float	GetTThCoM( float TTh, bool kinflag = false );
  float	GetBThLabT( float TTh, bool kinflag = false );
  float	GetTThLabB( float BTh, bool kinflag = false );
  float	GetBEnKin( float CoM );
  float	GetBEnKinB( float BTh, bool kinflag = false );
  float	GetBEnKinT( float TTh, bool kinflag = false );
  float	GetTEnKin( float CoM );
  float	GetTEnKinB( float BTh, bool kinflag = false );
  float	GetTEnKinT( float TTh, bool kinflag = false );
  float	GammaAng( float PTh, float PPhi, float GTh, float GPhi );
  float GetGTh(int cid, int sid); 
  float GetGPh(int cid, int sid); 
  float	Beta( float Ek, float m );
  float	DC( float PEn, float PTh, float PPhi, float GTh, float GPhi, float A );
  float	DC_elec( float een, float PEn, float PTh, float PPhi, float GTh, float GPhi, float A );
  bool	stoppingpowers( bool BT, bool TT, bool BS, bool TS, bool BC, bool TC );
  bool	stoppingpowers( string opt );
  void	reactionEnergy();
  void mbAngles(std::vector<Cluster> &clusters); // JP: re-define MB angles for tuning

  static string	convertInt( int number );
  static string	convertFloat( float number );
	
  inline float u_mass() { return 931494.0038; };
  inline bool	UseKin() { return usekin; };
	
private:
  static TRandom3	rand;
  TGraph *gSP[6];
  static string	gElName[110];
  static float	gates[64][3];
  float gamma_theta[8][3][7];
  float gamma_phi[8][3][7];
  int		Zb, Zt;
  int cluid, cid, sid; // for MB angle database
  float	Ab, At;
  float	Eb, Ex, thick, depth;
  float   Ereac; 			///< Energy at the reaction point (thick-depth)
  float	cddist, cdoffset;
  float	deadlayer;
  float	contaminant;
  float	spededist;
  float	bg_frac;
  TCutG	*Bcut, *Tcut;
  string  srimfile;
  bool	usekin;			///< flag to use two-body kinematics for velocity
  string calfile; // re-load MB geometry parameters
  Calibration *Cal;
  MBGeometry mbg;
  ClassDef(doppler,1);

};
#endif

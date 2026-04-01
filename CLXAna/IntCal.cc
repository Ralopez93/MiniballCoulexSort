#include "IntCal.hh"

#include <iostream>

#include <TFile.h>
#include <TMath.h>

using std::cout;
using std::endl;

IntCal::IntCal(string intcalfile_) {
  intcalfile = intcalfile_;
}

IntCal::~IntCal() {
  if (BCalGraph != nullptr) {
    free(BCalGraph);
    BCalGraph = nullptr;
  }
  
  if (TCalGraph != nullptr) {
    free(TCalGraph);
    BCalGraph = nullptr;
  }
}

bool IntCal::Setup() {

  size_t sep1 = intcalfile.find_first_of(":");
  size_t sep2 = intcalfile.find_last_of(":");

  if (sep1 <= 1 || sep2 <= 1 || sep1 > intcalfile.size() || sep2 > intcalfile.size()) {
    cout << "Format for the cutfile should be <cutfile.root>:<Bgraph>:<TGraph>\n";
    cout << "where <Bgraph> and <TGraph> are the TGraphG names of the beam and\n";
    cout << "target cuts in the root file, respectively." << endl;
    return false;
  }

  string str_file = intcalfile.substr(0, sep1);
  string str_bgraph = intcalfile.substr(sep1 + 1, sep2 - sep1 - 1);
  string str_tgraph = intcalfile.substr(sep2 + 1, intcalfile.size() - sep2 - 1);

  TFile *fgraph = new TFile(str_file.c_str());

  if (fgraph->IsZombie()) {
    cout << "Didn't open " << str_file << " correctly\n";
    cout << "Does it exist?\n";
    return false;
  }

  if (!fgraph->GetListOfKeys()->Contains(str_bgraph.c_str())) {

    cout << "Didn't find beam graph, "
         << str_bgraph
         << ", in "
         << str_file
         << endl;
    return false;
  }

  if (!fgraph->GetListOfKeys()->Contains(str_tgraph.c_str())) {
    cout << "Didn't find target graph, "
         << str_tgraph
         << ", in "
         << str_file
         << endl;
    return false;
  }

  BCalGraph = (TGraphErrors *)fgraph->Get(str_bgraph.c_str())->Clone();
  TCalGraph = (TGraphErrors *)fgraph->Get(str_tgraph.c_str())->Clone();

  fgraph->Close();

  free(fgraph);

  return true;
}

/**
 * @brief Get the internal calibration factor.
 * 
 * Factor for a given angle is defined as E(2B) - E(det),
 * where E(2B) is the estimated energy from 2-body kinematics and
 * E(det) is the mean detected energy (with attempted reconstruction).
 * 
 * @param BTh Beam scattering angle in lab frame.
 * 
 * @return Calibration factor.
 */
float IntCal::getBFactor(float BTh) {
  float angle = BTh * TMath::RadToDeg();
  return BCalGraph->Eval(angle);
}

/**
 * @brief Get the internal calibration factor.
 * 
 * Factor for a given angle is defined as E(2B) - E(det),
 * where E(2B) is the estimated energy from 2-body kinematics and
 * E(det) is the mean detected energy (with attempted reconstruction).
 * 
 * @param BTh Beam scattering angle in lab frame.
 * 
 * @return Calibration factor.
 */
float IntCal::getTFactor(float TTh) {
  float angle = TTh * TMath::RadToDeg();
  return TCalGraph->Eval(angle);
}
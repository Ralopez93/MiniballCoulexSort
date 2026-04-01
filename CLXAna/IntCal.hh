#ifndef INT_CAL_HH
#define INT_CAL_HH

#include <vector>
#include <string>

#include <TGraphErrors.h>

using std::pair;
using std::string;
using std::vector;

class IntCal {

public:

	IntCal(string intcalfile_);
	virtual ~IntCal();

  float getBFactor(float Th);
  float getTFactor(float Th);

  bool Setup();

private:

  TGraphErrors *BCalGraph = nullptr;
  TGraphErrors *TCalGraph = nullptr;

  string intcalfile;

};

#endif // INT_CAL_HH

// CLXAna: The main routine that calls g_clx
// Liam Gaffney (liam.gaffney@cern.ch) - 01/05/2017

#include "CLXAna.hh"

#ifndef __CINT__

void PrintInput() {

  cout << "Zb = " << Zb << endl;
  cout << "Ab = " << Ab << endl;
  cout << "Zt = " << Zt << endl;
  cout << "At = " << At << endl;
  cout << "Eb = " << Eb << " keV/u" << endl;
  cout << "Ex = " << Ex << " keV" << endl;
  cout << "thick = " << thick << " mg/cm2" << endl;
  cout << "depth = " << depth << " mg/cm2" << endl;
  cout << "cddist = " << cddist << " mm" << endl;
  cout << "cdoffset = " << cdoffset << " degrees" << endl;
  cout << "deadlayer = " << deadlayer << " mm" << endl;
  cout << "contaminant = " << contaminant << " mg/cm2" << endl;
  cout << "bg_frac = " << bg_frac << endl;
  cout << "srim = " << srim << endl;
  cout << "cutfile = " << cutfilename << endl;
  cout << "calfile = " << calfilename << endl;
  cout << "np_only = " << np_only << endl;
  if (print_angles)
    cout << "Priting crystal angles for each cluster." << endl;
  if (usekin)
    cout << "Using two-body kinematics for particle velocity" << endl;
  else
    cout << "Using detected particle energy for velocity calculation" << endl;

  if (usekinloss)
    cout << "Using SRIM to adjust velocity calculation from detected energy with two-body kinematics." << endl;

  if (emit_outside_targ)
    cout << "Assuming gamma emission taking place OUTSIDE target on average." << endl;
  else
    cout << "Assuming gamma emission taking place INSIDE target on average." << endl;

  cout << "\nOutputfile = " << outputfilename << endl << endl;

  return;
}

int main(int argc, char *argv[]) {

  CommandLineInterface *interface = new CommandLineInterface();

  interface->Add("-i", "Input file list", &inputfilenames);
  interface->Add("-o", "Output file name", &outputfilename);
  interface->Add("-c", "Configuration file", &configfilename);
  interface->Add("-cut", "Cutfile [cutfile.root:Bcut:Tcut]", &cutfilename);
  interface->Add("-Zb", "Zb", &Zb);
  interface->Add("-Ab", "Ab", &Ab);
  interface->Add("-Zt", "Zt", &Zt);
  interface->Add("-At", "At", &At);
  interface->Add("-Eb", "Eb", &Eb);
  interface->Add("-Ex", "Ex", &Ex);
  interface->Add("-thick", "Target thickness (mg/cm^2)", &thick);
  interface->Add("-depth", "Depth of interation in target (mg/cm^2)", &depth);
  interface->Add("-cddist", "Relative distance of CD and target (mm)", &cddist);
  interface->Add("-cdoffset", "Rotation of CD detector about phi from vertical (deg)", &cdoffset);
  interface->Add("-deadlayer", "Thickness of Si dead layer (mg/cm^2)", &deadlayer);
  interface->Add("-contaminant", "Thickness of contaminant layer on target (mg/cm^2)", &contaminant);
  interface->Add("-bg_frac", "Ratio of prompt and random for background subtraction", &bg_frac);
  interface->Add("-srim", "Directory containing the SRIM files", &srim);
  interface->Add("-emit_outside_targ", "Assume gamma emission outisde target.", &emit_outside_targ);
  interface->Add("-usekin", "Use two-body kinematics for particle velocity.", &usekin);
  interface->Add("-usekinloss", "Adjust velocity calculation from energy with two-body kinematics.", &usekinloss);
  interface->Add("-cal", "Calibration file", &calfilename);
  interface->Add("-print_angles", "Print crystal angles of each cluster", &print_angles);
  interface->Add("-cur_run_nbr", "Tag data with the given run number", &cur_run_nbr);
  interface->Add("-np_only", "Only sort entries with np number of particles, skipping rest.", &np_only);

  interface->CheckFlags(argc, argv);

  // Test if output file is there
  if (outputfilename.size() <= 0) {

    cout << "Did you specify your output file correctly?" << endl;
    return 0;
  }

  // Test if input files are there
  if (inputfilenames.size() <= 0) {

    cout << "Did you specify your input files correctly?" << endl;
    return 0;
  }

  if (calfilename.size() <= 0) {

    cout << "Did you specify your calibration file correctly?" << endl;
    return 0;
  }
  // Make chain for g_clx
  TChain *chain = new TChain("g_clx", "");
  for (unsigned int i = 0; i < inputfilenames.size(); i++) {

    chain->Add(inputfilenames[i].c_str());
  }

  // Initiate g_clx
  g_clx clx(chain);

  // Get the cut file
  if (cutfilename.size() > 0) {

    size_t sep1 = cutfilename.find_first_of(":");
    size_t sep2 = cutfilename.find_last_of(":");

    if (sep1 <= 1 || sep2 <= 1 || sep1 > cutfilename.size() || sep2 > cutfilename.size()) {

      cout << "Format for the cutfile should be <cutfile.root>:<Bcut>:<Tcut>\n";
      cout << "where <Bcut> and <Tcut> are the TCutG names of the beam and\n";
      cout << "target cuts in the root file, respectively." << endl;

      return 0;
    }

    string str_file = cutfilename.substr(0, sep1);
    string str_bcut = cutfilename.substr(sep1 + 1, sep2 - sep1 - 1);
    string str_tcut = cutfilename.substr(sep2 + 1, cutfilename.size() - sep2 - 1);

    TFile *fcut = new TFile(str_file.c_str());

    if (fcut->IsZombie()) {

      cout << "Didn't open " << str_file << " correctly\n";
      cout << "Does it exist?\n";

      return 0;
    }

    if (!fcut->GetListOfKeys()->Contains(str_bcut.c_str())) {

      cout << "Didn't find beam cut, " << str_bcut << ", in " << str_file
           << endl;
      return 0;
    }

    if (!fcut->GetListOfKeys()->Contains(str_tcut.c_str())) {

      cout << "Didn't find target cut, " << str_tcut << ", in " << str_file
           << endl;
      return 0;
    }

    clx.Bcut = (TCutG *)fcut->Get(str_bcut.c_str());
    clx.Tcut = (TCutG *)fcut->Get(str_tcut.c_str());

  } else {
    // if not cut file given, make empty cuts
    clx.Bcut = new TCutG();
    clx.Tcut = new TCutG();
  }

  // Test if we're using a config file (we overwrite the values if so)
  if (configfilename.size() > 0) {

    TEnv *config = new TEnv(configfilename.c_str());

    Zb = config->GetValue("Zb", -1);
    Ab = config->GetValue("Ab", -1.0);
    Zt = config->GetValue("Zt", -1);
    At = config->GetValue("At", -1.0);
    Eb = config->GetValue("Eb", -1.0);
    Ex = config->GetValue("Ex", -1.0);
    thick = config->GetValue("thick", -1.0);
    depth = config->GetValue("depth", -1.0);
    cddist = config->GetValue("cddist", -1.0);
    cdoffset = config->GetValue("cdoffset", 242.6);
    deadlayer = config->GetValue("deadlayer", 0.0007);
    contaminant = config->GetValue("contaminant", -1.0);
    bg_frac = config->GetValue("bg_frac", -1.0);
    srim = config->GetValue("srim", "./srim");
    emit_outside_targ = config->GetValue("emit_outside_targ", false);
    usekin = config->GetValue("usekin", false);
    usekinloss = config->GetValue("usekinloss", false);
    np_only = config->GetValue("np_only", 0);
    print_angles = config->GetValue("print_angles", false);
  }

  // Parameters are already read from the command line if not overwritten by config file.
  if (Zb > 0 && Zt > 0 && Ab > 0 && At > 0 && Eb > 0 && Ex > 0 && thick > 0 && depth > 0 && cddist > 0) {
    // Setup the experimental parameters
    // x.SetExpDefs( Zb, Ab, Zt, At, Eb, Ex, thick, depth, cddist, cdoffset );
    clx.Zb = Zb;
    clx.Ab = Ab;
    clx.Zt = Zt;
    clx.At = At;
    clx.Eb = Eb;
    clx.Ex = Ex;
    clx.thick = thick;
    clx.depth = depth;
    clx.cddist = cddist;
    clx.cdoffset = cdoffset;
    clx.deadlayer = deadlayer;
    clx.contaminant = contaminant;
    clx.bg_frac = bg_frac;
    clx.srim = srim;
    clx.emit_outside_targ = emit_outside_targ;
    clx.usekin = usekin;
    clx.usekinloss = usekinloss;
    clx.calfile = calfilename;
    clx.print_angles = print_angles;
    clx.cur_run_nbr = cur_run_nbr;
    clx.np_only = np_only;
    cout << "Input parameters:" << endl;
    PrintInput();
  } else {
    // In case something is missing, print out what we have and quit
    cout << "Some input is missing, please check:" << endl;
    PrintInput();
    cout << "Exiting..." << endl;

    return 0;
  }

  // Run sort
  cout << "Begin g_clx loop." << endl;
  clx.Loop(outputfilename);

  cout << "Finished." << endl;
  cout << "\a" << endl;

  return 0;
}
#endif

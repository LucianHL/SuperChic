// main11.cc is a part of the PYTHIA event generator.
// Copyright (C) 2023 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords:
//            Basic usage
//            LHE file

// This is a simple test program.
// It illustrates how Les Houches Event File input can be used in Pythia8.
// It uses the ttsample.lhe(.gz) input file, the latter only with 100 events.

// Modified for Superchic
#ifndef HEPMC2
#include "Pythia8Plugins/HepMC3.h"
#else
#include "Pythia8Plugins/HepMC2.h"
#endif
#include "Pythia8/Pythia.h"
using namespace Pythia8;
bool is_semi_exclusive_photon_initiated(int i) {
 return (( 48 <= i && i <=53 ) || ( 55 <= i && i <=83 ));

}
int main(int argc, char ** argv) {
  if (argc!=5) return 1;
 Pythia8::Pythia8ToHepMC topHepMC(argv[2]);
  // Generator. We here stick with default values, but changes
  // could be inserted with readString or readFile.
  Pythia pythia;

  // Initialize Les Houches Event File run. List initialization information.
  pythia.readString("Beams:frameType = 4");
  pythia.readString(std::string("Beams:LHEF = ")+argv[1]);
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("PartonLevel:Remnants = off");
  pythia.readString("Check:event = off");
  pythia.readString("LesHouches:matchInOut = off");
  
  if (std::string(argv[3]) != "dd" && std::string(argv[3]) != "sda" && std::string(argv[3]) != "sdb" && std::string(argv[3]) != "el") { return 7;}
    printf("Running in %s mode\n",argv[3]);
  
  int processnumber = atoi(argv[4]);
  if (is_semi_exclusive_photon_initiated(processnumber)) {
    pythia.readString("BeamRemnants:primordialKT = off");
    pythia.readString("PartonLevel:FSR = on");
    pythia.readString("SpaceShower:dipoleRecoil = on");
    pythia.readString("SpaceShower:pTmaxMatch = 2");
    pythia.readString("SpaceShower:QEDshowerByQ = off");
    pythia.readString("SpaceShower:pTdampMatch=1");
    if (std::string(argv[3]) == "dd") pythia.readString("BeamRemnants:unresolvedHadron = 0");
    if (std::string(argv[3]) == "sdb") pythia.readString("BeamRemnants:unresolvedHadron = 1");
    if (std::string(argv[3]) == "sda") pythia.readString("BeamRemnants:unresolvedHadron = 2");
    if (std::string(argv[3]) == "el") pythia.readString("BeamRemnants:unresolvedHadron = 3");
  }
  
//BeamRemnants:unresolvedHadron = 0 for double dissociation (dd), 1 for
//single dissociation (sdb), 2 for single dissociation (sda), 3 for elastic (el).

  pythia.init();
  size_t events=0;
  // Book histogram.
  Hist nCharged("charged particle multiplicity",100,-0.5,399.5);

  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;

  // Begin event loop; generate until none left in input file.
  while (iAbort < nAbort) {

    // Generate events, and check whether generation failed.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;
      ++iAbort;
      continue;
    }

    // Sum up final charged multiplicity and fill in histogram.
    int nChg = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
    if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
      ++nChg;
    nCharged.fill(nChg);
    topHepMC.writeNextEvent( pythia );
    events++;
  // End of event loop.
  }

  // Give statistics. Print histogram.
  pythia.stat();
  cout << nCharged;
  if (events==0) return 1;
  // Done.
  return 0;
}

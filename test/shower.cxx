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

std::vector<std::string> split(const std::string& str) {
    std::vector<std::string> tokens;
    std::string::size_type lastPos = str.find_first_not_of("\n", 0);
    std::string::size_type pos     = str.find_first_of("\n", lastPos);
    while (std::string::npos != pos || std::string::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of("\n", pos);
        pos = str.find_first_of("\n", lastPos);
    }
    return tokens;
}
std::string config_common =
R"""(Tune:pp = 5
Tune:ee = 3
TimeShower:pTminChgQ = 0.40000
TimeShower:pTmin = 0.40000
TimeShower:alphaSvalue = 0.13830
StringZ:rFactC = 1.00000
StringZ:rFactB = 0.67000
BeamRemnants:halfScaleForKT =  1.00000
BeamRemnants:primordialKThard = 2.00000
BeamRemnants:primordialKTsoft = 0.50000
ColourReconnection:range =  1.50000
Diffraction:largeMassSuppress =  2.00000
MultipartonInteractions:alphaSvalue = 0.13500
MultipartonInteractions:ecmPow = 0.19000
MultipartonInteractions:ecmRef = 1800.000
MultipartonInteractions:expPow = 2.00000
MultipartonInteractions:pT0Ref = 2.08500
Next:numberShowEvent = 0
PDF:pSet = LHAPDF6:cteq6l1
SigmaProcess:alphaSvalue = 0.13500
SpaceShower:alphaSvalue = 0.13700
SpaceShower:ecmRef = 1800.000
StringFlav:etaSup = 0.63000
StringFlav:mesonBvector = 3.00000
StringFlav:mesonCvector = 1.06000
StringFlav:mesonSvector = 0.72500
StringFlav:mesonUDvector = 0.62000
StringFlav:popcornSpair = 0.50000
StringFlav:probQQ1toQQ0 = 0.0270000
StringFlav:probQQtoQ = 0.0900000
StringFlav:probSQtoQQ = 1.00000
StringFlav:probStoUD = 0.19000
StringPT:sigma = 0.30400
StringZ:aExtraDiquark = 0.50000
StringZ:aLund = 0.30000
StringZ:bLund = 0.80000
)""";
std::string config_dd_pp = R"""(
PartonLevel:MPI = off
SpaceShower:dipoleRecoil = on
SpaceShower:pTmaxMatch = 2
SpaceShower:QEDshowerByQ = on
BeamRemnants:primordialKT = off
BeamRemnants:unresolvedHadron = 0
SpaceShower:pTdampMatch=1
PDF:pSet = LHAPDF6:MSHT20qed_nnlo
)""";

std::string config_el_pp = R"""(
PartonLevel:MPI = off
SpaceShower:pTmaxMatch = 2
BeamRemnants:primordialKT = off
BeamRemnants:unresolvedHadron = 3
SpaceShower:pTdampMatch=1
PartonLevel:ISR = off
LesHouches:matchInOut = off
PartonLevel:Remnants = off
)""";

std::string config_ds_pp = R"""(
PartonLevel:MPI = off
SpaceShower:dipoleRecoil = on
SpaceShower:pTmaxMatch = 2
BeamRemnants:primordialKT = off
SpaceShower:QEDshowerByQ = off
BeamRemnants:unresolvedHadron = 1
SpaceShower:pTdampMatch = 1
PDF:pSet = LHAPDF6:MSHT20qed_nnlo
)""";

std::string config_sd_pp = R"""(
PartonLevel:MPI = off
SpaceShower:dipoleRecoil = on
SpaceShower:pTmaxMatch = 2
BeamRemnants:primordialKT = off
SpaceShower:QEDshowerByQ = off
BeamRemnants:unresolvedHadron = 2
SpaceShower:pTdampMatch = 1
PDF:pSet = LHAPDF6:MSHT20qed_nnlo
)""";

class config {
public:
config(){}
void readString(const std::string& str) {m_strings.push_back(str);} 
void apply(Pythia& p)const { for ( auto s: m_strings ) p.readString(s);  }
void dump(const std::string& outname) const { 
  std::ofstream out(outname.c_str());
  for ( auto s: m_strings ) out<<s;  
  out.close();
}
private:
std::vector<std::string> m_strings;
};
/*
bool is_semi_exclusive_photon_initiated(int i) {
 return (( 48 <= i && i <=53 ) || ( 55 <= i && i <=83 ));

}*/
int main(int argc, char ** argv) {
  if (argc!=6) return 1;
  
  std::vector<std::string> c_c = split(config_common);
  std::vector<std::string> c_sd_pp = split(config_sd_pp);
  std::vector<std::string> c_ds_pp = split(config_ds_pp);
  std::vector<std::string> c_el_pp = split(config_el_pp);
  std::vector<std::string> c_dd_pp = split(config_dd_pp);
  
 Pythia8::Pythia8ToHepMC topHepMC(argv[2]);
  // Generator. We here stick with default values, but changes
  // could be inserted with readString or readFile.
  Pythia pythia;
  config pythiaconfig;
  // Initialize Les Houches Event File run. List initialization information.
  pythiaconfig.readString("Beams:frameType = 4");
  pythiaconfig.readString(std::string("Beams:LHEF = ")+argv[1]);
  if ( std::string(argv[5]) != "dummy" ) {

    std::string type = std::string(argv[3]);
    int processnumber = atoi(argv[4]);

    if (
      type != "dd_pp" && type != "sda_pp" && type != "sdb_pp" && type != "sd_pp" && type != "el_pp"
     && type != "dd_AA" && type != "sda_AA" && type != "sdb_AA" && type != "sd_AA" && type != "el_AA"
     && type != "dd_pA" && type != "sda_pA" && type != "sdb_pA" && type != "sd_pA" && type != "el_pA"
      && type != "dd_ee" && type != "sda_ee" && type != "sdb_ee" && type != "sd_ee" && type != "el_ee"
    ) { printf("Bad argument->%s<-\n",argv[3]); return 7;}
      printf("Running in %s mode\n",argv[3]);
      std::string beam ="pp";
      beam[0]=type[type.size()-2];
      beam[1]=type[type.size()-1];
      if (beam=="AA" || beam=="pA") { 
		  pythiaconfig.readString("Check:beams = off");
    }      
    if (type[type.size()-1] == 'A' ) type[type.size()-1] ='p';
    if (type[type.size()-2] == 'A' ) type[type.size()-2] ='p';
    bool showerconfigured=false;
  
    if (!showerconfigured  && (type == "dd_pp" || type == "sdb_pp" || type == "sda_pp" ||  type == "sd_pp" || type == "el_pp" ) )
    {
    for ( auto s: c_c) pythiaconfig.readString(s);
    pythiaconfig.readString("PartonLevel:MPI = off");
    pythiaconfig.readString("Check:event = off");
    pythiaconfig.readString("LesHouches:matchInOut = off");


    if (type == "dd_pp") for ( auto s: c_dd_pp) pythiaconfig.readString(s);
    if (type == "sdb_pp") for ( auto s: c_ds_pp) pythiaconfig.readString(s);
    if (type == "sda_pp" ) for ( auto s: c_sd_pp) pythiaconfig.readString(s);
    if (type == "el_pp") for ( auto s: c_el_pp) pythiaconfig.readString(s);
    showerconfigured=true;
    }
  
  

    if (!showerconfigured && (type=="el_ee"  || type=="dd_ee" ||type=="sda_ee"||type=="sdb_ee" ||type=="sd_ee")) {
      for ( auto s: c_c) pythiaconfig.readString(s);
      pythiaconfig.readString("PartonLevel:ISR = off");
      pythiaconfig.readString("PartonLevel:MPI = off");
      pythiaconfig.readString("PartonLevel:Remnants = off");
      pythiaconfig.readString("Check:event = off");
      pythiaconfig.readString("LesHouches:matchInOut = off");
      showerconfigured=true;
    }
    
  }
  pythiaconfig.apply(pythia);
  pythiaconfig.dump(std::string(argv[2])+"_config");
  
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

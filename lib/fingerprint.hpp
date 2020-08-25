
#ifndef AtomFingerprint_hpp
#define AtomFingerprint_hpp

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <sstream>
#include <math.h>
#include <algorithm>

#include "Coor.hpp"
#include "Atom.hpp"
#include "Analyze.hpp"
#include "Molecule.hpp"
#include "Mol2.hpp"
#include "Trajectory.hpp"
#include "Select.hpp"

#define PI 3.14159265

using namespace std;

// atom fingerprint
namespace AFP{
  typedef map<string, vector<double> > AtomDev;
  typedef map<string, map<int, double > > AtomFP; // {key: neighboring atom type; value: {key: eta, value: fingerprint vector} }
  bool getNeighDists(AtomDev& neighDists, Coor coor, Molecule* mol, double cutoff, vector<string> neighTypes);
  bool printNeighDists(AtomDev neighDists, bool dict=true);
  bool getFP(AtomFP& atomFP, Coor coor, Molecule* mol,
             double cutoff, vector<int> etas, vector<string> neighTypes);
  bool printFP(AtomFP atomFP);
}

// molecule fingerprint
namespace MFP{
  typedef map<Atom*,  AFP::AtomDev > MolDev;
  typedef map<Atom*, AFP::AtomFP > MolFP;
  typedef map<string, AFP::AtomFP > PoseMolFP;
  bool getNeighDists(MolDev& neighDists, vector<Atom*> atoms, Molecule* mol, double cutoff, vector<string> neighTypes);
  bool printNeighDists(MolDev neighDists, bool dict=true);
  bool getFP(MolFP& molFP, vector<Atom*> atoms, Molecule* mol,
             double cutoff, vector<int> etas, vector<string> neighTypes);
  bool printFP(MolFP molFP);
  bool assignMol2(vector<Atom*> atoms, Mol2* mol2);
  bool convert2PoseFP(PoseMolFP& poseFP, MolFP molFP);
  bool printPoseFP(PoseMolFP poseFP);

}

// WORK IN PROGRESS
// trajectory fingerprint (pose fingerprint, RNAPosers)
namespace TFP{
  // typedef map<int, MFP::MolFP > TrajFP; // map frameid to MolFP
  typedef map<int, MFP::PoseMolFP > PoseFP;
  // bool getTrajFP(TrajFP& trajFP, vector<Atom*> atoms, Trajectory* ftrjin, ifstream trjin,
  //           double cutoff, vector<int> etas, vector<string> neighTypes);
  bool getPoseFP(PoseFP& poseFP, vector<Atom*> atoms, Trajectory* ftrjin, ifstream trjin,
             double cutoff, vector<int> etas, vector<string> neighTypes,
             int start, int stop, int skip, bool startFlag);
  bool printPoseFP(PoseFP poseFP);
  bool getMedian(MFP::PoseMolFP& medianFP, PoseFP poseFP);
  bool normalizeFP(PoseFP& poseFP);
}


namespace FileProcess{
  Molecule* loadMol2(string mol2File);
  Trajectory* loadTraj(ifstream& trjin, string trajFile, Molecule* mol,
                       int skip, bool startFlag);
}
#endif

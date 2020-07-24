
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
#include "Select.hpp"

#define PI 3.14159265

using namespace std;

namespace AFP{
  typedef map<string, vector<double> > AtomDev;
  typedef map<string, map<int, double > > AtomFP;
  bool getNeighDists(AtomDev& neighDists, Atom* atom, Molecule* mol, double cutoff, vector<string> neighTypes);
  bool printNeighDists(AtomDev neighDists);
  bool getFP(AtomFP& atomFP, Atom* atom, Molecule* mol,
             double cutoff, vector<int> etas, vector<string> neighTypes);
  bool printFP(AtomFP atomFP);
}

namespace MFP{
  typedef map<Atom*,  AFP::AtomDev > MolDev;
  typedef map<Atom*, AFP::AtomFP > MolFP;
  bool getNeighDists(MolDev& neighDists, vector<Atom*> atoms, Molecule* mol, double cutoff, vector<string> neighTypes);
  bool printNeighDists(MolDev neighDists);
  bool getFP(MolFP& molFP, vector<Atom*> atoms, Molecule* mol,
             double cutoff, vector<int> etas, vector<string> neighTypes);
  bool printFP(MolFP molFP);
}

#endif

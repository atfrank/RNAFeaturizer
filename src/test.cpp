/*
  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.

  Author: Jingru Xie and Aaron T. Frank
*/

#include "Molecule.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Misc.hpp"
#include "Analyze.hpp"
#include "Trajectory.hpp"
#include "fingerprint.hpp"


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <time.h> // keep track of processing time of program

using namespace std;

void usage(){
  std::cerr << "====================================================" << std::endl;
  std::cerr << "For Test Atomic Fingerprint" << std::endl;
  std::cerr << "====================================================" << std::endl;
  exit(0);
}

int main (int argc, char **argv){
  int i;
  clock_t t_start, t_end;

  t_start = clock();

  std::string currArg;

  vector<string> pdbs;

  pdbs.clear();

  for (i = 1; i < argc; i++){
    currArg = argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0)
    {
      usage();
    }
    else if (currArg.compare(0,1,"-") == 0)
    {
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      pdbs.push_back(currArg);
    }
  }
  if (pdbs.size() == 0)
  {
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

  // initialize vars
  Molecule* mol=Molecule::readPDB(pdbs.at(0));
  mol->selAll();
  Atom* atom = mol->getAtom(1);
  double cutoff = 20;
  static const int arr[] = {2, 4, 6, 8, 16, 32};
  vector<int> etas (arr, arr + sizeof(arr) / sizeof(arr[0]) );
  static const string arr1[] = {":.C1'", ":.C4", ":.C8", ":.O1P", ":.O2P", ":.N4"};
  // static const string arr1[] = {":.heavy"};
  vector<string> neighTypes (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );

  // get dists
  // AFP::AtomDev neighDists;
  // AFP::getNeighDists(neighDists, atom, mol, cutoff, neighTypes);
  // AFP::printNeighDists(neighDists);
//
  // // get fingerprints
  // AFP::AtomFP atomFP;
  // AFP::getFP(atomFP, atom, mol, cutoff, etas, neighTypes);
  // AFP::printFP(atomFP);

  // molecular - get dists
  vector<Atom*> atoms = Select::makeSelVec(mol, ":.C8");
  MFP::MolDev dists;
  MFP::getNeighDists(dists, atoms, mol, cutoff, neighTypes);
  MFP::printNeighDists(dists);

  // molecular - get fingerprints
  MFP::MolFP molFP;
  MFP::getFP(molFP, atoms, mol, cutoff, etas, neighTypes);
  MFP::printFP(molFP);


  t_end = clock();
  float diff ((float)t_end - (float)t_start);
  cout << "Processing time: " << diff/CLOCKS_PER_SEC << " seconds." << endl << endl << endl;

  return 0;
}

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

  double x=0.0;
  double y=0.0;
  double z=0.0;
  vector<string> pdbs;


  for (i = 1; i < argc; i++){
    currArg = argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0)
    {
      usage();
    }
    else if (currArg.compare("-x") == 0){
      currArg=argv[++i];
      stringstream(currArg) >> x;
    }
    else if (currArg.compare("-y") == 0){
      currArg=argv[++i];
      stringstream(currArg) >> x;
    }
    else if (currArg.compare("-z") == 0){
      currArg=argv[++i];
      stringstream(currArg) >> x;
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
  Coor coor(x, y, z);
  Molecule* mol=Molecule::readPDB(pdbs.at(0));
  mol->selAll();
  double cutoff = 7.0;
  static const int arr[] = {2, 4, 8, 16};
  vector<int> etas (arr, arr + sizeof(arr) / sizeof(arr[0]) );
  static const string arr1[] = {":.C1'", ":.C4", ":.C8", ":.O1P", ":.O2P", ":.N4"};
  // static const string arr1[] = {":ADE.C1'", ":ADE.C2", ":ADE.C2'", ":ADE.C3'", ":ADE.C4", ":ADE.C4'", ":ADE.C5", ":ADE.C5'", ":ADE.C6", ":ADE.C8", ":ADE.N1", ":ADE.N3", ":ADE.N6", ":ADE.N7", ":ADE.N9", ":ADE.O2'", ":ADE.O3'", ":ADE.O4'", ":ADE.O5'", ":ADE.OP1", ":ADE.OP2", ":ADE.P", ":CYT.C1'", ":CYT.C2", ":CYT.C2'", ":CYT.C3'", ":CYT.C4", ":CYT.C4'", ":CYT.C5", ":CYT.C5'", ":CYT.C6", ":CYT.N1", ":CYT.N3", ":CYT.N4", ":CYT.O2", ":CYT.O2'", ":CYT.O3'", ":CYT.O4'", ":CYT.O5'", ":CYT.OP1", ":CYT.OP2", ":CYT.P", ":GUA.C1'", ":GUA.C2", ":GUA.C2'", ":GUA.C3'", ":GUA.C4", ":GUA.C4'", ":GUA.C5", ":GUA.C5'", ":GUA.C6", ":GUA.C8", ":GUA.N1", ":GUA.N2", ":GUA.N3", ":GUA.N7", ":GUA.N9", ":GUA.O2'", ":GUA.O3'", ":GUA.O4'", ":GUA.O5'", ":GUA.O6", ":GUA.OP1", ":GUA.OP2", ":GUA.P", ":URA.C1'", ":URA.C2", ":URA.C2'", ":URA.C3'", ":URA.C4", ":URA.C4'", ":URA.C5", ":URA.C5'", ":URA.C6", ":URA.N1", ":URA.N3", ":URA.O2", ":URA.O2'", ":URA.O3'", ":URA.O4", ":URA.O4'", ":URA.O5'", ":URA.OP1", ":URA.OP2", ":URA.P"};
  // static const string arr1[] = {":.heavy"};
  vector<string> neighTypes (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );

  // get dists
  AFP::AtomDev neighDists;
  AFP::getNeighDists(neighDists, coor, mol, cutoff, neighTypes);
  AFP::printNeighDists(neighDists);
//
  // // get fingerprints
  AFP::AtomFP atomFP;
  AFP::getFP(atomFP, coor, mol, cutoff, etas, neighTypes);
  AFP::printFP(atomFP);


  t_end = clock();
  float diff ((float)t_end - (float)t_start);
  cout << "Processing time: " << diff/CLOCKS_PER_SEC << " seconds." << endl;

  return 0;
}

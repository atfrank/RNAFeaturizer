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
#include "Grid.hpp"


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

void set_constrants(unordered_map<string, array<double, 2> >& dConstraints){
  array<double, 2> constr = {1.7, 7};
  dConstraints[":.N"] = constr;
  constr = {1.6, 6.2};
  dConstraints[":.OP1"] = constr;
  constr = {1.5, 5.25};
  dConstraints[":.OP2"] = constr;
  constr = {1.65, 4.4};
  dConstraints["!:.N+OP1+OP2"] = constr;
}

int main (int argc, char **argv){
  int i;
  clock_t t_start, t_end;

  t_start = clock();

  std::string currArg;

  vector<string> pdbs;
  string mol2file;

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
  double threshold = 8.0;
  double spacing = 2.8;
  // getGrids
  vector <Coor> grids;
  // get grids for molecule
  unordered_map<string, array<double, 2> > dConstraints;
  set_constrants(dConstraints);
  GRID::getGridsMol(grids, mol, threshold, spacing, dConstraints);
  // featurize

  // time analysis
  t_end = clock();
  float diff ((float)t_end - (float)t_start);
  cout << "Processing time: " << diff/CLOCKS_PER_SEC << " seconds." << endl;

  return 0;
}

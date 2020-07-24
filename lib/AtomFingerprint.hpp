//
//  AtomFingerprint.hpp
//
//
//  Created by Jingru Xie on 9/29/17.
//
//
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

class Molecule;

class AtomFingerprint{
public:
  AtomFingerprint();//default constructor
  AtomFingerprint(Molecule *mol=NULL);//value ctor
  bool getNeighDists();
  bool printNeighDists();
  bool getFP();
  bool printFP();

private:
  Atom* atom;
  Molecule* mol;
  double cutoff;
  vector<int> etas;
  vector<string> neighTypes;
  map<string, vector<double> > neighDists;
  map<string, map<int, double > > atomFP;

};


#endif /* AtomFingerprint_hpp */

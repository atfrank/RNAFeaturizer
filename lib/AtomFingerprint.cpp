//  AtomFingerprint.cpp
//  Created by J.X on 9/29/17
#include "AtomFingerprint.hpp"


// ctor
AtomFingerprint::AtomFingerprint(Molecule *mol){
  this->mol = mol;
  mol->selAll();
  atom = mol->getAtom(1);
  scalar = 1;
  cutoff = 20;
  etas.push_back(2);
  etas.push_back(4);
  etas.push_back(8);
  etas.push_back(16);
  etas.push_back(32);
  static const string arr[] = {":.C1'", ":.C4", ":.C8", ":.O1P", ":.O2P", ":.N4"};
  vector<string> temp (arr, arr + sizeof(arr) / sizeof(arr[0]) );
  neighTypes = temp;
}


bool AtomFingerprint::getNeighDists(){
  double dist;
  Coor coor;
  Coor coor_iter;
  string nType;
  Atom* aj;

  coor = atom->getCoor();
  // loop over neighbor atom types
  for (unsigned int k=0; k < neighTypes.size(); k ++){
    nType = neighTypes[k];
    neighDists.insert(make_pair(nType, vector<double>()));
    // mymap.insert(pair<int,vector<int> >(10, vector<int>()));
    // mol->select(nType, false);
    vector<Atom*> neighbors = Select::makeSelVec(mol, nType);
    // loop over each atom of the type
    for (unsigned int j=0; j < neighbors.size(); j ++){
      aj = neighbors[j];
      coor_iter = aj->getCoor();
      // get distance and keep distance if dist < cutoff
      dist = Analyze::distance(coor, coor_iter);
      if (!(coor == coor_iter) && (dist < cutoff)){
        neighDists[nType].push_back(dist);
      }
    }
  }
  return 0;
}

bool AtomFingerprint::printNeighDists(){
  vector<double> dists;
  string nType;
  for (unsigned int k=0; k < neighTypes.size(); k ++){
    nType = neighTypes[k];
    cout << nType << ":";
    dists = neighDists[nType];
    for (unsigned int j=0; j < dists.size(); j ++){
      cout << dists[j] << " ";
    }
    cout << endl;
  }
  return 0;
}


bool AtomFingerprint::getFP(){
  string nType;
  int eta;
  vector<double> dists;
  this->getNeighDists();
  for (unsigned int k=0; k < neighTypes.size(); k ++){
    nType = neighTypes[k];
    atomFP.insert(make_pair(nType, map<int, double >() ));
    dists = neighDists[nType];
    for (unsigned int i=0; i < etas.size(); i ++){
      eta = etas[i];
      atomFP[nType].insert(make_pair(eta, 0.));
      dists = neighDists[nType];
      for (unsigned int j=0; j < dists.size(); j ++){
        double dist = dists[j];
        atomFP[nType][eta] += exp((- dist/eta) * (dist/eta)) * (cos(PI * dist/ cutoff) + 1)/2;
      }
    }
  }
  return 0;
}

bool AtomFingerprint::printFP(){
  vector<double> dists;
  string nType;
  int eta;
  for (unsigned int k=0; k < neighTypes.size(); k ++){
    nType = neighTypes[k];
    for (unsigned int i=0; i < etas.size(); i ++){
      eta = etas[i];
      cout << nType << "_" << eta << "\t:" << atomFP[nType][eta] << endl;
    }
  }
  return 0;
}

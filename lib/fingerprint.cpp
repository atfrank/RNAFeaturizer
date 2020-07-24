#include "fingerprint.hpp"

namespace AFP{

  bool getNeighDists(AtomDev& neighDists, Atom* atom, Molecule* mol, double cutoff, vector<string> neighTypes){
    double dist;
    Coor coor;
    Coor coor_iter;
    string nType;
    Atom* aj;

    neighDists.clear();

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

  bool printNeighDists(AtomDev neighDists){
    vector<double> dists;
    string nType;

    for (AtomDev::iterator iter = neighDists.begin(); iter != neighDists.end(); ++iter){
      nType = iter->first;
    // for (unsigned int k=0; k < neighTypes.size(); k ++){
    //  nType = neighTypes[k];
      cout << nType << ":";
      dists = neighDists[nType];
      for (unsigned int j=0; j < dists.size(); j ++){
        cout << dists[j] << " ";
      }
      cout << endl;
    }
    return 0;
  }


  bool getFP(AtomFP& atomFP, Atom* atom, Molecule* mol,
             double cutoff, vector<int> etas, vector<string> neighTypes){
    string nType;
    int eta;
    vector<double> dists;
    AtomDev neighDists;

    atomFP.clear();
    // get distances
    getNeighDists(neighDists, atom, mol, cutoff, neighTypes);
    // convert distances into fingerprints
    for (AtomDev::iterator iter = neighDists.begin(); iter != neighDists.end(); ++iter){
      nType = iter->first;
      atomFP.insert(make_pair(nType, map<int, double >() ));
      dists = neighDists[nType];
      // multiple eta handling
      for (unsigned int i=0; i < etas.size(); i ++){
        eta = etas[i];
        atomFP[nType].insert(make_pair(eta, 0.));
        dists = neighDists[nType];
        // sum over all neighboring atoms of the same type
        for (unsigned int j=0; j < dists.size(); j ++){
          double dist = dists[j];
          atomFP[nType][eta] += exp((- dist/eta) * (dist/eta)) * (cos(PI * dist/ cutoff) + 1)/2;
        }
      }
    }
    return 0;
  }

  bool printFP(AtomFP atomFP){
    vector<double> dists;
    string nType;
    int eta;
    for (AtomFP::iterator iter = atomFP.begin(); iter != atomFP.end(); ++iter){
      nType = iter->first;
      map<int, double > fp = iter->second;
      for (map<int, double >::iterator iter2 = fp.begin(); iter2 != fp.end(); ++iter2){
        eta = iter2->first;
        cout << nType << "_" << eta << "\t:" << iter2->second << endl;
      }
    }
    return 0;
  }
}



namespace MFP{
  bool getNeighDists(MolDev& neighDists, vector<Atom*> atoms, Molecule* mol, double cutoff, vector<string> neighTypes){
    AFP::AtomDev atomDists;
    for (unsigned int i = 0; i < atoms.size(); i++){
      Atom* atom = atoms[i];
    // for (vector<Atom*>::iterator atom = atoms.begin(); atom != atoms.end(); ++atom){
      AFP::getNeighDists(atomDists, atom, mol, cutoff, neighTypes);
      neighDists[atom] = atomDists;
    }
    return 0;
  }

  bool printNeighDists(MolDev neighDists){
    for (MolDev::iterator iter = neighDists.begin(); iter != neighDists.end(); ++iter){
      Atom* atom = iter->first;
      cout << atom->getAtmNum() << " " << atom->getAtmName() << ":" << endl;
      AFP::printNeighDists(iter->second);
    }
    return 0;
  }

  bool getFP(MolFP& molFP, vector<Atom*> atoms, Molecule* mol, double cutoff, vector<int> etas,
             vector<string> neighTypes){
    AFP::AtomFP atomFP;
    for (unsigned int i = 0; i < atoms.size(); i++){
      Atom* atom = atoms[i];
    // for (vector<Atom*>::iterator atom = atoms.begin(); atom != atoms.end(); ++atom){
      AFP::getFP(atomFP, atom, mol, cutoff, etas, neighTypes);
      molFP[atom] = atomFP;
    }
    return 0;
  }

  bool printFP(MolFP molFP){
    for (MolFP::iterator iter = molFP.begin(); iter != molFP.end(); ++iter){
      Atom* atom = iter->first;
      cout << atom->getAtmNum() << " " << atom->getAtmName() << ":" << endl;
      AFP::printFP(iter->second);
    }
    return 0;
  }

}

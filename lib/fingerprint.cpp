#include "fingerprint.hpp"

// Helper functions
// get median in-place
bool getMedian(vector<double>& fp){
  double median;
  sort(fp.begin(), fp.end());
  int veclen = fp.size();
  // find median
  if (veclen % 2 == 0){
    median = (fp[veclen / 2 - 1] + fp[veclen / 2]) / 2;
  }
  else{
    median = fp[veclen / 2];
  }
  if (median != 0){
    for (unsigned int j = 0; j < veclen; ++j){
      // divide the whole column in the original vector matrix by column median
      fp[j] = fp[j] / median;
    }
    return 0;
  }
  else return 1;
}

// get list of map keys
template <class K, class V>
vector<K> getKeys(map<K, V> mymap){
  vector<K> keys;
  for (typename map<K, V>::iterator it = mymap.begin(); it != mymap.end(); ++it)
    keys.push_back(it->first);
  return keys;
}

// an alternative
// https://stackoverflow.com/a/2389291
template<class map_type>
class key_iterator : public map_type::iterator
{
public:
    typedef typename map_type::iterator map_iterator;
    typedef typename map_iterator::value_type::first_type key_type;

    key_iterator(const map_iterator& other) : map_type::iterator(other) {} ;

    key_type& operator *()
    {
        return map_type::iterator::operator*().first;
    }
};

namespace AFP{

  bool getNeighDists(AtomDev& neighDists, Coor coor, Molecule* mol, double cutoff, vector<string> neighTypes){
    double dist;
    Coor coor_iter;
    string nType;
    Atom* aj;

    neighDists.clear();

    // coor = atom->getCoor();
    // loop over neighbor atom types
    for (unsigned int k=0; k < neighTypes.size(); k ++){
      nType = neighTypes[k];
      neighDists.insert(make_pair(nType, vector<double>()));
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

  bool printNeighDists(AtomDev neighDists, bool dict){
    vector<double> dists;
    string nType;

    if (dict){
      cout << "{";
      for (AtomDev::iterator iter = neighDists.begin(); iter != neighDists.end(); ++iter){
        nType = iter->first;
        cout << '"' << nType << '"' << ": [";
        dists = neighDists[nType];
        for (unsigned int j=0; j < dists.size(); j ++){
          cout << dists[j] << ", ";
        }
        cout << "], ";
      }
      cout << "}" ;
    }
    else{
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
    }
    return 0;
  }

  bool getFP(AtomFP& atomFP, Coor coor, Molecule* mol,
             double cutoff, vector<int> etas, vector<string> neighTypes){
    string nType;
    int eta;
    vector<double> dists;
    AtomDev neighDists;

    atomFP.clear();
    // get distances
    getNeighDists(neighDists, coor, mol, cutoff, neighTypes);
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
      Coor coor = atoms[i]->getCoor();
    // for (vector<Atom*>::iterator atom = atoms.begin(); atom != atoms.end(); ++atom){
      AFP::getNeighDists(atomDists, coor, mol, cutoff, neighTypes);
      neighDists[atom] = atomDists;
    }
    return 0;
  }

  bool printNeighDists(MolDev neighDists, bool dict){
    if (dict) cout << "{";
    for (MolDev::iterator iter = neighDists.begin(); iter != neighDists.end(); ++iter){
      Atom* atom = iter->first;
      cout << atom->getAtmNum() << ":";
      AFP::printNeighDists(iter->second, dict);
      if (dict) cout << ",";
    }
    if (dict) cout << "}" << endl;
    return 0;
  }

  bool getFP(MolFP& molFP, vector<Atom*> atoms, Molecule* mol, double cutoff, vector<int> etas,
             vector<string> neighTypes){
    AFP::AtomFP atomFP;
    for (unsigned int i = 0; i < atoms.size(); i++){
      Atom* atom = atoms[i];
      Coor coor = atoms[i]->getCoor();
    // for (vector<Atom*>::iterator atom = atoms.begin(); atom != atoms.end(); ++atom){
      AFP::getFP(atomFP, coor, mol, cutoff, etas, neighTypes);
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

  bool assignMol2(vector<Atom*> atoms, Molecule* mol2){
    Atom* a1;
    Atom* a2;
    for (int k = 0; k < mol2->getNAtom(); k ++){
      a1 = atoms[k];
      a2 = mol2->getAtom(k);
      // sanity check - making sure Mol2 matches the selected ligand
      if ((a1->getResName() != a2->getResName()) || (a1->getAtmName() != a2->getAtmName())){
        cout << "Mol2 atoms do not match with selected molecule!" << endl;
        if ((a1->getResName() != a2->getResName())){
          cout << "Residue name mismatch at atom " << k << ": " << a1->getResName() << " and " << a2->getResName() << endl;
        }
        if ((a1->getAtmName() != a2->getAtmName())){
          cout << "Atom name mismatch at atom " << k << ": " << a1->getAtmName() << " and " << a2->getAtmName() << endl;
        }
        return 1;
      }
      atoms[k]->setAtmType(a2->getAtmType());
    }

    if (mol2->getNAtom() != atoms.size() ){
      cout << "mol2 size " << mol2->getNAtom() << " does not match mol size " << atoms.size() << endl;
      return 1;
    }
    return 0;
  }

  bool convert2PoseFP(PoseMolFP& poseFP, MolFP molFP){
    string atmType;
    for (MolFP::iterator iter = molFP.begin(); iter != molFP.end(); ++iter){
      Atom* atom = iter->first;
      AFP::AtomFP atomFP = iter->second;
      atmType = atom->getAtmType();
      if (poseFP.find(atmType) == poseFP.end() ) {
      // not found
        poseFP[atmType] = atomFP;
      }
      else{
      // found
        for(AFP::AtomFP::iterator it = atomFP.begin(); it != atomFP.end(); ++it){
          map<int, double > etaFP = it->second;
          for(map<int, double >::iterator it2 = etaFP.begin(); it2 != etaFP.end(); ++it2){
            poseFP[atmType][it->first][it2->first] += it2->second;
          }
        }
      }
    }
    return 0;
  }

  bool printPoseFP(PoseMolFP poseFP){
    for (PoseMolFP::iterator iter = poseFP.begin(); iter != poseFP.end(); ++iter){
      string atmType = iter->first;
      cout << atmType << ":" << endl;
      AFP::printFP(iter->second);
    }
    return 0;
  }
}


namespace TFP{
  bool getPoseFP(PoseFP& poseFP, vector<Atom*> atoms, Trajectory* ftrjin, ifstream trjin,
             double cutoff, vector<int> etas, vector<string> neighTypes,
             int start=0, int stop=-1, int skip=0, bool startFlag=false){
    MFP::MolFP molFP;
    MFP::PoseMolFP poseMolFP;
    unsigned int i;
    int frameid;
    /* Print out current trajectory info */
    cout << "Number of frames in traj: " << ftrjin->getNFrame() << endl;
    /* Loop through desired frames */
    if (skip > 0 && startFlag == false) start=skip;
    for (i = start; i < ftrjin->getNFrame() && i < stop; i = i + 1 + skip){
      if (ftrjin->readFrame(trjin, i) == false){ // read frame and setmol here
        std::cerr << "Warning: EOF found before the next frame could be read" << std::endl;
        break;
      }
      frameid = i + 1;
      cout << endl << "Frame " << frameid << ":" << endl;
      Molecule* mol = ftrjin->getMolecule();
      /* Print out molecule info at the first frame*/
      if (i == start){
        cout << "Number of Atoms: " << mol->getAtmVecSize() << endl;
        cout << "Number of Residues: " << mol->getResVecSize() << endl;
        cout << "Number of Chains: " << mol->getChnVecSize() << endl;
      }
      MFP::getFP(molFP, atoms, mol, cutoff, etas, neighTypes);
      MFP::convert2PoseFP(poseMolFP, molFP);
      poseFP[frameid] = poseMolFP;
    }
    return 0;
  }

  bool printPoseFP(PoseFP poseFP){
    for (PoseFP::iterator iter = poseFP.begin(); iter != poseFP.end(); ++iter){
      int frameid = iter->first;
      cout << frameid << ":" << endl;
      MFP::printPoseFP(iter->second);
    }
    return 0;
  }

  //
  bool normalizeFP(PoseFP& poseFP){
    // get keys
    vector<int> keys = getKeys(poseFP);
    for (unsigned int i = 0; i < keys.size(); i ++){
      cout << keys[i] << " ";
    }
    cout << endl;

  //   MFP::PoseMolFP poseMolFP = psoeFP[keys[0]];
  //     for (MFP::PoseMolFP::iterator poseMol = poseMolFP.begin(); poseMol != poseMolFP.end(); ++poseMol){
  //       for(AFP::AtomFP::iterator atom = poseMol->second.begin(); atom != poseMol->second.end(); ++atom){
  //         for(map<int, double >::iterator eta = atom->second.begin(); eta != atom->second.end(); ++eta){
  //           for k in keys:
  //             vector.push_back(poseFP[i][j][k][l]);
  //           median[i][j][k] = median(vector);
  //           vect0r.clear()
  //           eta->second;
  //           getMedian(atom->second);
  //       }
  //     }
  //   }
    return 0;
  }

}



namespace FileProcess{
  Molecule* loadMol2(string mol2file){
    Molecule *mol2 = Molecule::readMol2(mol2file);
    return mol2;
  }

  Trajectory* loadTraj(ifstream& trjin, string trajFile, Molecule* mol){
    trjin.open(trajFile.c_str(), std::ios::binary);
    if (trjin.is_open()){
      Trajectory* ftrjin=new Trajectory;
      ftrjin->setMolecule(mol);
      if (ftrjin->findFormat(trjin) == true){
        ftrjin->readHeader(trjin);
      }
      return ftrjin;
    }
    else{
      std::cerr << "Warning: Skipping unknown trajectory format \"";
      std::cerr << trajFile << "\"" << std::endl;
      return NULL;
    }
  }
}

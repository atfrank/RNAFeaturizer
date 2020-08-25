#include "Grid.hpp"

// hash function for array to be used with unordered_set
// ref: https://stackoverflow.com/questions/8026890/c-how-to-insert-array-into-hash-set/8026914#8026914
namespace std
{
    template<>
    struct hash<GRID::Index>
    {
        typedef GRID::Index argument_type;
        typedef size_t result_type;

        result_type operator()(const argument_type& a) const
        {
            hash<int> hasher;
            result_type h = 0;
            h = h * 31 + hasher(a.x);
            h = h * 31 + hasher(a.y);
            h = h * 31 + hasher(a.z);
            return h;
        }
    };
}

namespace GRID{
  // get min max boundary in each direction
  // boundary is a vector of 6 elements arranged as [xmin, xmax, ymin, ymax, zmin, zmax]
  bool getMinMax(Bound& boundary, Molecule *mol, double radius){
    for (unsigned int i = 0; i < mol->getAtmVecSize(); i++){
        double xcoor = mol->getAtom(i)->getX();
        double ycoor = mol->getAtom(i)->getY();
        double zcoor = mol->getAtom(i)->getZ();

        if (xcoor > boundary.xmax){
            boundary.xmax = xcoor;
        }
        if (ycoor > boundary.ymax){
            boundary.ymax = ycoor;
        }
        if (zcoor > boundary.zmax){
            boundary.zmax = zcoor;
        }
        if (xcoor < boundary.xmin){
            boundary.xmin = xcoor;
        }
        if (ycoor < boundary.ymin){
            boundary.ymin = ycoor;
        }
        if (zcoor < boundary.zmin){
            boundary.zmin = zcoor;
        }
    }
    boundary += radius;
    boundary.print();
    return 0;
  }

  // coor to index in 1D
  int c2i(double coor, double min, double spacing){
    return floor((coor - min) / spacing);
  }

  // index to coor in 1D
  double i2c(int index, double min, double spacing){
    return min + index * spacing;
  }

  Index coor2Index(const Coor& coor, const Bound& boundary, double spacing){
    double xcoor = coor.x();
    double ycoor = coor.y();
    double zcoor = coor.z();

    Index index;

    index.x = floor((xcoor - boundary.xmin) / spacing);
    index.y = floor((ycoor - boundary.ymin) / spacing);
    index.z = floor((zcoor - boundary.zmin) / spacing);

    return index;
  }


  Coor index2Coor(const Index& index, const Bound& boundary, double spacing){

    double xcoor = boundary.xmin + index.x * spacing;
    double ycoor = boundary.ymin + index.y * spacing;
    double zcoor = boundary.zmin + index.z * spacing;

    Coor coor(xcoor, ycoor, zcoor);
    return coor;
  }

  bool getBoundary(Bound& boundAtom, const Coor& center, double radius){
    double xcoor = center.x();
    double ycoor = center.y();
    double zcoor = center.z();

    boundAtom.xmin = xcoor - radius;
    boundAtom.xmax = xcoor + radius;
    boundAtom.ymin = ycoor - radius;
    boundAtom.ymax = ycoor + radius;
    boundAtom.zmin = zcoor - radius;
    boundAtom.zmax = zcoor + radius;

    return 0;
  }


  void getBoundIndex(Index& upperLeft, Index& lowerRight, const Bound& boundary,
                     const Coor& center, double radius, double spacing){
    Bound boundAtom;
    getBoundary(boundAtom, center, radius);
    // convert to index
    upperLeft.x = c2i(boundAtom.xmin, boundary.xmin, spacing);
    lowerRight.x = c2i(boundAtom.xmax, boundary.xmin, spacing);
    upperLeft.y = c2i(boundAtom.ymin, boundary.ymin, spacing);
    lowerRight.y = c2i(boundAtom.ymax, boundary.ymin, spacing);
    upperLeft.z = c2i(boundAtom.zmin, boundary.zmin, spacing);
    lowerRight.z = c2i(boundAtom.zmax, boundary.zmin, spacing);
  }


  bool getGrids(unordered_set<Index>& gridIndex, unordered_set<Index>& exclude,
                const vector<Atom*>& atoms, const Bound& boundary,
                double upLimit, double lowLimit, double spacing){
    // idea: First, loop over all RNA atoms to get the list of excluded grids
    // Next, loop over all RNA atoms again to get the list of grids
    // Finally, label the grids based on whether they are in the same cell as Mg2+
    Index ul, lr; // indices of the upperleft (minimum) and lowerright (maximum) corner of an atom
    // A grid is excluded if the distance < lowlimit between the grid center and an RNA atom
    // See Figure 1 of Yu, Zeyun, A List-Based Method for Fast Generation of Molecular Surfaces, 2009.
    for (unsigned int i = 0; i < atoms.size(); i++){
      Coor coor =  atoms[i]->getCoor();
      getBoundIndex(ul, lr, boundary, coor, lowLimit, spacing);
      for (int ix = ul.x; ix <= lr.x; ix++){
        for (int iy = ul.y; iy <= lr.y; iy++){
          for (int iz = ul.z; iz <= lr.z; iz++){
            Coor grid(i2c(ix, boundary.xmin, spacing), i2c(iy, boundary.ymin, spacing), i2c(iz, boundary.zmin, spacing));
            double dist = Analyze::distance(coor, grid);
            if (dist <= lowLimit){
              Index indices(ix, iy, iz);
              // Element is inserted only if it is not equivalent to any other element already in the container
              exclude.insert(indices);
            }
          }
        }
      }
    }
    // cout << "Grids excluded: " << endl;
    // for (auto const& ind : exclude){
    //  ind.print();
    // }

    // Get included indices
    for (unsigned int i = 0; i < atoms.size(); i++){
      Coor coor =  atoms[i]->getCoor();
      getBoundIndex(ul, lr, boundary, coor, upLimit, spacing);
      for (int ix = ul.x; ix <= lr.x; ix++){
        for (int iy = ul.y; iy <= lr.y; iy++){
          for (int iz = ul.z; iz <= lr.z; iz++){
            Coor grid(i2c(ix, boundary.xmin, spacing), i2c(iy, boundary.ymin, spacing), i2c(iz, boundary.zmin, spacing));
            double dist = Analyze::distance(coor, grid);
            Index indices(ix, iy, iz);
            if (dist <= upLimit && exclude.find(indices) == exclude.end()){
              gridIndex.insert(indices);
            }
          }
        }
      }
    }
    return 0;
  }

  bool getGridsMol(vector <Coor>& grids, Molecule* mol, double cutoff, double& spacing,
                const unordered_map<string, array<double, 2> >& dConstraints){
    Bound boundary;
    getMinMax(boundary, mol, cutoff);
    unordered_set<Index> gridIndex, exclude;
    for (auto const& lim: dConstraints){
      // select atoms
      vector<Atom*> atoms = Select::makeSelVec(mol, lim.first);
      // get grids based on distance constraints
      array<double, 2> limits  = lim.second;
      double lowLimit = limits[0];
      double upLimit = limits[1];
      getGrids(gridIndex, exclude, atoms, boundary, upLimit, lowLimit, spacing);
      // grids.insert(grids.end(), selGrids.begin(), selGrids.end());
      cout << "Number of grids placed so far: " << gridIndex.size() << endl;
    }

    // cout << "Grids Included: " << endl;
    for (auto const& ind : gridIndex){
      // ind.print();
      Coor coor = index2Coor(ind, boundary, spacing);
      grids.push_back(coor);
      cout << "(" << coor.x() << "," << coor.y() << "," << coor.z() << ") ";
    }
    cout << endl << "Total number of grids: " << gridIndex.size() << endl;
    return 0;
  }


}

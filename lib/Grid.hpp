// WORK IN PROGRESS
// Place grids on the surface of an RNA molecule, s.t. distance constraints
// Has been used for: Mg2+ binding site prediction
// Each grid location is represented by a Coor that contains its 3D coordinates
// Or an "Index" [n, m, l], a set of 3 integers that indexes the grid in x,y and z direction
// Coor can be computed from Index using limites and spacing
// and vice versa

#ifndef GRID_hpp
#define GRID_hpp

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <array>
#include <unordered_set>
#include <unordered_map>

#include "Coor.hpp"
#include "Atom.hpp"
#include "Analyze.hpp"
#include "Molecule.hpp"
#include "Mol2.hpp"
#include "Trajectory.hpp"
#include "Select.hpp"

#define PI 3.14159265

using namespace std;

namespace GRID{
  // Index: counterpart of Coor, a 3-array of grid indices
  struct Index {
    unsigned int x, y, z;
    Index() : x(0), y(0), z(0) { }
    Index(int ix, int iy, int iz) : x(ix), y(iy), z(iz) { }
    bool operator==(const Index& other) const{
      return x == other.x && y == other.y && z == other.z;
    }
    void print() const{
      cout << "Indices: [" << x << ", " << y << ", " << z << "]" << endl;
    }
  };

  // Boundary: minimum and maximum in each direction, a 6-array of coordinates
  struct Bound{
    double xmax, xmin, ymax, ymin, zmax, zmin;
    Bound() : xmax(-9999), xmin(9999), ymax(-9999), ymin(9999), zmax(-9999), zmin(9999)  { }
    void operator+=(const double& expand){
      xmax += expand;
      ymax += expand;
      zmax += expand;
      xmin -= expand;
      ymin -= expand;
      zmin -= expand;
    }
    void print() const{
      cout << "x: [" << xmin << ", " << xmax << "]" << endl;
      cout << "y: [" << ymin << ", " << ymax << "]" << endl;
      cout << "z: [" << zmin << ", " << zmax << "]" << endl;
    }
  };

  // get min max boundary in each direction
  // boundary is a vector of 6 elements arranged as [xmin, xmax, ymin, ymax, zmin, zmax]
  bool getMinMax(Bound& boundary, Molecule *mol, double radius);
  // coordinate to index and index to coordinate in 1D
  int c2i(double coor, double min, double spacing);
  double i2c(int index, double min, double spacing);
  // coordinate to index and index to coordinate in 3D
  Index coor2Index(const Coor& coor, const Bound& boundary, double spacing);
  Coor index2Coor(const Index& index, const Bound& boundary, double spacing);
  // get boundary of indices within radius around center
  bool getBoundary(Bound& boundary, const Coor& center, double radius);
  void getBoundIndex(Index& upperLeft, Index& lowerRight, const Bound& boundary,
                     const Coor& center, double radius, double spacing);
  // main function for getting grids
  bool getGrids(unordered_set<Index>& gridIndex, unordered_set<Index>& exclude, const vector<Atom*>& atoms, const Bound& boundary,
                double upLimit, double lowLimit, double spacing);
  bool getGridsMol(vector <Coor>& grids, Molecule* mol, double cutoff, double& spacing,
                   const unordered_map<string, array<double, 2> >& dConstraints);
}

#endif

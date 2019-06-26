/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_algorithms_HashGrid_cpp
#define bmpg_uncc_edu_algorithms_HashGrid_cpp

#include <iostream>
#include <vector>
#include <cmath>
#include <bmpg_uncc_edu/chemistry/PDBAtom.hpp>
#include <bmpg_uncc_edu/chemistry/Molecule.hpp>
#include <bmpg_uncc_edu/algorithms/HashGrid.hpp>

namespace bmpg_uncc_edu {
namespace algorithms {
using namespace bmpg_uncc_edu::chemistry;
using namespace bmpg_uncc_edu::util;

/**
	Constructor.
	@param mol molecule
	@param length grid length
*/
HashGrid::HashGrid(Molecule & mol,
		   double length) :
		   grid_length_(length)
{
	setup(mol.atoms());
	add_atoms(mol.atoms());
}

/**
	Constructor.
	@param aotms list of atoms
	@param length grid length
*/
HashGrid::HashGrid(const atom_list_t & atoms,
		   double length) :
		   grid_length_(length)
{
	setup(atoms);
	add_atoms(atoms);
}

/**
	Index of the grid specified by its x, y, and z positions.
	@param x coordinate
	@param y coordinate
	@param z coordinate
	@return index number
*/
int HashGrid::index(int x,
		    int y,
		    int z) const
{
	return z*nx_*ny_ + y*nx_ + x;
}

/**
	Convert a coordinate (a vector of 3 integers) to index.
	@param c coordinate
	@return index number
*/
int HashGrid::index(const coord_t & c) const
{
	return index(c.at(0),c.at(1),c.at(2));
}

/**
	Get the list of atoms in the neighboring grids and the current grid of a given atom.
	@param atom Atom
	@return list of atoms
*/
HashGrid::atom_list_t HashGrid::atom_neighbors(const PDBAtom* atom){
	atom_list_t atoms;
	grid_list_t grids = grid_neighbors(atom);
	
	grid_list_t::iterator ig;
	for(ig = grids.begin(); ig != grids.end(); ig++){
		atom_list_t& a = map_[*ig];					//returns empty list if a grid does not exist
		atoms.insert(atoms.begin(),a.begin(),a.end());			//append to the overall neighbor list
	}
	return atoms;
}

HashGrid::atom_list_t HashGrid::atom_neighbors(float x, float y, float z){
        coord_t co = coord(x, y, z);
	atom_list_t atoms;
	grid_list_t grids = grid_neighbors(co);
	
	grid_list_t::iterator ig;
	for(ig = grids.begin(); ig != grids.end(); ig++){
		atom_list_t& a = map_[*ig];					//returns empty list if a grid does not exist
		atoms.insert(atoms.begin(),a.begin(),a.end());			//append to the overall neighbor list
	}
	return atoms;
}

HashGrid::atom_list_t HashGrid::atom_neighbors_distance(float x, float y, float z, float probe, float boundary, int &returnCode){
        double lx = x - xmin_;
	double ly = y - ymin_;
	double lz = z - zmin_;
	int xi = (int)(lx/grid_length_) + 1;
	int yi = (int)(ly/grid_length_) + 1;
	int zi = (int)(lz/grid_length_) + 1;
    
	bool isOccupied = false;
        bool isCavity  =  false;
        bool isBoundary = false;
        bool isNeighbor = false;
        HashGrid::atom_list_t closeNeighbors;
        float distance = 0;
        
        int hash;
        
        for(int i = -1; i <= 1; i++){
            for(int j = -1; j <= 1; j++){
		for(int k = -1; k <= 1; k++){
                    hash = index(xi + i, yi + j, zi + k);
                    atom_list_t& a = map_[hash];					  
    
                    float vdw;
                    for (HashGrid::atom_list_t::const_iterator it = a.begin(); it != a.end(); it++){
                        isNeighbor = true;
                        const PDBAtom* atom = *it;  
                        vdw = atom->occupancy;
        
                        distance = (x - atom->x) * (x - atom->x) + (y - atom->y) * (y - atom->y) + (z - atom->z) * (z - atom->z);   
                        if (distance <= vdw*vdw) {
                            isOccupied = true;
                            break;
                        } 
                        else if (distance <= ((2 * probe + vdw) * (2 * probe + vdw) )) {
                            closeNeighbors.push_back(*it);
                            isCavity = true;
                        }
                        else if (distance < (boundary * boundary)) {
                            isBoundary = true;
                        }
                    }                    
		}
            }
	}
        if (isOccupied) {
            returnCode = 0;                                                     // Protein volume
        } else if (isCavity) {      
            returnCode = -1;                                                    // Void or microvoid
        } else if (isBoundary) {
            returnCode = 2;                                                     // Boundary Solvent
        } else if (!isNeighbor) {
            returnCode = 4;                                                     
        } else {
            returnCode = 3;   
        }
        return closeNeighbors;
}

/**
	Get the list of neighboring grids and the current grid of a given atom.
	@param atom Atom
	@return list of grids
*/
HashGrid::grid_list_t HashGrid::grid_neighbors(coord_t co)
{
	grid_list_t list;
	int x = co.at(0);
	int y = co.at(1);
	int z = co.at(2);
	for(int i = -1; i <= 1; i++){
		for(int j = -1; j <= 1; j++){
			for(int k = -1; k <= 1; k++){
				list.push_back(index(x + i, y + j, z + k));
			}
		}
	}
	return list;
}

HashGrid::grid_list_t HashGrid::grid_neighbors(const PDBAtom* atom)
{
	grid_list_t list;
	coord_t co = coord(atom);
	int x = co.at(0);
	int y = co.at(1);
	int z = co.at(2);
	for(int i = -1; i <= 1; i++){
		for(int j = -1; j <= 1; j++){
			for(int k = -1; k <= 1; k++){
				list.push_back(index(x + i, y + j, z + k));
			}
		}
	}
	return list;
}

/**
	Add a list of atoms to the hash grid.
	@param atoms list of atoms
*/
void HashGrid::add_atoms(const atom_list_t & atoms)
{
	atom_list_t::const_iterator it;
	for(it = atoms.begin(); it != atoms.end(); it++){
		int i = index(coord(*it));					//get the global index	
		map_[i].push_back(*it);		
	}
}

/**
	Get the coordinate of a given atom.
	@param atom Atom
	@return coordinate (vector of 3 integers)
*/
HashGrid::coord_t HashGrid::coord(const PDBAtom* atom)
{
	double lx = atom->x - xmin_;
	double ly = atom->y - ymin_;
	double lz = atom->z - zmin_;
	int i,j,k;
	i = (int)(lx/grid_length_) + 1;
	j = (int)(ly/grid_length_) + 1;
	k = (int)(lz/grid_length_) + 1;
	coord_t c;
	c.push_back(i);
	c.push_back(j);
	c.push_back(k);
	return c;
}

HashGrid::coord_t HashGrid::coord(float x, float y, float z)
{
	double lx = x - xmin_;
	double ly = y - ymin_;
	double lz = z - zmin_;
	int i,j,k;
	i = (int)(lx/grid_length_) + 1;
	j = (int)(ly/grid_length_) + 1;
	k = (int)(lz/grid_length_) + 1;
	coord_t c;
	c.push_back(i);
	c.push_back(j);
	c.push_back(k);
	return c;
}


/**
	Find the x, y, and z ranges of a given list of atoms.
	@param atoms list of atoms
	@param xmin min x
	@param xmax max x
	@param ymin min y
	@param ymax max y
	@param zmin min z
	@param zmax max z
*/
void HashGrid::find_ranges(const atom_list_t & atoms,
			   double & xmin,
			   double & xmax,
			   double & ymin,
			   double & ymax,
			   double & zmin,
			   double & zmax)
{
	double limit = numeric_limits<double>::max();
	xmax = ymax = zmax = -limit;
	xmin = ymin = zmin = limit;
	atom_list_t::const_iterator it;
	for(it = atoms.begin(); it != atoms.end(); it++){
		const PDBAtom* atom = *it;
		if(xmin > atom->x)xmin = atom->x;
		if(ymin > atom->y)ymin = atom->y;
		if(zmin > atom->z)zmin = atom->z;
		
		if(xmax < atom->x)xmax = atom->x;
		if(ymax < atom->y)ymax = atom->y;
		if(zmax < atom->z)zmax = atom->z;
	}
	
}
void HashGrid::setup(const atom_list_t & atoms)
{
	double xmin,xmax,ymin,ymax,zmin,zmax;
	find_ranges(atoms, xmin,xmax,ymin,ymax,zmin,zmax);
	
	xmin_ = xmin;
	ymin_ = ymin;
	zmin_ = zmin;
	
	double lx = xmax - xmin;
	double ly = ymax - ymin;
	double lz = zmax - zmin;
	nx_ = (int)(lx/grid_length_) + 1;
	ny_ = (int)(ly/grid_length_) + 1;
	nz_ = (int)(lz/grid_length_) + 1;
}

ostream& operator<<(ostream & out, HashGrid & rhs)
{
	out << "nx = " << rhs.nx_ << " " << "nx = " << rhs.ny_ << " " << "nx = " << rhs.nz_ << endl;
	for(int i = 0; i < rhs.nx_; i++){
		for(int j = 0; j < rhs.ny_; j++){
			for(int k = 0; k < rhs.nz_; k++){
				int ind = rhs.index(i,j,k);
				HashGrid::atom_list_t a = rhs.map_[ind];
				HashGrid::atom_list_t::const_iterator it;
				if(a.size() == 0)continue;
				out << "grid ("<<i << "," << j << "," << k << ")" << endl;
				for(it = a.begin(); it != a.end(); it++){
					string s = (*it)->atom_name;
					 out << s << " ";					 
				}
				out << "\n";
			}
		}
	}
	
	return out;
	
}
	


}	//namespace bmpg_uncc_edu::algorithms
}	//namespace bmpg_uncc_edu

#endif


/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_algorithms_HashGrid_hpp
#define bmpg_uncc_edu_algorithms_HashGrid_hpp

#include <limits>
#include <vector>
#include <map>

namespace bmpg_uncc_edu {
	namespace chemistry {
		class PDBAtom;
		class Molecule;
	}
}

namespace bmpg_uncc_edu {
namespace algorithms {
using namespace bmpg_uncc_edu::chemistry;
	
/**
	The HashGrid class sets up a hash grid for a given list of atoms. 
	If the grid size is chosen properly, to retrieve the neighbors
	of a given atom, only the atoms in the neighboring grids needs to be 
	checked for distance cutoff. Hence the HashGrid class
	provides a fast way for setting up neighbor lists.
*/
class HashGrid
{
public:
	typedef std::vector<PDBAtom*> atom_list_t;
	typedef std::vector<int> coord_t;					//x,y,z coordinates of a grid
	typedef std::vector<int> grid_list_t;					//list of global indices
	
	typedef map<int,atom_list_t> map_t;					//maps a grid index to the associated atoms

	HashGrid(const atom_list_t & list,
		 double length);
	
	HashGrid(Molecule & mol,
		 double length);

	size_t size() const;
	grid_list_t grid_neighbors(const PDBAtom* atom);
	atom_list_t atom_neighbors(const PDBAtom* atom);        
	grid_list_t grid_neighbors(coord_t co);
	atom_list_t atom_neighbors(float x, float y, float z);
	atom_list_t atom_neighbors_distance(float x, float y, float z, float probe, float boundary, int &returnCode);
	coord_t coord(const PDBAtom* atom);
        coord_t coord(float x, float y, float z);
	
	int index(int x,
		  int y,
		  int z) const;
	
	int index(const coord_t & c) const;
	friend ostream& operator<<(ostream & out, HashGrid & rhs);
	
private:
	void setup(const atom_list_t & atoms);
	void add_atoms(const atom_list_t & atoms);
	void find_ranges(const atom_list_t & atoms,
			 double & xmin,
			 double & xmax,
			 double & ymin,
			 double & ymax,
			 double & zmin,
			 double & zmax);
	
private:
	double grid_length_;
	atom_list_t atom_list_;							//list of all atoms in the grid
	map_t map_;
	int nx_, ny_, nz_;							//number of cells in the three directions
	double xmin_, ymin_, zmin_;
};

}	//namespace bmpg_uncc_edu::algorithms
}	//namespace bmpg_uncc_edu

#endif


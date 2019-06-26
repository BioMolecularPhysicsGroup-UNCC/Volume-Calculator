/* 
 * File:   HashGridVolume.hpp
 * Author: jenny
 *
 * Created on November 13, 2016, 7:23 AM
 */

#ifndef HASHGRIDVOLUME_HPP
#define	HASHGRIDVOLUME_HPP


#include <float.h>
#include <fstream>
#include <cmath>
#include <iostream>
#include "bmpg_uncc_edu/chemistry/PDBAtom.hpp"
#include "bmpg_uncc_edu/chemistry/PDBProtein.hpp"
#include "SpringModel.hpp"

using namespace std;
using namespace bmpg_uncc_edu::chemistry;

class HashGridVolume {

public:
    HashGridVolume(bool flagOpt, float grid, float probe, float shell, float cutoff, float ** atomCoords, int nAtoms);
    virtual ~HashGridVolume();
    
    void  CompleteSliceIn27Sections();
    void  volTypeSection(int ix, int jy, int kz, int c1, int c2, int * di, int * dj, int * dk);
    void  prepareInitialSlice();
    bool  prepareNextSlice();
    int   getVolumeType(int i, int j);
    float getHashSpacing();
    void  writeShiftedCoordinates(string PDBFileName);  
    
    int     nAtoms;
    int     closeAtom;
    int     veryCloseAtom;

    static const int maxMaxAtomsWithin = 1000;                                  // 1000 potential atom-probe clashes can be handled
    float   probe;
    float   grid;
    float   shell;
    float   cutoff;
    int     ni;
    int     nj;
    int     nk;
    
    float   xShift;
    float   yShift;
    float   zShift;
    
    float   xRange;
    float   yRange;
    float   zRange;
    
private: 
      
    bool    flagOpt_; 
    int     skipSpringMethod_;    
    float   maxRvdw_;                                                           //maximum van der Waal radius
    float   dBuffer_;                                                           //shifts coordinates to place atom hash grid in first quadrant in x, y, z directions 
    float   xmin_;
    float   xmax_;
    float   ymin_;
    float   ymax_;
    float   zmin_;
    float   zmax_;
// ------------------------------- setup all variables required for hash grid to efficiently locate atoms
    int     boxLx_; 
    int     boxLy_; 
    int     boxLz_; 
    int     n1st_;               //n1st_ = 26 => (0 -> 26) defines 27 grid cells used to check potential clashes with atoms
    int     nnxt_;               //nnxt_ = n1st_ + 1 = 27 = 1st index for 2nd layer of space to check if within boundary volume
    int     n2nd_;               //n2nd_ = 124 = last index for 2nd layer of space:  Note 0 to 124 gives 125 different cells
    int     * di_;               //displacement index for x-direction on discrete grid
    int     * dj_;               //displacement index for y-direction on discrete grid
    int     * dk_;               //displacement index for z-direction on discrete grid
    int     *gridChain_;
    int     *** hashGrid_; 

// ------------------------------------------------------ restricted displacement vector variables
    int    nnxtLLL_, n2ndLLL_, *diLLL_, *djLLL_, *dkLLL_; 
    int    nnxtLML_, n2ndLML_, *diLML_, *djLML_, *dkLML_; 
    int    nnxtLUL_, n2ndLUL_, *diLUL_, *djLUL_, *dkLUL_; 
    int    nnxtMLL_, n2ndMLL_, *diMLL_, *djMLL_, *dkMLL_; 
    int    nnxtMML_, n2ndMML_, *diMML_, *djMML_, *dkMML_; 
    int    nnxtMUL_, n2ndMUL_, *diMUL_, *djMUL_, *dkMUL_; 
    int    nnxtULL_, n2ndULL_, *diULL_, *djULL_, *dkULL_; 
    int    nnxtUML_, n2ndUML_, *diUML_, *djUML_, *dkUML_; 
    int    nnxtUUL_, n2ndUUL_, *diUUL_, *djUUL_, *dkUUL_; 

    int    nnxtLLM_, n2ndLLM_, *diLLM_, *djLLM_, *dkLLM_; 
    int    nnxtLMM_, n2ndLMM_, *diLMM_, *djLMM_, *dkLMM_; 
    int    nnxtLUM_, n2ndLUM_, *diLUM_, *djLUM_, *dkLUM_; 
    int    nnxtMLM_, n2ndMLM_, *diMLM_, *djMLM_, *dkMLM_; 
  //int    nnxtMMM_, n2ndMMM_, *diMMM_, *djMMM_, *dkMMM_;  <----- not needed
    int    nnxtMUM_, n2ndMUM_, *diMUM_, *djMUM_, *dkMUM_; 
    int    nnxtULM_, n2ndULM_, *diULM_, *djULM_, *dkULM_; 
    int    nnxtUMM_, n2ndUMM_, *diUMM_, *djUMM_, *dkUMM_; 
    int    nnxtUUM_, n2ndUUM_, *diUUM_, *djUUM_, *dkUUM_; 

    int    nnxtLLU_, n2ndLLU_, *diLLU_, *djLLU_, *dkLLU_; 
    int    nnxtLMU_, n2ndLMU_, *diLMU_, *djLMU_, *dkLMU_; 
    int    nnxtLUU_, n2ndLUU_, *diLUU_, *djLUU_, *dkLUU_; 
    int    nnxtMLU_, n2ndMLU_, *diMLU_, *djMLU_, *dkMLU_; 
    int    nnxtMMU_, n2ndMMU_, *diMMU_, *djMMU_, *dkMMU_; 
    int    nnxtMUU_, n2ndMUU_, *diMUU_, *djMUU_, *dkMUU_; 
    int    nnxtULU_, n2ndULU_, *diULU_, *djULU_, *dkULU_; 
    int    nnxtUMU_, n2ndUMU_, *diUMU_, *djUMU_, *dkUMU_; 
    int    nnxtUUU_, n2ndUUU_, *diUUU_, *djUUU_, *dkUUU_; 


// --------------------------- define range of fine grid in x,y,z directions for optimized algorithm
    int iL0_,iL1_,iM0_,iM1_,iU0_,iU1_;
    int jL0_,jL1_,jM0_,jM1_,jU0_,jU1_;
    int kL0_,kL1_,kM0_,kM1_,kU0_,kU1_;

// --------------------------- setup all variables required for layer to efficiently store type of grid points
    int     *** gridLayer_;  //stores a snap shot of volume type for all gridpoints in all slices comprising the layer 
    int     *** atomLayer_;  //stores closest atom labels for all gridpoints that correspond to gridLayer_
    int     *** closeAtomLayer_;  //stores closest atom labels for all gridpoints that correspond to gridLayer_
    bool    *** freeVolume;  //stores closest atom labels for all gridpoints that correspond to gridLayer_
    int     tag_;            //shift in identification labels that are active
    int     initTag_;         //initial tag
    static const int dTag_ = 4;   //includes variations of {-1,0,1,2} for {microvoid, protein-volume, void-not-bulk, void-bulk}
    int     minTag_;         //indicates old information in gridLayer whenever (tag_ - minTag_) < 0;
    int     refTag_;         //= tag_ + 1;
    int     kLayer_;         //width of layer
    int     * kLayerFromK_;  //array that stores precalculation of  k % kLayer_;
    int     * vTypeFromTag_; //array that stores precalculation of  tag_ % dTag_; 
    int     k_;              //actual interator in the z-direction, most forward position
    int     kLpos_;          //slice index in layer that corresponds to the forward  position  (pos => +)
    int     kLmid_;          //slice index in layer that corresponds to the current  position  (mid => 0)
    int     kLneg_;          //slice index in layer that corresponds to the previous position  (neg => -)
       
    float   maxExtension_;
    float   hashSpacing_;
    float   hashSpacing2_;  //squared hashSpacing_
    float   boundarySquared_;
    float   probeDiameter_;
    float   probeSquared_;
    float   ** atomCoords_;
    
    int     * atomsWithin_; 
    SpringModel * spring_;
    
    void generateHashGrid();
    void buildGridSearchPatterns();
    int  gridCount(int iStart, int iStop, int jStart, int jStop, int kStart, int kStop, int count);
    int  gridCount2(int iStart, int iStop, int jStart, int jStop, int kStart, int kStop, int count);
    void deleteOptimizedSearchPatterns();
    void generateGridLayer();
    
   
};

#endif	/* HASHGRIDVOLUME_HPP */


/* 
 * File:   HoshenKopelman.hpp
 * Author: jenny
 *
 * Created on November 18, 2016, 5:03 PM
 */

#ifndef HOSHENKOPELMAN_HPP
#define	HOSHENKOPELMAN_HPP

#include <assert.h>
#include <string>   
#include "HashGridVolume.hpp"
#include <fstream>
#include "bmpg_uncc_edu/chemistry/PDBAtom.hpp"

//#define VISUALIZE

using namespace std;
using namespace bmpg_uncc_edu::chemistry;

class HoshenKopelman {
public:
    HoshenKopelman(HashGridVolume * hashGrid, string fileBase);
    virtual ~HoshenKopelman();
    
    void process3DGrid();
    
    bool percolated = false;

    int * clusterSize;
    int   maxCluster;
    int   largestVoid;
    int   largestMicrovoid;
    
    int proteinVolume;
    int cavityVolume;
    int microvoidVolume;
    int boundaryVolume;
    int maximumClusterLabel;
    int numberOfMicrovoids;
    int numberOfVoids;
    int numOfErrorGdpts;
    
    int * microvoidPerAtom;
    int * boundaryPerAtom;
    int * cavityPerAtom;
    int * proteinPerAtom;
    
    int * cavityPerAtomClose;
    int * mvPerAtomClose;
    int * boundaryPerAtomClose;
    
    float rsm;
    
private:
    int ni_, nj_, nk_;
    int * clusterLabel_; 
    int   nLabels_; 
    int   minVoidSize_;
    int   voidBulkPerimeter_;
    
    int * cavityChain_;
    int * cavityLabel_;
    int * cavityCount_;
    int   cavityIndex_;
    
    int * cavityChainClose_;
    int * cavityLabelClose_;
    int * cavityCountClose_;
    int   cavityIndexClose_;
        
    float probe_;
    float grid_;
    int   nAtoms_;
    string fileBase_;

    int   nnHalf_;
    int * neighbor_;
    
    HashGridVolume * hashGrid_;
    ofstream    cavityVisual_;
    ofstream    mvVisual_;
    ofstream    mvCoord_;
    int mvAtomNum_ = 0;
    
    int uf_find(int x);
    int uf_union(int x, int y);
    int uf_make_set(int type);
    
    void processSlice(int ** priorSlice, int ** slice, int k);
    inline void getNeighbors(int ** priorSlice, int ** slice, int i, int j);
    void clusterFraction(int label, int atom);
    void clusterFractionClose(int label, int atom);
    void cluster();  
    void addVoidToFile(float x, float y, float z, int voidCluster, bool cavity);
    void clusterVoidFile(int falsePositives[], int nFalsePs, int cavities[], int nCavities, float largestCluster, bool cavity);


};

#endif	/* HOSHENKOPELMAN_HPP */


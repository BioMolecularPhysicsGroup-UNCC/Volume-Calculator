/* 
 * File:   HashGridVolume.cpp
 * Author: jenny
 * 
 * Created on November 13, 2016, 7:23 AM
 */

#include "HashGridVolume.hpp"

HashGridVolume::HashGridVolume(bool flagOpt, float grid, float probe, float shell, float cutoff, float ** atomCoords, int nAtoms) 
                             : flagOpt_(flagOpt),
                               grid(grid),
                               probe(probe),
                               shell(shell),
                               cutoff(cutoff),
                               atomCoords_(atomCoords), 
                               nAtoms(nAtoms) {
    if (flagOpt == 0) {                                                         //apply spring method for low number of constraints
       skipSpringMethod_ = 1;                                                   //in all cases a potential clash with 1 atom can ALWAYS be resolved exactly 
    } else {                                                                    //use empirical based results and make an assumption that voids will ALWAYS 
                                                                                //  occur below cutoff
        if (probe < 1.0) {
          skipSpringMethod_ = 1;                           
        } else if (probe < 2.4) {
          skipSpringMethod_ = floor(1.01 + (probe-0.8)*5);
        } else if (probe < 3.0) {
          skipSpringMethod_ = floor(9.01 + (probe-2.4)*10);
        } else if (probe < 4.0) {
          skipSpringMethod_ = floor(16.01 + (probe-3.0)*15);
        } else {                                              
          skipSpringMethod_ = floor(32.01 + (probe-4.0)*20);                    //unexplored after 3.8 but the pattern from before is repeated conservatively
        }

        if (grid < 0.3999999) {                                                 //smaller grid size implies user wants MORE ACCURACY over speed.
            skipSpringMethod_ -= floor( 0.019 + (0.4 - grid)*20 ); 
            if (skipSpringMethod_ < 1) { 
                skipSpringMethod_ = 1; 
            }
        }
    }

    maxRvdw_ = -1.0;                                                            //determine maximum van der Waals radius in system
    for (int n=1; n<=nAtoms; n++ ) {          
        if( atomCoords_[n][0] > maxRvdw_ ){
            maxRvdw_ = atomCoords_[n][0];                        
        }
    }
    probeDiameter_ = 2 * probe;
    probeSquared_ = probe * probe;  
    boundarySquared_ = (maxRvdw_ + shell)*(maxRvdw_ + shell);
    
    maxExtension_ = 0.0001 + maxRvdw_ + shell;                                  //0.0001 prevents dilemma of an exact equality
    float temp2 = maxExtension_/2;                                              //must have:  2*hashSpacing_ >= maxExtension_ 
    float temp1 = 0.0001 + maxRvdw_ + probeDiameter_;                           //this is the minimum possible hashSpacing_
    if (temp2 > temp1) {
        hashSpacing_ = temp2;
    } else {
        hashSpacing_ = temp1;
    }
    maxExtension_ = 2*hashSpacing_;                                             //Could be larger than (0.0001 + maxRvdw_ + shell)
    hashSpacing2_ = hashSpacing_*hashSpacing_;                                  //Could be larger than (0.0001 + maxRvdw_ + shell)
    
    
    generateHashGrid();
    buildGridSearchPatterns();
                                                                                //develop displacement indices to span level 1 and 2 neighborhoods
    generateGridLayer();    
   
    atomsWithin_      = (int *)   new   int[maxMaxAtomsWithin];
    spring_ = new SpringModel(maxMaxAtomsWithin, grid, probe, atomCoords_);
}

HashGridVolume::~HashGridVolume() {
    for(int i=0; i < boxLx_; i++){
        for(int j=0; j < boxLy_; j++){
            delete [] hashGrid_[i][j];
        }
        delete [] hashGrid_[i];
    }
    delete [] hashGrid_; 
    delete [] gridChain_;
    delete [] di_;
    delete [] dj_;
    delete [] dk_;

    delete [] kLayerFromK_;
    delete [] vTypeFromTag_;
    
    for(int i=0; i <= ni; i++){
        for(int j=0; j <= nj; j++){
            delete [] gridLayer_[i][j];
            delete [] atomLayer_[i][j];
            delete [] closeAtomLayer_[i][j];
        }
        delete [] gridLayer_[i];
        delete [] atomLayer_[i];
        delete [] closeAtomLayer_[i];
    }
    delete [] gridLayer_; 
    delete [] atomLayer_; 
    delete [] closeAtomLayer_; 

    if( flagOpt_ == 1 ) {
       deleteOptimizedSearchPatterns();
    }

    delete atomsWithin_;
    delete spring_;
}
 
float HashGridVolume::getHashSpacing() {
    return hashSpacing_;
}

void HashGridVolume::writeShiftedCoordinates(string PDBFileName) {
    ifstream fIn;
    ofstream fOut;
    
    fIn.open(PDBFileName.c_str());
    fOut.open((PDBFileName + "_adjusted.pdb").c_str());
    
    PDBAtom * voidAtom = new PDBAtom;
    
    ifstream fin(PDBFileName.c_str());
    PDBProtein *protein = 0;
    protein = new PDBProtein;
    protein->read(fin);
    fin.close();
    PDBProtein::AtomIterator iter;   
    
    int count = 1;
    for (iter = protein->atom_begin(); iter != protein->atom_end(); ++iter) {
        voidAtom = *iter;
        voidAtom->x = atomCoords_[count][1];
        voidAtom->y = atomCoords_[count][2];
        voidAtom->z = atomCoords_[count][3];
        count = count + 1;
        voidAtom->write(fOut, false);
        fOut << "\n";
    }    
    
    fIn.close();
    fOut.close();   
    delete voidAtom;
}

void HashGridVolume::generateHashGrid(){    
    
    float xmin = FLT_MAX;
    float xmax = -FLT_MAX;
    float ymin = FLT_MAX;
    float ymax = -FLT_MAX;
    float zmin = FLT_MAX;
    float zmax = -FLT_MAX;

    for(int i=1; i<=nAtoms; i++){
        if (atomCoords_[i][1] > xmax) { xmax = atomCoords_[i][1]; }
        if (atomCoords_[i][1] < xmin) { xmin = atomCoords_[i][1]; }
        if (atomCoords_[i][2] > ymax) { ymax = atomCoords_[i][2]; }
        if (atomCoords_[i][2] < ymin) { ymin = atomCoords_[i][2]; }
        if (atomCoords_[i][3] > zmax) { zmax = atomCoords_[i][3]; }
        if (atomCoords_[i][3] < zmin) { zmin = atomCoords_[i][3]; }
    }

    xRange = xmax - xmin;
    yRange = ymax - ymin;
    zRange = zmax - zmin;
            
                                                                                    //shift coordinates so closest atom is at (x,y,z) = (dBuffer_,dBuffer_,dBuffer_)
    dBuffer_ = maxExtension_ + 0.5 * grid;                                      //distance of buffer for large-box
    float shift = maxExtension_ + dBuffer_; 
    xShift = (shift - xmin) / 2;
    yShift = (shift - ymin) / 2;
    zShift = (shift - zmin) / 2;
    for(int i=1; i <= nAtoms; i++){
        atomCoords_[i][1] = (atomCoords_[i][1] - xmin) + shift;
        atomCoords_[i][2] = (atomCoords_[i][2] - ymin) + shift;
        atomCoords_[i][3] = (atomCoords_[i][3] - zmin) + shift;
    }

                                                                                //account for boundary volume associated with solvent exposed surface of the protein
    xmax += maxExtension_;   
    ymax += maxExtension_; 
    zmax += maxExtension_; 
    xmin -= maxExtension_; 
    ymin -= maxExtension_; 
    zmin -= maxExtension_; 
   
    float xRangeLocal = xmax - xmin; 
    float yRangeLocal = ymax - ymin; 
    float zRangeLocal = zmax - zmin; 

    ni = ceil(xRangeLocal/grid) + 1;                                            //the extra 1 will be split in half on both sides of the number line
    nj = ceil(yRangeLocal/grid) + 1;                                            //the extra 1 will be split in half on both sides of the number line
    nk = ceil(zRangeLocal/grid) + 1;                                            //the extra 1 will be split in half on both sides of the number line
                                                                                //must make nk even to handle processSlice along the z-direction
    int nk0 = nk + 1;                                                           //temporarily add 1 more and see what happens: (even => +0 , odd => +1)
    nk = 2*(nk0/2);                                                             //forces nk to be even, which may extend it by 1 step. 
        
    xmax = (ni - 1)*grid + dBuffer_;                                            //without left-side buffer (need only to reach to ni - 1)
    ymax = (nj - 1)*grid + dBuffer_;                                            //without left-side buffer (need only to reach to nj - 1)
    zmax = (nk - 1)*grid + dBuffer_;                                            //without left-side buffer (need only to reach to nk - 1)
    boxLx_ = 2 + (int)floor(xmax/hashSpacing_);                                 //adding 2 is the right-side buffer
    boxLy_ = 2 + (int)floor(ymax/hashSpacing_);                                 //adding 2 is the right-side buffer
    boxLz_ = 2 + (int)floor(zmax/hashSpacing_);                                 //adding 2 is the right-side buffer
    
                                                                                //These are the largest indices but we must include 0 and boxL?_
    boxLx_++;                                                                   //add 1 so the range becomes i=0 to i<boxLx_ instead of i<=boxLx_
    boxLy_++;
    boxLz_++;
    
    hashGrid_ = (int ***) new int** [boxLx_]; 
    for(int i = 0; i < boxLx_; i++){
        hashGrid_[i] = (int **) new int* [boxLy_];  
        for(int j = 0; j < boxLy_; j++){
            hashGrid_[i][j] = (int *) new int [boxLz_];
            for(int k = 0; k < boxLz_; k++){
                hashGrid_[i][j][k] = 0;
            }
        }
    }    
    
    gridChain_ = new int[nAtoms + 1];                                           //not using index 0, starting at atom 1, according to atomCoords array
    for(int iter = 0; iter <= nAtoms; iter++){
        gridChain_[iter] = 0;
    }
    
    int i; 
    int j; 
    int k; 
    int chainLoc; 

    xmin_ = FLT_MAX; 
    ymin_ = FLT_MAX; 
    zmin_ = FLT_MAX; 
    xmax_ = -FLT_MAX; 
    ymax_ = -FLT_MAX; 
    zmax_ = -FLT_MAX; 

    for(int iter=1; iter <= nAtoms; iter++){                                    //initialize hashGrid_ and chainLoc
        i = (int)floor(atomCoords_[iter][1] / hashSpacing_);
        j = (int)floor(atomCoords_[iter][2] / hashSpacing_);
        k = (int)floor(atomCoords_[iter][3] / hashSpacing_);
        if (hashGrid_[i][j][k] == 0) {
            hashGrid_[i][j][k] = iter;
        } else {
            chainLoc = hashGrid_[i][j][k];
            while(gridChain_[chainLoc] != 0) {
                chainLoc = gridChain_[chainLoc];
            }
            gridChain_[chainLoc] = iter;
        }

        if (atomCoords_[iter][1] > xmax_) { xmax_ = atomCoords_[iter][1]; }     //determine atom location limits
        if (atomCoords_[iter][1] < xmin_) { xmin_ = atomCoords_[iter][1]; } 
        if (atomCoords_[iter][2] > ymax_) { ymax_ = atomCoords_[iter][2]; } 
        if (atomCoords_[iter][2] < ymin_) { ymin_ = atomCoords_[iter][2]; } 
        if (atomCoords_[iter][3] > zmax_) { zmax_ = atomCoords_[iter][3]; } 
        if (atomCoords_[iter][3] < zmin_) { zmin_ = atomCoords_[iter][3]; } 
    }

    xmin = grid + dBuffer_;                                                     //calculate first and last fine grid points 
    ymin = grid + dBuffer_;
    zmin = grid + dBuffer_;
    int iFirst = (int) floor(xmin/hashSpacing_);
    int jFirst = (int) floor(ymin/hashSpacing_);
    int kFirst = (int) floor(zmin/hashSpacing_);
    xmax = (ni-1) * grid + dBuffer_;
    ymax = (nj-1) * grid + dBuffer_;
    zmax = (nk-1) * grid + dBuffer_;
    int iLast = (int) floor(xmax/hashSpacing_);
    int jLast = (int) floor(ymax/hashSpacing_);
    int kLast = (int) floor(zmax/hashSpacing_);

                                                                                //test for errors
    if (iFirst != 2) {                                                          //lower limits are well defined
        cout << "       " << endl;
        cout << "ABORT: " << endl;
        cout << "ERROR detected in HashGridVolume setup, with iFirst not equal to 2" << endl;
        exit(-1); 
    }
    if( jFirst != 2 ) {
        cout << "       " << endl;
        cout << "ABORT: " << endl;
        cout << "ERROR detected in HashGridVolume setup, with jFirst not equal to 2" << endl;
        exit(-1); 
    }
    if( kFirst != 2 ) {
        cout << "       " << endl;
        cout << "ABORT: " << endl;
        cout << "ERROR detected in HashGridVolume setup, with kFirst not equal to 2" << endl;
        exit(-1); 
    }
                                                                                //upper limits
    if( iLast > (boxLx_ - 3) ) {
        cout << "       " << endl;
        cout << "ABORT: " << endl;
        cout << "ERROR detected in HashGridVolume setup, with iLast > (boxLx_ - 3)" << endl;
        exit(-1); 
    }
    if( jLast > (boxLy_ - 3) ) {
        cout << "       " << endl;
        cout << "ABORT: " << endl;
        cout << "ERROR detected in HashGridVolume setup, with jLast > (boxLy_ - 3)" << endl;
        exit(-1); 
    }
    if( kLast > (boxLz_ - 3) ) {
        cout << "       " << endl;
        cout << "ABORT: " << endl;
        cout << "ERROR detected in HashGridVolume setup, with kLast > (boxLz_ - 3)" << endl;
        exit(-1); 
    }
                                                                                //determine limits for optimized search
                                                                                
    iL0_ = 1;                                                                   //define lower bins
    iL1_ = (int)floor( (xmin_ - dBuffer_)/grid );
    i = (int)floor( (dBuffer_ + iL1_*grid)/hashSpacing_ );
    while (i <= iFirst) {
        iL1_++; 
        i = (int) floor( (dBuffer_ + iL1_*grid)/hashSpacing_ );
        if (i > iFirst) {
            iL1_--;
            break;
        }
    }
    iM0_ = iL1_ + 1;

    jL0_ = 1;
    jL1_ = (int)floor( (ymin_ - dBuffer_)/grid );
    j = (int)floor( (dBuffer_ + jL1_*grid)/hashSpacing_ );
       while( j <= jFirst ) {
        jL1_++; 
        j = (int)floor( (dBuffer_ + jL1_*grid)/hashSpacing_ );
        if( j > jFirst ) {
            jL1_--;
            break;
        }
    }
    jM0_ = jL1_ + 1;

    kL0_ = 1;
    kL1_ = (int)floor( (zmin_ - dBuffer_)/grid );
    k = (int)floor( (dBuffer_ + kL1_*grid)/hashSpacing_ );
        while (k <= kFirst) {
        kL1_++; 
        k = (int)floor( (dBuffer_ + kL1_*grid)/hashSpacing_ );
        if( k > kFirst ) {
            kL1_--;
            break;
        }
    }
    kM0_ = kL1_ + 1;
                                                                                //define upper bins
    iU1_ = ni - 1;
    iU0_ = (int)ceil( (xmax_ - dBuffer_)/grid );
    i = (int)floor( (dBuffer_ + iU0_*grid)/hashSpacing_ );
    while (i >= iLast) {
        iU0_--; 
        i = (int)floor( (dBuffer_ + iU0_*grid)/hashSpacing_ );
        if( i < iLast ) {
            iU0_++;
            break;
        }
    }
    iM1_ = iU0_ - 1;

    jU1_ = nj - 1;
    jU0_ = (int)ceil( (ymax_ - dBuffer_)/grid );
    j = (int)floor( (dBuffer_ + jU0_*grid)/hashSpacing_ );
    while (j >= jLast) {
        jU0_--; 
        j = (int)floor( (dBuffer_ + jU0_*grid)/hashSpacing_ );
        if (j < jLast) {
            jU0_++;
            break;
        }
    }
    jM1_ = jU0_ - 1;

    kU1_ = nk - 1;
    kU0_ = (int)ceil( (zmax_ - dBuffer_)/grid );
    k = (int)floor( (dBuffer_ + kU0_*grid)/hashSpacing_ );
    while (k >= kLast) {
        kU0_--; 
        k = (int)floor( (dBuffer_ + kU0_*grid)/hashSpacing_ );
        if( k < kLast ) {
            kU0_++;
            break;
        }
    }
    kM1_ = kU0_ - 1;    
}

void HashGridVolume::deleteOptimizedSearchPatterns() {
    delete [] diLLL_, djLLL_, dkLLL_;
    delete [] diLML_, djLML_, dkLML_;
    delete [] diLUL_, djLUL_, dkLUL_;
    delete [] diMLL_, djMLL_, dkMLL_;
    delete [] diMML_, djMML_, dkMML_;
    delete [] diMUL_, djMUL_, dkMUL_;
    delete [] diULL_, djULL_, dkULL_;
    delete [] diUML_, djUML_, dkUML_;
    delete [] diUUL_, djUUL_, dkUUL_;

    delete [] diLLM_, djLLM_, dkLLM_;
    delete [] diLMM_, djLMM_, dkLMM_;
    delete [] diLUM_, djLUM_, dkLUM_;
    delete [] diMLM_, djMLM_, dkMLM_;
                                                                                //delete [] diMMM_, djMMM_, dkMMM_;  <-- no need for this case
    delete [] diMUM_, djMUM_, dkMUM_;
    delete [] diULM_, djULM_, dkULM_;
    delete [] diUMM_, djUMM_, dkUMM_;
    delete [] diUUM_, djUUM_, dkUUM_;

    delete [] diLLU_, djLLU_, dkLLU_;
    delete [] diLMU_, djLMU_, dkLMU_;
    delete [] diLUU_, djLUU_, dkLUU_;
    delete [] diMLU_, djMLU_, dkMLU_;
    delete [] diMMU_, djMMU_, dkMMU_;
    delete [] diMUU_, djMUU_, dkMUU_;
    delete [] diULU_, djULU_, dkULU_;
    delete [] diUMU_, djUMU_, dkUMU_;
    delete [] diUUU_, djUUU_, dkUUU_;

}

int HashGridVolume::gridCount(int iStart, int iStop, int jStart, int jStop, int kStart, int kStop, int count) {
    for (int i = iStart; i <= iStop; i++) {
        for (int j = jStart; j <= jStop; j++) {
            for (int k = kStart; k <= kStop; k++) {
                int flag = 0;
                if (i != 0 ) { flag++;} 
                else if (j != 0) { flag++; } 
                else if (k != 0) { flag++; } 
                  if( flag > 0 ) {
                    di_[count] = i;
                    dj_[count] = j;
                    dk_[count] = k;                         
                    count++;
                }
            }
        }
    }
    return count;
}

int HashGridVolume::gridCount2(int iStart, int iStop, int jStart, int jStop, int kStart, int kStop, int count) {
    for (int i = iStart; i <= iStop; i++) {
        for (int j = jStart; j <= jStop; j++) {
            for (int k = kStart; k <= kStop; k++) {
                int flag = 0;
                if (abs(i) >  1) { flag++; } 
                else if (abs(j) > 1) { flag++; } 
                else if (abs(k) > 1) { flag++; } 
                if (flag > 0)  {
                    di_[count] = i;
                    dj_[count] = j;
                    dk_[count] = k;                         
                    count++;
                }
            }
        }
    }
    return count;
}




// ------------------------------------------------- break triple sum into 27 sets of partial summations
// L = lower  range of index defined by iL0_ to iL1_, jL0_ to jL1_, kL0_ to kL1_ for x,y,z-directions  
// M = middle range of index defined by iM0_ to iM1_, jM0_ to jM1_, kM0_ to kM1_ for x,y,z-directions  
// U = upper  range of index defined by iU0_ to iU1_, jU0_ to jU1_, kU0_ to kU1_ for x,y,z-directions  
// sumX * sumY * sumZ = (sumXL + sumXM + sumXU) * (sumYL + sumYM + sumYU) * (sumZL + sumZM + sumZU) 
//
//                    = (sumXL + sumXM + sumXU) * (sumYL + sumYM + sumYU) * sumZL 
//                    + (sumXL + sumXM + sumXU) * (sumYL + sumYM + sumYU) * sumZM 
//                    + (sumXL + sumXM + sumXU) * (sumYL + sumYM + sumYU) * sumZU 
//
//                    = sumLLL + sumLML + sumLUL 
//                    + sumMLL + sumMML + sumMUL 
//                    + sumULL + sumUML + sumUUL 
//                      ------------------------
//                    + sumLLM + sumLMM + sumLUM 
//                    + sumMLM + sumMMM + sumMUM 
//                    + sumULM + sumUMM + sumUUM 
//                      ------------------------
//                    + sumLLU + sumLMU + sumLUU 
//                    + sumMLU + sumMMU + sumMUU 
//                    + sumULU + sumUMU + sumUUU 
// Recall number line: 
//
//  Example Number Line: All directions work the same. 
//  |<------ buffer ------>|
//                 i=0     1     2     3     4     5  ... ni-6  ni-5  ni-4  ni-3  ni-2  ni-1  ni
//                   |     |     |     |     |     |       |     |     |     |     |     |     |      
//  BBBBBBBBBBBBBBBBBBBBEEEEEEEEEEEEEEEEEEEE--------- ... -------------EEEEEEEEEEEEEEEEEEEEBBBBBBBBBBBBBBBBBBBB
//  |                   |<-------------------------- xRange ----------------------------->|
//  x=0                xMin
// ---------------------------------------------------------- determination of limits in x, y, z -directions
   
//  at x=0 defines the beginning of the first quadrant.  
//  xRange = inclusive of all the physics in the problem (atoms and boundary volume). 
//  E = extended range to account for boundary volume
//  B = buffer region that is added so that the hash grid never requires periodic translations.
//  xMin = maxExtension_ + 0.5*grid;    closest point to origin in 1st quadrant defines edge of large-box (system)
//  yMin = xMin;                        symmetry applies to x,y,z. Closest point is along body diagonal of large-box 
//  zMin = xMin;                        call xMin = yMin = zMin = dBuffer_ = distance of buffer for 1st quadrant
 

// ------------------------------------------------------------------------ z-direction is L
// part I:   sumXL * sumYL => nnxtLLL_  n2ndLLL_  diLLL_  djLLL_  dkLLL_            case LLL
//           sumXL * sumYM => nnxtLML_  n2ndLML_  diLML_  djLML_  dkLML_                 LML
//           sumXL * sumYU => nnxtLUL_  n2ndLUL_  diLUL_  djLUL_  dkLUL_                 LUL
//
// part II:  sumXM * sumYL => nnxtMLL_  n2ndMLL_  diMLL_  djMLL_  dkMLL_                 MLL
//           sumXM * sumYM => nnxtMML_  n2ndMML_  diMML_  djMML_  dkMML_                 MML
//           sumXM * sumYU => nnxtMUL_  n2ndMUL_  diMUL_  djMUL_  dkMUL_                 MUL
//
// part III: sumXU * sumYL => nnxtULL_  n2ndULL_  diULL_  djULL_  dkULL_                 ULL
//           sumXU * sumYM => nnxtUML_  n2ndUML_  diUML_  djUML_  dkUML_                 UML
//           sumXU * sumYU => nnxtUUL_  n2ndUUL_  diUUL_  djUUL_  dkUUL_                 UUL
// ------------------------------------------------------------------------ z-direction is M
// part I:   sumXL * sumYL => nnxtLLM_  n2ndLLM_  diLLM_  djLLM_  dkLLM_            case LLM
//           sumXL * sumYM => nnxtLMM_  n2ndLMM_  diLMM_  djLMM_  dkLMM_                 LMM
//           sumXL * sumYU => nnxtLUM_  n2ndLUM_  diLUM_  djLUM_  dkLUM_                 LUM
//
// part II:  sumXM * sumYL => nnxtMLM_  n2ndMLM_  diMLM_  djMLM_  dkMLM_                 MLM
//           sumXM * sumYM => nnxtMMM_  n2ndMMM_  diMMM_  djMMM_  dkMMM_                 MMM  <------ no need to do this one
//           sumXM * sumYU => nnxtMUM_  n2ndMUM_  diMUM_  djMUM_  dkMUM_                 MUM
//
// part III: sumXU * sumYL => nnxtULM_  n2ndULM_  diULM_  djULM_  dkULM_                 ULM
//           sumXU * sumYM => nnxtUMM_  n2ndUMM_  diUMM_  djUMM_  dkUMM_                 UMM
//           sumXU * sumYU => nnxtUUM_  n2ndUUM_  diUUM_  djUUM_  dkUUM_                 UUM
// ------------------------------------------------------------------------ z-direction is U
// part I:   sumXL * sumYL => nnxtLLU_  n2ndLLU_  diLLU_  djLLU_  dkLLU_            case LLU
//           sumXL * sumYM => nnxtLMU_  n2ndLMU_  diLMU_  djLMU_  dkLMU_                 LMU
//           sumXL * sumYU => nnxtLUU_  n2ndLUU_  diLUU_  djLUU_  dkLUU_                 LUU
//
// part II:  sumXM * sumYL => nnxtMLU_  n2ndMLU_  diMLU_  djMLU_  dkMLU_                 MLU
//           sumXM * sumYM => nnxtMMU_  n2ndMMU_  diMMU_  djMMU_  dkMMU_                 MMU
//           sumXM * sumYU => nnxtMUU_  n2ndMUU_  diMUU_  djMUU_  dkMUU_                 MUU
//
// part III: sumXU * sumYL => nnxtULU_  n2ndULU_  diULU_  djULU_  dkULU_                 ULU
//           sumXU * sumYM => nnxtUMU_  n2ndUMU_  diUMU_  djUMU_  dkUMU_                 UMU
//           sumXU * sumYU => nnxtUUU_  n2ndUUU_  diUUU_  djUUU_  dkUUU_                 UUU
//
// ------------------------------------------------------------------------ case LLL
    


void HashGridVolume::buildGridSearchPatterns(){    
                                                                                //allocate memory to hold complete displacement vector set
    di_ = new int[125];
    dj_ = new int[125];
    dk_ = new int[125];
    di_[0] = 0;
    dj_[0] = 0;
    dk_[0] = 0;        
    
    if (flagOpt_ == 1) {                                                        //build displacement vectors for optimized boundary conditions as subsets 
                                                                                //  of the complete set       
        nnxtLLL_ = gridCount(0, 1, 0, 1, 0, 1, 1);
        n2ndLLL_ = gridCount2(0, 2, 0, 2, 0, 2, nnxtLLL_);    
        diLLL_ = new int[n2ndLLL_];
        djLLL_ = new int[n2ndLLL_];
        dkLLL_ = new int[n2ndLLL_];
        for (int c=0; c < n2ndLLL_; c++) {
            diLLL_[c] = di_[c];
            djLLL_[c] = dj_[c];
            dkLLL_[c] = dk_[c];
        }
    
        nnxtLML_ = gridCount(0, 1, -1, 1, 0, 1, 1);
        n2ndLML_ = gridCount2(0, 2, -2, 2, 0, 2, nnxtLML_);    
        diLML_ = new int[n2ndLML_];
        djLML_ = new int[n2ndLML_];
        dkLML_ = new int[n2ndLML_];
        for (int c=0; c<n2ndLML_; c++) {
            diLML_[c] = di_[c];
            djLML_[c] = dj_[c];
            dkLML_[c] = dk_[c];
        }

        nnxtLUL_ = gridCount(0, 1, -1, 0, 0, 1, 1);
        n2ndLUL_ = gridCount2(0, 2, -2, 0, 0, 2, nnxtLUL_);  
        diLUL_ = new int[n2ndLUL_];
        djLUL_ = new int[n2ndLUL_];
        dkLUL_ = new int[n2ndLUL_];
        for (int c=0; c<n2ndLUL_; c++) {
            diLUL_[c] = di_[c];
            djLUL_[c] = dj_[c];
            dkLUL_[c] = dk_[c];
        }

        nnxtMLL_ = gridCount(-1, 1, 0, 1, 0, 1, 1);
        n2ndMLL_ = gridCount2(-2, 2, 0, 2, 0, 2, nnxtMLL_);
        diMLL_ = new int[n2ndMLL_];
        djMLL_ = new int[n2ndMLL_];
        dkMLL_ = new int[n2ndMLL_];
        for (int c=0; c<n2ndMLL_; c++) {
            diMLL_[c] = di_[c];
            djMLL_[c] = dj_[c];
            dkMLL_[c] = dk_[c];
        }

        nnxtMML_ = gridCount(-1, 1, -1, 1, 0, 1, 1);
        n2ndMML_ = gridCount2(-2, 2, -2, 2, 0, 2, nnxtMML_);
        diMML_ = new int[n2ndMML_];
        djMML_ = new int[n2ndMML_];
        dkMML_ = new int[n2ndMML_];

        for (int c=0; c<n2ndMML_; c++) {
            diMML_[c] = di_[c];
            djMML_[c] = dj_[c];
            dkMML_[c] = dk_[c];
        }

        nnxtMUL_ = gridCount(-1, 1, -1, 0, 0, 1, 1);
        n2ndMUL_ = gridCount2(-2, 2, -2, 0, 0, 2, nnxtMUL_);
        diMUL_ = new int[n2ndMUL_];
        djMUL_ = new int[n2ndMUL_];
        dkMUL_ = new int[n2ndMUL_];
        for(int c=0; c<n2ndMUL_; c++) {
            diMUL_[c] = di_[c];
            djMUL_[c] = dj_[c];
            dkMUL_[c] = dk_[c];
        }

        nnxtULL_ = gridCount(-1, 0, 0, 1, 0, 1, 1);
        n2ndULL_ = gridCount2(2, 0, 0, 2, 0, 2, nnxtULL_);
        diULL_ = new int[n2ndULL_];
        djULL_ = new int[n2ndULL_];
        dkULL_ = new int[n2ndULL_];
        for(int c=0; c<n2ndULL_; c++) {
            diULL_[c] = di_[c];
            djULL_[c] = dj_[c];
            dkULL_[c] = dk_[c];
        }

        nnxtUML_ = gridCount(-1, 0, -1, 1, 0, 1, 1);
        n2ndUML_ = gridCount2(-2, 0, -2, 2, 0, 2, nnxtUML_);
        diUML_ = new int[n2ndUML_];
        djUML_ = new int[n2ndUML_];
        dkUML_ = new int[n2ndUML_];
        for(int c=0; c<n2ndUML_; c++) {
            diUML_[c] = di_[c];
            djUML_[c] = dj_[c];
            dkUML_[c] = dk_[c];
        }

        nnxtUUL_ = gridCount(-1, 0, -1, 0, 0, 1, 1);
        n2ndUUL_ = gridCount2(-2, 0, -2, 0, 0, 2, nnxtUUL_);
        diUUL_ = new int[n2ndUUL_];
        djUUL_ = new int[n2ndUUL_];
        dkUUL_ = new int[n2ndUUL_];
        for(int c=0; c<n2ndUUL_; c++) {
            diUUL_[c] = di_[c];
            djUUL_[c] = dj_[c];
            dkUUL_[c] = dk_[c];
        }

        nnxtLLM_ = gridCount(0, 1, 0, 1, -1, 1, 1);
        n2ndLLM_ = gridCount2(0, 2, 0, 2, -2, 2, nnxtLLM_);
        diLLM_ = new int[n2ndLLM_];
        djLLM_ = new int[n2ndLLM_];
        dkLLM_ = new int[n2ndLLM_];
        for(int c=0; c<n2ndLLM_; c++) {
            diLLM_[c] = di_[c];
            djLLM_[c] = dj_[c];
            dkLLM_[c] = dk_[c];
        }

        nnxtLMM_ = gridCount(0, 1, -1, 1, -1, 1, 1);
        n2ndLMM_ = gridCount2(0, 2, -2, 2, -2, 2, nnxtLMM_);
        diLMM_ = new int[n2ndLMM_];
        djLMM_ = new int[n2ndLMM_];
        dkLMM_ = new int[n2ndLMM_];
        for(int c=0; c<n2ndLMM_; c++) {
            diLMM_[c] = di_[c];
            djLMM_[c] = dj_[c];
            dkLMM_[c] = dk_[c];
        }
   
        nnxtLUM_ = gridCount(0, 1, -1, 0, -1, 1, 1);
        n2ndLUM_ = gridCount2(0, 2, -2, 0, -2, 2, nnxtLUM_);
        diLUM_ = new int[n2ndLUM_];
        djLUM_ = new int[n2ndLUM_];
        dkLUM_ = new int[n2ndLUM_];
        for(int c=0; c<n2ndLUM_; c++) {
            diLUM_[c] = di_[c];
            djLUM_[c] = dj_[c];
            dkLUM_[c] = dk_[c];
        }

        nnxtMLM_ = gridCount(-1, 1, 0, 1, -1, 1, 1);
        n2ndMLM_ = gridCount2(-2, 2, 0, 2, -2, 2, nnxtMLM_);
        diMLM_ = new int[n2ndMLM_];
        djMLM_ = new int[n2ndMLM_];
        dkMLM_ = new int[n2ndMLM_];
        for(int c=0; c<n2ndMLM_; c++) {
            diMLM_[c] = di_[c];
            djMLM_[c] = dj_[c];
            dkMLM_[c] = dk_[c];
        }
                                                                                //skip MMM

        nnxtMUM_ = gridCount(-1, 1, -1, 0, -1, 1, 1);
        n2ndMUM_ = gridCount2(-2, 2, -2, 0, -2, 2, nnxtMUM_);
        diMUM_ = new int[n2ndMUM_];
        djMUM_ = new int[n2ndMUM_];
        dkMUM_ = new int[n2ndMUM_];
        for(int c=0; c<n2ndMUM_; c++) {
            diMUM_[c] = di_[c];
            djMUM_[c] = dj_[c];
            dkMUM_[c] = dk_[c];
        }

        nnxtULM_ = gridCount(-1, 0, 0, 1, -1, 1, 1);
        n2ndULM_ = gridCount2(-2, 0, 0, 2, -2, 2, nnxtULM_);
        diULM_ = new int[n2ndULM_];
        djULM_ = new int[n2ndULM_];
        dkULM_ = new int[n2ndULM_];
        for(int c=0; c<n2ndULM_; c++) {
            diULM_[c] = di_[c];
            djULM_[c] = dj_[c];
            dkULM_[c] = dk_[c];
        }

        nnxtUMM_ = gridCount(-1, 0, -1, 1, -1, 1, 1);
        n2ndUMM_ = gridCount2(-2, 0, -2, 2, -2, 2, nnxtUMM_);
        diUMM_ = new int[n2ndUMM_];
        djUMM_ = new int[n2ndUMM_];
        dkUMM_ = new int[n2ndUMM_];
        for(int c=0; c<n2ndUMM_; c++) {
            diUMM_[c] = di_[c];
            djUMM_[c] = dj_[c];
            dkUMM_[c] = dk_[c];
        }

        nnxtUUM_ = gridCount(-1, 0, -1, 1, -1, 1, 1);
        n2ndUUM_ = gridCount2(-2, 0, -2, 2, -2, 2, nnxtUUM_);
        diUUM_ = new int[n2ndUUM_];
        djUUM_ = new int[n2ndUUM_];
        dkUUM_ = new int[n2ndUUM_];
        for(int c=0; c<n2ndUUM_; c++) {
            diUUM_[c] = di_[c];
            djUUM_[c] = dj_[c];
            dkUUM_[c] = dk_[c];
        }

        nnxtLLU_ = gridCount(0, 1, 0, 1, -1, 0, 1);
        n2ndLLU_ = gridCount2(0, 2, 0, 2, -2, 0, nnxtLLU_);
        diLLU_ = new int[n2ndLLU_];
        djLLU_ = new int[n2ndLLU_];
        dkLLU_ = new int[n2ndLLU_];
        for(int c=0; c<n2ndLLU_; c++) {
            diLLU_[c] = di_[c];
            djLLU_[c] = dj_[c];
            dkLLU_[c] = dk_[c];
        }

        nnxtLMU_ = gridCount(0, 1, -1, 1, -1, 0, 1);
        n2ndLMU_ = gridCount2(0, 2, -2, 2, -2, 0, nnxtLMU_);
        diLMU_ = new int[n2ndLMU_];
        djLMU_ = new int[n2ndLMU_];
        dkLMU_ = new int[n2ndLMU_];
        for(int c=0; c<n2ndLMU_; c++) {
            diLMU_[c] = di_[c];
            djLMU_[c] = dj_[c];
            dkLMU_[c] = dk_[c];
        }

        nnxtLUU_ = gridCount(0, 1, -1, 0, -1, 0, 1);
        n2ndLUU_ = gridCount2(0, 2, -2, 0, -2, 0, nnxtLUU_);
        diLUU_ = new int[n2ndLUU_];
        djLUU_ = new int[n2ndLUU_];
        dkLUU_ = new int[n2ndLUU_];
        for(int c=0; c<n2ndLUU_; c++) {
            diLUU_[c] = di_[c];
            djLUU_[c] = dj_[c];
            dkLUU_[c] = dk_[c];
        }

        nnxtMLU_ = gridCount(-1, 1, 0, 1, -1, 0, 1);
        n2ndMLU_ = gridCount2(-2, 2, 0, 2, 2, 0, nnxtMLU_);
        diMLU_ = new int[n2ndMLU_];
        djMLU_ = new int[n2ndMLU_];
        dkMLU_ = new int[n2ndMLU_];
        for(int c=0; c<n2ndMLU_; c++) {
            diMLU_[c] = di_[c];
            djMLU_[c] = dj_[c];
            dkMLU_[c] = dk_[c];
        }

        nnxtMMU_ = gridCount(-1, 1, -1, 1, -1, 0, 1);
        n2ndMMU_ = gridCount2(-2, 2, -2, 2, -2, 0, nnxtMMU_);
        diMMU_ = new int[n2ndMMU_];
        djMMU_ = new int[n2ndMMU_];
        dkMMU_ = new int[n2ndMMU_];
        for(int c=0; c<n2ndMMU_; c++) {
            diMMU_[c] = di_[c];
            djMMU_[c] = dj_[c];
            dkMMU_[c] = dk_[c];
        }

        nnxtMUU_ = gridCount(-1, 1, -1, 0, -1, 0, 1);
        n2ndMUU_ = gridCount2(-2, 2, -2, 0, -2, 0, nnxtMUU_);
        diMUU_ = new int[n2ndMUU_];
        djMUU_ = new int[n2ndMUU_];
        dkMUU_ = new int[n2ndMUU_];
        for(int c=0; c<n2ndMUU_; c++) {
            diMUU_[c] = di_[c];
            djMUU_[c] = dj_[c];
            dkMUU_[c] = dk_[c];
        }

        nnxtULU_ = gridCount(-1, 0, 0, 1, -1, 0, 1);
        n2ndULU_ = gridCount2(-2, 0, 0, 2, -2, 0, nnxtULU_);
        diULU_ = new int[n2ndULU_];
        djULU_ = new int[n2ndULU_];
        dkULU_ = new int[n2ndULU_];
        for(int c=0; c<n2ndULU_; c++) {
            diULU_[c] = di_[c];
            djULU_[c] = dj_[c];
            dkULU_[c] = dk_[c];
        }

        nnxtUMU_ = gridCount(-1, 0, -1, 1, -1, 0, 1);
        n2ndUMU_ = gridCount2(-2, 0, -2, 2, -2, 0, nnxtUMU_);
        diUMU_ = new int[n2ndUMU_];
        djUMU_ = new int[n2ndUMU_];
        dkUMU_ = new int[n2ndUMU_];
        for(int c=0; c<n2ndUMU_; c++) {
            diUMU_[c] = di_[c];
            djUMU_[c] = dj_[c];
            dkUMU_[c] = dk_[c];
        }
        
        nnxtUUU_ = gridCount(-1, 0, -1, 0, -1, 0, 1);
        n2ndUUU_ = gridCount2(-2, 0, -2, 0, -2, 0, nnxtUUU_);
        diUUU_ = new int[n2ndUUU_];
        djUUU_ = new int[n2ndUUU_];
        dkUUU_ = new int[n2ndUUU_];
        for(int c=0; c<n2ndUUU_; c++) {
            diUUU_[c] = di_[c];
            djUUU_[c] = dj_[c];
            dkUUU_[c] = dk_[c];
        }

    }
 
    int c = gridCount(-1, 1, -1, 1, -1, 1, 1);
    n1st_ = 27;
    nnxt_ = 27;
    gridCount2(-2, 2, -2, 2, -2, 2, c);
    n2nd_ = 125;
}


void HashGridVolume::generateGridLayer(){    
    float width = probeDiameter_;                                               //width of layer is initially defined as maximum range of object of interest
    width += sqrt(3)*grid;                                                      //sqrt(3)*grid accounts for round off error to nearest gridpoint
    float f = 0.01 + width/grid;                                                //number of grids within width (0.01 adds good measure before rounding down)
    kLayer_ = (int)floor(f);                                                    //probe radius is used (current implementation)
    kLayer_ += 1;                                                               //ensures the entire range of width is covered inclusively along the z-direction
    kLayer_ += 1;                                                               //ensures a roundoff of 1 grid distance in the +z-direction is still covered
    int maxIndex = kLayer_;                                                     //this limit ensures the new information will never interfer with kLneg_
    kLayer_ += 1;                                                               //accounts for the unknown slice, defined by kLpos_
    kLayer_ += 1;                                                               //accounts for the current known slice, defined by kLmid_
    kLayer_ += 1;                                                               //accounts for the past slice defined by kLneg_
                                                                                //kLneg  kLmid  KLpos  +1  +2  +3  +4  +5  +6  +7  +8  (say when probeDiameter=2.8 and grid = 0.5)
    
    tag_ = 0;                                                                   //initial value defines final result: vType = gridLayer_[i][j][kL] - refTag
    minTag_ = -kLayer_* dTag_;                                                  //indicates old information in gridLayer whenever (tag_ - minTag_) < 0; 
    minTag_ += dTag_;                                                           //must add dTag_ so that only the kLpos_ is unknown while all other slices are known
    refTag_ = tag_ + 1;                                                         //refTag is always exactly 1 more more than tag_ 
                                                                                //gridpoint type: type = -1,0,1,2 => microvoid, protein-volume, void-not-bulk, void-bulk
  
/*       
 -------------------------------- the order of {ni,nj,kLayer_} versus {klayer_,nj,ni} could substantially affect CPU time
                                  klayer_ is typically very small compared to nj & ni
                                  use kLayer_ last so neighboring slices within the layer can be accessed quickly
                                  Is this assumption valid? Perhaps the reverse is better?
 */
    
    
    initTag_ = minTag_ - 1;                                                      //initial tag
    gridLayer_ = (int ***) new int** [ni+1];                                    // must add 1 in limit because we need to go from 0 to ni inclusive
    atomLayer_ = (int ***) new int** [ni+1];                                    // must add 1 in limit because we need to go from 0 to ni inclusive
    closeAtomLayer_ = (int ***) new int** [ni+1];        
    for (int i = 0; i <= ni; i++) {
        gridLayer_[i] = (int **) new int* [nj+1];                               // must add 1 in limit because we need to go from 0 to nj inclusive
        atomLayer_[i] = (int **) new int* [nj+1];                               // must add 1 in limit because we need to go from 0 to nj inclusive
        closeAtomLayer_[i] = (int **) new int* [nj+1]; 
        for(int j = 0; j <= nj; j++){
            gridLayer_[i][j] = (int *) new int [kLayer_];                       //the +1 was added above for the limit to go from 0 to < kLayer_
            atomLayer_[i][j] = (int *) new int [kLayer_];                       //the +1 was added above for the limit to go from 0 to < kLayer_
            closeAtomLayer_[i][j] = (int *) new int [kLayer_];  
            for(int kL = 0; kL < kLayer_; kL++){
                gridLayer_[i][j][kL] = initTag_;                                //intialize gridlayer not to be assigned anywhere
            }
        }
    }    
    int setTag = 3 + (nk + 5) * dTag_;                                          //This tag will never be reached because it is 4 more than last expected tag.    
                                                                                
    for(int i = 0; i <= ni; i++){                                               //bottom of box, initialize to completely bulk
        for(int j = 0; j <= nj; j++){
            gridLayer_[i][j][0] = 3;
            atomLayer_[i][j][0] = -1;                                           //void-bulk  => no atom is in range
            closeAtomLayer_[i][j][0] = -1; 
        }
     }    
                                                                                    
    int n0 = 0;                                                                 //initialize all boundary edges of box to bulk
    for(int i = 0; i <= ni; i++){
        for(int kL = 0; kL < kLayer_; kL++){
            gridLayer_[i][n0][kL] = setTag;
            atomLayer_[i][n0][kL] = -1;                                         //void-bulk => no atom is in range
            closeAtomLayer_[i][n0][kL] = -1; 
            gridLayer_[i][nj][kL] = setTag;
            atomLayer_[i][nj][kL] = -1;                                         //void-bulk => no atom is in range
            closeAtomLayer_[i][nj][kL] = -1;  
        }
     }
    for(int j = 0; j <= nj; j++){
        for(int kL = 0; kL < kLayer_; kL++){
            gridLayer_[n0][j][kL] = setTag;
            atomLayer_[n0][j][kL] = -1;                                         //void-bulk => no atom is in range
            closeAtomLayer_[n0][j][kL] = -1;    
            gridLayer_[ni][j][kL] = setTag;
            atomLayer_[ni][j][kL] = -1;                                         //void-bulk => no atom is in range
            closeAtomLayer_[ni][j][kL] = -1;   
        }
     }
                                                                                //precalculate MOD functions
    int   nk1 = nk + 1;                                                         //will need to go one more than last requested value  
    kLayerFromK_ = new int[nk1+1];
    vTypeFromTag_ = new int[setTag+1];
    for(int k=0; k <= nk1; k++) {
       int kL = k % kLayer_; 
       kLayerFromK_[k] = kL; 

       int vMicroVoid = tag_;
       int vProteinVol = tag_ + 1;                                              //vProteinVol - refTag => 0 = (tag % dtag_) - 1
       int vCavity = tag_ + 2;                                                  //vCavity - refTag     => 1 = (tag % dtag_) - 1
       int vBulk = tag_ + 3;                                                    //vBulk - refTag       => 2 = (tag % dtag_) - 1

       int vType = (vMicroVoid  % dTag_) - 1; 
       vTypeFromTag_[tag_] = vType;
       vType = (vProteinVol % dTag_) - 1; 
       vTypeFromTag_[tag_+1] = vType;
       vType = (vCavity     % dTag_) - 1; 
       vTypeFromTag_[tag_+2] = vType;
       vType = (vBulk       % dTag_) - 1; 
       vTypeFromTag_[tag_+3] = vType;
       tag_ += dTag_;   
    }
    tag_ = 0;      
    
   
}

void HashGridVolume::volTypeSection(int ix, int jy, int kz, int c1, int c2, int * di, int * dj, int * dk) {
    
        
    register float vdw;                                                         
    
    float x = ix*grid + dBuffer_;
    float y = jy*grid + dBuffer_;
    float z = kz*grid + dBuffer_;

    int i = (int) (x / hashSpacing_); 
    int j = (int) (y / hashSpacing_);
    int k = (int) (z / hashSpacing_);
    
    register int atom;                                                         
    int maxAtomsWithin = 0;                                                     // counts the number of atoms needed for the spring model
    int lastAtom = 0;    
    
    closeAtom = -1;                                                             // flags grid point as a void within bulk, not associated with an atom
    veryCloseAtom = -1;    
    float   closeDistance = 1000000;
    float   veryCloseDistance = 1000000;    
    float   distanceToSurface;
                                                                                
    for (int c = 0; c < c1 ; c++) {                                             // search 1st level of neighbors to check if spring model is necessary    
        int i2 = i + di[c];                                                     // Note i,j,k for hash grid,  ix, iy, kLpos_ for fine grid
        int j2 = j + dj[c];
        int k2 = k + dk[c];

        atom = hashGrid_[i2][j2][k2];            
        while (atom != 0){            
            vdw = atomCoords_[atom][0];                        
            float x1 = x - atomCoords_[atom][1];
            float y1 = y - atomCoords_[atom][2];
            float z1 = z - atomCoords_[atom][3];
            float sqDist2 = x1 * x1 + y1 * y1 + z1 * z1;        
            if (sqDist2 <= vdw * vdw) {                                         //gridpoint is in protein volume: neighborhood of atom
                gridLayer_[ix][jy][kLpos_] = refTag_;                           //vType = 0  => protein volume 
                closeAtom = atom;
                atomLayer_[ix][jy][kLpos_] = closeAtom;     
                return;                                                         //return, no other processing is necessary
            } else {         
                distanceToSurface = sqrt(sqDist2) - vdw;                        //must count all atoms within compressive spring range                                                        
                if (distanceToSurface <= shell) {                               //condition to determine if gridpoint is bulk or not-bulk  
                    if (distanceToSurface < closeDistance) {
                        closeDistance = distanceToSurface;
                        closeAtom = atom;
                    }                    
                    if (distanceToSurface <= probeDiameter_) {    
                        if (distanceToSurface <= cutoff) {
                            if (distanceToSurface < veryCloseDistance) {
                                veryCloseDistance = distanceToSurface;
                                veryCloseAtom = atom;
                            }
                        }                    
                        atomsWithin_[maxAtomsWithin++] = atom;                         
                        if (maxAtomsWithin > maxMaxAtomsWithin) {
                            cout << "ERROR: Require greater maximum number of atoms within" << endl;
                            cout << "       maxMaxAtomsWithin = " << maxMaxAtomsWithin << endl;
                            cout << "       increase maxMaxAtomsWithin in HashGridVolume.hpp " << endl;
                            exit(-1); 
                        }
                    } 
                }              
            }                                                                   //tracking closest atom requires checking all cases           
            lastAtom = atom;
            atom = gridChain_[atom];
        }
    }
                                                                                //check if closest atom has already been identified
    if (closeAtom < 0) {                                                        //search over 2nd level of neighbors for closest atom only
        for (int c = c1; c < c2; c++) {  
            int i2 = i + di[c];                                                 //Note i,j,k for hash grid,  ix, iy, kLpos_ for fine grid
            int j2 = j + dj[c];
            int k2 = k + dk[c];
            atom = hashGrid_[i2][j2][k2]; 
            while (atom != 0) {                
                float x1 = x - atomCoords_[atom][1];      
                float y1 = y - atomCoords_[atom][2];
                float z1 = z - atomCoords_[atom][3];                
                float sqDist2 = x1 * x1 + y1 * y1 + z1 * z1;     
                distanceToSurface = sqrt(sqDist2) - atomCoords_[atom][0]; 
                if (distanceToSurface <= shell){                                // checks if grid point is bulk or not-bulk
                    if (distanceToSurface < closeDistance) {
                        closeDistance = distanceToSurface;
                        closeAtom = atom;              
                    }
                    if (distanceToSurface <= cutoff) {
                        if (distanceToSurface < veryCloseDistance) {
                            veryCloseDistance = distanceToSurface;
                            veryCloseAtom = atom;
                        }
                    }
                }                
                lastAtom = atom;
                atom = gridChain_[atom];
            }
        }
    }
    
    if (maxAtomsWithin > skipSpringMethod_) {                                   //must invoke spring model to test grid point
        int type = spring_->testSpring(x, y, z, atomsWithin_, maxAtomsWithin);  //type = -1 or 1 for microvoid or void-not-bulk
        gridLayer_[ix][jy][kLpos_] = type + refTag_;                            //Note i,j,k is for hash grid, ix,jy,kLpos_ for fine grid
        atomLayer_[ix][jy][kLpos_] = closeAtom;       
        closeAtomLayer_[ix][jy][kLpos_] = veryCloseAtom;                
    } else if (closeAtom > 0) {                                                 // associated with an atom
        gridLayer_[ix][jy][kLpos_] = 1 + refTag_;                               // vType = 1   void-not-bulk
        atomLayer_[ix][jy][kLpos_] = closeAtom;        
        closeAtomLayer_[ix][jy][kLpos_] = veryCloseAtom;        
    } else {                                                                    //closestAtom = -1 => grid point is out of range from any neighboring atom
        gridLayer_[ix][jy][kLpos_] = 2 + refTag_;                               //vType = 2   void-bulk
    } 
}

void HashGridVolume::CompleteSliceIn27Sections() {
                                                                                //break into 27 cases divided as: (3 by k_) each of which has 9 cases 
                                                                                //bottom section of box (do not need kL0_ or kL1_)
    if (k_ < kM0_) {                                                            //part I: sumXL * (sumYL + sumYM + sumYU)
        for (int i = iL0_; i <= iL1_; i++) {                                    //part Ia: sumXL * sumYL                                                                                
            for (int j = jL0_; j <= jL1_; j++) {
                volTypeSection(i,j,k_,nnxtLLL_,n2ndLLL_,diLLL_,djLLL_,dkLLL_);  
            }
            for (int j = jM0_; j <= jM1_; j++) {                                //part Ib: sumXL * sumYM
                volTypeSection(i,j,k_,nnxtLML_,n2ndLML_,diLML_,djLML_,dkLML_); 
            }
            for (int j = jU0_; j <= jU1_; j++){                                 //part Ic: sumXL * sumYU
                volTypeSection(i,j,k_,nnxtLUL_,n2ndLUL_,diLUL_,djLUL_,dkLUL_); 
            }
        }
        for (int i = iM0_; i <= iM1_; i++) {                                    //part II: sumXM * (sumYL + sumYM + sumYU)
            for (int j = jL0_; j <= jL1_; j++) {                                //part IIa: sumXM * sumYL
                volTypeSection(i,j,k_,nnxtMLL_,n2ndMLL_,diMLL_,djMLL_,dkMLL_); 
            }
            for (int j = jM0_; j <= jM1_; j++){                                 //part IIb: sumXM * sumYM
                volTypeSection(i,j,k_,nnxtMML_,n2ndMML_,diMML_,djMML_,dkMML_); 
            }
            for (int j = jU0_; j <= jU1_; j++){                                 //part IIc: sumXM * sumYU
                volTypeSection(i,j,k_,nnxtMUL_,n2ndMUL_,diMUL_,djMUL_,dkMUL_); 
            }
        }
        for(int i = iU0_; i <= iU1_; i++){                                      //part III: sumXU * (sumYL + sumYM + sumYU)
            for(int j = jL0_; j <= jL1_; j++){                                  //part IIIa: sumXU * sumYL
                volTypeSection(i,j,k_,nnxtULL_,n2ndULL_,diULL_,djULL_,dkULL_); 
            }
            for(int j = jM0_; j <= jM1_; j++){                                  //part IIIb: sumXU * sumYM
                volTypeSection(i,j,k_,nnxtUML_,n2ndUML_,diUML_,djUML_,dkUML_); 
            }
            for(int j = jU0_; j <= jU1_; j++){                                  //part IIIc: sumXU * sumYU
                volTypeSection(i,j,k_,nnxtUUL_,n2ndUUL_,diUUL_,djUUL_,dkUUL_); 
            }
        }
    }
    else if (k_ < kU0_) {                                                       //middle section of box in z-direction (do not need kM1_)
        for (int i = iL0_; i <= iL1_; i++){                                     //part I: sumXL * (sumYL + sumYM + sumYU    
            for (int j = jL0_; j <= jL1_; j++){                                 //part Ia: sumXL * sumYL
                volTypeSection(i,j,k_,nnxtLLM_,n2ndLLM_,diLLM_,djLLM_,dkLLM_);  
            }
            for (int j = jM0_; j <= jM1_; j++){                                 //part Ib: sumXL * sumYM
                volTypeSection(i,j,k_,nnxtLMM_,n2ndLMM_,diLMM_,djLMM_,dkLMM_); 
            }
            for (int j = jU0_; j <= jU1_; j++){                                 //part Ic: sumXL * sumYU
                volTypeSection(i,j,k_,nnxtLUM_,n2ndLUM_,diLUM_,djLUM_,dkLUM_); 
            }
        }
        for (int i = iM0_; i <= iM1_; i++){                                     //part II: sumXM * (sumYL + sumYM + sumYU)
            for (int j = jL0_; j <= jL1_; j++){                                 //part IIa: sumXM * sumYL
                volTypeSection(i,j,k_,nnxtMLM_,n2ndMLM_,diMLM_,djMLM_,dkMLM_); 
            }
            for (int j = jM0_; j <= jM1_; j++){                                 //part IIb: sumXM * sumYM
                volTypeSection(i,j,k_,nnxt_,n2nd_,di_,dj_,dk_);                 //no need for a MMM case
            }
            for (int j = jU0_; j <= jU1_; j++){                                 //part IIc: sumXM * sumYU
                volTypeSection(i,j,k_,nnxtMUM_,n2ndMUM_,diMUM_,djMUM_,dkMUM_); 
            }
        }
        for (int i = iU0_; i <= iU1_; i++){                                     //part III: sumXU * (sumYL + sumYM + sumYU)    
            for (int j = jL0_; j <= jL1_; j++){                                 //part IIIa: sumXU * sumYL
                volTypeSection(i,j,k_,nnxtULM_,n2ndULM_,diULM_,djULM_,dkULM_); 
            }
            for (int j = jM0_; j <= jM1_; j++){                                 //part IIIb: sumXU * sumYM
                volTypeSection(i,j,k_,nnxtUMM_,n2ndUMM_,diUMM_,djUMM_,dkUMM_); 
            }
            for (int j = jU0_; j <= jU1_; j++){                                 //part IIIc: sumXU * sumYU
                volTypeSection(i,j,k_,nnxtUUM_,n2ndUUM_,diUUM_,djUUM_,dkUUM_); 
            }
        }
    }
    else {                                                                      //top section of box (do not need kU1_)
        for (int i = iL0_; i <= iL1_; i++){                                     //part I: sumXL * (sumYL + sumYM + sumYU)  
            for (int j = jL0_; j <= jL1_; j++){                                 //part Ia: sumXL * sumYL
                volTypeSection(i,j,k_,nnxtLLU_,n2ndLLU_,diLLU_,djLLU_,dkLLU_);  
            }
            for (int j = jM0_; j <= jM1_; j++){                                 //part Ib: sumXL * sumYM
                volTypeSection(i,j,k_,nnxtLMU_,n2ndLMU_,diLMU_,djLMU_,dkLMU_); 
            }
            for (int j = jU0_; j <= jU1_; j++){                                 //part Ic: sumXL * sumYU
                volTypeSection(i,j,k_,nnxtLUU_,n2ndLUU_,diLUU_,djLUU_,dkLUU_); 
            }
        }
        for (int i = iM0_; i <= iM1_; i++){                                     //part II: sumXM * (sumYL + sumYM + sumYU)
            for (int j = jL0_; j <= jL1_; j++){                                 //part IIa: sumXM * sumYL    
                volTypeSection(i,j,k_,nnxtMLU_,n2ndMLU_,diMLU_,djMLU_,dkMLU_); 
            }
            for (int j = jM0_; j <= jM1_; j++){                                 //part IIb: sumXM * sumYM
                volTypeSection(i,j,k_,nnxtMMU_,n2ndMMU_,diMMU_,djMMU_,dkMMU_); 
            }
            for (int j = jU0_; j <= jU1_; j++){                                 //part IIc: sumXM * sumYU
                volTypeSection(i,j,k_,nnxtMUU_,n2ndMUU_,diMUU_,djMUU_,dkMUU_); 
            }
        }
        for (int i = iU0_; i <= iU1_; i++){                                     //part III: sumXU * (sumYL + sumYM + sumYU)
            for (int j = jL0_; j <= jL1_; j++){                                 //part IIIa: sumXU * sumYL
                volTypeSection(i,j,k_,nnxtULU_,n2ndULU_,diULU_,djULU_,dkULU_); 
            }
            for (int j = jM0_; j <= jM1_; j++){                                 //part IIIb: sumXU * sumYM
                volTypeSection(i,j,k_,nnxtUMU_,n2ndUMU_,diUMU_,djUMU_,dkUMU_); 
            }
            for (int j = jU0_; j <= jU1_; j++){                                 //part IIIc: sumXU * sumYU
                volTypeSection(i,j,k_,nnxtUUU_,n2ndUUU_,diUUU_,djUUU_,dkUUU_); 
            }
        }
    }
}


void HashGridVolume::prepareInitialSlice() {
    tag_ += dTag_; 
    minTag_ += dTag_;
    refTag_ += dTag_;                                                                                
    k_ = 0;                                                                     //initial slice: Note that k=-1 is set by boundary conditions at bottom of box
    kLneg_ = kLayer_ - 2;                                                       //never used but accurate: The most right side of the layer with periodic boundary
    kLmid_ = kLayer_ - 1;                                                       //slice index in layer that corresponds to the current  position  (mid => 0)
    kLpos_ = 0;                                                                 //slice index in layer that corresponds to the forward  position  (pos => +)

                                                                                //determine volume types within the boarders of slice 1 = KL.
    if( flagOpt_ == 0 ) {
       for(int i = 1; i < ni; i++){
          for(int j = 1; j < nj; j++){
            volTypeSection(i, j, kLpos_, nnxt_, n2nd_, di_, dj_, dk_); 
          }
       }
    } else {                          
        CompleteSliceIn27Sections();      
    }
}

bool HashGridVolume::prepareNextSlice() {
    k_++;
    tag_ += dTag_; 
    minTag_ += dTag_;
    refTag_ += dTag_;
    kLneg_ = kLmid_;                                                            //middle slice shifts to the left
    kLmid_ = kLpos_;                                                            //right flanking slice shifts to the middle
    kLpos_ = kLayerFromK_[k_];                                                  //this direct look up is performing: kLpos_ = k_ % kLayer_;

                                                                                //kLpos_ defines new slice in periodic layer as the right flank
                                                                                //last plane is all void-bulk and the process should be terminated after 
                                                                                //  defining this plane
    
    if (k_ == nk) {
        for (int i = 1; i < ni; i++) {
            for (int j = 1; j < nj; j++) {
                gridLayer_[i][j][kLpos_] = 2 + refTag_;                         //vType = gridLayer_[i][j][kLpos_] - refTag   vType=2 => void-bulk
                atomLayer_[i][j][kLpos_] = -1;                                  //void-bulk => no atom is in range
                closeAtomLayer_[i][j][kLpos_] = -1;                             //void-bulk => no atom is in range
            }
        }
        return false;                                                           //process will be stopped, went far enough along the z-direction
    }

                                                                                //consider using function pointers: See  http://www.newty.de/fpt/fpt.html#chapter2
    if( flagOpt_ == 0 ) {                                                       //determine volume types within the boarders of slice form 1 to kL.                                                
       for(int i = 1; i < ni; i++){
          for(int j = 1; j < nj; j++){              
             volTypeSection(i, j, k_, nnxt_, n2nd_, di_, dj_, dk_); 
          }
       }
    } else {
       CompleteSliceIn27Sections();                                             //27 cases divided as: (3 by k_) each of which has 9 cases 
    }
    return true;                                                                //process will continue to propagate forward in the +z-direction
}


int HashGridVolume::getVolumeType(int i, int j) {          
    closeAtom = atomLayer_[i][j][kLmid_];                                     //for slice[][] cluster labeling always use current slice in layer
    veryCloseAtom = closeAtomLayer_[i][j][kLmid_];
    int type = vTypeFromTag_[ gridLayer_[i][j][kLmid_] ];
    if (closeAtom < 0) {
        if (type == 1) {
            cout << "shouldn't happen\n";
        }
    }
    
    return vTypeFromTag_[ gridLayer_[i][j][kLmid_] ];                           //vType = {-1, 0, 1, 2} 
}




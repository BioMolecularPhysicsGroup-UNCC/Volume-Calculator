/* 
 * File:   HoshenKopelman.cpp
 * Author: jenny
 * 
 * Created on November 18, 2016, 5:03 PM
 */

#include "HoshenKopelman.hpp"
#include "HashGridVolume.hpp"

HoshenKopelman::HoshenKopelman(HashGridVolume * hashGrid, string fileBase) : hashGrid_(hashGrid), fileBase_(fileBase){
    nAtoms_ = hashGrid->nAtoms;  
    grid_ = hashGrid->grid;
    probe_ = hashGrid->probe;

#ifdef VISUALIZE
    cavityVisual_.open((fileBase_ + "_cavityCoordinates.txt").c_str());
    mvVisual_.open((fileBase_ + "_mvCoordinates.txt").c_str());
#endif
    
    minVoidSize_ = ceil( 4.188790205* (probe_ * probe_ * probe_) / (grid_ * grid_ * grid_) );  
//    boundaryProbe_ = probe_ * 4;
    int nAtoms1 = nAtoms_ + 1;
    int mAtoms = nAtoms_ * 20;
    microvoidPerAtom = (int *) new int[nAtoms1]; 
    boundaryPerAtom  = (int *) new int[nAtoms1]; 
    cavityPerAtom    = (int *) new int[nAtoms1]; 
    proteinPerAtom   = (int *) new int[nAtoms1]; 
    cavityChain_ = (int *) new int[mAtoms]; 
    cavityLabel_ = (int *) new int[mAtoms]; 
    cavityCount_ = (int *) new int[mAtoms];       
    cavityIndex_ = 0;
    maxCluster = 1;
    
    
    boundaryPerAtomClose  = (int *) new int[nAtoms1]; 
    cavityPerAtomClose    = (int *) new int[nAtoms1]; 
    
    mvPerAtomClose    = (int *) new int[nAtoms1]; 
    cavityChainClose_ = (int *) new int[mAtoms]; 
    cavityLabelClose_ = (int *) new int[mAtoms]; 
    cavityCountClose_ = (int *) new int[mAtoms];       
    cavityIndexClose_ = 0;
    
    for (int i=1; i<=nAtoms_; i++) {                                            //atoms range from 1 to nAtoms_
        cavityPerAtom[i]    = 0;                                                //must be initialized now for start of linked list
        microvoidPerAtom[i] = 0;                                                //must be initialized now for counter
        boundaryPerAtom[i] = 0;       
        proteinPerAtom[i] = 0;
        
        cavityPerAtomClose[i]    = 0;     
        mvPerAtomClose[i]    = 0;     
        boundaryPerAtomClose[i]    = 0;     
    }
     for (int i=0; i<mAtoms; i++) {
        cavityChain_[i] = 0;
        cavityLabel_[i] = 0;
        cavityCount_[i] = 0;
        
        cavityChainClose_[i] = 0;
        cavityLabelClose_[i] = 0;
        cavityCountClose_[i] = 0;
    }
    nnHalf_ = 13;                                                               //see getNeighbor() inline function
    neighbor_ = (int *) new int[nnHalf_];    
}


HoshenKopelman::~HoshenKopelman() {
    delete [] clusterLabel_; 
    delete [] clusterSize;
     
    delete [] cavityPerAtom; 
    delete [] microvoidPerAtom; 
    delete [] boundaryPerAtom; 
    delete [] proteinPerAtom; 
    delete [] cavityChain_;  
    delete [] cavityLabel_; 
    delete [] cavityCount_;    
    
    delete [] cavityPerAtomClose; 
    delete [] mvPerAtomClose; 
    delete [] boundaryPerAtomClose; 

    delete [] neighbor_;       
}

int HoshenKopelman::uf_find(int x) {
    int y = x;
    while (clusterLabel_[y] != y){
        y = clusterLabel_[y];
    }
    while (clusterLabel_[x] != x) {
        int z = clusterLabel_[x];
        clusterLabel_[x] = y;
        x = z;
    }
    return y;
}

int HoshenKopelman::uf_union(int x, int y) {                                    // joins two clusters and returns lowest cluster label
    int x2 = uf_find(x);
    int y2 = uf_find(y);
    if( x2 < y2 ) {
       clusterLabel_[y2] = x2;
       return x2;
    } else {                                                                    //this function is not called when x = y
       clusterLabel_[x2] = y2;
       return y2;
    }
}

int HoshenKopelman::uf_make_set(int type) {                                     //create new cluster and return the new cluster label
    maxCluster++;                                                               //add a new cluster label
    assert(maxCluster < nLabels_);                                              //error check to verify not exceeding maximum cluster label
    clusterSize[maxCluster] = type;                                             //add first gridpoint into new cluster according to type
    return maxCluster;
}



void HoshenKopelman::addVoidToFile(float x, float y, float z, int voidCluster, bool cavity){
    if (cavity) {
        cavityVisual_ << voidCluster << " " << x << " " << y << " " << z <<  "\n";        
    } else {
        mvVisual_ << voidCluster << " " << x << " " << y << " " << z <<  "\n";  
    }    
}



void HoshenKopelman::clusterVoidFile(int falseP[], int p, int cavities[], int c, float largestCluster, bool cavity) {
    
    cavityVisual_.close();    
    mvVisual_.close();  
    ifstream fPDBin;
    ofstream fPDBout;
    ofstream fErrors;
    ofstream fPDB;
    ofstream fPDBLarge;
    ofstream fLarge;
    if (cavity) {
        fPDBin.open((fileBase_ + "_cavityCoordinates.txt").c_str());
        fPDBout.open((fileBase_ + "_cavityCoordinatesFinal.txt").c_str());
        fErrors.open((fileBase_ + "_cavityCoordinatesErrors.txt").c_str()); 
        fPDB.open((fileBase_ + "_cavityCoordinates.pdb").c_str());
    } else {
        fPDBin.open((fileBase_ + "_mvCoordinates.txt").c_str());
        fPDBout.open((fileBase_ + "_mvCoordinatesFinal.txt").c_str());
        fPDB.open((fileBase_ + "_mvCoordinates.pdb").c_str());
        fPDBLarge.open((fileBase_ + "_mvLargest.pdb").c_str());
        fLarge.open((fileBase_ + "_mvLargestRange.txt").c_str(), ios_base::app);
    }
        
    PDBAtom * voidAtom = new PDBAtom;
    int voidAtomNum = 0;
        
    string line;     
    int cluster;
    string clusterString;
    float x;
    float y;
    float z;
    float xmin = FLT_MAX;
    float xmax = -FLT_MAX;
    float ymin = FLT_MAX;
    float ymax = -FLT_MAX;
    float zmin = FLT_MAX;
    float zmax = -FLT_MAX;     
    float xminC = FLT_MAX;
    float xmaxC = -FLT_MAX;
    float yminC = FLT_MAX;
    float ymaxC = -FLT_MAX;
    float zminC = FLT_MAX;
    float zmaxC = -FLT_MAX;     
        
    
    while(getline(fPDBin, line)) {   
        std::stringstream   data(line);
        data >> cluster >> x >> y >> z;  
        int position = line.find(' ');
        clusterString = line.substr(0, position);
        bool foundFP = false;
        for (int count = 0; count < p; count++) {
            if (falseP[count] == cluster) {
                foundFP = true;
            }
        } 
        bool foundC = false;
        for (int count = 0; count < c; count++) {
            if (cavities[count] == cluster) {
                foundC = true;
            }
        }
        cluster = uf_find(cluster);
        string test = static_cast <ostringstream*> (&(ostringstream() << cluster))->str();
        if ((cluster != clusterLabel_[1]) && (cluster != 0)) {
            if (foundFP) {
                fErrors << line.replace(0, position, test) << "\n";
            } 
            if ((foundC) || (!cavity)) {
                if (x > xmax) xmax = x;
                if (x < xmin) xmin = x;
                if (y > ymax) ymax = y;
                if (y < ymin) ymin = y;
                if (z > zmax) zmax = z;
                if (z < zmin) zmin = z;
                fPDBout << line.replace(0, position, test) << "\n";  
                voidAtom->number = voidAtomNum++;
                char chain[2] = "A";
                voidAtom->chain_id = chain[0];
                strncpy(voidAtom->atom_name, "O  ", 3);     
//                strncpy(voidAtom->res_sname, "HOH", 3);     
                voidAtom->res_num = cluster; 
                float hash = hashGrid_->grid;
                voidAtom->x = x * hash;// - hashGrid_->xShift;
                voidAtom->y = y * hash;// - hashGrid_->yShift;
                voidAtom->z = z * hash;// - hashGrid_->zShift;
                voidAtom->write(fPDB, false);
                fPDB << "\n";
                if ((cluster == largestCluster) && (!cavity)) {
                    if (x > xmaxC) xmaxC = x;
                    if (x < xminC) xminC = x;
                    if (y > ymaxC) ymaxC = y;
                    if (y < yminC) yminC = y;
                    if (z > zmaxC) zmaxC = z;
                    if (z < zminC) zminC = z;
                    voidAtom->write(fPDBLarge, false);
                    fPDBLarge << "\n";
                }
            }
        }
    }    
    
    percolated = false;
    
    if ((xmax <= xmaxC) && (xmin >= xminC)) percolated = true;
    if ((ymax <= ymaxC) && (ymin >= yminC)) percolated = true;
    if ((zmax <= zmaxC) && (zmin >= zminC)) percolated = true;
    if (!cavity) {
        int xRange = xmaxC - xminC;
        int yRange = ymaxC - yminC;
        int zRange = zmaxC - zminC;        
        fLarge << xRange << " " << yRange << " " << zRange << "\n";
    }
    fPDBin.close();
    fPDBout.close();
    fPDB.close();
    if (cavity) {
        fErrors.close();
    } else {
        fLarge.close();
        fPDBLarge.close();
    }
    delete voidAtom;
}






void HoshenKopelman::process3DGrid() { 
    
// slice1[][] , slice2[][]  two parallel planes storing cluster labels that toggle successively in the forward direction to scan a 3D box. 
// slice1[][] , slice2[][]  = 0 => cluster label for the direct excluded volume of a molecule.  
// slice1[][] , slice2[][]  = 1 => cluster label for any solvent accessible gridpoint that is also within the boundary volume.
// slice1[][] , slice2[][]  > 1 => cluster label for microvoids or cavities that are separated from bulk solvent.   
    
    
    ni_ = hashGrid_->ni;                                                        //length of fine-grid along x-axis normal to the y-z plane. 
    nj_ = hashGrid_->nj;                                                        //length of fine-grid along y-axis normal to the z-x plane. 
    nk_ = hashGrid_->nk;                                                        //length of fine-grid along z-axis normal to the x-y plane. 

// ---------------------------------------------- create heap memory and initialize as appropriate
    int** slice1 = (int **) new int* [ni_+1]; 
    int** slice2 = (int **) new int* [ni_+1];
        for(int i=0; i<=ni_; i++){
        slice1[i] = (int *) new int [nj_+1];
        slice2[i] = (int *) new int [nj_+1];
            for(int j=0; j<=nj_; j++){         // one time initialization is performed at the time of creation of these arrays
            slice1[i][j] = 1;                  // bottom of box is known to be a slice of bulk solvent (outside of boundary volume)
            slice2[i][j] = 1;                  // 1 is required on perimeter, but within the slice parimeter any number > 0 will work.
            }                                  // Since slice2[i][j] < 0 => microvoid, slice2[i][j] must have an initial > 0 value.
        }
    int nBulkGridPoints = (ni_+1)*(nj_+1);     // IMPORTANT: Assumes bottom face of box is bulk solvent, outside of boundary volume
    voidBulkPerimeter_ = 2*(ni_ + nj_);        // # of bulk-void gridpoints residing on the perimeter of each slice added 
//    cout << "voidBulkPerimeter_ = " << voidBulkPerimeter_  << "\n";

    nLabels_ = (ni_ * nj_ * nk_) / 50;         // overestimate of number of cluster labels needed.  FIX ME: Determine a reduced/better/definitive range
    clusterLabel_ = (int *) new int[nLabels_]; // link list that holds which sets of labels are part of the same cluster
    clusterSize   = (int *) new int[nLabels_]; // collects total amount of gridpoints within in each cluster
    for (int i=0; i < nLabels_; i++){       // initialize cluster size counters and linked list for cluster labels
        clusterSize[i] = 0;                    // all gridpoints are counted and will be put into some cluster
        clusterLabel_[i] = i;                  // self-referencing => no clustering
    }                                      // Note:  void-bulk gridpoints are likely to represent bulk solvent. 
    clusterSize[1] = nBulkGridPoints;          // record initial size of cluster 1: bottom face of box (slice1) + perimeter of slice2
                                               // Notes: 1) The first slice (bottom of box) is (and must be) 100% bulk solvent.
                                               //        2) The perimeter for all other slices must also be 100% bulk solvent
// ---------------------------------------------- perform alternating-propagation starting with slice1 known and slice2 all > values. 
        
    
   
    hashGrid_->prepareInitialSlice();          //sets vType per gridpoint in slice 2 & 1 corresponding to k= -1 & 0 respectively
    bool processON = true;                     //using implicit counters; k increments are used in comments for understanding
    int k = 1;
    while (processON) {
        hashGrid_->prepareNextSlice();           //vType is known for slice1 at k, will find vType for slice2 at k+1
        processSlice(slice1, slice2, k++);            //cluster labels for slice1 known ---> will find cluster labels for slice2 
        processON = hashGrid_->prepareNextSlice(); //vType is known for slice2 at k, will find vType for slice1 at k+1
        processSlice(slice2, slice1, k++);             //cluster labels for slice2 known ---> will find cluster labels for slice1
    }

    
    
    
    
    
// ---------------------------------------------- delete heap memory related to slice information after scan is complete
    for(int i=0; i<=ni_; i++){
        delete [] slice1[i];
        delete [] slice2[i];
    }
    delete [] slice1; 
    delete [] slice2;
// ---------------------------------------------- cluster the information stored in linked lists and cluster size arrays
    cluster();
}



void HoshenKopelman::processSlice(int** priorSlice, int** slice, int k) {
    int hasJoined;                                                              //flag used to determine what to do for union-find each time through
    int closestAtom;                                                            //atom label for the atom closest to the grid point being considered
    int clabel;                                                                 //cluster label for either microvoid, cavity or boundary volume

                                                                                //initialize appropriate counts for the preset perimeter of a slice of box
    clusterSize[1] += voidBulkPerimeter_;                                       //accounts for the perimeter gridpoints in the current slice added to cluster 1
    int pVolume = clusterSize[0];                                               //replace clusterSize[0] to count protein volume gridpoints with pVolume
    clusterSize[0] = 0;                                                         //set clusterSize[0]=0 to serve as a CRITICAL flag to skip neighboring gridpoint

                                                                                //check every gridpoint WITHIN slice perimeter, but NOT the perimeter.
    int vType;                                                                  //specifies the type of gridpoint under consideration. 
    for(int i=1; i<ni_; i++){                                                   //scan along the x-direction WITHIN perimeter of slice. Skip i=0 & i=ni_.
        for(int j=1; j<nj_; j++){                                               //scan along the y-direction WITHIN perimeter of slice. Skip j=0 & j=nj_.            
            vType = hashGrid_->getVolumeType(i, j);                             //vType = -1, 0, 1, 2 for microvoid, protein-volume, void-not-bulk, void-bulk
            if (vType == 0) {                                                   //protein volume
                slice[i][j] = 0;                                                //assign to gridpoint its cluster label 
                closestAtom = hashGrid_->closeAtom;
                if (closestAtom < 0) {
                    cout << "atom=-1, protein\n";
                }
                proteinPerAtom[closestAtom]++;
                pVolume++;                                                      //count each gridpoint that falls within the excluded volume of atoms
                                                                                //cannot have slice[i][j] = 0 otherwise voids connect with protein volume
            } else {                                                            //a microvoid or void
                getNeighbors(priorSlice, slice, i, j);
                                                                                //check all neighbors and add gridpoint to an appropiate cluster 
                hasJoined = -1;                                                 //flag = 1 => gridpoint has been added into a cluster
                if (vType < 0) {                                                //gridpoint is a microvoid
                    for(int n=0; n < nnHalf_; n++){                             //test all neighbors to gridpoint                       
                        clabel = neighbor_[n];                                  //cluster label of the n-th neigbhor to gridpoint[i][j];                          
                        if (clusterSize[clabel] < 0){                           //neighbor gridpoint is a microvoid when clusterSize < 0
                            if (hasJoined < 0){                                 //gridpoint has not been assigned to a microvoid cluster shared by a neighbor
                                hasJoined = clabel;                             //flag the fact that gridpoint has joined a cluster by recording its label 
                                clusterSize[hasJoined]--;                       //add the microvoid gridpoint into the joined cluster counter (negative)
                            } else{                                               
                                hasJoined = uf_union(hasJoined, clabel);    //join two clusters with hasJoined updated to lowest cluster label                        //gridpoint already belongs to a cluster
                                if (hasJoined != clabel){                       //apply union find if cluster labels are different                                                                             
                                }                   
                            }
                        }                         
                    }                                                           //continue loop over all neighbors
                    closestAtom = hashGrid_->closeAtom;             
                    microvoidPerAtom[closestAtom]++;    
                    if (hashGrid_->veryCloseAtom > 0) {
                        mvPerAtomClose[hashGrid_->veryCloseAtom]++;//(hasJoined, hashGrid_->veryCloseAtom);
                    }
                    if (hasJoined < 0) {                                        //gridpoint will be assigned a new cluster label
                        hasJoined = uf_make_set(-1);                            //-1 => microvoid                      
                    } 
#ifdef VISUALIZE
                    addVoidToFile(i, j, k, hasJoined, false);                      
#endif
                } else {                                                        //gridpoint is a void
                    for (int n=0; n < nnHalf_; n++){                             //test all neighbors to gridpoint
                        clabel = neighbor_[n];                                  //cluster label of the n-th neigbhor to gridpoint[i][j];  
                        if (clusterSize[clabel] > 0){                           //neighbor gridpoint is void-not-bulk or void-bulk when clusterSize > 0 
                            if(hasJoined < 0){                                  //gridpoint has not been assigned to a void cluster shared by a neighbor
                                hasJoined = clabel;                             //flag the fact that gridpoint has joined a cluster by recording its label
                                clusterSize[hasJoined]++;                       //add the void gridpoint into the joined cluster counter (positive)
                            } else {                                            //gridpoint already belongs to a cluster
                                if (hasJoined != clabel){                        //apply union find if cluster labels are different
                                    hasJoined = uf_union(hasJoined,clabel);     //join two clusters with hasJoined updated to lowest cluster label
                                }                                        
                            }    
                        }                         
                    }                                                           //continue loop over all neighbors

                    if(hasJoined < 0) {                                         //gridpoint will be assigned a new cluster label
                        hasJoined = uf_make_set(1);                             //1 => void (cavity)                        
                    } 
                    if (vType == 1) {                                           //gridpoint is void-not-bulk: record partial volume information                        
                        closestAtom = hashGrid_->closeAtom;                
                        clusterFraction(hasJoined, closestAtom);
                        if (hashGrid_->veryCloseAtom > 0) {
                            clusterFractionClose(hasJoined, hashGrid_->veryCloseAtom);
                        }
#ifdef VISUALIZE
                        if (hasJoined > 1) {
                            addVoidToFile(i, j, k, hasJoined, true);
                        }
#endif
                    }
                }
                slice[i][j] = hasJoined;                                        // assign to gridpoint its cluster label 
            }                                                                   // finished working if-condition for known vType           
        }                                                                       // continue with loop over index j for slice 
    }                                                                           // continue with loop over index i for slice 
    clusterSize[0] = pVolume;                                                   // reset clusterSize[0] to size of cluster 0 as it should be 
}


inline void HoshenKopelman::getNeighbors(int ** priorSlice, int ** slice, int i, int j) {

// --------------------------------------------- populated nearest-neighbors (3 for cubic lattice)  
     neighbor_[0] = priorSlice[i][j];          //the distance between nn is grid_ 
     neighbor_[1] = slice[i][j-1];                       
     neighbor_[2] = slice[i-1][j]; 
// --------------------------------------------- populated next-nearest-neighbors (6 for cubic lattice)
     neighbor_[3] = slice[i-1][j-1];           //the distance between nnn is sqrt(2)*grid_
     neighbor_[4] = slice[i-1][j+1];
     neighbor_[5] = priorSlice[i-1][j];
     neighbor_[6] = priorSlice[i+1][j];
     neighbor_[7] = priorSlice[i][j-1];
     neighbor_[8] = priorSlice[i][j+1];
// --------------------------------------------- populated next-next-nearest-neighbors (4 for cubic lattice)
     neighbor_[9] = priorSlice[i-1][j-1];      //the distance between nnnn is sqrt(3)*grid_
     neighbor_[10]= priorSlice[i-1][j+1];
     neighbor_[11]= priorSlice[i+1][j-1];
     neighbor_[12]= priorSlice[i+1][j+1];
}


void HoshenKopelman::clusterFraction(int label, int atom) {   
    
    if (cavityPerAtom[atom] == 0) {  
        cavityIndex_++;
        cavityPerAtom[atom] = cavityIndex_;
        cavityLabel_[cavityIndex_] = label;
        cavityCount_[cavityIndex_]++;
    } else {                                                                    // must search for current cavity-labels
        int j = cavityPerAtom[atom];
        while (j > 0) {
            if (label == cavityLabel_[j]) {                                     // label not new
                cavityCount_[j]++;
                break;                          
            }
            j = cavityChain_[j];
        }
        if (j == 0) {                                                           // => label Lc is different than all members in the list
            cavityIndex_++;                                                     // extend the list for a new label
            cavityChain_[cavityIndex_] = cavityPerAtom[atom];                   // append to linked-list
            cavityPerAtom[atom] = cavityIndex_;                                 // start at most current cavity label in list
            cavityLabel_[cavityIndex_] = label;                                 // record the cluster label
            cavityCount_[cavityIndex_]++;                                       // count this new label
        }
    }      
}

void HoshenKopelman::clusterFractionClose(int label, int atom) {     
       
    if (cavityPerAtomClose[atom] == 0) {  
        cavityIndexClose_++;
        cavityPerAtomClose[atom] = cavityIndexClose_;
        cavityLabelClose_[cavityIndexClose_] = label;
        cavityCountClose_[cavityIndexClose_]++;
    } else {                                                                    // must search for current cavity-labels
        int j = cavityPerAtomClose[atom];
        while (j > 0) {
            if (label == cavityLabelClose_[j]) {                                // label not new
                cavityCountClose_[j]++;
                break;                          
            }
            j = cavityChainClose_[j];
        }
        if (j == 0) {                                                           //=> label Lc is different than all members in the list
            cavityIndexClose_++;                                                // extend the list for a new label
            cavityChainClose_[cavityIndexClose_] = cavityPerAtomClose[atom];    // append to linked-list
            cavityPerAtomClose[atom] = cavityIndexClose_;                       // start at most current cavity label in list
            cavityLabelClose_[cavityIndexClose_] = label;                       // record the cluster label
            cavityCountClose_[cavityIndexClose_]++;                             // count this new label
        }
    }      
}


void HoshenKopelman::cluster() {

// ----------------------------------------- consolidate all cluster labeling and cluster size information
  int j;
  int iBottom;                                                                  // lowest possible label for a cluster
  for(int i=2; i <= maxCluster; i++){                                           // work from lowest to greatest labels 
    j = clusterLabel_[i];
    if (i != j) {
        iBottom = clusterLabel_[j];
        clusterLabel_[i] = iBottom;                                             // two steps always gets to the bottom
                                                                                // count all gridpoints in cluster
        clusterSize[iBottom] += clusterSize[i]; 
        clusterSize[i] = 0; 
       }
    }
    for(int a=1; a<=nAtoms_; a++){
        boundaryPerAtom[a] = 0;                                                 // initialize boundaryPerAtom[]
        boundaryPerAtomClose[a] = 0;      
    }
    
    for (int a=1; a <= nAtoms_; a++) {
        int nGridpoints = 0;                                                    // temporary counter for # of cavity gridpoints
        int j = cavityPerAtom[a];
        while (j > 0) {  
            int label = cavityLabel_[j];                                        // recorded cluster label
            int bottomLabel = clusterLabel_[label];                             // due to consolidation of cluster labels due not need uf_find
            if (bottomLabel == 1) {                                             // belongs to boundary volume
               boundaryPerAtom[a] += cavityCount_[j];                           // add gridpoints to the boundary volume counter
            } else {                                                            // belongs to cavity
               nGridpoints += cavityCount_[j];                                  // add gridpoints to the cavity cluster counter
            }
            j = cavityChain_[j];
        }      
        cavityPerAtom[a] = nGridpoints;                                         // now using c_index as total count for cavity gridpoints.         
    }
  
    for (int a=1; a <= nAtoms_; a++) {
        int nGridpoints = 0;                                                    // temporary counter for # of cavity gridpoints
        int j = cavityPerAtomClose[a];
        while (j > 0) {  
            int label = cavityLabelClose_[j];                                   // recorded cluster label
            int bottomLabel = clusterLabel_[label];                             // due to consolidation of cluster labels due not need uf_find
            if (bottomLabel == 1) {                                             // belongs to boundary volume
               boundaryPerAtomClose[a] += cavityCountClose_[j];                 // add gridpoints to the boundary volume counter
            } else {                                                            // belongs to cavity
               nGridpoints += cavityCountClose_[j];                             // add gridpoints to the cavity cluster counter
            }
            j = cavityChainClose_[j];
        }      
        cavityPerAtomClose[a] = nGridpoints;                                    // now using c_index as total count for cavity gridpoints.         
    }

    proteinVolume = clusterSize[0];                                             // by definition
    int totalSolvent = clusterSize[1];                                          // by definition this is the percolating bulk-void cluster within the box. 
    int totalCavity = 0;
    int totalMicrovoid = 0;
    int totalBoundary = 0;    
    
    for(int atom=1; atom<=nAtoms_; atom++){
        totalCavity += cavityPerAtom[atom]; 
        totalMicrovoid += microvoidPerAtom[atom]; 
        totalBoundary += boundaryPerAtom[atom]; 
    }

    int bulkSolvent = totalSolvent - totalBoundary;    
    int expectedGridpoints = (ni_ + 1)*(nj_ + 1)*(nk_ + 1);                     //exact box size
    int totalGridpoints = totalSolvent + totalMicrovoid + proteinVolume + totalCavity;
    int gridpointDifference = totalGridpoints - expectedGridpoints;

    if (gridpointDifference != 0){
       cout << "     " << endl;
       cout << "ERROR DETECTED: Total number of gridpoints do not sum up to box size!" << endl;
       cout << "         BOX SIZE = " << expectedGridpoints << endl;
       cout << "total grid points = " << totalGridpoints << endl;
       cout << "       difference = " << gridpointDifference << endl;
       cout << "---------------------------------------------------------------------" << endl;
       cout << "     " << endl;
    }
    cout << " box size gridpoints = " << expectedGridpoints << "\n";
    cout << "bulkspace gridpoints = " << bulkSolvent << "\n";
    cout << " boundary gridpoints = " << totalBoundary << "\n";
    cout << "microvoid gridpoints = " << totalMicrovoid << "\n";
    cout << "  protein gridpoints = " << proteinVolume << "\n";
    cout << "   cavity gridpoints = " << totalCavity << "\n";
    cout << "------------------------------------------" << "\n";
    
    numberOfVoids = 0; 
    numberOfMicrovoids = 0;
    numOfErrorGdpts = 0;
    cavityVolume = 0;
    microvoidVolume = 0;
    largestVoid = 0;
    largestMicrovoid = 0;
    float largestLabel = 0;
   
#ifdef VISUALIZE
    int falsePositives[maxCluster];
    int nFalsePs = 0;
    int cavities[maxCluster];
    int nCavities = 0;
#endif
    
    maximumClusterLabel = maxCluster;
    int i = 2;     
    float sm = 0;
    while (i <= maxCluster){        
        if (clusterSize[i] > 0){
            if (clusterSize[i] >= minVoidSize_){                                // if the void is big enough (no false positive voids that are actually microvoids)
                cavityVolume += clusterSize[i];
                numberOfVoids++;
                if (clusterSize[i] > largestVoid){
                    largestVoid = clusterSize[i];
                }
#ifdef VISUALIZE
                cavities[nCavities++] = clusterLabel_[i];
#endif
            } else {                                                            // its a false positive
                numOfErrorGdpts = numOfErrorGdpts + clusterSize[i];             // add to number of error points
                clusterSize[i] = -1 * clusterSize[i];                           // put it into microvoids
                microvoidVolume -= clusterSize[i];                              // we add it to microvoid instead
                numberOfMicrovoids++;
#ifdef VISUALIZE
                falsePositives[nFalsePs++] = clusterLabel_[i];
#endif
            }
        }
        else if (clusterSize[i] < 0){
            microvoidVolume -= clusterSize[i]; //because negative
            numberOfMicrovoids++;
            if(clusterSize[i] < largestMicrovoid){
                largestMicrovoid = clusterSize[i];
                largestLabel = clusterLabel_[i];
            }
            sm += clusterSize[i] * clusterSize[i];
        }        
        i++;
    }
    
    rsm = (sm - largestMicrovoid * largestMicrovoid) / i;
    
   

#ifdef VISUALIZE
    clusterVoidFile(falsePositives, nFalsePs, cavities, nCavities, largestLabel, true);
    clusterVoidFile(falsePositives, nFalsePs, cavities, nCavities, largestLabel, false);
#endif
    
    largestMicrovoid = -1 * largestMicrovoid;    
    boundaryVolume = totalBoundary;
}

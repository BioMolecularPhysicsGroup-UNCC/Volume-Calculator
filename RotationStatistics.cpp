/* 
 * File:   RotationStatistics.cpp
 * Author: jenny
 * 
 * Created on November 26, 2016, 4:11 PM
 */

#include "RotationStatistics.hpp"
#include "AtomCoordinates.hpp"
#include "HoshenKopelman.hpp"
#include <iomanip> 


RotationStatistics::RotationStatistics(bool flagOpt, float grid, float probe, float shell, float cutoff) 
                                     : flagOpt_(flagOpt),
                                       grid_(grid),
                                       probe_(probe),
                                       shell_(shell),
                                       cutoff_(cutoff){
}


RotationStatistics::~RotationStatistics() {
    delete [] probeSizes_;
}

float erfInverse(float x){
   float tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0f : 1.0f;

   x = (1 - x)*(1 + x);        
   lnx = logf(x);

   tt1 = 2/(M_PI*0.147) + 0.5f * lnx;
   tt2 = 1/(0.147) * lnx;

   return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}

void RotationStatistics::calculateProbeSizes(int nRotations, int fixedProbe) {
   
    float sigma = 0;
    if (!fixedProbe) {
        sigma = 0.05 * probe_;
    }
    probeSizes_ = new float[nRotations];
    for (int i = 0; i < nRotations; i++) {
        float x = (1.0 * i + 1)/(nRotations + 1);
        probeSizes_[i] = probe_ + (sigma * sqrt(2)) * erfInverse(2 * x - 1);
    }
}

void RotationStatistics::rotateRun(int nRotations, string outputPath, string PDBFileName, int fixedProbe) {   
    
    
    calculateProbeSizes(nRotations, fixedProbe);
    
    /*********************************************************************************************************************************************** 
     * Get coordinates from the PDB file
     ***********************************************************************************************************************************************/
    
    AtomCoordinates * atoms = new AtomCoordinates(true);                        //true => use Bondi van der Waals radii
    atoms->getAtomCoordinates(PDBFileName);                                     //read from PDB file all atomic coordinates  
    
      
    float startProbe = probe_;
    float endProbe = probe_;
    if (fixedProbe == 2) {
        startProbe = 0.3;
        endProbe = 2.0;       
    }
    float r = startProbe;
    while (r <= endProbe) {
    
        probe_ = r;   
    
    
    /*********************************************************************************************************************************************** 
     * Variable initialization
     ***********************************************************************************************************************************************/
    
    float proteinVolume[nRotations];                                            // Average values per rotation
    float maxClusterLabel[nRotations]; 
    float cavityVolume[nRotations];
    float microvoidVolume[nRotations];
    float boundaryVolume[nRotations];
    float packingDensity[nRotations];
    float cpuTime[nRotations];
        
    float proteinVolumeSD = 0;                                                  // Standard deviation amongst rotations
    float maxClusterLabelSD = 0;    
    float cavityVolumeSD = 0;
    float microvoidVolumeSD = 0;
    float boundaryVolumeSD = 0;
    float packingDensitySD = 0;
    float cpuTimeSD = 0;
        
    float proteinVolumeTotal = 0;                                               // Averages for all rotations
    float maxClusterLabelTotal = 0; 
    float cavityVolumeTotal = 0;
    float microvoidVolumeTotal = 0;
    float boundaryVolumeTotal = 0;
    float packingDensityTotal = 0;
    float cpuTimeTotal = 0;
                     
    vector <int> voidHist;                                                      // Vectors to collect rotational statistics for clusters
    vector <int> microvoidHist;  
    
        
        
        
        
        
    
    /*********************************************************************************************************************************************** 
     * Output file initialization
     ***********************************************************************************************************************************************/
    
    ofstream voidHistFile;                                                     
    ofstream microvoidHistFile;                                             
    ofstream microvoidLargestClusterFile;
    ofstream summaryFile;
    ofstream rsmFile;
    
    ostringstream gridString; 
    gridString << fixed << setprecision(2) << grid_; 
    ostringstream probeString; 
    probeString << fixed << setprecision(3) << probe_;
    
    string fileBase = PDBFileName + "-"  + gridString.str() + "_" + probeString.str();    
    string fileBaseOut = outputPath + fileBase;

    summaryFile.open((fileBaseOut + "_summary.txt").c_str());
    summaryFile.is_open();
    voidHistFile.open((fileBaseOut  + "_voidHist.txt").c_str());
    microvoidHistFile.open ((fileBaseOut + "_microvoidHist.txt").c_str());
    rsmFile.open ((fileBaseOut + "rsm.txt").c_str());
    
  
    
 
    /*********************************************************************************************************************************************** 
     * Main rotation loop.  For each rotation, call Hoshen-Kopelman algorithm, collect volumes and clusters, then rotate protein and repeat.
     ***********************************************************************************************************************************************/
    
    for(int rotation=0; rotation < nRotations; rotation++){ 
//#ifndef VISUALIZE
        atoms->rotateCoordinates();
//#endif
        float probeActual = 0;
        if (fixedProbe) {
            probeActual = probe_;
        } else {
            probeActual = probeSizes_[rotation];
        }        
        clock_t algorithmTime;                                                  // begin clock for individual rotation
        algorithmTime = clock();    
                                                                                // create hash table and run Hoshen-Kopelman algorithm
        HashGridVolume * hashGrid = new HashGridVolume(flagOpt_, grid_, probeActual, shell_,  cutoff_, atoms->atomCoords, atoms->nAtoms);
#ifdef VISUALIZE
        hashGrid->writeShiftedCoordinates(outputPath + PDBFileName);
#endif
        HoshenKopelman * hk = new HoshenKopelman(hashGrid, fileBase);
        hk->process3DGrid();
        rsmFile << hk->rsm <<  " " << hk->proteinVolume << " " <<hk->microvoidVolume << " " << hk->cavityVolume << " " << hk->boundaryVolume << "\n";
        if (hk->percolated) {
            ofstream percFile;
            percFile.open(("perc" + fileBaseOut + ".txt").c_str());
            percFile << "P\n";
            percFile.close();
            break;
        }
        algorithmTime = clock() - algorithmTime;                                // stop clock - only timing HK algorithm, not post-processing steps    
        
        if (voidHist.size() <= hk->largestVoid) {                               // appropriately size vectors for storing clusters 
            voidHist.resize(hk->largestVoid + 1, 0);
        }
        if (microvoidHist.size() <= hk->largestMicrovoid) {
            microvoidHist.resize(hk->largestMicrovoid + 1, 0);
        }        
        
        double maxSize = 0;
        
        for (int i = 2; i <= hk->maxCluster; i++){                               // count up clusters for microvoid and cavity
            if (hk->clusterSize[i] > 0){
                voidHist[hk->clusterSize[i]] ++;
            }
            if (hk->clusterSize[i] < 0){
                if ((-1 * hk->clusterSize[i]) > maxSize) {
                    maxSize = -1 * hk->clusterSize[i];
                }
                microvoidHist[-1 * hk->clusterSize[i]] ++;
            }
        }          
        
        microvoidLargestClusterFile.open ((fileBaseOut + "_largestClusterSize.txt").c_str(), ios_base::app);
        microvoidLargestClusterFile <<  maxSize << "\n";
        microvoidLargestClusterFile.close();
        
        proteinVolume[rotation] = hk->proteinVolume;                            // store results for each rotation
        maxClusterLabel[rotation] = hk->maximumClusterLabel;
        cavityVolume[rotation] = hk->cavityVolume;
        microvoidVolume[rotation] = hk->microvoidVolume;
        boundaryVolume[rotation] = hk->boundaryVolume;
        cpuTime[rotation] = ((float) algorithmTime)/CLOCKS_PER_SEC;
        float totalVolume = hk->proteinVolume + hk->microvoidVolume + hk->cavityVolume; 
        packingDensity[rotation] = hk->proteinVolume / totalVolume;      
       
        
        cout << "\n\n\n";
        cout << "PDB:                    " << PDBFileName << "\n";
        cout << "Probe Radius:           " << probe_ << "\n";
        cout << "Grid:                   " << grid_ << "\n";
        cout << "Maximum Cluster Label   " << hk->maximumClusterLabel << "\n";
        cout << "Number of Voids:        " << hk->numberOfVoids << "\n";
        cout << "Number of Microvoids:   " << hk->numberOfMicrovoids << "\n";
        cout << "Protein Volume:         " << hk->proteinVolume << " gridpoints,   " << setprecision(3) << hk->proteinVolume * pow(grid_, 3.0) << " angstroms cubed.\n";
        cout << "Boundary Water Volume:  " << hk->boundaryVolume << " gridpoints   " << setprecision(3) << hk->boundaryVolume * pow(grid_, 3.0) << " angstroms cubed.\n";
        cout << "Void volume:            " << hk->cavityVolume << " gridpoints     " << setprecision(3) << hk->cavityVolume * pow(grid_, 3.0) << " angstroms cubed.\n";
        cout << "Microvoid volume:       " << hk->microvoidVolume << " gridpoints   " << setprecision(3) << hk->microvoidVolume * pow(grid_, 3.0) << " angstroms cubed.\n";
        cout << "Error gridpoints:       " << hk->numOfErrorGdpts << endl;
        cout << "Packing Density:        " << setprecision(3) << hk->proteinVolume / totalVolume << "\n";      
        cout << "CPU time:               " << cpuTime[rotation] << "\n\n";
        
        atoms->residueVolumeSum(hk->cavityPerAtom, hk->microvoidPerAtom, hk->boundaryPerAtom, hk->proteinPerAtom, hk->cavityPerAtomClose, hk->mvPerAtomClose, hk->boundaryPerAtomClose);      

        if (rotation == 0) {
            ofstream rangeFile;
            rangeFile.open((outputPath + PDBFileName + "Ranges.txt").c_str());
            rangeFile << hashGrid->xRange << "\n";
            rangeFile << hashGrid->yRange << "\n";
            rangeFile << hashGrid->zRange << "\n";
            rangeFile.close();
        }
        
        delete hashGrid;
        delete hk;
    }
    
    ostringstream shellString; 
    shellString << fixed << setprecision(1) << shell_;  
    ostringstream cutoffString; 
    cutoffString << fixed << setprecision(3) << cutoff_; 
    
    string c = (cutoffString.str()).substr(0, 4);
    
//    string file1 = (PDBFileName + "-" + gridString.str() + "_" + probeString.str() + "_" + shellString.str()+ "_fractionalVolumes.txt");
//    string file2 = (PDBFileName + "-" + gridString.str() + "_" + probeString.str() + "_" + shellString.str()+ "_fractionalVolumesAtoms.txt");
//    string file3 = (PDBFileName + "-" + gridString.str() + "_" + probeString.str() + "_" + shellString.str()+ "_fractionalVolumesClose.txt");
//    string file4 = (PDBFileName + "-" + gridString.str() + "_" + probeString.str() + "_" + shellString.str()+ "_fractionalVolumesAtomsClose.txt");
    
    string file1 = (outputPath + PDBFileName + "-" + gridString.str() + "_" + probeString.str() + "_" + c + "_fractionalVolumes.txt");
    string file2 = (outputPath + PDBFileName + "-" + gridString.str() + "_" + probeString.str() + "_" + c + "_fractionalVolumesAtoms.txt");
    string file3 = (outputPath + PDBFileName + "-" + gridString.str() + "_" + probeString.str() + "_" + c + "_fractionalVolumesClose.txt");
    string file4 = (outputPath + PDBFileName + "-" + gridString.str() + "_" + probeString.str() + "_" + c + "_fractionalVolumesAtomsClose.txt");
    
    atoms->writeFractionalVolumes(nRotations, (file1).c_str(), (file2).c_str(), (file3).c_str(), (file4).c_str());
    
    
    /*********************************************************************************************************************************************** 
     * Write cavity and microvoid clusters to the histogram files     
     ***********************************************************************************************************************************************/
        
    for(int j = 1; j < voidHist.size(); j++){                                   
        if(voidHist[j] > 0){
//            voidHistFile << j * grid_ * grid_ * grid_ << " " << voidHist[j] << "\n";
            voidHistFile << j  << " " << voidHist[j] << "\n";
        }
    }
    for(int j = 1; j < microvoidHist.size(); j++){
        if(microvoidHist[j] > 0){
            microvoidHistFile << j << " " << microvoidHist[j] << "\n";
        }
    }
  
     
    
    /*********************************************************************************************************************************************** 
     * Calculate rotation averages and standard deviation, write all results to summary file 
     ***********************************************************************************************************************************************/
    
    for(int i=0; i<nRotations; i++){                                            //sum results for all rotations
        proteinVolumeTotal += proteinVolume[i];                                          
        maxClusterLabelTotal += maxClusterLabel[i]; 
        cavityVolumeTotal += cavityVolume[i];
        microvoidVolumeTotal += microvoidVolume[i];
        boundaryVolumeTotal += boundaryVolume[i];
        packingDensityTotal += packingDensity[i];
        cpuTimeTotal += cpuTime[i];
    }    
    
    float nRot = ( (float) nRotations*1.0 ); 
    proteinVolumeTotal   /= nRot;                                               //calculate averages for rotations
    maxClusterLabelTotal /= nRot;
    cavityVolumeTotal    /= nRot;    
    microvoidVolumeTotal /= nRot;    
    boundaryVolumeTotal  /= nRot;    
    packingDensityTotal  /= nRot;    
    cpuTimeTotal         /= nRot;      
    
    for(int i=0; i<nRotations; i++){                                
        proteinVolumeSD += (proteinVolume[i] - proteinVolumeTotal) * (proteinVolume[i] - proteinVolumeTotal);
        maxClusterLabelSD += (maxClusterLabel[i] - maxClusterLabelTotal) * (maxClusterLabel[i] - maxClusterLabelTotal);
        cavityVolumeSD += (cavityVolume[i] - cavityVolumeTotal) * (cavityVolume[i] - cavityVolumeTotal);
        microvoidVolumeSD += (microvoidVolume[i] - microvoidVolumeTotal) * (microvoidVolume[i] - microvoidVolumeTotal);
        boundaryVolumeSD += (boundaryVolume[i] - boundaryVolumeTotal) * (boundaryVolume[i] - boundaryVolumeTotal);
        packingDensitySD += (packingDensity[i] - packingDensityTotal) * (packingDensity[i] - packingDensityTotal);
        cpuTimeSD += (cpuTime[i] - cpuTimeTotal) * (cpuTime[i] - cpuTimeTotal);       
    }
    
    float correction_factor = (nRotations - 1) * nRotations;                    //corrected sample standard deviation   
       
    proteinVolumeSD = pow(proteinVolumeSD / correction_factor, 0.5); 
    maxClusterLabelSD = pow(maxClusterLabelSD / correction_factor, 0.5); 
    cavityVolumeSD = pow(cavityVolumeSD / correction_factor, 0.5);
    microvoidVolumeSD = pow(microvoidVolumeSD / correction_factor, 0.5);
    boundaryVolumeSD = pow(boundaryVolumeSD / correction_factor, 0.5);
    cpuTimeSD = pow(cpuTimeSD / correction_factor, 0.5);
    packingDensitySD = pow(packingDensitySD / correction_factor, 0.5);   
   
    
    float v = 1;//grid_*grid_*grid_;                                           //write summary results to file
    summaryFile << "Grid:               " << grid_ << "\n";
    summaryFile << "Probe:              " << probe_ << "\n";
    summaryFile << "Protein volume:     " << v*proteinVolumeTotal   << " " << v*proteinVolumeSD << "\n";
    summaryFile << "Cavity volume:      " << v*cavityVolumeTotal    << " " << v*cavityVolumeSD << "\n";
    summaryFile << "Microvoid volume:   " << v*microvoidVolumeTotal << " " << v*microvoidVolumeSD << "\n";
    summaryFile << "Boundary volume:    " << v*boundaryVolumeTotal  << " " << v*boundaryVolumeSD << "\n";
    summaryFile << "Packing density:    " << packingDensityTotal    << " " << packingDensitySD << "\n";
    summaryFile << "CPU time:           " << cpuTimeTotal           << " " << cpuTimeSD << "\n";   
    summaryFile << "Number of residues  " << atoms->nResidues << "\n";
    summaryFile << "Number of atoms     " << atoms->nAtoms << "\n";
    summaryFile << "max cluster label:  " << maxClusterLabelTotal   << " " << maxClusterLabelSD << "\n";   
    summaryFile << "Number of rotations:" << nRotations << "\n";

   
    summaryFile.close();
    microvoidHistFile.close();
    voidHistFile.close();   
    rsmFile.close();

    if (r < 0.49) {
        r += 0.01;
    } else {
        r += 0.1;
    }

    
    }
    delete atoms;
}

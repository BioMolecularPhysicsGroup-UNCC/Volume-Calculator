/* 
 * File:   main.cpp
 * Author: jenny
 * version: jacobs
 *
 * Created on November 13, 2016, 7:11 AM
 */

#include <cstdlib>

#include "AtomCoordinates.hpp"
#include "HashGridVolume.hpp"
#include "HoshenKopelman.hpp"
#include "RotationStatistics.hpp"

using namespace std;

int main(int argc, char** argv) {
    
    /*********************************************************************************************************************************************** 
     * Variable initialization
     ***********************************************************************************************************************************************/
    
    float grid = 0.5;                                                             // default values for grid, probe, PDB, rotations, shell
    float probe = 2.0;
    float shell = -4;                                                             // boundary volume length or default = -(# of probe radii forming a length)
    float cutoff = 0;

    string PDBFileName = "1bp2.pdb";    
    string outputPath = "";//results/";
    bool oneRun = false;                                                        // set to true for a single run, false for multiple rotations with output file summaries
    int  rotations = 1;
    bool flagOpt = false;                                                        // true => use optimized search algorithm, false => benchmark using unoptimized search method
    int  iq;                                                                    // dummy variable -- not sure how to work getopt() function
    int  fixedProbe = 1;
    
    
    /*********************************************************************************************************************************************** 
     * Process command line arguments 
     ***********************************************************************************************************************************************/
        
    int c;                                                                      
    while ((c = getopt (argc, argv, "n:g:p:r:s:h:o:x:f:b:")) != -1)
    switch (c){
        case 'n':                                                               // PDB file name
            PDBFileName = optarg;
            cout << "PDB: " << PDBFileName << "\n";
            break;
        case 'b':                                                               // PDB file name
            outputPath = optarg;
            cout << "Output Path " << outputPath << "\n";
            break;
        case 'g':                                                               // grid size - angstroms
            grid = atof(optarg);
            cout << "grid: " << grid << "\n";
            break;
        case 'p':                                                               // probe radius - angstroms
            probe = atof(optarg);            
            cout << "probe: " << probe << "\n";
            break;   
        case 'x':
            fixedProbe = atoi(optarg);
            if (fixedProbe) {
                cout << "Using fixed probe for rotations\n";
            }
            break;
        case 'r':                                                               // number of rotations. 
            rotations = atoi(optarg);            
            cout << "rotations: " << rotations << "\n";
            if (rotations == 1) {
//                oneRun = true;                                                  // need at least 2 rotations to get statistics
            }
            break;      
        case 's':                                                               // shell length - angstroms
            shell = atof(optarg);            
            cout << "shell: " << shell << "\n";
            break;    
        case 'f':                                                               // shell length - angstroms
//            cutoff = probe / atof(optarg);            
            cutoff = atof(optarg);            
            cout << "cutoff: " << cutoff << "\n";
            break;   
        case 'h':                                                               // command line help, format
            cout << "-n filename -g grid_spacing -p probe_radius -r rotation_number -s shell_length -o 0 \n";
            exit(-1);
        case 'o':                                                               // command line help, format
            iq = atoi(optarg);            
               if( iq == 0 ) {
               flagOpt = false; 
               } 
            break;   
        default:
            abort();
    }
    if (fabs(shell + 4) < 1.0e-5) {
        shell = 4 * probe;
    }

//    cutoff = probe/2;
//    shell = 2.8 * probe;
    if (cutoff == 0) {
        cutoff = probe;
    }
    cutoff = 3.5;
    shell = 3.5;
    cout << "cutoff = " << cutoff << endl;
    cout << "shell = " << shell << endl;
    
    if (grid > probe){
        cout << "         " << endl;
        cout << "probe = " << probe << endl;
        cout << " grid = " << grid << endl;
        cout << "ERROR: grid cannot be larger than probe" << endl;
        cout << "         " << endl;
        exit(-1);
    }

    if (grid < 0.04999999  ){
        cout << "         " << endl;
        cout << " grid = " << grid << endl;
        cout << "NOT ALLOWED: It is not a good idea to set grid less than 0.05 Angstroms!" << endl;
        cout << "             A very small value that is probably too small for large systems is 0.1 Angstroms." << endl;
        cout << "         " << endl;
        exit(-1);
    }

    /*********************************************************************************************************************************************** 
     * Calculate volumes 
     ***********************************************************************************************************************************************/
        
    if (oneRun) {    
        AtomCoordinates * atoms = new AtomCoordinates(true);                    //true => use Bondi van der Waals radii
        atoms->getAtomCoordinates(PDBFileName);                                 //read from PDB file all atomic coordinates
        atoms->rotateCoordinates();
        
// ------------------------------------------------------------------------------------------------------------------------ FIX ME
// NEEDS TO BE ADDED
// 1. Calculate center of mass coordinates and shift coordinates so that the center of mass is the new origin. 
// 2. With the new coordinates calculate the moment of inertia tensor 
// 3. Rotate the coordinates so that the principle axes align with the x, y and z axes. 
// 4. Swap the x, y and z axes in such a way that the x-y plane has the smallest extent. 
//    This is determined by making sure that after the rotation (xmax - xmin)*(ymax - ymin) gives the mimimum area
// 5. Once this rotation is complete, send this structure into the 1-shot run
// 6. We want to time this OPTIMAL 1-shot run with the average rotation times and this result is part of the publication. 
//    NOTE: Sheridan programmed at least steps 1, 2, 3. I am not sure if he optimized the swapping, but he might have.
// ------------------------------------------------------------------------------------------------------------------------

        clock_t algorithmTime;                                    // clock grid generation and HK algorithm only
        algorithmTime = clock();    
        HashGridVolume * hashGrid = new HashGridVolume(flagOpt, 
        grid, probe, shell, cutoff, atoms->atomCoords, atoms->nAtoms);
        ofstream rangeFile;
        rangeFile.open((PDBFileName + "Ranges.txt").c_str());
        rangeFile << hashGrid ->xRange << "\n";
        rangeFile << hashGrid ->yRange << "\n";
        rangeFile << hashGrid ->zRange << "\n";
        rangeFile.close();
        HoshenKopelman * hk = new HoshenKopelman(hashGrid, "");
        hk->process3DGrid();
        algorithmTime = clock() - algorithmTime;       
        
        float totalVolume = hk->proteinVolume + hk->microvoidVolume + hk->cavityVolume; 
               
        cout << "\n\n\n";
        cout << "PDB:                    " << PDBFileName << "\n";
        cout << "Probe Radius:           " << probe << "\n";
        cout << "Grid:                   " << grid << "\n";
        cout << "Maximum Cluster Label:  " << hk->maximumClusterLabel << "\n";
        cout << "Number of Voids:        " << hk->numberOfVoids << "\n";
        cout << "Number of Microvoids:   " << hk->numberOfMicrovoids << "\n";
        cout << "Protein Volume:         " << hk->proteinVolume << " gridpoints,   " << setprecision(3) << hk->proteinVolume * pow(grid, 3.0) << " angstroms cubed.\n";
        cout << "Boundary Water Volume:  " << hk->boundaryVolume << " gridpoints   " << setprecision(3) << hk->boundaryVolume * pow(grid, 3.0) << " angstroms cubed.\n";
        cout << "Void volume:            " << hk->cavityVolume << " gridpoints     " << setprecision(3) << hk->cavityVolume * pow(grid, 3.0) << " angstroms cubed.\n";
        cout << "Microvoid volume:       " << hk->microvoidVolume << " gridpoints   " << setprecision(3) << hk->microvoidVolume * pow(grid, 3.0) << " angstroms cubed.\n";
        cout << "Error gridpoints:       " << hk->numOfErrorGdpts << endl;
        cout << "Packing Density:        " << setprecision(3) << hk->proteinVolume / totalVolume << "\n\n";            
        cout << "\n\n";       
        
        delete hashGrid;                                                        // cleanup
        delete hk;        
        delete atoms;       
    } else {      
        RotationStatistics * rotation = new RotationStatistics(flagOpt, grid, probe, shell, cutoff);
	rotation->rotateRun(rotations, outputPath, PDBFileName, fixedProbe);
	delete rotation;
    }
    return 0;
}


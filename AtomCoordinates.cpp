/* 
 * File:   AtomCoordinates.cpp
 * Author: jenny
 * 
 * Created on November 18, 2016, 12:50 PM
 */

#include "AtomCoordinates.hpp"
#define SQR(x) (x*x)

AtomCoordinates::AtomCoordinates(bool bondi) : bondi_(bondi) {
}


AtomCoordinates::~AtomCoordinates() {
    

    for (int i = 0; i < nChains_; i++) {
        delete [] cavityPerResidue1_[i];
        delete [] microvoidPerResidue1_[i];
        delete [] boundaryPerResidue1_[i];
        delete [] proteinPerResidue1_[i];

        delete [] cavityPerResidue2_[i];
        delete [] microvoidPerResidue2_[i];
        delete [] boundaryPerResidue2_[i];
        delete [] proteinPerResidue2_[i];     
        
        delete [] cavityPerResidueClose1_[i];
        delete [] mvPerResidueClose1_[i];
        delete [] boundaryPerResidueClose1_[i];

        delete [] cavityPerResidueClose2_[i];
        delete [] mvPerResidueClose2_[i];
        delete [] boundaryPerResidueClose2_[i];
    }
    
    
    delete [] cavityPerResidue1_;
    delete [] microvoidPerResidue1_;
    delete [] boundaryPerResidue1_;
    delete [] proteinPerResidue1_;

    delete [] cavityPerResidue2_;
    delete [] microvoidPerResidue2_;
    delete [] boundaryPerResidue2_;
    delete [] proteinPerResidue2_;
    
    delete [] cavityPerAtom1_;
    delete [] microvoidPerAtom1_;
    delete [] boundaryPerAtom1_;
    delete [] proteinPerAtom1_;

    delete [] cavityPerAtom2_;
    delete [] microvoidPerAtom2_;
    delete [] boundaryPerAtom2_;
    delete [] proteinPerAtom2_;
    
    delete [] cavityPerResidueClose1_;
    delete [] mvPerResidueClose1_;
    delete [] boundaryPerResidueClose1_;

    delete [] cavityPerResidueClose2_;
    delete [] mvPerResidueClose2_;
    delete [] boundaryPerResidueClose2_;
    
    delete [] cavityPerAtomClose1_;
    delete [] boundaryPerAtomClose1_;

    delete [] cavityPerAtomClose2_;
    delete [] boundaryPerAtomClose2_;

    for (int i = 0; i <= nAtoms; i++) {
        delete [] atomCoords[i];
    }

    delete [] atomCoords;
    delete [] atomTypes;
    delete [] residueTypes;

}

void AtomCoordinates::getAtomCoordinates(string PDBFileName) {
    
    Logger * logger = LoggerFactory::default_logger();
    logger->set_log_level(Logger::CRITICAL);
    
    FASTInitializer fi;
    fi.load_user_defined_libraries();

    ifstream fin(PDBFileName.c_str());
    PDBProtein *protein = 0;
    protein = new PDBProtein;
    try {
        protein->read(PDBFileName.c_str());
    } catch (Exception e) {
        throw Exception(e.what());
    }
    fin.close();
    PDBProtein::AtomIterator iter;   
    
    nAtoms = protein->num_atoms();
    nResidues = protein->num_residues();
    
    atomCoords = (float**) new float* [nAtoms + 1]; 
    atomTypes = new const char* [nAtoms + 1]; 
    for(int i=0; i <= nAtoms; i++) {
        atomCoords[i] = (float*) new float [7]; 
        atomTypes[i] = (char*) new char [4];
    }   
   
        
    int count = 1;
    const PDBAtom * a;
    
    chains_ = new char [maxChainsAllowed_];
    maxChainResidues_ = new int[maxChainsAllowed_];
    nChains_ = 0;
    for (int i = 0; i < maxChainsAllowed_; i++) {
        maxChainResidues_[i] = INT_MIN;
    }
    
    for (iter = protein->atom_begin(); iter != protein->atom_end(); ++iter) {
        a = *iter;
        atomCoords[count][1] = a->x;
        atomCoords[count][2] = a->y;
        atomCoords[count][3] = a->z;
        atomCoords[count][4] = a->res_num;
        atomCoords[count][5] = a->chain_id;
        atomCoords[count][6] = PDBHelper::is_back_bone(a);
        atomTypes[count] = a->res_sname;
        bool newChain = true;
        for (int chain = 0; chain < nChains_; chain++) {
            if (a->chain_id == chains_[chain]) {
                if (a->res_num > maxChainResidues_[chain]) {
                    maxChainResidues_[chain] = a->res_num;
                }
                newChain = false;
                break;
            }
        }
        if (newChain) {
            chains_[nChains_++] = a->chain_id;            
        }
        
        const char symbol = a->atomic_symbol[1]; 
        string name;
        if (bondi_) {
            switch (symbol){
                case 'N':
                    atomCoords[count][0] = 1.55;//1.3;//
                    break;
                case 'H':
                    atomCoords[count][0] = 1.2;//1.1;//
                    hydrogenVDW_ = 1.2;// 1.1;//
                    break;
                case 'C':
                    name = a->atom_name;
                    if (name.compare(" CA ") == 0){
                        atomCoords[count][0] = 1.7;//1.3;//
                    } else {
                        atomCoords[count][0] = 1.7;//1.5;//
                    }
                    break;
                case 'O':
                    atomCoords[count][0] = 1.52;//1.4;//
                    break;
                case 'S':
                    atomCoords[count][0] = 1.8;//1.75;//
                    break;
                default:
                    atomCoords[count][0] = 0;
                    cout << "";                    
            }
            
        }             
        count++;
    }
         
  
    
    maxResiduesInChain_ = 0;
    for (int chain = 0; chain < nChains_; chain++) {
        if (maxChainResidues_[chain] > maxResiduesInChain_) {
            maxResiduesInChain_ = maxChainResidues_[chain];
        }
    }
    
    residueTypes = new const char* [maxResiduesInChain_ + 1]; 
    for(int i=0; i <= maxResiduesInChain_; i++) {
        residueTypes[i] = (char*) new char [4];
    }
    
    cavityPerAtom1_ = (float *) new float [nAtoms + 1];                
    microvoidPerAtom1_ = (float *) new float [nAtoms + 1];
    boundaryPerAtom1_ = (float *) new float [nAtoms + 1];
    proteinPerAtom1_ = (float *) new float [nAtoms + 1];

    cavityPerAtom2_ = (float *) new float [nAtoms + 1];
    microvoidPerAtom2_ = (float *) new float [nAtoms + 1];
    boundaryPerAtom2_ = (float *) new float [nAtoms + 1];
    proteinPerAtom2_ = (float *) new float [nAtoms + 1];
    
    cavityPerAtomClose1_ = (float *) new float [nAtoms + 1];   
    mvPerAtomClose1_ = (float *) new float [nAtoms + 1]; 
    boundaryPerAtomClose1_ = (float *) new float [nAtoms + 1];

    cavityPerAtomClose2_ = (float *) new float [nAtoms + 1];
    mvPerAtomClose2_ = (float *) new float [nAtoms + 1];
    boundaryPerAtomClose2_ = (float *) new float [nAtoms + 1];
    
    for(int atom = 1; atom <= nAtoms; atom++){
        cavityPerAtom1_[atom] = 0.0;              
        microvoidPerAtom1_[atom] = 0.0;
        boundaryPerAtom1_[atom] = 0.0;
        proteinPerAtom1_[atom] = 0.0;

        cavityPerAtom2_[atom] = 0.0; 
        microvoidPerAtom2_[atom] = 0.0; 
        boundaryPerAtom2_[atom] = 0.0; 
        proteinPerAtom2_[atom] = 0.0;         
        
        cavityPerAtomClose1_[atom] = 0.0;           
        mvPerAtomClose1_[atom] = 0.0;      
        boundaryPerAtomClose1_[atom] = 0.0;

        cavityPerAtomClose2_[atom] = 0.0; 
        mvPerAtomClose2_[atom] = 0.0; 
        boundaryPerAtomClose2_[atom] = 0.0; 
    }

    cavityPerResidue1_ = (float **) new float [nChains_];                
    microvoidPerResidue1_ = (float **) new float [nChains_];
    boundaryPerResidue1_ = (float **) new float [nChains_];
    proteinPerResidue1_ = (float **) new float [nChains_];

    cavityPerResidue2_ = (float **) new float [nChains_];
    microvoidPerResidue2_ = (float **) new float [nChains_];
    boundaryPerResidue2_ = (float **) new float [nChains_];
    proteinPerResidue2_ = (float **) new float [nChains_];
    
    cavityPerResidueClose1_ = (float **) new float [nChains_];    
    mvPerResidueClose1_ = (float **) new float [nChains_];   
    boundaryPerResidueClose1_ = (float **) new float [nChains_];

    cavityPerResidueClose2_ = (float **) new float [nChains_];
    mvPerResidueClose2_ = (float **) new float [nChains_];
    boundaryPerResidueClose2_ = (float **) new float [nChains_];
    
    for (int i = 0; i < nChains_; i++) {
        cavityPerResidue1_[i] = (float *) new float[maxResiduesInChain_ + 1];
        microvoidPerResidue1_[i] = (float *) new float[maxResiduesInChain_ + 1];
        boundaryPerResidue1_[i] = (float *) new float[maxResiduesInChain_ + 1];
        proteinPerResidue1_[i] = (float *) new float[maxResiduesInChain_ + 1];
        cavityPerResidue2_[i] = (float *) new float[maxResiduesInChain_ + 1];
        microvoidPerResidue2_[i] = (float *) new float[maxResiduesInChain_ + 1];
        boundaryPerResidue2_[i] = (float *) new float[maxResiduesInChain_ + 1];
        proteinPerResidue2_[i] = (float *) new float[maxResiduesInChain_ + 1];
        
        cavityPerResidueClose1_[i] = (float *) new float[maxResiduesInChain_ + 1];
        mvPerResidueClose1_[i] = (float *) new float[maxResiduesInChain_ + 1];
        boundaryPerResidueClose1_[i] = (float *) new float[maxResiduesInChain_ + 1];
        
        cavityPerResidueClose2_[i] = (float *) new float[maxResiduesInChain_ + 1];
        mvPerResidueClose2_[i] = (float *) new float[maxResiduesInChain_ + 1];
        boundaryPerResidueClose2_[i] = (float *) new float[maxResiduesInChain_ + 1];
    }
    
    for (int chain = 0; chain < nChains_; chain++) {
        for (int i = 0; i <= maxResiduesInChain_; i++) {
            cavityPerResidue1_[chain][i] = 0.0;
            microvoidPerResidue1_[chain][i] = 0.0;
            boundaryPerResidue1_[chain][i] = 0.0;
            proteinPerResidue1_[chain][i] = 0.0;

            cavityPerResidue2_[chain][i] = 0.0;
            microvoidPerResidue2_[chain][i] = 0.0;
            boundaryPerResidue2_[chain][i] = 0.0;
            proteinPerResidue2_[chain][i] = 0.0;
            
            cavityPerResidueClose1_[chain][i] = 0.0;
            mvPerResidueClose1_[chain][i] = 0.0;
            boundaryPerResidueClose1_[chain][i] = 0.0;

            cavityPerResidueClose2_[chain][i] = 0.0;
            mvPerResidueClose2_[chain][i] = 0.0;
            boundaryPerResidueClose2_[chain][i] = 0.0;
        }
    }
//    alignAxes();
}

void AtomCoordinates::residueVolumeSum(int * cavityPerAtom, int * microvoidPerAtom, int * boundaryPerAtom, int * proteinPerAtom, int * cavityPerAtomClose, int * mvPerAtomClose, int * boundaryPerAtomClose){
    
    float  ** temp_cpr = (float **) new float [nChains_];
    float  ** temp_mpr = (float **) new float [nChains_];
    float  ** temp_bpr = (float **) new float [nChains_];
    float  ** temp_ppr = (float **) new float [nChains_];
    
    float  ** temp_cprClose = (float **) new float [nChains_];
    float  ** temp_mprClose = (float **) new float [nChains_];
    float  ** temp_bprClose = (float **) new float [nChains_];
    
    for (int i = 0; i < nChains_; i++) {         
        temp_cpr[i] = (float *) new float[maxResiduesInChain_ + 1];
        temp_mpr[i] = (float *) new float[maxResiduesInChain_ + 1];
        temp_bpr[i] = (float *) new float[maxResiduesInChain_ + 1];
        temp_ppr[i] = (float *) new float[maxResiduesInChain_ + 1];        
        
        temp_cprClose[i] = (float *) new float[maxResiduesInChain_ + 1];
        temp_mprClose[i] = (float *) new float[maxResiduesInChain_ + 1];
        temp_bprClose[i] = (float *) new float[maxResiduesInChain_ + 1];
    }
        
    for (int chain = 0; chain < nChains_; chain++) {
        for (int i = 0; i <= maxResiduesInChain_; i++) {
            temp_cpr[chain][i] = 0.0;
            temp_mpr[chain][i] = 0.0;
            temp_bpr[chain][i] = 0.0;
            temp_ppr[chain][i] = 0.0;            
            
            temp_cprClose[chain][i] = 0.0;
            temp_mprClose[chain][i] = 0.0;
            temp_bprClose[chain][i] = 0.0;
        }
    }
    
    int previousResidue = -1;
    int resCount = 1;
    for(int atom = 1; atom <= nAtoms; atom++){
        int currentResidue = atomCoords[atom][4];                               //residue of the atom
        if (currentResidue != previousResidue) {
            residueTypes[resCount++] = atomTypes[atom];
        }
        previousResidue = currentResidue;
        int currentChain = atomCoords[atom][5];
        int chainID;
        
        for (int chain = 0; chain < nChains_; chain++) {
            if (currentChain == chains_[chain]) {
                chainID = chain;
                break;
            }
        }       
        
        cavityPerAtom1_[atom]    += cavityPerAtom[atom];
        microvoidPerAtom1_[atom] += microvoidPerAtom[atom];
        boundaryPerAtom1_[atom]  += boundaryPerAtom[atom];
        proteinPerAtom1_[atom]   += proteinPerAtom[atom];        
         
        cavityPerAtom2_[atom]    += cavityPerAtom[atom] * cavityPerAtom[atom];
        microvoidPerAtom2_[atom] += microvoidPerAtom[atom] * microvoidPerAtom[atom];
        boundaryPerAtom2_[atom]  += boundaryPerAtom[atom] * boundaryPerAtom[atom];
        proteinPerAtom2_[atom]   += proteinPerAtom[atom] * proteinPerAtom[atom];
        
        
        cavityPerAtomClose1_[atom]    += cavityPerAtomClose[atom];
        mvPerAtomClose1_[atom]        += mvPerAtomClose[atom];   
        boundaryPerAtomClose1_[atom]  += boundaryPerAtomClose[atom];   
         
        cavityPerAtomClose2_[atom]    += cavityPerAtomClose[atom] * cavityPerAtomClose[atom];
        mvPerAtomClose2_[atom]  += mvPerAtomClose[atom] * mvPerAtomClose[atom];
        boundaryPerAtomClose2_[atom]  += boundaryPerAtomClose[atom] * boundaryPerAtomClose[atom];
                                                                                //sum over first moments per residue
        temp_cpr[chainID][currentResidue] += cavityPerAtom[atom];
        temp_mpr[chainID][currentResidue] += microvoidPerAtom[atom];
        temp_bpr[chainID][currentResidue] += boundaryPerAtom[atom];       
        temp_ppr[chainID][currentResidue] += proteinPerAtom[atom];      
        
        temp_cprClose[chainID][currentResidue] += cavityPerAtomClose[atom];
        temp_mprClose[chainID][currentResidue] += mvPerAtomClose[atom];   
        temp_bprClose[chainID][currentResidue] += boundaryPerAtomClose[atom];      
        

    }

    for (int chain = 0; chain < nChains_; chain++) {
        for (int i = 0; i <= maxResiduesInChain_; i++) {
            cavityPerResidue1_[chain][i]    += temp_cpr[chain][i];;
            microvoidPerResidue1_[chain][i] += temp_mpr[chain][i];
            boundaryPerResidue1_[chain][i]  += temp_bpr[chain][i];
            proteinPerResidue1_[chain][i]  += temp_ppr[chain][i];
                                                                                //sum over second moments per residue
            cavityPerResidue2_[chain][i]    += temp_cpr[chain][i] * temp_cpr[chain][i];
            microvoidPerResidue2_[chain][i] += temp_mpr[chain][i] * temp_mpr[chain][i];
            boundaryPerResidue2_[chain][i]  += temp_bpr[chain][i] * temp_bpr[chain][i];
            proteinPerResidue2_[chain][i]   += temp_ppr[chain][i] * temp_ppr[chain][i];
            
            cavityPerResidueClose1_[chain][i]    += temp_cprClose[chain][i];
            mvPerResidueClose1_[chain][i]  += temp_mprClose[chain][i];
            boundaryPerResidueClose1_[chain][i]  += temp_bprClose[chain][i];
                                                                                //sum over second moments per residue
            cavityPerResidueClose2_[chain][i]    += temp_cprClose[chain][i] * temp_cprClose[chain][i];
            mvPerResidueClose2_[chain][i]  += temp_mprClose[chain][i] * temp_mprClose[chain][i];
            boundaryPerResidueClose2_[chain][i]  += temp_bprClose[chain][i] * temp_bprClose[chain][i];
            
            
        }
    }
}

void AtomCoordinates::writeFractionalVolumes(int nRotations, const char * fileNameResidue, const char * fileNameAtom, const char * fileNameResidueClose, const char * fileNameAtomClose) {
    
    ofstream fractionalVolumesFile;    
    fractionalVolumesFile.open(fileNameResidue);  
    ofstream fractionalVolumesFileClose ;    
    fractionalVolumesFileClose.open(fileNameResidueClose);
    
    float temp = nRotations * 1.0;
    for (int chain = 0; chain < nChains_; chain++) {
        for (int residue = 1; residue <= maxResiduesInChain_; residue++) {
            float y1ave = cavityPerResidue1_[chain][residue]/temp;                     //first moments
            float y2ave = microvoidPerResidue1_[chain][residue]/temp; 
            float y3ave = boundaryPerResidue1_[chain][residue]/temp; 
            float y4ave = proteinPerResidue1_[chain][residue]/temp; 

            float y1asq = cavityPerResidue2_[chain][residue]/temp;              //second moments
            float y2asq = microvoidPerResidue2_[chain][residue]/temp; 
            float y3asq = boundaryPerResidue2_[chain][residue]/temp;                
            float y4asq = proteinPerResidue2_[chain][residue]/temp; 
            
            float dy1 = sqrt (y1asq - y1ave*y1ave);                  //error bars (standard deviation)
            float dy2 = sqrt (y2asq - y2ave*y2ave); 
            float dy3 = sqrt (y3asq - y3ave*y3ave); 
            float dy4 = sqrt (y4asq - y4ave*y4ave); 

//            fractionalVolumesFile << residue << "  " << y1ave << "  " << y2ave << "  " << y3ave << "  " << y4ave 
//                                      << "       " << dy1 << "  " << dy2 << "  " << dy3 << "  " << dy4 << "\n";
            
                     
            double fractionBuried = y2ave / (y2ave + y3ave);
            double saveMV = y2ave;
            double saveCavity = y1ave;
            
            y1ave = cavityPerResidueClose1_[chain][residue]/temp;                     //first moments
            y2ave = mvPerResidueClose1_[chain][residue]/temp; 
            y3ave = boundaryPerResidueClose1_[chain][residue]/temp; 

            y1asq = cavityPerResidueClose2_[chain][residue]/temp;              //second moments
            y2asq = mvPerResidueClose2_[chain][residue]/temp;      
            y3asq = boundaryPerResidueClose2_[chain][residue]/temp;      
            
            dy1 = sqrt (y1asq - y1ave*y1ave);                  //error bars (standard deviation)
            dy2 = sqrt (y2asq - y2ave*y2ave); 
            dy3 = sqrt (y3asq - y3ave*y3ave); 

            fractionalVolumesFileClose << residue << "  " << y1ave << "  " << y2ave << "  " << y3ave << "  " << y4ave 
                                      << "       " << dy1 << "  " << dy2 << "  " << dy3 << "  " << dy4 << "\n";
            
            double localPacking = y4ave / (y4ave + saveMV + saveCavity + y3ave);            
            fractionalVolumesFile << residueTypes[residue] << " " << fractionBuried << " " << localPacking <<  "\n";
            
        }
    }
    
    fractionalVolumesFile.close();
    fractionalVolumesFileClose.close();
    
    
    ofstream fractionalAtomsFile;    
    fractionalAtomsFile.open(fileNameAtom);
    ofstream fractionalAtomsFileClose;    
    fractionalAtomsFileClose.open(fileNameAtomClose);
    
    for(int atom = 1; atom <= nAtoms; atom++){
        float y1ave = cavityPerAtom1_[atom]/temp;                     //first moments
        float y2ave = microvoidPerAtom1_[atom]/temp; 
        float y3ave = boundaryPerAtom1_[atom]/temp; 
        float y4ave = proteinPerAtom1_[atom]/temp; 

        float y1asq = cavityPerAtom2_[atom]/temp;                     //second moments
        float y2asq = microvoidPerAtom2_[atom]/temp; 
        float y3asq = boundaryPerAtom2_[atom]/temp;                
        float y4asq = proteinPerAtom2_[atom]/temp; 
            
        float dy1 = sqrt (y1asq - y1ave*y1ave);                      //error bars (standard deviation)
        float dy2 = sqrt (y2asq - y2ave*y2ave); 
        float dy3 = sqrt (y3asq - y3ave*y3ave); 
        float dy4 = sqrt (y4asq - y4ave*y4ave); 

        fractionalAtomsFile << atom << "  " << atomCoords[atom][4] << "  " << atomCoords[atom][6] << "  " << y1ave << "  " << y2ave << "  " << y3ave << "  " << y4ave 
                                      << "       " << dy1 << "  " << dy2 << "  " << dy3 << "  " << dy4 << "\n";
        
        
        y1ave = cavityPerAtomClose1_[atom]/temp;                     //first moments
        y2ave = mvPerAtomClose1_[atom]/temp; 
        y3ave = boundaryPerAtomClose1_[atom]/temp; 

        y1asq = cavityPerAtomClose2_[atom]/temp;                     //second moments
        y2asq = mvPerAtomClose2_[atom]/temp;        
        y3asq = boundaryPerAtomClose2_[atom]/temp;        
            
        dy1 = sqrt (y1asq - y1ave*y1ave);                      //error bars (standard deviation)
        dy3 = sqrt (y3asq - y3ave*y3ave); 

        fractionalAtomsFileClose << atom << "  " << atomCoords[atom][4] << "  " << atomCoords[atom][6] << "  " << y1ave << "  " << y2ave << "  " << y3ave << "  " << y4ave 
                                      << "       " << dy1 << "  " << dy2 << "  " << dy3 << "  " << dy4 << "\n";
    }
    fractionalAtomsFile.close();
    fractionalAtomsFileClose.close();
       
}


void AtomCoordinates::rotateCoordinates() {
    float angles[3];
    float tempx; 
    float tempy; 
    float tempz;
    srand(1234567);
//
//    srand(time(NULL));
    angles[0] = 2 * M_PI * ((float) rand()/RAND_MAX);                           //M_PI makes the program run slower but not noticeable here
    angles[1] = asin(((float) rand()/RAND_MAX));
    angles[2] = M_PI * (2*((float) rand()/RAND_MAX)-1);                         //M_PI makes the program run slower but not noticeable here

    for(int i=1; i <= nAtoms; i++){
        tempx = rotatex(angles[0], angles[1], angles[2], atomCoords[i][1], atomCoords[i][2], atomCoords[i][3]);
        tempy = rotatey(angles[0], angles[1], angles[2], atomCoords[i][1], atomCoords[i][2], atomCoords[i][3]);
        tempz = rotatez(angles[0], angles[1], angles[2], atomCoords[i][1], atomCoords[i][2], atomCoords[i][3]);
        atomCoords[i][1] = tempx; 
        atomCoords[i][2] = tempy; 
        atomCoords[i][3] = tempz;
    }
}

float AtomCoordinates::rotatex(float phi, float theta, float psi, float x, float y, float z){
    return cos(psi)*cos(theta)*x + (cos(psi)*sin(theta)*sin(phi)- cos(phi)*sin(psi))*y + (sin(psi)*sin(phi) + cos(psi)*cos(phi)*sin(theta))*z;
}

float AtomCoordinates::rotatey(float phi, float theta, float psi, float x, float y, float z){
    return cos(theta)*sin(psi)*x + (cos(psi)*cos(phi) + sin(psi)*sin(theta)*sin(phi))*y + (cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi))*z;
}

float AtomCoordinates::rotatez(float phi, float theta, float psi, float x, float y, float z){
    return -sin(theta)*x + cos(theta)*sin(phi)*y + cos(theta)*cos(phi)*z;
}



int AtomCoordinates::dsyevj3(float A[3][3], float Q[3][3], float w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the Jacobi algorithm.
// The upper triangular part of A is destroyed during the calculation,
// the diagonal elements are read but not destroyed, and the lower
// triangular elements are not referenced at all.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
{
  const int n = 3;
  float sd, so;                                                                 //Sums of diagonal resp. off-diagonal elements
  float s, c, t;                                                                //sin(phi), cos(phi), tan(phi) and temporary storage
  float g, h, z, theta;                                                         //More temporary storage
  float thresh;
  
                                                                                //Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
  for (int i=0; i < n; i++)
  {
    Q[i][i] = 1.0;
    for (int j=0; j < i; j++)
      Q[i][j] = Q[j][i] = 0.0;
  }
#endif

  // Initialize w to diag(A)
  for (int i=0; i < n; i++)
    w[i] = A[i][i];

  // Calculate SQR(tr(A))  
  sd = 0.0;
  for (int i=0; i < n; i++)
    sd += fabs(w[i]);
  sd = SQR(sd);
 
  // Main iteration loop
  for (int nIter=0; nIter < 50; nIter++) {
    // Test for convergence 
    so = 0.0;
    for (int p=0; p < n; p++)
      for (int q=p+1; q < n; q++)
        so += fabs(A[p][q]);
    if (so == 0.0)
      return 0;

    if (nIter < 4)
      thresh = 0.2 * so / SQR(n);
    else
      thresh = 0.0;

    // Do sweep
    for (int p=0; p < n; p++)
      for (int q=p+1; q < n; q++) {
        g = 100.0 * fabs(A[p][q]);
        if (nIter > 4  &&  fabs(w[p]) + g == fabs(w[p])
                       &&  fabs(w[q]) + g == fabs(w[q])) {
          A[p][q] = 0.0;
        }
        else if (fabs(A[p][q]) > thresh) {
          // Calculate Jacobi transformation
          h = w[q] - w[p];
          if (fabs(h) + g == fabs(h)) {
            t = A[p][q] / h;
          }
          else {
            theta = 0.5 * h / A[p][q];
            if (theta < 0.0)
              t = -1.0 / (sqrt(1.0 + SQR(theta)) - theta);
            else
              t = 1.0 / (sqrt(1.0 + SQR(theta)) + theta);
          }
          c = 1.0/sqrt(1.0 + SQR(t));
          s = t * c;
          z = t * A[p][q];

          // Apply Jacobi transformation
          A[p][q] = 0.0;
          w[p] -= z;
          w[q] += z;
          for (int r=0; r < p; r++)
          {
            t = A[r][p];
            A[r][p] = c*t - s*A[r][q];
            A[r][q] = s*t + c*A[r][q];
          }
          for (int r=p+1; r < q; r++)
          {
            t = A[p][r];
            A[p][r] = c*t - s*A[r][q];
            A[r][q] = s*t + c*A[r][q];
          }
          for (int r=q+1; r < n; r++)
          {
            t = A[p][r];
            A[p][r] = c*t - s*A[q][r];
            A[q][r] = s*t + c*A[q][r];
          }

          // Update eigenvectors
#ifndef EVALS_ONLY          
          for (int r=0; r < n; r++)
          {
            t = Q[r][p];
            Q[r][p] = c*t - s*Q[r][q];
            Q[r][q] = s*t + c*Q[r][q];
          }
#endif
        }
      }
  }

  return -1;
}

void AtomCoordinates::alignAxes(){
        
                                                                                //calculate center of mass
    float hydrogenCount = 0; float cm[4] = {0};                                 //cm[0] wont be used, stores locations of center of mass
    for(int i = 1; i <= nAtoms; i++){
        if(atomCoords[i][0] == hydrogenVDW_){                                   //if it is Carbon == carbonVDW
            hydrogenCount++;
            for(int j=1; j<=3; j++){
                cm[j] = cm[j] + atomCoords[i][j];
            }
        }
    }
    for(int j=1; j<=3; j++){                                                    //get final values for center of mass
                cm[j] = cm[j] / hydrogenCount;
            }
                                                                                //translate to center of mass (generalize translate function),
    for(int i = 1; i <= nAtoms; i++){
        for(int j=1; j<=3; j++){
            atomCoords[i][j] = atomCoords[i][j] - cm[j];                        //translate to center of mass
        }
    }
                                                                                //calculate moment of inertia tensor
    float tensor[3][3] = {0};                                                   //inertia tensor
                                                                                //tensor[0][0] is Ixx, tensor[1][1] is Iyy, tensor[2][2] is Izz
    for(int i = 1; i <= nAtoms; i++){
        if(atomCoords[i][0] == hydrogenVDW_){
            tensor[0][0] = tensor[0][0] + (atomCoords[i][2]*atomCoords[i][2])+(atomCoords[i][3]*atomCoords[i][3]); //Ixx
            tensor[1][1] = tensor[1][1] + (atomCoords[i][1]*atomCoords[i][1])+(atomCoords[i][3]*atomCoords[i][3]); //Iyy
            tensor[2][2] = tensor[2][2] + (atomCoords[i][1]*atomCoords[i][1])+(atomCoords[i][2]*atomCoords[i][2]); //Izz
            tensor[0][1] = tensor[0][1] - (atomCoords[i][1]*atomCoords[i][2]); //Ixy = Iyx
            tensor[0][2] = tensor[0][2] - (atomCoords[i][1]*atomCoords[i][3]); //Ixz = Izx
            tensor[1][2] = tensor[1][2] - (atomCoords[i][2]*atomCoords[i][3]); //Iyz = Izy
        }
    }
    tensor[1][0] = tensor[0][1]; tensor[2][0] = tensor[0][2]; tensor[2][1] = tensor[1][2]; //corresponding products of inertia

                                                                                //diagonalize symmetric tensor
    float evector[3][3];
    float evalue[3];
    dsyevj3(tensor, evector, evalue);                                           //external algorithm
    //insertion sort
    float replacement; float evReplacement[3];
    for(int i=1; i<3; i++){                                                     //put the eigenvalues and eigenvectors in order from smallest to largest
        int j = i;
        while(j > 0 && evalue[j-1] > evalue[j]){
            for(int k = 0; k < 3; k++){
                evReplacement[k] = evector[k][j];
                evector[k][j] = evector[k][j-1];
                evector[k][j-1] = evReplacement[k];
            }
            replacement = evalue[j];
            evalue[j] = evalue[j-1];
            evalue[j-1] = replacement;
            j = j - 1;
        }
    }
    //iterate over the three possibilities that determine the z direction and find which one yields fastest times. compare to other times
    //initialize new algorithm to pick smallest footprint
    float temp[3] = {0};
    for(int permute = 0; permute < 1; permute++){                               //three possibilities
        
        for(int i = 1; i <= nAtoms; i++){
            for(int j=0; j<3; j++){
                temp[j] = evector[permute % 3][j]*atomCoords[i][1] + evector[(permute+1) % 3][j]*atomCoords[i][2] + evector[(permute + 2) % 3][j]*atomCoords[i][3];
            }
            atomCoords[i][1] = temp[0]; 
            atomCoords[i][2] = temp[1]; 
            atomCoords[i][3] = temp[2]; 
        }               
    }
}
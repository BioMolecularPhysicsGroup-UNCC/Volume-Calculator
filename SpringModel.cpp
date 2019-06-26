/* 
 * File:   SpringModel.cpp
 * Author: jenny
 * 
 * Created on November 17, 2016, 8:26 PM
 */

#include "SpringModel.hpp"

SpringModel::SpringModel(int maxAtomsWithin, float grid, float probe, float ** atomCoords) 
                             : grid_(grid),
                               probe_(probe), 
                               atomCoords_(atomCoords) {
// --------------------------------------------------- define geometrical tolerances that are applied to the underlying discrete grid
  float ratio = (grid/0.5);                          //keep length scales proportional to grid, using grid = 0.5 as a reference.  
  float maxDrift = ratio*ref_maxDrift_;              //should be no more than ~2% of grid length. For grid size of 0.5 maxDrift = 0.01
  drMinMagnitude_ = maxDrift/100;                    //maximum spatial resolution just before accuracy is unwarrented.
  maxEnergyAllowed_ = ref_maxDrift_*ref_maxDrift_;   //cutoff energy for selecting between stressed and not stressed (INDEPENDENT of grid) 
  float dEzero = maxEnergyAllowed_/1000;             //set effective zero energy to 1000 times smaller than the cutoff energy
  maxE_ = maxEnergyAllowed_ - dEzero;                //define a maximum energy reference point used to calculate energy differences
// --------------------------------------------------- sanity check: Using single precision variables (float not double)
     if( drMinMagnitude_ < 0.000005 ){
     cout << "ERROR: requested precision requires double precision calculations " << endl;
     cout << "       re-Write the code with double .... OR " << endl;
     cout << "       Re-compile with the flag that converts all single precision to double .... OR " << endl;
     cout << "       Do not attempt to use this small grid size! The gain in accuracy is probably not worth the effort " << endl;
     cout << "TERMINATED " << endl;
     exit(-1);
     }


// --------------------------------------------------- setup optimized relaxation stages
  int nR = ((int) ceil(probe*shell_rate_) );         //number of intermediate spring relaxation steps based on expanding probe size
  float dVolume = (probe * probe * probe)/nR; 
  float initialProbeRadius = 1.2*grid; 
  float temp = (initialProbeRadius * initialProbeRadius * initialProbeRadius); 
  float initialVolume = dVolume/2; 
     if( temp > initialVolume ) {
     initialVolume = temp;
     }
  nR++;                                              //this extra 1 is for the probe radius = probe_
  radiusSchedule_ = (float *) new float[nR]; 
           maxAW_ = (int *) new int[nR];
  numRadii_ = 0;
     for(int n = 0; n < nR; n++){
     float sphereVolume = initialVolume + n*dVolume;   //without the 4*pi/3 (put that in print statement)
     float radius = cbrt(sphereVolume);
        if( radius < probe ){
        radiusSchedule_[numRadii_] = radius;
        numRadii_++;
        }
     }
  radiusSchedule_[numRadii_] = probe;
  numRadii_++;
// --------------------------------------------------- vdw-radius, coordinates and equilibrium-length for each compressive spring
  springR_ = (float *) new float[maxAtomsWithin];             
  springX_ = (float *) new float[maxAtomsWithin];             
  springY_ = (float *) new float[maxAtomsWithin];             
  springZ_ = (float *) new float[maxAtomsWithin];             
// *NO*  radius2_ = (float *) new float[maxAtomsWithin];    //squared distance between position of spring post and gridpoint of interest
  springL_ = (float *) new float[maxAtomsWithin];    //natural equilibrium length of spring under compression only

//if( numRadii_ > -10000 ) { exit(-1); }             //a temporary stop condition for checking program
// --------------------------------------------------- define various counting bins for performance diagnostics
/*  springIterations  = (int *) new int[maxAtomsWithin];
  springOperations  = (int *) new int[maxAtomsWithin];
  springNeighbors   = (int *) new int[maxAtomsWithin];
  testParticleOut   = (int *) new int[maxAtomsWithin];
  springsStressed   = (int *) new int[maxAtomsWithin];
  cavityIdentified  = (int *) new int[maxAtomsWithin];

     for(int a = 0; a < maxAtomsWithin; a++){
     springIterations[a] = 0;
     springOperations[a] = 0;
     springNeighbors[a] = 0;
     testParticleOut[a] = 0;
     springsStressed[a] = 0;
     cavityIdentified[a] = 0;
     }
 */
}

SpringModel::~SpringModel() {

//  delete [] springIterations;
//  delete [] springOperations;
//  delete [] springNeighbors;
//  delete [] testParticleOut;
//  delete [] springsStressed;
//  delete [] cavityIdentified;
  delete [] springR_;
  delete [] springX_;
  delete [] springY_;
  delete [] springZ_;
// *NO*  delete [] radius2_;
  delete [] springL_;
}

void SpringModel::setAtomCoordinates(float ** atomCoords){
    atomCoords_ = atomCoords;
}

int SpringModel::testSpring(float x, float y, float z, int * atomsWithin, int maxAtomsWithin) {
   
   /********************************************************************************************************************************************** 
    * Calculate the magnitude of force exerted by all nearby atoms on a test particle grid point location. 
    * These forces are from compressive springs only. As they compress, energy is created. 
    * Find the spring configuration that minimizes this energy by moving the test particle within a small region of connected space. 
    * The range of motion of the test particle is constrained within a spherical shell of radius probe from the (x,y,z)-gridpoint. 
    * The geometrical constraint that restricts the range of motion to be within a spherical shell is enforced using contractive springs. 
    * Compressive springs exert forces only when they are compressed, while contractive springs exert a force only when they are stretched. 
    * If a non-stressed test location can be found, the gridpoint located at (x,y,z) is part of a void space (cavity). 
    * Stressed spring configurations correspond to microvoids. 
    * ADVANCED ALGORITHM
    * A small test sphere, if found to be a microvoid implies any larger sphere will be a microvoid. 
    * To test for a spherical probe of radius Rprobe, test a series of successively larger probe sizes, R1, R2, R3, ... until Rprobe. 
    * It is much easier and faster to test for the microvoid condition of smaller spheres. 
    * Whenever a probe of smaller radius identifies a gridpoint to be microvoid, the actual probe sphere will find the same result.   
    * Important to this series of tests, note that to find a microvoid using small spheres is much much faster than a large sphere. 
    * Each time a smaller sphere is not found to be a microvoid, it implies it can fit INSIDE the space of the next larger probe size. 
    * Upon each failed attempt to identify microvoid, the test particle location of a small probe defines a NEW geometical constraint. 
    * The NEW (or augmented) geometrical constraint requires all larger probes not to drift beyond this test particle point.
    * As the process continues, more inner geometrical constraints are added as contractive forces, creating bounds on test particle position.
    * The bounds on the test particle position can be used to discard compressive spring constrants from distance neighbor atoms. 
    * This is because the combined set of geometrical constraints from inner contractive forces maintains the stress conditions. 
    * Inner contractive forces combine to create large stresses when the probe (of any size) drifts too far, making microvoid detection easier.
    * Easier microvoid detection translates into a reduction of required iterations to reach converged equilibrium conditions (stressed or not). 
    *
    * FIX ME: Consider precalcuating (probe+vdw)  and  (probe+shell) and storing into  atomsCoords_[atom][n]  for n= 5 and 6 respectively. 
    **********************************************************************************************************************************************/
//    springNeighbors[maxAtomsWithin] += 1;                           //count number of times this many nearest neighbors appear


// ------------------------------------------------------------------ prepare system for series of evaluations of forces as probe size increases
    int sMin = -1;                                                  //lowest probe size to consider (-1 => not known at this time)
    int m = 0;                                 
       for(int s = 0; s < numRadii_; s++) {                         //pre-order atoms to follow the order of probe size interactions 
       float pD = 2*radiusSchedule_[s];                             //set probe diameter
          for(int n = m; n < maxAtomsWithin; n++){                  //check all atoms within range of maximum size probe
          int atom_n = atomsWithin[n];  
          float vdw = atomCoords_[atom_n][0];                       //vdw-radius for the n-th atom
          float span = pD + vdw;                                    //greatest distance atom can be away before becoming irrelvant
          float s2 = span*span;  
          float dx = atomCoords_[atom_n][1] - x; 
          float dy = atomCoords_[atom_n][2] - y; 
          float dz = atomCoords_[atom_n][3] - z; 
          float r2 = dx*dx + dy*dy + dz*dz;                         //squared distance between gridpoint & n-th atom from the ORIGINAL list
             if( r2 < s2 ) {                                        //atom is relevant for this size probe (swap order of atom listing)
             int atom_m = atomsWithin[m];    
             atomsWithin[m] = atom_n; 
             atomsWithin[n] = atom_m; 
             m++;
             }
          }
          if( (sMin < 0) && (m > 0) ){ 
          sMin = s; 
          } 
       maxAW_[s] = m;
       }

       for(int n = 0; n < maxAtomsWithin; n++){                     //store pertinent data about compressive springs in convenient list
       springR_[n] = atomCoords_[atomsWithin[n]][0];                //vdw-radius for the n-th atom from the SORTED list

       }


// ------------------------------------------------------------------- initialize variables for entire calculation over series of probe sizes
    float test_x = x;
    float test_y = y;
    float test_z = z;

    float wCos;                                       //cosine of the angle between current and previous direction of motion of the test particle
    float r_mag;                               
    float energy;
    int   nSwitches;                                             //counts # of times test particle switches between any of the three categories
    float energyBottom;                                          //for intermediate probe sizes energyBottom = 20*maxEnergyAllowed_ 
    int   nConsecutiveCrossings;                                 //counts consecutive changes from compressive to constraint force and vice versa
    int   nConsecutiveSingleBasin;                               //counts consecutive # of times force separation is not needed (=> single basin)
        
    int   nOps = 0;                                                          //monitoring purposes only
    int   iterTotal = 0;
    int   maxRepeats = 10; 
    float prevEnergy = 1.0e8; 
    float max_dr_mag = 0.25*probe_;                                          // 1/8 of probe diameter or 1/4 the probe radius sets maximum dr_mag
    int   sFinal = numRadii_ - 1;                                            //corresponds to the probe radius of the user-specified probe size

    float dr_mag = 0.0;
    int   numExclusions = 0;                                                 //count # of times geometrical constraints from atoms are used

// --------------------------------------------------------------------------- start series of progression in probe size
    for(int s = sMin; s < numRadii_; s++) {                          // LOOP-A
    pR_ = radiusSchedule_[s];                                                //set new probe radius (used extensively in forceSeparation)
    pR2_ = pR_*pR_;                                                          //probe radius squared (used extensively in forceSeparation)
    float pD = 2*pR_;                                                        //probe diameter
    int nAtoms = maxAW_[s]; 
    forceSeparation(test_x,test_y,test_z,atomsWithin,nAtoms);                //separate compressive forces from constraint forces 
//                                                   ^^^^^^------>             maximum number of atoms within range for probe size defined by s
//                                                                             this function defines springX_, springY_, springZ_, springL_ 
       if( s != sFinal ){                                                    //smaller probe size implies reaching a 100% converged energy for VOID is not necessary
       energyBottom = 20*maxEnergyAllowed_;   //avoid unneccessary calculations to identify VOID gridpoint by setting target 20x higher for intermediate probe sizes 
       } else {                                                              //on the final probe size that is of interest
       energyBottom = maxEnergyAllowed_;                                     //for final probe size pay cost of convergence by setting energyBottom what it must be
       }

// --------------------------------------------------------------------------- initialize all variables for relaxation on new probe size
    float prev_ux = 0.0;
    float prev_uy = 0.0;
    float prev_uz = 0.0;
    int   countConverged = 0;
    float scale = inverseViscousDrag_; 
    nConsecutiveSingleBasin = 0;                              //counts consecutive # of times force separation is not needed (=> single basin)
    nConsecutiveCrossings = 0;                                //counts consecutive changes from compressive to constraint force and vice versa
    nSwitches = 0;                                            //counts # of times test particle switches between any of the three categories

// --------------------------------------------------------------------------- apply same method used for each probe size
    int iter = 0;
       while( iter < maxIterations_ ) {                              // LOOP-B A maximum number of iterations is granted per probe size  
       iter++;
// -------------------------------- calculate net compressive force exerted by neighboring atoms on test particle representing sphere-center
       energy = 0;
       float Fnet_x = 0.0;
       float Fnet_y = 0.0;
       float Fnet_z = 0.0;
// --------------------------------------------------------------------------- sum over compressive forces from all directions from all atoms
          for(int n = 0; n < nF_; n++){
          float dx = springX_[n] - test_x; 
          float dy = springY_[n] - test_y; 
          float dz = springZ_[n] - test_z; 
          r_mag = magnitude(dx,dy,dz);  
          float k1PushD = r_mag - springL_[n];                               //k1 implies the compressive spring constant equals 1 in some units 
             if( k1PushD < 0 ){                                              //compressed spring 
             energy += k1PushD*k1PushD;                                      //two times the spring energy (the factor of 2 is on purpose)
             float temp = k1PushD/r_mag;                                     //r_mag will never be close to zero. Will always be > vdw
             Fnet_x += temp*dx;
             Fnet_y += temp*dy;
             Fnet_z += temp*dz;
             nOps++;
             }
          }
// --------------------------------------------------------------------------- test if finished
          if( energy <= energyBottom ) {                                     //sufficiently low energy -> consider this as no stress (or go to next probe size)
       
             if( s != sFinal ){                                              //a smaller probe size than specified
             //void identified for smaller probe size                        //contnue to consider this test particle but using bigger probe size
             break;                                                          //break out of LOOP-B (iterations) back into LOOP-A for next probe size
             // return 1; is not correct, must fetch next probe size
             }
       

//          cavityIdentified[maxAtomsWithin]++;
//          springIterations[maxAtomsWithin] += (iterTotal + iter);            //all previous probe sizes plus current probe size
//          springOperations[maxAtomsWithin] += nOps;
          return 1;                                                          // *void* (approximately not stressed => overestimate of cavity) 
          }

       if( iter > 200 ) {                                                    //200 iterations should have converged. Want to stop this.
       float Edrop = 100*maxIterations_*maxEnergyAllowed_/iter;              // will stop if not close to true final energy level
          if( energy > Edrop ) {

//          springIterations[maxAtomsWithin] += (iterTotal + iter);      //all previous probe sizes plus current probe size
//          springOperations[maxAtomsWithin] += nOps;
//          springsStressed[maxAtomsWithin]++;
          return -1;                                                   // *microvoid*   we already know energy > maxEnergyallowed_
          }
       }
        
          if( nConsecutiveSingleBasin > 1 ) {                                //to stop: require at least 2-steps within the same basin
          if( nConsecutiveSingleBasin > 5 ){ nSwitches = 0; }                //=> switching has stopped with 5 in a row
          float dEneg = prevEnergy - energy;                                 //equals -1* [ standard definition of change: (final - initial) ]
// --------------------------------------------------------------------------- criteria is based on precision/accuracy interplay
            if( dEneg > 0.0 ) {                                              //check for early convergence criteria
            float ratio = dEneg/(energy - maxE_);                            //simple metric for convergence: ratio --> 0 => beter convergence
            float Ecut = ratio*conversionFactor_;                            //Ecut must be proportional to convergence-ratio 
              if( energy > Ecut ){                                           //PROBABLY convergence is reached and gridpoint is a microvoid
              countConverged++;                                              //count consecutive # of times convergence criterion is met
                 if( countConverged > goodConvergence_ ){                    //presume energy error estimate is negligible compared to the energy

                
//                 springIterations[maxAtomsWithin] += (iterTotal + iter);      //all previous probe sizes plus current probe size
//                 springOperations[maxAtomsWithin] += nOps;
//                 springsStressed[maxAtomsWithin]++;
                 return -1;                                                   // *microvoid*   we already know energy > maxEnergyallowed_
                 } 
              } else {                                                        //energy < Ecut indicates relaxation is not converged yet
              countConverged = 0;                                             //=> restart countConverged for consecutive tests
              }
            } else {
            countConverged++;                                                //count consecutive # of times convergence criterion is met
            }
          } else {                                                           //=> crossings between basins indicates no convergence
          countConverged = 0;                                                //=> restart countConverged for consecutive tests
          }

          if( nConsecutiveCrossings > maxRepeats ){                          //exceeded maximum allowed number of oscillations 

           
//          springIterations[maxAtomsWithin] += (iterTotal + iter);         //all previous probe sizes plus current probe size
//          springOperations[maxAtomsWithin] += nOps;
//          springsStressed[maxAtomsWithin]++;
          return -1;                                                      // *microvoid*   we already know energy > maxEnergyallowed_
          }                                                              

// --------------------------------------------------------------------------- continue relaxation process because no stop condition was satisfied
       float Fnet = FLT_MIN + magnitude(Fnet_x, Fnet_y, Fnet_z);               //determine net force on particle
// --------------------------------------------------------------------------- find direction test particle moves from its current postion
       float ux = Fnet_x/Fnet; 
       float uy = Fnet_y/Fnet;                                               
       float uz = Fnet_z/Fnet;                                               

          if( nConsecutiveSingleBasin > 1 ){                                 //conditional upon the test particle moving within a single basin
          wCos = prev_ux*ux + prev_uy*uy + prev_uz*uz;                       //Note:  100% backwards == -1 <= wCos <= 1 == 100% forward
             if( wCos < -0.95 ){                                             //approximately a complete reversal
             scale *= 0.348678;                                              //= (0.9)^10 equivalent to 10 consecutive steps with wCos = -1
             } else if( wCos > 0.95 ) {                                      //approximately a complete forward direction
             scale *= 1.234568;                                              //Note: 0.348678*(1.234668)^5 = 1 => 5 to 1 ratio rule 
             } else {
             scale *= (1 + 0.1*wCos);                                        //continuously adapt step size based on persistence angle 
             }
          }
       prev_ux = ux;
       prev_uy = uy;
       prev_uz = uz;
       dr_mag = scale*Fnet;                                                  //particle displacement is proportional to net force on particle

          if( dr_mag > max_dr_mag ){                                         //must limit step size to maximum step size allowed. 
          scale = max_dr_mag/Fnet;                                           //if an adjustment is necessary, redefined scale to be consistent
          dr_mag = max_dr_mag; 
          }

     

// --------------------------------------------------------------------------- take a trial step assuming no geometrical constraints block move
       float trial_x = test_x + ux*dr_mag;
       float trial_y = test_y + uy*dr_mag;
       float trial_z = test_z + uz*dr_mag;

// --------------------------------------------------------------------------- apply strict geometrical constraint from (x,y,z)-gridpoint
       float rx = trial_x - x;                                               //calc. position vector from (x,y,z)-gridpoint to test particle
       float ry = trial_y - y;
       float rz = trial_z - z;
       float r2 = rx*rx + ry*ry + rz*rz;
          if( r2 > pR2_ ){                                                   //when test particle moves pass sphere radius, place it on the surface
          r_mag = sqrt(r2);
          float temp = pR_/r_mag; 
// --------------------------------------------------------------------------- location of test particle on surface of probe defines the new "current" position
          float dx = temp*rx;
          float dy = temp*ry;
          float dz = temp*rz;
          test_x = x + dx;
          test_y = y + dy;
          test_z = z + dz;
// --------------------------------------------------------------------------- calculate the current energy where the test particle is located on surface
          energy = 0;
             for(int n = 0; n < nAtoms; n++){                                //treat all forces as compressive forces (no geometrical constraints)
             dx = springX_[n] - test_x; 
             dy = springY_[n] - test_y; 
             dz = springZ_[n] - test_z; 
             r_mag = magnitude(dx,dy,dz);  
             float k1PushD = r_mag - springL_[n];                            //k1 implies the compressive spring constant equals 1 in some units 
                if( k1PushD < 0 ){                                           //compressed spring 
                energy += k1PushD*k1PushD;                                   //two times the spring energy (the factor of 2 is on purpose)
                nOps++;                                                      //count these evaluations as a full operations, albeit they are light
                }
             }
// =========================================================================== apply funnel diffusion to guard against slow down within a local minimum
          float maxStepSize = 0.062832*pR_;                                  //set maximum ds = (2*pi*radius)/100 => 1% of circumference 
          float stepSize = maxStepSize;                                      //start with maximum ds 
          bool  flagRepeat = false;
          float stopPoint = 250;
          float stoppingSum = 0;
             while(stepSize > drMinMagnitude_ ) {                    // LOOP-D performs funnel diffusion until step size reaches acceptable tolerance
            
                if( energy < energyBottom ) {                                //sufficiently low energy -> consider this as no stress (or go to next probe size)
                break;                                                       //break out of LOOP-D to stop funnel diffusion immediately    
                } else {                                                     //randomly generate another trial position for the test particle
                   if( flagRepeat == false ) {
                   dx = stepSize*( 0.1 + ((float)rand()/RAND_MAX) ); 
                   dy = stepSize*( 0.1 + ((float)rand()/RAND_MAX) ); 
                   dz = stepSize*( 0.1 + ((float)rand()/RAND_MAX) ); 
                   }
                trial_x = test_x + dx; 
                trial_y = test_y + dy; 
                trial_z = test_z + dz; 
                rx = trial_x - x;                                            //calc. position vector from (x,y,z)-gridpoint to test particle
                ry = trial_y - y;
                rz = trial_z - z;
                r2 = rx*rx + ry*ry + rz*rz;
                   if( r2 > pR2_ ){                                          //place particle on surface when it tries to escape
                   r_mag = sqrt(r2);
                   temp = pR_/r_mag;
                   trial_x = x + temp*rx;
                   trial_y = y + temp*ry;
                   trial_z = z + temp*rz;
                   } // else                                                 //trial position is within sphere
                }
// --------------------------------------------------------------------------- sum over compressive forces from all directions from all atoms
             float trial_energy = 0;
                for(int n = 0; n < nAtoms; n++){                             //treat all forces as compressive forces (no geometrical constraints)
                float dx = springX_[n] - trial_x; 
                float dy = springY_[n] - trial_y; 
                float dz = springZ_[n] - trial_z; 
                r_mag = magnitude(dx,dy,dz);  
                float k1PushD = r_mag - springL_[n];                         //k1 implies the compressive spring constant equals 1 in some units 
                   if( k1PushD < 0 ){                                        //compressed spring 
                   trial_energy += k1PushD*k1PushD;                          //two times the spring energy (the factor of 2 is on purpose)
                   nOps++;   
                   }
                }

                if( trial_energy < energy ){                                 //energy is reduced, so the trial is a sucessful step
                energy = trial_energy;
                flagRepeat = true;
                dx = trial_x - test_x;
                dy = trial_y - test_y;
                dz = trial_z - test_z;
                dx = 1.05*dx;
                dy = 1.05*dy;
                dz = 1.05*dz;
                stepSize = magnitude(dx,dy,dz);
                test_x = trial_x;
                test_y = trial_y;
                test_z = trial_z;
                stoppingSum = 0;
                } else {                                                     //failed attempt (must reduce step_size)
                   if( flagRepeat == true ) {
                   stepSize *= 0.666667;                                     //cut step size down by 1/3 when a correlated step fails
                   flagRepeat = false; 
                   } else {
                   stepSize *= 0.90;                                         //about 7 uncorrelated failed steps reduces stepSize by 1/2  0.9^7 = ~0.478 
                   }
                temp = energy/maxEnergyAllowed_;                             //this ratio of energies is based on the true final energy level we need to reach
                if( temp > 25 ) { temp = 25; }
                stoppingSum += temp;
                   if( stoppingSum > stopPoint ){
                   break;                                                    //break out of LOOP-D to stop funnel diffusion immediately    
                   }
                }
             }                                                       // LOOP-D completes the funnel diffusion process
// --------------------------------------------------------------------------- check if energy is low enough
             if( energy < energyBottom ){                                    //sufficiently low energy -> consider this as no stress (or go to next probe size)
       
       
                if( s != sFinal ){                                           //a smaller probe size than specified
                //void identified for smaller probe size                     //contnue to consider this test particle but using bigger probe size
                break;                                                       //break out of LOOP-B (iterations) back into LOOP-A for next probe size
                // return 1; is not correct, must fetch next probe size
                }

           

//             cavityIdentified[maxAtomsWithin]++;
//             springIterations[maxAtomsWithin] += (iterTotal + iter);         //all previous probe sizes plus current probe size
//             springOperations[maxAtomsWithin] += nOps;
             return 1;                                                       // *void* (approximately not stressed => overestimate of cavity) 
             } else {
                                                                            
//             springIterations[maxAtomsWithin] += (iterTotal + iter);         //all previous probe sizes plus current probe size
//             springOperations[maxAtomsWithin] += nOps;
//             springsStressed[maxAtomsWithin]++;                              //appears to be stressed, albeit not converged satisfactorily
             return -1;                                                      //microvoid --> energy too high and spatial resolution converaged
             } 
// ===========================================================================
          } else {                                                           //stays within accessible spherical region 
// --------------------------------------------------------------------------- check if test particle moves into an excluded region of space
//        Geometrical picture:     I = initial position of test particle       F = final position of test particle       A = postion of atom
//
//                      F                                NOTE: See https://en.wikipedia.org/wiki/Law_of_cosines
//                   */                                  ux,uy,uz is the unit vector that points from I to F
//                   /                                   B = angle between FIA
//                  /     *                              q = the distance between the atom at A to the test particle at initial position I.
//                 /                                     qCosB = (distance from I to A)*cos(B)
//                /             *                        L = length = springL_[n] is the distance from A to F (maximum interaction length)
//               /                                       rAF = distance between atom A and the final position of the test particle at F.
//              /                     *                  F is a trial position of the test particle, considered here as the final postion. 
//             /  B                                      Use law of cosines combined with Pythagorean Theorem to arrive at:
//            I_____________________________* A          x = distance between (I-F) = qCosB - sqrt(L*L - q2SinB2) 
//                          q
// --------------------------------------------------------------------------- apply geometrical constraints that exclude regions of space
          //cout << "intitial:  dr_mag= " << dr_mag << endl; 
          float d1 = dr_mag;                                                 //lowest step size before moving into an excluded region
          float d2 = dr_mag;                                                 //next lowest step size before moving into an excluded region
          numExclusions = 0;                                
             for(int n = nC_; n < nAtoms; n++){                      // LOOP-C check all constraint forces
             float dx = springX_[n] - trial_x; 
             float dy = springY_[n] - trial_y; 
             float dz = springZ_[n] - trial_z; 
             float rAF2 = dx*dx + dy*dy + dz*dz;                                //square of the distance AF.
             float L2 = springL_[n]*springL_[n];
                if( rAF2 < L2 ){                                                //test particle entered excluded region of space: Must stop this.
                dx = springX_[n] - test_x;                                      //working with distance q = sqrt(dx^2 + dy^2 + dz^2)
                dy = springY_[n] - test_y; 
                dz = springZ_[n] - test_z; 
                float q2 = dx*dx + dy*dy + dz*dz;
                float qCosB = dx*ux + dy*uy + dz*uz;                            //(distance from atom to test particle) * cos(angle)
                float q2SinB2 = q2 - qCosB*qCosB;                               //square of cos(angle) = 1 - sin(angle) * sin(angle)
                float sq_root = 0.0;
                   if( L2 > q2SinB2 ){                                          //two possible solution exists where lines I-F cross with A-F.
                   sq_root = sqrt( L2 - q2SinB2 );
                   } 
                numExclusions++;                                
                float temp = qCosB - sq_root;                                   //new dr_mag, smallest one requires negative square root
                   if( temp < 0 ){                                              //always take the positive solution (in direction of dr_mag)
                   temp = qCosB + sq_root; 
                  
                   }
                   if( temp < d2 ){                                             //among all restrictions, take the smallest case found
                   d2 = temp;                                                   //record the second smallest case too
                      if( d2 < d1 ) {                                           //enforce d1 <= d2
                      d2 = d1;
                      d1 = temp;
                      }
                   }
                }
             }                                                          // LOOP-C end checking over all constraint forces
             if( numExclusions > 0 ) {
             if( nConsecutiveCrossings == 0 ){ nSwitches++; }
             nConsecutiveCrossings++;                                           //allow test particle to cross into closest excluded region
             nConsecutiveSingleBasin = 0;                                       //resets to 0 whenever there is a crossing into another basin
             dr_mag = 0.95*d1 + 0.05*d2;                                        //put test particle just on other side of excluded boundary
             test_x += ux*dr_mag;
             test_y += uy*dr_mag;
             test_z += uz*dr_mag;
             forceSeparation(test_x,test_y,test_z,atomsWithin,nAtoms);          //must update separation of compressive and constraint forces 
             } else {                                                           //dr_mag is correct because test particle remains in same basin
             test_x = trial_x;
             test_y = trial_y;
             test_z = trial_z;
             if( nConsecutiveSingleBasin == 0 ){ nSwitches++; }
             nConsecutiveSingleBasin++;
             nConsecutiveCrossings = 0;      
             }
          }

          if( (dr_mag < drMinMagnitude_) && (nConsecutiveCrossings == 0) ) {    //=>spatial resolution achieved! => finish as a microvoid
                                                                                //do not measure dr_mag during basin jumps
         
//          springIterations[maxAtomsWithin] += (iterTotal + iter);            //all previous probe sizes plus current probe size
//          springOperations[maxAtomsWithin] += nOps;
//          springsStressed[maxAtomsWithin]++;                                 //appears to be stressed, albeit not converged satisfactorily
          return -1;                                                         //microvoid --> energy too high and spatial resolution converaged
          }

          if( nSwitches > 15 ){                                              //test particle can't decide if its coming or going.
                                                                                //do not measure dr_mag during basin jumps
       

//          springIterations[maxAtomsWithin] += (iterTotal + iter);            //all previous probe sizes plus current probe size
//          springOperations[maxAtomsWithin] += nOps;
//          springsStressed[maxAtomsWithin]++;                                 //appears to be stressed, albeit not converged satisfactorily
          return -1;                                                         //microvoid --> energy too high and spatial resolution converaged
          }


       nOps++;                                                               //moving test particle plus checking, together give 1 operation
       prevEnergy = energy;                                                  //prepare for next iteration by storing current energy
       }                                                              //LOOP-B maximum iterations for a given probe size exceeded 
    iterTotal += iter;                                                       //add all iterations from previous probe size to track total
    }                                                                 //LOOP-A finishes series of progressively larger probe radii
//    testParticleOut[maxAtomsWithin]++;                                       //counts out of loop decisions based on non-converged configuration

    

//    springIterations[maxAtomsWithin] += iterTotal;
//    springOperations[maxAtomsWithin] += nOps;
//    springsStressed[maxAtomsWithin]++;                                       //appears to be stressed, albeit not converged satisfactorily
    return -1;                                                               //assign gridpoint as microvoid due to lack of convergence 
}


inline void SpringModel::forceSeparation(float x, float y, float z, int * atomsWithin, int Na) {  
// Na = # of atoms that may produce a compressive or constraint force depending on probe size and its location, defined by the test particle.
// The array  atomsWithin[]  define the maximum set of atoms that can potentially create a force on the test particle. 
// x,y,z are the coordinates of the test particle.
// This program separates compressive forces from constraint forces. It also filters out all atoms that do not produce force on test particle. 
// The atoms that do not produce force are irrlevant and this happens when a radius smaller than the specified radius is considered. 
//
// --------------------------------------------------------------------------------------------------- explanation of how information is stored
//  OLD METHOD:
//  Consider the indexing for the set of arrays:  springX_, springY_, springZ_, springL_, etc
//  Example:    list of interactions: n= 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 = maxAtomsWithin
//                                       |              |                                 |                                |
//                    list of forces: n= 1  2  3  ...   ^=nF                              |                                |
//   list of geometrical constraints:                                                     ^=nC             ...             ^=maxAtomsWithin
//   remainder is currently irrelevant                     |<-- cannot exert force -->|
//   Calculate compressive forces using indices from 1 to nF.  If test particle is located in excluded space, use the net force to push it out.
//   Calculate geometrical constraints using indices from nc to maxAtomsWithin. If test particle is out of excluded space, do not let it enter.
//   int nC;                                           use to track geometrical constraints to prevent entry of test particle in excluded space
//   int nF;                                           use to track compressive spring-forces that push the test particle out of excluded space
//
//  NEW METHOD:
//  The new method is exactly the same as the OLD method, but because of presorting (at initialization) no irrelevant cases appear as probe
//  sizes increase. As such the new method does not need to peform the  *NO*  lines that have been commented out below. Moreover, the private
//  variable   radius2_[]   is now obsolete. Throughout this code, it is commented out with the note:  *NO*  added as a reference to here.
//  The  *NO*  lines have not been deleted, and the explanation above does not assume presorted so that the basic algorithm can be understood.
//  Now the advantage of presorting can be appreciated. Note that because irrelvant cases have been pushed to the far right, the new line is:
//
//                                                         nC                               Na = number of atoms within for this probe size
//                                                         |                                |
//  Example:    list of interactions: n= 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 = maxAtomsWithin
//                                       |              |                                                                 |
//                    list of forces: n= 1  2  3  ...   ^=nF                                |<- irrelavent atoms within ->|
// Note: nC = nF + 1 
// It is worth noting that before presorting was performed, this function refers to the ORIGINAL list of atoms within range to affect the probe. 
// In the new method the ORIGINAL list is actually the presorted list. 
// Another interesting fact is this new method would work on the ORIGINAL list that is not presorted, but the difference is there would 
// be more constraint forces. That is, all irrelavent atoms would be considered constraint forces. While this does not change the results in 
// any way, the more information kept that is irrelevant, the slower the program takes for NO REASON. As such, the goal of optimization 
// takes advantage of precaculations at the expense of complex program. 
// --------------------------------------------------------------------------- filter out springs that are out of reach from the get-go
    nF_ = 0;                                                                 //counts the number of compressive forces acting on the test particle
    nC_ = Na;                                                                //starting index for list of geometrical constraints on test particle
       for(int n = 0; n < Na; n++){                                          //search through all springs to identify their current status
// *NO*  float span = 2*pR_ + springR_[n];                                   //Note: springR_[n] refers to the n-th atom in the ORIGINAL list

// *NO*   if( radius2_[n] <= span*span ) {                                   //the n-th atom restricts location of a probe of radius pR_
          float length = pR_ + springR_[n];                                  //defines range of influence of the n-th atom in the ORIGINAL list
          float dx = atomCoords_[atomsWithin[n]][1] - x; 
          float dy = atomCoords_[atomsWithin[n]][2] - y; 
          float dz = atomCoords_[atomsWithin[n]][3] - z; 
          float r2 = dx*dx + dy*dy + dz*dz;                                  //squared distance between test particle & n-th atom of ORIGINAL list
          float r2min = length*length;                                       //A dividing distane between a compressive interaction or not
             if( r2 > r2min ) {                                              //test particle is currently outside of excluded region by atom
             nC_--;
             springX_[nC_] = atomCoords_[atomsWithin[n]][1];
             springY_[nC_] = atomCoords_[atomsWithin[n]][2];
             springZ_[nC_] = atomCoords_[atomsWithin[n]][3];
             springL_[nC_] = length;                                         //minimum distance tolerated before test particle enters excluded space 
             } else {                                                        //test particle currently finds itself in excluded space (bad place)
             springX_[nF_] = atomCoords_[atomsWithin[n]][1];
             springY_[nF_] = atomCoords_[atomsWithin[n]][2];
             springZ_[nF_] = atomCoords_[atomsWithin[n]][3];
             springL_[nF_] = length;                                         //equilibrium length of compressive spring used to push test particle out 
             nF_++;
             }
// *NO*   } // else -->                                                        this atom is irrelevant for this probe size because it is too far out
       }

}


inline float SpringModel::magnitude(float x, float y, float z) {
    return sqrt(x*x + y*y + z*z);
}


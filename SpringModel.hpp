/* 
 * File:   SpringModel.hpp
 * Author: jenny
 *
 * Created on November 17, 2016, 8:26 PM
 */

#ifndef SPRINGMODEL_HPP
#define	SPRINGMODEL_HPP

#include <boost/random.hpp>

using namespace std;

class SpringModel {
public:
/*    
    int * springIterations;
    int * springOperations;  
    int * springNeighbors;
    int * testParticleOut; 
    int * springsStressed;  
    int * cavityIdentified; 
*/       
    SpringModel(int maxAtomsWithin, float grid, float probe, float ** atomCoords);    
    virtual ~SpringModel();
    void setAtomCoordinates(float ** atomCoords);
    int testSpring(float x, float y, float z, int * atomsWithin, int nAtomsWithin);

    
private:    
    static const int maxIterations_ = 1000;         //plateau value of ~900 was determined by studying dependence from 45 to 2000
    static const float ref_maxDrift_ = 0.0125;      //reference tolerance in position (in Angstroms) for grid = 0.5 Angstroms  (2.5%)
    static const float inverseViscousDrag_ = 0.50;  //sets proportionality between force and displacement: largest value to use ~0.20
                                                    //larger values of inverseViscousDrag_ speeds up calculations, but decreases accuracy 
                                                    //lower accuracy at 0.75 is still good, but not advisable if accuracy is a concern. 
                                                    //a value of 0.125 works well for both speed and accuracy on one-shot probe sizes.
                                                    //gains in accuracy for lower values than 0.125 is not worth the extra calculation time
    static const float conversionFactor_ = 500;     //emperical experimentation shows 100 yields low errors, 500 has 1/5 as much error
                                                    //converts simple unitless convergence-ratio to Ecut used to truncate iteration loop. 
    static const float shell_rate_ = 5.0/2.5;       //use 5 intermediate (radii) before reaching Rprobe= 2.5 Angstroms
    static const float goodConvergence_ = 5;        //number of consecutive tests that says energy estimate is converged
    
    inline void forceSeparation(float x, float y, float z, int * atomsWithin, int Na);

    int     numRadii_;             //number of radii used in growing size of probe until desired final value
    float   * radiusSchedule_;     //stores list of probe radii for testing with growing the size of the probe
    int     * maxAW_;              //stores maximum number of atoms within sphere of influence for various probe size radii

    float drMinMagnitude_;         //based on set tolerance on MaxDallowed_
    float maxEnergyAllowed_;

    float maxE_;                   //a freqently used reference energy.   maxE_ = maxEnergyAllowed_ - dEzero_
    float grid_;  
    float probe_;
    float ** atomCoords_;    

// --------------------------------//variables to hold properties of all geometrical constraints and spring interactions
    int     nF_;
    int     nC_;
    float   pR_;
    float   pD_;
    float   pR2_;
    float * springR_;
    float * springX_;
    float * springY_;
    float * springZ_;
    float * springL_;
    float * radius2_;
    
    inline static float magnitude(float x, float y, float z);
};

#endif	/* SPRINGMODEL_HPP */


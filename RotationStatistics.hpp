/* 
 * File:   RotationStatistics.hpp
 * Author: jenny
 *
 * Created on November 26, 2016, 4:11 PM
 */

#ifndef ROTATIONSTATISTICS_HPP
#define	ROTATIONSTATISTICS_HPP

#include <fstream>
#include <string>
#include <cmath>

using namespace std;

class RotationStatistics {
public:
    RotationStatistics(bool flagOpt, float grid, float probe, float shell, float cutoff);
    virtual ~RotationStatistics();
    
    void rotateRun(int nRotations, string PDBFileName, string outputPath, int fixedProbe);
    
private:
 
    bool   flagOpt_;
    float  grid_;
    float  probe_;
    float  shell_;
    float  cutoff_;
  
    float * probeSizes_;
    void calculateProbeSizes(int nRotations, int fixedProbe);

};

#endif	/* ROTATIONSTATISTICS_HPP */


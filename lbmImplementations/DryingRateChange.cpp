
/************************************************************************************/
/*  Copyright 2025. Corning Incorporated. All rights reserved.                      */                                                                                     #
/*  This software may only be used in accordance with the identified license(s).    */  
/*                                                                                  */                                                                                      
/*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      */ 
/*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        */
/*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL         */
/*  CORNING BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN      */
/*  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN               */
/*  CONNECTION WITH THE SOFTWARE OR THE USE OF THE SOFTWARE.                        */
/************************************************************************************/
/*  Authors:                                                                        */
/* Hamed Haddadi Staff Scientist                                                    */
/*               haddadigh@corning.com                                              */
/* David Heine   Principal Scientist and Manager                                    */
/*               heinedr@corning.com                                                */
/************************************************************************************/

# include "../lbmDeclarations/DryingRateChange.h"


void DryingRateChange::setCohesionValues(const CohesionRangeParams<T> & cohesionRangeParams,
         const std::string & gChangeType) {
    T g00{0.}, g11{0.}, gMin{0.}, gMax{0.};
    g00 = cohesionRangeParams.g00;
    g11 = cohesionRangeParams.g11;
    gMin = cohesionRangeParams.gMin;
    gMax = cohesionRangeParams.gMax;

    if (gChangeType == "range") {
        T g{0}, gInc{0};
        numGSteps_ = cohesionRangeParams.gSteps;
        gValues_.resize(numGSteps_);
        gInc = (gMax - gMin)/(T) numGSteps_;
        for (plint n = 0; n < numGSteps_; n++) {
            g = gMin + n*gInc;
            gValues_.at(n) = g;
        }
    }

    else if (gChangeType == "step") {
        numGSteps_ = 2;
        gValues_.resize(numGSteps_);
        gValues_.at(0) = gMin;
        gValues_.at(1) = gMax;
        gChangeStep_ = cohesionRangeParams.gSteps;
    }

   
    // initialize spG
    spG_.resize(2, std::vector<T>(2));
    spG_.at(0).at(0) = g00;
    spG_.at(0).at(1) = 0.0;
    spG_.at(1).at(0) = 0.0;
    spG_.at(1).at(1) = g11;     

}

void DryingRateChange::setFluidsProperties(const CohesionRangeParams<T> & cohesionRangeParams, 
                        const FluidsParams<T> & fluidsParams, const std::string & gChangeType) {
    // method for changing G and constant omega values 
    T omegaF1{0}, omegaF2{0};
    setCohesionValues(cohesionRangeParams, gChangeType);
    omegaF1 = fluidsParams.omegaF1; 
    omegaF2 = fluidsParams.omegaF2;
    gF1S_ = fluidsParams.gF1S;

    for (plint numG = 0; numG < numGSteps_; ++numG) {
        std::vector<T> omegaVals{omegaF1, omegaF2};
        omegaValues_.push_back(omegaVals);
    }

}

void DryingRateChange::setFluidsProperties(const CohesionRangeParams<T> & cohesionRangeParams,
                 const FluidsParamsChangingOmega<T> & fluidsRangeParams, 
                    const std::string & gChangeType) {
    // method for changing G and Omega 
    T omegaMinF1{0}, omegaMaxF1{0}, omegaMinF2{0}, omegaMaxF2{0};
    T omegaF1Inc{0}, omegaF2Inc{0};
    T omegaF1{0.}, omegaF2{0.};
    
    setCohesionValues(cohesionRangeParams, gChangeType);

    omegaMinF1 = fluidsRangeParams.omegaMinF1;
    omegaMaxF1 = fluidsRangeParams.omegaMaxF1;
    omegaMinF2 = fluidsRangeParams.omegaMinF2; 
    omegaMaxF2 = fluidsRangeParams.omegaMaxF2;

    // change type of viscosity is similar to change type of cohesion 
    if (gChangeType == "range") {
        omegaF1Inc = (omegaMaxF1 - omegaMinF1)/(T) numGSteps_;
        omegaF2Inc = (omegaMaxF2 - omegaMinF2)/(T) numGSteps_;
        for (plint numG = 0; numG < numGSteps_; numG++) {
            omegaF1 = omegaMinF1 + omegaF1Inc; 
            omegaF2 = omegaMinF2 + omegaF2Inc; 
            std::vector<T> omegaVals{omegaF1, omegaF2};
            omegaValues_.push_back(omegaVals);
        }
    }

    else if (gChangeType == "step") {
        omegaValues_.resize(2, std::vector<T> (2, 0.));
        omegaValues_.at(0).at(0) = omegaMinF1;
        omegaValues_.at(0).at(1) = omegaMinF2;
        omegaValues_.at(1).at(0) = omegaMaxF1;
        omegaValues_.at(1).at(1) = omegaMaxF2;
    }

}


void DryingRateChange::setShanChen(T gValue, std::vector<T> omegaValues) {
    
    std::vector <MultiBlockLattice3D<T, MPDESCRIPTOR> *> blockLattices;    
    plint processorLevel = 1;
    blockLattices.push_back(& latticeFluidTwo_);
    blockLattices.push_back(& latticeFluidOne_);    
    
    spG_.at(0).at(1) = gValue; 
    spG_.at(1).at(0) = gValue;

    std::vector<T> constOmegaValues;
    constOmegaValues.assign({omegaValues.at(0), omegaValues.at(1)});

    integrateProcessingFunctional(new ShanChenMultiComponentProcessor3D <T, MPDESCRIPTOR> (spG_, constOmegaValues),
         Box3D(0, nx_ - 1, 0, ny_ - 1, 0, nz_ - 1), blockLattices, processorLevel);
}

void DryingRateChange::setUp() {
    // note that setShanChen is not called here
    readGeometry();
    setInletOutletDensities();
    initBoundaryPlanes();
    initPressureBC();
    // note that these values must be equal for the equilibrium stage
    setPressureBoundaryValues(rhoInitInlet_, rhoInitOutlet_);
    defineLatticeDynamics();
    initializeLatticeDensities();
    addExternalForces();
    initializeLattices();
}

void DryingRateChange::runPressureRamp(plint maxRampIter, plint outputFreq, plint checkFreq, T convCr) {

 //   pcout <<"performing the pressure ramp stage >>> "<<std::endl;
    plint iT{0}, gRampIter{0};
    T newAvgEnF1{}, newAvgEnF2{}, oldAvgEnF1{1.}, oldAvgEnF2{1.};
    T volume = (T)(nx_*ny_*ny_);
    std::vector<plint> gRampIters(numGSteps_);

    setPressureBoundaryValues(inletRhoValues_[1], outletRhoValues_[1]);           
    
    oldAvgEnF1 = 1.0;
    oldAvgEnF2 = 1.0;

    if (gChangeStep_ == 0) {
        gRampIter = maxRampIter/numGSteps_;
        for (plint it = 0; it < numGSteps_; ++it) {
            gRampIters.at(it) = gRampIter;
        }
    }
    else if (gChangeStep_ != 0) {
        gRampIters.at(0) = gChangeStep_;
        gRampIters.at(1) = maxRampIter - gChangeStep_;
    }

    for (plint numG = 0; numG < numGSteps_; ++numG) {
        
        T gValue = gValues_.at(numG);
        std::vector<T> omegaValues = omegaValues_.at(numG);
        setShanChen(gValue, omegaValues);
        gRampIter = gRampIters.at(numG);

        for (iT = 0; iT < gRampIter; ++iT) {
            latticeFluidOne_.collideAndStream();
            latticeFluidTwo_.collideAndStream();

            if (iT % outputFreq == 0) {
                writeRhoVTK(outCounter_);
                writeVelocityComponentsDAT(outCounter_);
                writeRhoDistributionDAT(outCounter_);
                ++outCounter_;   
                }

            if (iT % checkFreq == 0) {
                newAvgEnF1 = getStoredAverageDensity(latticeFluidOne_)*(volume);
                newAvgEnF2 = getStoredAverageDensity(latticeFluidTwo_)*(volume);
                oldAvgEnF1 = newAvgEnF1;
                oldAvgEnF2 = newAvgEnF2;
        }
    }
    }
}

void DryingRateChange::operator()(plint maxIter, plint maxRampIter, plint outputFreq, plint checkFreq, T convCr) {
    // uses a g value to run the equilibrium stage without evaporation
    setShanChen(1.0, omegaValues_.at(0));
    setUp();
    runEquilibrium(maxIter, outputFreq, checkFreq, convCr);
    runPressureRamp(maxRampIter, outputFreq, checkFreq, convCr);
}


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
// implementations of methods defined in MultiPhasePressure 
# include "../lbmDeclarations/MultiPhasePressure.h"


void MultiPhasePressure::setInletOutletDensities() {
    T cosTheta = std::abs(4.*gF1S_/(gc_*(rhoInitInlet_ - rhoNoFluid_)));
    T deltaRho = 6.0*sigma*cosTheta/minRadius_;
    T stepSize = deltaRho/totalNumRuns_; 

    for (plint runNum = 0; runNum <= totalNumRuns_; ++runNum) {
        outletRhoValues_.push_back(rhoInitOutlet_ - (T)runNum*stepSize);
        inletRhoValues_.push_back(rhoInitInlet_);

    }
}

void MultiPhasePressure::initPressureBC() {
    
    boundaryCondition_ = createLocalBoundaryCondition3D<T, MPDESCRIPTOR>();
    boundaryCondition_ -> addPressureBoundary0N(inlet_, latticeFluidOne_);
    boundaryCondition_ -> addPressureBoundary0N(inlet_, latticeFluidTwo_);

    boundaryCondition_ -> addPressureBoundary0P(outlet_, latticeFluidOne_);
    boundaryCondition_ -> addPressureBoundary0P(outlet_, latticeFluidTwo_);
}

void MultiPhasePressure::setPressureBoundaryValues(T rhoInlet, T rhoOutlet) {
    // inlets
    setBoundaryDensity(latticeFluidOne_, inlet_, rhoInlet);
    setBoundaryDensity(latticeFluidTwo_, inlet_, rhoNoFluid_);
    // outlets 
    setBoundaryDensity(latticeFluidOne_, outlet_, rhoNoFluid_);
    setBoundaryDensity(latticeFluidTwo_, outlet_, rhoOutlet);
}


void MultiPhasePressure::setUp() {
    setShanChen();
    readGeometry();
    setInletOutletDensities();
    initBoundaryPlanes();
    initPressureBC();
    setPressureBoundaryValues(rhoInitInlet_, rhoInitOutlet_);
    defineLatticeDynamics();
    initializeLatticeDensities();
    addExternalForces();
    initializeLattices();
}


void MultiPhasePressure::operator()(plint maxIter, plint checkFreq, plint outputFreq, T convCr) {
    // outputs are generated at the end of each converged step
    setUp();
    bool hasNotConverged{true};
    plint iT{0}, numOut{0}, totalNumIter{0}, numIter{0};
    T newAvgEnF1{}, newAvgEnF2{}, oldAvgEnF1{1.}, oldAvgEnF2{1.};
    T cyclePressure{0.};
    T volume = (T)(nx_*ny_*ny_);

  
    for (plint numRun = 0; numRun < totalNumRuns_; ++numRun) {
        if (numRun > 0) {
            setPressureBoundaryValues(inletRhoValues_[numRun], outletRhoValues_[numRun]);        
        }
        cyclePressure = (1./3.)*(inletRhoValues_.at(numRun) - outletRhoValues_.at(numRun));
        hasNotConverged = true;
        oldAvgEnF1 = 1.0;
        oldAvgEnF2 = 1.0;
        iT = 0;
        while (hasNotConverged) {

            latticeFluidOne_.collideAndStream();
            latticeFluidTwo_.collideAndStream();

            if (totalNumIter % outputFreq == 0) {
                writeRhoVTK(numOut);
                writeVelocityComponentsDAT(numOut);   
                writeRhoDistributionDAT(numOut);
                pressureValues_.push_back(cyclePressure);
                ++numOut;                 
            }

            if (iT % checkFreq == 0) {
                newAvgEnF1 = getStoredAverageDensity(latticeFluidOne_)*(volume);
                newAvgEnF2 = getStoredAverageDensity(latticeFluidTwo_)*(volume);
             
                if (simutils::hasConverged(oldAvgEnF1, oldAvgEnF2, newAvgEnF1, newAvgEnF2, (T) checkFreq, convCr)) {
                    hasNotConverged = false;
                }
                oldAvgEnF1 = newAvgEnF1;
                oldAvgEnF2 = newAvgEnF2;
            }

            if (iT >= maxIter) {
                hasNotConverged = false;
            }
            ++iT;
            ++totalNumIter;
        }
    }
    writeSimulationDatFile();
}


void MultiPhasePressure::writeSimulationDatFile() {
    std::string simFile = outputDir_ + "simulation.dat";
    plb_ofstream simInfo(simFile.c_str());
    plint pSize = pressureValues_.size();

    addSimulationGeneralInfo(simInfo);
    simInfo<<"delta_P: ";
    
    for (plint numP = 0; numP < pSize; ++numP) {
        simInfo<<" "<<pressureValues_.at(numP);
        if (numP == pSize - 1) {
            simInfo<<std::endl;
        }
    }
}

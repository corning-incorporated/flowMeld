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


# include "../lbmDeclarations/MultiPhaseRunOut.h"


void MultiPhaseRunOut::setInletOutletDensities() {
    T cosTheta = std::abs(4.*gF1S_/(gc_*(rhoInitInlet_ - rhoNoFluid_)));
    T deltaRho = 6.0*sigma*cosTheta/minRadius_;
    T stepSize = deltaRho/totalNumRuns_; 

    for (plint runNum = 0; runNum <= totalNumRuns_; ++runNum) {
        outletRhoValues_.push_back(rhoInitOutlet_ - (T)runNum*stepSize);
        inletRhoValues_.push_back(rhoInitInlet_);
    }    
}

void MultiPhaseRunOut::initPressureBC() {
    
    boundaryCondition_ = createLocalBoundaryCondition3D<T, MPDESCRIPTOR>();
  //  boundaryCondition_ -> addPressureBoundary0N(inlet_, latticeFluidOne_);
  //  boundaryCondition_ -> addPressureBoundary0P(outlet_, latticeFluidOne_);

    boundaryCondition_ -> addPressureBoundary0N(inlet_, latticeFluidTwo_);
    boundaryCondition_ -> addPressureBoundary0P(outlet_, latticeFluidTwo_);

}


void MultiPhaseRunOut::setPressureBoundaryValues(T rhoInlet, T rhoOutlet) {
    // tag 0: fluid Two; tag 3: fluid One 
    setBoundaryDensity(latticeFluidTwo_, inlet_, rhoInlet);
    setBoundaryDensity(latticeFluidTwo_, outlet_, rhoOutlet);

  //  setBoundaryDensity(latticeFluidOne_, inlet_, rhoNoFluid_);   
  //  setBoundaryDensity(latticeFluidOne_, outlet_, rhoNoFluid_);

}

void MultiPhaseRunOut::setUp() {
    // this is an initial setup for the imbibition stage 
  //  setShanChen();
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

void MultiPhaseRunOut::runEquilibrium(plint maxIter, plint outputFreq, plint checkFreq, T convCr) {
    // to simulate the initial imbibition stage 
    bool hasNotConverged{true};
    plint iT{0};
    T newAvgEnF1{}, newAvgEnF2{}, oldAvgEnF1{1.}, oldAvgEnF2{1.};

    pcout <<"performing the initial imbibition stage >>> "<<std::endl;

    for (iT = 0; iT < maxIter; ++iT) {
        latticeFluidOne_.collideAndStream();
        latticeFluidTwo_.collideAndStream();
        
        if ((iT % checkFreq == 0) && (hasNotConverged)) {
            newAvgEnF1 = getStoredAverageDensity(latticeFluidOne_);
            newAvgEnF2 = getStoredAverageDensity(latticeFluidTwo_);
           
            if (simutils::hasConverged(oldAvgEnF1, oldAvgEnF2, newAvgEnF1, newAvgEnF2, (T) checkFreq, convCr)) {
                hasNotConverged = false;
            }
            oldAvgEnF1 = newAvgEnF1;
            oldAvgEnF2 = newAvgEnF2;
        }

        if ((iT % outputFreq == 0)) {
            writeRhoVTK(outCounter_);
            writeVelocityComponentsDAT(outCounter_);
            writeRhoDistributionDAT(outCounter_);
            ++outCounter_;
        }

    }
    pcout <<"initial imbibition stage finished "<<std::endl;
}


void MultiPhaseRunOut::runPressureRamp(plint maxRampIter, plint outputFreq, plint checkFreq, T convCr) {
    // maxNumIter is different than maxIter above 
    pcout <<"performing the pressure ramp stage >>> "<<std::endl;
    bool hasNotConverged{true};
    plint iT{0}, totalNumIter{0};
    T newAvgEnF1{}, newAvgEnF2{}, oldAvgEnF1{1.}, oldAvgEnF2{1.};
    T volume = (T)(nx_*ny_*ny_);

    for (plint numRun = 0; numRun < totalNumRuns_; ++numRun) {
        if (numRun > 0) {
            setPressureBoundaryValues(inletRhoValues_[numRun], outletRhoValues_[numRun]);
        }
        
        hasNotConverged = true;
        oldAvgEnF1 = 1.0;
        oldAvgEnF2 = 1.0;
        iT = 0;
        while (hasNotConverged) {
            latticeFluidOne_.collideAndStream();
            latticeFluidTwo_.collideAndStream();

            if (totalNumIter % outputFreq == 0) {
                writeRhoVTK(outCounter_);
              //  writeVelocityComponentsDAT(numOut);   
                writeRhoDistributionDAT(outCounter_);
                ++outCounter_;                 
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

            if (iT >= maxRampIter) {
                hasNotConverged = false;
            }
            ++iT;
            ++totalNumIter;
        }
    }
}

void MultiPhaseRunOut::operator()(plint maxIter, plint maxRampIter, plint outputFreq, plint checkFreq, T convCr) {
    setShanChen();
    setUp();
    runEquilibrium(maxIter, outputFreq, checkFreq, convCr);
    runPressureRamp(maxRampIter, outputFreq, checkFreq, convCr);
}






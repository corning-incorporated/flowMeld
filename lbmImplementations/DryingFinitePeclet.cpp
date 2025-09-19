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
# include "../lbmDeclarations/DryingFinitePeclet.h"

void DryingFinitePeclet::setFluidsProperties(const CohesionParams<T> & cohesionParams, const FluidsParams<T> & fluidsParams) {
    
    T omegaF1{0}, omegaF2{0};

    spG_.resize(2, std::vector<T>(2));
    spG_.at(0).at(0) = cohesionParams.g00;
    // mutual values are set at 1 for the initial contact line equilibrium stage
    spG_.at(0).at(1) = 1.0;
    spG_.at(1).at(0) = 1.0;
    spG_.at(1).at(1) = cohesionParams.g11;
    // will be using this after contact line equilibrium 
    terminalG_ = cohesionParams.g01;

    omegaF1 = fluidsParams.omegaF1;
    omegaF2 = fluidsParams.omegaF2;
    gF1S_ = fluidsParams.gF1S;
    constOmegaValues_.assign({omegaF1, omegaF2});
    
}

// following methods will be run by setUp and do not need inputs 

void DryingFinitePeclet::setInletOutletDensities() {
  // for drying problem two inlet/outlet densities are added
  //    inlet = outlet and outlet = inlet - deltaRho*totalNumRuns 

  outletRhoValues_.push_back(rhoInitOutlet_);
  inletRhoValues_.push_back(rhoInitInlet_);

  T cosTheta  = std::abs(4.0*gF1S_/(spG_.at(0).at(1)*(rhoInitInlet_ - rhoNoFluid_)));
  T deltaRho = 6.0*sigma*cosTheta/minRadius_;

  if (totalNumRuns_ > 0) {
    outletRhoValues_.push_back(rhoInitOutlet_ - deltaRho*totalNumRuns_);
    inletRhoValues_.push_back(rhoInitInlet_);
  }
  else {
    pressureUpdate_ = false;
  }  

}


void DryingFinitePeclet::initPressureBC() {
    
    boundaryCondition_ = createLocalBoundaryCondition3D<T, MPDESCRIPTOR>();
    boundaryCondition_ -> addPressureBoundary0N(inlet_, latticeFluidTwo_);
    boundaryCondition_ -> addPressureBoundary0P(outlet_, latticeFluidTwo_);
    // add boundary conditions to the primary lattice as well 
    boundaryCondition_ -> addPressureBoundary0N(inlet_, latticeFluidOne_);
    boundaryCondition_ -> addPressureBoundary0P(outlet_, latticeFluidOne_);    
}

void DryingFinitePeclet::setPressureBoundaryValues(T rhoInlet, T rhoOutlet) {
    setBoundaryDensity(latticeFluidTwo_, inlet_, rhoInlet);
    setBoundaryDensity(latticeFluidTwo_, outlet_, rhoOutlet);

    setBoundaryDensity(latticeFluidOne_, inlet_, rhoNoFluid_);
    setBoundaryDensity(latticeFluidOne_, outlet_, rhoNoFluid_);

}


void DryingFinitePeclet::setShanChen(T gValue) {

    std::vector <MultiBlockLattice3D<T, MPDESCRIPTOR> *> blockLattices;
    plint processorLevel = 1;
    blockLattices.push_back(& latticeFluidTwo_);
    blockLattices.push_back(& latticeFluidOne_);

    spG_.at(0).at(1) = gValue;
    spG_.at(1).at(0) = gValue;

    integrateProcessingFunctional(new ShanChenMultiComponentProcessor3D <T, MPDESCRIPTOR> (spG_, constOmegaValues_),
         Box3D(0, nx_ - 1, 0, ny_ - 1, 0, nz_ - 1), blockLattices, processorLevel); 
}



void DryingFinitePeclet::runPressureRamp(plint maxRampIter, plint outputFreq, plint checkFreq, T convCr) {

    // maxNumIter is different than maxIter above 
    plint iT{0};
    T newAvgEnF1{}, newAvgEnF2{}, oldAvgEnF1{1.}, oldAvgEnF2{1.};
    T volume = (T)(nx_*ny_*ny_);

    
    setPressureBoundaryValues(inletRhoValues_[1], outletRhoValues_[1]);           

    for (iT = 0; iT < maxRampIter; ++iT) {
        latticeFluidOne_.collideAndStream();
        latticeFluidTwo_.collideAndStream();

        if (iT % outputFreq == 0) {
            writeRhoVTK(outCounter_);
            writeVelocityComponentsDAT(outCounter_);
            writeRhoDistributionDAT(outCounter_);
            ++outCounter_;   
            pcout  <<"generating output for the pressure flow and the drying "<<std::endl;             
        }

        if (iT % checkFreq == 0) {
            newAvgEnF1 = getStoredAverageDensity(latticeFluidOne_)*(volume);
            newAvgEnF2 = getStoredAverageDensity(latticeFluidTwo_)*(volume);
            oldAvgEnF1 = newAvgEnF1;
            oldAvgEnF2 = newAvgEnF2;
        }

    }
    
}

void DryingFinitePeclet::operator()(plint maxIter, plint maxRampIter, plint outputFreq, plint checkFreq, T convCr) {
    setShanChen(1.0);
    setUp();
    runEquilibrium(maxIter, outputFreq, checkFreq, convCr);
    setShanChen(terminalG_);
    if (pressureUpdate_) {
        runPressureRamp(maxRampIter, outputFreq, checkFreq, convCr);
    }
}



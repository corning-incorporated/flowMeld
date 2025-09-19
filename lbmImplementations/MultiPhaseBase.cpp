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

# include "../lbmDeclarations/MultiPhaseBase.h"

// helpers
void MultiPhaseBase::setDomainSize(const plint & nx, const plint & ny, const plint & nz) {
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
}

void MultiPhaseBase::setFileNames(const FileParams & fileParams) {
    geoFileName_ = fileParams.fNameIn;
    outputDir_ = fileParams.fNameOut;
    
}
void MultiPhaseBase::setDensities(const DensityParams<T> & densityParams) {
    rhoF1_ = densityParams.rhoF1;
    rhoF2_ = densityParams.rhoF2;
    rhoInitInlet_ = densityParams.rhoInitInlet;
    rhoInitOutlet_ = densityParams.rhoInitOutlet;
    rhoNoFluid_ = densityParams.rhoNoFluid;
}

void MultiPhaseBase::setPeriodicBCFlags(const PeriodicParams & periodicParams) {
    xPeriod_ = periodicParams.xPeriod;
    yPeriod_ = periodicParams.yPeriod;
    zPeriod_ = periodicParams.zPeriod;
    
    latticeFluidOne_.periodicity().toggle(0, xPeriod_);
    latticeFluidTwo_.periodicity().toggle(0, xPeriod_);

    latticeFluidOne_.periodicity().toggle(1, yPeriod_);
    latticeFluidTwo_.periodicity().toggle(1, yPeriod_);

    latticeFluidOne_.periodicity().toggle(2, zPeriod_);
    latticeFluidTwo_.periodicity().toggle(2, zPeriod_);

    // if periodic BC flags are off in y and z direction symmetry boundary
    // must be applied in subclasses that simulate the flow 
}

void MultiPhaseBase::setFluidsProperties(const FluidsParams<T> & fluidsParams) {
    T omegaF1{0}, omegaF2{0};
    omegaF1 = fluidsParams.omegaF1;
    omegaF2 = fluidsParams.omegaF2;
    constOmegaValues_.assign({omegaF1, omegaF2});
    gc_ = fluidsParams.gc;
    gF1S_ = fluidsParams.gF1S;
}


void MultiPhaseBase::setExternalForce(const ExternalForceParams<T> & externalForceParams) {
    forceDir_ = externalForceParams.forceDir;
    forceF1_ = externalForceParams.forceF1;
    forceF2_ = externalForceParams.forceF2;
}

void MultiPhaseBase::setShanChen() {        
    std::vector <MultiBlockLattice3D<T, MPDESCRIPTOR> *> blockLattices;
    plint processorLevel = 1;
    blockLattices.push_back(& latticeFluidTwo_);
    blockLattices.push_back(& latticeFluidOne_);

    integrateProcessingFunctional(new ShanChenMultiComponentProcessor3D <T, MPDESCRIPTOR> (gc_, constOmegaValues_),
         Box3D(0, nx_ - 1, 0, ny_ - 1, 0, nz_ - 1), blockLattices, processorLevel);      
}

void MultiPhaseBase::readGeometry() {
    Box3D slicebox(0,0, 0,ny_-1, 0,nz_-1);
    std::unique_ptr<MultiScalarField3D<int> > slice = generateMultiScalarField<int>(geometry_, slicebox);
    plb_ifstream geometryfile(geoFileName_.c_str());
    for (plint ix=0; ix<nx_; ++ix) {
        if (!geometryfile.is_open()) {
            pcout << "Error: could not open geometry file " << geoFileName_ << std::endl;
            exit(EXIT_FAILURE);
        }
        geometryfile >> *slice;
        copy(*slice, slice->getBoundingBox(), geometry_, Box3D(ix,ix, 0,ny_-1, 0,nz_-1));
    }
    geometryfile.close();
    {
        VtkImageOutput3D<T> vtkOut("porousMedium", 1.0);
        vtkOut.writeData<float>(*copyConvert<int, T> (geometry_, geometry_.getBoundingBox()), "tag", 1.0);
        std::unique_ptr<MultiScalarField3D<T> > floattags = copyConvert<int,T>(geometry_, geometry_.getBoundingBox());
        std::vector<T> isolevels;
        isolevels.push_back(0.5);
        typedef TriangleSet<T>::Triangle Triangle;
        std::vector<Triangle> triangles;
        Box3D domain = floattags->getBoundingBox().enlarge(-1);
        domain.x0++;
        domain.x1--;
        isoSurfaceMarchingCube(triangles, *floattags, isolevels, domain);
        TriangleSet<T> set(triangles);
        set.writeBinarySTL(outputDir_ + "porousMedium.stl");
    }        
}

void MultiPhaseBase::initBoundaryPlanes() {
    // note that boundary planes are defined internally
    // will be useful for setting up pressure boundary coinditions for drainage simulations
    // 0,0 and nx_ - 1 and nx_ - 1 will not work 
    inlet_ = Box3D(1, 2, 1, ny_-2, 1, nz_-2);
    outlet_ = Box3D(nx_-2, nx_-1, 1, ny_-2, 1, nz_-2);
    // side walls initialized directly: 
   // minYWall_ = Box3D(0, nx_-1, 0, 0, 0, nz_-1);
   // maxYWall_ = Box3D(0, nx_-1, ny_-1, ny_-1, 0, nz_-1);
   // minZWall_ = Box3D(0, nx_-1, 0, ny_-1, 0, 0);
   // maxZWall_ = Box3D(0, nx_-1, 0, ny_-1, nz_-1, nz_-1);
}

void MultiPhaseBase::defineLatticeDynamics() {
    // currently supports uniform wettability
    // 0: voids so no dynamics assigned
    // 1: surface nodes: bounce back with adhesion force
    // 2: interior solid nodes with bounce back or no dynamics for computational efficiency
    // interior solid nodes: no dynamics
    defineDynamics(latticeFluidOne_, geometry_, new NoDynamics<T,MPDESCRIPTOR>(), 2);
    defineDynamics(latticeFluidTwo_, geometry_, new NoDynamics<T,MPDESCRIPTOR>(), 2);

    //surface nodes with wettability: bounce back and adhesion 
    defineDynamics(latticeFluidOne_, geometry_, new BounceBack <T, MPDESCRIPTOR> (gF1S_), 1);
    defineDynamics(latticeFluidTwo_, geometry_, new BounceBack <T, MPDESCRIPTOR> (-1.0*gF1S_), 1);
}


void MultiPhaseBase::initializeLatticeDensities() {
    // 3: wetting fluid 0: non wetting fluid
    // for example: in contact angle measurements spreading fluid: -1, air: 0 
    Array<T, 3> zeroVelocity(0., 0., 0.);
    const plint nx = geometry_.getNx();
    const plint ny = geometry_.getNy();
    const plint nz = geometry_.getNz();

    for (plint iX = 0; iX < nx; iX++) {
        for (plint iY = 0; iY < ny; iY++) {
            for (plint iZ = 0; iZ < nz; iZ++) {
                plint tag = geometry_.get(iX, iY, iZ);
                // inert fluid (latticeFluidTwo_)
                if (tag == 0) {
                    initializeAtEquilibrium(latticeFluidTwo_, Box3D(iX,iX,iY,iY,iZ,iZ), rhoF2_, zeroVelocity);
                    initializeAtEquilibrium(latticeFluidOne_, Box3D(iX,iX,iY,iY,iZ,iZ), rhoNoFluid_, zeroVelocity);
                }
                // main fluid (latticeFluidOne_)
                // for drainage simulation tag == 3 is a fluid that invades the domain (it should be a non wetting fluid)
                else if (tag == 3) {
                    initializeAtEquilibrium(latticeFluidOne_, Box3D(iX,iX,iY,iY,iZ,iZ), rhoF1_, zeroVelocity);
                    initializeAtEquilibrium(latticeFluidTwo_, Box3D(iX,iX,iY,iY,iZ,iZ), rhoNoFluid_, zeroVelocity);
                }
            }
        }
    }
}

void MultiPhaseBase::addExternalForces() {
    Array<T, 3> forceF1; 
    Array<T, 3> forceF2;
    
    if (forceDir_ == "x") {
        forceF1 = Array<T, 3>(forceF1_, 0., 0.);
        forceF2 = Array<T, 3>(forceF2_, 0., 0.);
    }
    else if (forceDir_ == "y") {
        forceF1 = Array<T, 3>(0., forceF1_, 0.);
        forceF2 = Array<T, 3>(0., forceF2_, 0.); 
    }
    else if (forceDir_ == "z") {
        forceF1 = Array<T, 3>(0., 0., forceF1_);
        forceF2 = Array<T, 3>(0., 0., forceF2_); 
    }

    if (forceF1_ != 0.0) {
        setExternalVector(latticeFluidOne_, latticeFluidOne_.getBoundingBox(),
                    MPDESCRIPTOR<T>::ExternalField::forceBeginsAt, forceF1);
    }

    if (forceF2_ != 0.0) {
        setExternalVector(latticeFluidTwo_, latticeFluidTwo_.getBoundingBox(),
                MPDESCRIPTOR<T>::ExternalField::forceBeginsAt, forceF2);    
    }
}

void MultiPhaseBase::initializeLattices() {
    latticeFluidOne_.initialize();
    latticeFluidTwo_.initialize();
}

// main call() overriding operations
void MultiPhaseBase::setUp() {
    setShanChen();
    readGeometry();
    initBoundaryPlanes();
    defineLatticeDynamics();
    initializeLatticeDensities();
    addExternalForces();
    initializeLattices();
}

// checks convergence criteria
void MultiPhaseBase::operator()(plint checkFreq, plint outputFreq, plint maxIter, T convCr) {
    setUp();
    bool hasNotConverged{true};
    plint iT{0}, numOut{0};
    T newAvgEnF1{}, newAvgEnF2{}, oldAvgEnF1{1.}, oldAvgEnF2{1.};
    T relEF1{0}, relEF2{0};

    
    for (iT = 0; iT < maxIter; ++iT) {
        latticeFluidOne_.collideAndStream();
        latticeFluidTwo_.collideAndStream();
        
        if ((iT % checkFreq == 0) && (hasNotConverged)) {
            newAvgEnF1 = getStoredAverageDensity(latticeFluidOne_);
            newAvgEnF2 = getStoredAverageDensity(latticeFluidTwo_);
            relEF1 = std::fabs(oldAvgEnF1 - newAvgEnF1)*100.0/oldAvgEnF1/(T)checkFreq;
            relEF2 = std::fabs(oldAvgEnF2 - newAvgEnF2)*100.0/oldAvgEnF2/(T)checkFreq;
            pcout <<"the 1 energy value is "<<relEF1<<" cr: "<<convCr<<std::endl;
            pcout <<"the 2 energy value is "<<relEF2<<" cr: "<<convCr<<std::endl;
            if (simutils::hasConverged(oldAvgEnF1, oldAvgEnF2, newAvgEnF1, newAvgEnF2, (T) checkFreq, convCr)) {
                hasNotConverged = false;
                pcout <<"simulations converged at iteration "<<iT<<std::endl;
            }
            else {
                pcout <<"simulations has not converged yet at "<<iT<<std::endl;
            }
            oldAvgEnF1 = newAvgEnF1;
            oldAvgEnF2 = newAvgEnF2;
        }

      //  if ((iT % outputFreq == 0) && !(hasNotConverged)) {
        if (iT % outputFreq == 0) {
            pcout <<"generating output ... "<<iT<<std::endl;            
            writeRhoVTK(numOut);
        
        //    writeVelocityComponentsDAT(numOut);
            writeRhoDistributionDAT(numOut);
            ++numOut;
        }

    }
    writeSimulationDatFile();
}

// output methods 
void MultiPhaseBase::writeRhoVTK(plint it) {
    
    std::string rhoF1 = createFileName(outputDir_ + "f1_rho_step_", it, 6);
    std::string rhoF2 = createFileName(outputDir_ + "f2_rho_step_", it, 6);
    
    VtkImageOutput3D<T> vtkOutF1(rhoF1, 1.0);
    vtkOutF1.writeData<double> ((*computeDensity(latticeFluidOne_)), "density", 1.);

    VtkImageOutput3D<T> vtkOutF2(rhoF2, 1.0);
    vtkOutF2.writeData<double> ((*computeDensity(latticeFluidTwo_)), "density", 1.);
}


void MultiPhaseBase::writeVelocityComponentsDAT(plint it) {
    
    std::string fileSuffix = "_step_" + std::to_string(it) + ".dat";    
    std::string vxF1 = outputDir_ + "f1_vx" + fileSuffix;
    std::string vyF1 = outputDir_ + "f1_vy" + fileSuffix;
    std::string vzF1 = outputDir_ + "f1_vz" + fileSuffix;

    std::string vxF2 = outputDir_ + "f2_vx" + fileSuffix;
    std::string vyF2 = outputDir_ + "f2_vy" + fileSuffix;
    std::string vzF2 = outputDir_ + "f2_vz" + fileSuffix;

    plb_ofstream fOneX(vxF1.c_str());
    plb_ofstream fOneY(vyF1.c_str());
    plb_ofstream fOneZ(vzF1.c_str());

    plb_ofstream fTwoX(vxF2.c_str());
    plb_ofstream fTwoY(vyF2.c_str());
    plb_ofstream fTwoZ(vzF2.c_str());

    Box3D domain(0, nx_ - 1, 0, ny_ - 1, 0, nz_ - 1);
    fOneX<<*computeVelocityComponent(latticeFluidOne_, domain, 0)<<std::endl;
    fOneY<<*computeVelocityComponent(latticeFluidOne_, domain, 1)<<std::endl;
    fOneZ<<*computeVelocityComponent(latticeFluidOne_, domain, 2)<<std::endl;

    fTwoX<<*computeVelocityComponent(latticeFluidTwo_, domain, 0)<<std::endl;
    fTwoY<<*computeVelocityComponent(latticeFluidTwo_, domain, 1)<<std::endl;
    fTwoZ<<*computeVelocityComponent(latticeFluidTwo_, domain, 2)<<std::endl;    
}


void MultiPhaseBase::writeRhoDistributionDAT(plint it) {
    
    std::string fileSuffix = "_step_" + std::to_string(it) + ".dat";
    
    std::string rhoF1 = outputDir_ + "f1_rho_dist_" + fileSuffix;
    std::string rhoF2 = outputDir_ + "f2_rho_dist_" + fileSuffix;
    
    plb_ofstream rhoSF1(rhoF1.c_str());
    plb_ofstream rhoSF2(rhoF2.c_str());

    rhoSF1<<*computeDensity(latticeFluidOne_)<<std::endl;
    rhoSF2<<*computeDensity(latticeFluidTwo_)<<std::endl;

}

void MultiPhaseBase::addSimulationGeneralInfo(plb_ofstream & simInfo) {
    simInfo<<"f1_ads: "<<gF1S_<<std::endl;
    simInfo<<"diss_rho: "<<rhoNoFluid_<<std::endl;
}

void MultiPhaseBase::writeSimulationDatFile() {
    std::string simFile = outputDir_ + "simulation.dat";
    plb_ofstream simInfo(simFile.c_str());
    addSimulationGeneralInfo(simInfo);
    simInfo<<"delta_P: "<<(1./3.0)*(rhoInitInlet_ - rhoInitOutlet_)<<std::endl;
}




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

# include "../lbmDeclarations/SingleComponent.h"

void SingleComponent::setFileNames(const SingleCompFileParams & fileParams) {
    geoFileName_ = fileParams.geoName;
    outputDir_ = fileParams.outName;
    densityFileName_ = fileParams.rhoName;
}

void SingleComponent::setPeriodicBCFlags(const PeriodicParams & periodicParams) {
    xPeriod_ = periodicParams.xPeriod;
    yPeriod_ = periodicParams.yPeriod;
    zPeriod_ = periodicParams.zPeriod;
}

void SingleComponent::setFluidsProperties(const SingleCompFluidParams<T> & fluidParams) {
    omegaF_ = fluidParams.omegaF;
    gc_ = fluidParams.gc;
    gfs_ = fluidParams.gfs;
    nu_ = fluidParams.nu;
}

void SingleComponent::readGeometry() {
    Box3D slicebox(0,0, 0,ny_-1, 0,nz_-1);
    std::unique_ptr<MultiScalarField3D<int> > slice = generateMultiScalarField<int>(geometry_, slicebox);
    plb_ifstream geometryfile(geoFileName_.c_str());
    // note: in the previous version it was nx_ - 1
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



void SingleComponent::readDensity() {
    Box3D slicebox(0,0, 0,ny_-1, 0,nz_-1);
    std::unique_ptr<MultiScalarField3D<T> > slice = generateMultiScalarField<T>(density_, slicebox);
    plb_ifstream densityfile(densityFileName_.c_str());

    for (plint ix = 0; ix <nx_; ++ix) {
        if (!densityfile.is_open()) {
            pcout << "Error: could not open density file " << densityFileName_ << std::endl;
            exit(EXIT_FAILURE);
        }
        densityfile >> *slice;
        copy(*slice, slice->getBoundingBox(), density_, Box3D(ix,ix, 0,ny_-1, 0,nz_-1));
    }
}


void SingleComponent::setShanChenPotentialParameters(const T & rho_0, const T & psi_0) {
    rho_0_ = rho_0;
    psi_0_ = psi_0;
}

void SingleComponent::setExternalForce(const T & forceF) {
    forceF_ = forceF;
}

void SingleComponent::setUpShanChen() {
    plint processorLevel = 0;
    integrateProcessingFunctional(new ShanChenSingleComponentProcessor3D <T, MPDESCRIPTOR> (gc_, new interparticlePotential::PsiShanChen93<T>(rho_0_)), 
                lattice_.getBoundingBox(), lattice_, processorLevel);

}

void SingleComponent::defineLatticeDynamics() {
    defineDynamics(lattice_, geometry_, new NoDynamics<T, MPDESCRIPTOR>(), 2);
    defineDynamics(lattice_, geometry_, new BounceBack<T, MPDESCRIPTOR>(gfs_), 1);
}

void SingleComponent::initializeLatticeDensities() {
    // initialize densities using density_ which was read from a file before 
    Array<T, 3> zeroVelocity(0., 0., 0.);
    
    for (plint iX = 0; iX < nx_; iX++) {
        for (plint iY = 0; iY < ny_; iY++) {
            for (plint iZ = 0; iZ < nz_; iZ++) {
                plint tag = geometry_.get(iX, iY, iZ);
                if (tag == 0) {
                    T rho = density_.get(iX, iY, iZ);
                    initializeAtEquilibrium(lattice_, Box3D(iX, iX, iY, iY, iZ, iZ), rho, zeroVelocity);
                }
            }
        }
    }
}


void SingleComponent::addExternalForces() {
    if (forceF_ != 0.0) {
        setExternalVector(lattice_, lattice_.getBoundingBox(), 
                MPDESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(forceF_, 0., 0.));
    }
}

void SingleComponent::initializeLattices() {
    lattice_.initialize();
}

void SingleComponent::initializeSimulation() {
    setUpShanChen();
    readGeometry();
    readDensity();
    defineLatticeDynamics();
    initializeLatticeDensities();
    addExternalForces();
    initializeLattices();
}

void SingleComponent::operator()(plint checkFreq, plint outputFreq, plint maxIter, T convCr) {
    bool hasNotConverged{true};
    plint iT{0}, numOut{0};
    T newAvgEn{}, oldAvgEn{1.}; 
    T relEF{0};

    initializeSimulation();

    for (iT = 0; iT < maxIter; ++iT) {
        lattice_.collideAndStream();
        if ((iT % checkFreq) && (hasNotConverged)) {
            newAvgEn = getStoredAverageDensity(lattice_);
            if (simutils::hasConverged(oldAvgEn, newAvgEn, (T) checkFreq, convCr)) {
                hasNotConverged = false; 
                pcout <<"simulations converged at iteration "<<iT<<std::endl;
            }  
            else {
                pcout <<"simulations has not converged yet at "<<iT<<std::endl;
            }          
            oldAvgEn = newAvgEn;
        }

        if ((iT % outputFreq == 0) && !(hasNotConverged)) {
            pcout <<"generating output ... "<<iT <<std::endl;
            writeRhoVTK(numOut);
            writeRhoDistributionDAT(numOut);
            ++numOut;
        }
    }

}

/*** output methods ***/
void SingleComponent::writeRhoVTK(plint it) {
    std::string rhoF = createFileName(outputDir_ + "f_rho_step_", it, 6);
    VtkImageOutput3D<T> vtkOutF(rhoF, 1.0);
    vtkOutF.writeData<double> ((*computeDensity(lattice_)), "density", 1.);
}

void SingleComponent::writeVelocityComponentsDAT(plint it) {
    std::string fileSuffix = "_step_" + std::to_string(it) + ".dat";    
    std::string vxF = outputDir_ + "f_vx" + fileSuffix;
    std::string vyF = outputDir_ + "f_vy" + fileSuffix;
    std::string vzF = outputDir_ + "f_vz" + fileSuffix;

    plb_ofstream fOneX(vxF.c_str());
    plb_ofstream fOneY(vyF.c_str());
    plb_ofstream fOneZ(vzF.c_str());

    Box3D domain(0, nx_ - 1, 0, ny_ - 1, 0, nz_ - 1);
    fOneX<<*computeVelocityComponent(lattice_, domain, 0)<<std::endl;
    fOneY<<*computeVelocityComponent(lattice_, domain, 1)<<std::endl;
    fOneZ<<*computeVelocityComponent(lattice_, domain, 2)<<std::endl;
}

void SingleComponent::writeRhoDistributionDAT(plint it) {
    std::string fileSuffix = "_step_" + std::to_string(it) + ".dat";    
    std::string rhoF = outputDir_ + "f_rho_dist_" + fileSuffix;    
    plb_ofstream rhoSF(rhoF.c_str());
    rhoSF<<*computeDensity(lattice_)<<std::endl;
}


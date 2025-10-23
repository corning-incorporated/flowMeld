/************************************************************************************/
/*  Copyright 2025. Corning Incorporated. All rights reserved.                      */                                                                                     #
/**   This software may only be used in accordance with the identified license(s).    */  
/**                                                                                   */                                                                                      
/**   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      */ 
/**   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        */
/**   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL         */
/**   CORNING BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN      */
/**   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN               */
/**   CONNECTION WITH THE SOFTWARE OR THE USE OF THE SOFTWARE.                        */
/**************************************************************************************/
/**   Authors:                                                                        */
/**  Hamed Haddadi Staff Scientist                                                    */
/**                haddadigh@corning.com                                              */
/**  David Heine   Principal Scientist and Manager                                    */
/**                heinedr@corning.com                                                */
/**************************************************************************************/

# ifndef MULTIPHASEBASE_H_ 
# define MULTIPHASEBASE_H_ 

# include "../helpers/header.h"
# include "../helpers/mpParameterPacks.h"

class MultiPhaseBase {

    public:
        MultiPhaseBase(MultiBlockLattice3D<T, MPDESCRIPTOR> && latticeFluidOne,
                     MultiBlockLattice3D<T, MPDESCRIPTOR> && latticeFluidTwo, MultiScalarField3D<int> && geometry):
                        latticeFluidOne_{std::move(latticeFluidOne)},
                        latticeFluidTwo_{std::move(latticeFluidTwo)},
                        geometry_{std::move(geometry)}{};
        // class is not copyable
        MultiPhaseBase(const MultiPhaseBase&) = delete;
        MultiPhaseBase& operator=(const MultiPhaseBase&) = delete;

        // constructor and attribute flag helpers (called by client code)
        void setDomainSize(const plint &, const plint &, const plint &);
        void setFileNames(const FileParams &);
        void setPeriodicBCFlags(const PeriodicParams &); 
        void setExternalForce(const ExternalForceParams<T> &);
        // called by client code
        // computation methods
        void readGeometry();
        void defineLatticeDynamics();
        // if fluids are not loaded from a geometry file
        // invading and defending fluids are specified by initial coordinates
        void addExternalForces();
        void initializeLattices();
        void initializeLatticeDensities();
        // is used to initialize lattices from file, such as files for contact angle measurements
        // main call(): with and without checks for convergence 
        // output methods 
        void writeRhoVTK(plint);
        void writeVelocityComponentsDAT(plint);
        void writeRhoDistributionDAT(plint);
        void addSimulationGeneralInfo(plb_ofstream &) const;

        // virtual methods 
        // setShanChen can be different for drying problem 
        virtual void setShanChen();
        virtual void setFluidsProperties(const FluidsParams<T> &);
        virtual void setDensities(const DensityParams<T> &);
        virtual void setUp();
        virtual void initBoundaryPlanes();
        virtual void writeSimulationDatFile();
        virtual void operator()(plint, plint, plint, T);
        virtual ~MultiPhaseBase() = default;

    
    protected:
        // file names and directory paths
        std::string geoFileName_{}, outputDir_{}, forceDir_{};
        // domain size
        plint nx_{0}, ny_{0}, nz_{0};
        // fluid region bounds
        plint f1X1_{0}, f1X2_{0}, f1Y1_{0}, f1Y2_{0}, f1Z1_{0}, f1Z2_{0};
        plint f2X1_{0}, f2X2_{0}, f2Y1_{0}, f2Y2_{0}, f2Z1_{0}, f2Z2_{0};
        // periodic boundary conditions
        bool xPeriod_{true}, yPeriod_{true}, zPeriod_{true};
        // external forces
        T forceF1_{0}, forceF2_{0};
        // interaction parameters 
        T gF1S_{0}, gc_{0};
        // densities
        T rhoF1_{0}, rhoF2_{0}, rhoInitInlet_{0}, rhoInitOutlet_{0}, rhoNoFluid_{0}; 
        // relaxation times
        std::vector<T> constOmegaValues_;
        // core lattices
        MultiBlockLattice3D<T, MPDESCRIPTOR> latticeFluidOne_;
        MultiBlockLattice3D<T, MPDESCRIPTOR> latticeFluidTwo_;
        MultiScalarField3D<int> geometry_;
        // inlet and outlet boundaries
        Box3D inlet_, outlet_;
 
};

# endif 
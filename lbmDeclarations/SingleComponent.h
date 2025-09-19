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

# ifndef SINGLECOMPONENT_H_
# define SINGLECOMPONENT_H_ 

# include "../helpers/header.h"
# include "../helpers/mpParameterPacks.h"

class SingleComponent {

    public:
        SingleComponent(MultiBlockLattice3D<T, MPDESCRIPTOR> && lattice, 
                MultiScalarField3D<int> && geometry, 
                        MultiScalarField3D<T> && density):lattice_{std::move(lattice)}, 
                            geometry_{std::move(geometry)}, density_{std::move(density)}{
            nx_ = lattice_.getNx();
            ny_ = lattice_.getNy();
            nz_ = lattice_.getNz();
                            };
        
        void setShanChenPotentialParameters(const T &, const T &);
        void readDensity();
        void readGeometry();
        // overriden methods 
        // parameter set methods
        void setFileNames(const SingleCompFileParams &);
        void setPeriodicBCFlags(const PeriodicParams &);
        void setFluidsProperties(const SingleCompFluidParams<T> &);
        void setExternalForce(const T &);
        /* */
        void addExternalForces();
        void initializeLattices();
        void defineLatticeDynamics();
        void initializeLatticeDensities();
        void setUpShanChen();
        void initializeSimulation();
        void operator()(plint, plint, plint, T);
        // output methods
        void writeRhoVTK(plint); 
        void writeVelocityComponentsDAT(plint);
        void writeRhoDistributionDAT(plint); 
        // destructor 
        ~SingleComponent() = default;
        
    private:
        std::string geoFileName_{}, outputDir_{}, densityFileName_{};
        plint nx_{0}, ny_{0}, nz_{0};
        bool xPeriod_{true}, yPeriod_{true}, zPeriod_{true};       
        T omegaF_{}, gc_{}, gfs_{}, nu_{};
        T forceF_{};
        T rho_0_{0}, psi_0_{0};
        MultiBlockLattice3D<T, MPDESCRIPTOR> lattice_;
        // for perturbed/unpurturbed density field (in general: density initialization)
        MultiScalarField3D<int> geometry_;
        MultiScalarField3D<T> density_;
      //  interparticlePotential::PsiShanChen93<T> * psi_ptr_;     
};





# endif 
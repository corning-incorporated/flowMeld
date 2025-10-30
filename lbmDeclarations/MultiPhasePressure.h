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


// class for handling pressure boundary conditions 

# ifndef MULTIPHASEPRESSURE_H_ 
# define MULTIPHASEPRESSURE_H_ 

# include "./MultiPhaseBase.h"

class MultiPhasePressure: public MultiPhaseBase {
    
    public:
        MultiPhasePressure(MultiBlockLattice3D<T, MPDESCRIPTOR> && latticeFluidOne,
            MultiBlockLattice3D<T, MPDESCRIPTOR> && latticeFluidTwo,
            MultiScalarField3D<int> && geometry, plint numruns, T minradius):
            MultiPhaseBase(std::move(latticeFluidOne), std::move(latticeFluidTwo), std::move(geometry)),
            totalNumRuns_{numruns},
            minRadius_{minradius} {};
        
        MultiPhasePressure(const MultiPhasePressure &) = delete;
        MultiPhasePressure& operator=(const MultiPhasePressure &) = delete;

        MultiPhasePressure(MultiPhasePressure&& other) noexcept
            : MultiPhaseBase(std::move(other)), 
            totalNumRuns_(other.totalNumRuns_),
            minRadius_(other.minRadius_),
            inletRhoValues_(std::move(other.inletRhoValues_)),
            outletRhoValues_(std::move(other.outletRhoValues_)),
            pressureValues_(std::move(other.pressureValues_)),
            boundaryCondition_(other.boundaryCondition_)
            {
            other.boundaryCondition_ = nullptr;
            }

        MultiPhasePressure& operator=(MultiPhasePressure&& other) noexcept {
        if (this != &other) {
            MultiPhaseBase::operator=(std::move(other));
            totalNumRuns_ = other.totalNumRuns_;
            minRadius_ = other.minRadius_;
            inletRhoValues_ = std::move(other.inletRhoValues_);
            outletRhoValues_ = std::move(other.outletRhoValues_);
            pressureValues_ = std::move(other.pressureValues_);
            delete boundaryCondition_; // cleanup old pointer!
            boundaryCondition_ = other.boundaryCondition_;
            other.boundaryCondition_ = nullptr;
            }
            return *this;
            }

        
        // virtual methods
        virtual void initPressureBC();
        virtual void setInletOutletDensities(); 
        virtual void setPressureBoundaryValues(T, T);
        virtual void setUp();
        virtual void operator()(plint, plint, plint, T);
        virtual void writeSimulationDatFile();
        virtual ~MultiPhasePressure() {
            delete boundaryCondition_;
        }

    protected:
        // flow type: "drainage" or "runout" 
        plint totalNumRuns_{0};
        T minRadius_{};
        std::vector<T> inletRhoValues_; 
        std::vector<T> outletRhoValues_;
        std::vector<T> pressureValues_;
        OnLatticeBoundaryCondition3D<T, MPDESCRIPTOR>* boundaryCondition_{};

};


# endif 
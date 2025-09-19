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

# ifndef MULTIPHASERUNOUT_H_ 
# define MULTIPHASERUNOUT_H_ 

# include "./MultiPhasePressure.h"

class MultiPhaseRunOut : public MultiPhasePressure {
    public:
        MultiPhaseRunOut(MultiBlockLattice3D<T, MPDESCRIPTOR> && latticeFluidOne,
                     MultiBlockLattice3D<T, MPDESCRIPTOR> && latticeFluidTwo,
                      MultiScalarField3D<int> && geometry, plint numruns, T minradius):
                        MultiPhasePressure(std::move(latticeFluidOne), std::move(latticeFluidTwo),
                            std::move(geometry), numruns, minradius) {};
        
        MultiPhaseRunOut(const MultiPhaseRunOut &) = delete; 
        MultiPhaseRunOut& operator =(const MultiPhaseRunOut &) = delete;
        void runEquilibrium(plint, plint, plint, T);

        virtual void initPressureBC();
        virtual void setInletOutletDensities();
        virtual void setUp();
        virtual void setPressureBoundaryValues(T, T);
        virtual void runPressureRamp(plint, plint, plint, T); 
        virtual void operator()(plint, plint, plint, plint, T);

    
    protected:
        plint outCounter_{0};
        
};


# endif 
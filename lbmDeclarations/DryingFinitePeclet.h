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

// class for drying simulations
// main features are similar to the MultiPhaseBase class 
# ifndef DRYINGFINITEPECLET_H_ 
# define DRYINGFINITEPECLET_H_ 

# include "./MultiPhaseRunOut.h"

class DryingFinitePeclet: public MultiPhaseRunOut {
    // latticeFluidOne: water 
    // rhoF1_: water rhoF2_: air  
    // latticeFluidTwo: air 
    // pressure B.C
    // Note that Drying with Pressure BC is similar to drainage or runout  
    
    public:
        DryingFinitePeclet(MultiBlockLattice3D<T, MPDESCRIPTOR> && latticeFluidOne, 
            MultiBlockLattice3D<T, MPDESCRIPTOR> && latticeFluidTwo, 
                MultiScalarField3D<int> && geometry, plint deltapstrength, T minradius):MultiPhaseRunOut(std::move(latticeFluidOne), 
                    std::move(latticeFluidTwo), std::move(geometry), deltapstrength, minradius){};
        
        DryingFinitePeclet(const DryingFinitePeclet &) = delete;
        DryingFinitePeclet& operator=(const DryingFinitePeclet &) = delete;
        
        // all virtual functions 
        virtual void setFluidsProperties(const CohesionParams<T> &, const FluidsParams<T> &);
        virtual void setInletOutletDensities();
        virtual void setPressureBoundaryValues(T, T);
        virtual void initPressureBC();
        virtual void setShanChen(T);
       // virtual void setUp();
        virtual void runPressureRamp(plint, plint, plint, T); 
        virtual void operator()(plint, plint, plint, plint, T);
   
    protected:
        // vapor density
        bool pressureUpdate_{true};
        T terminalG_{0};
        std::vector<std::vector<T>> spG_;

};

# endif 
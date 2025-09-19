
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


# ifndef DRYINGRATECHANGE_H_ 
# define DRYINGRATECHANGE_H_ 

#include "./DryingFinitePeclet.h"

class DryingRateChange: public DryingFinitePeclet {
    public:
        DryingRateChange(MultiBlockLattice3D<T, MPDESCRIPTOR> && latticeFluidOne, 
            MultiBlockLattice3D<T, MPDESCRIPTOR> && latticeFluidTwo, 
                MultiScalarField3D<int> && geometry, plint deltapstrength, T minradius):
                    DryingFinitePeclet(std::move(latticeFluidOne), 
                        std::move(latticeFluidTwo), std::move(geometry), deltapstrength, minradius){};

        DryingRateChange(const DryingRateChange &) = delete;
        DryingRateChange& operator=(const DryingRateChange &) = delete;

        void setCohesionValues(const CohesionRangeParams<T> &, const std::string &);

        // methods
        virtual void setFluidsProperties(const CohesionRangeParams<T> &, const FluidsParams<T> &, const std::string &);
        virtual void setFluidsProperties(const CohesionRangeParams<T> &, const FluidsParamsChangingOmega<T> &, const std::string &);
        virtual void setShanChen(T, std::vector<T>);
        virtual void setUp();
        virtual void runPressureRamp(plint, plint, plint, T);
        virtual void operator()(plint, plint, plint, plint, T);
    

    protected:
        plint numGSteps_;
        plint gChangeStep_{0};
        std::vector<T> gValues_;
        // omega is a vector of [numGSteps_, 2] size and each row
        //  contains omegaf1, omegaf2
        std::vector<std::vector<T>> omegaValues_; 

};

# endif  
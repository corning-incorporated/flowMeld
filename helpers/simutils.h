
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


# ifndef SIMUTILS_H_
# define SIMUTILS_H_

# include <math.h>

namespace simutils {

template <typename U>
bool hasConverged(U oldAvgF1, U oldAvgF2, U newAvgF1, U newAvgF2, U checkFreq, U convCr) {
    bool converged{false};
    U relEF1 = std::fabs(oldAvgF1 - newAvgF1)*100.0/oldAvgF1/checkFreq;
    U relEF2 = std::fabs(oldAvgF2 - newAvgF2)*100.0/oldAvgF2/checkFreq;

    if (relEF1 < convCr && relEF2 < convCr) {
        converged = true;
    }
    else {
        converged = false;
    }
    
return converged;
}

/***** overloaded for single phase*/
template <typename U>
bool hasConverged(U oldAvg, U newAvg, U checkFreq, U convCr) {
    bool converged{false};
    U relEF = std::fabs(oldAvg - newAvg)*100.0/oldAvg/checkFreq;
    if (relEF < convCr) {
        converged = true;
    }
    else {
        converged = false;
    }
    return converged;
}
/*********************************/


}

# endif 
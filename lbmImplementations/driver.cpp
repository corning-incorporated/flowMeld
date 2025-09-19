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

// main driver file for multiphase flow

# include "../helpers/functionHeader.h"


int main(int argc, char ** argv) {
    
    plbInit(&argc, &argv);

    plint success{0};
    std::string modelName = argv[1];
    std::string xmlFileName = argv[2];

    global::timer("toma").restart();

    if (modelName == "multiphase") {
        success = runMultiPhaseMultiComponent(xmlFileName);
    }

    else if (modelName == "phasechange") {
        success = runMultiPhaseSingleComponent(xmlFileName);
    }

    if (success == 1) {
        T timeDuration = T();
        timeDuration = global::timer("toma").stop();
        pcout <<"LBM simulations successfully finished in "<<timeDuration<< std::endl; 
    }

   

    return 0;
}
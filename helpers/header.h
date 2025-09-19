
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

# ifndef HEADER_H_
# define HEADER_H_

# include "palabos2D.h"
# include "palabos2D.hh"   // include full template code
# include "palabos3D.h"
# include "palabos3D.hh"
# include "./simutils.h"

# include <iostream>
# include <vector>
# include <cmath>
# include <cstdlib>
# include <memory>
# include <algorithm>
# include <string>
# include <fstream>
# include <sstream>
# include <iterator>
# include <map>

using namespace plb;
typedef double T;
const double sigma = 0.15;

# define MPDESCRIPTOR descriptors::ForcedShanChenD3Q19Descriptor
# define SPDESCRIPTOR descriptors::D3Q19Descriptor 
# endif 
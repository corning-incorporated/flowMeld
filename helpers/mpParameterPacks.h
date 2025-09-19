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

// useful parameter packs for multiphase flow simulations
# ifndef MPPARAMETERPACKS_H_ 
# define MPPARAMETERPACKS_H_

# include "palabos3D.h"
# include "palabos3D.hh"

struct DomainParams {
    plint nx,ny,nz;
    DomainParams() = default; 
    DomainParams(plint x, plint y, plint z):nx{x},ny{y},nz{z}{};
};

struct PeriodicParams {
    bool xPeriod{true}, yPeriod{true}, zPeriod{true};
    PeriodicParams() = default;
    PeriodicParams(bool xp, bool yp, bool zp):xPeriod{xp},yPeriod{yp},zPeriod{zp}{};    
};

struct FileParams
{
    std::string fNameIn{}, fNameOut{};
    FileParams() = default;
    FileParams(std::string fnameIn, std::string fnameOut):fNameIn{fnameIn},fNameOut{fnameOut}{};
};


struct SingleCompFileParams
{
    std::string geoName{}, outName{}, rhoName{};
    SingleCompFileParams() = default;
    SingleCompFileParams(std::string geo, std::string out, std::string rho):geoName{geo}, outName{out}, rhoName{rho}{};
};

template <typename U>
struct FluidsParams {
    // note than Fluid 1 is the wetting fluid
    // currently only one type of solid surface is added
    U omegaF1{}, omegaF2{}, gc{}, gF1S{};
    FluidsParams() = default;
    FluidsParams(U omegaf1, U omegaf2, U g, U gF1S):omegaF1{omegaf1},omegaF2{omegaf2},
            gc{g}, gF1S{gF1S}{};
};

template <typename U>
struct ExternalForceParams {
    U forceF1{0.}, forceF2{0.};
    std::string forceDir{};
    ExternalForceParams() = default;
    ExternalForceParams(U forcef1, U forcef2, std::string forcedir):forceF1{forcef1},
                    forceF2{forcef2}, forceDir{forcedir}{};
};

// for changing viscosity during simulations 
template <typename U>
struct FluidsParamsChangingOmega {
    U omegaMinF1{}, omegaMaxF1{}, omegaMinF2{}, omegaMaxF2{};
    FluidsParamsChangingOmega() = default;
    FluidsParamsChangingOmega(U omegaminf1, U omegamaxf1, U omegaminf2, U omegamaxf2):omegaMinF1{omegaminf1},
             omegaMaxF1{omegamaxf1}, omegaMinF2{omegaminf2}, omegaMaxF2{omegamaxf2}{};
};


template <typename U>
struct SingleCompFluidParams {
    U omegaF{}, gc{}, gfs{}, nu{};
    SingleCompFluidParams() = default;
    SingleCompFluidParams(U omega, U g, U gs, U nu):omegaF{omega}, gc{g}, gfs{gs}, nu{nu}{};
};

// if cohesion is a vector
template <typename U>
struct CohesionParams {
    U g00{}, g01{}, g11{};
    CohesionParams() = default;
    CohesionParams(U g0, U g1, U g2):g00{g0}, g01{g1}, g11{g2}{};
};

template <typename U>
struct CohesionRangeParams {
    U g00{}, g11{}, gMin{}, gMax{};
    plint gSteps{0};
    CohesionRangeParams() = default;
    CohesionRangeParams(U g0, U g1, U gmin, U gmax, plint gsteps): g00{g0}, g11{g1}, gMin{gmin}, gMax{gmax}, gSteps{gsteps}{};
};


template <typename U>
struct DensityParams {
U rhoF1{}, rhoF2{}, rhoInitInlet{}, rhoInitOutlet{}, rhoNoFluid{};
DensityParams() = default; 
DensityParams(U rhof1, U rhof2, U rhoinitinlet, U rhoinitoutlet, U rhonofluid):rhoF1{rhof1},
    rhoF2{rhof2}, rhoInitInlet{rhoinitinlet}, rhoInitOutlet{rhoinitoutlet}, rhoNoFluid{rhonofluid}
     {};
};

template <typename U>
struct DryingDensityParams {
    U rhoF1{}, rhoF2{}, rhoNoFluid{};
    DryingDensityParams() = default; 
    DryingDensityParams(U rhof1, U rhof2, U rhonofluid):rhoF1{rhof1}, rhoF2{rhof2}, rhoNoFluid{rhonofluid}{};
};


struct CoordinateParams {
    plint fX1{0}, fX2{0}, fY1{0}, fY2{0}, fZ1{0}, fZ2{0};
    CoordinateParams() = default;
    CoordinateParams(plint fx1, plint fx2, plint fy1, plint fy2, plint fz1, plint fz2):
        fX1{fx1}, fX2{fx2}, fY1{fy1}, fY2{fy2}, fZ1{fz1}, fZ2{fz2}{};
};


# endif 
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


// runs the single component model (phase change)
# include "../lbmDeclarations/SingleComponent.h"
# include "../helpers/mpParameterPacks.h"

int runMultiPhaseSingleComponent(const std::string & xmlFileName) {

    std::string geomFileName{}, outputDir{}, densityFileName{};
    plint maxIter{0}, outputFreq{0};
    plint nx{0}, ny{0}, nz{0};
    plint convCheckFreq{0};
    bool xPeriod{true}, yPeriod{true}, zPeriod{true};
    T omegaF{}, gc{}, gfs{}, nu{};
    T rho0{}, psi0{};
    T forceF{};
    T convCr{};

    try {
        XMLreader document(xmlFileName);
        document["filenames"]["microstructure"].read(geomFileName);
        document["filenames"]["output_directory"].read(outputDir);  
        document["filenames"]["density_input"].read(densityFileName);
        // domain
        document["domain"]["resolution"]["x"].read(nx);
        document["domain"]["resolution"]["y"].read(ny);
        document["domain"]["resolution"]["z"].read(nz);
        // periodic boundary flags
        document["domain"]["periodic_bc"]["x"].read(xPeriod);
        document["domain"]["periodic_bc"]["x"].read(yPeriod);
        document["domain"]["periodic_bc"]["y"].read(zPeriod);
        // phase information 
        document["phase"]["relaxation_omega"].read(omegaF);
        document["phase"]["cohesion_gc"].read(gc);
        document["phase"]["adhesion_gfs"].read(gfs);
        document["phase"]["rho_0"].read(rho0);
        document["phase"]["psi_0"].read(psi0);
        document["phase"]["external_force"].read(forceF);
        //simulation info 
        document["simulations"]["max_iterations"].read(maxIter);
        document["simulations"]["output_frequency"].read(outputFreq);
        document["simulations"]["converge_check_frequency"].read(convCheckFreq);
        document["simulations"]["converge_criterion"].read(convCr);

    } catch (PlbIOException & exception) {
        pcout << exception.what() << std::endl; 
        return -1;
    }


    nu = ((T)1/omegaF - (T)0.5)/MPDESCRIPTOR<T>::invCs2;
    if (outputDir[outputDir.size() - 1] != '/') {
        outputDir += '/';
    }  

    global::directories().setOutputDir(outputDir);

    SingleCompFileParams fileParams(geomFileName, outputDir, densityFileName);
    PeriodicParams periodParams(xPeriod, yPeriod, zPeriod);
    SingleCompFluidParams<T> fluidParams(omegaF, gc, gfs, nu);

    // define the phase lattice
    MultiBlockLattice3D<T, MPDESCRIPTOR> lattice(nx, ny, nz, new ExternalMomentRegularizedBGKdynamics< T, MPDESCRIPTOR> (omegaF));
    MultiScalarField3D<int> geometry(nx, ny, nz);
    MultiScalarField3D<T> density(nx, ny, nz, T(0.0));

    SingleComponent singleComp(std::move(lattice), std::move(geometry), std::move(density));
    singleComp.setFileNames(fileParams);
    singleComp.setPeriodicBCFlags(periodParams);
    singleComp.setFluidsProperties(fluidParams);
    singleComp.setExternalForce(forceF);
    singleComp.setShanChenPotentialParameters(rho0, psi0);

    singleComp(convCheckFreq, outputFreq, maxIter, convCr);

    return 1;

}
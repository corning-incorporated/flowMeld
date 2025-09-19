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

// runs the multiphase model (3D imbibition/drainage model)
# include "../lbmDeclarations/MultiPhaseBase.h"
# include "../lbmDeclarations/MultiPhasePressure.h"
# include "../lbmDeclarations/MultiPhaseRunOut.h"
# include "../lbmDeclarations/DryingFinitePeclet.h"
# include "../lbmDeclarations/DryingRateChange.h"
# include "../helpers/mpParameterPacks.h"

int runMultiPhaseMultiComponent(const std::string & xmlFileName) {

    // declaration of variables (local to driver)
    // variables will be passed to other classes
    std::string tomaFileName{}, outputDir{}, simType{}, changeType{}, forceDir;
    plint maxIter{0}, maxRampIter{0}, outputFreq{0};
    plint nx{0}, ny{0}, nz{0};
    // applicable to pressure BC
    plint totalNumRuns{0}, minRadius{0};
    plint convCheckFreq{0};
    plint gsteps{0}, changeStep{0};
    bool xPeriod{true}, yPeriod{false}, zPeriod{false}, omegaChange{false};
    T omegaF1{}, omegaF2{}, gc{0.0}, gF1S{}, g00{0}, g01{0}, g11{0}, omegaMinF1{0.}, omegaMaxF1{0.}, omegaMinF2{0}, omegaMaxF2{0.};
    T gmin{0}, gmax{0};
    T rhoF1{}, rhoF2{}, rhoInitInlet{0.0}, rhoInitOutlet{0.0}, rhoNoFluid{};
    T forceF1{}, forceF2{};
    T convCr{};

    try {
        XMLreader document(xmlFileName);
        // input & output file names
        document["filenames"]["microstructure"].read(tomaFileName);
        document["filenames"]["output_directory"].read(outputDir);
        // domain information
        document["domain"]["resolution"]["x"].read(nx);
        document["domain"]["resolution"]["y"].read(ny);
        document["domain"]["resolution"]["z"].read(nz);
        // periodic boundary flags
        document["domain"]["periodic_bc"]["x"].read(xPeriod);
        document["domain"]["periodic_bc"]["x"].read(yPeriod);
        document["domain"]["periodic_bc"]["y"].read(zPeriod);
        // only for pressure BC
        document["flow"]["type"].read(simType);                
        document["flow"]["number_of_pressure_steps"].read(totalNumRuns);
        document["flow"]["min_throat_radius"].read(minRadius);   

        // fluid property parameters
        document["fluids"]["gc"].read(gc);
        // g11, g01, g11 should be zero for all models except drying
        // gChangeType: 'range' or 'step'
        document["fluids"]["change_type"].read(changeType);
        // constant rate drying
        document["fluids"]["g00"].read(g00);
        document["fluids"]["g01"].read(g01);
        document["fluids"]["g11"].read(g11);
        // gmin, gmax, gsteps: for drying rate change with increments 
        document["fluids"]["gmin"].read(gmin);
        document["fluids"]["gmax"].read(gmax);
        // currently: uniform wettability
        // how it is used: MultiPhaseBase::defineLatticeDynamics
        document["fluids"]["f1_fluid_surface_adhesion"].read(gF1S);
        document["fluids"]["omega_f1"].read(omegaF1);
        document["fluids"]["omega_f2"].read(omegaF2);
        // parameters related to change in omega 
        document["fluids"]["omega_change"].read(omegaChange);
        document["fluids"]["omega_min_f1"].read(omegaMinF1);
        document["fluids"]["omega_max_f1"].read(omegaMaxF1);
        document["fluids"]["omega_min_f2"].read(omegaMinF2);
        document["fluids"]["omega_max_f2"].read(omegaMaxF2);
        document["fluids"]["num_steps"].read(gsteps);
        document["fluids"]["change_step"].read(changeStep);
        document["fluids"]["density_f1"].read(rhoF1);
        document["fluids"]["density_f2"].read(rhoF2);
        document["fluids"]["density_no_fluid"].read(rhoNoFluid);
        document["fluids"]["force_f1"].read(forceF1);
        document["fluids"]["force_f2"].read(forceF2);
        document["fluids"]["force_direction"].read(forceDir);
        // invading-defending fluid initial coordinates
        // run and simulation parameters 
        document["simulations"]["max_iterations"].read(maxIter);
        document["simulations"]["max_pressure_iterations"].read(maxRampIter);
        document["simulations"]["output_frequency"].read(outputFreq);
        document["simulations"]["converge_check_frequency"].read(convCheckFreq);
        document["simulations"]["converge_criterion"].read(convCr);

    } catch (PlbIOException & exception) {
        pcout << exception.what() << std::endl; 
        return -1;
    }

    // set global variables
    // compute relaxation times 
   // nuF1 = ((T)1/omegaF1 - (T)0.5)/MPDESCRIPTOR<T>::invCs2;
   // nuF2 = ((T)1/omegaF2 - (T)0.5)/MPDESCRIPTOR<T>::invCs2;
    // output directory
    if (outputDir[outputDir.size() - 1] != '/') {
        outputDir += '/';
    }
    global::directories().setOutputDir(outputDir);
    if (rhoInitInlet == 0) {
        rhoInitInlet = rhoF1;
    }
    if (rhoInitOutlet == 0) {
        rhoInitOutlet = rhoF2;
    }

    // populate parameter packs
    FileParams fileParams(tomaFileName, outputDir);
    DomainParams domainParams(nx, ny, nz);
    PeriodicParams periodicParams(xPeriod, yPeriod, zPeriod);
    DensityParams<T> densityParams(rhoF1, rhoF2, rhoInitInlet, rhoInitOutlet, rhoNoFluid);
    FluidsParams<T> fluidsParams(omegaF1, omegaF2, gc, gF1S);
    CohesionParams<T> cohesionParams(g00, g01, g11);
    ExternalForceParams<T> externalForceParams(forceF1, forceF2, forceDir);
    CohesionRangeParams<T> cohesionRangeParams;
    FluidsParamsChangingOmega<T> fluidsRangeParams; 

    
    if (changeType == "range") {
        cohesionRangeParams = CohesionRangeParams<T>(g00, g11, gmin, gmax, gsteps);
    }
    else if (changeType == "step") {
        cohesionRangeParams = CohesionRangeParams<T>(g00, g11, gmin, gmax, changeStep);
    }

    if (omegaChange) {
        fluidsRangeParams = FluidsParamsChangingOmega<T>(omegaMinF1, omegaMaxF1, omegaMinF2, omegaMaxF2);
    }

    // define MultiBlock lattices and pass them to the class 
    MultiBlockLattice3D < T, MPDESCRIPTOR > latticeFluidOne(nx, ny, nz,
        new ExternalMomentRegularizedBGKdynamics < T, MPDESCRIPTOR > (omegaF1));

    MultiBlockLattice3D < T, MPDESCRIPTOR > latticeFluidTwo(nx, ny, nz, 
        new ExternalMomentRegularizedBGKdynamics < T, MPDESCRIPTOR > (omegaF2));
    
    MultiScalarField3D<int> geometry(nx, ny, nz);
    
    if (simType == "drainage") {
        MultiPhasePressure multiPressure(std::move(latticeFluidOne),
                    std::move(latticeFluidTwo), std::move(geometry), totalNumRuns, minRadius);       
        multiPressure.setDomainSize(nx, ny, nz);
        multiPressure.setFileNames(fileParams);
        multiPressure.setDensities(densityParams);
        multiPressure.setPeriodicBCFlags(periodicParams);
        multiPressure.setFluidsProperties(fluidsParams);
        multiPressure.setExternalForce(externalForceParams);
        multiPressure(maxIter, convCheckFreq, outputFreq, convCr);
    }

    else if (simType == "runout") {
        MultiPhaseRunOut multiRunOut(std::move(latticeFluidOne),
                    std::move(latticeFluidTwo), std::move(geometry), totalNumRuns, minRadius);
        multiRunOut.setDomainSize(nx, ny, nz);
        multiRunOut.setFileNames(fileParams);
        multiRunOut.setDensities(densityParams);
        multiRunOut.setPeriodicBCFlags(periodicParams);
        multiRunOut.setFluidsProperties(fluidsParams);
        multiRunOut.setExternalForce(externalForceParams);
        multiRunOut(maxIter, maxRampIter, outputFreq, convCheckFreq, convCr);

    }

    else if (simType == "imbibition") {
        MultiPhaseBase multiPhase(std::move(latticeFluidOne), std::move(latticeFluidTwo), std::move(geometry));
        multiPhase.setDomainSize(nx, ny, nz);
        multiPhase.setFileNames(fileParams);
        multiPhase.setDensities(densityParams);
        multiPhase.setPeriodicBCFlags(periodicParams);
        multiPhase.setFluidsProperties(fluidsParams);
        multiPhase.setExternalForce(externalForceParams);
        multiPhase(convCheckFreq, outputFreq, maxIter, convCr);
    }

    else if (simType == "drying") {
        DryingFinitePeclet drying(std::move(latticeFluidOne), std::move(latticeFluidTwo), std::move(geometry),
                        totalNumRuns, minRadius);
        drying.setDomainSize(nx, ny, nz);
        drying.setFileNames(fileParams);
        drying.setDensities(densityParams);
        drying.setPeriodicBCFlags(periodicParams);
        drying.setFluidsProperties(cohesionParams, fluidsParams);
        drying.setExternalForce(externalForceParams);
        drying(maxIter, maxRampIter, outputFreq, convCheckFreq, convCr);
    }

    else if (simType == "drying-rate") {
        DryingRateChange dryRate(std::move(latticeFluidOne), std::move(latticeFluidTwo), std::move(geometry),
                                totalNumRuns, minRadius);
        dryRate.setDomainSize(nx, ny, nz);
        dryRate.setFileNames(fileParams);
        dryRate.setDensities(densityParams);
        dryRate.setPeriodicBCFlags(periodicParams);
        dryRate.setPeriodicBCFlags(periodicParams);

        if (omegaChange) {
            dryRate.setFluidsProperties(cohesionRangeParams, fluidsRangeParams, changeType);
        }
        else {
            dryRate.setFluidsProperties(cohesionRangeParams, fluidsParams, changeType);
        }

        dryRate.setExternalForce(externalForceParams);
        dryRate(maxIter, maxRampIter, outputFreq, convCheckFreq, convCr);
    }

    return 1;

}
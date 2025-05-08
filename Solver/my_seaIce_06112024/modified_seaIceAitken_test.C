#include "fvCFD.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "EigenMatrix.H"

// NOTE NOTE: attempt to get isoAdvector to work
#include "isoAdvection.H"
//************************************

#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "Newtonian.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    //AITKEN PARAMETERS
    // scalar initial_relax_fact = 1.0;
    // scalar relax_fact;
    volVectorField U_old = U;
    // scalarField R(h.size()*2,0);
    // scalarField Ri(h.size()*2,0);

    while (runTime.run())
    {

        //NOTE NOTE: I changed it back to being based on U courante number, not alpha courante number
        //ALSO, if you would like to ensure that it outputs neatly to specific time steps instead of arbitrary
        //time steps then you should change "writeControl" to "adjustableRunTime" in system/ControlDict
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

       runTime++;

    Info<< "Time = " << runTime.timeName() << nl << endl;

    phi = fvc::interpolate(U) & mesh.Sf();

        #include "alphaControls.H"
        #include "alphaEqnSubCycle.H"

    interface.correct();

    h_Phi = fvc::interpolate(h)*phi;
    label iter_counter = 0;

    #include "hEqn.H"

    while (pimple.loop())
    {            
        // phi = fvc::interpolate(U) & mesh.Sf();
        // h_Phi = fvc::interpolate(h)*phi;

        // #include "hEqn.H" 

        iter_counter += 1;

        U_old = U;

	    #include "UEqn.H"

        //TODO Parrallelize. As master pull in U from all the domains

	    //==================================================================
	    //
	    // AITKEN DYNAMIC RELAXATION
            // TODO: Need to still pull all the information from all the CPUS
	    //
	    //==================================================================

        /// volVectorField Res = U - U_old;

        /// forAll(Res,i)
        /// {
        ///     R[i] = Res[i][0];
        ///     R[i+h.size()] = Res[i][1];
        /// }

        /// if (iter_counter ==1)
        /// {
        ///     Ri = R;
        ///     relax_fact = initial_relax_fact;
        /// }
        /// else
        /// {
        ///     scalar dot = 0.0;
        ///     scalar norm_square = 0.0;
        ///     forAll(R,i)
        ///     {
        ///         dot += Ri[i]*(R[i]-Ri[i]);
        ///         norm_square += sqr(Ri[i]-R[i]);
        ///     }
        ///     relax_fact = -relax_fact*dot/norm_square;
        ///     Ri = R;
        /// }
        /// U = relax_fact*U + (1-relax_fact)*U_old;
        /// Info << "Relax fact: " << relax_fact << endl;
        /// U.correctBoundaryConditions();

    } 

    runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;
}

Info<< "End\n" << endl;

return 0;
}

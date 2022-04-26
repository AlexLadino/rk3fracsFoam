/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
Application
    pisoFoam
Description
    Transient solver for incompressible, turbulent flow, using the PISO
    algorithm.
    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    turbulence->validate();
    
    #include "pEqn.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

            Uold=U; Uc=U;
            phi = (fvc::interpolate(U) & mesh.Sf());
            dU = runTime.deltaT()*(fvc::laplacian(turbulence->nuEff(),U) - fvc::div(phi,U) );
            Uc = Uc + (1.0/6.0)*dU; U = Uold + 0.5*dU;
            #include "pCorrection.H"
            phi = (fvc::interpolate(U) & mesh.Sf());
            dU = runTime.deltaT()*(fvc::laplacian(turbulence->nuEff(),U) - fvc::div(phi,U) );
            Uc = Uc + (2.0/3.0)*dU; U = Uold + 2.0*(Uold - U) + 2.0*dU;
            #include "pCorrection.H"
            phi = (fvc::interpolate(U) & mesh.Sf());
            dU = runTime.deltaT()*(fvc::laplacian(turbulence->nuEff(),U) - fvc::div(phi,U) );
            Uc = Uc + (1.0/6.0)*dU; U = Uc;
            #include "pCorrection.H"
            phi = (fvc::interpolate(U) & mesh.Sf());
            turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
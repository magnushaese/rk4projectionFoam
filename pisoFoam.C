/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    pisoFoam

Description
    Transient solver for incompressible flow.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pisoControl piso(mesh);

#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
	fvScalarMatrix pEqn
	( 
		fvm::laplacian(p)
	);
	pEqn.setReference(pRefCell, pRefValue);
  */                                         
scalarField betaRK=scalarField(4);
betaRK[0]=0.166666667;
betaRK[1]=0.333333333;
betaRK[2]=0.333333333;
betaRK[3]=0.166666667;

Info << "Time stepping coeficcients: ";
forAll(betaRK,bRK) {
	Info << "betaRK[" << bRK << "] = " << betaRK[bRK] << " ";
}
Info << endl;

scalarField alphaRK=scalarField(4);
alphaRK[0]=0.5;
alphaRK[1]=0.5;
alphaRK[2]=1.0;
alphaRK[3]=0.0;

Info << "aux coeficcients: ";
forAll(alphaRK,aRK) {
	Info << "alphaRK[" << aRK << "] = " << alphaRK[aRK] << " ";
}
Info << endl;


    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "CourantNo.H"
        runTime++;
	fvScalarMatrix pEqn
	( 
		fvm::laplacian(p)
	);
	pEqn.setReference(pRefCell, pRefValue);
 
	Uold = U; Uc = U; dt=runTime.deltaT();
	scalar ii=0;

	forAll(betaRK,bRK) 
	{
        	phi = (fvc::interpolate(U) & mesh.Sf());
        	dU = dt*(fvc::laplacian(turbulence->nuEff(), U) - fvc::div(phi,U));
       		Uc = Uc +betaRK[ii]*dU; U = Uold + alphaRK[ii]*dU;
		if(ii==3) U=Uc;
        	#include "pressureCorrection.H"
	ii++;
	}

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

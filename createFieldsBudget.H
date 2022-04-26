/*
UPrime = U - UMean;
pPrime = p - pMean;
// TKE instantanea
k = 0.5*magSqr(UPrime);
kMean = 0.5*tr(UPrime2Mean);
// TKE Budget terms
// convection (hacer fieldAverage convectionMean)
convection = -UMean & fvc::grad(k);
// Production (se calcula en fase de postProcesamiento)
productionMean = -UPrime2Mean && fvc::grad(UMean);
// viscous dissipation (hacer fieldAverage viscousDissipationMean)
viscousDissipation = -turbulence->nuEff()*tr(grad(UPrime)&grad(UPrime).T());
// transport Viscous Diffusion (en fase de postProcesamiento)
transportViscousDiffusionMean = turbulence->nuEff()*fvc::laplacian(kMean); 
// Transport Velocity-Pressure Gradient 
// 1. hacer fieldAverage transportVelocityPressureMean
// 2. en fase de postProcesamiento hacer -div(transportVelocityPressureMean);
transportVelocityPressure = pPrime*UPrime;
// Transport by turbulence
// 1. hacer fieldAverage transportTurbulenceMean
// 2. en fase de postProcesamiento hacer -div(transportTurbulenceMean);
transportTurbulence = k*UPrime;
*/

Info<< "Reading field UMean\n" << endl;
volVectorField UMean
(
    IOobject
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field pMean\n" << endl;
volScalarField pMean
(
    IOobject
    (
        "pMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field UPrime2Mean\n" << endl;
volSymmTensorField UPrime2Mean
(
    IOobject
    (
        "UPrime2Mean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volVectorField UPrime
(
    IOobject
    (
        "UPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,   
      dimensionedVector
      (
           "UPrime", 
           dimVelocity,
           vector(0.0,0.0,0.0)
      )
);

/*
volTensorField gradUPrime
(
    IOobject
    (
        "gradUPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,   
      dimensionedVector
      (
           "gradUPrime", 
           dimVelocity/dimLength,
           tensor(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
      )
);
*/

volTensorField gradUPrime(fvc::grad(UPrime));

volScalarField pPrime
(
    IOobject
    (
        "pPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    
      dimensionedScalar
      (
           "pPrime", 
           dimVelocity*dimVelocity,
           0.0
      )
);

volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    
      dimensionedScalar
      (
           "k", 
           dimVelocity*dimVelocity,
           0.0
      )
);

/*
volScalarField kMean
(
    IOobject
    (
        "kMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,   
      dimensionedScalar
      (
           "kMean", 
           dimVelocity*dimVelocity,
           0.0
      )
);
*/ 

volScalarField convection
(
    IOobject
    (
        "convection",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
      dimensionedScalar
      (
            "convection", 
            dimVelocity*dimVelocity/dimTime,
            0.0
      )
);

volScalarField viscousDissipation
(
    IOobject
    (
        "viscousDissipation",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
      dimensionedScalar
      (
            "viscousDissipation", 
            dimVelocity*dimVelocity/dimTime,
            0.0
      )
);

volVectorField transportVelocityPressureField
(
    IOobject
    (
        "transportVelocityPressureField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,   
      dimensionedVector
      (
           "transportVelocityPressureField", 
           dimVelocity*dimVelocity*dimLength/dimTime,
           vector(0.0,0.0,0.0)
      )
);

volVectorField transportTurbulenceField
(
    IOobject
    (
        "transportTurbulenceField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,   
      dimensionedVector
      (
           "transportTurbulenceField", 
           dimVelocity*dimVelocity*dimVelocity,
           vector(0.0,0.0,0.0)
      )
);


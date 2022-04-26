volScalarField e
(
    IOobject
    (
        "e",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "e", 
           dimVelocity*dimVelocity,
           0.0
      )
);

// unsteady term of kinetic energy
volScalarField dedt
(
    IOobject
    (
        "dedt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "dedt", 
           dimVelocity*dimVelocity/dimTime,
           0.0
      )
);

// kinetic energy flux
volScalarField deuidxi
(
    IOobject
    (
        "deuidxi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "deuidxi", 
           dimVelocity*dimVelocity/dimTime,
           0.0
      )
);

// pressure flux
volScalarField dpuidxi
(
    IOobject
    (
        "dpuidxi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "dpuidxi", 
           dimVelocity*dimVelocity/dimTime,
           0.0
      )
);

// viscous flux
volScalarField d2nue
(
    IOobject
    (
        "d2nue",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "d2nue", 
           dimVelocity*dimVelocity/dimTime,
           0.0
      )
);

// viscous dissipation
volScalarField eps_nu
(
    IOobject
    (
        "eps_nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "eps_nu", 
           dimVelocity*dimVelocity/dimTime,
           0.0
      )
);

// numerical dissipation at each node
volScalarField dissipation
(
    IOobject
    (
        "dissipation",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "dissipation", 
           dimVelocity*dimVelocity/dimTime,
           0.0
      )
);
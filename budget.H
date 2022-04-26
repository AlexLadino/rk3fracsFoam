      Info<< "Calculando UPrime" << nl << endl;
      UPrime = U - UMean;
      Info<< "Calculando gradUPrime" << nl << endl;      
      gradUPrime = fvc::grad(UPrime);
      Info<< "Calculando pPrime" << nl << endl;      
      pPrime = p - pMean;
      // TKE instantanea
      Info<< "Calculando k" << nl << endl;      
      k = 0.5*magSqr(UPrime);
      // Info<< "Calculando kMean" << nl << endl;      
      // kMean = 0.5*tr(UPrime2Mean);
 
      // TKE Budget terms
      // convection (hacer fieldAverage convectionMean)
      Info<< "Calculando convection" << nl << endl;
      convection = -UMean & fvc::grad(k);
      // Production (se calcula en fase de postProcesamiento)
            // productionMean = -UPrime2Mean && fvc::grad(UMean);
      // viscous dissipation (hacer fieldAverage viscousDissipationMean)
      Info<< "Calculando viscousDissipation" << nl << endl;      
      viscousDissipation = -turbulence->nuEff()*tr(gradUPrime & gradUPrime.T());
      // transport Viscous Diffusion (se calcula en fase de postProcesamiento)
            // transportViscousDiffusionMean = turbulence->nuEff()*fvc::laplacian(kMean); 
      // Transport Velocity-Pressure Gradient 
      // 1. hacer fieldAverage transportVelocityPressureMean
      // 2. en fase de postProcesamiento hacer -div(transportVelocityPressureMean);
      Info<< "Calculando transportVelocityPressure" << nl << endl;      
      transportVelocityPressureField = pPrime*UPrime;
      // Transport by turbulence
      // 1. hacer fieldAverage transportTurbulenceMean
      // 2. en fase de postProcesamiento hacer -div(transportTurbulenceMean);
      Info<< "Calculando transportTurbulence" << nl << endl;      
      transportTurbulenceField = k*UPrime;



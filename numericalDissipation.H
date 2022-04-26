// Begin numericl dissipation calculation

// kinetic energy
e = 0.5*magSqr(U);

// requires 2ns order centered time differentiation requires us to be on at
// least iteration 3
/*
if (runTime.timeIndex() >= 3)
{
      // velocity gradient tensor
      volTensorField gradU("gradU", fvc::grad(U));  
      // kinetic energy flux
      deuidxi = fvc::div(phi,e);     
      // unsteady term of kinetic energy
      dedt = fvc::DDt(phi,e) - deuidxi;
      // pressure flux
      dpuidxi = fvc::div(phi,p);
      // viscous flux
      d2nue = fvc::laplacian(turbulence->nuEff(),e);
      // viscous dissipation
      eps_nu = turbulence->nuEff()*(gradU && gradU);  
       
      // numerical dissipation at each node
      dissipation = -dedt - deuidxi - dpuidxi + d2nue - eps_nu; 
      // integral for numerical dissipation in the whole domain
      // Info<< "domainIntegrate=" << fvc::domainIntegrate(dissipation) << endl;
}
// End numerical dissipation calculation
*/

// ORIGINAL DEL PAPER
if (runTime.timeIndex() >= 3)
{
      // velocity gradient tensor
      const volTensorField gradU("gradU", fvc::grad(U.oldTime()));    
      // unsteady term of kinetic energy
      dedt = (e - e.oldTime().oldTime())/(2.0*runTime.deltaT());
      // kinetic energy flux
      deuidxi = fvc::div(phi,e.oldTime());
      // pressure flux
      dpuidxi = fvc::div(phi,p.oldTime());
      // viscous flux
      d2nue = fvc::laplacian(turbulence->nuEff(),e.oldTime());
      // viscous dissipation
      eps_nu = turbulence->nuEff()*(gradU && gradU);  
       
      // numerical dissipation at each node
      dissipation = -dedt - deuidxi - dpuidxi + d2nue - eps_nu;   
	       
      // integral for numerical dissipation in the whole domain
      // Info<< "domainIntegrate=" << fvc::domainIntegrate(dissipation) << endl;
}


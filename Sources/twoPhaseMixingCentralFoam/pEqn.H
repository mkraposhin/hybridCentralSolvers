if (pimple.corrPISO() > 1)
{
    thermo.correct();
}

volScalarField rhoT (thermo.rho());
rhoHat = rhoT - psi*p;
volVectorField HbyA ("HbyA", U);
volScalarField rAU ("rAU", 1.0 / UEqn.A());

HbyA = UEqn.H() * rAU;

{
    #include "updateCentralWeights.H"
    #include "rhoUpdateCentralFields.H"

    surfaceScalarField rho_nei = fvc::interpolate(rho, nei, "reconstruct(rho)");
    surfaceScalarField rho_own = fvc::interpolate(rho, own, "reconstruct(rho)");
    surfaceScalarField rhoLiq_nei = fvc::interpolate(thermo.thermoLiq().rho(), nei, "reconstruct(rho)");
    surfaceScalarField rhoLiq_own = fvc::interpolate(thermo.thermoLiq().rho(), own, "reconstruct(rho)");
    surfaceScalarField rhoGas_nei = fvc::interpolate(thermo.thermoGas().rho(), nei, "reconstruct(rho)");
    surfaceScalarField rhoGas_own = fvc::interpolate(thermo.thermoGas().rho(), own, "reconstruct(rho)");
    surfaceScalarField YbarLiqf   = min(max(fvc::interpolate(YbarLiq, "reconstruct(YbarLiq)"),0.0),1.0);
    surfaceScalarField YbarGasf   = 1.0 - YbarLiqf;

    surfaceScalarField deltaRho_nei = ((YbarLiqf*rhoLiq_nei + YbarGasf*rhoGas_nei) - rho_nei)*aSf;
    surfaceScalarField deltaRho_own = (rho_own - (YbarLiqf*rhoLiq_own + YbarGasf*rhoGas_own))*aSf;
    //surfaceScalarField deltaRho_nei = -rhoHat_nei*aSf;
    //surfaceScalarField deltaRho_own =  rhoHat_own*aSf;
    surfaceScalarField deltaRho_tot = deltaRho_nei + deltaRho_own;
    //volScalarField divDeltaRho = fvc::div(deltaRho_tot);
    //divDeltaRho.rename("divDeltaRho");
    //divDeltaRho.write();
    //volScalarField divPhiRhoHat = fvc::div(phiRhoHat_own) + fvc::div(phiRhoHat_nei);
    //divPhiRhoHat.rename("divPhiRhoHat");
    //divPhiRhoHat.write();
    phiRhoHat_own += deltaRho_own;
    phiRhoHat_nei += deltaRho_nei;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn_own
        (
            fvm::div(phid_own,p) + fvc::div(phiRhoHat_own) - fvm::laplacian(Dp_own,p)
        );
    
        fvScalarMatrix pEqn_nei
        (
            fvm::div(phid_nei,p) + fvc::div(phiRhoHat_nei) - fvm::laplacian(Dp_nei,p)
        );
    
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi,p)
            +
            fvc::ddt(rhoHat)
            +
            //psi.oldTime()*fvm::ddt(p)+
            pEqn_own
            +
            pEqn_nei
            //+
            //fvc::div(deltaRho_tot)
            //-
            //thermo.thermoLiq().rho()*(fvc::ddt(YbarLiq) + (U & fvc::grad(YbarLiq)))
            //-
            //thermo.thermoGas().rho()*(fvc::ddt(YbarGas) + (U & fvc::grad(YbarGas)))
            //-
            //rhoHat * (U & fvc::grad(YLiq))
            //-
            //rhoHat * fvc::ddt(YLiq)
            ==
            fvOptions(psi,p,rho.name())
        );
    
        fvOptions.constrain(pEqn);
    
        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    
        if (pimple.finalNonOrthogonalIter())
        {
            phi_own = pEqn_own.flux() + phiRhoHat_own;
            phi_nei = pEqn_nei.flux() + phiRhoHat_nei;
            phi = phi_own + phi_nei;
        }
    }
    
    p_own = fvc::interpolate(p, own, "reconstruct(p)");
    p_nei = fvc::interpolate(p, nei, "reconstruct(p)");
    
    gradp = fvc::div((alpha_own*p_own + alpha_nei*p_nei)*mesh.Sf());

    #include "limitPressureGradient.H"
}

U = HbyA - rAU * gradp;
U.correctBoundaryConditions();
fvOptions.correct(U);

Info << "max(mag(U)): " << max(mag(U)).value() << endl;
Info << "max/min p: " << max(p).value() << "/" << min(p).value() << endl;

//p = max(p, pMin);
thermo.correctRealDensities();
rhoT = 1.0 / (YLiq / thermo.thermoLiq().rho() + YGas / thermo.thermoGas().rho());
//rhoT = rhoHat + psi*p;
#include "massEqn.H"

#include "centralContinuityErrs.H"
rho = rhoT;

//volScalarField ddtRho(fvc::ddt(rho));
//ddtRho.rename("ddtRho");
//ddtRho.write();

rho_own = fvc::interpolate(rho, own, "reconstruct(rho)");
rho_nei = fvc::interpolate(rho, nei, "reconstruct(rho)");
correctCentralACMIInterpolation(rho_nei);

//
//END-OF-FILE
//


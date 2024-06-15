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

    const volScalarField& psiLiqEff = thermo.psiLiqEff();
    const volScalarField& psiGasEff = thermo.psiGasEff();
    volScalarField OneByRhoSqr ("OneByRhoSqr",Foam::pow(YLiq / thermo.thermoLiq().rho() + YGas / thermo.thermoGas().rho(), -2.0));
    volScalarField YLiqCoef (YLiq*OneByRhoSqr);
    volScalarField YGasCoef (YGas*OneByRhoSqr);
    volScalarField rhoLiqHat = YLiq*rhoT - thermo.psiLiqEff()*YLiqCoef*p;
    volScalarField rhoGasHat = YGas*rhoT - thermo.psiGasEff()*YLiqCoef*p;



    volScalarField rhoLiqEff = YLiq*rhoT;

    surfaceScalarField rhoLiqEff_nei = fvc::interpolate(rhoLiqEff, nei, "reconstruct(rho)");
    surfaceScalarField rhoLiqEff_own = fvc::interpolate(rhoLiqEff, own, "reconstruct(rho)");

    surfaceScalarField rho_nei = fvc::interpolate(rho, nei, "reconstruct(rho)");
    surfaceScalarField rho_own = fvc::interpolate(rho, own, "reconstruct(rho)");
    surfaceScalarField rhoLiq_nei = fvc::interpolate(thermo.thermoLiq().rho(), nei, "reconstruct(rho)");
    surfaceScalarField rhoLiq_own = fvc::interpolate(thermo.thermoLiq().rho(), own, "reconstruct(rho)");
    surfaceScalarField rhoGas_nei = fvc::interpolate(thermo.thermoGas().rho(), nei, "reconstruct(rho)");
    surfaceScalarField rhoGas_own = fvc::interpolate(thermo.thermoGas().rho(), own, "reconstruct(rho)");

    surfaceScalarField rhoLiqHat_nei = fvc::interpolate(rhoLiqHat, nei, "reconstruct(rho)");
    surfaceScalarField rhoLiqHat_own = fvc::interpolate(rhoLiqHat, own, "reconstruct(rho)");
    surfaceScalarField rhoGasHat_nei = fvc::interpolate(rhoGasHat, nei, "reconstruct(rho)");
    surfaceScalarField rhoGasHat_own = fvc::interpolate(rhoGasHat, own, "reconstruct(rho)");

    surfaceScalarField psiLiqEff_nei = fvc::interpolate(psiLiqEff, nei, "reconstruct(rho)");
    surfaceScalarField psiLiqEff_own = fvc::interpolate(psiLiqEff, own, "reconstruct(rho)");
    surfaceScalarField psiGasEff_nei = fvc::interpolate(psiGasEff, nei, "reconstruct(rho)");
    surfaceScalarField psiGasEff_own = fvc::interpolate(psiGasEff, own, "reconstruct(rho)");
    
    surfaceScalarField YbarLiqf   = min(max(fvc::interpolate(YbarLiq, "reconstruct(YbarLiq)"),0.0),1.0);
    surfaceScalarField YbarGasf   = 1.0 - YbarLiqf;

    surfaceScalarField InvRho2f   = fvc::interpolate(OneByRhoSqr, "reconstruct(OneByRhoSqr)");
    surfaceScalarField YLiqf      = min(max(fvc::interpolate(YLiq, "reconstruct(YLiq)"),0.0),1.0);
    surfaceScalarField YGasf      = 1.0 - YLiqf;

    //surfaceScalarField YbarLiqf   = (alpha_own*rhoLiqEff_own + alpha_nei*rhoLiqEff_nei) / (alpha_own*rhoLiq_own + alpha_nei*rhoLiq_nei);
    //surfaceScalarField YbarGasf   = 1.0 - YbarLiqf;


    surfaceScalarField deltaRho_nei = ((YbarLiqf*rhoLiq_nei + YbarGasf*rhoGas_nei) - rho_nei)*aSf;
    surfaceScalarField deltaRho_own = (rho_own - (YbarLiqf*rhoLiq_own + YbarGasf*rhoGas_own))*aSf;
    //surfaceScalarField deltaRho_nei = -rhoHat_nei*aSf;
    //surfaceScalarField deltaRho_own =  rhoHat_own*aSf;
    //surfaceScalarField deltaRho_tot = deltaRho_nei + deltaRho_own;
    //volScalarField divDeltaRho = fvc::div(deltaRho_tot);
    //divDeltaRho.rename("divDeltaRho");
    //divDeltaRho.write();
    //volScalarField divPhiRhoHat = fvc::div(phiRhoHat_own) + fvc::div(phiRhoHat_nei);
    //divPhiRhoHat.rename("divPhiRhoHat");
    //divPhiRhoHat.write();

    //surfaceScalarField deltaPsi_nei = ((YbarLiqf*psiLiqEff_nei + YbarGasf*psiGasEff_nei) - psi_nei)*aSf;
    //surfaceScalarField deltaPsi_own = (psi_own - (YbarLiqf*psiLiqEff_own + YbarGasf*psiGasEff_own))*aSf;
    surfaceScalarField deltaPsi_nei = ((YLiqf*InvRho2f*psiLiqEff_nei + YGasf*InvRho2f*psiGasEff_nei) - psi_nei)*aSf;
    surfaceScalarField deltaPsi_own = (psi_own - (YLiqf*InvRho2f*psiLiqEff_own + YGasf*InvRho2f*psiGasEff_own))*aSf;

    p_nei = fvc::interpolate(p, nei, "reconstruct(p)");
    p_own = fvc::interpolate(p, own, "reconstruct(p)");
    
    rhoLiqHat_nei = YbarLiqf*rhoLiq_nei - YLiqf*InvRho2f*psiLiqEff_nei*p_nei;
    rhoLiqHat_own = YbarLiqf*rhoLiq_own - YLiqf*InvRho2f*psiLiqEff_own*p_own;
    rhoGasHat_nei = YbarGasf*rhoGas_nei - YGasf*InvRho2f*psiGasEff_nei*p_nei;
    rhoGasHat_own = YbarGasf*rhoGas_own - YGasf*InvRho2f*psiGasEff_own*p_own;

    surfaceScalarField deltaRhoHat_nei = ((rhoLiqHat_nei + rhoGasHat_nei) - rhoHat_nei)*aSf;
    surfaceScalarField deltaRhoHat_own = (rhoHat_own - (rhoLiqHat_own + rhoGasHat_own))*aSf;

    //looks like we need better function for kappa in two-phase mode

    if (runTime.outputTime())
    {
        phiRhoHat_nei.write();
        phiRhoHat_own.write();
        deltaRho_nei.rename("deltaRho_nei");
        deltaRho_own.rename("deltaRho_own");
        deltaRho_nei.write();
        deltaRho_own.write();
        alpha_own.write();
        alpha_nei.write();
        aSf.write();
    }
    //phiRhoHat_nei += deltaRho_nei;
    //phiRhoHat_own += deltaRho_own;
    //phiRhoHat_nei = (YbarLiqf*rhoLiqHat_nei + YbarGasf*rhoGasHat_nei)*aphiv_nei;
    //phiRhoHat_own = (YbarLiqf*rhoLiqHat_own + YbarGasf*rhoGasHat_own)*aphiv_own;


    phid_nei += deltaPsi_nei;
    phid_own += deltaPsi_own;
    phiRhoHat_nei += deltaRhoHat_nei;
    phiRhoHat_own += deltaRhoHat_own;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn_nei
        (
            fvm::div(phid_nei,p) + fvc::div(phiRhoHat_nei) - fvm::laplacian(Dp_nei,p)
        );
        fvScalarMatrix pEqn_own
        (
            fvm::div(phid_own,p) + fvc::div(phiRhoHat_own) - fvm::laplacian(Dp_own,p)
        );

        fvScalarMatrix pEqn
        (
            fvm::ddt(psi,p)
            +
            fvc::ddt(rhoHat)
            +
            pEqn_nei
            +
            pEqn_own
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
    //gradp = fvc::grad(p);

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


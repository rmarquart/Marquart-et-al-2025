const scalar pi = constant::mathematical::pi;

gradU = fvc::grad(U);
gradU.correctBoundaryConditions();

epsilon = 0.5*(gradU + gradU.T());
epsilon.correctBoundaryConditions();

forAll(Delta,i)
{
	Delta[i] = (Foam::sqrt((sqr(epsilon[i].xx()) +
		sqr(epsilon[i].yy()))*(scalar(1)+pow(e.value(),-2)) +
	4*pow(e.value(),-2)*sqr(epsilon[i].xy()) +
	2*epsilon[i].xx()*epsilon[i].yy()*(scalar(1)-pow(e.value(),-2))));
}

Delta = max(Delta,DeltaMin);
Delta.correctBoundaryConditions();

pF = pstar*h;
pF.correctBoundaryConditions();

zetaF = pF/(2*Delta);
zetaF.correctBoundaryConditions();

etaF = zetaF/(e*e);
etaF.correctBoundaryConditions();

SigmaF = etaF*twoSymm(gradU) + (zetaF-etaF)*Ei*tr(gradU) - (pF*Ei)/2.0;
SigmaF.correctBoundaryConditions();

epsilonMag = Foam::sqrt(2.0) * mag(symm(gradU));
epsilonMag.correctBoundaryConditions();

etaInf = etaInf0*pow((h/hG0),alphaInf);
etaInf.correctBoundaryConditions();
etaNil = etaNil0*pow((h/hG0),alphaNil);
etaNil.correctBoundaryConditions();

etaG = (etaInf + (etaNil - etaInf)/(1 + pow(K*epsilonMag,m)));
etaG.correctBoundaryConditions();

surfaceScalarField etaF_zetaF_ = fvc::interpolate((alpha1*etaF)/rho + (alpha1*zetaF)/rho);
surfaceScalarField etaG_ = fvc::interpolate((alpha2*etaG)/rho);

SigmaG = etaG*(gradU + gradU.T());
SigmaG.correctBoundaryConditions();

eta = alpha1*etaF + alpha2*etaG;
eta.correctBoundaryConditions();

Sigma = alpha1*SigmaF + alpha2*SigmaG;
Sigma.correctBoundaryConditions();

volTensorField LHSPart1 = alpha1*((etaF*gradU.T() + (zetaF - etaF)*Ei*tr(gradU) - zetaF*gradU)/rho);
LHSPart1.correctBoundaryConditions();

volTensorField LHSPart2 = alpha2*((etaG*gradU.T())/rho);
LHSPart2.correctBoundaryConditions();

// waves:
dimensionedScalar Tp = Foam::sqrt(((2.0*pi)*lwave)/g);
// Pout<< "Wave period Tp = " << Tp << nl << endl;

dimensionedScalar omega_w = (2.0*pi)/Tp;

forAll(k,cellI)
{
	if (mag(Foam::sqrt(sqr(Uw[cellI][0])+sqr(Uw[cellI][1]))) < 1e-8)
	{
		k = sqr(omega_w)/g;
	}
	else
	{
		k[cellI] = (2.0*Foam::sqrt(sqr(Uw[cellI][0])+sqr(Uw[cellI][1]))*omega_w.value()+g.value()-Foam::sqrt(4.0*Foam::sqrt(sqr(Uw[cellI][0])+sqr(Uw[cellI][1]))*g.value()*omega_w.value()+sqr(g.value())))/(2.0*sqr(Foam::sqrt(sqr(Uw[cellI][0])+sqr(Uw[cellI][1]))));
	}
}

k1 = k*Foam::cos(thetaWa);
k2 = k*Foam::sin(thetaWa);

theta = omega_w*runTime - k1*mesh.C().component(0) - k2*mesh.C().component(1);
theta.correctBoundaryConditions();

z_sur = a*Foam::sin(theta);
z_sur.correctBoundaryConditions();

gradAlpha = fvc::grad(alpha1);

volScalarField magGradAlpha = mag(gradAlpha);
scalar maxGrad = gMax(magGradAlpha);

forAll(interfaceRegion,i)
{
	if (magGradAlpha[i] > maxGrad/2.0)
	{
		interfaceRegion[i] = 1;
	}
	else
	{
		interfaceRegion[i] = 0;
	}
}

nHat = -(gradAlpha/(mag(gradAlpha) + deltaN_)*interfaceRegion);

hw = (rho/rhoW)*h;

tauWI = nHat*sqr(omega_w)*hw*a*Foam::cos(theta);

volScalarField Uwtx = Uw.component(0) + (omega_w-k1*Uw.component(0))*a*(Foam::sin(theta)*Foam::cos(thetaWa));
volScalarField Uwty = Uw.component(1) + (omega_w-k2*Uw.component(1))*a*(Foam::sin(theta)*Foam::sin(thetaWa));
volScalarField Uwtz = (omega_w-k1*Uw.component(0)-k2*Uw.component(1))*a*Foam::cos(theta);

forAll(Uwt, cellI)
{
	Uwt[cellI] = vector(Uwtx[cellI], Uwty[cellI], Uwtz[cellI]);
}

volScalarField twt_mag = (rhoW/rho)*Cw*(mag(Uwt-U));
twt_mag.correctBoundaryConditions();

twt_explicit = twt_mag*(cos(thetaW)*Uwt + sin(thetaW)*(vector(0,0,1)^(Uwt-U)));
twt_explicit.correctBoundaryConditions();

// // wind velocity profile
// Ustar = u10*pow(C10,0.5);
// z0 = (alphaW*sqr(Foam::sqrt(sqr(Ustar.component(0)) + sqr(Ustar.component(1)))))/g;

// // volScalarField h_emerged1 = h1-(rho1/rhoW)*h1;
// volScalarField h_emerged = alpha1*h-((alpha1*rho)/rhoW)*(alpha1*h);

// hz = z0 + h_emerged;

// forAll(z_wind, cellI)
// {
// 	z_wind[cellI] = hz[cellI] + z_sur[cellI] + a.value();
// }
// z_wind.correctBoundaryConditions();

// forAll(Ca,i)
// {
// 	Ca[i] = Ca_mult.value()*(sqr(kappa.value())/sqr(Foam::log(z_wind[i])-Foam::log(z0[i])));
// }
// Ca.correctBoundaryConditions();

// forAll(Ua_component,i)
// {
// 	Ua_component[i] = Foam::log(z_wind[i])-Foam::log(z0[i]);
// }
// Ua_component.correctBoundaryConditions();

// forAll(Ua,i)
// {
// 	Ua[i][0] = (Ustar[i][0]/kappa.value())*Ua_component[i];
// 	Ua[i][1] = (Ustar[i][1]/kappa.value())*Ua_component[i];
// 	Ua[i][2] = 0;
// }
// Ua.correctBoundaryConditions();

// volScalarField ta_mag = rhoA/rho*Ca*(mag(Ua-U));
// ta_mag.correctBoundaryConditions();

// tauAir_explicit = ta_mag*(cos(thetaA)*Ua + sin(thetaA)*(vector(0,0,1)^(Ua-U)));
// tauAir_explicit.correctBoundaryConditions();

// volScalarField tw_mag = rhoW/rho*Cw*(mag(Uw-U));
// tw_mag.correctBoundaryConditions();

// volVectorField tauWater_explicit = tw_mag*(cos(thetaW)*Uw + sin(thetaW)*(vector(0,0,1)^(Uw-U)));
// tauWater_explicit.correctBoundaryConditions();

// fvVectorMatrix UEqn
// (
// 	fvm::ddt(h,U)
// 	+ fvm::div(h_Phi,U)
// 	- fvm::laplacian(etaF_zetaF_,U)
// 	- fvm::laplacian(etaG_,U)
// 	+ fvm::Sp(cos(thetaW)*twt_mag,U)
// 	// + fvm::Sp(cos(thetaA)*ta_mag,U)
// 	==
// 	fvc::div(LHSPart1)
// 	+ fvc::div(LHSPart2)
// 	// + tauAir_explicit
// 	+ twt_explicit
// 	+ tauWI
// 	- alpha1*fvc::grad(pF/(2*rho))
// 	);

// UEqn.relax();
// solve(UEqn);
// U.correctBoundaryConditions();

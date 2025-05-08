
volScalarField alpha3
 (
     IOobject
     (
         "alpha3",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
 );

Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
        ),
    mesh
    );

#   include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties 
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
        )
    );

const dimensionedScalar rhoA = dimensionedScalar("rhoA", transportProperties);
const dimensionedScalar rhoW = dimensionedScalar("rhoW", transportProperties);
const dimensionedScalar thetaA = dimensionedScalar("thetaA", transportProperties);
const dimensionedScalar thetaW = dimensionedScalar("thetaW", transportProperties);
const dimensionedScalar thetaWa = dimensionedScalar("thetaWa", transportProperties);
const dimensionedScalar g = dimensionedScalar("g", transportProperties);
const dimensionedScalar a = dimensionedScalar("a", transportProperties);
const dimensionedScalar lwave = dimensionedScalar("lwave", transportProperties);
const dimensionedScalar kappa = dimensionedScalar("kappa", transportProperties);
const dimensionedScalar alphaW = dimensionedScalar("alphaW", transportProperties);
const dimensionedScalar alphaInf = dimensionedScalar("alphaInf", transportProperties);
const dimensionedScalar alphaNil = dimensionedScalar("alphaNil", transportProperties);
const dimensionedScalar C10 = dimensionedScalar("C10", transportProperties);
const dimensionedScalar Ca_mult = dimensionedScalar("Ca_mult", transportProperties);
const dimensionedScalar etaInf0 = dimensionedScalar("etaInf0", transportProperties);
const dimensionedScalar etaNil0 = dimensionedScalar("etaNil0", transportProperties);
const dimensionedScalar hG0 = dimensionedScalar("hG0", transportProperties);
const dimensionedScalar K = dimensionedScalar("K", transportProperties);
const dimensionedScalar m = dimensionedScalar("m", transportProperties);
const dimensionedScalar pstar = dimensionedScalar("pstar", transportProperties);
const dimensionedScalar e = dimensionedScalar("e", transportProperties);
const dimensionedScalar DeltaMin = dimensionedScalar("DeltaMin", transportProperties);

//NOTE: constructing here what I believe to be the bare minimum for MULES with 2 phases
twoPhaseMixture mixture(U.mesh(), transportProperties);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

interfaceProperties interface(alpha1, U, transportProperties);

Info << mixture.phase1Name() << endl;
Info << mixture.phase2Name() << endl;

wordList phaseList = wordList(transportProperties.lookup("phases"));

Info<< "Reading field h\n" << endl;

volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
        ),
    mesh
    );

const dimensionedScalar Cw1 = dimensionedScalar("Cw", transportProperties.subDict(phaseList[0]));
const dimensionedScalar Cw2 = dimensionedScalar("Cw", transportProperties.subDict(phaseList[1]));

volScalarField Cw
(
 IOobject
 (
    "Cw",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT
    ),
 alpha1*Cw1 + alpha2*Cw2
 );

const dimensionedScalar rho1 = dimensionedScalar("rho", transportProperties.subDict(phaseList[0]));
const dimensionedScalar rho2 = dimensionedScalar("rho", transportProperties.subDict(phaseList[1]));

volScalarField rho
(
 IOobject
 (
    "rho",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
    ),
 alpha1*rho1 + alpha2*rho2
 );

volVectorField Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedVector("Ua", dimensionSet(0,1,-1,0,0,0,0),vector(0,0,0))
    );

volVectorField Uw
(
    IOobject
    (
        "Uw",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedVector("Uw", dimensionSet(0,1,-1,0,0,0,0),vector(0,0,0))
    );

volTensorField epsilon
(
    IOobject
    (
       "epsilon",
       runTime.timeName(),
       mesh,
       IOobject::MUST_READ,
       IOobject::AUTO_WRITE
       ),
    mesh,
    dimensionedTensor("epsilon", dimensionSet(0,0,-1,0,0,0,0), tensor::zero)
    );

volScalarField epsilonMag
(
    IOobject
    (
        "epsilonMag",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar("epsilonMag", dimensionSet(0,0,-1,0,0,0,0), scalar(0))
    );

volTensorField SigmaF
(
    IOobject
    (
        "SigmaF",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedTensor("SigmaF", dimensionSet(1,0,-2,0,0,0,0), tensor::zero)
    );

volTensorField SigmaG
(
    IOobject
    (
        "SigmaG",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedTensor("SigmaG", dimensionSet(1,0,-2,0,0,0,0), tensor::zero)
    );

volTensorField Sigma
(
    IOobject
    (
        "Sigma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedTensor("Sigma", dimensionSet(1,0,-2,0,0,0,0), tensor::zero)
    );

volTensorField gradU
(
    IOobject
    (
        "gradU",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedTensor("gradU", dimensionSet(0,0,-1,0,0,0,0), tensor::zero)
    );

volScalarField Delta
(
    IOobject
    (
        "Delta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar("Delta", dimensionSet(0,0,-1,0,0,0,0), scalar(0))
    );

volScalarField zetaF
(
    IOobject
    (
        "zetaF",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("zetaF", dimensionSet(1,0,-1,0,0,0,0), scalar(0))
    );

volScalarField etaInf
(
    IOobject
    (
        "etaInf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("etaInf", dimensionSet(1,0,-1,0,0,0,0), scalar(0))
    );

volScalarField etaNil
(
    IOobject
    (
        "etaNil",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("etaNil", dimensionSet(1,0,-1,0,0,0,0), scalar(0))
    );

volScalarField eta
(
    IOobject
    (
        "eta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar("eta", dimensionSet(1,0,-1,0,0,0,0), scalar(0))
    );

volScalarField etaF
(
    IOobject
    (
        "etaF",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("etaF", dimensionSet(1,0,-1,0,0,0,0), scalar(0))
    );

volScalarField etaG
(
    IOobject
    (
        "etaG",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("etaG", dimensionSet(1,0,-1,0,0,0,0), scalar(0))
    );

volScalarField pF
(
    IOobject
    (
        "pF",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("pF", dimensionSet(1,0,-2,0,0,0,0), scalar(0))
    );

surfaceScalarField h_Phi = fvc::interpolate(h)*phi;
tensor Ei(1,0,0,0,1,0,0,0,1);

volVectorField tauAir_explicit
(
    IOobject
    (
        "tauAir_explicit",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedVector("tauAir_explicit", dimensionSet(0,2,-2,0,0,0,0),vector(0,0,0))
    );

// waves:
volScalarField hw
(
    IOobject
    (
        "hw",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("hw", dimensionSet(0,1,0,0,0,0,0), scalar(0))
    );

volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("k", dimensionSet(0,-1,0,0,0,0,0),scalar(0))
    );

volScalarField k1
(
    IOobject
    (
        "k1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("k1", dimensionSet(0,-1,0,0,0,0,0),scalar(0))
    );

volScalarField k2
(
    IOobject
    (
        "k2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("k2", dimensionSet(0,-1,0,0,0,0,0),scalar(0))
    );

volVectorField tauWI
(
    IOobject
    (
        "tauWI",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedVector("tauWI", dimensionSet(0,2,-2,0,0,0,0),vector(0,0,0))
    );

volVectorField twt_explicit
(
    IOobject
    (
        "twt_explicit",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedVector("twt_explicit", dimensionSet(0,2,-2,0,0,0,0),vector(0,0,0))
    );

volScalarField z_sur
(
    IOobject
    (
        "z_sur",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("z_sur", dimensionSet(0,1,0,0,0,0,0),scalar(0))
    );

volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("theta", dimensionSet(0,0,0,0,0,0,0),scalar(0))
    );

volVectorField Uwt
(
    IOobject
    (
        "Uwt",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedVector("Uwt", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0))
    );

dimensionedScalar deltaN_
(
    "deltaN_",
    dimensionSet(0,-1,0,0,0,0,0),
    scalar(1e-8)
    );

volScalarField interfaceRegion
(
    IOobject
    (
        "interfaceRegion",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("interfaceRegion", dimensionSet(0,0,0,0,0,0,0), scalar(0))
    );

volVectorField gradAlpha
(
    IOobject
    (
        "gradAlpha",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedVector("gradAlpha", dimensionSet(0,-1,0,0,0,0,0), vector(0,0,0))
    );

volVectorField nHat
(
    IOobject
    (
        "nHat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedVector("nHat", dimensionSet(0,0,0,0,0,0,0), vector(0,0,0))
    );

// wind profile
volScalarField Ca
(
    IOobject
    (
        "Ca",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("Ca", dimensionSet(0,0,0,0,0,0,0), scalar(0))
    );

volScalarField z0
(
    IOobject
    (
        "z0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("z0", dimensionSet(0,1,0,0,0,0,0), scalar(0))
    );

volVectorField u10
(
    IOobject
    (
        "u10",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedVector("u10", dimensionSet(0,1,-1,0,0,0,0),vector(1e-10,1e-10,0))
    );

volVectorField Ustar
(
    IOobject
    (
        "Ustar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedVector("Ustar", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0))
    );

volScalarField hz
(
    IOobject
    (
        "hz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("hz", dimensionSet(0,1,0,0,0,0,0), scalar(0))
    );

volScalarField z_wind
(
    IOobject
    (
        "z_wind",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("z_wind", dimensionSet(0,1,0,0,0,0,0), scalar(0))
    );

volScalarField Ua_component
(
    IOobject
    (
        "Ua_component",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
    mesh,
    dimensionedScalar("Ua_component", dimensionSet(0,0,0,0,0,0,0), scalar(0))
    );

mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"

//NOTE: constructing here isoAdvector:
isoAdvection advector(alpha1, phi, U);

while (pimple.correctNonOrthogonal())
{
    solve
    (
        fvm::ddt(h)
        + fvm::div(phi,h)
        );
}
h.correctBoundaryConditions();

forAll(h,i)
{
  if (h[i] < 1e-8)
  {
      h[i] = 1e-8;
  }
}

while (pimple.correctNonOrthogonal())
{
    solve
    (
        fvm::ddt(alpha3)
        + fvm::div(phi,alpha3)
        );
}
alpha3.correctBoundaryConditions();

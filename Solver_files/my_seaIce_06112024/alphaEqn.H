// WIce-FOAM 1.0: Coupled dynamic and thermodynamic modelling of heterogeneous sea ice
// Copyright (C) 2025 Marquart et al.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// For questions or inquiries, contact: rutger.marquart@uct.ac.za

{
    // Temporarily making U relative to mesh motion
    if (mesh.moving())
    {
        U -= fvc::reconstruct(mesh.phi());
    }

    // Updating alpha1
    #include "alphaSuSp.H"
    advector.advect(Sp, Su);

    // Making U absolute again after advection step
    if (mesh.moving())
    {
        U += fvc::reconstruct(mesh.phi());
    }

    // #include "rhofs.H"
    // rhoPhi = advector.getRhoPhi(rho1f, rho2f);

    alpha2 = 1.0 - alpha1;
    interface.correct();
}

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;

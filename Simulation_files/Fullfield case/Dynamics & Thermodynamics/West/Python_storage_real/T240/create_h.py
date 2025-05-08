# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 09:35:07 2023

@author: Ma-re Admin
"""

import os
import numpy as np

# Load the exported variable
ii = float(os.environ.get("ii"))
num = str(os.environ.get("num"))

file_path = f'/scratch/mrqrut001/simulation/postdoc/Python_storage_real/T{num}/h_updated.npy'
h = np.load(file_path)

l1 = '/*--------------------------------*- C++ -*----------------------------------*'
l2 = '| =========                |                                                 |'
l3 = '| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'
l4 = '|  \\    /   O peration     | Version:  2106                                  |'
l5 = '|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |'
l6 = '|    \\/     M anipulation  |                                                 |'
l7 = '*---------------------------------------------------------------------------*/'
l8 = 'FoamFile'
l9 = '{'
l10= 'version     2.0;'
l11= 'format      ascii;'
l12= 'class       volScalarField;'
l13= 'location    "{}";'.format(ii)
l14= 'object      h;'
l15= '}'
l16= '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'
l17= ' '
l18= 'dimensions      [0 1 0 0 0 0 0];'
l19= ' '
l20= ' '
l21= 'internalField   nonuniform List<scalar>'
l22= '254016'
l23= '('
l24 = ')'
l25 = ';'
l26 = ' '
l27 = ' '
l28 = 'boundaryField'
l29 = '{'
l30 = '    left'
l31 = '    {'
l32 = '        type            cyclic;'
l33 = '    }'
l34 = '    right'
l35 = '    {'
l36 = '        type            cyclic;'
l37 = '    }'
l38 = '    top'
l39 = '    {'
l40 = '        type            cyclic;'
l41 = '    }'
l42 = '    bottom'
l43 = '    {'
l44 = '        type            cyclic;'
l45 = '    }'
l46 = '    frontAndBack'
l47 = '    {'
l48 = '        type            empty;'
l49 = '    }'
l50 = '}'

# Create and write to h.txt
with open('h', 'w') as fid:
    fid.write(l1 + '\n' + l2 + '\n' + l3 + '\n' + l4 + '\n' + l5 + '\n' + l6 + '\n' + l7 + '\n' +
              l8 + '\n' + l9 + '\n' + l10 + '\n' + l11 + '\n' + l12 + '\n' + l13 + '\n' + l14 + '\n' +
              l15 + '\n' + l16 + '\n' + l17 + '\n' + l18 + '\n' + l19 + '\n' + l20 + '\n' + l21 + '\n' +
              l22 + '\n' + l23 + '\n')

# Write h_dummy.txt
np.savetxt('h_dummy.txt', h, delimiter=' ')

# Append h_dummy.txt to h.txt
h_dummy = np.loadtxt('h_dummy.txt')
with open('h', 'a') as fid:
    np.savetxt(fid, h_dummy, delimiter=' ', fmt='%f')

# Append remaining lines to h.txt
with open('h', 'a') as fid:
    fid.write(l24 + '\n' + l25 + '\n' + l26 + '\n' + l27 + '\n' + l28 + '\n' + l29 + '\n' +
              l30 + '\n' + l31 + '\n' + l32 + '\n' + l33 + '\n' + l34 + '\n' + l35 + '\n' +
              l36 + '\n' + l37 + '\n' + l38 + '\n' + l39 + '\n' + l40 + '\n' + l41 + '\n' +
              l42 + '\n' + l43 + '\n' + l44 + '\n' + l45 + '\n' + l46 + '\n' + l47 + '\n' +
              l48 + '\n' + l49 + '\n' + l50 + '\n')
    
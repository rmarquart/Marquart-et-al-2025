# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 09:35:07 2023

@author: Ma-re Admin
"""

import os

# Load the exported variable
ii = float(os.environ.get("ii"))
num = str(os.environ.get("num"))

output_dir = "/scratch/mrqrut001/simulation/postdoc/Python_real/T{}/system".format(num)

controlDict_lines = [
        '/*--------------------------------*- C++ -*----------------------------------*',
        '| =========                |                                                 |',
        '| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |',
        '|  \\    /   O peration     | Version:  2106                                  |',
        '|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |',
        '|    \\/     M anipulation  |                                                 |',
        '*---------------------------------------------------------------------------*/',
        'FoamFile',
        '{',
        '    version       2.0;',
        '    format        ascii;',
        '    class         dictionary;',
        '    location      "system";',
        '    object        controlDict;',
        '}',
        '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //',
        '',
        'application       my_seaIce_06112024_v2;',
        '',
        'startFrom         startTime;',
        '',
        'startTime         {};'.format(ii),
        '',
        'stopAt            endTime;',
        '',
        'endTime           {};'.format(round(ii + 15.2, 1)),
        '',
        'deltaT            0.1;',
        '',
        'writeControl      runTime;',
        '',
        'writeInterval     7.6;',
        '',
        'purgeWrite        0;',
        '',
        'writeFormat       ascii;',
        '',
        'writePrecision    6;',
        '',
        'writeCompression  off;',
        '',
        'timeFormat        general;',
        '',
        'timePrecision     6;',
        '',
        'runTimeModifiable true;',
        '',
        '// ************************************************************************* //'
    ]

# Create and write the controlDict file
controlDict_path = os.path.join(output_dir, 'controlDict'.format(ii))
with open(controlDict_path, 'w') as controlDict_file:
    controlDict_file.write('\n'.join(controlDict_lines))

    ## Update i for the next timestep
    ii += 15.2

print("controlDict file updated successfully.")

#!/usr/bin/env python
"""
Calculates the radial shape distribution for interface residues

An interface residue is defined as a residue that has at least one interface
atom.  An interface atom is defined as an atom that is within distance <rad>
of any atom in the opposing molecule.

The radial shape distribution is a sphere subdivided into <shell> concentric
shells with width <width>.  The total and average shape score, and the
number of FADE points are reported for each shell.

Defaults: shell = 10, width = 1

Sample: A  10 GLN  -0.119 -0.029750   4   0.173  0.010812  16 ...
        -residue- -------shell 1------- -------shell 2-------

"""

# Steven Darnell / sdarnell@biochem.wisc.edu / UW-Madison
#
# Code History
#
# 01/11/2006 - First version.
# 01/27/2006 - Changed input to match FEATURE syntax. [rad_sh width] replaced
#              by [shell width]. Documentation edited.
#
# 2011 - gw, fixed some bugs, including making output go to file and
#        warning/error messages go to stdout to avoid crashing
#        FormatData.py when it tried to parse warning messages as data.
#        
import sys, ProteinInterface

def main(argv=None):
    """
    Main method

    Usage from command-line (brackets optional):
    python fnearby.py mol1 mol2 dad [rad] [shell width]

    """
    if argv == None:
	argv = sys.argv
    argc = len(argv)
    if argc != 4 and argc != 5 and argc != 7:
	print 'Usage: python %s mol1 mol2 dad [rad] [shell width]' % argv[0]
	sys.exit()

    if argc == 4:
	prog, mol1, mol2, dad = argv
	rad, shell, width = 4, 10, 1
    elif argc == 5:
	prog, mol1, mol2, dad, rad = argv
	rad, shell, width = int(rad), 10, 1
    else:
	prog, mol1, mol2, dad, rad, shell, width = argv
	rad, shell, width = int(rad), int(shell), int(width)

    face  = ProteinInterface.ProteinInterface(mol1,mol2)
    sites = face.interfaceResidues(rad)
    shape = face.shapeDistribution(sites,dad,shell,width)

    try:
       fnbfh = open("fnearby.fnb","w")
    except:
       print 'ERROR: Unable to open file "fnearby.fnb" for writing.'
       sys.exit(1)

    for i,id in enumerate(sites):
	print >> fnbfh, '%s %3d %s' % id,
	for stat in shape[i]:
	    if stat[1] == 0:
		    nstat=0.000
	    else:
		    nstat=stat[0]/stat[1]
		
	    print  >> fnbfh, '%7.3f %9.6f %3d' % tuple([stat[0], nstat, stat[1]]),
	print >> fnbfh
    fnbfh.close()   


if __name__ == '__main__':
    main()

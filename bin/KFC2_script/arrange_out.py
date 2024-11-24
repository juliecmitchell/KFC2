"""
Formats data for training decision tree

Formats the KFC data so it can be read by the c5.0 decisiontree method.  
Shape data file lists all interface residues. 
Contact data file lists all contacts identified by WHATIF.
Position data file lists interface position of all residues in train set

"""

# File:    arrange_out.py
# Sources: None
#
# Authors:
#
# Steven Darnell / sdarnell@biochem.wisc.edu / UW-Madison
# Julie Mitchell / mitchell@math.wisc.edu    / UW-Madison
# Xiaolei Zhu    /xzhu32@wisc.edu            / UW-Madison
#
# Reference:
#
# Darnell, S.J., Page, D., Mitchell, J.C. (2006) "An Automated Decision-Tree
# Approach to Predicting Protein Interaction Hot Spots." Proteins, submitted.
#
# Copyright (c) 2006 University of Wisconsin System Board of Regents
#
# All rights reserved. This software may not be redistributed in any form
# without permission from the authors. The software is distributed "as is"
# with no warranty of any kind, expressed or implied.

import sys


def kfc2a_result(inFile):
    """
    kfc2a_result(inFile) -> (float)

    inFile (str)  KFC2A result file name

    result  (list)

    """
    try:
	data = open(inFile,'r').readlines()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))
    result = []
    for line in data:
	x=line.strip()
	result.append(x)
    return result

def kfc2b_result(inFile):
    """
    kfc2b_result(inFile) -> (float)

    inFile  (str)  KFC2B result file name

	result  (list)

    """
    try:
	data = open(inFile,'r').readlines()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))
    result = []
    for line in data:
	x=line.strip()
	result.append(x)
    return result

def initSite(inFile):
    """
    initSite(inFile) -> residue's position at interface

    inFile  (str)  residue's position file name

    sites   (list)

    """
    try:
	data = open(inFile, 'r').readlines()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))
    sites = []
    for line in data:
	resname,ch, nid = line.split()
	sites.append((resname,ch,nid))
    return sites



def toFileResult(outFile,result_a,result_b,sites):
    """
    Write kfc2 formated results

    outFile  (str)   Output file name
    result_a (list)  score by kfc2a
    result_b (list)  score by kfc2b
	sites    (list)  Residues to evaluate

    """
    try:
	output = open(outFile,'w')
	for i in range(len(sites)):
	    output.write('%-4s %2s %4d %8.3f%8.3f\n' % (sites[i][0],sites[i][1],int(sites[i][2]),float(result_a[i]),float(result_b[i])))
	output.close()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))


def toFileSites(outFile,sites):
    """
    Writes list of interface residues

    outFile (str)  Output file name
    sites   (list) Interface residues

    """
    try:
	output = open(outFile,'w')
	for site in sites:
	    output.write('%s %s %4d\n' % (site[0],site[1],int(site[2])))
	output.close()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))

def main(argv=None):
    """
    Main method

    Usage from command-line:
	python arrange_out.py  result_a.file result_b.file sites.file all_result.file

    """
    if argv == None:
	argv = sys.argv
    if len(argv) != 5:
	sys.exit('Usage: python %s result_a.txt result_b.txt sites.txt KFC2_result.txt' \
			% argv[0])

    result_a = kfc2a_result(argv[1])
    result_b = kfc2b_result(argv[2])
    sites = initSite(argv[3])
    toFileResult(argv[4],result_a,result_b,sites)


if __name__ == '__main__':
    main()

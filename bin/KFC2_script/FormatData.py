"""
Formats data for training decision tree

Formats the KFC data so it can be read by the c5.0 decisiontree method.  
Shape data file lists all interface residues. 
Contact data file lists all contacts identified by WHATIF.
Position data file lists interface position of all residues in train set

"""

# File:    FormatData.py
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

def residueTrans(amino_s):
    """
    residueTrans(amino_s) -> 'A->ALA', 'C->CYS', 'R->ARG', or 'like this'
    
    amino_s (str) Amino acid (1-letter)
    
    """
    return {'A':'ALA', 'G':'GLY', 'S':'SER', 'C':'CYS', 'D':'ASP', 'N':'ASN',
            'P':'PRO', 'T':'THR', 'V':'VAL', 'E':'GLU', 'F':'PHE', 'H':'HIS',
            'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'R':'ARG', 'Q':'GLN',
            'W':'TRP', 'Y':'TYR'}[amino_s]

def residueHydrobicity(amino_t):
    """
    residueHydrobicity(amino_t) -> 'ALA->0.31', 'ARG->-1.01', or 'like this'
    
    amino_t (str) Amino acid (3-letter)
    
    """
    return {'ALA':'0.31', 'GLY':'0.00', 'SER':'-0.04', 'CYS':'1.54', 'ASP':'-0.77', 'ASN':'-0.60',
            'PRO':'0.72', 'THR':'0.26', 'VAL':'1.22', 'GLU':'-0.64', 'PHE':'1.79', 'HIS':'0.13',
            'ILE':'1.80', 'LYS':'-0.99', 'LEU':'1.70', 'MET':'1.23', 'ARG':'-1.01', 'GLN':'-0.22',
            'TRP':'2.25', 'TYR':'0.96'}[amino_t]

def residueSize(amino):
    """
    residueSize(amino) -> 'Small', 'Medium', or 'Large'

    amino (str) Amino acid (3-letter)

    """
    return {'ZN2':'Small' ,'ALA':'Small' ,'GLY':'Small' ,'SER':'Small' ,
	    'CYS':'Medium','CYX':'Medium','ASP':'Medium','ASN':'Medium','PRO':'Medium',
	    'THR':'Medium','VAL':'Medium','GLU':'Large' ,'PHE':'Large' ,
	    'HIS':'Large' ,'ILE':'Large' ,'LYS':'Large' ,'LEU':'Large' ,
	    'MET':'Large' ,'MSE':'Large' ,'ARG':'Large' ,'GLN':'Large' ,'TRP':'Large' ,
	    'TYR':'Large' ,'IIP':'Small' , 'DA' :'Large' , 'DC' :'Large' ,
	    'DT'  :'Large' , 'DG' :'Large' , 'A' :'Large' , 'C' :'Large' ,
	    'T'  :'Large' , 'G' :'Large' , 'U' :'Large' ,'ADE':'Large' ,
	    'CYT':'Large' ,'THY':'Large' ,'GUA':'Large' ,'URA':'Large' }[amino]

def nresidueSize(amino):
    """
	nresidueSize(amino) -> solvent accessible area of that residue

    amino (str) Amino acid (3-letter)

    """
    return {'ALA':'102.68' ,'GLY':'70.27' ,'SER':'111.97' ,
	    'CYS':'127.72','ASP':'141.61','ASN':'155.22','PRO':'126.78',
	    'THR':'134.28','VAL':'148.51','GLU':'173.46' ,'PHE':'209.64' ,
	    'HIS':'184.79' ,'ILE':'176.22' ,'LYS':'193.73' ,'LEU':'184.33' ,
	    'MET':'196.87' ,'ARG':'235.45' ,'GLN':'175.42' ,'TRP':'251.78' ,
	    'TYR':'224.68' }[amino]
																		

def residueType(amino):
    """
    residueType(amino) -> 'Nonpolar', 'Polar', or 'Charged'

    amino (str) Amino acid (3-letter)

    """
    return {'ZN2':'Charged' ,'ALA':'Nonpolar','GLY':'Nonpolar','ILE':'Nonpolar',
	    'LEU':'Nonpolar','MET':'Nonpolar','MSE':'Nonpolar','PHE':'Nonpolar',
	    'PRO':'Nonpolar','TRP':'Nonpolar','VAL':'Nonpolar',
	    'ASN':'Polar'   ,'CYS':'Polar'   ,'CYX':'Polar',  'GLN':'Polar'   ,
	    'SER':'Polar'   ,'THR':'Polar'   ,'TYR':'Polar'   ,
	    'ARG':'Charged' ,'ASP':'Charged' ,'GLU':'Charged' ,
	    'HIS':'Charged' ,'LYS':'Charged' ,
            'IIP':'Polar'   , 'DA' :'Charged' ,
            'DC'  :'Charged' , 'DT' :'Charged' , 'DG' :'Charged', 'A' :'Charged' ,
            'C'  :'Charged' , 'T' :'Charged' , 'G' :'Charged', 'U' :'Charged' ,
            'ADE':'Charged' ,'CYT':'Charged' ,'THY':'Charged','GUA':'Charged' ,
	    'URA':'Charged' }[amino]

def initShapeEnvironment(inFile):
    """
    initShapeEnvironment(inFile) -> (shape,fpoint)

    inFile (str)  K-FADE data file name

    shape  (dict) Key: residue  Value: shape scores in 10 shells (list)
    fpoint (dict) Key: residue  Value: FADE points in 10 shells (list) 

    """
    try:
	data = open(inFile,'r').readlines()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))
#    shape  = {}
    fpoint = {}
    fpoint1 = {}
    for line in data:
#	pdbcode,ch,nid,sh0,fp0,sh1,fp1,sh2,fp2,sh3,fp3,sh4,fp4,sh5,fp5,sh6,fp6, \
#	    sh7,fp7,sh8,fp8,sh9,fp9,fp_f5,fp_b5,sh_f5,sh_b5 = line.split()
	ch,nid,resname,sh0,avesh0,fp0,sh1,avesh1,fp1,sh2,avesh2,fp2,sh3,avesh3,fp3,sh4,avesh4,fp4,sh5,avesh5,fp5,sh6,avesh6,fp6,\
		sh7,avesh7,fp7,sh8,avesh8,fp8,sh9,avesh9,fp9 = line.split()
#	key = (pdbcode,ch,nid)
	key = (resname,ch,nid)
#	shape[key]  = [sh0,sh1,sh2,sh3,sh4,sh5,sh6,sh7,sh8,sh9,sh_f5,sh_b5]
#	fpoint[key] = [fp0,fp1,fp2,fp3,fp4,fp5,fp6,fp7,fp8,fp9,fp_f5,fp_b5]
	fpoint[key] = [fp9,]
	fpoint1[key] = [fp8,]
    return fpoint,fpoint1

def initBiochemContacts(inFile):
    """
    initBiochemContacts(inFile) -> contact

    inFile  (str)  K-CON data file name

    contact (dict) Key: residue  Value: contacts [polar,generic,h-bond]

    """
    try:
	data = open(inFile,'r').readlines()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))
    contact = {}
    for line in data:
	pdbcode,ch,nid,p,g,hb = line.split()
	key = (pdbcode,ch,nid)
	contact[key] = [p,g,hb]
    return contact

def initCoreRim(inFile):
    """
    initCoreRim(inFile) -> residue's position at interface

    inFile  (str)  residue's position file name

    corerim (dict) Key: residue  Value: [delta_tot, per, pos_per]

    """
    try:
	data = open(inFile, 'r').readlines()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))
    corerim = {}
    for line in data:
#	pdbcode,ch,nid,delta_tot,per,pos_per = line.split()
	resname,ch,nid,com_abs_tot,mon_abs_tot,delta_tot,per,pos_per = line.split()
	key = (resname,ch,nid)
	corerim[key] = [pos_per,]
    return corerim

def initNewSas(inFile):
    """
    initNewSas(inFile) -> all features related with SAS except pos_per

    inFile  (str)  residue's SAS related file name

    newsas  (dict) Key: residue  Value: [tot_delta, per, sc_delta, com_tot/maxASA, tot_delta/maxASA, com_sc/maxASAsc, sc_delta/maxASAsc, \
	tot_delta*maxASA/MABS_tot]

    """
    try:
	data = open(inFile, 'r').readlines()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))
    newsas = {}
    newsas1 = {}
    for line in data:
#	pdbcode,ch,nid,delta1,per,delta2,ratio1,ratio2,ratio3,ratio4,ratio5,x1,x2,x3,x4,x5,x6,x7,x8 = line.split()
	resname,ch,nid,c_tot,m_tot,delta1,per,c_sc,m_sc,delta2,ratio1,ratio2,ratio3,ratio4,ratio5,x1,x2,x3,x4,x5,x6,x7,x8 = line.split()
#	key = (pdbcode,ch,nid)
	key = (resname,ch,nid)
#	newsas[key] = [delta1,per,delta2,ratio1,ratio2,ratio3,ratio4,ratio5]
	newsas[key] = [delta1,per]
	newsas1[key] = [ratio5,]
    return newsas, newsas1
	
def initMutCom(inFile):
    """
    initMutCom(inFile) -> the possibility to compensate a mutated residue to ALA

    inFile   (str)  residue's mut_compensate file name

    mut_com  (dict) Key: residue  Value: [rot, nrot, resn, atmn, wt_sc, hp, rot1, nrot1, resn1, atmn1, wt_sc1, hp1, score, score1]
    score -> nrot/atmn/maxASAsc
    score1 -> nrot1/atmn1/maxASAsc

    """
    try:
	data = open(inFile, 'r').readlines()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))
    mut_com = {}
    mut_com1 = {}
    for line in data:
#	pdbcode,ch,nid,rot, nrot, resn, atmn, wt_sc,rot1, nrot1, resn1, atmn1, wt_sc1,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,\
#	x16,x17,x18,x19,x20,x21,x22,hp,hp1,hp2 = line.split()
	resname,ch,nid,rot,nrot,resn,atmn,wt_sc,hp,rot1,nrot1,resn1,atmn1,wt_sc1,hp1,rot2,nrot2,resn2,atmn2,wt_sc2,hp2,\
		maxASAsc,ResMWsc,ratio3,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14 = line.split()
	if maxASAsc == '0.00':
		score = 100.00
		score1 = 100.00
	else:
		if atmn == '0.00':
			score = 100.00
		else:
			score = 100.0*float(nrot)/float(atmn)/float(maxASAsc)
		if atmn1 == '0.00':
			score1 = 100.00
		else:
			score1 = 100.0*float(nrot1)/float(atmn1)/float(maxASAsc)
	nscore = str(score)
	nscore1 = str(score1)
#	key = (pdbcode,ch,nid)
	key = (resname,ch,nid)
#	mut_com[key] = [rot, nrot, resn, atmn, wt_sc,hp, rot1, nrot1, resn1, atmn1, wt_sc1,hp1,nscore,nscore1]
	mut_com[key] = [rot1, nscore, nscore1]
	mut_com1[key] = [rot, rot1, hp1]
    return mut_com,mut_com1


def initNewMsf(inFile):
    """
    initNewMsf(inFile) -> residues' side chain msf divided by average

	inFile    (str)   residues' new msf file

	newmsf    (dict)  Key: residue  Value: newmsf[msf,msf2,msf3,pos1,pos2,pos3]

	"""
    try:
	data = open(inFile, 'r').readlines()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))
    newmsf = {}
    for line in data:
	pdbcode,ch,nid,msf,msf2,msf3,pos1,pos2,pos3 = line.split()
	key = (pdbcode,ch,nid)
	newmsf[key] = [msf,msf2,msf3,pos1,pos2,pos3]
    return newmsf
	

def initFeatures(inFile):
    """
    initFeatures(inFile) -> residue's position at interface
	
    inFile  (str)  residue's position file name
	
    features (dict) Key: residue  Value: [pos,msf1,msf2.msf3,ss]

    """
    try:
	data = open(inFile, 'r').readlines()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))
    features = {}
    for line in data:
	pdbcode,ch,nid,pos,msf1,msf2,msf3,ss = line.split()
	key = (pdbcode,ch,nid)		
	features[key] = [pos,msf1,msf2,msf3,ss]
    return features

def initStack(inFile):
    """
	initStack(inFile) -> pi-pi interaction of residues at interface
	
	inFile  (str) residues' pi-pi interaction feature file name
	
	stack   (dict) Key: residue   Value: continue

	"""
    try:
	data = open(inFile, "r").readlines()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))
    stack = {}
    for line in data:
	pdbcode,ch,nid,pipi = line.split()
	key = (pdbcode,ch,nid)
	stack[key] = [pipi,]
    return stack

def initExper(inFile):
    """
    initExper(inFile) -> if the residue is a hot spot

    inFile  (str)  experiment result file

    hotdef  (dict) Key: residue Value: hot [size, hot, resname, id]

    """
    try:
	data = open(inFile,'r').readlines()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))
    hotdef = {}
    for line in data:
	pdbcode,ch,nid,hot = line.split()
	key = (pdbcode,ch,nid)
	s_resname = nid[0]
	id = nid[1:]
	resname = residueTrans(s_resname)
	size = nresidueSize(resname)
	hotdef[key] = [size,hot,resname,id,pdbcode]
    return hotdef

#def toFileCoreRim(outFile,hotdef,corerim,sites):
#    """
#    Writes CoreRim formated data

#    outFile  (str)  Output file name
#    hotdef   (dict) Key: residue  Value: hot [size,hot,resname,id]
#    corerim  (dict) Key: residue  Value: pos

#   """
#    try:
#	output = open(outFile,'w')
#	for site in sites:
#	    sitedata = []
#	    sitedata.extend(hotdef[site])
#	    sitedata.extend(corerim[site])
#	    output.write(','.join(sitedata)+'.\n')
#	    output.write(','.join(sitedata))
#	    output.write(',')
#	    output.write(corerim[site])
#	    output.write('.\n')
#	output.close()
#    except IOError,io:
#	sys.exit('IOError: %s' % str(io))


def toFileKFADE(outFile,shape,fpoint,hotdef,corerim,stack,newsas,mut_com,contact,sites):
    """
    Writes K-FADE formated data (skips unknown residues)

    outFile  (str)  Output file name
    shape    (dict) Key: residue  Value: shape scores in 10 shells (list)
    fpoint   (dict) Key: residue  Value: FADE points in 10 shells (list) 
    hotdef   (dict) Key: residue  Value: hot [size,hot,resname,id]
    corerim  (dict) Key: residue  Value: features [pos_per,]
	stack    (dict) key: residue  Value: pipi
	newsas   (dict) Key: residue  Value: [delta1,per,delta2,ratio1,ratio2,ratio3,ratio4,ratio5]
	mut_com  (dict) Key: residue  Value: [rot, nrot, resn, atmn, wt_sc,hp, rot1, nrot1, resn1, atmn1, wt_sc1,hp1,score,score1]
	contact  (dict) Key: residue  Value: contacts [polar,generic,h-bond]
    sites    (list) Residues to evaluate

    """
    try:
	output = open(outFile,'w')
	for site in sites:
#	    c,id,aa = sitei
	    size,hot,resname,id,pdbcode = hotdef[site]
	    try:
#		sitedata = ['?',residueSize(aa)]
		sitedata = [residueHydrobicity(resname),]
	    except:
		print 'Note: unrecognized residue (%s %s %d), skipping' % \
		    (resname,id,pdbcode)
		continue
#	    sitedata = []
	    sitedata.extend(hotdef[site])
	    sitedata.extend(corerim[site])
	    sitedata.extend(newsas[site])
	    sitedata.extend(mut_com[site])
#	    sitedata.extend(newmsf[site])
	    sitedata.extend(shape[site])
	    sitedata.extend(fpoint[site])
	    sitedata.extend(stack[site])
	    try:
		sitedata.extend(contact[site])
	    except:
		sitedata.extend(['0','0','0.000'])
	    output.write(','.join(sitedata)+'.\n')
	output.close()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))

def toFileKCON(outFile,contact,fpoint,hotdef,corerim,stack,newsas,mut_com,sites):
    """
    Writes K-CON formated data (skips unknown residues)

    outFile (str)  Output file name
    contact (dict) Key: residue  Value: contacts [polar,generic,h-bond]
    fpoint  (dict) Key: residue  Value: FADE points in 10 shells (list) 
    hotdef   (dict) Key: residue  Value: hot [size,hot,resname,id,pdbcode]
    corerim  (dict) Key: residue  Value: features[pos_per,]
	stack dict key :residue Value: pipi[yes/no,]
	newsas   (dict) Key: residue  Value: [delta1,per,delta2,ratio1,ratio2,ratio3,ratio4,ratio5]
	mut_com  (dict) Key: residue  Value: [rot, nrot, resn, atmn, wt_sc,hp, rot1, nrot1, resn1, atmn1, wt_sc1,hp1,score,score1]
    sites   (list) Residues to evaluate

    """
    try:
	output = open(outFile,'w')
	for site in sites:
#	    c,id,aa = site
	    size,hot,resname,id,pdbcode = hotdef[site]
	    try:
#		sitedata = [residueType(resname),]
		sitedata = [residueHydrobicity(resname),]
	    except:
		print 'Note: unrecognized residue (%s %s %s), skipping' % \
		    (resname,id,pdbcode)
		continue
#	    sitedata = []
	    sitedata.extend(hotdef[site])
	    sitedata.extend(corerim[site])
#	    sitedata.extend(newmsf[site])
	    sitedata.extend(newsas[site])
	    sitedata.extend(mut_com[site])
	    try:
		sitedata.extend(contact[site])
	    except:
		sitedata.extend(['0','0','0.000'])
	    sitedata.extend(fpoint[site])
	    sitedata.extend(stack[site])
	    output.write(','.join(sitedata)+'.\n')
	output.close()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))

def toFileKfc2a(outFile,core_rim,newsas,mut_com,fpoint,sites):
    """
    Write kfc2a formated data

    outFile  (str)   Output file name
    core_rim (dict)  Key: residue   Value: [posper,]
    newsas   (dict)  Key: residue   Value: [delta1,per]
	mut_com  (dict)  Key: residue   Value: [rot1,score,score1]
	fpoint   (dict)  Key: residue   Value: [fp9,]
	sites    (list)  Residues to evaluate

    """
    try:
	output = open(outFile,'w')
	for site in sites:
	    resname,ch,nid = site
	    try:
		sitedata = [residueHydrobicity(resname),]
	    except:
		print 'Note: unrecognized residue (%s %s %s), skipping' % \
			(resname,ch,id)
		continue
	    sitedata.extend(core_rim[site])
	    sitedata.extend(newsas[site])
	    try:
		sitedata.extend(mut_com[site])
	    except:
		sitedata.extend(['10.0','100.0','100.0'])
	    try:
		sitedata.extend(fpoint[site])
	    except:
		sitedata.extend(['0',])
	    output.write(','.join(sitedata)+'.\n')
	output.close()
    except IOError,io:
	sys.exit('IOError: %s' % str(io))

def toFileKfc2b(outFile,core_rim,newsas1,mut_com1,fpoint1,sites):
    """
    Write kfc2a formated data

    outFile  (str)   Output file name
    core_rim (dict)  Key: residue   Value: [posper,]
    newsas1  (dict)  Key: residue   Value: [delta1,per]
    mut_com1 (dict)  Key: residue   Value: [rot1,score,score1]
    fpoint1  (dict)  Key: residue   Value: [fp9,]
    sites    (list)  Residues to evaluate

    """
    try:
	output = open(outFile,'w')
	for site in sites:
	    resname,ch,nid = site
	    try:
		sitedata = [nresidueSize(resname),]
	    except:
		print 'Note: unrecognized residue (%s %s %s), skipping' % \
			(resname,ch,id)
		continue
	    sitedata.extend(core_rim[site])
	    sitedata.extend(newsas1[site])
	    try:
		sitedata.extend(mut_com1[site])
	    except:
		sitedata.extend(['10.0','10.0','-10.0'])
	    try:
		sitedata.extend(fpoint1[site])
	    except:
		sitedata.extend(['0',])
	    output.write(','.join(sitedata)+'.\n')
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
#    python FormatData.py hotfile featurefile shapefile contactfile stackfile newmsffile out-nkfade out-nkcon
#	python FormatData.py hotfile shapefile contactfile stackfile newsasfile core_rimfile mut_comfile out-nkfade out-nkcon
	python FormatData.py core_rim.txt corerim.txt compensation.txt file.fnb outForKfc2a outForKfc2b

    """
    if argv == None:
	argv = sys.argv
#    if len(argv) != 10:
#	sys.exit('Usage: python %s hotfile shape contact stackfile newmsf core_rim out-kfade out-kcon' \
#		 % argv[0])
    if len(argv) != 8:
	sys.exit('Usage: python %s core_rim.txt corerim.txt compensation.txt file.fnb outForKfc2a outForKfc2b outForInter' \
			% argv[0])

    core_rim=initCoreRim(argv[1])
    newsas,newsas1 = initNewSas(argv[2])
    mut_com,mut_com1 = initMutCom(argv[3])
    fpoint,fpoint1 = initShapeEnvironment(argv[4])
    sites = core_rim.keys()
    sites.sort()
    toFileKfc2a(argv[5],core_rim,newsas,mut_com,fpoint,sites)
    toFileKfc2b(argv[6],core_rim,newsas1,mut_com1,fpoint1,sites)
    toFileSites(argv[7],sites)

#    shape,fpoint = initShapeEnvironment(argv[2])
#    contact = initBiochemContacts(argv[3])
#    hotdef = initExper(argv[1])
##    nfeatures = initFeatures(argv[2])
#    nstack = initStack(argv[4])
##    newmsf = initNewMsf(argv[5])
#    newsas = initNewSas(argv[5])
#    mut_com = initMutCom(argv[7])
#    core_rim=initCoreRim(argv[6])
#    sites = hotdef.keys()
#    sites.sort()
##    toFileCoreRim(argv[3],hotdef,corerim,sites)
##    toFileKFADE(argv[7],shape,fpoint,hotdef,core_rim,nstack,newmsf,sites)
#    toFileKFADE(argv[8],shape,fpoint,hotdef,core_rim,nstack,newsas,mut_com,contact,sites)
##    toFileKCON(argv[9],contact,fpoint,hotdef,core_rim,nstack,newmsf,sites)
#    toFileKCON(argv[9],contact,fpoint,hotdef,core_rim,nstack,newsas,mut_com,sites)
##    toFileSites(argv[5],sites)

if __name__ == '__main__':
    main()

"""
gw 3/2011 - keep data from auxillary sources even if the chain has no KFC hits

Formats predictions from multiple sources
Formats the predictions from KFC, ConSurf, Rosetta, and
experimental data into a final results file.

"""

# File:    CombineResults.py
# Sources: None
#
# Authors:
#
# Steven Darnell / sdarnell@biochem.wisc.edu / UW-Madison
# Julie Mitchell / mitchell@math.wisc.edu    / UW-Madison
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

def makeKey(chn,res,num):
    try:
        key = (chn,int(num),res.upper())
    except ValueError:
        key = None
    return key

def inv_Key( key ):
    if key != None:
       ( chn,num,res) = key
       return ( chn, str(num), res )
    else:
       return 0


def loadResults(file,type):
    d = {}
    try:
        if type.lower() == 'consurf':
            data = open(file,'r').readlines()

            for ln in data:
                try:
                    fields = ln.split()
                    l3atom = fields[2]
                    res = l3atom[:3]
                    num = l3atom[3:-2]
                    chn = l3atom[-1]
                    key = makeKey(chn,res,num)

                    color = fields[4]
                    if not color.isdigit():
                        color = color[:-1]
                    if int(color) >= 7:
                        d[key] = ('Conserv',color)
                    else:
                        d[key] = ('-------',color)

                except:
                    continue

        elif type.lower() == 'rosetta':
            iToRes = {1:'ALA', 2:'CYS', 3:'ASP', 4:'GLU', 5:'PHE', 6:'GLY',
                      7:'HIS', 8:'ILE', 9:'LYS', 10:'LEU', 11:'MET', 12:'ASN',
                      13:'PRO', 14:'GLN', 15:'ARG', 16:'SER', 17:'THR',
                      18:'VAL', 19:'TRP', 20:'TYR'}
            data = open(file,'r').readlines()[3:]

            for ln in data:
                try:
                    fields = ln.split()
                    num = fields[0]
                    chn = fields[1]
                    res = iToRes[int(fields[4])]
                    key = makeKey(chn,res,num)

                    ddg = fields[5]
                    if float(ddg) >= 2.0:
                        d[key] = ('Hotspot',ddg)
                    else:
                        d[key] = ('-------',ddg)

                except:
                    continue

        elif type.lower() == 'experimental':
            data = open(file,'r').readlines()
            for ln in data:
                try:
                    fields = ln.split()
                    chn = fields[0]
                    res = fields[1]
                    num = fields[2]
                    key = makeKey(chn,res,num)

                    hot,val = int(fields[3]),' '.join(fields[4:])
                    if hot:
                        d[key] = ('Hotspot',val)
                    else:
                        d[key] = ('-------',val)

                except:
                    continue
                    
        else:
            print('Warning: file type \'%s\' not recognized' % type)
            return 0

    except IOError:
        print('Note: %s file not found, skipping file' % type)
        return 0

    return d

def readData(data,key):
    if key in data:
        return data[key]
    else:
        return ('-------','---')

def main(argv=None):
    if argv == None:
        argv = sys.argv
    if len(argv) != 2:
        sys.exit('Usage: %s jobid' % argv[0])

    jobid = argv[1]


    combine = open(jobid+'.kfc.results.temp','w')
    results = open(jobid+'.kfc.results','r').readlines()

    # write header to combined file...
    combine.writelines(results[:3])
    ln = results[3][:-1]
    combine.write(ln+'  ConSurf ConSu  Rosetta  Roset   Exper  Exper\n')
    ln = results[4][:-1]
    combine.write(ln+'   Class  Value   Class    DDG    Class  Value\n')
    ln = results[5][:-1]
    combine.write(ln+'------------------------------------------------ \n')

    chainkeyset = {}

    # read through kfc2 results and get all keys.
    fileid = 0
    kfchitlines = {}
    for ln in results[6:]:
        ln = ln[:-1]
        chn,res,num = ln.split()[:3]
        num = int(num)
        if chn not in chainkeyset:
           chainkeyset[chn] = {}
        if num not in chainkeyset[chn]:
           chainkeyset[chn][num] = [(res,fileid)]
        else:
           #check for case when same chain/resnum occurs more than once for different resnames!
           for hits in chainkeyset[chn][num]:
              ( rs, fid ) = hits 
              if res != rs:
                  print "WARNING, same chain/resnum occurs more than once for different resnames!"
                  print "   check:",  chn,rs,num, "fileid",fid
           chainkeyset[chn][num].append( (res,fileid) )

        key = makeKey(chn,res,num)
        kfchitlines[key] = ln

    consurf = loadResults(jobid+'.consurf','Consurf')
    if consurf != 0:
       # read through consurf results, if present, and get all keys.
       fileid = 1
       for key in consurf:
           ikey = inv_Key( key )
           if ikey == 0: continue
           (chn,num,res) = ikey
           num = int(num)
           if chn not in chainkeyset:
              chainkeyset[chn] = {}
           if num not in chainkeyset[chn]:
              chainkeyset[chn][num] = [(res,fileid)]
           else:
              #check for case when same chain/resnum occurs more than once for different resnames!
              for hits in chainkeyset[chn][num]:
                 ( rs, fid ) = hits 
                 if res != rs:
                     print "WARNING, same chain/resnum occurs more than once for different resnames!"
                     print "   check:",  chn,rs,num, "fileid",fid
              chainkeyset[chn][num].append( (res,fileid) )
    else:
       consurf = [None]
           

    rosetta = loadResults(jobid+'.rosetta','Rosetta')
    if rosetta != 0:
       # read through rosetta results, if present, and get all keys.
       fileid = 2
       for key in rosetta:
           ikey = inv_Key( key )
           if ikey == 0: continue
           (chn,num,res) = ikey
           num = int(num)
           if chn not in chainkeyset:
              chainkeyset[chn] = {}
           if num not in chainkeyset[chn]:
              chainkeyset[chn][num] = [(res,fileid)]
           else:
              #check for case when same chain/resnum occurs more than once for different resnames!
              for hits in chainkeyset[chn][num]:
                 ( rs, fid ) = hits 
                 if res != rs:
                     print "WARNING, same chain/resnum occurs more than once for different resnames!"
                     print "   check:",  chn,rs,num, "fileid",fid
              chainkeyset[chn][num].append( (res,fileid) )
    else:
       rosetta = [None]


    expermt = loadResults(jobid+'.experimental','Experimental')
    if expermt != 0:
       # read through expermt results, if present, and get all keys.
       fileid = 3
       for key in expermt:
           ikey = inv_Key( key )
           if ikey == 0: continue
           (chn,num,res) = ikey
           num = int(num)
           if chn not in chainkeyset:
              chainkeyset[chn] = {}
           if num not in chainkeyset[chn]:
              chainkeyset[chn][num] = [(res,fileid)]
           else:
              #check for case when same chain/resnum occurs more than once for different resnames!
              for hits in chainkeyset[chn][num]:
                 ( rs, fid ) = hits 
                 if res != rs:
                     print "WARNING, same chain/resnum occurs more than once for different resnames!"
                     print "   check:",  chn,rs,num, "fileid",fid
              chainkeyset[chn][num].append( (res,fileid) )
    else:
        expermt = [None]



#    C      ALA    419  -------  -2.84  -------   -0.78
#012345678901234567890123456789012345678901234567890123456789

    chainkeyset_keys = chainkeyset.keys()
    chainkeyset_keys.sort()

    for k in chainkeyset_keys:
        resn = chainkeyset[k]
        resnkeys = resn.keys()
        resnkeys.sort()
        for rn in resnkeys:
           print k, rn
           if len( chainkeyset[k][rn] ) == 1:
              if chainkeyset[k][rn][0][1] == 1:
                 print 'SKIP consurf only residue:',rn,k
                 continue

           for (resname, fileid) in  chainkeyset[k][rn]:
              print '    ',resname, fileid

           key = makeKey(k,resname,str(rn))

           ln = "    %1s      %3s %6d" %(k, resname,rn )
           if key in kfchitlines: 
              ln += kfchitlines[key][21:].rstrip() + '  '
           else:
              ln += '  -------  -----  -------   -----  '

           chot,cval = readData(consurf,key)
           rhot,rval = readData(rosetta,key)
           ehot,eval = readData(expermt,key)
           print >> combine, ln + ' %7s %5s  %7s %5s  %7s  %s ' %(chot,cval,rhot,rval,ehot,eval)

    combine.close()


if __name__ == '__main__':
    main()

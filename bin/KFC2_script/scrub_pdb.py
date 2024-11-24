#!/usr/bin/env python

# gw 2011


# Example 1:  From job 534 pdb -----------------------------------------

# before: ATOM      6  CG  ASP L   1       6.253  12.635  14.396  1.00 99.99
# before: ATOM      7 OD1  ASP L   1       6.498  12.782  15.604  1.00 99.99
# before: ATOM      8 OD2  ASP L   1       5.176  12.173  13.936  1.00 99.99
# before: ATOM      9   N  ILE L   2      10.370  13.020  12.156  1.00 99.99
#         012345678901234567890123456789012345678901234567890123456789012345678
# after:  ATOM      6  CG  ASP L   1       6.253  12.635  14.396  1.00 99.99
# after:  ATOM      7  OD1 ASP L   1       6.498  12.782  15.604  1.00 99.99 
# after:  ATOM      8  OD2 ASP L   1       5.176  12.173  13.936  1.00 99.99 
# after:  ATOM      9  N   ILE L   2      10.370  13.020  12.156  1.00 99.99 

# Example 2:  From job 929 pdb -----------------------------------------

# before: ATOM      1  N   MET A     1     -26.005  -1.070 -12.128  1.00 27.19           N
# before: ATOM      2  CA  MET A     1     -25.700   0.377 -12.066  1.00 26.69           C
# before: ATOM      3  C   MET A     1     -24.536   0.629 -11.118  1.00 25.53           C
# before: ATOM      4  O   MET A     1     -24.034  -0.298 -10.486  1.00 25.11           O
#         012345678901234567890123456789012345678901234567890123456789012345678
# after:  ATOM      1  N   MET A   1     -26.005  -1.070 -12.128  1.00 27.19
# after:  ATOM      2  CA  MET A   1     -25.700   0.377 -12.066  1.00 26.69
# after:  ATOM      3  C   MET A   1     -24.536   0.629 -11.118  1.00 25.53
# after:  ATOM      4  O   MET A   1     -24.034  -0.298 -10.486  1.00 25.11

import sys, os
filename = sys.argv[1]
if not os.path.isfile( filename ):
   print "File not found in scrub_pdb.py:", filename
   sys.exit(1)

fh = open(filename, 'r')
lines = fh.readlines()
fh.close()
if len(sys.argv) < 3:
   print "Not emough command line args given scrub_pdb.py."
   sys.exit(1)
  
try:
   fhnew = open(sys.argv[2], 'w' )
except:
   print "Can not open write file in scrub_pdb.py:", sys.argv[2]
   sys.exit(1)
 
std_aa = ["SER","LEU","GLY","ALA","LYS","VAL","THR","ASP","GLU",
          "PRO","ASN","ARG","PHE","ILE","GLN","TYR","CYS","HIS",
          "MET","TRP",'DG', 'DT', 'DC', 'DA', 'DI', 'A', 'G', 'C', 'T', 'U', 'I']
sub_aa = { "MSE":"MET", "HID":"HIS", "HIE":"HIS", "HSE":"HIS", "HSD":"HIS", "CYX":"CYS" }
sub_aa_know = sub_aa.keys()
renamed_aa = []
skipped_aa = []

lnum = 1

# must assume chain label is always in l[21] but we want to allow
# for empty chain labels...  (5/6/11)


# 5/13/2011 (Stop After First MODEL...)
HitModel = False

for l in lines:
   l = l.replace(chr(0),'').rstrip()
   if l[0:5] == 'ATOM ':

      lchain_to_end = l[21:]
      w = lchain_to_end.split()
      # 4 digit resnum will result in no space between chain and num, like A2001
      # with iCode possible to get A2001T
      # also, might have negative resnumbers!

      if l[21] == ' ':
         chain = ' '
         resnum = w[0]
      else:
         chain = w[0]
         if len(chain) > 1:
            chain = w[0][0]
            resnum = w[0][1:]
         else:
            resnum = w[1]

      if resnum.isdigit():
         iCode = " "
      else: # check for - res number, leave as neg but logic for icode needs check
         if resnum[0] == '-':
            if resnum[1:].isdigit():
               iCode = " "
            else:
               if resnum[1:-1].isdigit():
                  iCode = resnum[-1]
                  resnum = resnum[:-1]
               else:
                  print "ERROR: Bad Atom record in PDB, -resNum/icode on line number:", lnum
                  print l
                  sys.exit(1)
         else:
            if resnum[0:-1].isdigit():
               iCode = resnum[-1]
               resnum = resnum[:-1]
            else:
               print "ERROR: Bad Atom record in PDB, +resNum/icode on line number:", lnum
               print l
               sys.exit(1)


      # Assume coords on out are formatted OK, but in some cases might be shifted to the
      # right (or left?) as in case ii in top comment above...
      # Look for decimal pt in x at column 35, should be first decimal pt lchain_to_end,
      # use this position to reegister rest of line...
      # (go back 4 to left of decimal

      icoords = lchain_to_end.find('.') - 4

      linepart = " %1s%4s%1s   %s" %(chain,resnum,iCode,lchain_to_end[icoords:])

      anam = l[12:16].strip()
      len_anam = len(anam)
      if len_anam < 4:
         a = " %-3s" %anam
      else:
         a = "%-4s" %anam

      rnam = l[17:20]
      rnamcheck = l[17:20].strip()
      if rnamcheck not in std_aa:
         if rnamcheck not in sub_aa_know:
            if rnamcheck not in skipped_aa:
               skipped_aa.append(rnam)
            continue
         else:
            if rnamcheck not in renamed_aa:
               renamed_aa.append(rnam)
            rnam = sub_aa[rnam]

      lnew = l[0:12] + a + l[16] + rnam + linepart
      print >> fhnew, lnew

   elif l[0:6] == 'HETATM':
      a = "%-4s" %l[12:16].strip()
      lnew = l[0:12] + a + l[16:]
      print >> fhnew, lnew

   elif l[0:6] == 'MODEL ':      
      if HitModel: break
      HitModel = True

   else:
      print >> fhnew, l

   lnum += 1

fhnew.close()

if len(renamed_aa) > 0:
   print "WARNING: Some residue names were changed to a standard name:", renamed_aa
if len(skipped_aa) > 0:
   print "WARNING: Some residues with unrecognized names were SKIPPED:", skipped_aa


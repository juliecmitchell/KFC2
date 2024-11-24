#!/usr/bin/env python

import sys
import os

stripped_pdb = sys.argv[1]

# PREFERABLY, THIS PDB FILE SHOULD ONLY CONSIST OF THE CHAINS BEING USED IN CURRENT KFC JOB.
# That is, the sum of the KFC script .mol1 and .mol2. files.
# (I did add chain check below, just in case..)

if not os.path.isfile(stripped_pdb):
   print "ERROR: File not found:", stripped_pdb
   sys.exit(1)

hetatms_pdb = sys.argv[2]
if not os.path.isfile(hetatms_pdb):
   print "ERROR: File not found:", hetatms_pdb
   sys.exit(1)

fh = open(stripped_pdb, 'r')
alines = fh.readlines()
fh.close()

fh = open(hetatms_pdb, 'r')
hlines = fh.readlines()
fh.close()

chainsA = sys.argv[3]
chainsB = sys.argv[4]

chains = chainsA + chainsB

try:
   fh = open(hetatms_pdb+'.pdb', 'w')
except:
   print "ERROR: cannot open output HETATMS keeper file."
   sys.exit(1)

# Use chain labels for HETATMS if available, keep no-matchers for additional testing...
hetsnotinchains = []
for l in hlines:
   if l[0:6] != 'HETATM': continue
   if l[21] in chains:
      print >> fh, l.strip()
   else: 
      hetsnotinchains.append(l)
del hlines
if len( hetsnotinchains ) == 0:
   fh.close()
   sys.exit(0)

# Store all coords for polymer chains we're using in this job...
ax = []
ay = []
az = []
for l in alines:
   if l[0:5] != 'ATOM ': continue
   if l[21] not in chains: continue
   ax.append( float( l[30:38] ) )
   ay.append( float( l[38:46] ) )
   az.append( float( l[46:54] ) )
del alines

dcut = 4.00
dcut2 = dcut**2
keephet = []
for l in hetsnotinchains:
   x = float( l[30:38] )
   y = float( l[38:46] )
   z = float( l[46:54] )
   ia = -1
   keep = False
   for xa in ax:
      ia += 1   
      dx = abs(xa - x) 
      if dx < dcut:
         dy = abs(ay[ia] - y)
         if dy < dcut:
            dz = abs(az[ia] - z)
            if dz < dcut:
               d2 = dx**2 + dy**2
               if d2 < dcut2:
                  d2 += dz**2
                  if d2 <= dcut2:
                     keep = True
                     break

   if keep:
      hetresnamechainnumber = l[17:26]
      print "keeper!", hetresnamechainnumber 
      if hetresnamechainnumber not in keephet:
         keephet.append( hetresnamechainnumber )


for l in hetsnotinchains:
   hetresnamechainnumber = l[17:26]
   if hetresnamechainnumber in keephet:
      print >> fh, l.strip()

fh.close()

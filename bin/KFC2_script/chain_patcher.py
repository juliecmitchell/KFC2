#!/usr/bin/env python

# gw 2011, use this for case when no chain sets were given by user,
#          need to define unique labels for both sets
#          allow cases like:   A TER A B C,  need to relabel second "A" to "a" or some
#                              other label if "a" is already in use
#          another case, unlabled chains:   " " TER " " 
#                              assign label to first chain and another label to everything after
#          another case: A B TER unlabled chain(s)
#                            assign unique label to unlabeled

import sys, os, string


jmolchaincolor = { "A": "C0D0FF", "B": "B0FFB0", "C": "FFC0C8", "D": "FFFF80", "E": "FFC0FF", "F": "B0F0F0",
                   "G": "FFD070", "H": "F08080", "I": "F5DEB3", "J": "00BFFF", "K": "CD5C5C", "L": "66CDAA",
                   "M": "9ACD32", "N": "EE82EE", "O": "00CED1", "P": "00FF7F", "Q": "3CB371", "R": "00008B",
                   "S": "BDB76B", "T": "006400", "U": "800000", "V": "808000", "W": "800080", "X": "008080",
                   "Y": "B8860B", "Z": "B22222", "0": "00FF7F", "1": "3CB371", "2": "00008B", "3": "BDB76B",
                   "4": "006400", "5": "800000", "6": "808000", "7": "800080", "8": "008080", "9": "B8860B"  } 

defaultcol = "FFFFFF"	 


def jmolcolspanhtml( chain ):
   # Jmol colors lower/upper the same.
   upchain = chain.upper()
   if upchain in jmolchaincolor:
      col = jmolchaincolor[upchain]
   else:
      col = defaultcol
   colspanhtmlle = '<span style="background-color:#' + col + '">       </span>'
   return  colspanhtmlle

def html_top( jobid, fhtml ):
   print >> fhtml, '<head>'
   print >> fhtml, '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />'
   print >> fhtml, '<title>Job ' + str(jobid) + ' Chain Names</title>'
   print >> fhtml, '</head>'
   print >> fhtml, '<body>'
   print >> fhtml, '<pre>'

def html_bot( fhtml ):
   print >> fhtml, '</body>'
   print >> fhtml, '</pre><br><br>'
   print >> fhtml, '</html>'

def get_unique_label_reverse( chains ):
   # first, try uppercase from Z-->A:
   for c in string.uppercase[::-1]:
      if c not in chains: return c

   # next, try lowercase from z-->a:
   for c in string.lowercase[::-1]:
      if c not in chains: return c

   # next, try digits from 9-->0:
   for c in string.digits[::-1]:
      if c not in chains: return c

   # ok, we're getting desperate...
   for c in string.punctuation:
      if c not in chains: return c

   print "Error in chain_patcher.py: ran out of unique chain labels!"
   sys.exit(1)

def get_unique_label( chains ):
   # first, try uppercase from A-->Z:
   for c in string.uppercase:
      if c not in chains: return c

   # next, try lowercase from a-->z:
   for c in string.lowercase:
      if c not in chains: return c

   # next, try digits from 0-->9:
   for c in string.digits:
      if c not in chains: return c

   # ok, we're getting desperate...
   for c in string.punctuation:
      if c not in chains: return c

   print "Error in chain_patcher.py: ran out of unique chain labels!"
   sys.exit(1)


filename = sys.argv[1]
if not os.path.isfile( filename ):
   print "File not found in chain_patcher.py:", filename
   sys.exit(1)

fh = open(filename, 'r')
lines = fh.readlines()
fh.close()

if len(sys.argv) < 4:
   print "Not enough command line args given chain_patcher.py."
   sys.exit(1)
  
try:
   fhnew = open(sys.argv[2], 'w' )
except:
   print "Can not open write file in chain_patcher.py:", sys.argv[2]
   sys.exit(1)

jobid = sys.argv[3]

chains = []
#
# Get list list of chain names already in use in initial PDB file.
# ( if needed, will generate new names that are not laready in use, must update/check this list.) 

Nters = 0
LastTer = False
for l in lines:
   l = l.rstrip()
   if l[0:5] == 'ATOM ' or l[0:7] == 'HETATM ':
      LastTer = False
      chain = l[21]
      if chain not in chains:
         chains.append(chain)
   else:
      if l[0:3] == 'TER':
         if LastTer: continue  # skip multiple TERs in a row...
         Nters += 1
         LastTer = True

if len(chains) < 1:
   print "ERROR: No Atoms or HetAtms found in PDB file."
   sys.exit(1)

if (Nters == 0) or (LastTer and Nters == 1):
   if len(chains) == 1 and chain == ' ':
      print "ERROR: No chain labels and no TER between Atom found in PDB file, cannot separate."
      sys.exit(1)

   print "USING FIRST CHAIN (NOT TERS)."
   UseFirstChain = True
else:
   UseFirstChain = False
   print "FOUND NON-TERMINAL TER."


print "Chains used in original PDB file:", chains
resranges = {}
hitFirstTER = False
haveNewChainLabel = False
blank_chains_newnames = [] #
chain_set1 = []
chain_set2 = []
secondset_remapping = {}  # for case like A TER A TER A TER,  map both A's after first TER to "a" (OK?)
original = {}
lastchain = ''
chainprevious = '00'
rawchainprevious = '00'
resnumlast = '-99999'
chainset1 = ''
chainset2 = ''
chain_has_hetatms = {}
chain_has_atoms = {}
NRENAMED = 0

for l in lines:
   l = l.rstrip()
   if l[0:5] == 'ATOM ' or l[0:7] == 'HETATM ':

      if l[0] == 'H':
         HET = True
         ATM = False
      else:
         HET = False
         ATM = True

      resnum = l[22:26]
      chain = l[21]

      if chain != rawchainprevious and rawchainprevious != '00':
          if UseFirstChain:
             #resranges[ chainprevious ].append( (resnum1,resnumlast) )
             #resnum1 = resnum
             if chainprevious not in chain_set1:
                chain_set1.append( chainprevious )
                ###resranges[ chain ] = []
             haveNewChainLabel = False
             hitFirstTER = True
             print >> fhnew, 'TER'
             UseFirstChain = False

      rawchainprevious = chain

      if chain == ' ':
         if not haveNewChainLabel:
            # Create a label for this unlabled chain...
            chain = get_unique_label( chains )
            NRENAMED += 1
            chains.append( chain ) 
            original[chain] = ' '
            resranges[ chain ] = []
            resnum1 = resnum
            blank_chains_newnames.append( chain )
            chainlast = chain
            haveNewChainLabel = True
            # Initialize...
            chain_has_hetatms[chain] = False
            chain_has_atoms[chain] = False

         else:
            # Already creates a label for this unlabled chain...
            chain = chainlast
      else:
         haveNewChainLabel = False
         if hitFirstTER:

            if chain in chain_set1:

               # Chain label in set2 is already in use by set1. Gotta change label.

               if chain in secondset_remapping:

                  # Already have a new label for this labeled chain, re-use it. ( AB_TER_ABABAB OK?)

                  chain = secondset_remapping[chain]

               else:

                  # Need to generate new unique label...

                  if chain.isupper():
                     chainuse = chain.lower()
                     if chainuse in chains:
                        chainuse = get_unique_label_reverse( chains )
                  else:
                     chainuse = get_unique_label_reverse( chains )

                  # Save a bunch of stuff...
                  chains.append( chainuse )
                  secondset_remapping[chain] = chainuse
                  original[chainuse] = chain
                  chain = chainuse
                  NRENAMED += 1
                  # Initialize...
                  chain_has_hetatms[chain] = False
                  chain_has_atoms[chain] = False

            else:
               # Chain is already unique...
               if chain not in original:
                  original[chain] = chain
                  chain_has_hetatms[chain] = False
                  chain_has_atoms[chain] = False

         else:
            # labeled chains before any TER, leave labels alone?  What about ABAB-TER-CDEF? OK?  
            # Note that this chain in already in chains since we grabbed all
            # chains in first pass thru file above
            if chain not in original:
               original[chain] = chain
               chain_has_hetatms[chain] = False
               chain_has_atoms[chain] = False

      if HET:
         chain_has_hetatms[chain] = True
      if ATM:
         chain_has_atoms[chain] = True
 
      if hitFirstTER:
         if chain not in chain_set2:
            chain_set2.append( chain )
            resranges[ chain ] = []
            resnum1 = resnum
      else:
         if chain not in chain_set1:
            chain_set1.append( chain )
            resranges[ chain ] = []
            resnum1 = resnum

      if chain != chainprevious and chainprevious != '00':
         resranges[ chainprevious ].append( (resnum1,resnumlast) )
         resnum1 = resnum

      chainprevious = chain
      resnumlast = resnum

      print >> fhnew, ''.join( ( l[0:21], chain, l[22:] ) )

   elif l[0:3] == 'TER':
      resranges[ chain ].append( (resnum1,resnumlast) )
      resnumlast = '-99999'
      chainprevious = '00'
      haveNewChainLabel = False
      hitFirstTER = True
      print >> fhnew, 'TER'

print >> fhnew, 'END'
fhnew.close()

if resnumlast != '-99999':
    resranges[ chain ].append( (resnum1,resnumlast) )

try:
   fhtml = open(jobid + "_chains.html", "w" )
except:
   print "Warning: Can not open write file in chain_patcher.py: reset_chain_label_sets.html"
   sys.exit(0)

nth_blnk = 0

chain_set1_not_all_het = []
chain_set2_not_all_het = []

print 'chainorig:', original
print 'chain_set1:',  chain_set1
print 'chain_set2:',  chain_set2

print 'Number Renamed:', NRENAMED 

html_top( jobid, fhtml )
print >> fhtml, "Chain Set 1:"
print >> fhtml, 112*'-'
for chain in chain_set1:

   if chain_has_hetatms[chain]:
      if chain_has_atoms[chain]:
         extra = "ATOMS and HETATMS"
         chain_set1_not_all_het.append(chain)
      else:
         extra = "HETATMS only"
   else:
      if chain_has_atoms[chain]:
         extra = "ATOMS only"
         chain_set1_not_all_het.append(chain)
      else:
         extra = "NOTHING?"

   chainorig = original[chain]
   '''
   if chainorig != ' ':
      newchain = chain
   else:
      newchain = blank_chains_newnames[nth_blnk]
      nth_blnk += 1
   '''

   colblock = jmolcolspanhtml( chain )

   for respair in resranges[chain]:
      (res1, res2) = respair
      res1 = res1.rstrip()
      res2 = res2.lstrip()
      if chain != chainorig:
         print >> fhtml,  colblock + '  Original Chain: "' + chainorig + '"   Renamed To: "' + chain + '"   Residue Range:', \
                          res1 + ' - ' + res2, \
                         '  Contains:', extra
      else:
         print >> fhtml,  colblock + '  Original Chain: "' + chainorig + '"   Not Renamed.      Residue Range:', \
                          res1 + ' - ' + res2, \
                         '  Contains:', extra



print >> fhtml, ' '
print >> fhtml, ' '
print >> fhtml, "Chain Set 2:"
print >> fhtml, 112*'-'
for chain in chain_set2:

   if chain_has_hetatms[chain]:
      if chain_has_atoms[chain]:
         extra = "ATOMS and HETATMS"
         chain_set2_not_all_het.append(chain)
      else:
         extra = "HETATMS only"
   else:
      if chain_has_atoms[chain]:
         extra = "ATOMS only"
         chain_set2_not_all_het.append(chain)
      else:
         extra = "NOTHING?"


   colblock = jmolcolspanhtml( chain )

   chainorig = original[chain]

   for respair in resranges[chain]:
      (res1, res2) = respair
      res1 = res1.rstrip()
      res2 = res2.lstrip()
      if chain != chainorig:
         print >> fhtml, colblock + '  Original Chain: "' + chainorig + '"   Renamed To: "' + chain + '"   Residue Range:', \
                          res1 + ' - ' + res2, \
                         '  Contains:', extra
      else:
         print >> fhtml,  colblock + '  Original Chain: "' + chainorig + '"   Not Renamed.      Residue Range:', \
                          res1 + ' - ' + res2, \
                         '  Contains:', extra
html_bot( fhtml )
fhtml.close()

try:
   fhsh = open("reset_chain_label_sets.sh", "w" )
except:
   print "ERROR: Can not open write file in chain_patcher.py: reset_chain_label_sets.sh"
   sys.exit(1)

print >> fhsh, 'C1=' + ''.join( chain_set1_not_all_het )
print >> fhsh, 'C2=' + ''.join( chain_set2_not_all_het )

fhsh.close()

sys.exit(0)





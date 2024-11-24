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


def resint(s):
   if s[-1] not in string.digits:
      return int(s[0:-1])
   return int(s)      


def rescheck( resnum1, resnum2 ):
   print ' resnum1, resnum2:',  resnum1, resnum2
   if resnum1[-1] not in string.digits:
      # Assume that if second residue in a series is a specially labelled one and that
      # (should be checked in calling test statement) same sequence label, then this new
      # residue is in same chain...
      return True
   else:
      ir1 = int(resnum1)

   if resnum2[-1] not in string.digits:
      # Assume that if second residue in a series is a specially labelled one and that
      # (should be checked in calling test statement) same sequence label, then this new
      # residue is in same chain...
      return True
   else:
      ir2 = int(resnum2)

   if ir2 >= ir1:
      return True
   else:
      return False


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

def html_bot( fhtml ):
   print >> fhtml, '</body>'
   print >> fhtml, '</pre><br><br>'
   print >> fhtml, '</html>'

def get_unique_label_reverse_strict( chains ):
   # first, try lowercase from z-->a:
   for c in string.lowercase[::-1]:
      if c not in chains and c.upper() not in chains: return c

   # next, try uppercase from Z-->A:
   for c in string.uppercase[::-1]:
      if c not in chains and c.lower() not in chains: return c

   # next, try digits from 9-->0:
   for c in string.digits[::-1]:
      if c not in chains: return c

   # ok, we're getting desperate...
   for c in string.punctuation:
      if c not in chains: return c

   return ' '

def get_unique_label_reverse( chains ):
   # first, try lowercase from z-->a:
   for c in string.lowercase[::-1]:
      if c not in chains: return c

   # next, try uppercase from Z-->A:
   for c in string.uppercase[::-1]:
      if c not in chains: return c

   # next, try digits from 9-->0:
   for c in string.digits[::-1]:
      if c not in chains: return c

   # ok, we're getting desperate...
   for c in string.punctuation:
      if c not in chains: return c

   return ' '

def get_unique_label_strict( chains ):

   # first, try lowercase from a-->z:
   for c in string.lowercase:
      if c not in chains and c.upper() not in chains: return c

   # next, try uppercase from A-->Z:
   for c in string.uppercase:
      if c not in chains and c.lower() not in chains: return c

   # next, try digits from 0-->9:
   for c in string.digits:
      if c not in chains: return c

   # ok, we're getting desperate...
   for c in string.punctuation:
      if c not in chains: return c

   return ' '

def get_unique_label( chains ):

   # first, try lowercase from a-->z:
   for c in string.lowercase:
      if c not in chains: return c

   # next, try uppercase from A-->Z:
   for c in string.uppercase:
      if c not in chains: return c

   # next, try digits from 0-->9:
   for c in string.digits:
      if c not in chains: return c

   # ok, we're getting desperate...
   for c in string.punctuation:
      if c not in chains: return c

   return ' '

def unique_label_reverse( chains):
   chain = get_unique_label_reverse_strict( chains )
   if chain == ' ':
      chain = get_unique_label_reverse( chains )
      if chain == ' ':
         print "Unable to generate unique chain label for chain."
         chain = string.punctuation[-1]
   return chain

def unique_label( chains):
   chain = get_unique_label_strict( chains )
   if chain == ' ':
      chain = get_unique_label( chains )
      if chain == ' ':
         print "Unable to generate unique chain label for chain."
         chain = string.punctuation[-1]
   return chain

filename = sys.argv[1]
if not os.path.isfile( filename ):
   print "ERROR: Input PDB file not found in chain_patcher_general.py:", filename
   sys.exit(1)

fh = open(filename, 'r')
lines = fh.readlines()
fh.close()

if len(sys.argv) < 6:
   print "ERROR: Not enough command line args given chain_patcher_general.py."
   print "ERROR: Usage: chain_patcher_general.py PDBfilein PDBfileout chain_set_1 chain_set_0 JobID"
   sys.exit(1)
  
try:
   fhnew = open(sys.argv[2], 'w' )
except:
   print "ERROR: Can not open write file in chain_patcher_general.py:", sys.argv[2]
   sys.exit(1)

in_chain_set_1 = sys.argv[3]
in_chain_set_2 = sys.argv[4]

if in_chain_set_1 == '_' and in_chain_set_2 == '_':
   do_not_rename_these_chains = ''
else:
   do_not_rename_these_chains = in_chain_set_1 + in_chain_set_2

jobid = sys.argv[5]

chains = []
#
# Get list list of chain names already in use in initial PDB file.
# ( if needed, will generate new names that are not laready in use, must update/check this list.) 

Nters = 0
LastTer = False
chainorder = []
ith_chain = 0
lastchain = '999'
for l in lines:
   l = l.rstrip()
   if l[0:4] == 'ATOM' or l[0:6] == 'HETATM':
      chain = l[21]
      if chain not in chains:
         chains.append(chain)
         chainorder.append([chain])
         ith_chain += 1
      else:
         if LastTer or lastchain != chain:
            chainorder.append([chain])
            ith_chain += 1
      LastTer = False
      lastchain = chain
   else:
      if l[0:3] == 'TER':
         if LastTer: continue  # skip multiple TERs in a row...
         Nters += 1
         LastTer = True

if len(chains) < 1:
   print "ERROR: No Atoms or HetAtms found in PDB file."
   sys.exit(1)


if do_not_rename_these_chains == '':
   if (Nters == 0) or (LastTer and Nters == 1):
      if len(chains) == 1 and chain == ' ':
         print "ERROR: No chain labels and no TER between Atom found in PDB file, cannot separate."
         sys.exit(1)
      print "USING FIRST CHAIN (NOT TERS)."
      UseFirstChain = True
   else:
      UseFirstChain = False
      print "FOUND NON-TERMINAL TER."
else:
   UseFirstChain = False
   

print "Chains labels used in original PDB file:", chains
print "Nchains, Order of chains:", len(chainorder), chainorder
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
chain_has_atoms = False
chain_has_hetatms = False
NRENAMED = 0
Nthchain = 0
justhitTER = False
last_chain_set_in = 1
HitModel = False

for l in lines:

   l = l.rstrip()
   print l

   if l[0:4] == 'ATOM' or l[0:6] == 'HETATM':

      if l[0] == 'H':
         HET = True
         ATM = False
      else:
         HET = False
         ATM = True

      resnum = l[22:27]
      chainraw = l[21]
      chain = chainraw  # Beware! code below uses chain for both read-in and relabled

      if chain != rawchainprevious and rawchainprevious != '00':
         if UseFirstChain:
            if chainprevious not in chain_set1:
               chain_set1.append( chainprevious )
            haveNewChainLabel = False
            hitFirstTER = True
            print >> fhnew, 'TER'
            UseFirstChain = False

      if chain == ' ':
         if not haveNewChainLabel:

            # Create a label for this unlabled chain...
            chain = unique_label( chains )
            NRENAMED += 1
            chains.append( chain ) 
            original[chain] = ' '
            resranges[ chain ] = []
            resnum1 = resnum
            blank_chains_newnames.append( chain )
            chainlast = chain
            haveNewChainLabel = True

         else:
            # Already creates a label for this unlabled chain...
            chain = chainlast
      else:

         haveNewChainLabel = False   # in case we hit another blank chain.

         if hitFirstTER and do_not_rename_these_chains == '': # start second set:

            # OLD: Already have a new label for this labeled chain, re-use it. ( AB_TER_ABABAB OK?)
            # NEW: Revision: to rename chain in case where ABA or AA if reslastA1 > resfirstA2

            # if ATM and last_chain_set_in == 2 and chain == rawchainprevious and resint(resnumlast) <= resint(resnum):
          
            if last_chain_set_in == 2 and chain == rawchainprevious and ((ATM and rescheck(resnumlast,resnum)) or HET):
               # Assuming either a chain break here or the same chain as previous...
               if justhitTER:
                  resnum1 = resnum
                  print " continue second set:", chain, chainprevious, resnum

               chain = chainprevious

            else:
               last_chain_set_in = 2
               if (chain not in chain_set1 and chain not in chain_set2) or chain in do_not_rename_these_chains:
                  # just use original chain...
                  if chain not in original: 
                     original[chain] = chain
               else:
                  # Need to generate new unique label...
                  if chain.isupper():
                     chainuse = chain.lower()
                     if chainuse in chains:
                        chainuse = unique_label_reverse( chains )
                  else:
                     chainuse = unique_label_reverse( chains )

                  # Save a bunch of stuff...
                  chains.append( chainuse )
                  original[chainuse] = chain
                  chain = chainuse
                  NRENAMED += 1

         else: # first set

            # OLD: labeled chains before any TER, leave labels alone?  What about ABAB-TER-CDEF? OK?  
            # NEW: Revision: to rename chain in case where ABA or AA if reslastA1 > resfirstA2

            # Note that this chain in already in chains since we grabbed all
            # chains in first pass thru file above

            if chain == rawchainprevious and ((ATM and rescheck(resnumlast,resnum)) or HET):
               # Assuming either a chain break here or the same chain as previous...

               if justhitTER:
                  resnum1 = resnum
                  print " continue first set:", chain, chainprevious, resnum

               chain = chainprevious

            else:
               if (chain not in chain_set1 and chain not in chain_set2) or chain in do_not_rename_these_chains:
                  # just use original chain...
                  if chain not in original: 
                     original[chain] = chain
               else:
                  # Need to generate new unique label...
                  if chain.isupper():
                     chainuse = chain.lower()
                     if chainuse in chains:
                        chainuse = unique_label_reverse( chains )
                  else:
                     chainuse = unique_label_reverse( chains )

                  # Save a bunch of stuff...
                  chains.append( chainuse )
                  original[chainuse] = chain
                  chain = chainuse
                  NRENAMED += 1

      if HET:
         chain_has_hetatms = True
      if ATM:
         chain_has_atoms = True

      if not justhitTER and chain != chainprevious and chainprevious != '00':
         print "here2:", Nthchain, chainprevious, chain, resnum,  (resnum1,resnumlast), 'raw:',  rawchainprevious, chainraw
         resranges[ chainprevious ].append( (resnum1,resnumlast, chain_has_atoms, chain_has_hetatms ) )
         chainorder[Nthchain].append(chainprevious)
         chainorder[Nthchain].append(last_chain_set_in)
         Nthchain += 1
         resnum1 = resnum

 
      if hitFirstTER and do_not_rename_these_chains == '':   # using second set:
         if chain not in chain_set2:
            chain_set2.append( chain )
            resranges[ chain ] = []
            resnum1 = resnum
            print "here1B:", chain, resnum
      else:
         if chain not in chain_set1:
            chain_set1.append( chain )
            resranges[ chain ] = []
            resnum1 = resnum
            print "here1A:", chain, resnum

      chainprevious = chain
      rawchainprevious = chainraw
      resnumlast = resnum

      print >> fhnew, ''.join( ( l[0:21], chain, l[22:] ) )
      justhitTER = False

   elif l[0:3] == 'TER' and not justhitTER:

      print "here3:", Nthchain, chainprevious, chain, resnum,  (resnum1,resnumlast), 'raw:',  rawchainprevious, chainraw
      resranges[ chainprevious ].append( (resnum1,resnumlast, chain_has_atoms, chain_has_hetatms ) )

      chainorder[Nthchain].append(chainprevious)
      chainorder[Nthchain].append(last_chain_set_in)

      Nthchain += 1
      #resnumlast = '-99999'
      #chainprevious = '00'
      haveNewChainLabel = False
      hitFirstTER = True
      justhitTER = True
      chain_has_atoms = False
      chain_has_hetatms = False

      print >> fhnew, l.strip()

   elif l[0:6] == 'MODEL':      
      print 'M O D E L'
      if HitModel: break
      HitModel = True

print >> fhnew, 'END'
fhnew.close()

if not justhitTER:
   resranges[ chain ].append( (resnum1,resnumlast, chain_has_atoms, chain_has_hetatms) )
   if Nthchain < len(chainorder):
      chainorder[Nthchain].append(chain)
      chainorder[Nthchain].append(last_chain_set_in)
      Nthchain += 1

print "Chainorder:    Nthchain:", Nthchain
for i in range(0, Nthchain):
   print i, "chain order:", chainorder[i]

print "Resranges:"
for chain in resranges:
   print 'chain:', chain
   for rr in resranges[chain]:
      print '    ',rr

print "Chain Set 1 (autoselect):", chain_set1
print "Chain Set 2 (autoselect):", chain_set2

try:
   fhtml = open(jobid + "_chains.html", "w" )
except:
   print "WARNING: Can not open write file in chain_patcher_general.py: reset_chain_label_sets.html"
   sys.exit(0)

html_top( jobid, fhtml )

print >> fhtml, '<H3> Chains Used in Job', jobid, '</H3>'

if do_not_rename_these_chains == '':
   set1 = chain_set1
   set2 = chain_set2
   print >> fhtml, '''
   This page provides you with information about the chains contained in the PDB file for this job.
   Since you did not specify chain labels in the submission form for this job, the KFC2 server attempts to automatically
   define the interface for you.  If a <b>TER</b> record was encountered in the PDB file, then all atoms before the <b>TER</b>
   are used for the first set, and all the atoms after it for the second set.  If there were no <b>TER</b> records found separating atoms,
   then the chain name (column 22 of the <b>ATOM</b> records) is used to define the interface.  The first interface set is taken to be
   all atoms before a second unique chain label is encountered.  All the remaining atoms are then used for the second set.
   HETATM records are listed but are not used in the hotspot predictions.
   Non-water HETATMS will be included when using the viewer to display results.
   </p><p>
   Please check the chain lists below to be sure that the proper chains were automatically selected.  If the chains selected are not correct,
   you will need to edit the PDB file and re-submit a job and consider directly specifying the interface chains in the submission form.
   If redundant chain names were found in this file, they have been renamed as indicated below.
   </p>
   '''

   print >> fhtml, '<pre>'

else:
   set1 = in_chain_set_1
   set2 = in_chain_set_2
   print >> fhtml, '''
   This page provides you with information about the chains contained in the PDB file for this job.
   The interface chain sets used are based on the chain labels you specified in the submission form for this job.
   If redundant chain names were found in this file, they have been renamed as indicated below.
   HETATM records are listed but are not used in the hotspot prediction calculations.
   Non-water HETATMS will be included when using the viewer to display results.
   </p>
   '''

   print >> fhtml, '<pre>'
   print >> fhtml, 'Interface Chain Set 1 from job submission form: ',  set1
   print >> fhtml, 'Interface Chain Set 2 from job submission form: ',  set2
   print >> fhtml, '  '
   print >> fhtml, '  '



setsunion = chain_set1 + chain_set2

chain_set1_not_all_het = []
chain_set2_not_all_het = []

print 'chainorig:', original
print 'chain_set1:',  set1
print 'chain_set2:',  set2

print 'Number Renamed:', NRENAMED 


print >> fhtml, "Interface Chain Set 1:"
print >> fhtml, 106*'-'
for chain in set1:

   chainorig = original[chain]
   colblock = jmolcolspanhtml( chain )

   for respair in resranges[chain]:
      (res1, res2, chain_has_atoms, chain_has_hetatms) = respair

      colblockuse = colblock

      if chain_has_hetatms:
         if chain_has_atoms:
            extra = "ATOMS and HETATMS"
            chain_set1_not_all_het.append(chain)
         else:
            colblockuse = '       '
            extra = "HETATMS only"
      else:
         if chain_has_atoms:
            extra = "ATOMS only"
            chain_set1_not_all_het.append(chain)
         else:
            extra = "NOTHING?"

      res1 = res1.rstrip()
      res2 = res2.lstrip()
      if chain != chainorig:
         print >> fhtml,  colblockuse + '  Original Chain: "' + chainorig + '"   Renamed To: "' + chain + '"   Residue Range:', \
                          res1 + ' - ' + res2, \
                         '  Contains:', extra
      else:
         print >> fhtml,  colblockuse + '  Original Chain: "' + chainorig + '"   Not Renamed.      Residue Range:', \
                          res1 + ' - ' + res2, \
                         '  Contains:', extra

print >> fhtml, ' '
print >> fhtml, ' '
print >> fhtml, "Interface Chain Set 2:"
print >> fhtml, 106*'-'

for chain in set2:


   colblock = jmolcolspanhtml( chain )
   chainorig = original[chain]

   for respair in resranges[chain]:
      (res1, res2, chain_has_atoms, chain_has_hetatms) = respair

      colblockuse = colblock

      if chain_has_hetatms:
         if chain_has_atoms:
            extra = "ATOMS and HETATMS"
            chain_set2_not_all_het.append(chain)
         else:
            colblockuse = '       '
            extra = "HETATMS only"
      else:
         if chain_has_atoms:
            extra = "ATOMS only"
            chain_set2_not_all_het.append(chain)
         else:
            extra = "NOTHING?"

      res1 = res1.rstrip()
      res2 = res2.lstrip()
      if chain != chainorig:
         print >> fhtml, colblockuse + '  Original Chain: "' + chainorig + '"   Renamed To: "' + chain + '"   Residue Range:', \
                          res1 + ' - ' + res2, \
                         '  Contains:', extra
      else:
         print >> fhtml,  colblockuse + '  Original Chain: "' + chainorig + '"   Not Renamed.      Residue Range:', \
                          res1 + ' - ' + res2, \
                         '  Contains:', extra



if len(setsunion) > len(set1) + len(set2):

   print >> fhtml, ' '
   print >> fhtml, ' '
   print >> fhtml, "Not Selected for Interface:"
   print >> fhtml, 106*'-'

   for chain in setsunion:

      if chain not in set1 and chain not in set2:

         colblock = jmolcolspanhtml( chain )

         chainorig = original[chain]

         for respair in resranges[chain]:
            (res1, res2, chain_has_atoms, chain_has_hetatms ) = respair

            colblockuse = colblock

            if chain_has_hetatms:
               if chain_has_atoms:
                  extra = "ATOMS and HETATMS"
               else:
                  extra = "HETATMS only"
                  colblockuse = '       '
            else:
               if chain_has_atoms:
                  extra = "ATOMS only"
               else:
                  extra = "NOTHING?"

            res1 = res1.rstrip()
            res2 = res2.lstrip()
            if chain != chainorig:
               print >> fhtml, colblockuse + '  Original Chain: "' + chainorig + '"   Renamed To: "' + chain + '"   Residue Range:', \
                                res1 + ' - ' + res2, \
                               '  Contains:', extra
            else:
               print >> fhtml,  colblockuse + '  Original Chain: "' + chainorig + '"   Not Renamed.      Residue Range:', \
                                res1 + ' - ' + res2, \
                               '  Contains:', extra


html_bot( fhtml )
fhtml.close()

try:
   fhsh = open("reset_chain_label_sets.sh", "w" )
except:
   print "ERROR: Can not open write file in chain_patcher_general.py: reset_chain_label_sets.sh"
   sys.exit(1)

if do_not_rename_these_chains == '':
   print >> fhsh, 'C1=' + ''.join( chain_set1_not_all_het )
   print >> fhsh, 'C2=' + ''.join( chain_set2_not_all_het )
else:
   print >> fhsh, 'C1=' + in_chain_set_1
   print >> fhsh, 'C2=' + in_chain_set_2

fhsh.close()

sys.exit(0)

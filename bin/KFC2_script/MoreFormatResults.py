#!/usr/bin/env python
#
# gw, 2/2011
#
import sys
import os
from time import strftime

if len(sys.argv) < 3:
   print 'Usage: %s scoresfilename  outputfilename jobname(optional)' % sys.argv[0]
   sys.exit(1)

if not os.path.isfile(sys.argv[1]):
   print "No such file:", sys.argv[1]   
   sys.exit(2)

if len(sys.argv) == 7:
   jobid = sys.argv[3]
   jobname = sys.argv[4]
   chainset1 = sys.argv[5]
   chainset2 = sys.argv[6]

elif len(sys.argv) == 5:
   jobid = sys.argv[3]
   jobname = sys.argv[4]
   chainset1 = '(not given)'
   chainset2 = '(not given)'
else:
   jobid = '(not given)'
   jobname = '(not given)'
   chainset1 = '(not given)'
   chainset2 = '(not given)'

#TESTHEAD=False
#if TESTHEAD:
#   print 'KFC2 Hot Spot Prediction Server @mitchell-lab.org from '+ strftime("%a, %d %b %Y %H:%M:%S %Z\n")
#   print 'JobId: ' + jobid + '   JobName: '+ jobname+'   ChainSet1: '+chainset1+'   ChainSet2: '+chainset2+'\n\n'
#   print '                       KFC2-A  KFC2-A  KFC2-B   KFC2-B \n'
#   print '  Chain    Res    Num  Class    Conf   Class     Conf  \n'
#   print 56*'-' + '\n'
#   sys.exit()

scores = open(sys.argv[1],'r')
lines = scores.readlines()
scores.close()

resname = []
chain = []
resnum = []
score1 = []
score2 = []
chanum = {}

i = 0
for l in lines:
   w = l.strip().split()
   resname.append(w[0])
   chain.append(w[1])
   resnum.append(w[2])
   key = w[1]+w[2]
   chanum[key] = i
   score1.append(w[3])
   score2.append(w[4])
   i += 1

output = open(sys.argv[2],'w')

output.write('KFC2 Hot Spot Prediction Server @mitchell-lab.org from '+ strftime("%a, %d %b %Y %H:%M:%S %Z\n"))
output.write('JobId: ' + jobid + '   JobName: '+ jobname+'   ChainSet1: '+chainset1+'   ChainSet2: '+chainset2+'\n\n')
output.write('                       KFC2-A  KFC2-A  KFC2-B   KFC2-B \n')
output.write('  Chain    Res    Num  Class    Conf   Class     Conf  \n')
output.write(56*'-' + '\n')

keys = chanum.keys()
keys.sort()

for key in keys:
   i = chanum[key]
   s1 = float(score1[i])
   s2 = float(score2[i])
   ## print resname[i], chain[i], resnum[i], score1[i], score2[i]
   if s1 > 0.0:
       class1 = 'Hotspot'
   else:
       class1 = '-------'
   if s2 > 0.0:
       class2 = 'Hotspot'
   else:
       class2 = '-------'

   frmt = '%5s  %7s  %5s  %7s %6.2f  %7s  %6.2f  '

   print >> output, frmt %( chain[i], resname[i], resnum[i], class1, s1, class2, s2 )

output.close()
sys.exit(0)

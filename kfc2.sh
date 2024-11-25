# MIT License
#
# Copyright (c) 2011 Julie C. Mitchell and the University of Wisconsin
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#!/usr/bin/env bash

# Script to run KFC 2 job

###
### Please look for "important installation note" in this script below to make
### sure the KFC2 package installation path is correctly set for your site.
###

###
### TO OBTAIN "naccess", see this Web page:
### http://www.bioinf.manchester.ac.uk/naccess/
###

###
### TO OBTAIN "svm_light", see this Web page:
### http://svmlight.joachims.org/
###

# Initialize variables based on command line input parameters...
PDBIN=$1                           # PDB file name
C1=$2                              # Chains in molecule 1
C2=$3                              # Chains in molecule 2
JOB=$4                             # Job Name

# for some reason, the next command is needed on a mac
export MallocNanoZone=0

echo " "
if [ ! -f $PDBIN ]; then
    echo "Error. PDB file not found:" $PDBIN
    echo "Check command line arg 1 for valid PDB file name."
    exit 1
else
    echo "Using PDB file: $PDBIN"
fi

if [ "$C1" == "" ]; then
   echo "You must specify Chain 1 for interface."
   echo "Use command line arg 2 to specify Chain 1 for interface."
   exit 1
else
   echo "Chain 1 of interface: $C1"
fi

if [ "$C2" == "" ]; then
   echo "You must specify Chain 2 for interface."
   echo "Use command line arg 3 to specify Chain 2 for interface."
   exit 1
else
   echo "Chain 2 of interface: $C2"
fi

if [ "$JOB" == "" ]; then
   echo "You must specify a Job Name for this run."
   echo "Use command line arg 4 to specify a Job Name for this run."
   exit 1
else
   echo "Job Name for this KFC2 run: $JOB"
fi

PARAMSFILE=logfile
if [ -f $PARAMSFILE ]; then
   rm -f $PARAMSFILE
fi

echo " "
echo 'Check file "'${PARAMSFILE}'" for log information regarding this KFC2 run.'
echo " "

echo " " > $PARAMSFILE
echo "Begin KFC2 Job on: "`date` >> $PARAMSFILE
echo " " >> $PARAMSFILE
echo "   Using PDB file: $PDBIN" >> $PARAMSFILE
echo "   Chain 1 of interface: $C1" >> $PARAMSFILE
echo "   Chain 2 of interface: $C2" >> $PARAMSFILE
echo "   Job Name for this KFC2 run: $JOB" >> $PARAMSFILE
echo " " >> $PARAMSFILE
echo "   If you use results from this program, please cite:" >> $PARAMSFILE
echo '        "KFC2: A knowledge-based hot spot prediction method based on' >> $PARAMSFILE
echo '        interface solvation, atomic density and plasticity features",' >> $PARAMSFILE
echo '        X. Zhu and JC Mitchell (2011) Proteins, 79:2671-2683.' >> $PARAMSFILE
echo " " >> $PARAMSFILE
echo "   For KFC2 documentation please see this Web page:" >> $PARAMSFILE
echo '       http://kfc.mitchell-lab.org/kfc2manual.html' >> $PARAMSFILE
echo " " >> $PARAMSFILE


### ---------------important installation note------------------
###
### EDIT NEXT LINE TO POINT TO YOUR LOCAL KFC2 INSTALLATION DIR:
KFC2=/Users/juliemitchell/Desktop/KFC2
KFCBIN=$KFC2/bin                    # KFC2 binaries (splitpdb, calculation, FADE)
KFCDAT=$KFC2/data                   # KFC2 template file

### ---------------important installation note------------------
###
### MAYBE EDIT NEXT LINES TO POINT TO YOUR LOCAL INSTALLATION DIRS FOR NACCS AND SVM_LIGHT:
###
SVMBIN=$KFCBIN/svm_light            # svm_light binary path
export NACCS=$KFCBIN/naccess2.1.1   # naccess binary path (must configure EXE_PATH in naccess)
SCRPT=$KFCBIN/KFC2_script           # KFC2 python scripts (Python, Bash)

function clean_up {
# Clean-up temporary files, set to 1 for clean-up
cleanup_option=1
if [ $cleanup_option == 1 ];
then
   rm -f ${JOB}_$JOB.pdb
   rm -f ${JOB}_$JOB.fad.pdb
   rm -f $JOB.remark.pdb
   rm -f $JOB.mol1.pdb
   rm -f $JOB.mol2.pdb
   rm -f $JOB.mol1mol2.pdb
   rm -f $JOB.fad.pdb
   rm -f $JOB.pdb.nohets
   rm -f $JOB.pars
   rm -f ${JOB}_chains.html
   rm -f ${JOB}_$JOB.*
   rm -f ${JOB}.log
   rm -f input_*
   rm -f reset_chain_label_sets.sh
   rm -f *.asa
   rm -f *.rsa
   rm -f *.fnb
   rm -f *.txt
   rm -f *.result
   rm -f *.temp*
   rm -f temp.out
fi
}

PDB=`basename $PDBIN`

python $SCRPT/scrub_pdb.py $PDBIN $PDB.tmp > temp.out
if [ $? -ne 0 ]; then
    echo " " >> $PARAMSFILE
    echo "Error running scrub_pdb.py." >> $PARAMSFILE
    echo " " >> $PARAMSFILE
    grep -i -e "err" -e "warn" temp.out >> $PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi

grep -i "warn" temp.out >> $PARAMSFILE
rm -f temp.out
mv $PDB.tmp $PDB

python $SCRPT/chain_patcher_general.py $PDB $PDB.tmp $C1 $C2 $JOB > temp.out
if [ $? -ne 0 ]; then
    echo " " >> $PARAMSFILE
    echo "Error running chain_patcher_general.py." >> $PARAMSFILE
    echo " " >> $PARAMSFILE
    grep -i -e "err" -e "warn" temp.out >> $PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi


grep -i "warn" temp.out >> $PARAMSFILE
rm -f temp.out
mv $PDB.tmp $PDB
# REASSIGN CHAIN LABEL SETS TO NEW NAMES...
source ./reset_chain_label_sets.sh


# Remove hydrogens and almost never seen WATs...
egrep -v  -e '^ATOM[ ]+[0-9]+[ ]+H|^ATOM[ ]+[0-9]+[ ]+[1-9]+H| HOH | HOH | WAT ' $PDB > $PDB.temp_nohyds
if [ $? -ne 0 ]; then
    echo " " >> $PARAMSFILE
    echo "Error removing hydrogens and waters with egrep." >> $PARAMSFILE
    echo " " >> $PARAMSFILE
    clean_up
    exit 1
fi

egrep -v '^HETATM ' $PDB.temp_nohyds > $PDB.temp
if [ $? -ne 0 ]; then
    echo " " >> $PARAMSFILE
    echo "Error removing hetatms with egrep." >> $PARAMSFILE
    echo " " >> $PARAMSFILE
    clean_up
    exit 1
fi

mv -f $PDB.temp $PDB

# Don't check for status of the next command since if there are no
# HETATM's found, no hit gives non-zero status. Skip WATERS..
egrep '^HETATM ' $PDB.temp_nohyds | egrep -v ' HOH ' > $PDB.temp_hets
rm -f $PDB.temp_nohyds

if [ $C1 == "_" ] &&  [ $C2 == "_" ]; then
   python $SCRPT/chain_patcher.py $PDB $PDB.tmp > temp.out
   if [ $? -ne 0 ]; then
       echo " " >> $PARAMSFILE
       echo "Error running chain_patcher.py." >> $PARAMSFILE
       echo " " >> $PARAMSFILE
       grep -i -e "err" -e "warn" temp.out >> $PARAMSFILE
       rm -f temp.out
       echo " " >>$PARAMSFILE
       clean_up
       exit 1
   fi
   grep -i "warn" temp.out >> $PARAMSFILE
   rm -f temp.out
   mv $PDB.tmp $PDB 
   # REASSIGN CHAIN LABEL SETS TO NEW NAMES...
   source reset_chain_label_sets.sh
fi


$KFCBIN/KFC2_exec/splitpdb $PDB $JOB.mol1.pdb $C1 $JOB.mol2.pdb $C2 > temp.out
if [ $? -ne 0 ]; then
    echo " " >> $PARAMSFILE
    echo '    Error running "splitpdb" program.  Check PDB file format.' >> $PARAMSFILE
    echo " " >> $PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

# Check to see mol1 file is empty, if so, end here.
if [ ! -s $JOB.mol1.pdb ]; then
    echo " " >>$PARAMSFILE
    echo "    Error. File $JOB.mol1.pdb is empty." >> $PARAMSFILE
    echo "    Please check the chain labels and TER records in your PDB file." >> $PARAMSFILE
    echo " " >> $PARAMSFILE
    clean_up
    exit 1
fi

# Check to see mol2 file is empty, if so, end here.
if [ ! -s $JOB.mol2.pdb ]; then
    echo " " >>$PARAMSFILE
    echo "    Error. File $JOB.mol2.pdb is empty." >>$PARAMSFILE
    echo "    Please check the chain labels and TER records in your PDB file." >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

cp $NACCS/vdw.radii .
cp $NACCS/standard.data .

$NACCS/naccess -r vdw.radii $JOB.mol1.pdb > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error in "naccess" calculation for chain set 1 atoms.' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
#    clean_up
#    exit 1
fi
rm -f temp.out

mv $JOB.rsa  $JOB.mol1.rsa
mv $JOB.asa  $JOB.mol1.asa

$NACCS/naccess -r vdw.radii $JOB.mol2.pdb > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error in "naccess" calculation for chain set 2 atoms.' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    rm -f temp.out
    echo " " >>$PARAMSFILE
#    clean_up
#    exit 1
fi
rm -f temp.out

mv $JOB.rsa  $JOB.mol2.rsa
mv $JOB.asa  $JOB.mol2.asa

cat $JOB.mol1.pdb $JOB.mol2.pdb >$JOB.pdb
$NACCS/naccess -r vdw.radii $JOB.pdb > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error in "naccess" calculation for full complex atoms.' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
#    clean_up
#    exit 1
fi
rm -f temp.out

rm -f vdw.radii
rm -f standard.data
$KFCBIN/KFC2_exec/core_rim_2 $JOB.rsa $JOB.mol1.rsa $JOB.mol2.rsa > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error in "core_rim_2" calculation.' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    rm -f temp.out
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

$KFCBIN/KFC2_exec/core_rim_1 $JOB.rsa $JOB.mol1.rsa $JOB.mol2.rsa > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error in "core_rim_1" calculation.' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

grep -v HETATM $JOB.pdb > $JOB.pdb.nohets
$KFCBIN/KFC2_exec/mut_compensation_1 $JOB.pdb.nohets $JOB.rsa $JOB.mol1.rsa $JOB.mol2.rsa > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error in "mut_compensation_1" calculation.' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f $JOB.pdb.nohets
rm -f temp.out

cat core_inter.txt rim_inter.txt > corerim.txt

$KFCBIN/FADE/FADE -q -p $JOB.mol1.pdb $JOB.mol2.pdb > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error in "FADE" calculation.' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

python $SCRPT/microenv_network/fnearby.py $JOB.mol1.pdb $JOB.mol2.pdb ${JOB}_$JOB.dad 6 10 1 > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error in "microenv fnearby" calculation.' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    clean_up
    exit 1
fi
grep -i "warn"  temp.out  >>$PARAMSFILE
rm -f temp.out


python $SCRPT/FormatData.py core_rim.txt corerim.txt compensation.txt fnearby.fnb input_a.txt input_b.txt ${JOB}_inter.txt > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error running "FormatData.py." program.'>>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

sh $SCRPT/gen_input_kfc2a.sh input_a.txt input_kfc2a.svm > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error running script "gen_input_kfc2a.sh."' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

sh $SCRPT/gen_input_kfc2b.sh input_b.txt input_kfc2b.svm > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error running script "gen_input_kfc2b.sh."' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

$SVMBIN/svm_classify input_kfc2a.svm $KFCDAT/kfc2a.model input_a.result > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error in running "svm_classify" on chain set 1.' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

$SVMBIN/svm_classify input_kfc2b.svm $KFCDAT/kfc2b.model input_b.result > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error in running "svm_classify" on chain set 2.' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

python $SCRPT/arrange_out.py input_a.result input_b.result ${JOB}_inter.txt ${JOB}_kfc2.result > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error running "arrange_out.py."' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

python $SCRPT/MoreFormatResults.py  ${JOB}_kfc2.result $JOB.kfc.results $JOB $JOB $C1 $C2 > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error running "MoreFormatResults.py."' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

# ADD IN AUXILLARY DATA ( Consurf, Rosetta, Experimental ), IF PRESENT...

python $SCRPT/CombineResults.py $JOB > temp.out
if [ $? -ne 0 ]; then
    echo " " >>$PARAMSFILE
    echo '    Error running "CombineResults.py."' >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    grep -i "err"  temp.out  >>$PARAMSFILE
    echo " " >>$PARAMSFILE
    clean_up
    exit 1
fi
rm -f temp.out

mv $JOB.kfc.results.temp $JOB.kfc.results

# Assemble the PDB file which will be fed to the Jmol viewer...

#sed "s/XXXX/$JOB/" $KFCDAT/REMARK.TMP > $JOB.remark.pdb
#grep -v REMARK ${JOB}_$JOB.fad.pdb > $JOB.fad.pdb

if [ -s $PDB.temp_hets ]; then
   cat $JOB.mol1.pdb $JOB.mol2.pdb>  $JOB.mol1mol2.pdb
   python $SCRPT/find_close_hetams.py $JOB.mol1mol2.pdb $PDB.temp_hets  $C1  $C2 > temp.out
   if [ $? -ne 0 ]; then
       echo " " >>$PARAMSFILE
       echo '    Non-fatal error running "find_close_hetams.py"' >>$PARAMSFILE
       echo " " >>$PARAMSFILE
       grep -i "err"  temp.out  >>$PARAMSFILE
       echo " " >>$PARAMSFILE
       cat  $JOB.remark.pdb $JOB.mol1.pdb $JOB.mol2.pdb $PDB.temp_hets.pdb $JOB.fad.pdb > $JOB.fade.pdb
   else
       cat  $JOB.remark.pdb $JOB.mol1mol2.pdb $PDB.temp_hets.pdb $JOB.fad.pdb > $JOB.fade.pdb
   fi 

else
   cat  $JOB.remark.pdb $JOB.mol1.pdb $JOB.mol2.pdb $JOB.fad.pdb > $JOB.fade.pdb
fi

clean_up

echo " " >>$PARAMSFILE
echo "   Job completed successfully." >>$PARAMSFILE
echo " " >>$PARAMSFILE
echo '   See file "'${JOB}'.kfc.results" for KFC2 Hotspot predicitons.' >>$PARAMSFILE
echo " " >>$PARAMSFILE
echo "End KFC2 Job on: "`date` >> $PARAMSFILE
echo " " >>$PARAMSFILE

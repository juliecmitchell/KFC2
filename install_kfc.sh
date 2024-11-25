#! /bin/sh

cd src

rm -f ../bin/KFC2_exec/core_rim_1
rm -f ../bin/KFC2_exec/core_rim_2
rm -f ../bin/KFC2_exec/mut_compensation_1
rm -f ../bin/KFC2_exec/splitpdb


gcc -o ../bin/KFC2_exec/core_rim_1 core_rim_1.c
gcc -o ../bin/KFC2_exec/core_rim_2 core_rim_2.c
gcc -o ../bin/KFC2_exec/mut_compensation_1 mut_compensation_1.c
gcc -o ../bin/KFC2_exec/splitpdb splitpdb.c

echo ""
echo ""
echo "*** DEPENDENCIES ***"
echo ""
echo "Install FADE as ./bin/KFC2_exec/FADE"
echo "https://jcmitchelllab.wordpress.com/2024/07/28/fade-and-padre-molecular-shape-analysis/"
echo ""
echo "Install naccess as ./bin/naccess2.1.1/naccess"
echo "http://www.bioinf.manchester.ac.uk/naccess/nacdownload.html"
echo ""
echo "Install svm_classify into ./bin/svm_light"
echo "https://www.cs.cornell.edu/people/tj/svm_light/"
echo ""
echo ""
echo "*** EXAMPLE ***"
echo ""
echo "sh kfc2.sh data/1dfj.pdb A B 1dfj"

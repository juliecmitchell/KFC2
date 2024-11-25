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

#! /bin/sh

cd src

rm ../bin/KFC2_exec/core_rim_1
rm ../bin/KFC2_exec/core_rim_2
rm ../bin/KFC2_exec/mut_compensation_1
rm ../bin/KFC2_exec/splitpdb


gcc -o ../bin/KFC2_exec/core_rim_1 core_rim_1.c
gcc -o ../bin/KFC2_exec/core_rim_2 core_rim_2.c
gcc -c -ggdb mut_compensation_1.c
gcc -o ../bin/KFC2_exec/mut_compensation_1 mut_compensation_1.o
gcc -o ../bin/KFC2_exec/splitpdb splitpdb.c

echo "Install FADE <------ https://linktofade"
echo "Install naccess2.1.1 <------ https://linktona"
echo "Install svm_light <------ https://svm"

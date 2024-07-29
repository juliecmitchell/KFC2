# KFC2
This is the code needed to run the KFC2 model for predicting alanine scanning mutagenesis hotspots with protein-protein interactions. Please note the need to install NACCESS and svm-light. We are unable to include them with our distribution. 

1) Obtain and install program "naccess" :

      http://www.bioinf.manchester.ac.uk/naccess/

   Make sure to put this inside the "KFC2/bin" directory inside a directory called "naccess2.1.1"

2) Obtain and install program "svm_light" :
 
      http://svmlight.joachims.org/

   Make sure to put this inside the "KFC2/bin" directory, inside a directory called "svm_light" and make sure that this folder contains the "svm_classify" executable.


3) Obtain and install FADE:

      http://mitchell-lab.org

   Again, make sure you add the "naccess2.1.1" and "svm_light" inside the "bin" directory of this distribution. 

5) Edit main script "kfc2.sh" as needed to change any paths (see comments in script):

   Look for comment "--important installation note--" in script.


6) To run a test of the installation:

   If you run test within KFC2 distribution directory, do this:

      mkdir test
      cd test
      ../kfc2.sh ../data/1bsl.pdb  A  B  test
                          ^        ^  ^  ^
                          |        |  |  |
                          |        |  |  |
                          |        |  |  a name for this particular job
                          |        |  |  
                          |        |  PDB chain label for second chain of interface
                          |        |   
                          |        PDB chain label for first chain of interface
                          |
                          full path and name of a PDB file containing at least two chains


   File with hotspot predictions will be called:  "test.kfc.results"
      
   Sample output is "sampleout" directory.

// MIT License
//
// Copyright (c) 2011 Julie C. Mitchell and the University of Wisconsin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define ss 1000

/* GW 5/2011, allow lowercase chainlabels, removed all toupper */

/* global vars */

int main(int argc, char **argv){

	char helpStr[] = "Usage:  splitpdb  infile.pdb  outfile1.pdb  chains1  outfile2.pdb  chains2\n";
	FILE *infile, *outfile1, *outfile2;
	char inName[ss], outName1[ss], outName2[ss], chains1[ss], chains2[ss];
	char inStr[ss], thisChain[1], userChain[1];
	int  i,numChain1, numChain2, goodChain1, goodChain2, TERsplit=0, postTER=0, postENDMDL=0;		

	/* read user input */
	
	if (argc < 6){
		fprintf(stderr,"%s",helpStr);
		exit(1);
	} else {	
		strcpy(inName,argv[1]);
		strcpy(outName1,argv[2]);
		strcpy(chains1,argv[3]);
		strcpy(outName2,argv[4]);
		strcpy(chains2,argv[5]);
	}
	
	/* check if splitting by first TER */
	if ( (strcmp(chains1,".")==0) && (strcmp(chains2,".")==0) )
		TERsplit=1;

        if ( (strcmp(chains1,"_")==0) && (strcmp(chains2,"_")==0) )
                TERsplit=1;
	
	/* open input/output files */
	
	infile = fopen(inName,"r");
	
	if (infile == NULL){
		fprintf(stderr,"ERROR:  could not find file %s\n",inName);
		exit(1);
	}
	
	outfile1 = fopen(outName1,"w");
	outfile2 = fopen(outName2,"w");
	
	/* get number of user-input chains */
	
	numChain1 = strlen(chains1);	
	numChain2 = strlen(chains2);	

	/* read lines of input file */
	
	while ( (fgets(inStr,ss,infile) != NULL) && (!postENDMDL)){
		
		/* check to see if we are past the first model in a multi-model file */
		
		if (strncmp(inStr,"ENDMDL",6) == 0)
			postENDMDL = 1;
					
		/* check to see if we are past the first TER */

//		if (strncmp(inStr,"TER",3) == 0)
		if (inStr[0]=='T' && inStr[1]=='E' && inStr[2]=='R')
			postTER = 1;
		

			/* check to see if line is an atom or TER record */
		
		if  ( (strncmp (inStr,"ATOM",4) == 0) || (strncmp(inStr,"HETATM",6) == 0) )  {
		
			/* check if chain is desirable */
			
			goodChain1 = 0;
			goodChain2 = 0;

			if (TERsplit) {
								
				if (postTER) {
					goodChain2 = 1;
				} else {
					goodChain1 = 1;
				}
				
			} else {
				thisChain[0] = inStr[21];
							
				for (i=0;i<numChain1;i++) {
					userChain[0] = chains1[i];
					if ( strncmp(userChain,thisChain,1) == 0) {
						goodChain1 = 1;
						i = numChain1;
					}
				}

				for (i=0;i<numChain2;i++) {
					userChain[0] = chains2[i];
					if ( strncmp(userChain,thisChain,1) == 0) {
						goodChain2 = 1;
						i = numChain2;
					}
				}
			}
			
			/* write out pdb record */
		
			if (goodChain1) {
				fprintf(outfile1,"%s",inStr);
			}

			if (goodChain2) {
				fprintf(outfile2,"%s",inStr);
			}
			
		}
	}
	
	/* close io files */
	fclose(infile);
	fclose(outfile1);
	fclose(outfile2);
	
return 0;}


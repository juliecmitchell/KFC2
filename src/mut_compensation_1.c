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

#include "mut_com.h"
/***define interface residue distance cutof***/
#define      DCUTOFF 6.0

/***************************************************************************
  This program will calculate the possibility a mutated residue be compensated 
  by residues and waters around it
  *****************************************************************************/

/****************************************************************************
 * Function: main
 ****************************************************************************/
int main(int argc, char **argv) {

    PDBsas_t* com_p;
	PDBsas_t* mono_p1;
	PDBsas_t* mono_p2;
	PDBData_t* mol_com;
	FILE *fp;
	FILE *fp1;
	if (argc<4) {
		fprintf(stderr,"Usage:  mut_compensation_1  complex.pdb complex.rsa monomer1.rsa  monomer2.rsa\n\n");
		exit(1);
	}

	//read in the protein complex
	mol_com = ReadPDB(argv[1]);

	// read in the residues' solvent accessible surface of complex 
	com_p = Read_SAS(argv[2]);

	//read in the residues' solvent accessible surface of monomer
	mono_p1 = Read_SAS(argv[3]);
	mono_p2 = Read_SAS(argv[4]);

	core_rim(com_p, mono_p1, mol_com);
	core_rim(com_p, mono_p2, mol_com);
	if(com_p) 
	{
		free (com_p->resolacs);
		free (com_p);
	}
	if(mono_p1)
	{
		free (mono_p1->resolacs);
		free (mono_p1);
	}
	if(mono_p2)
	{
		free (mono_p2->resolacs);
		free (mono_p2);
	}
	if(mol_com)
	{
		free (mol_com->atoms);
		free (mol_com->residues);
		free (mol_com);
	}


	return 0;
}

/****************************************************************************
 * Function: msf_calc
 ****************************************************************************/
int msf_calc(PDBData_t **mol_p, int nmodel) {
	
	int i,j,k,mod_j;	
	double xyzref[3];
	double total[3];
	double msf;

	// calculate atom-wise msf
	
	for (j=0; j<mol_p[0]->natoms; j++) {
		
		// calculate average coordinate (note, we could use this as a reference point instead of the orig molecule, but it is commented out for now

		mod_j=atomID(mol_p[0],j,mol_p[1]);
		total[0]=0.0;
		total[1]=0.0;
		total[2]=0.0;
		
		for (/*k=0*/k=1; k<nmodel+1; k++) {
			
			total[0] += mol_p[k]->atoms[mod_j].xyz[0];  
			total[1] += mol_p[k]->atoms[mod_j].xyz[1];		        
			total[2] += mol_p[k]->atoms[mod_j].xyz[2];
			
		}
		
		xyzref[0]=(total[0])/(double) nmodel;
		xyzref[1]=(total[1])/(double) nmodel;
		xyzref[2]=(total[2])/(double) nmodel;
				
		// calculate mean square fluctuation (store in mol_p[0])		
		msf = 0.0;
		
		for (k=1; k<nmodel+1; k++) {		
			msf += (mol_p[k]->atoms[mod_j].xyz[0]-xyzref[0])*(mol_p[k]->atoms[mod_j].xyz[0]-xyzref[0]);
			msf += (mol_p[k]->atoms[mod_j].xyz[1]-xyzref[1])*(mol_p[k]->atoms[mod_j].xyz[1]-xyzref[1]);
			msf += (mol_p[k]->atoms[mod_j].xyz[2]-xyzref[2])*(mol_p[k]->atoms[mod_j].xyz[2]-xyzref[2]);
		}
		
		mol_p[0]->atoms[j].msf = msf/(double) nmodel;
		if (mol_p[0]->atoms[j].msf > 999.99) mol_p[0]->atoms[j].msf = 999.99;
				
	}
	
	return 0;
}

/****************************************************************************
 * Function: myfopen
 ****************************************************************************/
FILE* myfopen(char *filename, char *iotype) {
	FILE *fp;

	/* open file for writing */
	if (strncmp(iotype,"w",1)==0) {
		fp = fopen(filename,"w");		
		if (fp == NULL){
			fprintf(stderr,"ERROR:  There is a problem writing the file %s.  Exiting. \n",filename);
			exit(EXIT_FAILURE);
		} 
	} else if (strncmp(iotype,"r",1)==0) {
		fp = fopen(filename,"r");

		if (fp == NULL){
			fprintf(stderr,"ERROR: There is a problem reading the file %s.  Exiting.  \n",filename);
			exit(EXIT_FAILURE);
		} 
	} else {
		fp = NULL;
	}

	return fp;
}

/****************************************************************************
 * Function: ExtractString
 ****************************************************************************/
int ExtractString(int len, char *src, char *dst ) {
    char *ptr;
    char ch;
    int i;

    ptr = dst;
    for( i=0; i<len; i++ ){   
		if( *src ) {   
			ch = *src++;
			*dst++ = ch;
			if(ch != ' ') 
				ptr = dst;
		} else {
			break;
		}
    }
    *ptr = 0;
    
return 0;}

/****************************************************************************
 * Function: deblank
 ****************************************************************************/
int deblank(char *origStr) {
	char tempStr[132];
	
	strcpy(tempStr,origStr);
	sscanf(tempStr,"%s",origStr);
		
return 0;}

/****************************************************************************
 * Function: isDigit
 ****************************************************************************/
int isDigit(char A) {

	if ((A == '0') || (A == '1') || (A == '2') || (A == '3') || (A == '4') ||
		(A == '5') || (A == '6') || (A == '7') || (A == '8') || (A == '9'))
			
	return 1;

return 0;}

/****************************************************************************
 * Function: isH
 ****************************************************************************/
int isH(char *sAtom) {

	if ((sAtom[0] == 'H') || ((sAtom[0] == ' ') && (sAtom[1] == 'H') ) ||
		((isDigit(sAtom[0]) == 1) && (sAtom[1] == 'H') ) )
		
	return 1;

return 0;}

/****************************************************************************
 * Function: isNuc
 ****************************************************************************/
int isNuc(char *sResName) {
	int len, nucleic;
	nucleic = 0; 
	
	deblank(sResName);
	
	len = strlen(sResName);
	if (len <= 2) {
		if ((sResName[len-1] == 'A') || (sResName[len-1] == 'C') || (sResName[len-1] == 'G') ||
			(sResName[len-1] == 'T') || (sResName[len-1] == 'U')) {
			
			nucleic = 1;
			
			if (len == 2) 
				if ((sResName[0] != 'D') && (sResName[0] != 'R'))
					nucleic = 0;
		}
	} 
	    				
	return nucleic;
}

/****************************************************************************
 * Function: fixresname
 ****************************************************************************/
 int fixresname(char* str) {

	if (strcmp(str,"HSD") == 0) strcpy(str,"HIS\0");
	if (strcmp(str,"HSE") == 0) strcpy(str,"HIS\0");
	if (strcmp(str,"HID") == 0) strcpy(str,"HIS\0");
	if (strcmp(str,"HIE") == 0) strcpy(str,"HIS\0");
	if (strcmp(str,"CYX") == 0) strcpy(str,"CYS\0");
	if (strcmp(str,"MSE") == 0) strcpy(str,"MET\0");
	
	return 0;
}

/****************************************************************************
 * Function: ReadPDB
 ****************************************************************************/
PDBData_t* ReadPDB (char *pdbfile) {
	int i,j,k;
	int acount=0;                       /* atom counter */
	int rcount=0;                       /* residue counter */
	int lastresnum = 0;                /* last residue number stored */
	char inStr[200];
	PDBData_t *mol_p=NULL;
	FILE *fp;
    char *blank = " ";
	int atomnum, resnum;
	char s_atomnum[6], s_resnum[5];
	char type[7], atomname[5], alt_conformation, resname[5], tmp_resname[4], chain[2];
	double coordx, coordy, coordz, q, r, sigma_2, epsilon,sufix1,sufix2;
	int s;
	char s_coordx[10], s_coordy[10], s_coordz[10], s_q[10], s_r[10], s_s[10], s_sigma_2[10], s_epsilon[10],s_sufix1[10],
		 s_sufix2[10];
	char   *ptok, *rtok, *ctok;
	atom_t *thisatom;

	/* allocate and initialize new pdb struct  */
	mol_p = (PDBData_t *) malloc(sizeof(PDBData_t));
	memset(mol_p, 0, sizeof(PDBData_t));

	/* Fill the structure with data read in from the .pdb file */
	fp = myfopen(pdbfile, "r");


	while (fgets(inStr, 200, fp)) {
						
		ExtractString(6, &inStr[0], type);
		ExtractString(5, &inStr[6], s_atomnum);
		ExtractString(4, &inStr[12], atomname);
		alt_conformation = inStr[16];
		ExtractString(3, &inStr[17], tmp_resname);
		ExtractString(1, &inStr[21], chain);
		ExtractString(4, &inStr[22], s_resnum);
		ExtractString(8, &inStr[30], s_coordx);
		ExtractString(8, &inStr[38], s_coordy);
		ExtractString(8, &inStr[46], s_coordz);
		ExtractString(6, &inStr[54], s_sufix1);
		ExtractString(6, &inStr[60], s_sufix2);

		atomnum = atoi(s_atomnum);
		resnum = atoi(s_resnum);
		coordx = atof(s_coordx);
		coordy = atof(s_coordy);
		coordz = atof(s_coordz);
		sufix1 = atof(s_sufix1);
		sufix2 = atof(s_sufix2);

		if ( (strncmp(type, "ATOM", 4) == 0) || (strncmp(type, "HETATM", 6) == 0)) {
		    //fprintf(stderr,"Here 0: %s N%cN \n", atomname, alt_conformation );
                    // allow A conformations to be included, reject others...  (gw, 3/2011)
     		    if ( (alt_conformation == ' ') || (alt_conformation == 'A') ) {
		        //printf("Here: %s %s %s %s \n", atomname, tmp_resname, s_resnum, chain );
		        /* to keep consistency with existing data structures, copy into char [5] (gw, 3/2011) */
                        strcpy(resname, blank);
                        strcat(resname, tmp_resname );

			/* Allocate, initialize, then store an atom */
			mol_p->atoms = (atom_t*) realloc (mol_p->atoms, sizeof(atom_t) * (acount+1));

			thisatom = &mol_p->atoms[acount];

			ptok = strtok (atomname, " \n\t");
			strcpy(thisatom->aname, ptok);
			strcat (thisatom->aname, "\0");

			ctok = strtok (chain, " \n\t");
			strcpy(thisatom->chain, chain);
			strcat (thisatom->chain, " \0");

			/* figure out if we need to allocate a new residue_t for this guy */
			if ((rcount == 0) || (resnum != lastresnum)) {

				mol_p->residues = (residue_t*) realloc (mol_p->residues, sizeof(residue_t) * (rcount+1));
				mol_p->residues[rcount].firstatom = acount;

				rtok = strtok (resname, " \n\t");
				strcpy(mol_p->residues[rcount].rname, rtok);
				strcat(mol_p->residues[rcount].rname, "\0");
				fixresname(mol_p->residues[rcount].rname);
				
				mol_p->residues[rcount].rnum = resnum;
				mol_p->residues[rcount].ca = 0;
				//mol_p->residues[rcount].surface = 1;

				mol_p->residues[rcount].phi = 0.0;
				mol_p->residues[rcount].psi = 0.0;

				lastresnum = resnum;

				rcount++;
			} else if (strcmp(mol_p->residues[0].rname,"NTER")==0) {
				rtok = strtok (resname, " \n\t");
				if (strcmp(rtok,"NTER")!=0) {
					strcpy(mol_p->residues[0].rname, rtok);
					strcat(mol_p->residues[0].rname, "\0");
					fixresname(mol_p->residues[0].rname);
				}
			}
                    
			/* get the atom coordinates, charge & radius */
			thisatom->xyz[0] = coordx;
			thisatom->xyz[1] = coordy;
			thisatom->xyz[2] = coordz;
			thisatom->occ = sufix1;
			thisatom->bfact = sufix2;
			thisatom->msf = 0.0;
			thisatom->myres = rcount-1;
			thisatom->resnum = resnum;

			/* store coordinates of the CA atom */

			if ((strncmp(atomname, "CA", 2) == 0) || (strncmp(atomname, " CA", 3) == 0) ) 
				mol_p->residues[rcount-1].ca = acount;
			
			acount++;
                 }
	     } 
	}

	mol_p->natoms = acount;
	mol_p->nres = rcount;

	mol_p->residues = (residue_t*) realloc (mol_p->residues, sizeof(residue_t) * (rcount+1));
	mol_p->residues[rcount].firstatom=acount; /* add extra dummy residue for loop upper bounds */
	
	fclose(fp);
			
	return mol_p;
}


/****************************************************************************
 * Function: WriteModel
 ****************************************************************************/
int WriteModel(PDBData_t* mol_p, char* filename) {

    int i, j, ii, jj, acount;
	FILE* fp;
	int myres;

	fp = myfopen (filename, "w");

	acount = 1;
	
	for (i=0; i< mol_p->natoms; i++) {

		myres = mol_p->atoms[i].myres;
		
		if (!isH(mol_p->atoms[i].aname)){
			if (strlen(mol_p->atoms[i].aname) == 4){
				fprintf(fp,"ATOM %5i   %-4s%3s %c%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n",
					 acount,mol_p->atoms[i].aname,mol_p->residues[myres].rname,mol_p->atoms[i].chain[0],mol_p->residues[myres].rnum,
					 mol_p->atoms[i].xyz[0],mol_p->atoms[i].xyz[1],mol_p->atoms[i].xyz[2],mol_p->atoms[i].msf);
			} else {
				fprintf(fp,"ATOM  %5i  %-4s%3s %c%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n",
					 acount,mol_p->atoms[i].aname,mol_p->residues[myres].rname,mol_p->atoms[i].chain[0],mol_p->residues[myres].rnum,
					 mol_p->atoms[i].xyz[0],mol_p->atoms[i].xyz[1],mol_p->atoms[i].xyz[2],mol_p->atoms[i].msf);
			}
			acount++;
		}
		
	}

        fprintf(fp,"TER    %4i      %3s  %4i\n",
               acount,mol_p->residues[myres].rname,mol_p->residues[myres].rnum);
	
	fclose (fp);

	return 0;
}


/**************************************************************************
 * Function: atomID
 **************************************************************************/
int atomID(PDBData_t *mol_p, int atomid_old,PDBData_t *modller_p1)
{
	int i;
	int atom_id_new=0,res_count_n=0;
	res_count_n=mol_p->atoms[atomid_old].myres;
	for(i=modller_p1->residues[res_count_n].firstatom;i<modller_p1->residues[res_count_n+1].firstatom;i++)
	{
		if(strcmp(modller_p1->atoms[i].aname,mol_p->atoms[atomid_old].aname)==0) break;
	}
	if(i<modller_p1->residues[res_count_n+1].firstatom) atom_id_new=i;
	else
	{
		printf("Error:Could not find atom %s in residue %d in modeller structure!\n",mol_p->atoms[atomid_old].aname,
				modller_p1->residues[res_count_n].rnum);
		exit(1);
	}
	return atom_id_new;
}

/***************************************************************************
 * Function: interRes
 * to find all the interface residues of a complex
 ***************************************************************************/
int interRes(PDBData_t *mol_p, PDBData_t *mol_q)
{
	int i,j,k;
	FILE *fp,*fp1;
	double main_msf,side_msf;
	int main_n,side_n,rotsigbond;
	//find interface residues
	for(i=0;i<mol_p->nres;i++)
	{
		for(j=mol_p->residues[i].firstatom;j<mol_p->residues[i+1].firstatom;j++)
		{
			for(k=0;k<mol_q->natoms;k++)
			{
				if(dist(mol_p->atoms[j].xyz,mol_q->atoms[k].xyz)<DCUTOFF) break;
			}
			if(k<mol_q->natoms)break;
			else
				continue;
		}
		if(j<mol_p->residues[i+1].firstatom)
			mol_p->residues[i].interface=1;
		else
			mol_p->residues[i].interface=0;
	}
	for(i=0;i<mol_q->nres;i++)
	{
		for(j=mol_q->residues[i].firstatom;j<mol_q->residues[i+1].firstatom;j++)
		{
			for(k=0;k<mol_p->natoms;k++)
			{
				if(dist(mol_q->atoms[j].xyz,mol_p->atoms[k].xyz)<DCUTOFF) break;
			}
			if(k<mol_p->natoms)break;
			else
				continue;
		}
		if(j<mol_q->residues[i+1].firstatom)
			mol_q->residues[i].interface=1;
		else
			mol_q->residues[i].interface=0;
	}
	fp=myfopen("interface_6.pdb","w");
	fp1=myfopen("interface_6.rmsf","w");
	//output interface residues' coordinate
	for(i=0;i<mol_p->nres;i++)
		if(mol_p->residues[i].interface==1)
			for(j=mol_p->residues[i].firstatom;j<mol_p->residues[i+1].firstatom;j++)
			{
				if(strlen(mol_p->atoms[j].aname) == 4)
					fprintf(fp,"ATOM  %5i %-4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n",j,mol_p->atoms[j].aname,
							mol_p->residues[i].rname,mol_p->atoms[j].chain[0],mol_p->residues[i].rnum,mol_p->atoms[j].xyz[0],
							mol_p->atoms[j].xyz[1],mol_p->atoms[j].xyz[2],mol_p->atoms[j].occ,mol_p->atoms[j].bfact);
				else
					fprintf(fp,"ATOM  %5i  %-4s%3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n",j,mol_p->atoms[j].aname,
							mol_p->residues[i].rname,mol_p->atoms[j].chain[0],mol_p->residues[i].rnum,mol_p->atoms[j].xyz[0],
							mol_p->atoms[j].xyz[1],mol_p->atoms[j].xyz[2],mol_p->atoms[j].occ,mol_p->atoms[j].bfact);
			}
	fprintf(fp,"TER\n");
	for(i=0;i<mol_q->nres;i++)
		if(mol_q->residues[i].interface==1)
			for(j=mol_q->residues[i].firstatom;j<mol_q->residues[i+1].firstatom;j++)
			{
				if(strlen(mol_q->atoms[j].aname) == 4)
					fprintf(fp,"ATOM  %5i %-4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n",2000+j,mol_q->atoms[j].aname,
							mol_q->residues[i].rname,mol_q->atoms[j].chain[0],mol_q->residues[i].rnum,mol_q->atoms[j].xyz[0],
							mol_q->atoms[j].xyz[1],mol_q->atoms[j].xyz[2],mol_q->atoms[j].occ,mol_q->atoms[j].bfact);
				else
					fprintf(fp,"ATOM  %5i  %-4s%3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n",2000+j,mol_q->atoms[j].aname,
							mol_q->residues[i].rname,mol_q->atoms[j].chain[0],mol_q->residues[i].rnum,mol_q->atoms[j].xyz[0],
							mol_q->atoms[j].xyz[1],mol_q->atoms[j].xyz[2],mol_q->atoms[j].occ,mol_q->atoms[j].bfact);
			}
	fprintf(fp,"END\n");
	fclose(fp);
	//calculate the interface residues' side chain and mainchain+cb average msf
	for(i=0;i<mol_p->nres;i++)
	{
		if(mol_p->residues[i].interface==1)
		{
			main_msf=0.0;
			side_msf=0.0;
			main_n=0;
			side_n=0;
			for(j=mol_p->residues[i].firstatom;j<mol_p->residues[i+1].firstatom;j++)
			{
				if(strcmp(mol_p->atoms[j].aname,"N")==0 || strcmp(mol_p->atoms[j].aname,"CA")==0 || 
						strcmp(mol_p->atoms[j].aname,"C")==0 || strcmp(mol_p->atoms[j].aname,"O")==0 || 
						strcmp(mol_p->atoms[j].aname,"CB")==0 || strcmp(mol_p->atoms[j].aname,"OXT")==0)
				{
					main_n++;
					main_msf+=mol_p->atoms[j].bfact;
				}
				else
				{
					side_n++;
					side_msf+=mol_p->atoms[j].bfact;
				}
			}
			main_msf=main_msf/main_n;
			if(side_n==0) side_msf=0.0;
			else
				side_msf=side_msf/side_n;
			rotsigbond=rot_sig_bond(mol_p->residues[i].rname);
			if(rotsigbond<0)
			{
				printf("ERROR:the residue %4s is not a suitable residue!\n",mol_p->residues[i].rname);
				exit(1);
			}
			fprintf(fp1,"%-7s%4i  %c   %8.3f%8.3f%8.3f%8.3f\n",mol_p->residues[i].rname,mol_p->residues[i].rnum,
					mol_p->atoms[mol_p->residues[i].firstatom].chain[0],main_msf,side_msf,side_msf/rotsigbond,
					side_msf/sqrt(1.0*rotsigbond));
		}
		else
			continue;
	}
	for(i=0;i<mol_q->nres;i++)
	{
		if(mol_q->residues[i].interface==1)
		{
			main_msf=0.0;
			side_msf=0.0;
			main_n=0;
			side_n=0;
			for(j=mol_q->residues[i].firstatom;j<mol_q->residues[i+1].firstatom;j++)
			{
				if(strcmp(mol_q->atoms[j].aname,"N")==0 || strcmp(mol_q->atoms[j].aname,"CA")==0 || 
						strcmp(mol_q->atoms[j].aname,"C")==0 || strcmp(mol_q->atoms[j].aname,"O")==0 || 
						strcmp(mol_q->atoms[j].aname,"CB")==0 || strcmp(mol_q->atoms[j].aname,"OXT")==0)
				{
					main_n++;
					main_msf+=mol_q->atoms[j].bfact;
				}
				else
				{
					side_n++;
					side_msf+=mol_q->atoms[j].bfact;
				}
			}
			main_msf=main_msf/main_n;
			if(side_n==0) side_msf=0.0;
			else
				side_msf=side_msf/side_n;
			if(rotsigbond<0)
			{
				printf("ERROR:the residue %4s is not a suitable residue!\n",mol_q->residues[i].rname);
				exit(1);
			}
			rotsigbond=rot_sig_bond(mol_q->residues[i].rname);
			fprintf(fp1,"%-7s%4i  %c   %8.3f%8.3f%8.3f%8.3f\n",mol_q->residues[i].rname,mol_q->residues[i].rnum,
					mol_q->atoms[mol_q->residues[i].firstatom].chain[0],main_msf,side_msf,side_msf/rotsigbond,
					side_msf/sqrt(1.0*rotsigbond));
		}
		else
			continue;
	}
	fprintf(fp1,"END\n");
	fclose(fp1);
	return 0;
}

double dist(double *coor,double *coor1)
{
	int i;
	double z;
	z=0.0;
	for(i=0;i<3;i++)
		z=z+(coor[i]-coor1[i])*(coor[i]-coor1[i]);
	z=sqrt(z);
	return z;
}

int rot_sig_bond(char *resname)
{
	int rot_bond_num=-100;
	if(strcmp(resname,"GLY")==0 || strcmp(resname,"ALA")==0 || strcmp(resname,"PRO")==0)
		rot_bond_num=0;
	else if(strcmp(resname,"LYS")==0)
		rot_bond_num=4;
	else if(strcmp(resname,"ARG")==0)
		rot_bond_num=5;
	else if(strcmp(resname,"ASN")==0 || strcmp(resname,"ASP")==0 || strcmp(resname,"HIS")==0 || strcmp(resname,"ILE")==0 || 
			strcmp(resname,"LEU")==0 || strcmp(resname,"PHE")==0 || strcmp(resname,"TRP")==0 || strcmp(resname,"TYR")==0)
		rot_bond_num=2;
	else if(strcmp(resname,"CYS")==0 || strcmp(resname,"SER")==0 || strcmp(resname,"THR")==0 ||strcmp(resname,"VAL")==0) 
		rot_bond_num=1;
	else if(strcmp(resname,"GLN")==0 || strcmp(resname,"GLU")==0 || strcmp(resname,"MET")==0)
		rot_bond_num=3;
	else
	{
		printf("ERROR:the residue %4s is not a natural residue\n",resname);
		exit(1);
	}
	return rot_bond_num;
}

double wt_rot_sig_bond(char *resname)
{
	int rot_bond_num=0;
	double wt_rot_num=0.0;
	rot_bond_num=rot_sig_bond(resname);
	if(strcmp(resname,"GLY")==0 || strcmp(resname,"ALA")==0 || strcmp(resname,"PRO")==0)
		wt_rot_num=0.0;
	else if(strcmp(resname,"ASN")==0 || strcmp(resname,"ASP")==0 || strcmp(resname,"ILE")==0 || strcmp(resname,"LEU")==0 || 
			strcmp(resname,"MET")==0 )
		wt_rot_num=rot_bond_num/4.0;
	else if(strcmp(resname,"LYS")==0 || strcmp(resname,"GLN")==0 || strcmp(resname,"GLU")==0 )
		wt_rot_num=rot_bond_num/5.0;
	else if(strcmp(resname,"PHE")==0 || strcmp(resname,"ARG")==0)
		wt_rot_num=rot_bond_num/7.0;
	else if(strcmp(resname,"THR")==0 ||strcmp(resname,"VAL")==0)
		wt_rot_num=rot_bond_num/3.0;
	else if(strcmp(resname,"TRP")==0)
		wt_rot_num=rot_bond_num/10.0;
	else if(strcmp(resname,"HIS")==0)
		wt_rot_num=rot_bond_num/6.0;
	else if(strcmp(resname,"TYR")==0)
		wt_rot_num=rot_bond_num/8.0;
	else if(strcmp(resname,"CYS")==0 || strcmp(resname,"SER")==0)
		wt_rot_num=rot_bond_num/2.0;
	else
	{
		printf("ERROR:the residue %4s is not a natural residue\n",resname);
		exit(1);
	}
	return wt_rot_num;
}
/***************************************************
  **** Function: Read pdb.rsa  *********************
  **************************************************/
PDBsas_t* Read_SAS (char *sasfile)
{
	int rcount=0;
	char inStr[200],head[6];
	PDBsas_t *sas_p=NULL;
	FILE *fp;
	rsa_t *thisres;
	/* allocate and initialize new pdb struct  */
	sas_p = (PDBsas_t *) malloc(sizeof(PDBsas_t));
	memset(sas_p, 0, sizeof(PDBsas_t));
	/* Fill the structure with data read in from the pdb.rsa file */
	fp = myfopen(sasfile, "r");
	while(!feof(fp))
	{
		fgets(inStr,200,fp);
		sscanf(inStr,"%s",head);
		if(strcmp(head,"RES")==0)
		{
			sas_p->resolacs = (rsa_t*) realloc (sas_p->resolacs, sizeof(rsa_t) * (rcount+1));
			thisres = &sas_p->resolacs[rcount];
			if(inStr[8]!=' ')
			{
				thisres->chain[0]=inStr[8];
				inStr[8]=' ';
				sscanf(inStr, "%*s%s%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", thisres->rname, &thisres->rnum,
						&thisres->abs_tot, &thisres->rel_tot, &thisres->abs_sc, &thisres->rel_sc, &thisres->abs_mc, &thisres->rel_mc,
						&thisres->abs_npol, &thisres->rel_npol, &thisres->abs_pol, &thisres->rel_pol);
				rcount++;
			}
			else
			{
				thisres->chain[0] = '_';
				sscanf(inStr, "%*s%s%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", thisres->rname, &thisres->rnum,&thisres->abs_tot, 
						&thisres->rel_tot, &thisres->abs_sc, &thisres->rel_sc, &thisres->abs_mc, &thisres->rel_mc,
						&thisres->abs_npol, &thisres->rel_npol, &thisres->abs_pol, &thisres->rel_pol);
				rcount++;
			}
		}
	}
	fclose(fp);
	sas_p->nres = rcount;
	return sas_p;
}

int  core_rim(PDBsas_t *com_p, PDBsas_t *mono_p, PDBData_t *mol_com)
{
	int i,j,count;
	double deltarsa,per,sc_delta,maxASA,maxASAsc,ratio1,ratio2,ratio3,ratio4, ratio5,delta_np,delta_p,np_per,p_per;
	double score[14], ResMWsc;
	FILE *fp;
	FILE *fp1;
	//res_feature *res_f;
	around_res_t around_f, around_f_1, around_f_2;
	if((fp=fopen("compensation.txt","a"))==NULL) 
	{
		printf("could not open file :%s\n","compensation.txt");
		exit(1);
	}
	/*if((fp1=fopen("rim_inter.txt","a"))==NULL)
	{
		printf("could not open file :%s\n","rim_inter.txt");
		exit(1);
	}*/

	for(i=0;i<com_p->nres;i++)
	{
		for(j=0;j<mono_p->nres;j++)
		{
			if(strcmp(com_p->resolacs[i].rname,mono_p->resolacs[j].rname)==0 && 
					com_p->resolacs[i].chain[0] == mono_p->resolacs[j].chain[0] && 
					com_p->resolacs[i].rnum == mono_p->resolacs[j].rnum)
			{
				if(mono_p->resolacs[j].abs_tot > EPSMIN )
				{
					if(strcmp(com_p->resolacs[i].rname,"GLY")==0)
					{
						/*maxASA = MaxASA(com_p->resolacs[i].rname);
						deltarsa = fabs(com_p->resolacs[i].abs_tot - mono_p->resolacs[j].abs_tot);
						per = com_p->resolacs[i].abs_tot/mono_p->resolacs[j].abs_tot;
						sc_delta = 0.00;
						ratio1 = 100*com_p->resolacs[i].abs_tot/maxASA;
						ratio2 = 100*deltarsa/maxASA;
						ratio3 = 100.00;
						ratio4 = 0.00;
						ratio5 = deltarsa*maxASA/mono_p->resolacs[j].abs_tot;
						delta_np = fabs(com_p->resolacs[i].abs_npol - mono_p->resolacs[j].abs_npol);
						if(mono_p->resolacs[j].abs_npol < EPSMIN)
							np_per = 1.0;
						else
							np_per = com_p->resolacs[i].abs_npol/mono_p->resolacs[j].abs_npol;
						delta_p = fabs(com_p->resolacs[i].abs_pol - mono_p->resolacs[j].abs_pol);
						if(mono_p->resolacs[j].abs_pol < EPSMIN)
							p_per = 1.0;
						else
							p_per = com_p->resolacs[i].abs_pol/mono_p->resolacs[j].abs_pol;
						if(deltarsa > EPSMIN && (per - 0.3)<= EPSMIN)
						{
							fprintf(fp,"%-4s%c%5d%9.2f%8.2f%8.2f%8.2f%9.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%9.2f%9.2f%9.2f%6.2f%9.2f%9.2f%9.2f%6.2f\n", com_p->resolacs[i].rname, com_p->resolacs[i].chain[0],com_p->resolacs[i].rnum, 
									com_p->resolacs[i].abs_tot, mono_p->resolacs[j].abs_tot, deltarsa, per, com_p->resolacs[i].abs_sc, 
									mono_p->resolacs[j].abs_sc, sc_delta, ratio1, ratio2, ratio3, ratio4, ratio5, com_p->resolacs[i].abs_npol,
									mono_p->resolacs[j].abs_npol,delta_np, np_per, com_p->resolacs[i].abs_pol, mono_p->resolacs[j].abs_pol,
									delta_p, p_per);
						}
						else if(deltarsa > EPSMIN && (per - 0.3) > EPSMIN)
						{
							fprintf(fp1,"%-4s%c%5d%9.2f%8.2f%8.2f%8.2f%9.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%9.2f%9.2f%9.2f%6.2f%9.2f%9.2f%9.2f%6.2f\n", com_p->resolacs[i].rname, com_p->resolacs[i].chain[0], com_p->resolacs[i].rnum, 
									com_p->resolacs[i].abs_tot, mono_p->resolacs[j].abs_tot, deltarsa, per, com_p->resolacs[i].abs_sc, 
									mono_p->resolacs[j].abs_sc, sc_delta, ratio1, ratio2, ratio3, ratio4, ratio5, com_p->resolacs[i].abs_npol,
									mono_p->resolacs[j].abs_npol,delta_np, np_per, com_p->resolacs[i].abs_pol, mono_p->resolacs[j].abs_pol,
									delta_p, p_per);
						}
						else
							continue;*/
						deltarsa = fabs(com_p->resolacs[i].abs_tot - mono_p->resolacs[j].abs_tot);
						if(deltarsa > EPSMIN)
						fprintf(fp,"%-4s%c%6d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %9.2f %9.2f %7.2f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %8.3f\n",
								com_p->resolacs[i].rname, com_p->resolacs[i].chain[0], com_p->resolacs[i].rnum, 0.0, 0.0, 0.0, 0.0, 0.0, 
								0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 100.0, 100.0, 100.0, 
								100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0);
						else
							continue;
					}
					else
					{
						sc_delta = fabs(com_p->resolacs[i].abs_sc - mono_p->resolacs[j].abs_sc);
						if(sc_delta > EPSMIN)
						{
							feature_around(com_p->resolacs[i], mol_com, &around_f, 4.0);
							feature_around(com_p->resolacs[i], mol_com, &around_f_1, 5.0);
							feature_around(com_p->resolacs[i], mol_com, &around_f_2, 5.5);
							maxASAsc = MaxASAsc(com_p->resolacs[i].rname);
							ResMWsc = resMW_sc(com_p->resolacs[i].rname);
							ratio3 = 100*com_p->resolacs[i].abs_sc/maxASAsc;
							score[0] = around_f.flexibility * around_f.res_tot * around_f.res_tot_wt * ratio3 / maxASAsc;
							score[1] = around_f_1.flexibility * around_f_1.res_tot * around_f_1.res_tot_wt * ratio3 / maxASAsc;
							score[2] = around_f_2.flexibility * around_f_2.res_tot * around_f_2.res_tot_wt * ratio3 / maxASAsc;
							score[3] = around_f.flexibility * around_f.res_tot * around_f.res_tot_wt * ratio3 / ResMWsc;
							score[4] = around_f_1.flexibility * around_f_1.res_tot * around_f_1.res_tot_wt * ratio3 / ResMWsc;
							score[5] = around_f_2.flexibility * around_f_2.res_tot * around_f_2.res_tot_wt * ratio3 / ResMWsc;
							score[6] = around_f.flexibility * around_f.res_tot * around_f.res_tot_wt / maxASAsc;
							score[7] = around_f_1.flexibility * around_f_1.res_tot * around_f_1.res_tot_wt / maxASAsc;
							score[8] = around_f_2.flexibility * around_f_2.res_tot * around_f_2.res_tot_wt / maxASAsc;
							score[9] = around_f.flexibility * around_f.res_tot * around_f.res_tot_wt / ResMWsc;
							score[10] = around_f_1.flexibility * around_f_1.res_tot * around_f_1.res_tot_wt / ResMWsc;
							score[11] = around_f_2.flexibility * around_f_2.res_tot * around_f_2.res_tot_wt / ResMWsc;
							score[12] = 100.0*around_f_1.flexibility * around_f_1.res_tot_wt * ratio3 / pow(maxASAsc,1.5);
							score[13] = 100.0*around_f_1.flexibility * around_f_1.res_tot_wt / pow(maxASAsc,1.5);
							fprintf(fp,"%-4s%c%6d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %9.2f %9.2f %7.2f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %8.3f\n",
								com_p->resolacs[i].rname, com_p->resolacs[i].chain[0], com_p->resolacs[i].rnum, around_f.flexibility, 
								around_f.nflex, around_f.res_tot, around_f.atm_tot, around_f.res_tot_wt, around_f.hydrophobicity, 
								around_f_1.flexibility, around_f_1.nflex, around_f_1.res_tot, around_f_1.atm_tot, around_f_1.res_tot_wt, 
								around_f_1.hydrophobicity, around_f_2.flexibility, around_f_2.nflex, 
								around_f_2.res_tot, around_f_2.atm_tot, around_f_2.res_tot_wt, around_f_2.hydrophobicity, maxASAsc, 
								ResMWsc, ratio3, score[0], score[1], score[2], score[3], score[4], score[5], score[6], score[7], score[8], 
								score[9], score[10], score[11], score[12], score[13]);
						}
						else
							continue;
					}
				}
				break;
			}
		}
	}
	fclose(fp);
	return 0;
}

/*********************************************************************
Function: MaxASA
used to calculate the maxASA of the whole residue
the values come from the article "Quantifying the accessible surface area of protein residues in their local environment"
*********************************************************************/
double MaxASA(char *resname)
{
	double res_asa=0.0;
	if(strcmp(resname, "GLY")==0)
		res_asa = 70.27;
	else if(strcmp(resname, "ALA")==0)
		res_asa = 102.68;
	else if(strcmp(resname, "SER")==0)
		res_asa = 111.97;
	else if(strcmp(resname, "CYS")==0)
		res_asa = 127.72;
	else if(strcmp(resname, "THR")==0)
		res_asa = 134.28;
	else if(strcmp(resname, "ASP")==0)
		res_asa = 141.61;
	else if(strcmp(resname, "PRO")==0)
		res_asa = 126.78;
	else if(strcmp(resname, "ASN")==0)
		res_asa = 155.22;
	else if(strcmp(resname, "VAL")==0)
		res_asa = 148.51;
	else if(strcmp(resname, "GLU")==0)
		res_asa = 173.46;
	else if(strcmp(resname, "GLN")==0)
		res_asa = 175.42;
	else if(strcmp(resname, "HIS")==0)
		res_asa = 184.79;
	else if(strcmp(resname, "LEU")==0)
		res_asa = 184.33;
	else if(strcmp(resname, "ILE")==0)
		res_asa = 176.22;
	else if(strcmp(resname, "MET")==0)
		res_asa = 196.87;
	else if(strcmp(resname, "LYS")==0)
		res_asa = 193.73;
	else if(strcmp(resname, "PHE")==0)
		res_asa = 209.64;
	else if(strcmp(resname, "TYR")==0)
		res_asa = 224.68;
	else if(strcmp(resname, "ARG")==0)
		res_asa = 235.45;
	else if(strcmp(resname, "TRP")==0)
		res_asa = 251.78;
	else
	{
		printf("ERROR: this is not a natural residue!\n");
		exit(1);
	}
	return res_asa;
}

/***********************************************************
  *** Function: MaxASAsc
  ***** similar with MaxASA but used for side chain********/
double MaxASAsc(char *resname)
{
	double sc_asa=0.0;
    if(strcmp(resname, "GLY")==0)
		sc_asa = 0.00;
	else if(strcmp(resname, "ALA")==0)
		sc_asa = 55.28;
	else if(strcmp(resname, "SER")==0)
		sc_asa = 68.97;
	else if(strcmp(resname, "CYS")==0)
		sc_asa = 81.91;
	else if(strcmp(resname, "THR")==0)
		sc_asa = 88.45;
	else if(strcmp(resname, "ASP")==0)
		sc_asa = 97.66;
	else if(strcmp(resname, "PRO")==0)
		sc_asa = 98.23;
	else if(strcmp(resname, "ASN")==0)
		sc_asa = 109.87;
	else if(strcmp(resname, "VAL")==0)
		sc_asa = 103.00;
	else if(strcmp(resname, "GLU")==0)
		sc_asa = 132.42;
	else if(strcmp(resname, "GLN")==0)
		sc_asa = 129.52;
	else if(strcmp(resname, "HIS")==0)
		sc_asa = 141.17;
	else if(strcmp(resname, "LEU")==0)
		sc_asa = 141.47;
	else if(strcmp(resname, "ILE")==0)
		sc_asa = 130.58;
	else if(strcmp(resname, "MET")==0)
		sc_asa = 150.32;
	else if(strcmp(resname, "LYS")==0)
		sc_asa = 147.83;
	else if(strcmp(resname, "PHE")==0)
		sc_asa = 164.14;
	else if(strcmp(resname, "TYR")==0)
		sc_asa = 180.01;
	else if(strcmp(resname, "ARG")==0)
		sc_asa = 190.04;
	else if(strcmp(resname, "TRP")==0)
		sc_asa = 209.57;
	else
	{
		printf("ERROR: this is not a natural residue!\n");
		exit(1);
	}
	return sc_asa;
}

/***********************************************************
 *** Function: resMW
 ***** residues' molecular weight********/
double resMW(char *resname)
{
	double mol_w=0.0;
	if(strcmp(resname, "GLY")==0)
		mol_w = 75.0669 - 18;
	else if(strcmp(resname, "ALA")==0)
		mol_w = 89.0935 - 18;
	else if(strcmp(resname, "SER")==0)
		mol_w = 105.0930 - 18;
	else if(strcmp(resname, "CYS")==0)
		mol_w = 121.1590 - 18;
	else if(strcmp(resname, "THR")==0)
		mol_w = 119.1197 - 18;
	else if(strcmp(resname, "ASP")==0)
		mol_w = 133.1032 - 18;
	else if(strcmp(resname, "PRO")==0)
		mol_w = 115.1310 - 18;
	else if(strcmp(resname, "ASN")==0)
		mol_w = 132.1184 - 18;
	else if(strcmp(resname, "VAL")==0)
		mol_w = 117.1469 - 18;
	else if(strcmp(resname, "GLU")==0)
		mol_w = 147.1299 - 18;
	else if(strcmp(resname, "GLN")==0)
		mol_w = 146.1451 - 18;
	else if(strcmp(resname, "HIS")==0)
		mol_w = 155.1552 - 18;
	else if(strcmp(resname, "LEU")==0)
		mol_w = 131.1736 - 18;
	else if(strcmp(resname, "ILE")==0)
		mol_w = 131.1736 - 18;
	else if(strcmp(resname, "MET")==0)
		mol_w = 149.2124 - 18;
	else if(strcmp(resname, "LYS")==0)
		mol_w = 146.1882 - 18;
	else if(strcmp(resname, "PHE")==0)
		mol_w = 165.1900 - 18;
	else if(strcmp(resname, "TYR")==0)
		mol_w = 181.1894 - 18;
	else if(strcmp(resname, "ARG")==0)
		mol_w = 174.2017 - 18;
	else if(strcmp(resname, "TRP")==0)
		mol_w = 204.2262 - 18;
	else
	{
		printf("ERROR: this is not a natural residue!\n");
		exit(1);
	}
	return mol_w;
}

/***********************************************************
 *** Function: resMW_sc
 ***** residues' side chain molecular weight********/
double resMW_sc(char *resname)
{
	double mol_w_sc=0.0;
	if(strcmp(resname, "GLY")==0)
		mol_w_sc = 75.0669 - 18 -57.0669;
	else if(strcmp(resname, "ALA")==0)
		mol_w_sc = 89.0935 - 18 - 56.0669;
	else if(strcmp(resname, "SER")==0)
		mol_w_sc = 105.0930 - 18 - 56.0669;
	else if(strcmp(resname, "CYS")==0)
		mol_w_sc = 121.1590 - 18 - 56.0669;
	else if(strcmp(resname, "THR")==0)
		mol_w_sc = 119.1197 - 18 - 56.0669;
	else if(strcmp(resname, "ASP")==0)
		mol_w_sc = 133.1032 - 18 - 56.0669;
	else if(strcmp(resname, "PRO")==0)
		mol_w_sc = 115.1310 - 18 - 56.0669;
	else if(strcmp(resname, "ASN")==0)
		mol_w_sc = 132.1184 - 18 - 56.0669;
	else if(strcmp(resname, "VAL")==0)
		mol_w_sc = 117.1469 - 18 - 56.0669;
	else if(strcmp(resname, "GLU")==0)
		mol_w_sc = 147.1299 - 18 - 56.0669;
	else if(strcmp(resname, "GLN")==0)
		mol_w_sc = 146.1451 - 18 - 56.0669;
	else if(strcmp(resname, "HIS")==0)
		mol_w_sc = 155.1552 - 18 - 56.0669;
	else if(strcmp(resname, "LEU")==0)
		mol_w_sc = 131.1736 - 18 - 56.0669;
	else if(strcmp(resname, "ILE")==0)
		mol_w_sc = 131.1736 - 18 - 56.0669;
	else if(strcmp(resname, "MET")==0)
		mol_w_sc = 149.2124 - 18 - 56.0669;
	else if(strcmp(resname, "LYS")==0)
		mol_w_sc = 146.1882 - 18 - 56.0669;
	else if(strcmp(resname, "PHE")==0)
		mol_w_sc = 165.1900 - 18 - 56.0669;
	else if(strcmp(resname, "TYR")==0)
		mol_w_sc = 181.1894 - 18 - 56.0669;
	else if(strcmp(resname, "ARG")==0)
		mol_w_sc = 174.2017 - 18 - 56.0669;
	else if(strcmp(resname, "TRP")==0)
		mol_w_sc = 204.2262 - 18 - 56.0669;
	else
	{
		printf("ERROR: this is not a natural residue!\n");
		exit(1);
	}
	return mol_w_sc;
}

int feature_around(rsa_t r_rsa, PDBData_t * mol_com, around_res_t *around_f, double ndist)
{
	int i,j,k;
	int num1, res_n, num2;
	char resnames[30][6];
	for(i=0;i<mol_com->nres;i++)
	{
		if(strcmp(r_rsa.rname,mol_com->residues[i].rname)==0 && r_rsa.rnum == mol_com->residues[i].rnum && r_rsa.chain[0] == mol_com->atoms[mol_com->residues[i].firstatom].chain[0])
		{
			res_n=i;
			break;
		}
	}
	num2 = 0;
	for(i=0; i<mol_com->natoms; i++)
	{
		if((i>=mol_com->residues[res_n].firstatom && i<mol_com->residues[res_n+1].firstatom) || strcmp(mol_com->atoms[i].aname,"N")==0 || 
				strcmp(mol_com->atoms[i].aname,"CA")==0 || strcmp(mol_com->atoms[i].aname,"C")==0 || strcmp(mol_com->atoms[i].aname,"O")==0 || 
				strcmp(mol_com->atoms[i].aname,"OXT")==0)
			continue;
		for(j=mol_com->residues[res_n].firstatom; j<mol_com->residues[res_n+1].firstatom; j++)
		{
			if(strcmp(mol_com->atoms[j].aname,"N")==0 || strcmp(mol_com->atoms[j].aname,"CA")==0 || 
					strcmp(mol_com->atoms[j].aname,"C")==0 || strcmp(mol_com->atoms[j].aname,"O")==0 || 
					strcmp(mol_com->atoms[j].aname,"OXT")==0)
				continue;
			else if(dist(mol_com->atoms[i].xyz, mol_com->atoms[j].xyz)<=ndist)  break;
			else
				continue;
		}
		if(j<mol_com->residues[res_n+1].firstatom)
			num2++;
		else
			continue;
	}
	around_f->atm_tot = 1.0 * num2;
	num1=0;
	for(i=0;i<mol_com->nres;i++)
	{
		if(i==res_n) continue;
		for(j=mol_com->residues[i].firstatom; j<mol_com->residues[i+1].firstatom; j++)
		{
			if(strcmp(mol_com->atoms[j].aname,"N")==0 || strcmp(mol_com->atoms[j].aname,"CA")==0 || 
					strcmp(mol_com->atoms[j].aname,"C")==0 || strcmp(mol_com->atoms[j].aname,"O")==0 ||
					strcmp(mol_com->atoms[j].aname,"OXT")==0)
				continue;
			for(k=mol_com->residues[res_n].firstatom;k<mol_com->residues[res_n+1].firstatom;k++)
			{
				if(strcmp(mol_com->atoms[k].aname,"N")==0 || strcmp(mol_com->atoms[k].aname,"CA")==0 || 
						strcmp(mol_com->atoms[k].aname,"C")==0 || strcmp(mol_com->atoms[k].aname,"O")==0 || 
						strcmp(mol_com->atoms[k].aname,"OXT")==0)
					continue;
				if(dist(mol_com->atoms[j].xyz,mol_com->atoms[k].xyz)<=ndist) break;
			}
			if(k<mol_com->residues[res_n+1].firstatom) break;
			else
				continue;
		}
		if(j<mol_com->residues[i+1].firstatom)
		{
			strcpy(resnames[num1],mol_com->residues[i].rname);
			num1++;
			if(num1==30)
			{
				printf("Error: the residues around residue %-5s%4d%3s is more than 30\n",r_rsa.rname, r_rsa.rnum, r_rsa.chain);
				exit(1);
			}
		}
		else
			continue;
	}
	around_f->res_tot = 1.0*num1;
	around_f->flexibility = 0.0;
	around_f->nflex = 0.0;
	around_f->res_tot_wt = 0.0;
	around_f->hydrophobicity = 0.0;
	for(i=0; i<num1; i++)
	{
		around_f->flexibility = around_f->flexibility + rot_sig_bond(resnames[i]); //number of rotatable single bonds could be substitued by number of rotamers
		around_f->nflex = around_f->nflex + wt_rot_sig_bond(resnames[i]);
		around_f->res_tot_wt = around_f->res_tot_wt + simi_character(r_rsa.rname, resnames[i]);
		around_f->hydrophobicity = around_f->hydrophobicity + resHP(resnames[i]);
	}
	return 0;
}

// similarity between residues
double simi_character(char *res1name, char *res2name)
{
	int i, j;
	double simi_score=0.0;
	static double simi_matrix[11][11] = {{0.8,0.4,0.0,0.0,0.1,0.0,0.0,0.0,0.0,0.1,0.0},
										 {0.5,0.8,0.2,0.1,0.3,0.0,0.0,0.0,0.1,0.1,0.7},
										 {0.2,0.6,0.4,0.1,0.1,0.0,0.0,0.0,0.1,0.1,0.4},
										 {0.1,0.5,0.2,0.6,0.1,0.0,0.0,0.0,0.4,0.1,0.1},
										 {0.5,0.4,0.1,0.1,0.6,0.0,0.0,0.0,0.2,0.2,0.0},
										 {0.8,0.7,0.4,0.3,0.4,0.0,0.0,0.0,0.2,0.3,0.4},
										 {0.8,0.8,0.4,0.4,0.6,0.0,0.0,0.0,0.2,0.6,0.4},
										 {0.8,0.6,0.1,0.1,0.4,0.0,0.0,0.0,0.8,0.6,0.2},
										 {0.4,0.2,0.1,0.1,0.2,0.0,0.0,0.1,0.8,0.6,0.0},
										 {0.4,0.1,0.0,0.0,0.2,0.0,0.0,0.0,0.3,0.8,0.0},
										 {0.0,0.7,0.2,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.8}};
	i = classify(res1name);
	j = classify(res2name);
	simi_score = simi_matrix[i][j];
	return simi_score;
}

//classify residues according to their character
int classify(char *resname)
{
	int type=0;
	if(strcmp(resname, "ARG")==0 || strcmp(resname,"LYS")==0)
		type=0;
	else if(strcmp(resname, "GLN")==0 || strcmp(resname,"ASN")==0)
		type=1;
	else if(strcmp(resname, "SER")==0 || strcmp(resname,"THR")==0)
		type=2;
	else if(strcmp(resname, "CYS")==0)
		type=3;
	else if(strcmp(resname, "HIS")==0)
		type=4;
	else if(strcmp(resname, "GLY")==0)
		type=5;
	else if(strcmp(resname, "PRO")==0)
		type=6;
	else if(strcmp(resname, "ALA")==0 || strcmp(resname, "VAL")==0)
		type=7;
	else if(strcmp(resname, "ILE")==0 || strcmp(resname, "LEU")==0 || strcmp(resname, "MET")==0)
		type=8;
	else if(strcmp(resname, "TRP")==0 || strcmp(resname, "PHE")==0 || strcmp(resname, "TYR")==0)
		type=9;
	else 
		type=10;
	return type;
}

//residues' hydrophobicity
double resHP(char *resname)
{
	double hp=0.0;
	if(strcmp(resname, "ALA") == 0) hp = 0.31;
	if(strcmp(resname, "LEU") == 0) hp = 1.70;
	if(strcmp(resname, "ARG") == 0) hp = -1.01;
	if(strcmp(resname, "LYS") == 0) hp = -0.99;
	if(strcmp(resname, "ASN") == 0) hp = -0.60;
	if(strcmp(resname, "MET") == 0) hp = 1.23;
	if(strcmp(resname, "ASP") == 0) hp = -0.77;
	if(strcmp(resname, "PHE") == 0) hp = 1.79;
	if(strcmp(resname, "CYS") == 0) hp = 1.54;
	if(strcmp(resname, "PRO") == 0) hp = 0.72;
	if(strcmp(resname, "GLN") == 0) hp = -0.22;
	if(strcmp(resname, "SER") == 0) hp = -0.04;
	if(strcmp(resname, "GLU") == 0) hp = -0.64;
	if(strcmp(resname, "THR") == 0) hp = 0.26;
	if(strcmp(resname, "GLY") == 0) hp = 0.00;
	if(strcmp(resname, "TRP") == 0) hp = 2.25;
	if(strcmp(resname, "HIS") == 0) hp = 0.13;
	if(strcmp(resname, "TYR") == 0) hp = 0.96;
	if(strcmp(resname, "ILE") == 0) hp = 1.80;
	if(strcmp(resname, "VAL") == 0) hp = 1.22;
	return hp;
}

#include "core_rim_1.h"
/***define interface residue distance cutof***/
#define      DCUTOFF 6.0
/***************************************************************************
  Compared with core_rim.c, I will calculated the absolutely solvent accessible
  surface area for residue's sidechain and delta of that residue and these value
  compared with the maxASAi
  *****************************************************************************/

/****************************************************************************
 * Function: main
 ****************************************************************************/
//PDBData_t* mol_p[20];
int main(int argc, char **argv) {
	
	PDBsas_t* com_p;
	PDBsas_t* mono_p1;
	PDBsas_t* mono_p2;
	FILE *fp;
	FILE *fp1;
	if (argc<4) {
		fprintf(stderr,"Usage:  core_rim  complex.rsa monomer1.rsa  monomer2.rsa\n\n");
		exit(1);
	}
	
	// read in the residues' solvent accessible surface of complex 
	com_p = Read_SAS(argv[1]);

	//read in the residues' solvent accessible surface of monomer
	mono_p1 = Read_SAS(argv[2]);
	mono_p2 = Read_SAS(argv[3]);

	core_rim(com_p,mono_p1);
	core_rim(com_p,mono_p2);
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


	return 0;
}

/****************************************************************************
 * Function: msf_calc
 ****************************************************************************/
double msf_calc(PDBData_t **mol_p, int nmodel) {
	
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
void ExtractString(int len, char *src, char *dst ) {
    register char *ptr;
    register char ch;
    register int i;

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
    
}

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
 void fixresname(char* str) {

	if (strcmp(str,"HSD") == 0) strcpy(str,"HIS\0");
	if (strcmp(str,"HSE") == 0) strcpy(str,"HIS\0");
	if (strcmp(str,"HID") == 0) strcpy(str,"HIS\0");
	if (strcmp(str,"HIE") == 0) strcpy(str,"HIS\0");
	if (strcmp(str,"CYX") == 0) strcpy(str,"CYS\0");
	if (strcmp(str,"MSE") == 0) strcpy(str,"MET\0");
	
	return;
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

	int atomnum, resnum;
	char s_atomnum[5], s_resnum[5];
	char type[6], atomname[5], resname[5], chain[2];
	double coordx, coordy, coordz, q, r, sigma_2, epsilon,sufix1,sufix2;
	int  s;
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
		ExtractString(4, &inStr[16], resname);
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

		if ((strncmp(type, "ATOM", 4) == 0) || (strncmp(type, "HETATM", 6) == 0)) {

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
void WriteModel(PDBData_t* mol_p, char* filename) {

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

	return;
}


/**************************************************************************
 * Function: atomID
 **************************************************************************/
int atomID(PDBData_t *mol_p, int atomid_old,PDBData_t *modller_p1)
{
	int i;
	int atom_id_new,res_count_n;
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
	return 1;
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
	int rot_bond_num;
	rot_bond_num=-100;
	if(strcmp(resname,"GLY")==0 || strcmp(resname,"ALA")==0 || strcmp(resname,"PRO")==0)
		rot_bond_num=1;
	else if(strcmp(resname,"ARG")==0 || strcmp(resname,"LYS")==0)
		rot_bond_num=4;
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

int  core_rim(PDBsas_t *com_p, PDBsas_t *mono_p)
{
	int i,j;
	double deltarsa,per,sc_delta,maxASA,maxASAsc,ratio1,ratio2,ratio3,ratio4, ratio5,delta_np,delta_p,np_per,p_per;
	FILE *fp;
	FILE *fp1;
	if((fp=fopen("core_inter.txt","a"))==NULL)
	{
		printf("could not open file :%s\n","core_inter.txt");
		exit(1);
	}
	if((fp1=fopen("rim_inter.txt","a"))==NULL)
	{
		printf("could not open file :%s\n","rim_inter.txt");
		exit(1);
	}

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
						maxASA = MaxASA(com_p->resolacs[i].rname);
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
							fprintf(fp,"%-4s%c%5d %9.2f %8.2f %8.2f %8.2f %9.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %9.2f %9.2f %9.2f %6.2f %9.2f %9.2f %9.2f %6.2f\n", com_p->resolacs[i].rname, com_p->resolacs[i].chain[0],com_p->resolacs[i].rnum, 
									com_p->resolacs[i].abs_tot, mono_p->resolacs[j].abs_tot, deltarsa, per, com_p->resolacs[i].abs_sc, 
									mono_p->resolacs[j].abs_sc, sc_delta, ratio1, ratio2, ratio3, ratio4, ratio5, com_p->resolacs[i].abs_npol,
									mono_p->resolacs[j].abs_npol,delta_np, np_per, com_p->resolacs[i].abs_pol, mono_p->resolacs[j].abs_pol,
									delta_p, p_per);
						}
						else if(deltarsa > EPSMIN && (per - 0.3) > EPSMIN)
						{
							fprintf(fp1,"%-4s%c%5d %9.2f %8.2f %8.2f %8.2f %9.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %9.2f %9.2f %9.2f %6.2f %9.2f %9.2f %9.2f %6.2f\n", com_p->resolacs[i].rname, com_p->resolacs[i].chain[0], com_p->resolacs[i].rnum, 
									com_p->resolacs[i].abs_tot, mono_p->resolacs[j].abs_tot, deltarsa, per, com_p->resolacs[i].abs_sc, 
									mono_p->resolacs[j].abs_sc, sc_delta, ratio1, ratio2, ratio3, ratio4, ratio5, com_p->resolacs[i].abs_npol,
									mono_p->resolacs[j].abs_npol,delta_np, np_per, com_p->resolacs[i].abs_pol, mono_p->resolacs[j].abs_pol,
									delta_p, p_per);
						}
						else
							continue;
					}
					else
					{
						maxASA = MaxASA(com_p->resolacs[i].rname);
						maxASAsc = MaxASAsc(com_p->resolacs[i].rname);
						deltarsa = fabs(com_p->resolacs[i].abs_tot - mono_p->resolacs[j].abs_tot);
						per = com_p->resolacs[i].abs_tot/mono_p->resolacs[j].abs_tot;
						sc_delta = fabs(com_p->resolacs[i].abs_sc - mono_p->resolacs[j].abs_sc);
						ratio1 = 100*com_p->resolacs[i].abs_tot/maxASA;
						ratio2 = 100*deltarsa/maxASA;
						ratio3 = 100*com_p->resolacs[i].abs_sc/maxASAsc;
						ratio4 = 100*sc_delta/maxASAsc;
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
							fprintf(fp,"%-4s%c%5d %9.2f %8.2f %8.2f %8.2f %9.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %9.2f %9.2f %9.2f %6.2f %9.2f %9.2f %9.2f %6.2f\n", com_p->resolacs[i].rname, com_p->resolacs[i].chain[0],com_p->resolacs[i].rnum, 
									com_p->resolacs[i].abs_tot, mono_p->resolacs[j].abs_tot, deltarsa, per, com_p->resolacs[i].abs_sc, 
									mono_p->resolacs[j].abs_sc, sc_delta, ratio1, ratio2, ratio3, ratio4, ratio5, com_p->resolacs[i].abs_npol,
									mono_p->resolacs[j].abs_npol,delta_np, np_per, com_p->resolacs[i].abs_pol, mono_p->resolacs[j].abs_pol,
									delta_p, p_per);
						}
						else if(deltarsa > EPSMIN && (per - 0.3) > EPSMIN)
						{
							fprintf(fp1,"%-4s%c%5d %9.2f %8.2f %8.2f %8.2f %9.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %9.2f %9.2f %9.2f %6.2f %9.2f %9.2f %9.2f %6.2f\n", com_p->resolacs[i].rname, com_p->resolacs[i].chain[0], com_p->resolacs[i].rnum, 
									com_p->resolacs[i].abs_tot, mono_p->resolacs[j].abs_tot, deltarsa, per, com_p->resolacs[i].abs_sc, 
									mono_p->resolacs[j].abs_sc, sc_delta, ratio1, ratio2, ratio3, ratio4, ratio5, com_p->resolacs[i].abs_npol,
									mono_p->resolacs[j].abs_npol,delta_np, np_per, com_p->resolacs[i].abs_pol, mono_p->resolacs[j].abs_pol,
									delta_p, p_per);
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
	fclose(fp1);
	return 1;
}

/*********************************************************************
Function: MaxASA
used to calculate the maxASA of the whole residue
the values come from the article "Quantifying the accessible surface area of protein residues in their local environment"
*********************************************************************/
double MaxASA(char *resname)
{
	double res_asa;
	if(strcmp(resname, "GLY")==0)
		res_asa = 70.27;
	else if(strcmp(resname, "ALA")==0)
		res_asa = 102.68;
	else if(strcmp(resname, "SER")==0)
		res_asa = 111.97;
	else if(strcmp(resname, "CYS")==0)
		res_asa = 127.72;
	else if(strcmp(resname, "CYX")==0)
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
    double sc_asa;
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
	else if(strcmp(resname, "MSE")==0)
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
	else if(strcmp(resname, "MSE")==0)
		sc_asa = 150.32;
	else if(strcmp(resname, "CYX")==0)
		sc_asa = 81.91;
	else if(strcmp(resname, "HIE")==0)
		sc_asa = 141.17;
	else if(strcmp(resname, "HID")==0)
		sc_asa = 141.17;
	else if(strcmp(resname, "HSD")==0)
		sc_asa = 141.17;
	else if(strcmp(resname, "HSE")==0)
		sc_asa = 141.17;
	else
	{
		printf("ERROR: this is not a natural residue!\n");
		exit(1);
	}
	return sc_asa;
}

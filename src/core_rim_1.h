#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include <errno.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.142
#endif

#define      VARS		(250)	/* max # of variables	*/     
#define      RHO_BEGIN	(0.5)	/* stepsize geometric shrink */
#define      EPSMIN		(1E-6)	/* ending value of stepsize  */
#define      IMAX		(5000)	/* max # of iterations	*/    

#define      MAX_REPLICAS 40 // Max number of replicas
#define      MAX_MODELS_PER_REPLICA 40 // Max number of models per replica
#define      MAX_SEGMENTS_PER_REPLICA 20 // Number of non-contiguous segments per replica

#define      MAX_ATOMS_PER_REPLICA 5000 // This is an overall number for all segments (but not all models)

#define		 MAX_NUM_CLUST 500

#define		 NVAR 6

/**** ENERGY ORDERING DEFS ****/
#define		 LJP 0
#define		 LJN 1
#define		 DS  2
#define		 EL  3
#define		 HB  4
#define		 PI  5
#define		 PC  6
#define		 IT  7
#define		 ETERMS  8


/**** RESIDUES ****/

#define		 ALA 0
#define		 ARG 1
#define		 ASN 2
#define		 ASP 3
#define		 CYS 4
#define		 GLN 5
#define		 GLU 6
#define		 GLY 7
#define		 HIS 8
#define		 ILE 9
#define		 LEU 10
#define		 LYS 11
#define		 MET 12
#define		 PHE 13
#define		 PRO 14
#define		 SER 15
#define		 THR 16
#define		 TRP 17
#define		 TYR 18
#define		 VAL 19



/**** TYPE DEFINITIONS ****/
typedef		double			mat[3][3];
typedef		double			matoms[MAX_MODELS_PER_REPLICA][3];

typedef struct {
	char		rname[5];	/* residue name */
	int			resid;		/* residue id */
	int			rnum;		/* residue number */
	int			firstatom;	/* index into atom array */
	int			ca;			/* position of CA atom */
	int			surface;	/* surface residue */
	int			interface;	/* interface residue */
	int			hbond;		/* hbond with clash */
	int			canchor[5];	/* chi angle anchor */
	double 		phi; 		/* phi angle */
	double		psi; 		/* psi angle */
	double		omg; 		/* omega angle */
	double		msf;
} residue_t;

typedef struct {
	int		myres;
	int		resnum;
	char	aname[5];	/* atom name */
	char	chain[2];	/* chain ID */
	double	xyz[3];		/* atom coordinates */
	double  msf;
	double  occ;        /* atom occupancy */
	double  bfact;      /* atom b-factor */
} atom_t;

typedef struct {
	int			natoms;			/* number of atoms in the molecule */
	int			nres;			/* number of residues in the molecule */
	atom_t*		atoms;			/* array of atom_t's */
	residue_t*	residues;		/* array of residue_t's */
} PDBData_t;

typedef struct {
	char rname[5];      /* residues' name */
	int  rnum;          /* residues' sequence number */
	char chain[2];      /* chain name */
	double abs_tot;     /* all atom absolute solvent accessable surface of the residue */
	double rel_tot;     /* all atom relative solvent accessable surface of the residue */
	double abs_sc;      /* side chain absolute solvent accessable surface of the residue */
	double rel_sc;      /* side chain relative solvent accessable surface of the residue */
	double abs_mc;      /* main chain absolute solvent accessable surface of the residue */
	double rel_mc;      /* main chain relative solvent accessable surface of the residue */
	double abs_npol;    /* non-polar atoms' absolute solvent accessable surface of the residue */
	double rel_npol;    /* non-polar atoms' relative solvent accessable surface of the residue */
	double abs_pol;     /* polar atoms' absolute solvent accessable surface of the residue */
	double rel_pol;     /* polar atoms' relative solvent accessable surface of the residue */
} rsa_t;

typedef struct {
	int nres;           /* number of residues */
	rsa_t* resolacs;    /* array of residues's solvent access surfaces */
} PDBsas_t;

typedef struct {
	double flexibility;  /* total rotatable single bonds number around the residue*/
	double nflex;        /* rotatable single bonds number averaged by side chain atoms number */
	double res_tot;      /* total number of residues around the residue */
	double atm_tot;      /* total number of atoms around the residue */
	double res_tot_wt;   /* weighted total number of residues around the residue */
	double hydrophobicity; /* total hydrophobicity of residues around the target residue */
} around_res_t;

/**** FUNCTION PROTOTYPES ****/                  

/**** replicopt_io ****/
FILE*			myfopen(char* filename, char* iotype);
PDBData_t*		ReadPDB(char* filename);
void			WriteModel(PDBData_t* receptor_p, char* filename);
void			ExtractString(int len, char *src, char *dst );
int				deblank(char *origStr);
int				isDigit(char A);
int				isH(char *A);
int				isNuc(char *sResName);
void			fixresname(char *str);
double			msf_calc(PDBData_t **mol_p, int nmodel);
int atomID(PDBData_t *mol_p, int atomid_old,PDBData_t *modller_p1);
double dist(double *coor,double *coor1);
int interRes(PDBData_t *mol_p, PDBData_t *mol_q);
int rot_sig_bond(char *resname);
double wt_rot_sig_bond(char *resname);
PDBsas_t* Read_SAS (char *sasfile);
double MaxASA(char *resname);
double MaxASAsc(char *resname);
//int  core_rim(PDBsas_t *com_p,PDBsas_t *mono_p,PDBData_t *mol_com);
int  core_rim(PDBsas_t *com_p,PDBsas_t *mono_p);
double resMW(char *resname);
double resMW_sc(char *resname);
int feature_around(rsa_t r_rsa, PDBData_t * mol_com, around_res_t *around_f, double dist);
double simi_character(char *res1name, char *res2name);
int classify(char *resname);
double resHP(char *resname);


/*
____________________________________________________________________________________
Title:		CloneIt (trade mark)
____________________________________________________________________________________
Version		2.1 (CGI version)
Date		1998
Langage:	ANSI-C
Author: 	Pierre LINDENBAUM
Adress:		Domaine de Vilvert INRA CRJ VIM
			Laboratoire de Biologie Moléculaire des rotavirus.
			78352 JOUY EN JOSAS FRANCE
E-Mail:		lindenb@biotec.jouy.inra.fr
____________________________________________________________________________________
Usage:		CloneIt is a molecular biology program finding strategies
			to sub-clone a DNA sequence (INSERT) into a plasmid (VECTOR) using
			restriction enzymes .
			This program handle parameters such as
				- usage of a phosphatase (CIP)
				- usage of modifying polymerases that blunt overhanged ends
				- partial digestions
				- 5' or/and 3' in-frame ligation
				- digestion post ligation
				- digestion to direct insert
				- detection of stop codon after ligation
			Another capabilities of this program are:
				- find excision of a DNA sequence by restriction enzyme that
				  will produce an in-frame deletion
				- find digestion, fill-in and ligation that produces a frameshift
					into sequence
				- compare 2 sequences from the point of view of restriction sites
				- display a restriction map.
				- Translation of restriction enzymes databases.
				- Management of large cloning projects.
____________________________________________________________________________________
Compilation:This program has been successfully tested on the GNU compiler gcc but it
			should be compatible on any ANSI-C compiler.
			To compile it, copy this file ( that should be named "CloneItCGI.c") on
			your UNIX directory.
			
			
			Please before compiling, look at this
	
			
			To compile ype:
			
			gcc -oCloneIt CloneIt.c
			
	
			this creates a program called "CloneIt".
			Type "CloneIt" to launch the program
____________________________________________________________________________________
Trade Mark:	National Number 97/704582
			Institut National de la Propriété Industrielle
			26, rue de St Petersbourg 75800 Paris Cedex 08 FRANCE
			Order Number:18.NOV1997
____________________________________________________________________________________
Copyright Notice: Permission to use, copy, modify, and distribute this software and
	its documentation is hereby granted, subject to the following restrictions
	and understandings:

		1) Any copy of this software or any copy of software derived
		from it must include this copyright notice in full.

		2) All materials or software developed as a consequence of the
		use of this software or software derived requires the express,
  		written author permission.

		3) The software may be used by anyone for any purpose, except
		that its redistribution for profit or inclusion in other
		software sold for profit requires the express, written
		permission of the author.
 
		4) This software is provided AS IS with no warranties of any
 		kind. In no event will the author be liable for any lost revenue 
 		or profits or other special, indirect and consequential damages. 
____________________________________________________________________________________
Revision history:
		Apr		1997 I learn C and i write a small piece of code
		Sept	1997 Final version submited to CABIOS
		Feb		1998 CABIOS revision (projects management...)
		Jun		1998 Cabios publication
		Aug		1998 Strategies are linked list; output to HTML, use GAP, Bestfit, 
					silent mutation,edit,DNA Strider Sequences,
		Jn		2002 define sequence_t
____________________________________________________________________________________
I would be very glad to receive suggestions and criticisms about this program from the
users...Pierre LINDENBAUM November 1998
____________________________________________________________________________________*/


	
	
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<errno.h>
#include<time.h>
#include <stdarg.h>

#define MAX_NPB					10000				/** Number max of bp within the sequence    **/
#define NOM_MAX_ENZYME			15					/** Length max of a string such as "EcoR1"  **/
#define PALINDROME_MAX_ENZYME	16					/** Length max of a string such as "GAATTC" **/
#define SITE_MAX_ENZYME			20					/** Length max of a string such as "G^AATTC" **/
#define MAX_SITE				100					/** Number max of sites for an enzyme       **/
#define MAX_LENGHT_POLY			30
#define LARGEUR_ECRAN			80					/** Length of the screen                    **/
#define HAUTEUR_ECRAN			50
#define BEEP 					
#define INKEY               	
#define BACK
#define DRAW_LINE	printf("<BR>\n");		
#define MAX_NOM_FICHIER			FILENAME_MAX					/** Length max of a string such as "MyFile" **/
#define VAR_VERSION				"CloneIt V2.2" 	/** Name of the application Please don't change this**/
#define POLYLINKER_FILE			"Polylinkers.Set"	/** Polylinkers File **/
#define VAR_CITATION			"LINDENBAUM Pierre (1998) CloneIt (tm):Finding Cloning strategies. BioInformatics Vol 14, 5 pp 465-466." 
#define VAR_URL					"http://locus.jouy.inra.fr/soft/cloneit/cloneit.html"
#define MAX(a,b)				((a)>(b)?(a):(b))
#define MIN(a,b)				((a)<(b)?(a):(b))
#define LOWER(c)				(((c)<=90 && (c)>=65)? (c)+32:(c))
#define UPPER(c)				(((c)>=97 && (c)<=122)?(c)-32:(c))
#define FORMAT_HTML				3
#define EOS						'\0'
#ifndef TRUE
	#define TRUE				1
	#define FALSE				0
#endif
#define AMBIGOUS				2					/* Nucleotide is N or Y or R or ....*/
#define ARE_IDENTIC				0
#define TYPE_PALINDROMIC		1					/* ex: G^AATTC      */
#define TYPE_ASYMETRIC			2					/* ex: GAATTC(1/-1) */
#define IS_IN_FRAME				0
#define VAR_OCHRE 				'*'
#define VAR_AMBRE 				'$'
#define VAR_OPALE 				'!'
#define TYPE_BLUNT				1 /*** SmaI  ***/
#define TYPE_3_OVER				2 /*** EcoR1 ***/
#define TYPE_5_OVER				3 /*** Pst1  ***/
#define NO_TREATMENT			1
#define KLENOW_TREATMENT		2
#define T4_TREATMENT			2
#define DNTP_TREATMENT			3
#define DELTA_TREATMENT			4
#define MAX_VAR_PARTIAL			3
#define VECTOR					0
#define INSERT					1
#define UNSELECT				2
#define IN_FRAME				1
#define NO_FRAME				0
#define PARAGRAF				printf("\n")
#define REBASE_CLASSIC			"classic_allenz"

#ifndef REBASE_FILE
	#define REBASE_FILE	 "allenz"
#endif

#ifndef DOC_URL
	#define DOC_URL "http://topaze.jouy.inra.fr/doc/"
#endif
	
#define _ABORT					2
#define GRAPHIC					50	/* length of graphical representation of insert */
#define SMALL_SITE				4 	/* length recognized by an enzyme that will be discarded (see optimizing memory in Main menu)*/
#define VEC_5	0
#define INS_5	1
#define VEC_3	0
#define INS_3	1
#define SIDE_5	0
#define SIDE_3	1
#define MAX_SAVE	10
#define MAX_LINE	200 /* max lenght line from file input */
#define RETURN_NULL				{err_printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);exit(0);}
#define ERROR_USER				{err_printf("Error from the user (file %s line %d)\n",__FILE__,__LINE__);exit(0);}

#define MAX_CORP	30 /* max number of society */
#define SUB_CLONING 		0
#define DELETION_FRAME_NO	1
#define DELETION_FRAME_T4	2
#define DELETION_CARBOXY	3
#define FRAME_SHIFT			4
#define BLUE_COLOR			"0000FF"
#define RED_COLOR			"FF0000"
#define YELLOW_COLOR		"00FF00"
#define DEFAULT_MAX_NUM_STGY 100
#define LINE_LEN 1024

 
/*
	define MENU : used in main(),
	it will be easer to modify the menu,
	to insert a new item,
	in the future with this system
*/



#define CGI_FILE_ERROR "Error.html"
#define PRINT_CGI_HEADER 	printf("Content-type: text/html%c%c",10,10)
#define SERVER_ERROR {  printf("<HTML><BODY><H1>SERVER ERROR</H1></BODY></HTML>");\
						exit(0);}
#define GET_NEXT_FIELD '&'
#define GET_EQUAL '='
#define HEAD_PAGE "head.html"
#define TAIL_PAGE "tail.html"
#define HOME_PAGE "cloneit_cgi.html"


#ifndef DATAS_PATH
	#define DATAS_PATH "./"
	#warning "env var DATAS_PATH undefined default is ./"
#endif

typedef int nbp_t;
typedef char bool_t;

/* this structure describes a sequence */
typedef struct _sequence_t
	{
	nbp_t	pos_ATG;
	nbp_t npb;
	nbp_t var_max;
	nbp_t var_min;
	int degenerate;
	nbp_t var_min_int,var_max_int;
	char sequence[MAX_NPB]; /* num bases sequence */
	char FICHIER_ADN[MAX_NOM_FICHIER]; /* names of sequence */
	} sequence_t;

static sequence_t seq[2];

/* this structure describes an enzyme */
typedef struct StructEnz
	{
	char 	nom[NOM_MAX_ENZYME];				/* name */
	char 	site[PALINDROME_MAX_ENZYME];		/* site GAATTC */
	char 	site_complet[SITE_MAX_ENZYME];		/* raw site G^AATTC */
	nbp_t 	taille_site;						/* lentgh of site: 6 */
	bool_t	select[2];		/* is this enzyme is selected for 0:VECTOR or 1:INSERT */
	nbp_t 	pos5_3;			/* distance on strand 5'->3' to the cuting site */
	nbp_t 	pos3_5;			/* distance on strand 3'->5' to the cuting site */
	bool_t	palindromic;	/* Is this Enzyme palindromic or not ? */
	nbp_t		Nbr_N5;			/* number of aNy base on side 5' of site ex: NNATG (used in function GetSite)*/
	nbp_t		Nbr_N3;			/* number of aNy base on side 3' of site ex: GTANN */
	short	corp;/* number of customers */
	} STRUCT_ENZYME;

/* this structure describes the preferences */
typedef struct Prefs
	{
	int			partial;		/* max number of partial digestions */
	int			var_pct_min;	/* percentage of insert that will be considered in the cloning box by default  */
	int			DeltaMax;		/* max percentage of insert  that will be deleted see:finding in frame deletion */
	int 		DeltaMin;		/* min percentage of insert  that will be deleted see:finding in frame deletion */
	bool_t 	anti_V; /* anti paralle VECTOR */
	bool_t 	anti_I; /* anti paralle INSERT */
	bool_t 	allow_T4;		/* allow modifying polymerase */
	bool_t 	allow_CIP;		/* allow phosphatase */
	bool_t		allow_all_sol;	/* show all solutions for a couple of enzymes */
	bool_t		allow_part_overhang; /* allow ligation between 2 sites that are compatibles but don't have the same overhang lenght */
	bool_t		display_messages; 
	bool_t 	side_5;	/*clone INSERT in frame in 5' */
	bool_t		side_3; /*clone INSERT in frame in 3' */
	bool_t		memory; /* discard short ,too degenerate, or ubiquitous enzymes */
	bool_t		partial_only_blunt; /* use partial digestion only if the enzyme cuts blunt */
	bool_t		search_C_term; /* search for C terminal deletion */
	int 		numMaxStgys;
	short 		strategy;
	bool_t		allenz_classic;
	char RebasePref[MAX_NOM_FICHIER];
	}STRUCT_PREFS;

/* this structure describes a site on a sequence */
typedef struct StructSite
	{
	int		NumEnz; /* Enzyme number see: STRUCT_ENZYME*/
	int		NumSeq; /* Sequence number 0:VECTOR or 1:INSERT*/
	int		Loc;	/* localisation on the sequence */
	} STRUCT_SITE;


/* this structure describes a cloning strategy on 5' side OR 3' side */
typedef struct StructStrategy_2
	{
	bool_t			Test_Trans;	/* TRUE: VECTOR and INSERT are considered FALSE: only INSERT */
	int				NumSite[2];	/* site number see:STRUCT_SITE*/
	int 			Treatment[2]; /*polymerase used*/
	int 			Partial[2];	/* number of partial sites */
	int				Pos5_3[2];	/* distance on strand 5'->3' from base 1 to the cuting site */
	int				Pos3_5[2];	/* distance on strand 3'->5' from base 1 to the cuting site */
	}STRUCT_STRATEGY_2;


/* this structure describes a cloning strategy  */
typedef struct StructStrategy
	{		
	STRUCT_STRATEGY_2	Couple[2]; /* STRUCT_STRATEGY_2 on side 5' and side 3' */
	bool_t				use_CIP;	/* is phosphatase used */
	struct StructStrategy      *next;
	struct StructStrategy      *prev;
	short type;
	}STRUCT_STRATEGY;
	
/* this structure describes the position of a site in function ResMap  */
typedef struct StructEnzxy
		{
		int		LocX; /* column number */
		int		LocY; /* line number */
		int		NumSite; /* site number see:STRUCT_SITE*/
		} STRUCT_ENZXY;
/* this structure describes the fragment released after enzymatic digestion  */
typedef struct fragment_struct
	{
	int Start;/* site number see:STRUCT_SITE*/
	int End;/* site number see:STRUCT_SITE*/
	}STRUCT_FRAGMENT;



/*____________________________________________________________________________________*/
  /*****************************************/
 /* Functions and procedures declarations */
/*****************************************/
bool_t AntiParallele(char*s ,nbp_t len);
void	err_printf(char *format, ... );
void	DRAW_HR(FILE *out,char car);
void	EchangeStgy(STRUCT_STRATEGY *A, STRUCT_STRATEGY *B);
void	SortStgys(void);
int 	Sign(int vara);
int 	Fct_Frame(int NumSeq,int nombre);
int		Fct_Pos(int NumSeq,int position);
int 	Fct_type(STRUCT_ENZYME Enzyme);
int 	Fct_est_ADN(char lettre);
char	Fct_Complementaire(char lettre);
int 	Fct_Identique(char lettre1,char lettre2);
void	Fct_Reverse(char *mot);
int 	Is_Codon_Stop(char AminoAcid);
char	Fct_Traduction(char c1,char c2,char c3);
char	Translation_at(short NumSeq,int pos);
void	ShowSeq(short mode,FILE *out,STRUCT_SITE site_5, STRUCT_SITE site_3);
int cgi_val(const char *string_number);
void cgi_get_sequence(char *the_seq,int NumSeq);
int 	ResMap(FILE *out,short mode,int var_seq);
bool_t	Intersections(FILE *out, short mode);
bool_t Are_Same_Asymetric(int NumSite_3,int NumSite_5);
int 	Cmd_Get_Info2(short mode, FILE *out,int NumEnz);
void	Cmd_Get_Info(short mode, FILE *out,int NumEnz);
void	Discard_Small(void);
void	Cmd_Get_Rebase(void);
int 	Add_Site(int k,int _NumSeq,int _Loc);
void	Cmd_Get_Site(void);
int 	Get_ATG(int NumSeq);
int 	Cmd_POLYLINKER2(int NumSeq,char mot1[],char mot2[],char mot3[]);
void	Cmd_POLYLINKER1(int NumSeq);
bool_t Are_in_Frame(STRUCT_STRATEGY_2 *Stg_2);
int 	Fct_Test(int NumSiteA,int NumSiteB,int Treatment ,int _Frame, STRUCT_STRATEGY_2 *Stg_2);
void	Find_stop_codon(short mode,FILE *out,int NumSite,int var_limit,int var_side);
int 	Fct_N_sites(int Numenz,int NumSeq,int _min,int _max,int Boolean);
void	Init_Enz(void);
void	Digest(short mode,FILE *out, int NumSeq,int NumEnz1, int NumEnz2, int AlertPartial);
bool_t Fct_Compatible2( int SeqA, int varA_5_3, int varA_3_5, int SeqB, int varB_5_3, int varB_3_5);
int 	Fct_Compatible_No_Treatment (int NumSiteA,int NumSiteB, STRUCT_STRATEGY_2 *Stg_2);
int 	Fct_Compatible(int NumSiteA,int NumSiteB,int Treatment,STRUCT_STRATEGY_2 *Stg_2);
void	Orient_Insert_After_CIP(STRUCT_STRATEGY *Strategy,short mode,FILE *out);
bool_t	Detect_Stop(FILE *out,short mode,int NumSiteA, int NumSiteB,STRUCT_STRATEGY_2 *Stg_2);
int 	Test_Use_CIP(int NumSeq,STRUCT_STRATEGY Stg);
void	Treatments(int *Treat_5, int *Treat_3,int Type5,int Type3);
int		Sub_Cloning(void);
int 	Fct_BaseShift(int NumSite, int Position);
int 	FrameShift(void);
void	Save_alignment( int NumSiteI_5, int partialI_5,int NumSiteI_3,int partialI_3,bool_t _TREATMENT,bool_t it_is_a_C_term_del ,FILE *out, short mode);
int 	Solution_C_Terminal(short mode, FILE *out,int NumSiteI_5,STRUCT_STRATEGY *Strategy);
int 	Solution_Frame(short mode, FILE *out,int NumSiteI_5, int partialI_5,int NumSiteI_3,int partialI_3, int _TREATMENT);
int		DeltaFrame(void);
void	Cmd_init_preferences(STRUCT_PREFS *Pref);


/* 21/06/98 */
bool_t TestStgy(void);
STRUCT_STRATEGY *DeleteStgy(STRUCT_STRATEGY *Stgy);
STRUCT_STRATEGY *NewStrategy(STRUCT_STRATEGY *TheNewStrategy);
STRUCT_STRATEGY *AddStgy(STRUCT_STRATEGY *Stgy);
STRUCT_STRATEGY *SetStgy(STRUCT_STRATEGY *Stgy);
STRUCT_STRATEGY *SetFirstStgy(void);
STRUCT_STRATEGY *SetLastStgy(void);
STRUCT_STRATEGY *GetCurrStgy(void);
STRUCT_STRATEGY *NextStgy(void);
STRUCT_STRATEGY *PrevStgy(void);
int CountStgy(int *pos);
int Make_New_enzyme(char *EnzymeName, char *EnzymeSite,int Enzyme_Available);
int Discard_iso(void);
void InitStgy(void);
void DeleteAllStgys(void);
void Show_All_Alignment(FILE *out,short mode);
void Show_All_Shift(FILE *out,short mode);
void Show_cloning_solution(FILE *out,STRUCT_STRATEGY *Strategy,short mode);
void SemiSolution(FILE *out,int NumSeq,STRUCT_STRATEGY *Stg,short mode);
void printHeader(FILE *out,STRUCT_STRATEGY *Strategy,short mode);
void printEnzymeName(int NumEnz);
void SolutionShift(short mode, FILE *out,STRUCT_STRATEGY *Strategy);
void printWebEnzyme(char *name,FILE *out);
void printTitle(char *title,short decal, bool_t begin);
void write_date(FILE *out);
void printMainHeader(FILE *out,short mode,bool_t begin);
void printEntreprise(FILE *out,char c);
bool_t SaveCloningSolutions(short mode);
bool_t SaveDeltaSolutions(short mode);
bool_t SaveShiftSolutions(short mode);
bool_t freadRebase(FILE *Buffer, char *mot1, char *mot2, char *mot3);
char x2c(char *what);
void cgi_print_page(char *filename);
void get_CloneIt_datas(char *word);
void get_CloneIt_fields(char *query);
bool_t cgi_bool(const char *string_number);


char			FICHIER_ENZYME[MAX_NOM_FICHIER]; /* name of rebase file */
STRUCT_PREFS	Preference; /* user's preferences */
int 	nbr_enzyme; /* number of enzymes in the  list of Enzymes*/
int 			true_nbr_enzyme; /* the true number of enzymes ( non palindromic count twice)*/
STRUCT_ENZYME   *Enzymes; /* a pointer to the list of Enzymes */
int 			nbr_sites; /* number of sites in memory */
STRUCT_SITE 	*Sites; /* a pointer to the list of sites */
char 			*string_seq[2]={"VECTOR","INSERT"};
STRUCT_ENZXY	*EnzXY; /* used in function Resmap */
int				var_nbr_xy=0; 	/* used in function Resmap */
int				NumFragment;/* used in Digestion */
STRUCT_FRAGMENT *TheFragment;/* a pointer used in Digestion */
bool_t			Search_Done=FALSE;
int				TotalSolutions=0;
/* 21/06/98 */
STRUCT_STRATEGY *FirstStgy=NULL;
STRUCT_STRATEGY *LastStgy=NULL;
STRUCT_STRATEGY *CurrStgy=NULL;
STRUCT_STRATEGY *NewStgy=NULL;
/* 9/11/98*/

#define NUM_CGI_FIELDS 27

const char *cgi_field[NUM_CGI_FIELDS]={
						"V_name=", /* 0 */
						"V_ATG=",/* 1 */
						"V_min=",/* 2 */
						"V_max=",/* 3 */
						"I_name=",/* 4 */
						"I_ATG=",/* 5 */
						"I_min=",/* 6 */
						"I_min_int=",/* 7 */
						"I_max_int=",/* 8 */
						"I_max=",/* 9 */
						"V_seq=",/* 10 */
						"I_seq=",/* 11 */
						"cip=",/* 12 */
						"frame5=",/* 13 */
						"frame3=",/* 14 */
						"partial=",/* 15 */
						"partial_blunt=",/* 16 */
						"short=",/* 17 */
						"klenow=",/* 18 */
						"carboxy=",/* 19 */
						"max_deletion=",/* 20 */
						"min_deletion=",/* 21 */
						"max_num=",/* 22 */
						"strategy=",/* 23 */
						"V_anti=",/* 24 */
						"I_anti=",/* 25 */
						"allenz_classic="/* 26 */
						};

char * local_query=NULL;
int len_local_query=0;
#define _GETENV_(a)  (char*)(a)
char *THE_QUERY_STRING="";






/*********************************************************************************/
char x2c(char *what)
	{
    char digit;
    
    digit = (what[0] >= 'A' ? ((what[0] & 0xdf) - 'A')+10 : (what[0] - '0'));
    digit *= 16;
    digit += (what[1] >= 'A' ? ((what[1] & 0xdf) - 'A')+10 : (what[1] - '0'));
    return(digit);
	}

/*********************************************************************************/
void get_CloneIt_datas(char *word)
	{
	int i=0,item=0;
	char *begin=NULL;
	bool_t found=FALSE;

	if(strcmp(word,"")==ARE_IDENTIC) return;
	if(strrchr(word,GET_EQUAL)==NULL)
		{
		err_printf("Bad input can't find '=' in the CGI query");
		//err_printf("Field Error: \"%s\" does not contain '%c'",word,GET_EQUAL);
		}
 	for(i=0;i<NUM_CGI_FIELDS;i++)
 		if(strncmp(cgi_field[i],word,strlen(cgi_field[i]))==ARE_IDENTIC)
 			{
 			begin=&word[strlen(cgi_field[i])];
 			found=TRUE;
 			item=i;
 			break;
 			}
 	if(found==FALSE)
	err_printf("Unknown field :'%s'",word);
 	switch(item)
 		{
 		case(0):strncpy(seq[VECTOR].FICHIER_ADN,begin,MAX_NOM_FICHIER);break;
 		case(1):seq[VECTOR].pos_ATG=cgi_val(begin);break;
 		case(2):seq[VECTOR].var_min=cgi_val(begin);break;
 		case(3):seq[VECTOR].var_max=cgi_val(begin);break;
 		case(4):strncpy(seq[INSERT].FICHIER_ADN,begin,MAX_NOM_FICHIER);break;
		case(5):seq[INSERT].pos_ATG=cgi_val(begin);break;
 		case(6):seq[INSERT].var_min=cgi_val(begin);break;
 		case(7):seq[INSERT].var_min_int=cgi_val(begin);break;
 		case(8):seq[INSERT].var_max_int=cgi_val(begin);break;
		case(9):seq[INSERT].var_max=cgi_val(begin);break;
		case(10):cgi_get_sequence(begin,VECTOR);break;
		case(11):cgi_get_sequence(begin,INSERT);break;
		case(12):Preference.allow_CIP=cgi_bool(begin);break;
		case(13):Preference.side_5=cgi_bool(begin);break;
		case(14):Preference.side_3=cgi_bool(begin);break;
		case(15):Preference.partial=cgi_bool(begin);break;
		case(16):Preference.partial_only_blunt=cgi_bool(begin);break;
		case(17):Preference.memory=cgi_bool(begin)==TRUE;break;
		case(18):Preference.allow_T4=cgi_bool(begin);break;
		case(19):Preference.search_C_term=cgi_bool(begin);break;
	    case(20):Preference.DeltaMax=cgi_val(begin);break;
		case(21):Preference.DeltaMin=cgi_val(begin);break;
		case(22):Preference.numMaxStgys=cgi_val(begin);break;
		case(23):Preference.strategy=cgi_val(begin);break;
		case(24):Preference.anti_V=cgi_val(begin);break; 
		case(25):Preference.anti_I=cgi_val(begin);break;
		case(26):Preference.allenz_classic=cgi_val(begin);break;
    	default:err_printf("Error :item=%d",item);break;
    	}
	}

/*********************************************************************************/
bool_t cgi_bool(const char *string_number)
	{
	int i=0;
	if(strncmp(string_number,"TRUE",4)==ARE_IDENTIC) return(TRUE);
	if(strncmp(string_number,"FALSE",5)==ARE_IDENTIC) return(FALSE);
	i=cgi_val(string_number);
	if(i!=FALSE && i!=TRUE)
		err_printf("Error:\"%s\" is not bool_t",string_number);
	return(i);
	}
/*********************************************************************************/
int cgi_val(const char *string_number)
	{
	register int i;
	int var_return=0,thelen=0;
	if(strspn(string_number," \n\r\t0123456789") != strlen(string_number))
		err_printf("Error:\"%s\" should be a base number",string_number);
	for(i=0;i<(int)strlen(string_number);i++)
		{
		if(string_number[i]==' ' || string_number[i]=='\n' || string_number[i]=='\r' || string_number[i]=='\t')
			continue;
		if(++thelen>6)
			err_printf("Error:\"%s\" is too large",string_number);
		var_return = var_return*10+string_number[i]-'0';
		}
	return(var_return);
	}
/*********************************************************************************/
void get_CloneIt_fields(char *query)
	{
 	int i=0,x=0,y=0;
 	char *word;
	char *begin=NULL;
	int the_len=0,query_len=0;
 	if(strrchr(query,GET_NEXT_FIELD)==NULL)
 		err_printf("Query Error: query does not contain '%c'",GET_NEXT_FIELD);
	query_len=strlen(query);
	begin=&query[0];

	for(i=0;i<=strlen(query);i++)
		{
		the_len++;
		if(query[i]==EOS && i==0) break;
		if(query[i]==GET_NEXT_FIELD || query[i]==EOS)
			{
	    	if((word=(char*)malloc((the_len)*sizeof(char)))==NULL)
	    		err_printf("Server Error: Out of memory (%d)",__LINE__);
	    	strncpy(word,begin,the_len);
	  		word[the_len]=EOS;
			/*enlever les '+' */
	    	for(x=0;x<strlen(word);x++)
	    		{
	    		if(word[x] == '+') word[x] = ' ';
	    		}
	    	/* enlever les signes ASCII */
	    	x=0;y=0;
	    	for(y=0;y<strlen(word);y++)
	    		{
	    		if(word[y]=='%')
	    			{
	    			word[x]=x2c(&word[y+1]);
	    			if(word[++y]==EOS) err_printf("Unespected end-of line in %s",word);
	    			if(word[++y]==EOS) err_printf("Unespected end-of line in %s",word);
	    			}
	    		else
	    			{
	    			word[x]=word[y];
	    			}
	    		x++;
	    		}
		    if(query[i]==GET_NEXT_FIELD)
		    	word[x-1]=EOS;
		    for(x=0;x<strlen(word);x++)
	    		{
	    		if(isspace(word[x])) word[x] = ' ';
	    		}
	    	/* decoder le champs */


	    	get_CloneIt_datas(word);
	    	free(word);
			begin=&query[i+1];
			the_len=0;  
			}
}
}

/*********************************************************************************/
void cgi_print_page(char *filename)
	{
	FILE *in=NULL;
	char line[LINE_LEN];
	if((in=fopen(filename,"r"))==NULL) SERVER_ERROR;
	while(fgets(line,LINE_LEN,in)!=NULL) puts(line);
	fclose(in);
	}
/*********************************************************************************/
void cgi_get_sequence(char *the_seq,int NumSeq)
	{
	char c;
	int index=0;
	/* intit sequence parameters */
	seq[NumSeq].npb=0;
	seq[NumSeq].sequence[0]='\0';
	seq[NumSeq].sequence[1]='\0';

	while((c=the_seq[index++])!= EOS )
		{
		if(!isspace(c) && c!=' ')
		if(c<48 || c>57)
			{
			if(Fct_est_ADN(c)>=TRUE)
				{
				seq[NumSeq].npb++;
				if(seq[NumSeq].npb<MAX_NPB)
					{
					seq[NumSeq].sequence[seq[NumSeq].npb]=UPPER(c);
					if(Fct_est_ADN(c)>=AMBIGOUS)
						{
						err_printf("Error: %s sequence contains degenerate base '%c' at position %d",
						 (NumSeq==VECTOR?"VECTOR":"INSERT"),c,seq[NumSeq].npb);	
						}
					}
				else
					err_printf("Error: %s sequence is too large\n",(NumSeq==VECTOR?"VECTOR":"INSERT"));	
				}
			else
				{
				err_printf("Error: In %s sequence: base '%c' at position %d is not a base symbol",
						 (NumSeq==VECTOR?"VECTOR":"INSERT"),c,seq[NumSeq].npb);	
				}
			}
		seq[NumSeq].sequence[seq[NumSeq].npb+1]='\0';
		}
	}

bool_t
AntiParallele(char*s, nbp_t len)
        {
        char* start=s;
        char* end=&s[len-1];
        char c;
        while(1)
        	{
        	if(start>end)
        		{
        		break;
        		}
        	else if(start==end)
        		{
        		*start=Fct_Complementaire(*start);
        		break;
        		}
        	else
        		{
	            c=*start; /* i=1->c=npb */ /* i=npb c=0 */
	            *start=Fct_Complementaire(*end);
	            *end=Fct_Complementaire(c);
        		}
        	++start;
        	--end;
        	}
        
        return(TRUE);
        }


/*********************************************************************************/
int main (int argc, char *argv[])
	{
	short i;
	FILE *in;
   	char *buff=NULL;
   	char *len1;
  	size_t contentlength;
  	char *endptr;
  	
	/*****************************************************/
	PRINT_CGI_HEADER;
	/*****************************************************/
	errno=FALSE;
	nbr_enzyme=0;
	nbr_sites=0;
	Search_Done=FALSE;
	strcpy(FICHIER_ENZYME,"");
	Cmd_init_preferences(&Preference);

	for(i=VECTOR;i<=INSERT;++i)
		{
		seq[i].degenerate=FALSE;
		seq[i].npb = 0;
		seq[i].var_max = 0;
		seq[i].var_min = 0;
		seq[i].pos_ATG = FALSE;
		seq[i].FICHIER_ADN[0]='\0';
		}
	InitStgy();
	/*****************************************************/
	cgi_print_page(DATAS_PATH HEAD_PAGE);
	/* check REQUEST_METHOD */
	


	if(getenv("CONTENT_LENGTH")==NULL || getenv("REQUEST_METHOD")==NULL)
            {
 			cgi_print_page(DATAS_PATH HOME_PAGE);
            cgi_print_page(DATAS_PATH TAIL_PAGE);
			exit(0);
			}
	len1 = (char*)getenv("CONTENT_LENGTH");
	contentlength=strtol(len1, &endptr, 10);
	if((buff=(char*)malloc((contentlength+2)*sizeof(char)))==NULL)
        {err_printf("Sorry:Out of memory.");}
	fread(buff, contentlength, 1, stdin);
	buff[contentlength]=EOS;
//err_printf("DEBUG buff %s\n",buff);
//DEBUG pos_atg, var_min, var_max, npb all zero here


	get_CloneIt_fields(buff);
	free(buff);
//DEBUG pos_atg 50K, var_min, var_max scaled between 0 and 100K, npb set here, scale them all
#if 0
err_printf("DEBUG atg, min, max, min_int, max_int, npb, name values: <br>seq0 %d %d %d %d %d %d %s <br>seq1 %d %d %d %d %d %d %s\n",
seq[0].pos_ATG, seq[0].var_min, seq[0].var_max, seq[0].var_min_int, seq[0].var_max_int, seq[0].npb, seq[0].FICHIER_ADN,
seq[1].pos_ATG, seq[1].var_min, seq[1].var_max, seq[1].var_min_int, seq[1].var_max_int, seq[1].npb, seq[1].FICHIER_ADN);
#endif
        /* the interface uses the "range" type and different browsers send different results if
           the user doesn't enter anything.  Internet Explorer sends back a 0, which is what the
           original code was expecting.  Firefox and Opera sends back the middle value of the range.
           First convert any IE style responses to Firefox style. The head.html has been modified
           to tell them to only use a browser that has a slider. */
        

        seq[0].pos_ATG     =      ((float) seq[0].pos_ATG    /100000.0) * seq[0].npb;
        seq[0].var_min     =  MAX(((float) seq[0].var_min    /100000.0) * seq[0].npb,1);
        seq[0].var_max     =      ((float) seq[0].var_max    /100000.0) * seq[0].npb;
        seq[0].var_min_int =      ((float) seq[0].var_min_int/100000.0) * seq[0].npb;
        seq[0].var_max_int =      ((float) seq[0].var_max_int/100000.0) * seq[0].npb;
        seq[1].pos_ATG     =      ((float) seq[1].pos_ATG    /100000.0) * seq[1].npb;
        seq[1].var_min     =  MAX(((float) seq[1].var_min    /100000.0) * seq[1].npb,1);
        seq[1].var_max     =      ((float) seq[1].var_max    /100000.0) * seq[1].npb;
        seq[1].var_min_int =      ((float) seq[1].var_min_int/100000.0) * seq[1].npb;
        seq[1].var_max_int =      ((float) seq[1].var_max_int/100000.0) * seq[1].npb;
        

	/* get rebase */
	
	if(Preference.allenz_classic==TRUE)
		{
		strcpy(FICHIER_ENZYME,DATAS_PATH REBASE_CLASSIC);
		}
	else
		{
		strcpy(FICHIER_ENZYME,DATAS_PATH REBASE_FILE);
		}

	if((in=fopen(FICHIER_ENZYME,"r"))!=NULL)
		{
		fclose(in);
		Cmd_Get_Rebase();
		}
	else
		err_printf("Server Error:Can't get REBASE \"%s\".",FICHIER_ENZYME);



	/* check datas */
	if(Preference.DeltaMax>100) err_printf("Error:deletion maximum overflow %d",Preference.DeltaMax);
	if(Preference.DeltaMin<1  ) err_printf("Error:deletion minimum overflow %d",Preference.DeltaMin);
	if(Preference.DeltaMin>=Preference.DeltaMax) err_printf("Error:deletion minimum/maximum overflow %d/%d",Preference.DeltaMin,Preference.DeltaMax);
	Preference.numMaxStgys=(Preference.numMaxStgys>100?100:Preference.numMaxStgys);
	Preference.numMaxStgys=(Preference.numMaxStgys<=0?1:Preference.numMaxStgys);
	if(seq[INSERT].npb<20 && (Preference.strategy!=4 && Preference.strategy!=5))
		err_printf("Error:INSERT sequence is not defined or is too short");

	if(seq[VECTOR].npb>=20)
		if(Preference.anti_V==TRUE) AntiParallele(&seq[VECTOR].sequence[1], seq[VECTOR].npb);
	if(seq[INSERT].npb>=20)
		if(Preference.anti_I==TRUE) AntiParallele(&seq[INSERT].sequence[1], seq[INSERT].npb);

	for(i=VECTOR;i<=INSERT;i++)
		{
		if(seq[i].npb>=20)
			{
			if(strcmp(seq[i].FICHIER_ADN,"")==ARE_IDENTIC)
				strcpy(seq[INSERT].FICHIER_ADN,(i==VECTOR?"VECTOR":"INSERT"));
// user set nothing, have to calculate values, was like this, but it doesn't work
//			if(seq[i].var_min==0 || seq[i].var_max==0) Cmd_POLYLINKER1(i);
//			if(seq[i].pos_ATG==0 ) Get_ATG(i);
			if(seq[i].var_min== seq[i].var_max) Cmd_POLYLINKER1(i);
                        int displacement = seq[i].pos_ATG - seq[i].npb/2;
			if(displacement >= -5 && displacement <= 5) Get_ATG(i);
			}
		}

#if 0
err_printf("DEBUG atg, min, max, min_int, max_int, npb, name values: <br>seq0 %d %d %d %d %d %d %s <br>seq1 %d %d %d %d %d %d %s\n",
seq[0].pos_ATG, seq[0].var_min, seq[0].var_max, seq[0].var_min_int, seq[0].var_max_int, seq[0].npb, seq[0].FICHIER_ADN,
seq[1].pos_ATG, seq[1].var_min, seq[1].var_max, seq[1].var_min_int, seq[1].var_max_int, seq[1].npb, seq[1].FICHIER_ADN);
#endif

	for(i=VECTOR;i<=INSERT;i++)
		{
		if(seq[i].npb>=20)
			{
		if(seq[i].pos_ATG<=0 || seq[i].pos_ATG>seq[i].npb) err_printf("Error: %s's ATG valor overflow (%d)",(i==0?"VECTOR":"INSERT"),seq[i].pos_ATG);
		if(seq[i].var_min<=0 || seq[i].var_min>seq[i].npb) err_printf("Error: %s's var_min valor overflow (%d)",(i==0?"VECTOR":"INSERT"),seq[i].var_min);
		if(seq[i].var_max<=0 || seq[i].var_max>seq[i].npb) err_printf("Error: %s's var_max valor overflow (%d)",(i==0?"VECTOR":"INSERT"),seq[i].var_max);
		if(seq[i].var_min>seq[i].var_max) err_printf("Error: %s's var_max &lt; var_min",(i==0?"VECTOR":"INSERT"));
			}
		}

        if(Preference.strategy!=4 && Preference.strategy!=5){
           if(seq[INSERT].var_min_int <= seq[INSERT].var_min){
               err_printf("Error: 5' internal cloning box limits are reversed: %d, %d (lower limit, upper limit)<BR>",seq[INSERT].var_min,seq[INSERT].var_min_int);
           }
           if(seq[INSERT].var_min_int >= seq[INSERT].var_max_int){
               err_printf("Error: 5' internal cloning box overlaps 3' cloning box:  %d, %d (upper of 5', lower of 3')<BR>",seq[INSERT].var_min_int,seq[INSERT].var_max_int);
           }
           if(seq[INSERT].var_max_int >= seq[INSERT].var_max){
               err_printf("Error: 3' internal cloning box limits are reversed: %d, %d (lower limit, upper limit)<BR>",seq[INSERT].var_max_int,seq[INSERT].var_max);
           }
        }


	if(Preference.memory==TRUE) Discard_Small();
	
	if(Preference.strategy<4)
		printf("<BR><CENTER><H2>Please wait, searching for cloning strategies, this may take a few minutes</H2></CENTER><BR>");
	
	switch(Preference.strategy)
		{
		case(1):
			{
			if(seq[VECTOR].npb<20)
				err_printf("Error:VECTOR sequence is not defined or is too short");
			Sub_Cloning();
			}
			break;
		case(2):
			{
			FrameShift();
			}break;
			
		case(3):
			{
			DeltaFrame();
			}break;
		case(4):
			{
			printMainHeader(stdout,FORMAT_HTML,TRUE);
			Intersections(stdout,FORMAT_HTML);
			printMainHeader(stdout,FORMAT_HTML,FALSE);
			}break;
		case(5):
			{
			if(Search_Done==FALSE) Cmd_Get_Site();
			printMainHeader(stdout,FORMAT_HTML,TRUE);
			printf("<P><CENTER><H2>Restriction Maps</H2></CENTER><P>");
			ResMap(stdout,FORMAT_HTML,VECTOR);
			printf("<HR>");
			ResMap(stdout,FORMAT_HTML,INSERT);
			printMainHeader(stdout,FORMAT_HTML,FALSE);
			}break;
		default:err_printf("Error: type of research %d unknown",Preference.strategy);break;
		}
	cgi_print_page(DATAS_PATH TAIL_PAGE);
	exit(0);
	return(0);
}

/***********************************************************************/
void Cmd_init_preferences(STRUCT_PREFS *Pref)
	{
	Pref->anti_V=FALSE;
	Pref->anti_I=FALSE;
	Pref->side_5=FALSE;
	Pref->side_3=FALSE;
	Pref->partial=FALSE;
	Pref->partial_only_blunt=FALSE;
	Pref->allow_CIP=FALSE;
	Pref->var_pct_min= 10;
	Pref->DeltaMax=80;
	Pref->DeltaMin=20;
	Pref->allow_T4=FALSE;
	Pref->memory=FALSE;
	Pref->allow_all_sol=FALSE;
	Pref->allow_part_overhang=FALSE;
	Pref->display_messages=FALSE;
	Pref->search_C_term=FALSE;
	Pref->allenz_classic=FALSE;
	Pref->numMaxStgys=DEFAULT_MAX_NUM_STGY;
	strcpy(Pref->RebasePref,REBASE_FILE);	
	}


void err_printf(char *format, ... )	
	{
	va_list ap;
	int result;
	FILE *in=NULL;
	char line[LINE_LEN];
	
	
	if((in=fopen(DATAS_PATH CGI_FILE_ERROR,"r"))==NULL)
		{
		printf("CloneIt: Server ERROR Can't open:"DATAS_PATH CGI_FILE_ERROR);
		}
	else
		{
		while(fgets(line,LINE_LEN,in)!=NULL)
			{
			if(strncmp(line,"£",1)==0)
				{
				va_start( ap, format );
				result = vprintf(format, ap );
				va_end( ap );
				}
			else
				puts(line);
			}
		fclose(in);
		}
	cgi_print_page(DATAS_PATH TAIL_PAGE);
	exit(0);
	}


/******************************************************************************/
void printEnzymeName(int NumEnz)
	{
	printWebEnzyme(Enzymes[NumEnz].nom,stdout);
	}

/******************************************************************************/
void printTitle(char *title,short decal, bool_t begin)
	{
	if(begin==TRUE) printf("<UL><LI><H%d>%s</H%d><P>",decal,title,decal);
		else printf("</UL><P>");
	}

/******************************************************************************/
void printMainHeader(FILE *out,short mode,bool_t begin)
	{
	char line[100];
	FILE *in=NULL;
	register int j;
	if(begin==TRUE)
		{
		fprintf(out,"<HR>Sequence file(s) used:<UL>");			
		if(seq[VECTOR].npb>0) fprintf(out,"<LI>%s %s (%d pb) [%d-%d] ATG frame at:%d <CODE>",
				"VECTOR",seq[VECTOR].FICHIER_ADN,seq[VECTOR].npb,seq[VECTOR].var_min,seq[VECTOR].var_max,seq[VECTOR].pos_ATG);
		if(seq[VECTOR].npb>0)	
		for(j=seq[VECTOR].pos_ATG;(j<=seq[VECTOR].var_max && (j-seq[VECTOR].pos_ATG)<=150);j+=3)
			{
			printf("%c",Translation_at(VECTOR,j));
			}
		printf("...</CODE>");
		if(seq[INSERT].npb>0) fprintf(out,"<LI>%s %s (%d pb) [%d-%d]-[%d-%d] ATG frame at:%d <CODE>",
				"INSERT",seq[INSERT].FICHIER_ADN,seq[INSERT].npb,seq[INSERT].var_min,seq[INSERT].var_min_int,seq[INSERT].var_max_int,seq[INSERT].var_max,seq[INSERT].pos_ATG);
		if(seq[VECTOR].npb>0)
		for(j=seq[INSERT].pos_ATG;(j<=seq[INSERT].var_max && (j-seq[INSERT].pos_ATG)<=150);j+=3)
			{
			printf("%c",Translation_at(INSERT,j));
			}
		printf("...</CODE>");
		/* rebase */
		fprintf(out,"</UL>REBASE file used:%s (%d enzymes)\n",FICHIER_ENZYME,true_nbr_enzyme);
	
		if((in=fopen(FICHIER_ENZYME,"r"))!=NULL)
			{
			printf("<BR><SMALL><PRE><TABLE><tr bgcolor=\"#e1e0ee\"><td>");
			while(fgets(line,99,in)!=NULL)
				{
				if(strstr(line,"allenz")!=NULL)
					{
					printf("\n%s",line);
					break;
					}
				}
			printf("</td></tr></table></PRE></SMALL><HR>");
			fclose(in);
			}
		else
			err_printf("Error:can't find %s.",FICHIER_ENZYME);
		/* parameters */
		fputs("Parameters:<UL>",out);
		
		if(Preference.memory!=FALSE) fputs("<LI>Speed optimization (Discard short sites).\n",out);
		if(Preference.allow_part_overhang!=FALSE) fputs("<LI>Allow non-overlapping overhangs.\n",out);
		if(Preference.side_5!=FALSE) fputs("<LI>Clone in frame at 5' of insert (NH2).\n",out);
		if(Preference.side_3!=FALSE) fputs("<LI>Clone in frame at 3' of insert (COOH).\n",out);
		fprintf(out,"<LI>Number of partial digestion allowed = %d.\n",Preference.partial);
		if(Preference.partial_only_blunt!=FALSE) fputs("<LI>Allow partial digestion only if enzyme blunt.\n",out);
		fprintf(out,"<LI>%sllow usage of C.I.A.P.\n",(Preference.allow_CIP==TRUE?"A":"Do NOT a"));
		fprintf(out,"<LI>%sllow usage of modifying polymerase.\n",(Preference.allow_T4==TRUE?"A":"Do NOT a"));
		fprintf(out,"<LI>Deletion length comprised between %d%% and %d%% of the insert.\n",Preference.DeltaMin,Preference.DeltaMax);
		fprintf(out,"<LI>%sook for Carboxy terminal deletions.\n",(Preference.search_C_term==TRUE?"L":"Do NOT l"));
		fprintf(out,"<LI>Maximum number of strategies displayed=%d.\n",Preference.numMaxStgys);
		fputs("</UL><HR><PRE>",out);
		}
	else
		{
		printTitle("Miscellanous",1,TRUE);
		fprintf(out,"</PRE>");
		cgi_print_page(DATAS_PATH "misc.html");
		printTitle("",1,FALSE);	
		}
	}







/*******************************************************************************
 get the most probable atg frame (0,1 or 2) in the sequence
		that is to say the length max of sequence where there is no stop codon
********************************************************************************/
int Get_ATG(int NumSeq)
	{
	int i,j,vara,_max=0;
	
	/* no boundaries are defined */
	if(seq[NumSeq].var_min==seq[NumSeq].var_max) err_printf("Boundaries not defined in seq %d.",NumSeq);
	seq[NumSeq].pos_ATG=FALSE;

		if(NumSeq==INSERT)
			{
			/* if the sequence is insert, the program is checking the frame, where there is
				the longest ORF without stop codon IN the insert */
			for(i=0;i<=2;i++)
				{
				vara=0;
				for(j=seq[INSERT].var_min+i;j<=seq[INSERT].var_max;j=j+3)
					{
					if(Is_Codon_Stop(Fct_Traduction(seq[INSERT].sequence[Fct_Pos(INSERT,j)],
													seq[INSERT].sequence[Fct_Pos(INSERT,j+1)],
													seq[INSERT].sequence[Fct_Pos(INSERT,j+2)]))!=TRUE)
						{
						vara++;
					
						if(vara>_max)
							{
							_max=vara;
							seq[INSERT].pos_ATG=seq[INSERT].var_min+i;
							}
						}
					else
						vara=0;
					}
				}
			}
		else if(NumSeq==VECTOR)
			{
			/* if the sequence is vector, the program is checking the frame, where there is
				the longest ORF without stop codon from each side (but not in ) of the cloning box */
			for(i=0;i<=2;i++)
				{
				/* search on the left of cloning box */
				vara=0;
				for(j=seq[VECTOR].var_min+i;j>=1;j=j-3)
					{
					if(Is_Codon_Stop(Fct_Traduction(seq[VECTOR].sequence[Fct_Pos(VECTOR,j)],
													seq[VECTOR].sequence[Fct_Pos(VECTOR,j+1)],
													seq[VECTOR].sequence[Fct_Pos(VECTOR,j+2)]))==TRUE)
						break;
					else
						vara++;
					}
				if(vara>_max)
					{
					
					_max=vara;
					seq[VECTOR].pos_ATG=seq[VECTOR].var_min+i;
					}
				/* search on the right of cloning box */
				vara=0;
				for(j=seq[VECTOR].var_max+i;j<=seq[VECTOR].npb;j=j+3)
					{
					if(Is_Codon_Stop(Fct_Traduction(seq[VECTOR].sequence[Fct_Pos(VECTOR,j)],
													seq[VECTOR].sequence[Fct_Pos(VECTOR,j+1)],
													seq[VECTOR].sequence[Fct_Pos(VECTOR,j+2)]))==TRUE)
						break;
					else
						vara++;
					}
				if(vara>_max)
					{
					_max=vara;
					seq[VECTOR].pos_ATG=seq[VECTOR].var_max+i;
					}
				/* search IN the cloning box */
				vara=0;
				for(j=seq[VECTOR].var_min+i;j<=seq[VECTOR].var_max;j=j+3)
					{
					if(Is_Codon_Stop(Fct_Traduction(seq[VECTOR].sequence[Fct_Pos(VECTOR,j)],
													seq[VECTOR].sequence[Fct_Pos(VECTOR,j+1)],
													seq[VECTOR].sequence[Fct_Pos(VECTOR,j+2)]))==TRUE)
						break;
					else
						vara++;
					}
				if(vara>_max)
					{
					
					_max=vara;
					seq[VECTOR].pos_ATG=seq[VECTOR].var_min+i;
					}
				}
			}
		
	return(TRUE);
	}

/********************************  Cmd_POLYLINKER2  ******************/
/* this function try to find the words mot2 and mot3 in the sequence */
/* this will be the boundaries of the cloning box                    */
/*********************************************************************/
int Cmd_POLYLINKER2(int NumSeq,char mot1[],char mot2[],char mot3[])
	{
	int var_return=FALSE,vara=TRUE,poly_left=0,poly_right=0; 
	char *result;
	char *p_seq;
		
	p_seq=&seq[NumSeq].sequence[1];
	result=strstr(p_seq, mot2);
	if(result!=NULL)
		{
		vara=FALSE;
		poly_left=seq[NumSeq].npb-strlen(result)+1+strlen(mot2);
		}
	if(vara==FALSE)
		{
		p_seq=&seq[NumSeq].sequence[poly_left];
		result=strstr(p_seq, mot3);
		if(result!=NULL)
			{
			poly_right=seq[NumSeq].npb-strlen(result)+1;
			}
		}
	if(poly_left<poly_right && poly_left!=0 && poly_right!=0)
		{
		var_return=TRUE;
		seq[NumSeq].var_min=poly_left;
		seq[NumSeq].var_max=poly_right;
		
		printf("<HR>The %s sequence seems to be a <B>\"%s\"</B> type plasmid. ",(NumSeq==VECTOR?"VECTOR":"INSERT"),mot1);
		
		if(NumSeq==VECTOR)
			{
			printf("The <A HREF=\"%sCloneIt_Help.html#cloning_box\">cloning box</A> may be [<B>%d</B> - <B>%d</B>]",DOC_URL,seq[NumSeq].var_min,seq[NumSeq].var_max);
			}
		else
			{
			seq[INSERT].var_min_int=seq[INSERT].var_min+(int)((seq[INSERT].var_max-seq[INSERT].var_min)/Preference.var_pct_min);
			seq[INSERT].var_max_int=seq[INSERT].var_max-(int)((seq[INSERT].var_max-seq[INSERT].var_min)/Preference.var_pct_min);
			printf("The <A HREF=\"%sCloneIt_Help.html#cloning_box\">cloning boxes</A> may be  <B>[%d-%d]</B> and <B>[%d-%d]</B>",DOC_URL,seq[NumSeq].var_min,seq[INSERT].var_min_int,seq[INSERT].var_max_int,seq[NumSeq].var_max);
			}
		printf("<BR>");
		fflush(stdout);
		}
	return(var_return);
	}
	
/********************************  Cmd_POLYLINKER1  ***/
/* Look in the polylinker.set file for cloning box boundaries */
/*** POLYLINKER_FILE example:
		;
		; List of Polylinkers.
		
		MyPlasmid,	gacgtc,gaagaggagag,
		
****/
void Cmd_POLYLINKER1(int NumSeq)
	{
	FILE *tampon;
	/* mot 1: name of plasmid        mot2 & mot3:oligo used */
	char mot1[NOM_MAX_ENZYME],mot2[MAX_LENGHT_POLY],mot3[MAX_LENGHT_POLY];
	char c='\0';
	int vara=0,varb=0;
	char ThePolyfFile[FILENAME_MAX];
	

	sprintf(ThePolyfFile,"%s",DATAS_PATH POLYLINKER_FILE);
		
	
	Search_Done=FALSE;
	/* open POLYLINKER_FILE */
	tampon = fopen(ThePolyfFile, "r" );
	if (tampon==NULL) 
		{
		err_printf("SERVER ERROR: Can't find polinker file");
		}
	else
		{
		while(c != EOF )
			{
			/* read file */
			c=fgetc(tampon);
			if(c==';')
				{
				while((c=fgetc(tampon))!='\n' && c!='\r' && c!=21 && c!=8) {};
				vara=0;
				varb=0;
				} /* discard commentary */
			if(c != EOF )
				{
				if(c==',') {vara++;varb=0;} /* next field is found */
				else if(c=='\n' || c=='\r' || c==21 || c==8) /* end of line: polylinker can be searched */
					{
					if(vara==3) /* if mot1 & mot 2 & mot3 found */
						{
						Fct_Reverse(mot3); /* inverse mot3 for its search in Cmd_POLYLINKER1*/
						if (Cmd_POLYLINKER2(NumSeq,mot1,mot2,mot3)==TRUE)
							break;
						else /* if those oligos are not found in direct , try to find them in reverse */
							{
							Fct_Reverse(mot3);
							Fct_Reverse(mot2);
							if (Cmd_POLYLINKER2(NumSeq,mot1,mot3,mot2)==TRUE) break;
							}
							
						}
					vara=0;varb=0;
					}
				else if(varb<MAX_LENGHT_POLY)
					{
					if(Fct_est_ADN(c)>=AMBIGOUS && (vara==1 || vara==2))
						{
						err_printf("Server error: polylinker %s should not contain degenerate bases.",mot1);
						}
					if (vara==0) {mot1[varb++]=c;mot1[varb]='\0';} /* complete name of plasmid */
					if (vara==1 && Fct_est_ADN(c)>0 && varb<MAX_LENGHT_POLY)
						{mot2[varb++]=UPPER(c);mot2[varb]='\0';} /* complete oligo1 */
					if (vara==2 && Fct_est_ADN(c)>0 && varb<MAX_LENGHT_POLY)
						{mot3[varb++]=UPPER(c);mot3[varb]='\0';} /* complete oligo1 */
					}
				else if(varb>=MAX_LENGHT_POLY)
					{
					err_printf("Reduce your oligonucleotide length (%s).\n",mot1);
					}
				}
			}
		fclose(tampon);
		}
	if(seq[NumSeq].var_min==seq[NumSeq].var_max) /* if polylinker is not defined set boundaries to all the sequence*/
			{
			seq[NumSeq].var_min=1;
			seq[NumSeq].var_max=seq[NumSeq].npb;
			printf("<H1>%s can't find the cloning box(es). Please define the boundaries</H1><BR>",VAR_VERSION); 
			printf("Here is the polylinkers file used.<P><table border=0 cellpadding=2 cellspacing=1><tr bgcolor=\"#e1e0ee\"><td><PRE>\n");
			cgi_print_page(DATAS_PATH POLYLINKER_FILE);
			printf("</PRE></td></tr></table><BR><A HREF=\"mailto:lindenb@biotec.jouy.inra.fr\">You can ask for us for adding a new item</A><HR>\n");
			err_printf("Error:can't find the cloning box(es).");
			}
	}


void printEntreprise(FILE *out,char c)
	{
	switch(c)
		{
		case('A'):fprintf(out,"<LI><A HREF=\"http://www.apbiotech.com/\">Amersham Life Sciences-USB</A>");break;
		case('B'):fprintf(out,"<LI><A HREF=\"http://www.lifetech.com\">Life Technologies Inc, Gibco-BRL</A>");break;
		case('C'):fprintf(out,"<LI>Minotech Molecular Biology Products");break;
		case('D'):fprintf(out,"<LI>Angewandte Gentechnologie Systeme");break;
		case('E'):fprintf(out,"<LI><A HREF=\"http://www.strategene.com\">Stratagene</A>");break;
		case('F'):fprintf(out,"<LI><A HREF=\"http://www.fermentas.com\">Fermentas AB</A>");break;
		case('G'):fprintf(out,"<LI><A HREF=\"http://www.oncor.com/prod-app.htm\">Appligene Oncor</A>");break;
		case('H'):fprintf(out,"<LI>American Allied Biochemical, Inc.");break;
		case('I'):fprintf(out,"<LI>SibEnzyme Ltd.");break;
		case('J'):fprintf(out,"<LI>Nippon Gene Co., Ltd.");break;
		case('K'):fprintf(out,"<LI><A HREF=\"http://www.takara.co.jp/\">Takara Shuzo Co. Ltd.</A>");break;
		case('L'):fprintf(out,"<LI>Kramel Biotech");break;
		case('M'):fprintf(out,"<LI><A HREF=\"http://www.boehringer-mannheim.com\">Boehringer-Mannheim</A>");break;
		case('N'):fprintf(out,"<LI><A HREF=\"http://www.neb.com\">New England BioLabs</A>");break;
		case('O'):fprintf(out,"<LI>Toyobo Biochemicals");break;
		case('P'):fprintf(out,"<LI><A HREF=\"http://www.apbiotech.com/\">Pharmacia Biotech Inc.</A>");break;
		case('Q'):fprintf(out,"<LI><A HREF=\"http://www.chimerx.com\">CHIMERx</A>");break;
		case('R'):fprintf(out,"<LI><A HREF=\"http://www.promega.com\">Promega Corporation</A>");break;
		case('S'):fprintf(out,"<LI><A HREF=\"http://www.sigma.sial.com/\">Sigma</A>");break;
		case('T'):fprintf(out,"<LI><A HREF=\"http://www.adbio.co.uk/\">Advanced Biotechnologies Ltd.</A>");break;
		default:fprintf(out,"<LI>? Unknow '%c'.",c);break;	
		}
	}


/************************************************************************************/
void printWebEnzyme(char *name,FILE *out)
	{
	int i=0;
	fprintf(out,"<A HREF=\"http://rebase.neb.com/rebase/enz/");
	while(name[i]!='\0')
		{
		if(name[i]!=' ')
			fprintf(out,"%c",name[i]);
		i++;
		}
	fprintf(out,".html\">%s</A>\n",name);
	}


/************************************************************************************/
void Cmd_Get_Info(short mode, FILE *out,int NumEnz)
	{
	/************************************************************************************
	This procedure gets informations about an Enzyme from the Rebase File...
	____________________________________________________________________________________*/
	fprintf(out,"<LI><H2><A NAME=\"TAG%d\"> Informations about %s [%s].</H2>\n",
			Sites[NumEnz].NumEnz,
			Enzymes[Sites[NumEnz].NumEnz].nom,Enzymes[Sites[NumEnz].NumEnz].site_complet);
	
	
	/**Create a link to the rebase www site***/
	fprintf(out,"Internet WWW link:");
	 printWebEnzyme(Enzymes[Sites[NumEnz].NumEnz].nom,out);
	if(seq[VECTOR].npb>0 && Fct_N_sites(Sites[NumEnz].NumEnz,VECTOR,1,seq[VECTOR].npb,TRUE)) {fprintf(out,"VECTOR Pattern ");Digest(mode,out,VECTOR,NumEnz,NumEnz,FALSE);}
	if(seq[INSERT].npb>0 && Fct_N_sites(Sites[NumEnz].NumEnz,INSERT,1,seq[INSERT].npb,TRUE)) {fprintf(out,"INSERT Pattern ");Digest(mode,out,INSERT,NumEnz,NumEnz,FALSE);}
	Cmd_Get_Info2(mode,out,Sites[NumEnz].NumEnz);
	}


int Cmd_Get_Info2(short mode, FILE *out,int NumEnz)
{
/************************************************************************************
This procedure gets informations about an Enzyme from the Rebase File...
____________________________________________________________________________________*/
/******************************************************/
FILE	*Buffer;
char	c2,mot1[NOM_MAX_ENZYME],mot2[NOM_MAX_ENZYME],mot3[MAX_CORP]="";
int		var_find=0,fric=0,find_iso=FALSE,k;
char line[MAX_LINE+1];
  /*******************************************/
 /* first find information about the enzyme */
/*******************************************/

if ((Buffer = fopen(FICHIER_ENZYME, "r" ))==NULL)
	{
	err_printf("Error: Server can't acces RE file %s.\n",FICHIER_ENZYME);
	}
else
	{
	while(fgets(line,MAX_LINE,Buffer) != NULL )
		{
		fric=FALSE;/* set to non commercialy available */
		/* search the tag < and > that delimit Rebase file fields */
		if (line[0]!='<' || line[2]!='>')
			continue;
		if(line[strlen(line)-1]=='\n' || line[strlen(line)-1]=='\r')
			line[strlen(line)-1]='\0';
		c2=line[1];
		if (c2=='1') /* Get the Enzyme's name */
				{
				var_find=1; /* name  found */
				memmove(&mot1[0],&line[3],MIN(strlen(line)-2,SITE_MAX_ENZYME));
				mot1[MIN(strlen(line)-3+1,SITE_MAX_ENZYME)]=EOS;
				if(strcmp(mot1,Enzymes[NumEnz].nom)!=ARE_IDENTIC)
					var_find=FALSE;
				}
				/* prototype */
		else if (var_find==1)
			{
			switch(c2)
				{
				case('2'):fprintf(out," Prototype                :%s\n",&line[3]);break;
				case('3'):fprintf(out," Microorganism            :%s\n",&line[3]);break;
				case('4'):fprintf(out," Source                   :%s\n",&line[3]);break;
				case('6'):fprintf(out," Methylation site         :%s\n",&line[3]);break;
				case('7'):fprintf(out," Commercial availability  :\n");break;
				case('8'):fprintf(out," New England Biolabs Refs :%s\n",&line[3]);break;
				default:break;
				}
			if(c2=='7')
				{
				fprintf(out,"<UL>");
				for(k=3;k<strlen(line);k++)
					printEntreprise(out,line[k]);
				fprintf(out,"</UL>");
				fprintf(out,"<BR>");
				}
			if(c2=='8') break;
			}/* end var_find */
		}/*end while*/
	/**************************************
	Now search isoschyzomers in Rebase File
	***************************************/
	fprintf(out,"Looking for Isoschizomers.\n");
	fprintf(out,"<UL>");
	find_iso=FALSE;
	fseek(Buffer,0L,SEEK_SET);/* redo from start */
	while(freadRebase(Buffer,mot1,mot2,mot3)==TRUE)
		{
		if(strcmp(mot1,Enzymes[NumEnz].nom)==ARE_IDENTIC) continue;
		if(strcmp(mot2,Enzymes[NumEnz].site_complet)!=ARE_IDENTIC)  continue;
		if(strlen(mot3)<=0) continue;
		find_iso=TRUE;
		fprintf(out,"<LI>%s  %s.\n",mot1,mot2);

		if(strlen(mot3)>0)
			{
			fprintf(out,"<LI>salso available at:\n");
			fprintf(out,"<UL>");

			for(fric=0;fric<strlen(mot3);fric++)
				{
				
				printEntreprise(out,mot3[fric]);
				}
			fprintf(out,"</UL>");
			}
		}
	if(find_iso==FALSE)
		fprintf(out,"\tNo other enzyme was found.\n");
	fprintf(out,"</UL>");
	fclose(Buffer);
	}/* end if file */
return(TRUE);
}



/************************************************************************************/
bool_t Are_Same_Asymetric(int NumSite_3,int NumSite_5)
	{
	/* this function is useful to detect two enzyme that are non plaindromic but
		have the same name, such enzymes are generated by Cmd_Get_Rebase */

	if(Sites[NumSite_3].NumEnz == Sites[NumSite_5].NumEnz)
		return(TRUE);
	
	if(strcmp(	Enzymes[Sites[NumSite_3].NumEnz].site_complet,
				Enzymes[Sites[NumSite_5].NumEnz].site_complet)==ARE_IDENTIC)
		return(TRUE);
	return(FALSE);
	}

/************************************************************************************
	as many short sites have a very low probabilty to be used in a clonage, because
	they cut sequence  everywhere, this function will be used to discard
	enzymes like  Taq I [t^cga] or NlaIV[ggn^ncc].
	This function also always discard isoschizomers, the less commercialy available will be discarded
____________________________________________________________________________________*/

void Discard_Small(void)
{
int i,j;
double _num,vara=FALSE;/* number of bases cuted by an enzyme */


/* scan all enzymes */
for(i=0;i<nbr_enzyme;i++)
	{
	vara=(double)(Enzymes[i].taille_site);
	if(Preference.memory==TRUE)
		{
		/* take in charge in _num the degenerate bases in the site */
		for(j=0;j<Enzymes[i].taille_site;j++)
			if( (_num=(double)Fct_est_ADN(Enzymes[i].site[j] )) >= AMBIGOUS)
				vara=vara-(_num/4.0);
		}
		
	if(vara<=SMALL_SITE) /* if an enzyme was discarded then move the array */
		{
		if(Enzymes[i].palindromic!=TRUE && i+1<nbr_enzyme)
			{
			if(strcmp(Enzymes[i].site_complet,Enzymes[i+1].site_complet)==ARE_IDENTIC)
				true_nbr_enzyme++;/* else true_number will be sub-estimated */
			}
			
		for(j=i;j<nbr_enzyme-1;j++)
			{
			memcpy(&Enzymes[j],&Enzymes[j+1],sizeof(STRUCT_ENZYME));
			}
		i--;
		true_nbr_enzyme--;
		nbr_enzyme--; /* deleted last enzyme and realloc memory */
		Enzymes=(STRUCT_ENZYME*)realloc(Enzymes,(nbr_enzyme)*sizeof(STRUCT_ENZYME));
		}
	}
}
/************************************************************************************/
bool_t freadRebase(FILE *Buffer, char *mot1, char *mot2, char *mot3)
	{
	/*int Enzyme_Available=0;*/
	char line[MAX_LINE];
	int var_find=0;
	int i;
	
	while(fgets(line,MAX_LINE,Buffer) != NULL )
		{
		/* Enzyme_Available=FALSE;set to non commercialy available */
		/* search the tag < and > that delimit Rebase file fields */
		if (line[0]!='<' || line[2]!='>')
			continue;
		
		if(line[strlen(line)-1]=='\n' || line[strlen(line)-1]=='\r')
			line[strlen(line)-1]=EOS;
			
		/*purge non authorized caracters */
		for(i=0;i<strlen(line);i++) if(line[i]=='\t') line[i]=' ';
			
		if (line[1]=='1') /* Get the Enzyme's name */
			{
			var_find=1; /* enzyme name  found */
			memmove(&mot1[0],&line[3],MIN(strlen(line)-2,SITE_MAX_ENZYME));
			mot1[MIN(strlen(line)-3+1,SITE_MAX_ENZYME)]=EOS;
			}

		else if (line[1]=='5') /* Get the Enzyme's NewEnzyme.site of cleavage */
			{
			if(strchr(line,'?')!=NULL || strpbrk(line,"(/)^")==NULL)
				{ var_find=FALSE;}
			var_find++; /* enzyme name  found */
			memmove(&mot2[0],&line[3],MIN(strlen(line)-2,SITE_MAX_ENZYME));
			mot2[MIN(strlen(line)-3+1,SITE_MAX_ENZYME)]=EOS;
			}
		else if (line[1]=='7') /** find if the Enzyme is commercialy available **/
			{
			var_find++; /* commercial availability  found */
			memmove(&mot3[0],&line[3],MIN(strlen(line)-2,MAX_CORP));
			mot3[MIN(strlen(line)-3+1,MAX_CORP)]=EOS;
			if(var_find==3) 
				{
				/*printf("mot1=%s mot4=%d mot5=%c mot6=%d/%d/%d/%d\n",
									mot1,*mot4,*mot5,mot6[0],mot6[1],mot6[2],mot6[3]);
				fflush(stdout);*/
				
				return(TRUE);
				}
			else
				var_find=0;
			}
		}
	return(FALSE);
	}

/************************************************************************************/
int Make_New_enzyme(char *EnzymeName, char *EnzymeSite,int Enzyme_Available)
	{
	int i;
	int j,vara,varb;
	int var_parenthese_droite=0,var_parenthese_gauche=0,var_slash=0,nbr_parenthese=0;
	int pos5_3,pos3_5,taille_site=0;
	char site[PALINDROME_MAX_ENZYME];
	bool_t var_non_palindromic=0;
	/* realloc memory for the new enzyme */	Enzymes=(STRUCT_ENZYME*)realloc(Enzymes,(nbr_enzyme+1)*sizeof(STRUCT_ENZYME));
	if(Enzymes==NULL)
		{
		err_printf("Server OUT OF MEMORY !  (too many enzymes ! %d )****\n",nbr_enzyme);
		}
	else
		{
		/**** Analyse the site *****************************************/
		i=j=0;
		vara=FALSE;
		nbr_parenthese=0;
		
		site[0]=EOS;
		var_non_palindromic=(strpbrk(EnzymeSite,"(/)")!=NULL?TRUE:FALSE);/*** It's an asymetric site not palindromic ***/
		
		/*** It's a palindromic site ***********************************************************/
		if(var_non_palindromic==FALSE)	 
			{
			if(strchr(EnzymeSite,'^')==NULL)
				{
				err_printf("Server error %s-%s",EnzymeName,EnzymeSite);
				}
			/***Find the cleavage site  ^ *******/
			i=j=0;
			pos5_3=-1;
			while(EnzymeSite[i]!='\0')
				{
				if (Fct_est_ADN(EnzymeSite[i])>=TRUE) /*** write site without '/' ***/
					{
					if(j<PALINDROME_MAX_ENZYME-1)
						{
						site[j]=UPPER(EnzymeSite[i]);
						site[++j]=EOS;
						}
					else
						{
						return(FALSE);
						}	
					taille_site=strlen(site);
					}
				if (EnzymeSite[i]=='^')
					pos5_3=i; /*** pos5_3 is localized at the cleavage site  *******/
				i++;
				}
			pos3_5=taille_site-pos5_3; /*** pos3_5 calculation *******/
			if(pos5_3!=-1)
				{
				/* write the new enzyme */
				strcpy(Enzymes[nbr_enzyme].nom,EnzymeName);
				strcpy(Enzymes[nbr_enzyme].site_complet,EnzymeSite);
				strcpy(Enzymes[nbr_enzyme].site,site);
				Enzymes[nbr_enzyme].pos5_3=pos5_3;
				Enzymes[nbr_enzyme].pos3_5=pos3_5;
				Enzymes[nbr_enzyme].taille_site=taille_site;
				Enzymes[nbr_enzyme].palindromic=TRUE;
				Enzymes[nbr_enzyme].Nbr_N5=0;
				Enzymes[nbr_enzyme].Nbr_N3=0;
				Enzymes[nbr_enzyme].select[VECTOR]=Enzyme_Available;
				Enzymes[nbr_enzyme].corp=(short)Enzyme_Available;
				
				
				nbr_enzyme++;
				true_nbr_enzyme++;
				}
			
			}/* end if vara */
		else
			{
			/*** Find symbols  ()/   ***/
			i=0;
			while(EnzymeSite[i]!='\0')
				{
				if (EnzymeSite[i]=='(') var_parenthese_gauche=i;
				if (EnzymeSite[i]==')') var_parenthese_droite=i;
				if (EnzymeSite[i]=='/') var_slash=i;	
				i++;
				}
			/**Detect negative number ***/
			vara=(EnzymeSite[var_parenthese_gauche+1]=='-'?TRUE:FALSE);
			varb=(EnzymeSite[var_slash+1]			 =='-'?TRUE:FALSE);
				
			/*** write site without symbols ***/
			if(var_parenthese_gauche-1<PALINDROME_MAX_ENZYME-1)
				{
				for(i=0;i<=(var_parenthese_gauche-1);i++)
					{
					site[i]=UPPER(EnzymeSite[i]);
					taille_site=i+1;
					site[i+1]='\0';
					}
				taille_site=strlen(site);
				}
			else
				{
				return(FALSE);
				}
			/*** Find pos5_3 and pos3_5 ****/
			pos5_3=0;
			for(i=(var_parenthese_gauche+1+vara);i<=(var_slash-1);i++)
				pos5_3=pos5_3*10+(EnzymeSite[i]-'0');
			pos5_3=(vara==FALSE?pos5_3+taille_site:taille_site-pos5_3); 					/*** Positive number ***/
			
				
			pos3_5=0;
			for(i=(var_slash+1+varb);i<=(var_parenthese_droite-1);i++)
				pos3_5 = pos3_5*10+(EnzymeSite[i]-'0');
			pos3_5=(varb==FALSE?pos3_5+taille_site:taille_site - pos3_5);
			
			/*** fill the complement of restriction site with 'N' ******************/
			vara=0;
			for(i=taille_site;i<MAX(pos3_5,pos5_3);i++)
				vara++;
			/*taille_site = MAX( pos3_5 , MAX( taille_site , pos5_3) );*/
			/* write the new enzyme */
			strcpy(Enzymes[nbr_enzyme].nom,EnzymeName);
			strcpy(Enzymes[nbr_enzyme].site_complet,EnzymeSite);
			strcpy(Enzymes[nbr_enzyme].site,site);
			Enzymes[nbr_enzyme].pos5_3=pos5_3;
			Enzymes[nbr_enzyme].pos3_5=pos3_5;
			Enzymes[nbr_enzyme].taille_site=taille_site;
			Enzymes[nbr_enzyme].palindromic=FALSE;
			Enzymes[nbr_enzyme].Nbr_N5=0;
			Enzymes[nbr_enzyme].Nbr_N3=vara;
			Enzymes[nbr_enzyme].select[VECTOR]=Enzyme_Available;
			Enzymes[nbr_enzyme].corp=(short)Enzyme_Available;

			/* I decided to use this field in order to use Enzyme_Available in Discard_Small procedure */		

			nbr_enzyme++;
			true_nbr_enzyme++;
			/* add an enzyme with the same name BUT with the anti-parallele site */
			Enzymes=(STRUCT_ENZYME*)realloc(Enzymes,(nbr_enzyme+1)*sizeof(STRUCT_ENZYME));
			if(Enzymes==NULL)
				{
				err_printf("OUT OF MEMORY !  (too Many ER ! %d )",nbr_enzyme);
				}
			strcpy(Enzymes[nbr_enzyme].nom,EnzymeName);
			strcpy(Enzymes[nbr_enzyme].site_complet,EnzymeSite);
			strcpy(Enzymes[nbr_enzyme].site,site);
			Enzymes[nbr_enzyme].pos5_3=taille_site-pos3_5+1;
			Enzymes[nbr_enzyme].pos3_5=taille_site-pos5_3+1;
			Enzymes[nbr_enzyme].taille_site=taille_site;
			Enzymes[nbr_enzyme].palindromic=FALSE;
			Enzymes[nbr_enzyme].Nbr_N5=vara;
			Enzymes[nbr_enzyme].Nbr_N3=0;
			Enzymes[nbr_enzyme].select[VECTOR]=Enzyme_Available;
			Enzymes[nbr_enzyme].corp=(short)Enzyme_Available;


			for(i=0;i<strlen(Enzymes[nbr_enzyme].site);i++)
				Enzymes[nbr_enzyme].site[i]=Fct_Complementaire(Enzymes[nbr_enzyme-1].site[strlen(Enzymes[nbr_enzyme-1].site)-1-i]);
			Enzymes[nbr_enzyme].site[taille_site]='\0';
			nbr_enzyme++;
			/* true_nbr_enzyme++;*/
			}
		}/* end else */
	return(nbr_enzyme);
}


/************************************************************************************/

void Cmd_Get_Rebase(void)
{
/************************************************************************************
This procedure read the Rebase File and create an array with all the enzymes and their
caracteristics.
____________________________________________________________________________________*/

FILE *Buffer;
char 	mot1[NOM_MAX_ENZYME],mot2[SITE_MAX_ENZYME],mot3[MAX_CORP];
/*int Enzyme_Available=0;*/



true_nbr_enzyme=0;
nbr_enzyme=0;


	if ((Buffer = fopen(FICHIER_ENZYME, "r" ))==NULL)
		{
		err_printf("Server can't find the REBASE file \"%s\"",FICHIER_ENZYME);
		}
	else
		{
		while(freadRebase(Buffer,mot1,mot2,mot3)==TRUE)
			{
			if(strlen(mot3)>0)
				{
				Make_New_enzyme(mot1, mot2,(int)strlen(mot3));
				}
			}/*end while*/
		fclose(Buffer);

		/* now, discard isoschizomeres of enzyme i, note that procedure
			Cmd_Get_Rebase stored the number of companies selling the
			enzymes in enzyme field Enzymes[i].select[VECTOR] */
		Discard_iso();

		if(nbr_enzyme==0)
			err_printf("No enzyme was found in %s",FICHIER_ENZYME);

		}/* end if file */
Search_Done=FALSE;
}/* end void */



int Discard_iso(void)
	{
	int i,j,vara;
	
	
	for(i=0;i<nbr_enzyme;i++)
		{
		vara=TRUE;
		for(j=0;j<nbr_enzyme;j++)
			{
			if(j==i) continue;
			if( strcmp(Enzymes[i].site_complet,Enzymes[j].site_complet)==ARE_IDENTIC &&
				strcmp(Enzymes[i].site,Enzymes[j].site)==ARE_IDENTIC &&
				Enzymes[i].corp<=Enzymes[j].corp)
				{
					{
					vara=FALSE;
					break;
					}
				}
			}
		if(vara==FALSE) /* if an enzyme was discarded then move the array */
			{
			if(Enzymes[i].palindromic!=TRUE && i+1<nbr_enzyme)
				{
				if(strcmp(Enzymes[i].site_complet,Enzymes[i+1].site_complet)==ARE_IDENTIC)
				true_nbr_enzyme++;/* else true_number will be sub-estimated */
				}
			for(j=i;j<nbr_enzyme-1;j++)
				{
				memcpy(&Enzymes[j],&Enzymes[j+1],sizeof(STRUCT_ENZYME));
				}
			i--;
			true_nbr_enzyme--;
			nbr_enzyme--; /* deleted last enzyme and realloc memory */
			Enzymes=(STRUCT_ENZYME*)realloc(Enzymes,(nbr_enzyme)*sizeof(STRUCT_ENZYME));
			}
		}
	return(nbr_enzyme);
	}



/************************************************************************************/
void ShowSeq(short mode,FILE *out,STRUCT_SITE site_5, STRUCT_SITE site_3)
	{/******************************************************************************
	Display the sequence from site_5 to site_3
	where the cuting site of an enzyme will be found a sign '/' will be printed
	if the sequence position is in FRAME a sign '.' will be printed
	if the sequence position is in FRAME the translation codon will be printed
	******************************************************************************/
	register int i,N_Seq;
	/*****************************************************************************/
	N_Seq = site_5.NumSeq; /* number of the sequence VECTOR or INSERT */
	/* if the two sites are too FAR from each other, write the sequence in two times */
	if(site_5.Loc +PALINDROME_MAX_ENZYME< site_3.Loc - 4)
		{
		/* write first part */
		fprintf(out,"\n  5'  --");
		
		fprintf(out,"<FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));
		
		for(i= site_5.Loc - 4; i<= site_5.Loc +PALINDROME_MAX_ENZYME;i++)
			{
			if (Fct_Frame(N_Seq,i)==IS_IN_FRAME) fprintf(out,".");
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos5_3) {fprintf(out,"/");
			fprintf(out,"</FONT><FONT COLOR=\"#%s\">",(N_Seq==VECTOR?RED_COLOR:BLUE_COLOR));}
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos3_5) fprintf(out," ");
			fprintf(out,"%c",seq[N_Seq].sequence[Fct_Pos(N_Seq,i)]);
			}
		/* write second part */
		fprintf(out,"--  --");
		for(i= site_3.Loc - 4; i<= site_3.Loc + PALINDROME_MAX_ENZYME;i++)
			{
			if (Fct_Frame(N_Seq,i)==IS_IN_FRAME) fprintf(out,".");
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos5_3) {fprintf(out,"/");
			fprintf(out,"</FONT><FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));}
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos3_5) fprintf(out," ");
			fprintf(out,"%c",seq[N_Seq].sequence[Fct_Pos(N_Seq,i)]);
			}
		
		fprintf(out,"</FONT>");
		
		/* write anti-first part */
		fprintf(out,"--  3'\n  3'  --");
		
		fprintf(out,"<FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));
		
		for(i= site_5.Loc - 4; i<= site_5.Loc +PALINDROME_MAX_ENZYME;i++)
			{
			if (Fct_Frame(N_Seq,i)==IS_IN_FRAME) fprintf(out,".");
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos3_5) {fprintf(out,"/");
			fprintf(out,"</FONT><FONT COLOR=\"#%s\">",(N_Seq==VECTOR?RED_COLOR:BLUE_COLOR));}
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos5_3) fprintf(out," ");
			fprintf(out,"%c",Fct_Complementaire(seq[N_Seq].sequence[Fct_Pos(N_Seq,i)]));
			}
		/* write anti-second part */
		fprintf(out,"--  --");
		for(i= site_3.Loc - 4; i<= site_3.Loc + PALINDROME_MAX_ENZYME;i++)
			{
			if (Fct_Frame(N_Seq,i)==IS_IN_FRAME) fprintf(out,".");
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos3_5) {fprintf(out,"/");
			fprintf(out,"</FONT><FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));}
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos5_3) fprintf(out," ");
			fprintf(out,"%c",Fct_Complementaire(seq[N_Seq].sequence[Fct_Pos(N_Seq,i)]));
			}
			
		fprintf(out,"</FONT>");
			
		fprintf(out,"--  5'\n");
		/* write the translation */
		
		

		
		if(seq[INSERT].pos_ATG!=FALSE || seq[VECTOR].pos_ATG!=FALSE)
			{
			/* write first part */
			fprintf(out,"  NH2    ");
			fprintf(out,"<FONT COLOR=\"#%s\">",YELLOW_COLOR);
			for(i= site_5.Loc - 4; i<= site_5.Loc +PALINDROME_MAX_ENZYME;i++)
				{
				if(Fct_Frame(N_Seq,i)==IS_IN_FRAME)
					fprintf(out,"%c ",Fct_Traduction(seq[N_Seq].sequence[Fct_Pos(N_Seq,i)]
												,seq[N_Seq].sequence[Fct_Pos(N_Seq,i+1)]
												,seq[N_Seq].sequence[Fct_Pos(N_Seq,i+2)]));
				else
					fprintf(out," ");
				if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos5_3) fprintf(out," ");
				if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos3_5) fprintf(out," ");

				}
			/* write second part */
			fprintf(out,"--  --");
			for(i= site_3.Loc - 4; i<= site_3.Loc + PALINDROME_MAX_ENZYME;i++)
				{
				if(Fct_Frame(N_Seq,i)==IS_IN_FRAME)
					fprintf(out,"%c ",Fct_Traduction(seq[N_Seq].sequence[Fct_Pos(N_Seq,i)]
											,seq[N_Seq].sequence[Fct_Pos(N_Seq,i+1)]
											,seq[N_Seq].sequence[Fct_Pos(N_Seq,i+2)]));
				else
					fprintf(out," ");
				if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos5_3) fprintf(out," ");
				if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos3_5) fprintf(out," ");
				}
			fprintf(out,"</FONT>");
			fprintf(out,".COOH\n");
			}
		}
	else /* if the two sites are close from each other, write the sequence in one time */
		{
		/* write direct strand */
		fprintf(out,"\n  5'  --");
		
		fprintf(out,"<FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));

		
		for(i= site_5.Loc - 4; i<= site_3.Loc +PALINDROME_MAX_ENZYME;i++)
			{
			if (Fct_Frame(N_Seq,i)==IS_IN_FRAME) fprintf(out,".");
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos5_3) {fprintf(out,"/");fprintf(out,"</FONT><FONT COLOR=\"#%s\">",(N_Seq==VECTOR?RED_COLOR:BLUE_COLOR));}
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos3_5) fprintf(out," ");
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos5_3) {fprintf(out,"/");
			fprintf(out,"</FONT><FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));}
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos3_5) fprintf(out," ");
			fprintf(out,"%c",seq[N_Seq].sequence[Fct_Pos(N_Seq,i)]);
			}
		/* write anti-strand */
		fprintf(out,"</FONT>");
		fprintf(out,"--  3'\n  3'  --");
		fprintf(out,"<FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));
		for(i= site_5.Loc - 4; i<= site_3.Loc +PALINDROME_MAX_ENZYME;i++)
			{
			if (Fct_Frame(N_Seq,i)==IS_IN_FRAME) fprintf(out,".");
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos5_3) fprintf(out," ");
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos3_5) {fprintf(out,"/");
			fprintf(out,"</FONT><FONT COLOR=\"#%s\">",(N_Seq==VECTOR?RED_COLOR:BLUE_COLOR));}
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos5_3) fprintf(out," ");
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos3_5) {fprintf(out,"/");
			fprintf(out,"</FONT><FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));}
			fprintf(out,"%c",Fct_Complementaire(seq[N_Seq].sequence[Fct_Pos(N_Seq,i)]));
			}
		fprintf(out,"</FONT>");
		fprintf(out,"--  5'\n");
		
		/* write the translation */
		if(seq[INSERT].pos_ATG!=FALSE || seq[VECTOR].pos_ATG!=FALSE)
			{
			fprintf(out,"  NH2    ");
			fprintf(out,"<FONT COLOR=\"#%s\">",YELLOW_COLOR);
			for(i= site_5.Loc - 4; i<= site_3.Loc +PALINDROME_MAX_ENZYME;i++)
				{
				if(Fct_Frame(N_Seq,i)==IS_IN_FRAME)
					fprintf(out,"%c ",Fct_Traduction(seq[N_Seq].sequence[Fct_Pos(N_Seq,i)]
												,seq[N_Seq].sequence[Fct_Pos(N_Seq,i+1)]
												,seq[N_Seq].sequence[Fct_Pos(N_Seq,i+2)]));
				else
					fprintf(out," ");
				if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos5_3) fprintf(out," ");
				if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos3_5) fprintf(out," ");
				if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos5_3) fprintf(out," ");
				if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos3_5) fprintf(out," ");
					
				}
			fprintf(out,"</FONT>");
			fprintf(out,".COOH\n");
			}
		}
	}



/************************************************************************************/




/*** Add a Site in memory see : Cmd_Get_Site  **/
int Add_Site(int k,int _NumSeq,int _Loc)
	{
	
	Sites=(STRUCT_SITE*)realloc(Sites,(nbr_sites+1)*sizeof(STRUCT_SITE));
	if(Sites==NULL)
		{
		err_printf("SERVER OUT OF MEMORY (TOO MUCH SITES %d)",nbr_sites);
		Sites=(STRUCT_SITE*)realloc(Sites,(1)*sizeof(STRUCT_SITE));
		return(FALSE);
		}
	else
		{
		Sites[nbr_sites].NumEnz=k;
		Sites[nbr_sites].NumSeq=_NumSeq;
		Sites[nbr_sites].Loc=_Loc;
		nbr_sites++;
		}
	return(_Loc);
	}


void Cmd_Get_Site(void)
{
/***********************************
	Get all sites that are in VECTOR
	and in INSERT
************************************/
int j,vara=0,NumSeq,Nbr_N5=0,start;
register int i,k;
char *result;
char *p_seq;
/**************************************/
nbr_sites=0;
/* scan VECTOR and INSERT */
for(NumSeq=VECTOR;NumSeq<=INSERT;NumSeq++)
	{
	if(seq[NumSeq].npb>0)
	/* scan all enzymes */
	for(k=0;k<nbr_enzyme;k++)
		{
		/* if the enzyme site does not contain degenerate base, strpbrk function can be used */
		if(strpbrk(Enzymes[k].site,"YRMKSWBDHVN")==NULL) /* Does not contain degenerate bases */
			{
			p_seq=&seq[NumSeq].sequence[1];
			while((result=strstr(p_seq, Enzymes[k].site ))!=NULL)
				{
				/* if the site has been found, add a new site */
				i=seq[NumSeq].npb-strlen(result)+1;
				if(Add_Site(k,NumSeq,seq[NumSeq].npb-strlen(result)+1)==FALSE)
						{NumSeq=INSERT+1;nbr_enzyme=0;break;} /* abort if problem */
				p_seq=result+1;
				}
			start=seq[NumSeq].npb-Enzymes[k].taille_site+2;
			}
		else
			{
			start=1; /* else scan from the begin of the sequence */
			}

		Nbr_N5 = Enzymes[k].Nbr_N5; /* there is no use to scan the 'N' at the begining of the site see: Get_Rebase */
		for(i=start;i<=seq[NumSeq].npb;i++) /*** Scanning sequence *******/ 
			{
			/* site not found by default */
			vara=FALSE;
			for(j=0;j<Enzymes[k].taille_site;j++)
				if(Fct_Identique(seq[NumSeq].sequence[Fct_Pos(NumSeq,i+j)],Enzymes[k].site[j])!=1)
					{vara=TRUE;break;}
			/* if site was found, add a site in memory */
			if (vara==FALSE)
				{
				if(Add_Site(k,NumSeq,i-Nbr_N5)==FALSE)
					{NumSeq=INSERT+1;i=seq[NumSeq].npb+1;k=nbr_enzyme;NumSeq=3;nbr_enzyme=0;break;}/* abort if problem */
				}
			}
		}
	}
Search_Done=TRUE;
}


  /******************************************************/
 /* display a linar  restriction map with the 6 frames */
/******************************************************/


int ResMap(FILE *out,short mode,int var_seq)
	{
	/* there was a bug in this function till those functions where set as internal ???!!*/ 
	
	int var_h[LARGEUR_ECRAN+LARGEUR_ECRAN]; /* define a line on screen if TRUE: there is a column to fill*/
	char line[LARGEUR_ECRAN+LARGEUR_ECRAN+LARGEUR_ECRAN];
	STRUCT_ENZXY Topaze; /* define all sites that are present on a DNA sequence segment */
	/* i is the beginning of the sequence fragment displayed on screen */
	int j,k,m,p,q,nbr_line=0,vara=0,var_len_screen=0,linelen=0;
	int i,n;
	int var_begin=0,var_end=0; /* begin and end of the sequence description */
	div_t	r;
	/**********************************************************************************/
	var_len_screen=LARGEUR_ECRAN-NOM_MAX_ENZYME; /* length of a on the screen */
	line[0]=EOS;linelen=0;
	if(seq[var_seq].npb==0) return(FALSE);
	
	
	var_begin = seq[var_seq].var_min;
	var_end   = seq[var_seq].var_max;
	
	printf("<H3>%s [%d-%d]</H3><BR>",seq[var_seq].FICHIER_ADN,var_begin,var_end);
		
	var_nbr_xy=0; /*no site defined on fragment */
	
	/* here, I use the field: 'select' to tag the UNIQUE sites in var_seq*/
	for(i=0;i<nbr_enzyme;i++)
		{
		Enzymes[i].select[var_seq]=FALSE; /* this enzyme is not unique by default */
		if(Enzymes[i].palindromic==TRUE)
			{
			if(Fct_N_sites(i,var_seq,1,seq[var_seq].npb,TRUE)==1)
				Enzymes[i].select[var_seq]=TRUE; /* this enzyme is unique */
			}
		else if(i<nbr_enzyme-1)
			if(strcmp(Enzymes[i].site_complet,Enzymes[i+1].site_complet)==ARE_IDENTIC)
			{
			Enzymes[i+1].select[var_seq]=FALSE;
			if(	Fct_N_sites(i,var_seq,1,seq[var_seq].npb,TRUE) + Fct_N_sites(i+1,var_seq,1,seq[var_seq].npb,TRUE)==1 )
				{
				Enzymes[i].select[var_seq]=TRUE; /* this enzyme is unique */
				Enzymes[i+1].select[var_seq]=TRUE; /* this enzyme is unique */
				}
			i++;
			}
		}
	/**********************************************************************************/

/* scan all the sequence fragment */
for(i=var_begin;i<= var_end;i=i+var_len_screen)
	{
	/* find all sites that are present on this portion of DNA */
	var_nbr_xy=0;
	nbr_line=0;

	for(n=0;n<nbr_sites;n++)
		{
		/* is this site is on var_seq ?*/
		if(Sites[n].NumSeq!=var_seq)
			continue;
		/* is this site present on this portion ? */
		if((Sites[n].Loc >= i) && (Sites[n].Loc< i+var_len_screen) )
			{
			/* memorise site */
			if((EnzXY=(STRUCT_ENZXY*)realloc(EnzXY,(var_nbr_xy+1)*sizeof(STRUCT_ENZXY)))==NULL)
				{
				printf("Memory Full ! var_nbr_xy=%d et n=%d\n",var_nbr_xy,n);
				printf("*** ERROR ***\n");
				INKEY;BEEP;
				break;
				}
			else
				{
				EnzXY[var_nbr_xy].LocX= (Sites[n].Loc)-i;
				EnzXY[var_nbr_xy].LocY= 0;
				EnzXY[var_nbr_xy].NumSite=n;
				var_nbr_xy++;
				}
			}
		}

	/**********************************************************************/
	/* sort sites from  5' to 3' */
	vara=TRUE;
	while(vara==TRUE)
		{
		vara=FALSE;
		for(p=0;p<var_nbr_xy-1;p++)
			if(EnzXY[p].LocX>EnzXY[p+1].LocX)
				{
				vara=TRUE;
				memcpy(&Topaze,&EnzXY[p],sizeof(STRUCT_ENZXY));
				memcpy(&EnzXY[p],&EnzXY[p+1],sizeof(STRUCT_ENZXY));
				memcpy(&EnzXY[p+1],&Topaze,sizeof(STRUCT_ENZXY));
				}
		}
	/* if  a site displayed on the line n°:LocY will be write OVER a second site, then add 1 to Loc Y , add the number of line*/
	vara=TRUE;
	while(vara==TRUE)
		{
		vara=FALSE;
		for(m=0;m<var_nbr_xy-1;m++)
			for(p=m+1;p<var_nbr_xy;p++)
		 		{
		 		if(EnzXY[m].LocX+2+strlen(Enzymes[Sites[EnzXY[m].NumSite].NumEnz].nom)>=EnzXY[p].LocX && 
		 		   EnzXY[m].LocY == EnzXY[p].LocY)
		 			{
		 			EnzXY[p].LocY=EnzXY[m].LocY+1;
		 			vara=TRUE;
		 			}
		 		if(EnzXY[p].LocY>nbr_line)
		 			nbr_line=EnzXY[p].LocY;
		 		}
		 }
	/* init var_h to FALSE */
	 for(m=0;m<(LARGEUR_ECRAN+LARGEUR_ECRAN);m++)
	 	{var_h[m]=FALSE;}
	/* display all lines with the name of enzymes */
	 for(m=nbr_line;m>=0;m--)
	 	{
	 	linelen=0;
	 	/* scan all column in the line */
	 	for(n=0;n<var_len_screen;n++)
	 		{
	 		vara=FALSE;
	 		/* scan all the sites of EnzXY */
		 	for(p=0;p<var_nbr_xy;p++)
		 		{
		 		/* if this site is localised in this column */
		 		if(EnzXY[p].LocY==m && EnzXY[p].LocX==n)
		 			{
		 			vara=TRUE;
		 			var_h[n]=TRUE;
		 			/* if enzyme is unique displays it in UPPERCASE */
		 			if(Enzymes[ Sites[EnzXY[p].NumSite].NumEnz ].select[var_seq]==TRUE)
		 				{
		 				for(q=0;q<strlen(Enzymes[ Sites[EnzXY[p].NumSite].NumEnz ].nom);q++)
		 					/*fprintf(out,"%c",UPPER(Enzymes[ Sites[EnzXY[p].NumSite].NumEnz ].nom[q]));*/
		 					line[linelen++]=UPPER(Enzymes[ Sites[EnzXY[p].NumSite].NumEnz ].nom[q]);
		 				}
		 			else /* displays it in LOWERCASE */
		 				{
		 				for(q=0;q<strlen(Enzymes[ Sites[EnzXY[p].NumSite].NumEnz ].nom);q++)
		 					/*fprintf(out,"%c",LOWER(Enzymes[ Sites[EnzXY[p].NumSite].NumEnz ].nom[q]));*/
		 					line[linelen++]=LOWER(Enzymes[ Sites[EnzXY[p].NumSite].NumEnz ].nom[q]);
		 				}
		 			/* add the length of the word that as been displayed to n */
		 			n = n + strlen(Enzymes[Sites[EnzXY[p].NumSite].NumEnz].nom)-1;
		 			}
		 		}
		 	if(vara==FALSE)
		 	 	{
		 	 	/* if there was nothing in the upper line print ' ' else ':' */
		 		/*fprintf(out,"%c",((var_h[n]==TRUE)?':':' '));*/
		 		line[linelen++]=((var_h[n]==TRUE)?':':' ');
		 		}
		 	}
	 	line[linelen]=EOS;
	 	fprintf(out,"%s\n",line);
	 	}
	 linelen=0;
	 /* display one line */
	 for(n=0;n<var_len_screen;n++)
		 	/*fprintf(out,"%c",((var_h[n]==TRUE)?':':' '));*/
		 	line[linelen++]=((var_h[n]==TRUE)?':':' ');
	line[linelen]=EOS;
	 	fprintf(out,"%s\n",line);
	 linelen=0;
	 /**********************************************************************/
	 /* display one strand */
	for(j=i;(j< i+var_len_screen && j<= seq[var_seq].npb);j++)
	 	/*fprintf(out,"%c",seq[var_seq].sequence[Fct_Pos(var_seq,j)]);*/
	 	line[linelen++]=seq[var_seq].sequence[Fct_Pos(var_seq,j)];
	line[linelen]=EOS;
	 	fprintf(out,"%s",line);
	fprintf(out," \\ %d\n",Fct_Pos(var_seq,i));
	/* display € each 10 bases */
	linelen=0;
	for(j=i;(j< i+var_len_screen && j<= seq[var_seq].npb);j++)
		{
		r = div(j,10);
		/*fprintf(out,"%c",(((int)r.rem==0)?'€':' '));*/
		line[linelen++]=(((int)r.rem==0)?'€':' ');
		}
	line[linelen]=EOS;
	 fprintf(out,"%s  \\\n",line);
	/* display anti strand */
	linelen=0;
	for(j=i;(j< i+var_len_screen && j<= seq[var_seq].npb);j++)
		/*fprintf(out,"%c",Fct_Complementaire(seq[var_seq].sequence[Fct_Pos(var_seq,j)]));*/
		line[linelen++]=Fct_Complementaire(seq[var_seq].sequence[Fct_Pos(var_seq,j)]);
	line[linelen]=EOS;
	fprintf(out,"%s   \\ %d\n",line,Fct_Pos(var_seq,i+var_len_screen-1));
	/* for each 3 frames 5'->3' show the resulting sequence ******************/
	
	for(k=0;k<=2;k++)
	 	{
	 	linelen=0;
	 	for(j=i+k;(j< i+var_len_screen+k && j<= seq[var_seq].npb);j++)
	 		{
	 		if(Fct_Frame(var_seq,Fct_Pos(var_seq,k+j))==0)
	 			line[linelen++]=Translation_at(var_seq,j);
	 			/*fprintf(out,"%c",Fct_Traduction( seq[var_seq].sequence[Fct_Pos(var_seq,j)],
	 										seq[var_seq].sequence[Fct_Pos(var_seq,j+1)],
	 										seq[var_seq].sequence[Fct_Pos(var_seq,j+2)]));*/
	 		else
	 			line[linelen++]=((var_h[j-i-k]==TRUE)?':':' ');
	 			/*fprintf(out,"%c",((var_h[j-i-k]==TRUE)?':':' '));*/
	 		}
	 	line[linelen]=EOS;
	 	fprintf(out,"%s ->\n",line);
	 	}
	 /* for each 3 frames 3'->5' show the resulting sequence ******************/
	for(k=2;k>=0;k--)
	 	{
	 	linelen=0;
	 	for(j=i+k;(j< i+var_len_screen+k && j<= seq[var_seq].npb);j++)
	 		{
	 		if(Fct_Frame(var_seq,Fct_Pos(var_seq,k+j))==0)
	 			line[linelen++]=Fct_Traduction( seq[var_seq].sequence[Fct_Pos(var_seq,j)],
	 										seq[var_seq].sequence[Fct_Pos(var_seq,j-1)],
	 										seq[var_seq].sequence[Fct_Pos(var_seq,j-2)]);
	 			/*fprintf(out,"%c",Fct_Traduction( seq[var_seq].sequence[Fct_Pos(var_seq,j)],
	 										seq[var_seq].sequence[Fct_Pos(var_seq,j-1)],
	 										seq[var_seq].sequence[Fct_Pos(var_seq,j-2)]));*/
	 		else
	 			line[linelen++]=((var_h[j-i-k]==TRUE)?':':' ');
	 			/*fprintf(out,"%c",((var_h[j-i-k]==TRUE)?':':' '));*/
	 		}
	 	line[linelen]=EOS;
	 	fprintf(out,"%s <-\n",line);
	 	}
	 
	/**********************************************************************/
	 /* display one line */
	 linelen=0;
	for(n=0;n<var_len_screen;n++)
		 /*fprintf(out,"%c",((var_h[n]==TRUE)?':':' '));*/
		 line[linelen++]=((var_h[n]==TRUE)?':':' ');
	line[linelen]=EOS;
	fprintf(out,"%s\n",line);
	/* do the same thing as previously, but this time, displays the loc */
	for(m=0;m<=nbr_line;m++)
	 	{
	 	for(n=0;n<var_len_screen;n++)
	 		{
	 		vara=FALSE;
		 	for(p=0;p<var_nbr_xy;p++)
		 		{
		 		if(EnzXY[p].LocY==m && EnzXY[p].LocX==n)
		 			{
		 			vara=TRUE;
		 			var_h[n]=FALSE; /* remove the information: there was something in the column */
		 			fprintf(out,"%-5d",EnzXY[p].LocX+i); /* number displayed over 5 column */
		 			n = n + 5 -1;
		 			}
		 		}
		 	if(vara==FALSE)
		 	 	{
		 		fprintf(out,"%c",((var_h[n]==TRUE)?':':' '));
		 		}
		 	}
	 	fprintf(out,"\n");
	 	}
	 fprintf(out,"\n\n");
	 }
	return(TRUE);
	}
/************************************************************************/

  /************************************************************************/
 /* compare the number of site in INSERT/VECTOR and in the cloning boxes */
/************************************************************************/
bool_t Intersections(FILE *out, short mode)
	{
	int IN_I,OUT_I,IN_V,OUT_V;
	register int i;
	
	if(Search_Done==FALSE) Cmd_Get_Site();
	
	printf("<HR><CENTER><H1>INTERSECTIONS</H1></CENTER><BR>");
	printf("<B>IN</B>:inside the region delimited by the external cloning boxes.<BR>");
	printf("<B>OUT</B>:outside the the region delimited by the external cloning boxes.<BR>");
	fprintf(out,"<TABLE BORDER=2 CELLSPACING=2 CELLPADDING=2  ALIGN=\"CENTER\">");
	fprintf(out,"<TR>");
	fprintf(out,"<TD COLSPAN=\"75\" ALIGN=\"LEFT\"></TD>");
	fprintf(out,"<TD COLSPAN=\"30\"><CENTER>VECTOR</CENTER></CENTER>");
	fprintf(out,"</TD>");
	fprintf(out,"<TD COLSPAN=\"30\"><CENTER>INSERT</CENTER></TD>");
	fprintf(out,"</TR>");
	fprintf(out,"<TR>");
	fprintf(out,"<TD COLSPAN=\"75\"></TD>");
	fprintf(out,"<TD COLSPAN=\"15\"><CENTER>IN</CENTER></TD>");
	fprintf(out,"<TD COLSPAN=\"15\"><CENTER>OUT</CENTER></TD>");
	fprintf(out,"<TD COLSPAN=\"15\" ><CENTER>IN</CENTER></TD>");
	fprintf(out,"<TD COLSPAN=\"15\"><CENTER>OUT</CENTER></TD>");
	fprintf(out,"</TR>");
	
	for(i=0;i<nbr_enzyme;i++)
		{
		/* get the number of site in INSERT, IN the cloning box */
		IN_I  = Fct_N_sites(i,INSERT,seq[INSERT].var_min,seq[INSERT].var_max,TRUE);
		/* get the number of site in INSERT, out of the cloning box */
		OUT_I = Fct_N_sites(i,INSERT,seq[INSERT].var_min-1,seq[INSERT].var_max+1,FALSE);
		/* get the number of site in VECTOR, IN the cloning box */
		IN_V  = Fct_N_sites(i,VECTOR,seq[VECTOR].var_min,seq[VECTOR].var_max,TRUE);
		/* get the number of site in VECTOR, out of the cloning box */
		OUT_V = Fct_N_sites(i,VECTOR,seq[VECTOR].var_min-1,seq[VECTOR].var_max+1,FALSE);
		
		
		
		fprintf(out,"<TR>");
		fprintf(out,"<TD COLSPAN=\"75\">");
		printWebEnzyme(Enzymes[i].nom,out);
		fprintf(out," [%s]",Enzymes[i].site_complet);
		fprintf(out,"</TD>");
		if(IN_V==1 && OUT_V==1)
			{
			fprintf(out,"<TD COLSPAN=\"15\" BGCOLOR=\"#00FF00\">%d</TD>",IN_V);
			fprintf(out,"<TD COLSPAN=\"15\" BGCOLOR=\"#00FF00\">%d</TD>",OUT_V);
			}
		else if((IN_V==1 && OUT_V==0) || (IN_V==0 && OUT_V==1))
			{
			fprintf(out,"<TD COLSPAN=\"15\" BGCOLOR=\"#00FFFF\">%d</TD>",IN_V);
			fprintf(out,"<TD COLSPAN=\"15\" BGCOLOR=\"#00FFFF\">%d</TD>",OUT_V);
			}
		else 
			{
			fprintf(out,"<TD COLSPAN=\"15\">%d</TD>",IN_V);
			fprintf(out,"<TD COLSPAN=\"15\">%d</TD>",OUT_V);
			}
			
		if(IN_I==1 && OUT_I==1)
			{
			fprintf(out,"<TD COLSPAN=\"15\" BGCOLOR=\"#00FF00\">%d</TD>",IN_I);
			fprintf(out,"<TD COLSPAN=\"15\" BGCOLOR=\"#00FF00\">%d</TD>",OUT_I);
			}
		else if((IN_I==1 && OUT_I==0) || (IN_I==0 && OUT_I==1))
			{
			fprintf(out,"<TD COLSPAN=\"15\" BGCOLOR=\"#00FFFF\">%d</TD>",IN_I);
			fprintf(out,"<TD COLSPAN=\"15\" BGCOLOR=\"#00FFFF\">%d</TD>",OUT_I);
			}
		else
			{
			fprintf(out,"<TD COLSPAN=\"15\">%d</TD>",IN_I);
			fprintf(out,"<TD COLSPAN=\"15\">%d</TD>",OUT_I);
			}
		fprintf(out,"</TR>");

		
		
		
		if(i<nbr_enzyme-1)
			{
			if(Enzymes[i].palindromic!=TRUE) i++;
			}
		}
	
		fprintf(out,"<TR>");
		fprintf(out,"<TD COLSPAN=\"75\"></TD>");
		fprintf(out,"<TD COLSPAN=\"15\"><CENTER>IN</CENTER></TD>");
		fprintf(out,"<TD COLSPAN=\"15\"><CENTER>OUT</CENTER></TD>");
		fprintf(out,"<TD COLSPAN=\"15\" ><CENTER>IN</CENTER></TD>");
		fprintf(out,"<TD COLSPAN=\"15\"><CENTER>OUT</CENTER></TD>");
		fprintf(out,"</TR>");
		fprintf(out,"<TR>");
		fprintf(out,"<TD COLSPAN=\"75\" ALIGN=\"LEFT\"></TD>");
		fprintf(out,"<TD COLSPAN=\"30\"><CENTER>VECTOR</CENTER></CENTER>");
		fprintf(out,"</TD>");
		fprintf(out,"<TD COLSPAN=\"30\"><CENTER>INSERT</CENTER></TD>");
		fprintf(out,"</TR>");
		fprintf(out,"</TABLE>");
	return(TRUE);
	}







void Orient_Insert_After_CIP(STRUCT_STRATEGY *Strategy,short mode,FILE *out)
	{
	/*****************************************************************************
	this function find a enzyme that cut one time in insert and one time in vector 
		to direct insert in dephoso-phorylated vector
	******************************************************************************/
	int vara=FALSE,i;
	/******************************************************************************/
	if(Strategy->use_CIP==TRUE)
		fprintf(out,"Enzyme(s) that could be used to direct insert:");
	else
		fprintf(out,"Enzyme(s) that could be used to verify your construction:<BR><UL>");
	for(i=0;i<nbr_enzyme;i++)
		{
		/* discard this enzyme if it is used in the strategy */
		if( (i == Sites[Strategy->Couple[SIDE_5].NumSite[VECTOR]].NumEnz) ||
			(i == Sites[Strategy->Couple[SIDE_3].NumSite[VECTOR]].NumEnz) ||
			(i == Sites[Strategy->Couple[SIDE_5].NumSite[INSERT]].NumEnz) ||
			(i == Sites[Strategy->Couple[SIDE_3].NumSite[INSERT]].NumEnz) )
				continue;
		
		/* one time out of vector  and one time in insert */
		if( Fct_N_sites(i,VECTOR,
			Sites[Strategy->Couple[SIDE_5].NumSite[VECTOR]].Loc-1,
			Sites[Strategy->Couple[SIDE_3].NumSite[VECTOR]].Loc+1,
			FALSE)==1 &&
			Fct_N_sites(i,INSERT,
			Sites[Strategy->Couple[SIDE_5].NumSite[INSERT]].Loc+1,
			Sites[Strategy->Couple[SIDE_3].NumSite[INSERT]].Loc-1,
			TRUE)==1)
				{
				vara=TRUE;
				fprintf(out,"<LI>");
				printEnzymeName(i);
				}
		/* if enzyme is not palindromic do not use the next one that will be the same */
		if(Enzymes[i].palindromic!=TRUE) i++;
		}
	if(vara==FALSE)
		{
		fprintf(out,"<LI><I>no enzyme was found.</I>");
		}
	fprintf(out,"</UL><BR>\n");
	}



/* draw a line **********************************/
void DRAW_HR(FILE *out,char car)
	{
	int i;
	for(i=0;i<LARGEUR_ECRAN;i++) fprintf(out,"%c",car);
	fprintf(out,"<BR>\n");
	}
/***********************************/
int Sign(int vara)
	{
	return( (abs(vara)==vara) ? 1 : -1 );
	}
/* return the frame position ******/
int Fct_Frame(int NumSeq,int nombre)
	{
	/***********************
	This function returns a number (0,1 or 2)
	that is the remainder of nombre/3
	by default the ATG localisation return 0
	
	123	456	789	012	345	678 <- localisation
	ATG CAT GCA TGC ATG CAT
	012	012	012	012	012	012 <- frame
	***********************/
	div_t	r,s;
	int i,j;

	r = div(nombre,3);
	if(seq[NumSeq].pos_ATG==FALSE)
		{
		return(r.rem);
		}
	else
		{
		s = div(seq[NumSeq].pos_ATG,3);
		i=(int)r.rem; /* remainder of the localisation of nombre */
		j=(int)s.rem; /* remainder of the localisation of ATG    */
		
		switch(j)
			{
			case(0):return(i);break;
			case(1):switch(i)
				{
				case(0):return(2);break;
				case(1):return(0);break;
				case(2):return(1);break;
				}break;
			case(2):switch(i)
				{
				case(0):return(1);break;
				case(1):return(2);break;
				case(2):return(0);break;
				}break;
			}
		}
	return(FALSE);
	}
/************************************ Fct_Pos ***/	
int	Fct_Pos(int NumSeq,int position)
/*** I use this function to find restriction site that are between the end and the beginning of the cirular sequence...**/
	{
	while( position>seq[NumSeq].npb)	position=  position-seq[NumSeq].npb;
	while( position<1)				position=1+position+seq[NumSeq].npb;
	return( position);
	}



void write_date(FILE *out)
	{
	time_t now;
	char *s;
	now = time( NULL );
	s = ctime( &now );
	fprintf(out,"%s", s );
	}



/***returns enzyme type *****/
int Fct_type(STRUCT_ENZYME Enzyme)
	{
	if(Enzyme.pos3_5==Enzyme.pos5_3) return(TYPE_BLUNT);
	else if(Enzyme.pos3_5<Enzyme.pos5_3) return(TYPE_3_OVER);
	else return(TYPE_5_OVER);
	}

/****returns if a char is coding a nucleotide  ***/
int Fct_est_ADN(char lettre)
	{
	int i=0;
	switch(LOWER(lettre))
		{
		case 'a':
		case 't':
		case 'g':
		case 'c': i=TRUE;break;
		
		case 'y':
		case 'r':
		case 'm':
		case 'k':
		case 's':
		case 'w': i=AMBIGOUS;break;
		case 'b':
		case 'd':
		case 'h':
		case 'v': i=AMBIGOUS+1;break;
		case 'n': i=AMBIGOUS+2;break;
		default: i=FALSE;break;
		}
	return(i);
	}
/****** returns complementary nucleotide  ***/
char Fct_Complementaire(char lettre)
	{
	char mot;
	bool_t var_font=FALSE;
	
	if(isupper(lettre)!=0) var_font=TRUE;
	lettre=LOWER(lettre);
	switch(lettre)
		{
		case 'a': mot='t';break;
		case 't': mot='a';break;
		case 'g': mot='c';break;
		case 'c': mot='g';break;		
		case 'n': mot='n';break;
				
		case 'y': mot='r';break;
		case 'r': mot='y';break;
		
		case 'm': mot='k';break;
		case 'k': mot='m';break;
		case 's': mot='s';break;
		case 'w': mot='w';break;
		
		case 'b': mot='v';break;
		case 'd': mot='h';break;
		case 'h': mot='d';break;
		case 'v': mot='b';break;
		
		
		default: mot='?';break;
		}
		if (var_font==TRUE) mot=UPPER(mot);
	return(mot);
	}
/**** returns if lettre1 is compatible with lettre2 ***/
int Fct_Identique(char lettre1,char lettre2)
	{
	int retour=0;
	lettre2=LOWER(lettre2);
	
	switch(LOWER(lettre1))
		{
		case 'a':
				if(lettre2=='a')
					retour=TRUE;
				else if(lettre2=='n' || lettre2=='d'|| lettre2=='h'|| lettre2=='v'
								|| lettre2=='m'|| lettre2=='w'|| lettre2=='r')
					retour=TRUE;
				break;
		case 't':
				if(lettre2=='t')
					retour=TRUE;
				else if(lettre2=='n' || lettre2=='y'|| lettre2=='b'|| lettre2=='d'
								|| lettre2=='h'|| lettre2=='k'|| lettre2=='w')
					retour=TRUE;
				break;
		case 'g':
				if(lettre2=='g')
					retour=TRUE;
				else if(lettre2=='n' || lettre2=='r'|| lettre2=='b'|| lettre2=='d'
								|| lettre2=='v'|| lettre2=='k'|| lettre2=='s')
					retour=TRUE;
				break;
		case 'c':
				if(lettre2=='c')
					retour=TRUE;
				else if(lettre2=='n' || lettre2=='y'|| lettre2=='b'|| lettre2=='h'
								|| lettre2=='v'|| lettre2=='m'|| lettre2=='s')
					retour=TRUE;
				break;
		case 'n':
				retour=TRUE;break;
		case 'y':
				if(lettre2=='n' || lettre2=='y' || lettre2=='c' || lettre2=='t')
					retour=TRUE;break;
		case 'r':
				if(lettre2=='n' || lettre2=='r' || lettre2=='a' || lettre2=='g')
					retour=TRUE;break;
		case 'm':
				if(lettre2=='n' || lettre2=='m' || lettre2=='a' || lettre2=='c')
					retour=TRUE;break;
		case 'k':
				if(lettre2=='n' || lettre2=='k' || lettre2=='g' || lettre2=='t')
					retour=TRUE;break;
		case 's':
				if(lettre2=='n' || lettre2=='s' || lettre2=='g' || lettre2=='c')
					retour=TRUE;break;
		case 'w':
				if(lettre2=='n' || lettre2=='w' || lettre2=='a' || lettre2=='t')
					retour=TRUE;break;
		case 'b':
				if(lettre2=='n' || lettre2=='b' || lettre2=='c' || lettre2=='g' || lettre2=='t')
					retour=TRUE;break;
		case 'd':
				if(lettre2=='n' || lettre2=='d' || lettre2=='a' || lettre2=='g' || lettre2=='t')
					retour=TRUE;break;
		case 'h':
				if(lettre2=='n' || lettre2=='h' || lettre2=='a' || lettre2=='c' || lettre2=='t')
					retour=TRUE;break;
		case 'v':
				if(lettre2=='n' || lettre2=='v' || lettre2=='a' || lettre2=='c' || lettre2=='g')
					retour=TRUE;break;
		default: retour=FALSE;break;
		
		}
	return(retour);
	}
/*****Gives the anti-parallele sequence **************************/
void Fct_Reverse(char *mot)
	{
	nbp_t len=strlen(mot);
	AntiParallele(mot,len);
	}

/***** returns if AminoAcid Is a stop Codon **************************/
int Is_Codon_Stop(char AminoAcid)
	{
	if(AminoAcid==VAR_OCHRE || AminoAcid==VAR_AMBRE || AminoAcid==VAR_OPALE)
		return(TRUE);
	else
		return(FALSE);
	}

/***** returns the translation at position pos ****/
char Translation_at(short NumSeq,int pos)
	{
	return(Fct_Traduction(
		seq[NumSeq].sequence[Fct_Pos(NumSeq,pos  )],
		seq[NumSeq].sequence[Fct_Pos(NumSeq,pos+1)],
		seq[NumSeq].sequence[Fct_Pos(NumSeq,pos+2)]));
	}

/***** returns the translation of codon c1,c2,c3 ****/
char Fct_Traduction(char c1,char c2,char c3)
	{
	char var_return;
	
	c1=LOWER(c1);
	c2=LOWER(c2);
	c3=LOWER(c3);
	switch(c1)
		{
		case 'a':switch(c2)
			{
			case 'a':
				switch(c3)
					{
					case 'a':var_return='K';break;
					case 't':var_return='N';break;
					case 'g':var_return='K';break;
					case 'c':var_return='N';break;
					default: var_return='?';break;
					}break;

			case 't':
				switch(c3)
					{
					case 'a':var_return='I';break;
					case 't':var_return='I';break;
					case 'g':var_return='M';break;
					case 'c':var_return='I';break;
					default: var_return='?';break;
					}break;

			case 'g':
				switch(c3)
					{
					case 'a':var_return='R';break;
					case 't':var_return='S';break;
					case 'g':var_return='R';break;
					case 'c':var_return='S';break;
					default: var_return='?';break;
					}break;

			case 'c':
				switch(c3)
					{
					case 'a':
					case 't':
					case 'g':
					case 'c':var_return='T';break;
					default: var_return='?';break;
					}break;
			default: var_return='?';break;
			}break;
		case 't':switch(c2)
			{
			case 'a':
				switch(c3)
					{
					case 'a':var_return=VAR_OCHRE;break;
					case 't':var_return='Y';break;
					case 'g':var_return=VAR_AMBRE;break;
					case 'c':var_return='Y';break;
					default: var_return='?';break;
					}break;

			case 't':
				switch(c3)
					{
					case 'a':var_return='L';break;
					case 't':var_return='F';break;
					case 'g':var_return='L';break;
					case 'c':var_return='F';break;
					default: var_return='?';break;
					}break;

			case 'g':
				switch(c3)
					{
					case 'a':var_return=VAR_OPALE;break;
					case 't':var_return='C';break;
					case 'g':var_return='T';break;
					case 'c':var_return='C';break;
					default: var_return='?';break;
					}break;

			case 'c':
				switch(c3)
					{
					case 'a':
					case 't':
					case 'g':
					case 'c':var_return='S';break;
					default: var_return='?';break;
					}break;

			default: var_return='?';break;
			}break;
		case 'g':switch(c2)
			{
			case 'a':
				switch(c3)
					{
					case 'a':var_return='E';break;
					case 't':var_return='D';break;
					case 'g':var_return='E';break;
					case 'c':var_return='D';break;
					default: var_return='?';break;
					}break;

			case 't':
				switch(c3)
					{
					case 'a':
					case 't':
					case 'g':
					case 'c':var_return='V';break;
					default: var_return='?';break;
					}break;

			case 'g':
				switch(c3)
					{
					case 'a':
					case 't':
					case 'g':
					case 'c':var_return='G';break;
					default: var_return='?';break;
					}break;
			case 'c':
				switch(c3)
					{
					case 'a':
					case 't':
					case 'g':
					case 'c':var_return='A';break;
					default: var_return='?';break;
					}break;
			default: var_return='?';break;
			}break;
		case 'c':switch(c2)
			{
			case 'a':
				switch(c3)
					{
					case 'a':var_return='Q';break;
					case 't':var_return='H';break;
					case 'g':var_return='Q';break;
					case 'c':var_return='H';break;
					default: var_return='?';break;
					}break;
			case 't':
				switch(c3)
					{
					case 'a':
					case 't':
					case 'g':
					case 'c':var_return='L';break;
					default: var_return='?';break;
					}break;

			case 'g':
				switch(c3)
					{
					case 'a':
					case 't':
					case 'g':
					case 'c':var_return='R';break;
					default: var_return='?';break;
					}break;
			case 'c':
				switch(c3)
					{
					case 'a':
					case 't':
					case 'g':
					case 'c':var_return='P';break;
					default: var_return='?';break;
					}break;
			default: var_return='?';break;
			}break;
		default: var_return='?';break;
		}
	return(var_return);
	}


/***********************************************************************/
void InitStgy(void)
	{FirstStgy=NULL;LastStgy=NULL;CurrStgy=NULL;TotalSolutions=0;}
/***********************************************************************/
STRUCT_STRATEGY *DeleteStgy(STRUCT_STRATEGY *Stgy)
	{
	if(Stgy==NULL) return(NULL);
	if(Stgy==FirstStgy || Stgy->prev==NULL)
		{
		if(FirstStgy==LastStgy)
			{FirstStgy=NULL;LastStgy=NULL;CurrStgy=NULL;}/*plus rien*/
		else
			{
			Stgy->next->prev=NULL;
			FirstStgy=Stgy->next;/*on enleve le premier*/
			CurrStgy=FirstStgy;
			}
		}
	else
		{
		if(Stgy==LastStgy || Stgy->next==NULL)
			{
			if(Stgy->prev==NULL) {err_printf("Error %s %d\n",__FILE__,__LINE__);}
			Stgy->prev->next=NULL;
			CurrStgy=Stgy->prev;
			LastStgy=Stgy->prev;/*on enleve le dernier*/
			}
		else
			{
			/*circuit normal*/
			Stgy->next->prev=Stgy->prev;
			Stgy->prev->next=Stgy->next;
			CurrStgy=Stgy->prev;
			}
		}
	if(Stgy==NULL) {err_printf("Error %s %d\n",__FILE__,__LINE__);}
	free(Stgy);
	TotalSolutions--;
	return(CurrStgy);
	}
/***********************************************************************/
STRUCT_STRATEGY *NewStrategy(STRUCT_STRATEGY *TheNewStrategy)
	{
	if((NewStgy=(STRUCT_STRATEGY*)malloc(sizeof(STRUCT_STRATEGY)))==NULL)
		return(NULL);
	else
		{
		memcpy(NewStgy,TheNewStrategy,sizeof(STRUCT_STRATEGY));
		if(AddStgy(NewStgy)==NULL)
			return(NULL);
		}
	return(CurrStgy);
	}
/***********************************************************************/
bool_t TestStgy(void)
	{
	if(FirstStgy==NULL || LastStgy==NULL  || CurrStgy==NULL)
		return(FALSE);
	else
		return(TRUE);
	}
/***********************************************************************/
STRUCT_STRATEGY *AddStgy(STRUCT_STRATEGY *Stgy)
	{	
	if(Stgy==NULL) return(NULL);
	if((LastStgy==NULL || FirstStgy==NULL) && TotalSolutions !=0)
		{
		err_printf("Internal error %s %d !!!\n",__FILE__,__LINE__);
		}
	if(LastStgy==NULL && FirstStgy==NULL)
		{
		TotalSolutions=1;
		Stgy->prev=NULL;
		Stgy->next=NULL;
		FirstStgy=Stgy;
		LastStgy=Stgy;
		CurrStgy=Stgy;
		
		}
	else
		{
		LastStgy->next=Stgy;
		Stgy->prev=LastStgy;
		Stgy->next=NULL;
		LastStgy=Stgy;
		CurrStgy=Stgy;
		TotalSolutions++;
		if(TotalSolutions>Preference.numMaxStgys)
			{
			SortStgys();
			SetStgy(LastStgy);
			DeleteStgy(CurrStgy);
			SetStgy(LastStgy);
			}
	/*	if(CurrStgy->prev->next==NULL) {err_printf("Total =.next est null\n");}*/
		}
	/*if(CurrStgy->prev->next==NULL) {printf("Total =%d.next est null\n");}*/
	return(Stgy);
	}
/***********************************************************************/
STRUCT_STRATEGY *SetStgy(STRUCT_STRATEGY *Stgy)	{CurrStgy=Stgy; return(CurrStgy);}
STRUCT_STRATEGY *SetFirstStgy(void)
	{
	/*if(FirstStgy==NULL && TotalSolutions!=0) printf("Error (file %s line %d)\n",__FILE__,__LINE__);*/
	CurrStgy=FirstStgy;
	return CurrStgy;
	}
STRUCT_STRATEGY *SetLastStgy(void)
	{
	if(LastStgy==NULL) err_printf("Error (file %s line %d)\n",__FILE__,__LINE__);
	CurrStgy=LastStgy;
	return LastStgy;
	}
STRUCT_STRATEGY *GetCurrStgy(void)
	{
	return CurrStgy;
	}
STRUCT_STRATEGY *NextStgy(void)
	{
	if(CurrStgy==NULL) return(NULL);
	CurrStgy=CurrStgy->next;
	return(CurrStgy);
	}
/***********************************************************************/
STRUCT_STRATEGY *PrevStgy(void)
	{
	if(CurrStgy==NULL) return(NULL);
	CurrStgy=CurrStgy->prev;
	return(CurrStgy);
	}
/***********************************************************************/
void DeleteAllStgys(void)
	{
	if(SetFirstStgy()!=NULL)
		while(DeleteStgy(GetCurrStgy())!=NULL) /* delete */;
	if(TestStgy()==TRUE)
		ERROR_USER;
	InitStgy();
	}
/***********************************************************************/
int CountStgy(int *pos)
	{
	STRUCT_STRATEGY *save;
	register int i=0;
	*pos=0;
	save=GetCurrStgy();
	
	SetFirstStgy();
	while(GetCurrStgy()!=NULL)
		{
		i++;
		if(GetCurrStgy()==save) (*pos)=i;
		NextStgy();
		}
	SetStgy(save);
	if(i!=TotalSolutions) 
		{
		err_printf("Problem:= %s %d i=%d et TotalSolutions=%d et pos=%d\n",__FILE__,__LINE__,i,TotalSolutions,*pos);
		INKEY;
		}
	return(TotalSolutions);
	}
/***********************************************************************/
void EchangeStgy(STRUCT_STRATEGY *A, STRUCT_STRATEGY *B)
	{
	STRUCT_STRATEGY C;
	STRUCT_STRATEGY *save;
	
	save=A->next;
	A->next=B->next;
	B->next=save;
	
	save=A->prev;
	A->prev=B->prev;
	B->prev=save;
		
	memcpy(&C,A,sizeof(STRUCT_STRATEGY));
	memcpy(A,B,sizeof(STRUCT_STRATEGY));
	memcpy(B,&C,sizeof(STRUCT_STRATEGY));
	}
/***********************************************************************/
void SortStgys(void)
	{
	bool_t change=TRUE;
	

if(TotalSolutions>0)
	{
	
	/* partial sites */
	change=TRUE;
	change=TRUE;
        SetFirstStgy();
        if(CurrStgy!=NULL)
	while(change==TRUE)
		{
		change=FALSE;
		SetFirstStgy();
		while(GetCurrStgy()!=NULL)
			{
			if(CurrStgy!=LastStgy)
			if(CurrStgy->Couple[SIDE_5].Partial[VECTOR] + CurrStgy->Couple[SIDE_5].Partial[INSERT] +
			   CurrStgy->Couple[SIDE_3].Partial[VECTOR] + CurrStgy->Couple[SIDE_3].Partial[INSERT] >
			   CurrStgy->next->Couple[SIDE_5].Partial[VECTOR] + CurrStgy->next->Couple[SIDE_5].Partial[INSERT] +
			   CurrStgy->next->Couple[SIDE_3].Partial[VECTOR] + CurrStgy->next->Couple[SIDE_3].Partial[INSERT])
				{
				EchangeStgy(CurrStgy, CurrStgy->next);
				change=TRUE;
				}
			NextStgy();
			}
		}

	/* phosphatase */
	change=TRUE;
	SetFirstStgy();	
	if(CurrStgy!=NULL)
	if(CurrStgy->type==SUB_CLONING)
		{
		while(change==TRUE)
			{
			change=FALSE;
			SetFirstStgy();
			while(GetCurrStgy()!=NULL)
				{
				if(CurrStgy!=LastStgy)
				if(CurrStgy->use_CIP==TRUE && CurrStgy->next->use_CIP==FALSE)
					{
					EchangeStgy(CurrStgy,CurrStgy->next);
					change=TRUE;
					}
				NextStgy();
				}
			}
		}

	/* treatment */
	change=TRUE;
        SetFirstStgy();
        if(CurrStgy!=NULL)	
	while(change==TRUE)
		{
		change=FALSE;
		SetFirstStgy();
		while(GetCurrStgy()!=NULL)
			{
			if(CurrStgy!=LastStgy)
				{
				if(CurrStgy->Couple[SIDE_5].Treatment[VECTOR] + CurrStgy->Couple[SIDE_5].Treatment[INSERT] +
	  			   CurrStgy->Couple[SIDE_3].Treatment[VECTOR] + CurrStgy->Couple[SIDE_3].Treatment[INSERT] >
				   CurrStgy->next->Couple[SIDE_5].Treatment[VECTOR] + CurrStgy->next->Couple[SIDE_5].Treatment[INSERT] +
				   CurrStgy->next->Couple[SIDE_3].Treatment[VECTOR] + CurrStgy->next->Couple[SIDE_3].Treatment[INSERT])
					{
					EchangeStgy(CurrStgy, CurrStgy->next);
					change=TRUE;
					}
				}
			NextStgy();
			}
		}
	
	/* commercial */
	change=TRUE;
        SetFirstStgy();
        if(CurrStgy!=NULL)
	while(change==TRUE)
		{
		change=FALSE;
		SetFirstStgy();
		while(GetCurrStgy()!=NULL)
			{
			if(CurrStgy!=LastStgy)
				{
				if(Enzymes[Sites[CurrStgy->Couple[SIDE_5].NumSite[VECTOR]].NumEnz].corp + 
				   Enzymes[Sites[CurrStgy->Couple[SIDE_5].NumSite[INSERT]].NumEnz].corp +
	  			   Enzymes[Sites[CurrStgy->Couple[SIDE_3].NumSite[VECTOR]].NumEnz].corp +
	  			   Enzymes[Sites[CurrStgy->Couple[SIDE_3].NumSite[INSERT]].NumEnz].corp <
				   Enzymes[Sites[CurrStgy->next->Couple[SIDE_5].NumSite[VECTOR]].NumEnz].corp +
				   Enzymes[Sites[CurrStgy->next->Couple[SIDE_5].NumSite[INSERT]].NumEnz].corp +
				   Enzymes[Sites[CurrStgy->next->Couple[SIDE_3].NumSite[VECTOR]].NumEnz].corp +
				   Enzymes[Sites[CurrStgy->next->Couple[SIDE_3].NumSite[INSERT]].NumEnz].corp)
					{
					EchangeStgy(CurrStgy, CurrStgy->next);
					change=TRUE;
					}
				}
			NextStgy();
			}
		}
	}	
}





/****************************************************************
      this procedure finds the restriction pattern of a plasmid
                   and diplays it on screen
*****************************************************************/
void Digest(short mode,FILE *out,int NumSeq,int NumEnz1, int NumEnz2, int AlertPartial)
	{

	/***********************************************************/
	
	STRUCT_FRAGMENT Topaze;
	
	int i,j,L1,L2,var_change=FALSE;
	/***********************************************************/
	NumFragment=0;
	/* diplay enzymes used */

	fprintf(out,"%s [%s] (%d) ", Enzymes[Sites[NumEnz1].NumEnz].nom,Enzymes[Sites[NumEnz1].NumEnz].site_complet,Sites[NumEnz1].Loc);
	if(Sites[NumEnz1].NumEnz != Sites[NumEnz2].NumEnz)
		fprintf(out,"and %s [%s] (%d) ", Enzymes[Sites[NumEnz2].NumEnz].nom,Enzymes[Sites[NumEnz2].NumEnz].site_complet,Sites[NumEnz2].Loc);
	fprintf(out,"%s digestion.\n", seq[NumSeq].FICHIER_ADN);
	fprintf(out,"<BR>\n");
	/* scan all sites */
	for(i=0;i<nbr_sites;i++)
	{
	/*
	
	0-----------------------------------npb
	      |                |        |
	      start1           end1     |
                           start2   end2
                                    start0
          end0
	*/
	/* if site i cuts the plasmid */
	if(  Sites[i].NumSeq==NumSeq && (
		(strcmp(Enzymes[Sites[NumEnz1].NumEnz].nom,Enzymes[Sites[i].NumEnz].nom)==ARE_IDENTIC && Are_Same_Asymetric(NumEnz1,i)==TRUE )|| 
		(strcmp(Enzymes[Sites[NumEnz2].NumEnz].nom,Enzymes[Sites[i].NumEnz].nom)==ARE_IDENTIC && Are_Same_Asymetric(NumEnz2,i)==TRUE)))
		{
		/* alloc a new fragment */
		if((TheFragment=(STRUCT_FRAGMENT *)realloc(TheFragment,(NumFragment+1)*sizeof(STRUCT_FRAGMENT )))==NULL)
			{
			err_printf("Memory Full ! NumFragment=%d\n\n",NumFragment);
			break;
			}
		else
			{
			/* if there no fragment has been defined alloc the first fragment...*/
			if(NumFragment==0)
				{
				TheFragment[0].Start=i;
				TheFragment[0].End=i;
				NumFragment=1;
				}
			/* if there is just one fragment defined*/
			else if(NumFragment==1)
				{
				if(Sites[i].Loc > Sites[TheFragment[0].Start].Loc)
					{
					TheFragment[1].Start=TheFragment[0].End;
					TheFragment[1].End=i;
					TheFragment[0].Start=i;
				/*  TheFragment[0].End doesn't change*/
					}
				else if(Sites[i].Loc < Sites[TheFragment[0].Start].Loc)
					{
					TheFragment[1].Start=i;
					TheFragment[1].End=TheFragment[0].Start;
					/*TheFragment[0].Start */
					TheFragment[0].End=i;
					}
				/* note that element[0] is overlapping */
				NumFragment=2;
				}/* end else frag==1 */
			/* else scan all fragments and find where is localized site i */
			else for(j=0;j<NumFragment;j++)
				{
				/* if the site is in a fragment overlapping the boundaries*/
				if(j==0)
					{
					if(Sites[TheFragment[0].Start].Loc<Sites[TheFragment[0].End  ].Loc)
						err_printf("ERROR !\n");
					if(Sites[i].Loc > Sites[TheFragment[0].Start].Loc &&
					   Sites[i].Loc > Sites[TheFragment[0].End  ].Loc)
						{
						TheFragment[NumFragment].Start=TheFragment[0].Start;
						TheFragment[NumFragment].End=i;
						TheFragment[0].Start=i;
					 /* TheFragment[0].End Doesn't change*/					
						NumFragment++;
						break;
						}
					else
					if (Sites[i].Loc < Sites[TheFragment[0].Start].Loc &&
					    Sites[i].Loc < Sites[TheFragment[0].End  ].Loc)
						{
						TheFragment[NumFragment].Start=i;
					    TheFragment[NumFragment].End=TheFragment[0].End;
					 /* TheFragment[0].start Doesn't change*/
						TheFragment[0].End=i;			
						NumFragment++;
						break;
						}
					}/* end j==0 */
				else if(Sites[i].Loc > Sites[TheFragment[j].Start].Loc &&
					    Sites[i].Loc < Sites[TheFragment[j].End  ].Loc)
					{
					/*
					-------[j].Start---[NumFrag].Loc--[j].End----------
					*/
				  /*TheFragment[j].Start don't change */
				  	TheFragment[NumFragment].Start=i;
				 	TheFragment[NumFragment].End=TheFragment[j].End;
				  	TheFragment[j].End=i;
					NumFragment++;
					break;
					}/* end else */
				}/* end for j=0 */
			}/* end else alloc */
		}/* end if */
	}/* end for i*/
	/* sort by length */
	var_change=FALSE;
	
	
	while(var_change==FALSE)
		{
		var_change=TRUE;
		for(i=0;i<NumFragment-1;i++)
			{
			L1=(Sites[TheFragment[i].End].Loc- Sites[TheFragment[i].Start].Loc);
			L1=(L1>0?L1:seq[NumSeq].npb+L1);
			L2=(Sites[TheFragment[i+1].End].Loc- Sites[TheFragment[i+1].Start].Loc);
			L2=(L2>0?L2:seq[NumSeq].npb+L2);
			if( L1 < L2 )
				{
				var_change=FALSE;
				memcpy(&Topaze,&TheFragment[i],sizeof(STRUCT_FRAGMENT));
				memcpy(&TheFragment[i],&TheFragment[i+1],sizeof(STRUCT_FRAGMENT));
				memcpy(&TheFragment[i+1],&Topaze,sizeof(STRUCT_FRAGMENT));
				}
			}
		}
	
	fprintf(out,"<TABLE ALIGN=\"CENTER\" BORDER=1>");
	for(j=0;j<NumFragment;j++)
		{
		fprintf(out,"<TR><TD>");
		fprintf(out," %d<TD>",j+1);
		var_change=(Sites[TheFragment[j].End].Loc- Sites[TheFragment[j].Start].Loc);
		fprintf(out,"%d pb<TD>",(var_change>0?var_change:seq[NumSeq].npb+var_change));
		fprintf(out,"%s  %d - %s  %d",
			Enzymes[Sites[TheFragment[j].Start].NumEnz].nom,
			Sites[TheFragment[j].Start].Loc,
			Enzymes[Sites[TheFragment[j].End].NumEnz].nom,
			Sites[TheFragment[j].End].Loc);
		}
	fprintf(out,"</TABLE>");
fprintf(out,"<BR>\n");
if(AlertPartial>0)
	fprintf(out,"<BLINK>Beware this strategy needs partial digestion(s)</BLINK>");
if((TheFragment=(STRUCT_FRAGMENT *)realloc(TheFragment,(1)*sizeof(STRUCT_FRAGMENT )))==NULL)
	exit(0);
}/* end void */


/********************************************************
 general function to see if two sites are compatible
	- with or without treatment
	- with or without frame pression
********************************************************/
int Fct_Test(int NumSiteA,int NumSiteB,int Treatment ,int _Frame, STRUCT_STRATEGY_2 *Stg_2)
	{
	/* test if compatible with treatment */
	if(Fct_Compatible(NumSiteA,NumSiteB,Treatment,Stg_2)==TRUE)
		{
		if(_Frame==FALSE)
			{
			return(TRUE);
			}
		else
			{
			/* test if compatible with frame pression */
			return(Are_in_Frame(Stg_2));
			}
		}
	else
		return(FALSE);
	}

void Find_stop_codon(short mode, FILE *out,int NumSite,int var_limit,int var_side)
	{
/********************************************************
 find the first stop codon from a site , on the 5' side or
	on the 3' side, to var_limit. This is just a displayed information that
	may be an important information for the biologist:
	"What is the use of cloning in frame if there is a stop
	codon a few bases after the ligation site ?"
********************************************************/
	int i,var_Seq;
	if(Preference.side_3==TRUE || Preference.side_5==TRUE)
		{
		var_Seq=Sites[NumSite].NumSeq;
		/* seach stop codon on the left */
		if(var_side==SIDE_5)
			{
			/* download sequence */
			for(i=Sites[NumSite].Loc-3;i>=var_limit;i--)
			   if(Fct_Frame(var_Seq,i)==IS_IN_FRAME) /* if this position is in frame */
				if(Is_Codon_Stop(Fct_Traduction(seq[var_Seq].sequence[i],
												seq[var_Seq].sequence[i+1],
												seq[var_Seq].sequence[i+2]))==TRUE)
					{
	fprintf(out,"The first stop codon detected <B>BEFORE</B> the %s site (%d) is localized at position <B>%d</B> on %s.<BR>",
						Enzymes[Sites[NumSite].NumEnz].nom,Sites[NumSite].Loc,i,((var_Seq==0)?"vector":"insert"));
					break;
					}
			}
		else if(var_side==SIDE_3) /* seach stop codon on the right */
			{
			/* upload sequence */
			for(i=Sites[NumSite].Loc+Enzymes[Sites[NumSite].NumEnz].taille_site-1;i<=var_limit;i++)
			   if(Fct_Frame(var_Seq,i)==IS_IN_FRAME) /* if this position is in frame */
				if(Is_Codon_Stop(Fct_Traduction(seq[var_Seq].sequence[i],
												seq[var_Seq].sequence[i+1],
												seq[var_Seq].sequence[i+2]))==TRUE)
					{
fprintf(out,"The first stop codon detected <B>AFTER</B> the %s site (%d) is localized at position <B>%d</B> on %s.<BR>",
						Enzymes[Sites[NumSite].NumEnz].nom,Sites[NumSite].Loc,i,((var_Seq==VECTOR)?"vector":"insert"));
					break;
					}
			}
		}
	}

			
int Fct_N_sites(int Numenz,int NumSeq,int _min,int _max,int Boolean)
	{
/************************************************************************
 this is one of the most used function in this program
	if Boolean==TRUE
		it returns the number of sites cutes by Enzyme n°Numenz in the 
			sequence fragment between _min et _max.
	if Boolean==FALSE
		it returns the number of sites cutes by Enzyme n°Numenz OUT OF the 
			sequence fragment between _min et _max.
************************************************************************/
	register int i=0,j=0;
	int var_result=0;
	
	/* scan all sites */
	for(i=0;i<nbr_sites;i++)
		{
		/* if the site is cuted by NumEnz and is one NumSeq */
		if(Sites[i].NumEnz==Numenz && Sites[i].NumSeq==NumSeq)
			{
			if(Sites[i].Loc>=_min && Sites[i].Loc<=_max && Boolean==TRUE)
				var_result++; /* site is in the fragment */
			else if((Sites[i].Loc<=_min || Sites[i].Loc>=_max) && Boolean==FALSE)
				var_result++; /* site is out of the fragment */
			}
		}
	/* then there is the problem of the non palindromic enzymes that are coded
		twice in the array *Enzymes (See: Get_Rebase) */
	if(Enzymes[Numenz].palindromic!=TRUE)
		{
		/* is 'i' is the same  enzyme than previous i-1 */
		if(Numenz+1<nbr_enzyme)
			{
			if(strcmp(Enzymes[Numenz].site_complet,Enzymes[Numenz+1].site_complet)==ARE_IDENTIC)
				j=Numenz+1;
			}
		if(Numenz-1>=0) /* or is 'i' is the same  enzyme than next i+1 */
			{
			if(strcmp(Enzymes[Numenz].site_complet,Enzymes[Numenz-1].site_complet)==ARE_IDENTIC)
				j=Numenz-1;
			}
		/* add the sites cuted by the second "same" enzyme */
		for(i=0;i<nbr_sites;i++)
			{
			if(Sites[i].NumEnz==j && Sites[i].NumSeq==NumSeq)
				{
				if(Sites[i].Loc>=_min && Sites[i].Loc<=_max && Boolean==TRUE)
					var_result++;/* site is in the fragment */
				else if((Sites[i].Loc<=_min || Sites[i].Loc>=_max) && Boolean==FALSE)
					var_result++;/* site is out of the fragment */
				}
			}
		}
	return(var_result);
	}

/*********************************************************************************/
void Init_Enz(void)
	{
	/***************************************************************************
	 this command discard all enzymes that don't have any site in the VECTOR or the INSERT
		cloning boxes (It accelerate the search for cloning strategies)
		It sets field *Enzymes[i].select[sequence] to FALSE
	*****************************************************************************/
	register int i;
	
	/* scan all enzymes */
	for(i=0;i<nbr_enzyme;i++)
		{
		Enzymes[i].select[VECTOR]=TRUE; /* TRUE by default */
		Enzymes[i].select[INSERT]=TRUE; /* TRUE by default */
		/* discard if no enz in vector */
		if(Enzymes[i].palindromic==TRUE)
			{
			/* discard if no Enz in vector cloning box*/
			if(Fct_N_sites(i,VECTOR,seq[VECTOR].var_min,seq[VECTOR].var_max,TRUE)==0)
				Enzymes[i].select[VECTOR]=FALSE;
			/* discard if no Enz in insert cloning box */
			if( Fct_N_sites(i,INSERT,seq[INSERT].var_min,seq[INSERT].var_min_int,TRUE)==0	&&
				Fct_N_sites(i,INSERT,seq[INSERT].var_max_int,seq[INSERT].var_max,TRUE)==0)
				Enzymes[i].select[INSERT]=FALSE;

			}
	/* then there is the problem of the non palindromic enzymes that are coded
		twice in the array *Enzymes (See: Get_Rebase) */

		else if(strcmp(Enzymes[i].site_complet,Enzymes[i+1].site_complet)==ARE_IDENTIC)
			{
			Enzymes[i+1].select[VECTOR]=TRUE; /* TRUE by default */
			Enzymes[i+1].select[INSERT]=TRUE; /* TRUE by default */
			/* discard if no Enz in vector cloning box*/
			if(	Fct_N_sites(i,VECTOR,seq[VECTOR].var_min,seq[VECTOR].var_max,TRUE)==0 &&
				Fct_N_sites(i+1,VECTOR,seq[VECTOR].var_min,seq[VECTOR].var_max,TRUE)==0)
					{
					Enzymes[i].select[VECTOR]=FALSE;
					Enzymes[i+1].select[VECTOR]=FALSE;
					}
			/* discard if no Enz in insert cloning box*/
			if( Fct_N_sites(i,INSERT,seq[INSERT].var_min,seq[INSERT].var_min_int,TRUE)==0	&&
				Fct_N_sites(i,INSERT,seq[INSERT].var_max_int,seq[INSERT].var_max,TRUE)==0	&&
				Fct_N_sites(i+1,INSERT,seq[INSERT].var_min,seq[INSERT].var_min_int,TRUE)==0	&&
				Fct_N_sites(i+1,INSERT,seq[INSERT].var_max_int,seq[INSERT].var_max,TRUE)==0)
					{
					Enzymes[i].select[INSERT]=FALSE;
					Enzymes[i+1].select[INSERT]=FALSE;
					}
			i++; /* next enzyme have just been done */
			}
		}
	}
/*********************************************************************************/
bool_t Detect_Stop(FILE *out,short mode,int NumSiteA, int NumSiteB,STRUCT_STRATEGY_2 *Stg_2)
	{
	/***************************************************************
	This function detects if a stop codon was created after ligation
	of site NumSiteA with NumSiteB according to Strategy 
	***************************************************************/

	int SeqA,SeqB,stop_was_detected=FALSE;
	int Loc5_3_5,Loc5_5_3,Loc3_3_5,Loc3_5_3,overhang_A,overhang_B;
	char c1,c2,c3,c4; /* bases implicated in the ligation */
	char chara='?',charb='?',charc='?';/* bases displayed at the end  if a stop codon was found */
	
	/* c1 to c4 are the four bases that can be implicated in the creation of a 
	stop codon after ligation
	         cccc
	         1234
	AAAAAAAAAAABBBBBBB
	         --- <-stop codon if c1 is in frame ?
	          --- <-stop codon if c2 is in frame ?
	*/
	
	/***************************************************************************/
	/* no matter to detect a stop if frame is not important */
	if(Preference.side_3==FALSE && Preference.side_5==FALSE) return(FALSE);

	/* if Test_Trans==FALSE, then the function is called by FrameShift, or Delta Frame
		else it is called by Sub_Cloning */
	SeqA=((Stg_2->Test_Trans==FALSE) ? VECTOR :Sites[NumSiteA].NumSeq);
	SeqB=((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[NumSiteB].NumSeq);

	/* copy valors , it is easier to understand the program */
	Loc5_5_3= Stg_2->Pos5_3[SeqA];
	Loc5_3_5= Stg_2->Pos3_5[SeqA];
	Loc3_5_3= Stg_2->Pos5_3[SeqB];
	Loc3_3_5= Stg_2->Pos3_5[SeqB];

	/* if Test_Trans==FALSE, then the function works only with INSERT, previous line
	 was just used to extract the valors of Strategy  */
	SeqA=((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[NumSiteA].NumSeq);
	SeqB=((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[NumSiteB].NumSeq);

	/* overhang length */
	overhang_A= Loc5_5_3-Loc5_3_5;
	overhang_B= Loc3_5_3-Loc3_3_5;
	/***************************************************************************/
	if( overhang_A == overhang_B && Sign(overhang_A)==Sign(overhang_B))
		{
		/**********************************
		5-------- ------------------3
		3--------- -----------------5
		or
		5--------- -----------------3
		3-------- ------------------5
		or
		5--------- -----------------3
		3--------- -----------------5
		**********************************/

		c1= seq[SeqA].sequence[Loc5_5_3-2];
		c2= seq[SeqA].sequence[Loc5_5_3-1];
		c3= seq[SeqB].sequence[Loc3_5_3];
		c4= seq[SeqB].sequence[Loc3_5_3+1];
		
		/* if c2 is in frame */
		if(Fct_Frame(SeqA,Loc5_5_3-1)==IS_IN_FRAME)
			{
			if(Is_Codon_Stop(Fct_Traduction(c2,c3,c4))==TRUE)
				{stop_was_detected=TRUE;
				chara=LOWER(c2);charb=UPPER(c3);charc=UPPER(c4);}
			}
		/* if c1 is in frame */
		else if(Fct_Frame(SeqA,Loc5_5_3-2)==IS_IN_FRAME)
			{
			if(Is_Codon_Stop(Fct_Traduction(c1,c2,c3))==TRUE)
				{stop_was_detected=TRUE;
				chara=LOWER(c1);charb=LOWER(c2);charc=UPPER(c3);}
			}
		}
	/***************************************************************************/
	else if((overhang_A < overhang_B) && Sign(overhang_A)==Sign(overhang_B))
		{
		/**********************************
		5------        -------------------3
		3--------                       --5
		or
		5--------                       --3
		3------        -------------------5
		**********************************/
		
		c1= seq[SeqA].sequence[MIN(Loc5_5_3,Loc5_3_5)-2];
		c2= seq[SeqA].sequence[MIN(Loc5_5_3,Loc5_3_5)-1];
		c3= seq[SeqB].sequence[MIN(Loc3_5_3,Loc3_3_5)  ];
		c4= seq[SeqB].sequence[MIN(Loc3_5_3,Loc3_3_5)+1];
		/* if c2 is in frame */
		if(Fct_Frame(SeqA,MIN(Loc5_5_3,Loc5_3_5)-1)==IS_IN_FRAME)
			{
			if(Is_Codon_Stop(Fct_Traduction(c2,c3,c4))==TRUE)
				{stop_was_detected=TRUE;
				chara=LOWER(c2);charb=UPPER(c3);charc=UPPER(c4);}
			}
		/* if c1 is in frame */
		else if(Fct_Frame(SeqA,MIN(Loc5_5_3,Loc5_3_5)-2)==IS_IN_FRAME)
			{
			if(Is_Codon_Stop(Fct_Traduction(c1,c2,c3))==TRUE)
				{stop_was_detected=TRUE;
				chara=LOWER(c1);charb=LOWER(c2);charc=UPPER(c3);}
			}
		}
	else if((overhang_A > overhang_B) && Sign(overhang_A)==Sign(overhang_B))
		{
		/**********************************
		5-------------                   -3
		3-                             ---5
		or
		5-                             ---3
		3-------------                   -5
		**********************************/
		c1= seq[SeqA].sequence[MAX(Loc5_5_3,Loc5_3_5)-2];
		c2= seq[SeqA].sequence[MAX(Loc5_5_3,Loc5_3_5)-1];
		c3= seq[SeqB].sequence[MAX(Loc3_5_3,Loc3_3_5)  ];
		c4= seq[SeqB].sequence[MAX(Loc3_5_3,Loc3_3_5)+1];
		
		/* if c2 is in frame */
		if(Fct_Frame(SeqA,MAX(Loc5_5_3,Loc5_3_5)-1)==IS_IN_FRAME)
			{
			if(Is_Codon_Stop(Fct_Traduction(c2,c3,c4))==TRUE)
				{stop_was_detected=TRUE;
				chara=LOWER(c2);charb=UPPER(c3);charc=UPPER(c4);}
			}
		/* if c1 is in frame */
		else if(Fct_Frame(SeqA,MAX(Loc5_5_3,Loc5_3_5)-2)==IS_IN_FRAME)
			{
			if(Is_Codon_Stop(Fct_Traduction(c1,c2,c3))==TRUE)
				{stop_was_detected=TRUE;
				chara=LOWER(c1);charb=LOWER(c2);charc=UPPER(c3);}
			}
		}
	else return(FALSE);
	
	if(stop_was_detected==TRUE)
		{
		fprintf(out,"<BLINK>BEWARE: A stop codon is created after the ligation [%c%c%c]!..</BLINK>",chara,charb,charc);
		return(TRUE);
		}
	return(FALSE);
	}
	
/**************************************************************************************/	
int Test_Use_CIP(int NumSeq,STRUCT_STRATEGY Stg)
	{
	
	/* Test to see if vector can be self ligated  */
	if(NumSeq!=VECTOR) /* should always be used with the vector */
		return(FALSE);
	return(Fct_Compatible2( VECTOR, Stg.Couple[SIDE_5].Pos5_3[VEC_5],
								 	Stg.Couple[SIDE_5].Pos3_5[VEC_5],
						 	VECTOR, Stg.Couple[SIDE_3].Pos5_3[VEC_3],
								 	Stg.Couple[SIDE_3].Pos3_5[VEC_3]));
	/* TRUE if enzyme on SIDE_5 is compatible with enzyme on SIDE_3 */
	}
	
/**************************************************************************************/

void Treatments(int *Treat_5, int *Treat_3,int Type5,int Type3)
{
	/*********************************************************************
	This function look the two treatments of the strategy used on one
	sequence. It simplifies the treatment used.
		For example:
		-there is no use to use modifying polymerase
		if the 2 endonucleases used generate blunt ends.
		- there is no use to make a 2 time treatment without/with
		 modifying polymerase if one of the 2 endonucleases used
		 generate blunt ends.
	*********************************************************************/



if(*Treat_5==NO_TREATMENT)
	{
	switch(Type5)
		{
		case(TYPE_BLUNT):
			switch(*Treat_5)
				{
				case(NO_TREATMENT):
					switch(Type3)
						{
						case(TYPE_BLUNT):	*Treat_5=NO_TREATMENT;*Treat_3=NO_TREATMENT;break;
						default:			*Treat_5=T4_TREATMENT;*Treat_3=T4_TREATMENT;break;								}break;
				case(T4_TREATMENT):
					switch(Type3)
						{
						case(TYPE_BLUNT):	*Treat_5=NO_TREATMENT;*Treat_3=NO_TREATMENT;break;
						default:			*Treat_5=T4_TREATMENT;*Treat_3=T4_TREATMENT;break;								}break;
				}break;
		default:
			switch(*Treat_5)
				{
				case(NO_TREATMENT):
					switch(Type3)
						{
						case(TYPE_BLUNT):	*Treat_5=NO_TREATMENT;*Treat_3=NO_TREATMENT;break;
						default:			*Treat_5=NO_TREATMENT;*Treat_3=T4_TREATMENT;break;								}break;
				case(T4_TREATMENT):
					switch(Type3)
						{
						case(TYPE_BLUNT):	*Treat_5=T4_TREATMENT;*Treat_3=T4_TREATMENT;break;
						default:			*Treat_5=T4_TREATMENT;*Treat_3=T4_TREATMENT;break;								}break;
				}break;
		}
	}
else /*treat 5 = T4 treatment */
	{
	switch(Type5)
		{
		case(TYPE_BLUNT):
			switch(*Treat_3)
				{
				case(NO_TREATMENT):
					switch(Type3)
						{
						case(TYPE_BLUNT):	*Treat_5=NO_TREATMENT;*Treat_3=NO_TREATMENT;break;
						default:			*Treat_5=NO_TREATMENT;*Treat_3=NO_TREATMENT;break;								}break;
				case(T4_TREATMENT):
					switch(Type3)
						{
						case(TYPE_BLUNT):	*Treat_5=NO_TREATMENT;*Treat_3=NO_TREATMENT;break;
						default:			*Treat_5=T4_TREATMENT;*Treat_3=T4_TREATMENT;break;								}break;
				}break;
		default:
			switch(*Treat_3)
				{
				case(NO_TREATMENT):
					switch(Type3)
						{
						case(TYPE_BLUNT):	*Treat_5=T4_TREATMENT;*Treat_3=T4_TREATMENT;break;
						default:			*Treat_5=T4_TREATMENT;*Treat_3=NO_TREATMENT;break;								}break;
				case(T4_TREATMENT):
					switch(Type3)
						{
						case(TYPE_BLUNT):	*Treat_5=T4_TREATMENT;*Treat_3=T4_TREATMENT;break;
						default:			*Treat_5=T4_TREATMENT;*Treat_3=T4_TREATMENT;break;								}break;
				}break;
		}
	}

}




int Fct_Compatible_No_Treatment (int NumSiteA,int NumSiteB, STRUCT_STRATEGY_2 *Stg_2)
	{
	/***************************************************************************/
	/* this function returns TRUE if two sites are compatible without treatment*/
	/***************************************************************************/
	int var_return=FALSE,typ1,typ2,NumEnzA,NumEnzB,SeqA,SeqB;
	bool_t vara=FALSE;
	
	NumEnzA	= Sites[NumSiteA].NumEnz;
	NumEnzB = Sites[NumSiteB].NumEnz;
	
	/* if Test_Trans==FALSE, then the function works only with INSERT 
	if Test_Trans==FALSE, then the function is called by FrameShift,
	or Delta Frame else it is called by Sub_Cloning */
	
	SeqA=((Stg_2->Test_Trans==FALSE) ? VECTOR :Sites[NumSiteA].NumSeq);
	SeqB=((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[NumSiteB].NumSeq);
	
	typ1 = Fct_type(Enzymes[NumEnzA]);
	typ2 = Fct_type(Enzymes[NumEnzB]);
	
	/* each type combinaison is scanned */
	switch(typ1)
		{
		case(TYPE_BLUNT):
			switch(typ2)
				{
				case(TYPE_BLUNT):return(TRUE);break;
				case(TYPE_3_OVER):return(FALSE);break;
				case(TYPE_5_OVER):return(FALSE);break;
				}break;
		case(TYPE_3_OVER):
			switch(typ2)
				{
				case(TYPE_BLUNT):return(FALSE);break;
				case(TYPE_3_OVER):
					{
					/* 3' overhang and 3' overhang may be compatible, let's
						check it on the level of the sequence with Fct_Compatible2 */
					vara=Fct_Compatible2( ((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[NumSiteA].NumSeq),Stg_2->Pos5_3[SeqA],Stg_2->Pos3_5[SeqA],
						 				  ((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[NumSiteB].NumSeq),Stg_2->Pos5_3[SeqB],Stg_2->Pos3_5[SeqB]);
						 				 
					return(vara);
					}break;
				case(TYPE_5_OVER):return(FALSE);break;
				}
			break;
		case(TYPE_5_OVER):
			switch(typ2)
				{
				case(TYPE_BLUNT):return(FALSE);break;
				case(TYPE_3_OVER):return(FALSE);break;
				case(TYPE_5_OVER):
					{
					/* 5' overhang and 5' overhang may be compatible, let's
						check it on the level of the sequence with Fct_Compatible2 */
					vara=Fct_Compatible2( ((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[NumSiteA].NumSeq),Stg_2->Pos5_3[SeqA],Stg_2->Pos3_5[SeqA],
						 				  ((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[NumSiteB].NumSeq),Stg_2->Pos5_3[SeqB],Stg_2->Pos3_5[SeqB]);
					return(vara);
					}break;
				}
			break;
		}
	return(var_return);
	}

/********************************************************/
int Fct_Compatible(int NumSiteA,int NumSiteB,int Treatment,STRUCT_STRATEGY_2 *Stg_2)
	{
	/*________________________________________________________________________
	this function returns TRUE if two sites are compatible
	- polymerase treatment is explored
	- in this function STRUCT_STRATEGY_2.Pos5_3 and others are calculated
	
	if there is no treatment:
		Pos5_3= Site.Loc+ Enzyme.pos5_3
		Pos3_5= Site.Loc+ Enzyme.pos3_5
	
	if overhang are filled by polymerase.
		(here is an exemple with an EcoR1 site but you can check that it 
		also works with Pst (5' overhanged))
		
		on the left (5'):
			pos5_3--------->
			5'NNNNNNNNNNNN_G....
			3'NNNNNNNNNNNN_CTTAA
			pos3_5------------->
	
			Pos5_3= Site.Loc+ Enzyme.pos3_5
			Pos3_5= Site.Loc+ Enzyme.pos3_5
			
		on the left (3'):
			pos5_3-->
					AATTC_NNNNNNN 3'
					....G_NNNNNNN 5'
			pos3_5------>
				
			Pos5_3= Site.Loc+ Enzyme.pos5_3
			Pos3_5= Site.Loc+ Enzyme.pos5_3
	 
	________________________________________________________________________*/

	if(Treatment==NO_TREATMENT)
		{
		Stg_2->Pos5_3[((Stg_2->Test_Trans==FALSE) ? VECTOR :Sites[NumSiteA].NumSeq)] = Sites[NumSiteA].Loc + Enzymes[Sites[NumSiteA].NumEnz].pos5_3;
		Stg_2->Pos5_3[((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[NumSiteB].NumSeq)] = Sites[NumSiteB].Loc + Enzymes[Sites[NumSiteB].NumEnz].pos5_3;
		Stg_2->Pos3_5[((Stg_2->Test_Trans==FALSE) ? VECTOR :Sites[NumSiteA].NumSeq)] = Sites[NumSiteA].Loc + Enzymes[Sites[NumSiteA].NumEnz].pos3_5;
		Stg_2->Pos3_5[((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[NumSiteB].NumSeq)] = Sites[NumSiteB].Loc + Enzymes[Sites[NumSiteB].NumEnz].pos3_5;		
		return(Fct_Compatible_No_Treatment (NumSiteA,NumSiteB,Stg_2));
		}
	else if(Treatment==T4_TREATMENT)
		{
		Stg_2->Pos5_3[((Stg_2->Test_Trans==FALSE) ? VECTOR :Sites[NumSiteA].NumSeq)] = Sites[NumSiteA].Loc + Enzymes[Sites[NumSiteA].NumEnz].pos3_5;
		Stg_2->Pos5_3[((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[NumSiteB].NumSeq)] = Sites[NumSiteB].Loc + Enzymes[Sites[NumSiteB].NumEnz].pos5_3;
		Stg_2->Pos3_5[((Stg_2->Test_Trans==FALSE) ? VECTOR :Sites[NumSiteA].NumSeq)] = Sites[NumSiteA].Loc + Enzymes[Sites[NumSiteA].NumEnz].pos3_5;
		Stg_2->Pos3_5[((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[NumSiteB].NumSeq)] = Sites[NumSiteB].Loc + Enzymes[Sites[NumSiteB].NumEnz].pos5_3;
		return(TRUE);
		}
	else 
		return(FALSE);
	}



/********************************************************/
bool_t Are_in_Frame(STRUCT_STRATEGY_2 *Stg_2)
	{
	/****************************************************/
	/* return if two sites , when ligated, are in frame */
	/****************************************************/
	
	int SeqA,SeqB;
	int Loc5_3_5,Loc5_5_3,Loc3_3_5,Loc3_5_3,overhang_A,overhang_B;
	/* frame is conserved if user don't have selected the frame selection pressure */
	if(Preference.side_3==FALSE && Preference.side_5==FALSE) return(TRUE);

	/* if Test_Trans==FALSE, then the function is called by FrameShift, or Delta Frame
		else it is called by Sub_Cloning */
		
	SeqA=((Stg_2->Test_Trans==FALSE) ? VECTOR :Sites[Stg_2->NumSite[VECTOR]].NumSeq);
	SeqB=((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[Stg_2->NumSite[INSERT]].NumSeq);

	/* copy valors , it is easier to understand the program */
	Loc5_5_3= Stg_2->Pos5_3[SeqA]; /* distance 5->3 on side 5' of SeqA */
	Loc5_3_5= Stg_2->Pos3_5[SeqA]; /* etc ... */
	Loc3_5_3= Stg_2->Pos5_3[SeqB];
	Loc3_3_5= Stg_2->Pos3_5[SeqB];

	/* if Test_Trans==FALSE, then the function works only with INSERT, previous line
	 was just used to extract the valors of Strategy  */
	SeqA=((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[Stg_2->NumSite[VECTOR]].NumSeq);
	SeqB=((Stg_2->Test_Trans==FALSE) ? INSERT :Sites[Stg_2->NumSite[INSERT]].NumSeq);

	/* overhang length */
	overhang_A= Loc5_5_3-Loc5_3_5;
	overhang_B= Loc3_5_3-Loc3_3_5;

		
	if( overhang_A == overhang_B && Sign(overhang_A)==Sign(overhang_B))
		{
		/**********************************
		"an image count for 1000 words" [Confucius]
		
		5------    -----------------------3
		3--------    ---------------------5
		or
		5--------    ---------------------3
		3------    -----------------------5
		**********************************/
		if( Fct_Frame(SeqA,Fct_Pos(SeqA,Loc5_5_3))==
			Fct_Frame(SeqB,Fct_Pos(SeqB,Loc3_5_3)))
				return(TRUE);
		else
			return(FALSE);
		}
	else if((overhang_A < overhang_B) && Sign(overhang_A)==Sign(overhang_B))
		{
		/**********************************
		5------        -------------------3
		3--------                       --5
		or
		5--------                       --3
		3------        -------------------5
		**********************************/
		
		if( Fct_Frame(SeqA,MIN(Loc5_5_3,Loc5_3_5))==
			Fct_Frame(SeqB,MIN(Loc3_5_3,Loc3_3_5)))
				return(TRUE);
		else
			return(FALSE);
		}
	else if((overhang_A > overhang_B) && Sign(overhang_A)==Sign(overhang_B))
		{
		/**********************************
		5-------------                   -3
		3-                             ---5
		or
		5-                             ---3
		3-------------                   -5
		**********************************/
		if( Fct_Frame(SeqA,MAX(Loc5_5_3,Loc5_3_5))==
			Fct_Frame(SeqB,MAX(Loc3_5_3,Loc3_3_5)))
				return(TRUE);
		else
			return(FALSE);
		}
	else return(FALSE);	
	}





bool_t Fct_Compatible2( int SeqA, int varA_5_3, int varA_3_5,
						 int SeqB, int varB_5_3, int varB_3_5)
	{
	/****************************************************************************************
	This function check if bases from varA_5_3 to varA_3_5 on Sequence SeqA are compatible
	with bases from varB_5_3 to varB_3_5 on Sequence SeqB
	*****************************************************************************************/

	int vara=TRUE;
	int overhang,i;
	
	/**************************************/
	/* non overlapping overhang allowed ? */
	/**************************************/
	if(Preference.allow_part_overhang==FALSE && ((varA_5_3-varA_3_5)!=(varB_5_3-varB_3_5)))
		{
		return(FALSE);
		}
	/****************/
	/* Enzyme blunt */
	/****************/
	/*
		SeqA			SeqB
	5-----------	----------3
	3-----------    ----------5
	
	*/
	else if((varA_5_3-varA_3_5)==0 &&  (varB_5_3-varB_3_5)==0)
		return(TRUE);
	/****************/
	/*  3' overhang */
	/****************/
	else if((varA_5_3-varA_3_5)>0 &&  (varB_5_3-varB_3_5)>0)
		{
		/**************************************
		* if len_overhang(A)<=len_overhang(B) *
		**************************************/
		/*
			SeqA			SeqB
		5-----------	       --3
		3----------    ----------5
		
		*/

		if((varA_5_3-varA_3_5)<=(varB_5_3-varB_3_5))
			{
			overhang = (varA_5_3-varA_3_5);
			for(i=0;i< overhang ;i++)
				{
				if(	Fct_Identique(	seq[SeqA].sequence[Fct_Pos(SeqA,varA_3_5+i)],
									seq[SeqB].sequence[Fct_Pos(SeqB,varB_3_5+i)])!=TRUE)
					{vara=FALSE;break;}
				}
			return(vara);	
			}

		/**************************************
		* if len_overhang(A) >len_overhang(B) *
		**************************************/
		/*
			SeqA			SeqB
		5-----------	 --------3
		3---           ----------5
		
		*/

		else if((varA_5_3-varA_3_5) > (varB_5_3-varB_3_5))
			{
			overhang = (varB_5_3-varB_3_5);
			for(i=0;i< overhang ;i++)
				{
				if(	Fct_Identique(	seq[SeqA].sequence[Fct_Pos(SeqA,varA_5_3-overhang+i)],
									seq[SeqB].sequence[Fct_Pos(SeqB,varB_3_5+i)])!=TRUE)
					{vara=FALSE;break;}
				}
			return(vara);	
			}
		else return(FALSE);
		}
	
	/****************/
	/*  5' overhang */
	/****************/
	else if((varA_5_3-varA_3_5)<0 &&  (varB_5_3-varB_3_5)<0)
		{
		/**************************************
		* if len_overhang(A)<=len_overhang(B) *
		**************************************/
		/*
			SeqA			SeqB
		5----------     ----------3
		3-----------            --5
		
		*/

		if((varA_3_5-varA_5_3)<=(varB_3_5-varB_5_3))
			{
			overhang = (varA_3_5-varA_5_3);
			for(i=0;i< overhang ;i++)
				{
				if(	Fct_Identique(	seq[SeqA].sequence[Fct_Pos(SeqA,varA_5_3+i)],
									seq[SeqB].sequence[Fct_Pos(SeqB,varB_5_3+i)])!=TRUE)
					{vara=FALSE;break;}
				}
			return(vara);	
			}
		/**************************************
		* if len_overhang(A) >len_overhang(B) *
		***************************************/
		/*
			SeqA			SeqB
		5--              ----------3
		3-----------      ---------5
		
		*/

		else if((varA_3_5-varA_5_3) > (varB_3_5-varB_5_3))
			{
			overhang = (varB_3_5-varB_5_3);
			for(i=0;i< overhang ;i++)
				{
				if(	Fct_Identique(	seq[SeqA].sequence[Fct_Pos(SeqA,varA_3_5-overhang+i)],
									seq[SeqB].sequence[Fct_Pos(SeqB,varB_5_3+i)])!=TRUE)
					{vara=FALSE;break;}
				}
			return(vara);	
			}
		else return(FALSE);
		}
	else return(FALSE);
	}


	/*************************************************************************************/                 
	/*     Sub_Cloning find all the solution to clone INSERT (delimited in the sequence
			by the sites NumSiteI5 and NumSiteI3) into VECTOR (delimited in the sequence
		by the sites NumSiteV5 and NumSiteV3)											 */
	/*************************************************************************************/
	
int Sub_Cloning(void)
	{
	/*****local variables ************************************/
	bool_t IsFind=FALSE;
	int				partialV_5=0,partialV_3=0,partialI_5=0,partialI_3=0;	/* number of partial digestions for each site */
	register int	NumSiteV_5,NumSiteV_3,NumSiteI_5,NumSiteI_3;			/* sites used during the search */
	int				TreatV_5,TreatV_3,TreatI_5,TreatI_3;					/* treatment used during the search */
	STRUCT_STRATEGY Strategy;												/* cloning strategy used */
	/**********************************************************/
	Strategy.type=SUB_CLONING;/* sub cloning strategie */
	if(seq[VECTOR].npb==0)
		{err_printf("No VECTOR sequence defined !.\n");BEEP;return(FALSE);}
	if(seq[INSERT].npb==0)
		{err_printf("No INSERT sequence defined !.\n");BEEP;return(FALSE);}
	if(Preference.side_5==TRUE || Preference.side_3==TRUE)
		{
		if(seq[VECTOR].pos_ATG==FALSE)
			{err_printf("No ATG pos. defined for VECTOR !.\n");BEEP;return(FALSE);}
		if(seq[INSERT].pos_ATG==FALSE)
			{err_printf("No ATG pos. defined for INSERT!.\n");BEEP;return(FALSE);}
		}
	if(Search_Done==FALSE) Cmd_Get_Site();
	if(nbr_sites<=0)
		{
		err_printf("No site was found with %s !!!.",FICHIER_ENZYME);
		return(FALSE);
		}
	Init_Enz();
	/**********************************************************/
	
	
	/* 21/07/98 */
	DeleteAllStgys();
	/* This statement is useful for other functions to know if there is one or two sequence to handle */
	Strategy.Couple[SIDE_5].Test_Trans=TRUE;
	Strategy.Couple[SIDE_3].Test_Trans=TRUE;
	
	
	/*************************************************/
	for(NumSiteV_5=0;NumSiteV_5<nbr_sites;NumSiteV_5++)
	/*************************************************/
		{
		
		
			 
		/* Is Enzyme V5 selected ? (cf InitEnz) */
		if(Enzymes[Sites[NumSiteV_5].NumEnz].select[VECTOR]==FALSE)
			continue;
		/** site V5 is on VECTOR ? **/
		if(Sites[NumSiteV_5].NumSeq != VECTOR)
			continue;
		/** site in box **/
		if(Sites[NumSiteV_5].Loc<seq[VECTOR].var_min || Sites[NumSiteV_5].Loc>seq[VECTOR].var_max)
			continue;
		/** site on VECTOR everywhere in box between V5 and var_max, test partial else **/
		partialV_5=Fct_N_sites(Sites[NumSiteV_5].NumEnz,VECTOR,Sites[NumSiteV_5].Loc-1,seq[VECTOR].var_max+1,FALSE);
		if(partialV_5 > Preference.partial)
			continue;
		/*************************************************/
		for(NumSiteI_5=0;NumSiteI_5<nbr_sites;NumSiteI_5++)
		/*************************************************/
			{
			/* Enzyme selected */
			if(Enzymes[Sites[NumSiteI_5].NumEnz].select[INSERT]==FALSE)
				continue;
			/** site on INSERT **/
			if(Sites[NumSiteI_5].NumSeq != INSERT)
				continue;
			/** this site must be in the box **/
			if(Sites[NumSiteI_5].Loc<seq[INSERT].var_min || Sites[NumSiteI_5].Loc>seq[INSERT].var_max)
				continue;
			/*** Site on INSERT nowhere in [var_min-sigma5][sigma3-var_max] everywhere else **/
			if( Fct_N_sites(Sites[NumSiteI_5].NumEnz,INSERT,seq[INSERT].var_min,seq[INSERT].var_max,TRUE)>0 	
			&&  (partialI_5=Fct_N_sites(Sites[NumSiteI_5].NumEnz,INSERT,seq[INSERT].var_min_int,seq[INSERT].var_max_int,TRUE))>Preference.partial)
				continue;			
			/** Must be before sigma5 % */
			if(Sites[NumSiteI_5].Loc>seq[INSERT].var_min_int)
					continue;
			/*************************************************/
			for(NumSiteV_3=0;NumSiteV_3<nbr_sites;NumSiteV_3++)
			/*************************************************/
				{
				/* Enzyme selected */
				if(Enzymes[Sites[NumSiteV_3].NumEnz].select[VECTOR]==FALSE)
					continue;
				/** site on VECTOR **/
				if(Sites[NumSiteV_3].NumSeq != VECTOR)
					continue;
				
				/** Must be on the right of V5 ONLY if V3 != V5**/
				if(NumSiteV_5!=NumSiteV_3)
					if((Sites[NumSiteV_5].Loc+Enzymes[Sites[NumSiteV_5].NumEnz].taille_site) > Sites[NumSiteV_3].Loc)
						continue;
				/** site on VECTOR everywhere in box from V5 **/
				if(Sites[NumSiteV_3].Loc<Sites[NumSiteV_5].Loc || Sites[NumSiteV_3].Loc>seq[VECTOR].var_max)
					continue;
				/* V3 nowhere else from V5 to V3 **/
				partialV_3 = Fct_N_sites(Sites[NumSiteV_3].NumEnz,VECTOR,Sites[NumSiteV_5].Loc-1,Sites[NumSiteV_3].Loc+1,FALSE);
				if(partialV_3>Preference.partial)
					continue;
				/* V5 nowhere else from V5 to V3 **/
				partialV_5 = Fct_N_sites(Sites[NumSiteV_5].NumEnz,VECTOR,Sites[NumSiteV_5].Loc-1,Sites[NumSiteV_3].Loc+1,FALSE);
				if(partialV_5>Preference.partial)
					continue;
				
				/* partial problem with neoschyzomeres G^AATTC - R^AATTY */
				if(Are_Same_Asymetric(NumSiteV_3,NumSiteV_5)==FALSE)
						{
						partialV_3+=Fct_N_sites(Sites[NumSiteV_3].NumEnz,INSERT,Sites[NumSiteV_5].Loc,Sites[NumSiteV_5].Loc,TRUE);
						partialV_5+=Fct_N_sites(Sites[NumSiteV_5].NumEnz,INSERT,Sites[NumSiteV_3].Loc,Sites[NumSiteV_3].Loc,TRUE);
						}

				
				
				/* check number of partial sites , take care if enz are non-palidromic */
					if(Are_Same_Asymetric(NumSiteV_3,NumSiteV_5)==TRUE)
						{
						/*if it is the same ENZYME that cut at 5' and at 3' (two cuts)
						then discard one partial digestion*/
						partialV_5=0;						
						 }
				/* if there is a partial site, is non-blunt allowed ? */
					if( partialV_5>0 && Preference.partial_only_blunt==TRUE &&
						Fct_type(Enzymes[Sites[NumSiteV_5].NumEnz])!=TYPE_BLUNT)
							continue;						
					if( partialV_3>0 && Preference.partial_only_blunt==TRUE &&
						Fct_type(Enzymes[Sites[NumSiteV_3].NumEnz])!=TYPE_BLUNT)
							continue;
				/* may be there can have 2 enzymes with degenerate sites that are not
				the same (GAATTC and GAANTC) but cutting the same site... So the
				number of partial digestion will be overestimated. This case
				have not been explored. So a (very, very....) few solutions involving 
				partial digestions won't be find */
					if(partialV_3+partialV_5>Preference.partial)
						 	continue;
				/*************************************************/
				for(NumSiteI_3=0;NumSiteI_3<nbr_sites;NumSiteI_3++)
				/*************************************************/
					{
					/* Enzyme selected */
					if(Enzymes[Sites[NumSiteI_3].NumEnz].select[INSERT]==FALSE)
						continue;
					/** site is not I5 */
					if(NumSiteI_3 == NumSiteI_5)
						continue;
					/** site on INSERT **/
					if(Sites[NumSiteI_3].NumSeq != INSERT)
						continue;
					/** this site must be in the box and in seq[INSERT].var_max_int**/
					if(Sites[NumSiteI_3].Loc<(seq[INSERT].var_max_int) || Sites[NumSiteI_3].Loc>seq[INSERT].var_max)
						continue;
					/** Must be on the right of I5 **/
					if( (Sites[NumSiteI_5].Loc+Enzymes[Sites[NumSiteI_5].NumEnz].taille_site) > Sites[NumSiteI_3].Loc)
						continue;
					/** site I3 on INSERT nowhere in box from I5, everywhere else **/
					
					partialI_3=Fct_N_sites(Sites[NumSiteI_3].NumEnz,INSERT,Sites[NumSiteI_5].Loc+1,Sites[NumSiteI_3].Loc-1,TRUE);
					partialI_5=Fct_N_sites(Sites[NumSiteI_5].NumEnz,INSERT,Sites[NumSiteI_5].Loc+1,Sites[NumSiteI_3].Loc-1,TRUE);

					/* partial problem with neoschyzomeres G^AATTC - R^AATTY */
					if(Are_Same_Asymetric(NumSiteI_3,NumSiteI_5)==FALSE)
						{
						partialI_3+=Fct_N_sites(Sites[NumSiteI_3].NumEnz,INSERT,Sites[NumSiteI_5].Loc,Sites[NumSiteI_5].Loc,TRUE);
						partialI_5+=Fct_N_sites(Sites[NumSiteI_5].NumEnz,INSERT,Sites[NumSiteI_3].Loc,Sites[NumSiteI_3].Loc,TRUE);
						}
					
					if(Are_Same_Asymetric(NumSiteI_3,NumSiteI_5)==TRUE)
						{
						/*if it is the same ENZYME that cut at 5' and at 3' (two cuts)
						then discard one partial digestion*/
						partialI_5=0;
						 }
					/* if there is a partial site, is non-blunt allowed ? */
					if( partialI_5>0 && Preference.partial_only_blunt==TRUE &&
						Fct_type(Enzymes[Sites[NumSiteI_5].NumEnz])!=TYPE_BLUNT)
							continue;						
					if( partialI_3>0 && Preference.partial_only_blunt==TRUE &&
						Fct_type(Enzymes[Sites[NumSiteI_3].NumEnz])!=TYPE_BLUNT)
							continue;
					if(partialI_3+partialI_5>Preference.partial)
						 	continue;
/***************************************************************************************************/
/* init the strategy with the sites numbers and the number of partial digestion */
Strategy.Couple[SIDE_5].NumSite[VEC_5]		= NumSiteV_5;
Strategy.Couple[SIDE_5].Partial[VEC_5]		= partialV_5;
Strategy.Couple[SIDE_5].NumSite[INS_5]		= NumSiteI_5;
Strategy.Couple[SIDE_5].Partial[INS_5]		= partialI_5;
Strategy.Couple[SIDE_3].NumSite[VEC_3]		= NumSiteV_3;
Strategy.Couple[SIDE_3].Partial[VEC_3]		= partialV_3;
Strategy.Couple[SIDE_3].NumSite[INS_3]		= NumSiteI_3;
Strategy.Couple[SIDE_3].Partial[INS_3]		= partialI_3;
IsFind=FALSE;
/**************************************/
/*Ligation 5' is OK without treatment */
/**************************************/

if(Fct_Test(NumSiteV_5,NumSiteI_5,NO_TREATMENT,Preference.side_5,&Strategy.Couple[SIDE_5])==TRUE)
	{
	/****************************************/
	/*  Ligation 3' is OK without treatment */
	/****************************************/
	
	if(Fct_Test(NumSiteI_3,NumSiteV_3,NO_TREATMENT,Preference.side_3,&Strategy.Couple[SIDE_3])==TRUE)
		{
		/*----------------------------------------*/
		/* init strategy treatments */
		Strategy.Couple[SIDE_5].Treatment[VEC_5]	= NO_TREATMENT;
		Strategy.Couple[SIDE_5].Treatment[INS_5]	= NO_TREATMENT;
		Strategy.Couple[SIDE_3].Treatment[VEC_3]	= NO_TREATMENT;
		Strategy.Couple[SIDE_3].Treatment[INS_3]	= NO_TREATMENT;
		Strategy.use_CIP = Test_Use_CIP(VECTOR, Strategy);
		/*----------------------------------------*/
		if(Strategy.use_CIP==TRUE && Preference.allow_CIP==FALSE)
			continue;
		/*Pouet ajouté le 21/06/98 */
		
		if(NewStrategy(&Strategy)==NULL)
			{
			err_printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
			INKEY;
			DeleteAllStgys();
			return(FALSE);
			}
		IsFind=TRUE;
		}
	/******************************************/
	/*  else ligation 3' is OK with treatment */
	/******************************************/
	if(Preference.allow_all_sol==TRUE || IsFind==FALSE)
	if(Preference.allow_T4!=FALSE) /* if usage of modifying polymerases is allowed */
	  if(Fct_Test(NumSiteI_3,NumSiteV_3,T4_TREATMENT,Preference.side_3,&Strategy.Couple[SIDE_3])==TRUE)
		{
		/* impossible treatment if I3 is the same enzyme than I5 and I5 is not blunt
		   because there is no treatment on side 5' */
		if(Are_Same_Asymetric(NumSiteI_3,NumSiteI_5)==TRUE)
			if(Fct_type(Enzymes[Sites[NumSiteI_5].NumEnz])!=TYPE_BLUNT)
				continue;
		/* impossible treatment if V3 is the same enzyme than V5 and V5 is not blunt
		   because there is no treatment on side 5' */
		if(Are_Same_Asymetric(NumSiteV_3,NumSiteV_5)==TRUE)
			if(Fct_type(Enzymes[Sites[NumSiteV_5].NumEnz])!=TYPE_BLUNT)
				continue;
		/*----------------------------------------*/
		/* simplifies treatments */
		TreatV_5=NO_TREATMENT;TreatV_3=T4_TREATMENT;
		Treatments(&TreatV_5,&TreatV_3,
			Fct_type(Enzymes[Sites[NumSiteV_5].NumEnz]),
			Fct_type(Enzymes[Sites[NumSiteV_3].NumEnz]));
		
		TreatI_5=NO_TREATMENT;TreatI_3=T4_TREATMENT;
		Treatments(&TreatI_5,&TreatI_3,
			Fct_type(Enzymes[Sites[NumSiteI_5].NumEnz]),
			Fct_type(Enzymes[Sites[NumSiteI_3].NumEnz]));
		/*----------------------------------------*/
		/* init strategy treatments */
		Strategy.Couple[SIDE_5].Treatment[VECTOR]	= TreatV_5;
		Strategy.Couple[SIDE_5].Treatment[INSERT]	= TreatI_5;
		Strategy.Couple[SIDE_3].Treatment[VECTOR]	= TreatV_3;
		Strategy.Couple[SIDE_3].Treatment[INSERT]	= TreatI_3;
		Strategy.use_CIP		= Test_Use_CIP(VECTOR, Strategy);
		/*----------------------------------------*/
		if(Strategy.use_CIP==TRUE && Preference.allow_CIP==FALSE)
			continue;

		/*Pouet ajouté le 21/06/98 */
		
		if(NewStrategy(&Strategy)==NULL)
			{err_printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
			INKEY;
			DeleteAllStgys();
			return(FALSE);
			}
		IsFind=TRUE;
		}
	}/* end ligation 5' */
/*****************************************/
/* Else Ligation 5' is OK with treatment */
/*****************************************/
  if(Preference.allow_all_sol==TRUE || IsFind==FALSE)

  if(Fct_Test(NumSiteV_5,NumSiteI_5,T4_TREATMENT,Preference.side_5,&Strategy.Couple[SIDE_5])==TRUE)
	{
	if(Preference.allow_T4==FALSE) /* if usage of modifying polymerases is not allowed */
		continue;
	/****************************************/
	/*  Ligation 3' is OK without treatment */
	/****************************************/
	
	if(Fct_Test(NumSiteI_3,NumSiteV_3,NO_TREATMENT,Preference.side_3,&Strategy.Couple[SIDE_3])==TRUE)
		{
		/* impossible if I3 is the same enzyme than I5 and I5 is not blunt
		   because there is no treatment on side 3' */
		if(Are_Same_Asymetric(NumSiteI_3,NumSiteI_5)==TRUE)
			if(Fct_type(Enzymes[Sites[NumSiteI_5].NumEnz])!=TYPE_BLUNT)
				continue;
		/* impossible if V5 is the same enzyme than V3 and V3 is not blunt
		   because there is no treatment on side 3' */
		if(Are_Same_Asymetric(NumSiteV_3,NumSiteV_5)==TRUE)
			if(Fct_type(Enzymes[Sites[NumSiteV_5].NumEnz])!=TYPE_BLUNT)
				continue;
		/*----------------------------------------*/
		/* simplifies treatments */
		TreatV_5=T4_TREATMENT;TreatV_3=NO_TREATMENT;
		Treatments(&TreatV_5,&TreatV_3,
			Fct_type(Enzymes[Sites[NumSiteV_5].NumEnz]),
			Fct_type(Enzymes[Sites[NumSiteV_3].NumEnz]));

		
		TreatI_5=T4_TREATMENT;TreatI_3=NO_TREATMENT;
		Treatments(&TreatI_5,&TreatI_3,
			Fct_type(Enzymes[Sites[NumSiteI_5].NumEnz]),
			Fct_type(Enzymes[Sites[NumSiteI_3].NumEnz]));
		/*----------------------------------------*/
		/* init strategy treatments */
		Strategy.Couple[SIDE_5].Treatment[VECTOR]	= TreatV_5;
		Strategy.Couple[SIDE_5].Treatment[INSERT]	= TreatI_5;
		Strategy.Couple[SIDE_3].Treatment[VECTOR]	= TreatV_3;
		Strategy.Couple[SIDE_3].Treatment[INSERT]	= TreatI_3;
		Strategy.use_CIP		= Test_Use_CIP(VECTOR, Strategy);		
		/*----------------------------------------*/
		if(Strategy.use_CIP==TRUE && Preference.allow_CIP==FALSE)
			continue;
		
		
		/*Pouet ajouté le 21/06/98 */
		if(NewStrategy(&Strategy)==NULL)
			{err_printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
			INKEY;
			DeleteAllStgys();
			return(FALSE);}
		IsFind=TRUE;
		}
	/******************************************/
	/*  else ligation 3' is OK with treatment */
	/******************************************/
	if(Preference.allow_all_sol==TRUE  || IsFind==FALSE)
	if(Fct_Test(NumSiteI_3,NumSiteV_3,T4_TREATMENT,Preference.side_3,&Strategy.Couple[SIDE_3])==TRUE)
		{
		/*----------------------------------------*/
		/* simplifies treatments */
		TreatV_5=T4_TREATMENT;TreatV_3=T4_TREATMENT;
		Treatments(&TreatV_5,&TreatV_3,
			Fct_type(Enzymes[Sites[NumSiteV_5].NumEnz]),
			Fct_type(Enzymes[Sites[NumSiteV_3].NumEnz]));
		
		TreatI_5=T4_TREATMENT;TreatI_3=T4_TREATMENT;
		Treatments(&TreatI_5,&TreatI_3,
			Fct_type(Enzymes[Sites[NumSiteI_5].NumEnz]),
			Fct_type(Enzymes[Sites[NumSiteI_3].NumEnz]));
		/*----------------------------------------*/
		/* init strategy treatments */
		Strategy.Couple[SIDE_5].Treatment[VECTOR]	= TreatV_5;
		Strategy.Couple[SIDE_5].Treatment[INSERT]	= TreatI_5;
		Strategy.Couple[SIDE_3].Treatment[VECTOR]	= TreatV_3;
		Strategy.Couple[SIDE_3].Treatment[INSERT]	= TreatI_3;
		Strategy.use_CIP		= Test_Use_CIP(VECTOR, Strategy);		
		/*----------------------------------------*/
		if(Strategy.use_CIP==TRUE && Preference.allow_CIP==FALSE)
			continue;
		
		/*Pouet ajouté le 21/06/98 */
		
		if(NewStrategy(&Strategy)==NULL)
			{err_printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
			INKEY;
			DeleteAllStgys();
			return(FALSE);
			}
		IsFind=TRUE;
		}
	}/* end liagtion 5'*/
/***************************************************************************************************/
					
					}/*for(NumSiteI_3=0;NumSiteI_3<nbr_sites;NumSiteI_3++)*/
				}
			}
		}/* end for(NumSiteV_5=0;NumSiteV_5<nbr_sites;NumSiteV_5++) */
	
		
	if(TotalSolutions==0)
		{
		printf("<H2>No solution was found.</H2><UL>\n");
		if(Preference.allenz_classic==TRUE)
			printf("<LI>May be could you use a complete REBASE file ?");
		if(Preference.partial==0)
			printf("<LI>May be could you allow partial digestions ?");
		if(Preference.allow_CIP==FALSE)
			printf("<LI> May be could you allow usage of C.I.P. ?");
		if(Preference.allow_T4==FALSE)
			printf("<LI> May be could you allow usage of modifying polymerases ?");
		if(Preference.memory==TRUE)
			printf("<LI> May be could you avoid to discard short sites ?");
		if(Preference.allow_part_overhang==FALSE)
			printf("<LI> May be could you allow partial overhangs ligation ?");
		printf("</UL><BR>\n");
		return(FALSE);
		}
	else
		{
		SaveCloningSolutions(FORMAT_HTML);
		return(TRUE);
		}
	}


/*******************************************************************************************/
bool_t SaveCloningSolutions(short mode)
	{
	STRUCT_STRATEGY *Stgy=NULL;
	FILE *out=NULL;
	int i;
	
	
	
	out=(FILE*)stdout;
	
	printf("<HR><CENTER><H1>SubCloning Solutions</H1></CENTER><BR>");
	/*write header */
	printMainHeader(out,mode,TRUE);
	
	
	fprintf(out,"</PRE><UL>\n");
	fprintf(out,"<LI><A HREF=\"#CLONING_SOLUTIONS\">Solutions</A>\n");
	fprintf(out,"<LI><A HREF=\"#ENZYMES_LIST\">Enzymes</A>\n");
	fprintf(out,"<LI><A HREF=\"#MISC\">Miscellanous</A>\n");
	fprintf(out,"</UL><PRE>\n");
			
	/* write title 1: solutions */
	fprintf(out,"<A NAME=\"CLONING_SOLUTIONS\">\n");

		printTitle("Cloning Solutions",1,TRUE);
				SetFirstStgy();
		while((Stgy=GetCurrStgy())!=NULL)
			{
			Show_cloning_solution(out,Stgy,mode);
			Digest(mode,out,VECTOR,Stgy->Couple[SIDE_5].NumSite[VECTOR],Stgy->Couple[SIDE_3].NumSite[VECTOR],Stgy->Couple[SIDE_5].Partial[VECTOR]+Stgy->Couple[SIDE_3].Partial[VECTOR]);
			Digest(mode,out,INSERT,Stgy->Couple[SIDE_5].NumSite[INSERT],Stgy->Couple[SIDE_3].NumSite[INSERT],Stgy->Couple[SIDE_5].Partial[INSERT]+Stgy->Couple[SIDE_3].Partial[INSERT]);
			NextStgy();
			}
		printTitle("",1,FALSE);
	/* write title 1:enzymes */
	fprintf(out,"<A NAME=\"ENZYMES_LIST\">\n");
		printTitle("Enzymes Used",1,TRUE);
		fprintf(out,"<UL>");
		for(i=0;i<nbr_enzyme;i++) Enzymes[i].select[VECTOR]=FALSE;
				SetFirstStgy();
		while((Stgy=GetCurrStgy())!=NULL)
			{
			Enzymes[Sites[Stgy->Couple[SIDE_5].NumSite[VECTOR]].NumEnz].select[VECTOR]=TRUE;
			Enzymes[Sites[Stgy->Couple[SIDE_3].NumSite[VECTOR]].NumEnz].select[VECTOR]=TRUE;
			Enzymes[Sites[Stgy->Couple[SIDE_5].NumSite[INSERT]].NumEnz].select[VECTOR]=TRUE;
			Enzymes[Sites[Stgy->Couple[SIDE_3].NumSite[INSERT]].NumEnz].select[VECTOR]=TRUE;
			NextStgy();
			}
		SetFirstStgy();	
		for(i=0;i<nbr_sites;i++)
			{
			if(Enzymes[Sites[i].NumEnz].select[VECTOR]==TRUE)
				{Cmd_Get_Info(mode,out,i);}
			Enzymes[Sites[i].NumEnz].select[VECTOR]=FALSE;
			}
		fprintf(out,"</UL>");
		printTitle("",1,FALSE);
	
	/* write Footer */
	fprintf(out,"<A NAME=\"MISC\">\n");
	printMainHeader(out,mode,FALSE);

	return(TRUE);
	}


/***********************************************************************/
void printHeader(FILE *out,STRUCT_STRATEGY *Strategy,short mode)
	{
	char s[40];
	int i;
	CountStgy(&i);
	sprintf(s,"Solution %d/%d",i,TotalSolutions);
	printTitle(s,2,TRUE);
	}
/***********************************************************************/

void Show_cloning_solution(FILE *out,STRUCT_STRATEGY *Strategy,short mode)
	{
	
	int vara,i;
	char s[40];
	CountStgy(&i);
	sprintf(s,"Solution %d/%d",i,TotalSolutions);
	printTitle(s,2,TRUE);
	
	
	/**********************************************************/
	/* display strategy for vector */
	SemiSolution(out,VECTOR,Strategy,mode);
	fprintf(out,"\n\n");
	/* display strategy for insert */
	SemiSolution(out,INSERT,Strategy,mode);
	/**********************************************************/
	if(Preference.side_5==TRUE && Preference.side_3==FALSE)
		fprintf(out,"Sites will be in frame ligated in 5'.\n\n");
	if(Preference.side_5==FALSE && Preference.side_3==TRUE)
		fprintf(out,"Sites will be in frame ligated in 3'.\n\n");
	if(Preference.side_5==TRUE && Preference.side_3==TRUE)
		fprintf(out,"Sites will be in frame on both sides.\n\n");
	/**********************************************************/
	/* if CIP used, is there a digestion that can help to direct insert ?*/
	Orient_Insert_After_CIP(Strategy,mode,out);
	/**********************************************************/
	/* if frame dependant is there a stop codon ?*/
	if(Preference.side_5==TRUE || Preference.side_3==TRUE)
		{
		/** Is there a stop codon created by the ligation ? **/
		Detect_Stop( out,mode,Strategy->Couple[SIDE_5].NumSite[VECTOR],
						Strategy->Couple[SIDE_5].NumSite[INSERT],
						&Strategy->Couple[SIDE_5]);
		Detect_Stop( out,mode,Strategy->Couple[SIDE_3].NumSite[INSERT],
						Strategy->Couple[SIDE_3].NumSite[VECTOR],
						&Strategy->Couple[SIDE_3]);
		/** where are the other stop codons ? */
		if(Preference.side_5==TRUE)
			{
			Find_stop_codon(mode,out,Strategy->Couple[SIDE_5].NumSite[VECTOR],1,SIDE_5);
			Find_stop_codon(mode,out,Strategy->Couple[SIDE_5].NumSite[INSERT],Sites[Strategy->Couple[SIDE_3].NumSite[INSERT]].Loc,SIDE_3);
			}
		if(Preference.side_3==TRUE)
			{
			Find_stop_codon(mode,out,Strategy->Couple[SIDE_3].NumSite[INSERT],Sites[Strategy->Couple[SIDE_5].NumSite[INSERT]].Loc,SIDE_5);
			Find_stop_codon(mode,out,Strategy->Couple[SIDE_3].NumSite[VECTOR],seq[VECTOR].npb,SIDE_3);
			}
		fprintf(out,"\n");
		}
	/**********************************************************/
	/* find sites that are localised between the 2 vector site and that could be
	used for post-ligation digestion */
	if(Strategy->Couple[SIDE_5].NumSite[VECTOR] != Strategy->Couple[SIDE_3].NumSite[VECTOR])
		{
		vara=FALSE;
		fprintf(out,"Discard non recombinant molecules by a digestion post-ligation with:<BR><UL>");
		for(i=0;i<nbr_enzyme;i++)
			{
			/*
				at least one site between V5 and V3 (vector)
				and no site out of V5 and V3 (vector )
				and no site between I5 and I3 (insert)
			*/
			if( Fct_N_sites(i,VECTOR,
							Sites[Strategy->Couple[SIDE_5].NumSite[VECTOR]].Loc+1,
							Sites[Strategy->Couple[SIDE_3].NumSite[VECTOR]].Loc-1,
							TRUE) != 0 &&
				Fct_N_sites(i,VECTOR,
							Sites[Strategy->Couple[SIDE_5].NumSite[VECTOR]].Loc,
							Sites[Strategy->Couple[SIDE_3].NumSite[VECTOR]].Loc,
							FALSE)== 0 && 
				Fct_N_sites(i,INSERT,
							Sites[Strategy->Couple[SIDE_5].NumSite[INSERT]].Loc,
							Sites[Strategy->Couple[SIDE_3].NumSite[INSERT]].Loc+Enzymes[Sites[Strategy->Couple[SIDE_3].NumSite[INSERT]].NumEnz].taille_site,
							TRUE)== 0)
				{vara=TRUE;fprintf(out,"<LI>%s",Enzymes[i].nom);}
			}
		if(vara!=TRUE)
			fprintf(out,"<LI>no enzyme was found.\n");
		}
	printf("</UL><P>");
	printTitle("",2,FALSE);
	}




void SemiSolution(FILE *out,int NumSeq,STRUCT_STRATEGY *Stg,short mode)
	{
	/*************************************************************/
	/* This procedure displays a strategy for vector OR insert   */
	/*************************************************************/
	char *string_pol[4]={"?","?","T4 DNA polymerase","Klenow DNA polymerase"};
	/***********************************************************/

	fprintf(out,"  Digest %s %s with ",string_seq[NumSeq],seq[NumSeq].FICHIER_ADN);
	fprintf(out,"<A HREF=\"#TAG%d\">%s</A>",
			Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz,
			Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ].nom);
	
	fprintf(out,"[%s] (%d",
		Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ].site_complet,
		Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].Loc);
	
	/* if site SIDE_3 isn't the same enzyme at SIDE5 (take care of asymetric enz)*/
	if(Are_Same_Asymetric(Stg->Couple[SIDE_5].NumSite[NumSeq],Stg->Couple[SIDE_3].NumSite[NumSeq])==FALSE)
		{
		fprintf(out,") and ");
		
		
		fprintf(out,"<A HREF=\"#TAG%d\">%s</A>",
				Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ,
				Enzymes[ Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ].nom);

		fprintf(out," [%s] (%d).\n",
			Enzymes[ Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ].site_complet,
			Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].Loc);
		}
	else	/* if sites are the same enzyme */
		{
		/* if this isn't the same Loc*/
		if(Stg->Couple[SIDE_3].NumSite[NumSeq]!=Stg->Couple[SIDE_5].NumSite[NumSeq])
			fprintf(out,"-%d).\n",Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].Loc);
		/* if this is the same Loc*/
		else
			fprintf(out,").\n");
		}
	ShowSeq(mode,out,Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]], Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]]);
	/**************************************/
	/* modifying    polymerase treatment */
	/************************************/
	if((Stg->Couple[SIDE_5].Treatment[NumSeq])==DELTA_TREATMENT) /* used in function: find frameshift */
		fprintf(out,"Polymerase treatment may be needed in function of the second enzyme.\n");
	else
	switch(Stg->Couple[SIDE_5].Treatment[NumSeq])
		{
		case(NO_TREATMENT  ):
			{
			switch(Stg->Couple[SIDE_3].Treatment[NumSeq])
				{
				case(NO_TREATMENT  ):fprintf(out,"No polymerase treatment is needed.\n");break;
				case(T4_TREATMENT  ):
					{
					fprintf(out,"Digest first with %s. ",Enzymes[ Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ].nom);
					fprintf(out,"Then treat with %s.",string_pol[Fct_type(Enzymes[ Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ])]);
					fprintf(out,"Finally digest with %s .\n",Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ].nom);
					}break;
				default:break;
				}
			}break;
		case(T4_TREATMENT  ):
			{
			switch(Stg->Couple[SIDE_3].Treatment[NumSeq])
				{
				case(NO_TREATMENT  ):
					{
					fprintf(out,"Digest first with %s .",Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ].nom);
					fprintf(out,"Then treat with %s.",string_pol[Fct_type(Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ])]);
					fprintf(out,"Finally digest with %s .\n",Enzymes[ Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ].nom);
					}break;
				case(T4_TREATMENT  ):
					{
					if(	Fct_type(Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ])==TYPE_5_OVER && 
						Fct_type(Enzymes[ Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ])==TYPE_5_OVER)
							fprintf(out,"Treat with Klenow DNA polymerase.\n");
					else if(	Fct_type(Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ])==TYPE_5_OVER && 
								Fct_type(Enzymes[ Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ])==TYPE_BLUNT)
							fprintf(out," Treat with Klenow polymerase.\n");
					else if(	Fct_type(Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ])==TYPE_BLUNT && 
								Fct_type(Enzymes[ Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ])==TYPE_5_OVER)
							fprintf(out," Treat with Klenow polymerase.\n");

					else
							fprintf(out," Treat with T4 DNA polymerase.\n");
					}break;
				default:
					INKEY;
					break;
				}
			}break;
		default:break;
		}
	  /**********************************************/
	 /* Test to see if vector can be self ligated  */
	/**********************************************/
	if(NumSeq==VECTOR && (Stg->use_CIP)==TRUE)
				fprintf(out,"<BLINK>You will have to dephosphorylate your vector.</BLINK>.");
	
	 /*********************************/	
	/* Test for partials digestions  */
   /*********************************/
	if(Stg->Couple[SIDE_5].Partial[NumSeq] > 0 ||  Stg->Couple[SIDE_3].Partial[NumSeq]>0 )
		{
		/* if it is the same side on side 5 and on side 3 */
		if(Are_Same_Asymetric(Stg->Couple[SIDE_5].NumSite[NumSeq],Stg->Couple[SIDE_3].NumSite[NumSeq])==TRUE)
 				fprintf(out,"<BLINK>Beware : %s [%d partial site%c]</BLINK> .\n",
 					
					Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ].nom,
					MAX(Stg->Couple[SIDE_3].Partial[NumSeq],Stg->Couple[SIDE_5].Partial[NumSeq]),
					(MAX(Stg->Couple[SIDE_3].Partial[NumSeq],Stg->Couple[SIDE_5].Partial[NumSeq])>1?'s':'\0'));
		else
			{
			fprintf(out,"<BLINK>Beware :");
			if(Stg->Couple[SIDE_5].Partial[NumSeq] > 0)
				fprintf(out," %s [%d partial site%c]   ",
				Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ].nom ,
				Stg->Couple[SIDE_5].Partial[NumSeq],
				(Stg->Couple[SIDE_5].Partial[NumSeq]>1?'s':'\0'));
			if(Stg->Couple[SIDE_3].Partial[NumSeq] > 0)
				fprintf(out," %s [%d partial site%c]",
				Enzymes[ Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ].nom ,
				Stg->Couple[SIDE_3].Partial[NumSeq],
				(Stg->Couple[SIDE_3].Partial[NumSeq]>1?'s':'\0'));
			fprintf(out,"</BLINK>.\n");
			}
		}
	fprintf(out,"\n");		
	}





int Fct_BaseShift(int NumSite, int Position)
	{
	/**************************************************************/
	/*  this function returns the virtual
		position of a base after a digestion fill-in and ligation
		see: Function Frameshift
		ex: G^AATTC --> GAATTAATTC
			0 12345     0123456789								*/
	/***********************************************************/

	int gain;

	/* if it is a 5' overhang enzyme */
	if(Fct_type(Enzymes[Sites[NumSite].NumEnz])==TYPE_5_OVER)
		{
		gain= Enzymes[Sites[NumSite].NumEnz].pos3_5 - Enzymes[Sites[NumSite].NumEnz].pos5_3;
		
		if(Position < Sites[NumSite].Loc + Enzymes[Sites[NumSite].NumEnz].pos3_5)
			return(Fct_Pos(Sites[NumSite].NumSeq,Position));
		else if(Position>=(Sites[NumSite].Loc+Enzymes[Sites[NumSite].NumEnz].pos3_5) &&
				Position<=(Sites[NumSite].Loc+Enzymes[Sites[NumSite].NumEnz].taille_site+ gain))
			return(Fct_Pos(Sites[NumSite].NumSeq,Position-gain));
		else if(Position>(Sites[NumSite].Loc+Enzymes[Sites[NumSite].NumEnz].taille_site + gain))
			return(Fct_Pos(Sites[NumSite].NumSeq,Position-gain));
		else return('?');
		}
	/* if it is a 3' overhang enzyme */
	else if(Fct_type(Enzymes[Sites[NumSite].NumEnz])==TYPE_3_OVER)
		{
		gain=Enzymes[Sites[NumSite].NumEnz].pos5_3-Enzymes[Sites[NumSite].NumEnz].pos3_5;
		if(Position < Sites[NumSite].Loc + Enzymes[Sites[NumSite].NumEnz].pos3_5)
			return(Fct_Pos(Sites[NumSite].NumSeq,Position));
		else if(Position>=(Sites[NumSite].Loc+Enzymes[Sites[NumSite].NumEnz].pos3_5) &&
				Position<=(Sites[NumSite].Loc+Enzymes[Sites[NumSite].NumEnz].pos5_3))
			return(Fct_Pos(Sites[NumSite].NumSeq,Position+gain));
		else if(Position>(Sites[NumSite].Loc+Enzymes[Sites[NumSite].NumEnz].pos5_3))
			return(Fct_Pos(Sites[NumSite].NumSeq,Position+gain));
		else return('?');
		}
	else return ('?');
	}

int FrameShift(void)
	{
	/*************************************************************************************/
	/* This function finds site cuted by 3' or 5' overhang enzymes filled-in and ligated */
	/* that will produce a frameshift whithin the sequence */
	/*************************************************************************************/


	int				partial=0;
	register int	NumSite;

	bool_t 		stock5,stock3,Mem_T4;
	STRUCT_STRATEGY Strategy;
/*******************************************/
Strategy.type=FRAME_SHIFT;
if(seq[INSERT].npb==0)
	{err_printf("No INSERT sequence defined !.\n");return(FALSE);}
if(seq[INSERT].pos_ATG==FALSE)
	{
	Get_ATG(INSERT);
	if(seq[INSERT].pos_ATG==FALSE)
		{err_printf("No frame can be defined !.\n");return(FALSE);}
	}
if(Search_Done==FALSE) Cmd_Get_Site();
if(nbr_sites<=0)
	{
	err_printf("No site was found !!! with %s",FICHIER_ENZYME);
	return(FALSE);
	}
/*******************************************/
Mem_T4=Preference.allow_T4; /* this procedure implies to use modifying polymerase, so let's memorise the old */
Preference.allow_T4=TRUE;	/* allow_t4 and set it to TRUE */
stock5=Preference.side_5;	/* Frame will be important in this function so let's */
stock3=Preference.side_3;	/* memorise the old valor and set them to TRUE */
Preference.side_5=Preference.side_3=TRUE;
Strategy.Couple[SIDE_5].Test_Trans=Strategy.Couple[SIDE_3].Test_Trans=FALSE; /* the strategy is only used with INSERT */
/* as a fill-in is recquired, a modifying polymerase will always be used */
Strategy.Couple[SIDE_5].Treatment[VECTOR]	= T4_TREATMENT;
Strategy.Couple[SIDE_3].Treatment[VECTOR]	= T4_TREATMENT;
Strategy.Couple[SIDE_5].Treatment[INSERT]	= T4_TREATMENT;
Strategy.Couple[SIDE_3].Treatment[INSERT]	= T4_TREATMENT;


DeleteAllStgys();
/*******************************************/
for(NumSite=0;NumSite<nbr_sites;NumSite++)
/*******************************************/
	{

	/** site on INSERT **/
	if(Sites[NumSite].NumSeq != INSERT)
		continue;
	/** Discard blunt enzymes */
	if(Fct_type(Enzymes[Sites[NumSite].NumEnz])==TYPE_BLUNT)
		continue;
	/** Site on INSERT must be localized  in the box **/
	if(	Sites[NumSite].Loc<seq[INSERT].var_min || 
		Sites[NumSite].Loc>seq[INSERT].var_max)
		continue;
	/** Site must be in the right tolerance percentage **/
	if( ((int)(100.0-100.0*((double)(seq[INSERT].var_max-Sites[NumSite].Loc)/(double)(seq[INSERT].var_max-seq[INSERT].var_min)))) < Preference.DeltaMin)
		continue;
	if( ((int)(100.0-100.0*((double)(seq[INSERT].var_max-Sites[NumSite].Loc)/(double)(seq[INSERT].var_max-seq[INSERT].var_min)))) > Preference.DeltaMax)
		continue;
	/** manage partial siteq **/
	partial=Fct_N_sites(Sites[NumSite].NumEnz,INSERT,
						Sites[NumSite].Loc-1,
						Sites[NumSite].Loc+1,FALSE);
	if(partial>Preference.partial)
		 	continue;
	/********************************************************************/

	/* init the strategy */
	Strategy.Couple[SIDE_5].NumSite[VECTOR]		= NumSite;
	Strategy.Couple[SIDE_5].Partial[VECTOR]		= partial;
	Strategy.Couple[SIDE_3].NumSite[VECTOR]		= NumSite;
	Strategy.Couple[SIDE_3].Partial[VECTOR]		= partial;
	Strategy.Couple[SIDE_5].NumSite[INSERT]		= NumSite;
	Strategy.Couple[SIDE_5].Partial[INSERT]		= partial;
	Strategy.Couple[SIDE_3].NumSite[INSERT]		= NumSite;
	Strategy.Couple[SIDE_3].Partial[INSERT]		= partial;
	

	
	if(NewStrategy(&Strategy)==NULL)
		{
		err_printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
		DeleteAllStgys();
		return(FALSE);
		}


	}
/***************************************************************************************************/


Preference.side_5=stock5;/* set those valors to original */
Preference.side_3=stock3;
Preference.allow_T4=Mem_T4;


if(TotalSolutions==FALSE)
	{
	printf("<H2>No solution was found</H2><UL>");
	if(Preference.allenz_classic==TRUE)
		printf("<LI>May be could you use a complete REBASE file ?");
	if(Preference.partial==0)
		printf("<LI>May be could you allow partial digestions ?");
	if(Preference.memory==TRUE)
		printf("<LI>May be could you avoid to discard short sites ?");
	printf("</UL><BR>");
	return(FALSE);
	}
else 
	{
	SaveShiftSolutions(FORMAT_HTML);
	return(TRUE);
	}
}





/*******************************************************************/
void SolutionShift(short mode, FILE *out,STRUCT_STRATEGY *Strategy)
	{
	int 			i,j,varb=0,val_shift=0,overhang;
	int				partial=0,site_reconstitued=0;
	int				NumSite;
	double 			pctage=0.0;
	char c1;
	div_t	r;
	char s[40];
	/********************************************************************/
	CountStgy(&i);
	fprintf(out,"<A NAME=\"SOL_%d\">",i);
	sprintf(s,"Solution %d/%d",i,TotalSolutions);
	printTitle(s,2,TRUE);
		
	NumSite=Strategy->Couple[SIDE_5].NumSite[VECTOR];
	partial=Strategy->Couple[SIDE_5].Partial[VECTOR];

	varb =(int)(((double)(Sites[NumSite].Loc-seq[INSERT].var_min)/(double)(seq[INSERT].var_max-seq[INSERT].var_min))*(double)GRAPHIC);
	r = div(abs(Enzymes[Sites[NumSite].NumEnz].pos5_3-Enzymes[Sites[NumSite].NumEnz].pos3_5),3);
	val_shift=r.rem;
	pctage=((100.0*((double)(seq[INSERT].var_max-Sites[NumSite].Loc)/(double)(seq[INSERT].var_max-seq[INSERT].var_min))));
	overhang=abs(Enzymes[Sites[NumSite].NumEnz].pos5_3-Enzymes[Sites[NumSite].NumEnz].pos3_5);

	/* check if site is reconstitued after ligation */
	for(i=Sites[NumSite].Loc;i<=Sites[NumSite].Loc+ Enzymes[Sites[NumSite].NumEnz].taille_site+MAX(Enzymes[Sites[NumSite].NumEnz].pos5_3,Enzymes[Sites[NumSite].NumEnz].pos3_5);i++)
		{
		site_reconstitued=TRUE;
		for(j=0;j<Enzymes[Sites[NumSite].NumEnz].taille_site;j++)
			if(Fct_Identique(Fct_BaseShift(NumSite,Fct_Pos(INSERT,i+j)),Enzymes[Sites[NumSite].NumEnz].site[j])!=TRUE)
				{
				site_reconstitued=FALSE;
				break;
				}
		if(site_reconstitued==TRUE) break;
		}
	
	
	
	fprintf(out,"\n%s has found an Enzyme that could induce frameshift.\n\n",VAR_VERSION);

	/* display sequence */
	SemiSolution(out,INSERT,Strategy,mode);
	/* display sequence after fill-in and liagation */
	fprintf(out,"\nAfter digestion, fill-in and ligation:\n\n");
	/* write the direct strand */
	fprintf(out,"<FONT COLOR=\"#%s\">",RED_COLOR);

	fprintf(out,"  5'  ");
	for(i= Sites[NumSite].Loc - 4; i<= Sites[NumSite].Loc +PALINDROME_MAX_ENZYME;i++)
		{
		if (Fct_Frame(INSERT,i)==IS_IN_FRAME) fprintf(out,".");
		fprintf(out,"%c",seq[INSERT].sequence[Fct_BaseShift(NumSite,Fct_Pos(INSERT,i))]);
		}
	/* write the anti strand */

	fprintf(out," 3'</FONT>\n  3'  <FONT COLOR=\"#%s\">",RED_COLOR);
	
	for(i= Sites[NumSite].Loc - 4; i<= Sites[NumSite].Loc +PALINDROME_MAX_ENZYME;i++)
		{
		if (Fct_Frame(INSERT,i)==IS_IN_FRAME) fprintf(out,".");
		fprintf(out,"%c",Fct_Complementaire(seq[INSERT].sequence[Fct_BaseShift(NumSite,Fct_Pos(INSERT,i))]));
		}
	fprintf(out,"</FONT>");
	fprintf(out," 5'\n");
	/* write the proteic sequence */
	if(seq[INSERT].pos_ATG!=FALSE)
		{			
		fprintf(out,"  NH2 ");
		fprintf(out,"<FONT COLOR=\"#%s\">",YELLOW_COLOR);
		for(i= Sites[NumSite].Loc - 4; i<= Sites[NumSite].Loc +PALINDROME_MAX_ENZYME;i++)
			if(Fct_Frame(INSERT,i)==IS_IN_FRAME)
				{
				fprintf(out," %c",(Fct_Traduction(	seq[INSERT].sequence[Fct_BaseShift(NumSite,i)],
									seq[INSERT].sequence[Fct_BaseShift(NumSite,i+1)],
									seq[INSERT].sequence[Fct_BaseShift(NumSite,i+2)])));
				}
			else
				fprintf(out," ");
		fprintf(out,"  \tCOOH\n");
		fprintf(out,"</FONT>");
		}
	/********************************************************************/
	fprintf(out,"<UL>");
	/* write parameters */
	if(Enzymes[Sites[NumSite].NumEnz].pos3_5-Enzymes[Sites[NumSite].NumEnz].pos5_3>=0)
		{
		fprintf(out,"<LI> FrameShift (+ %d)\n",val_shift);
		fprintf(out,"<LI> %d base%c Added.\n",
		(Enzymes[Sites[NumSite].NumEnz].pos3_5-Enzymes[Sites[NumSite].NumEnz].pos5_3),
		((Enzymes[Sites[NumSite].NumEnz].pos3_5-Enzymes[Sites[NumSite].NumEnz].pos5_3)>1?'s':'\0'));
		}
	else
		{
		fprintf(out,"<LI> FrameShift (- %d)\n",val_shift);
		fprintf(out,"<LI> %d base%c Deleted.\n",
		(Enzymes[Sites[NumSite].NumEnz].taille_site-(Enzymes[Sites[NumSite].NumEnz].pos5_3-Enzymes[Sites[NumSite].NumEnz].pos3_5)),
		((Enzymes[Sites[NumSite].NumEnz].taille_site-(Enzymes[Sites[NumSite].NumEnz].pos5_3-Enzymes[Sites[NumSite].NumEnz].pos3_5))>1?'s':'\0'));
		}
	fprintf(out,"<LI> Site is %sreconstitued after ligation.\n",((site_reconstitued==TRUE)?"":"NOT "));
	fprintf(out,"<LI> [%d %%] percentage of Insert.\n\n",(int)pctage);
	fprintf(out,"</UL>");
	/********************************************************************/
	/* display a figure */
	fprintf(out,"\n\tFIGURE:\n\t\t");
	for(i=0;i<=varb;i++) fprintf(out,"=");
	fprintf(out,"\n\t\t");
	for(i=0;i<varb;i++) fprintf(out," ");
	fprintf(out,"| (%c%d)\n\t\t",((Enzymes[Sites[NumSite].NumEnz].pos3_5-Enzymes[Sites[NumSite].NumEnz].pos5_3>=0)?'+':'-'),val_shift);
	for(i=0;i<varb;i++) fprintf(out," ");
	for(i=varb;i<=GRAPHIC;i++) fprintf(out,"=");
	fprintf(out,"\n\n");
	/********************************************************************/
	fprintf(out,"Translated sequence:...\n");
	/* Scan all the fragment of INSERT */
	j=0;
	for(i=seq[INSERT].var_min;i<=seq[INSERT].var_max;i++)
		{
		if(Fct_Frame(INSERT,i)==IS_IN_FRAME)
			{
			c1=Fct_Traduction(	seq[INSERT].sequence[Fct_Pos(INSERT,Fct_BaseShift(NumSite,i  ))],
								seq[INSERT].sequence[Fct_Pos(INSERT,Fct_BaseShift(NumSite,i+1))],
								seq[INSERT].sequence[Fct_Pos(INSERT,Fct_BaseShift(NumSite,i+2))]);
			if(i>= Sites[NumSite].Loc && i<=(Sites[NumSite].Loc+Enzymes[Sites[NumSite].NumEnz].taille_site))
				c1=LOWER(c1);
			fprintf(out,"%c",c1);
			if(++j>LARGEUR_ECRAN-2) {j=0;fprintf(out,"\n");}
			/* if there is a stop after the site, break the loop */
			if(Is_Codon_Stop(c1)==TRUE && i>=Sites[NumSite].Loc) 
				break;
			}
		}
	fprintf(out,"...\n\n");
	/********************************************************************/
	printTitle(s,2,FALSE);
	}




  /*=================================================================*/
 /*========= show all alignments of truncated protein           ====*/
/*=================================================================*/
	
	
	
void Show_All_Alignment(FILE *out,short mode)
	{
	int m=0,stop_codon=FALSE,var5_5_3,var3_5_3;
	char c1;
	int pos,r, partialI_5,partialI_3; 
	register int NumSiteI_5,NumSiteI_3;
	STRUCT_STRATEGY *save=NULL,*Stgy=NULL;
	bool_t _TREATMENT, it_is_a_C_term_del ;
	
	save=GetCurrStgy();
	
	fprintf(out,"<UL>");
	
	
		/* Scan all the fragment of INSERT */
	for(r=seq[INSERT].var_min;r<=seq[INSERT].var_max;r=r+(LARGEUR_ECRAN+LARGEUR_ECRAN))
		{
		/*fprintf(out,"r=%d<->%d\n",r,r+(LARGEUR_ECRAN+LARGEUR_ECRAN));*/
		fprintf(out,"<LI>");
		SetFirstStgy();
		while(GetCurrStgy()!=NULL)
			{
			Stgy=CurrStgy;
			NumSiteI_5=Stgy->Couple[SIDE_5].NumSite[VECTOR];
			partialI_5=Stgy->Couple[SIDE_5].Partial[VECTOR];
			NumSiteI_3=Stgy->Couple[SIDE_3].NumSite[INSERT];
			partialI_3=Stgy->Couple[SIDE_3].Partial[INSERT];
			_TREATMENT=(Stgy->type==DELETION_FRAME_NO?NO_TREATMENT:T4_TREATMENT);
			it_is_a_C_term_del=(Stgy->type==DELETION_CARBOXY?TRUE:FALSE);
			/* adjust length extremities in function of treatment */
						
			
			if(NumSiteI_5>=nbr_sites || NumSiteI_3>=nbr_sites || Sites[NumSiteI_3].NumEnz>=nbr_enzyme || Sites[NumSiteI_5].NumEnz>=nbr_enzyme)
				while(1) printf("Erreur !!!\n");
			
			
			var5_5_3=(_TREATMENT==NO_TREATMENT ? Enzymes[Sites[NumSiteI_5].NumEnz].pos5_3: Enzymes[Sites[NumSiteI_5].NumEnz].pos3_5);
			var3_5_3=(_TREATMENT==NO_TREATMENT ? Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3: Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3);
			CountStgy(&pos);
			fprintf(out,"<LI><A HREF=\"#SOL_%d\">",pos);
			fprintf(out,"%03d:",pos);
			fprintf(out,"</A>");
			stop_codon=FALSE;
			for(m=seq[INSERT].var_min;m<=seq[INSERT].var_max;m++)
				{
				if(Fct_Frame(INSERT,m)==IS_IN_FRAME)
					{
					/* c1 will be the amino acid displayed */
					/* if pos m is just at the level of the junction of the ligation */
					if(m==(Sites[NumSiteI_5].Loc+var5_5_3-1))
						c1=Fct_Traduction(seq[INSERT].sequence[Sites[NumSiteI_5].Loc+var5_5_3-1],
										  seq[INSERT].sequence[Sites[NumSiteI_3].Loc+var3_5_3],
										  seq[INSERT].sequence[Sites[NumSiteI_3].Loc+var3_5_3+1]);
					else if(m==(Sites[NumSiteI_5].Loc+var5_5_3-2))
						c1=Fct_Traduction(seq[INSERT].sequence[Sites[NumSiteI_5].Loc+var5_5_3-2],
										  seq[INSERT].sequence[Sites[NumSiteI_5].Loc+var5_5_3-1],
										  seq[INSERT].sequence[Sites[NumSiteI_3].Loc+var3_5_3]);
					/* else if m is BEFORE the junction */
					else if (m <= (Sites[NumSiteI_5].Loc+var5_5_3-1) )
						c1=Fct_Traduction(seq[INSERT].sequence[m],seq[INSERT].sequence[m+1],seq[INSERT].sequence[m+2]);
					/* else if m is AFTER the junction */
					else if (m >= (Sites[NumSiteI_3].Loc+var3_5_3)   )
						c1=Fct_Traduction(seq[INSERT].sequence[m],seq[INSERT].sequence[m+1],seq[INSERT].sequence[m+2]);
					else c1='-';
					
					if(it_is_a_C_term_del==TRUE && m > (Sites[NumSiteI_5].Loc))
						c1='-';
					/* if there is a stop, display in lower case */
					if(Is_Codon_Stop(c1)==TRUE) stop_codon=TRUE;
					
					if(m<=r) continue;
					if(m>=r+(LARGEUR_ECRAN+LARGEUR_ECRAN)) break;
					
					c1=(stop_codon==FALSE?c1:LOWER(c1));
					fprintf(out,"%c",c1);
					
					}
				
				}
			if(m>=seq[INSERT].var_max &&  stop_codon==FALSE)
				fprintf(out,"[stop]...");
			NextStgy();
			fprintf(out,"\n");
			}
		}
	fprintf(out,"</UL>");
	
	SetStgy(save);
	fprintf(out,"\n");
	
	}
	


  /*=================================================================*/
 /*========= save alignment of truncated protein           =========*/
/*=================================================================*/
	
void Save_alignment( int NumSiteI_5, int partialI_5,
					int NumSiteI_3,int partialI_3, 
					bool_t _TREATMENT,bool_t it_is_a_C_term_del ,FILE *out, short mode)
	{

	int m=0,stop_codon=FALSE,var5_5_3,var3_5_3;
	char c1;
	int i=0;
	
	/* adjust length extremities in function of treatment */
	var5_5_3=(_TREATMENT==NO_TREATMENT ? Enzymes[Sites[NumSiteI_5].NumEnz].pos5_3: Enzymes[Sites[NumSiteI_5].NumEnz].pos3_5);
	var3_5_3=(_TREATMENT==NO_TREATMENT ? Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3: Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3);
	
	fprintf(out,"<DL><DT>Translated truncated sequence:<DD>");
	

	/* Scan all the fragment of INSERT */
	for(m=seq[INSERT].var_min;m<=seq[INSERT].var_max;m++)
		{
		if(m>=seq[INSERT].var_max &&  stop_codon==FALSE)
			fprintf(out,"[stop]...\n");
		else if(Fct_Frame(INSERT,m)==IS_IN_FRAME)
			{
			/* c1 will be the amino acid displayed */
			/* if pos m is just at the level of the junction of the ligation */
			if(m==(Sites[NumSiteI_5].Loc+var5_5_3-1))
				c1=Fct_Traduction(seq[INSERT].sequence[Sites[NumSiteI_5].Loc+var5_5_3-1],
								  seq[INSERT].sequence[Sites[NumSiteI_3].Loc+var3_5_3],
								  seq[INSERT].sequence[Sites[NumSiteI_3].Loc+var3_5_3+1]);
			else if(m==(Sites[NumSiteI_5].Loc+var5_5_3-2))
				c1=Fct_Traduction(seq[INSERT].sequence[Sites[NumSiteI_5].Loc+var5_5_3-2],
								  seq[INSERT].sequence[Sites[NumSiteI_5].Loc+var5_5_3-1],
								  seq[INSERT].sequence[Sites[NumSiteI_3].Loc+var3_5_3]);
			/* else if m is BEFORE the junction */
			else if (m <= (Sites[NumSiteI_5].Loc+var5_5_3-1) )
				c1=Fct_Traduction(seq[INSERT].sequence[m],seq[INSERT].sequence[m+1],seq[INSERT].sequence[m+2]);
			/* else if m is AFTER the junction */
			else if (m >= (Sites[NumSiteI_3].Loc+var3_5_3)   )
				c1=Fct_Traduction(seq[INSERT].sequence[m],seq[INSERT].sequence[m+1],seq[INSERT].sequence[m+2]);
			else c1='-';
			
			if(it_is_a_C_term_del==TRUE && m > (Sites[NumSiteI_5].Loc))
				c1='-';
			
			/* if there is a stop, display in lower case */
			if(Is_Codon_Stop(c1)==TRUE) stop_codon=TRUE;
			c1=(stop_codon==FALSE?c1:LOWER(c1));
			fprintf(out,"%c",c1);
			if(++i>=LARGEUR_ECRAN-1)
				{i=0;fprintf(out,"\n");}
			}
		}
	fprintf(out,"</DD></DL>");
	}



/********************************************************************************/
	/*____________________________________________________________
	    This procedure look at in frame deletion and
		carboxy terminal deletion among all the sites in INSERT  
	____________________________________________________________*/
	
int DeltaFrame(void)
	{
	int 			i,vara;
	double			pctage=0.0;
	int				partialI_5=0,partialI_3=0;
	register int	NumSiteI_5,NumSiteI_3;
	register int	Loop_I_5,Loop_I_3;
	bool_t 		stock5,stock3,IsFind=FALSE;

	STRUCT_STRATEGY Strategy;
	/****************************************************/
	Strategy.type=DELETION_FRAME_NO;
	if(seq[INSERT].npb==0)
		{err_printf("No INSERT sequence defined !.");BEEP;return(FALSE);}
	if(seq[INSERT].pos_ATG==FALSE)
		{
		Get_ATG(INSERT);
		if(seq[INSERT].pos_ATG==FALSE)
			{err_printf("No frame can be defined !.");BEEP;return(FALSE);}
		}
	if(Search_Done==FALSE) Cmd_Get_Site();
	if(nbr_sites<=0)
		{
		err_printf("No site was found with %s!!!.",FICHIER_ENZYME);
		return(FALSE);
		}

	/****************************************************/
	DeleteAllStgys();
	stock5=Preference.side_5;   /* Frame will be important in this function so let's */
	stock3=Preference.side_3;	/* memorise the old valor and set them to TRUE */
	Preference.side_5=Preference.side_3=TRUE;
	Strategy.Couple[SIDE_5].Test_Trans=Strategy.Couple[SIDE_3].Test_Trans=FALSE;
/****************************************************/
for(Loop_I_5=0;Loop_I_5<nbr_sites;Loop_I_5++)
/****************************************************/
	{

	
	/** site on INSERT **/
	if(Sites[Loop_I_5].NumSeq != INSERT)
		continue;
	/** Site on INSERT must be localized  in the box **/
	if(	Sites[Loop_I_5].Loc<seq[INSERT].var_min || 
		Sites[Loop_I_5].Loc>seq[INSERT].var_max)
		continue;
	/***************************************************************/
	for(Loop_I_3=Loop_I_5+1;Loop_I_3<nbr_sites;Loop_I_3++)
	/***************************************************************/
		{
		  /********************/
		 /** site on INSERT **/
		/********************/
		if(Sites[Loop_I_3].NumSeq != INSERT)
			continue;
		/* site in cloning box */
		if(	Sites[Loop_I_3].Loc<seq[INSERT].var_min || 
			Sites[Loop_I_3].Loc>seq[INSERT].var_max)
				continue;
		/********************************/
		/** Must be on the right of I5 **/
		/********************************/
		/* orient the two sites defined by NumSite */
		if(Sites[Loop_I_5].Loc<Sites[Loop_I_3].Loc)
			{NumSiteI_5=Loop_I_5;NumSiteI_3=Loop_I_3;}
		else
			{NumSiteI_5=Loop_I_3;NumSiteI_3=Loop_I_5;}
		/* Numsite5 is not OVER Numsite3 */
		if( (Sites[NumSiteI_5].Loc+Enzymes[Sites[NumSiteI_5].NumEnz].taille_site) > Sites[NumSiteI_3].Loc)
			 continue;
				
		pctage=(double)((double)(Sites[NumSiteI_3].Loc-Sites[NumSiteI_5].Loc)/(double)(seq[INSERT].var_max-seq[INSERT].var_min))*100.0;
		if(pctage>(double)Preference.DeltaMax || pctage<(double)Preference.DeltaMin)
			continue;
		  /************************/
		 /** partial digestions **/
		/************************/
		partialI_3=Fct_N_sites(Sites[NumSiteI_3].NumEnz,INSERT,Sites[NumSiteI_5].Loc-1,Sites[NumSiteI_3].Loc+1,FALSE);
		partialI_5=Fct_N_sites(Sites[NumSiteI_5].NumEnz,INSERT,Sites[NumSiteI_5].Loc-1,Sites[NumSiteI_3].Loc+1,FALSE);
		
		/* partial problem with neoschyzomeres G^AATTC - R^AATTY */
		if(Are_Same_Asymetric(NumSiteI_3,NumSiteI_5)==FALSE)
			{
			partialI_3+=Fct_N_sites(Sites[NumSiteI_3].NumEnz,INSERT,Sites[NumSiteI_5].Loc,Sites[NumSiteI_5].Loc,TRUE);
			partialI_5+=Fct_N_sites(Sites[NumSiteI_5].NumEnz,INSERT,Sites[NumSiteI_3].Loc,Sites[NumSiteI_3].Loc,TRUE);
			}

		
		if(Are_Same_Asymetric(NumSiteI_3,NumSiteI_5)==TRUE)
			{
			partialI_5=0;
			 }
		if(partialI_3+partialI_5>Preference.partial)
			 	continue;
		/* if there is a partial site, is it blunt ? */
		if( partialI_5>0 && Preference.partial_only_blunt==TRUE && Fct_type(Enzymes[Sites[NumSiteI_5].NumEnz])!=TYPE_BLUNT)
				continue;						
		if( partialI_3>0 && Preference.partial_only_blunt==TRUE && Fct_type(Enzymes[Sites[NumSiteI_3].NumEnz])!=TYPE_BLUNT)
				continue;
		IsFind=FALSE;
		/********************************************************/
		Strategy.Couple[SIDE_5].NumSite[VECTOR]		= NumSiteI_5;
		Strategy.Couple[SIDE_5].Treatment[VECTOR]	= NO_TREATMENT;
		Strategy.Couple[SIDE_5].Partial[VECTOR]		= partialI_5;
		
		Strategy.Couple[SIDE_3].NumSite[VECTOR]		= NumSiteI_5;
		Strategy.Couple[SIDE_3].Treatment[VECTOR]	= NO_TREATMENT;
		Strategy.Couple[SIDE_3].Partial[VECTOR]		= partialI_5;
		
		Strategy.Couple[SIDE_5].NumSite[INSERT]		= NumSiteI_3;
		Strategy.Couple[SIDE_5].Treatment[INSERT]	= NO_TREATMENT;
		Strategy.Couple[SIDE_5].Partial[INSERT]		= partialI_3;
		
		Strategy.Couple[SIDE_3].NumSite[INSERT]		= NumSiteI_3;
		Strategy.Couple[SIDE_3].Treatment[INSERT]	= NO_TREATMENT;
		Strategy.Couple[SIDE_3].Partial[INSERT]		= partialI_3;
		/* DELETION_FRAME strategie */
		/* are I_5 and I_3 compatible without treatment ?*/
		if(Fct_Test(NumSiteI_5,NumSiteI_3,NO_TREATMENT,TRUE,&Strategy.Couple[SIDE_5])==TRUE)
			{
			Strategy.type=DELETION_FRAME_NO;

			if(NewStrategy(&Strategy)==NULL)
					{
					printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
					INKEY;
					DeleteAllStgys();
					return(FALSE);}
			IsFind=TRUE;
			}
		/* are I_5 and I_3 compatible WITH treatment ?*/
		if(Preference.allow_all_sol==TRUE || IsFind==FALSE)
		if(Preference.allow_T4!=FALSE)
			{
			Strategy.Couple[SIDE_5].Treatment[VECTOR]	= T4_TREATMENT;
			Strategy.Couple[SIDE_3].Treatment[VECTOR]	= T4_TREATMENT;
			Strategy.Couple[SIDE_5].Treatment[INSERT]	= T4_TREATMENT;
			Strategy.Couple[SIDE_3].Treatment[INSERT]	= T4_TREATMENT;
			Strategy.type=DELETION_FRAME_T4;
			if(Fct_Test(NumSiteI_5,NumSiteI_3,T4_TREATMENT,TRUE,&Strategy.Couple[SIDE_5])==TRUE)
				{
				if(NewStrategy(&Strategy)==NULL)
					{
					printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
					INKEY;
					DeleteAllStgys();
					return(FALSE);
					}
				}
			}
		}
	
/**** DELETIONS C TERMINAL **************************************************/
	if(Preference.search_C_term==TRUE)
		{
		/* check frameshift length */
		NumSiteI_5=Loop_I_5;
		pctage=(double)((double)(seq[INSERT].var_max-Sites[NumSiteI_5].Loc)/(double)(seq[INSERT].var_max-seq[INSERT].var_min))*100.0;
		if(pctage>(double)Preference.DeltaMax || pctage<(double)Preference.DeltaMin)
				continue;
		/** check partial site**/
		partialI_5=Fct_N_sites(Sites[NumSiteI_5].NumEnz,INSERT,Sites[NumSiteI_5].Loc-1,seq[INSERT].var_max+1,FALSE);
		if(partialI_5>Preference.partial)
			 continue;
		/* if there is a partial site, is it blunt ? */
		if( partialI_5>0 && Preference.partial_only_blunt==TRUE && Fct_type(Enzymes[Sites[NumSiteI_5].NumEnz])!=TYPE_BLUNT)
				continue;
		/* is there a site between I5 and var_max  that could be used as a second enzyme ? **********/
		vara=FALSE; /* no enzyme found */
		for(i=0;i<nbr_sites;i++)
			{

			if(Preference.allow_T4==FALSE)
				if(Fct_Test(NumSiteI_5,i,NO_TREATMENT,TRUE,&Strategy.Couple[SIDE_5])==FALSE)
					continue;
			if(Fct_N_sites(Sites[i].NumEnz,INSERT,
							Sites[NumSiteI_5].Loc+Enzymes[Sites[NumSiteI_5].NumEnz].taille_site,
							seq[INSERT].var_max,
							TRUE)>0 &&
			  Fct_N_sites(Sites[i].NumEnz,INSERT,
							Sites[NumSiteI_5].Loc-1,
							seq[INSERT].var_max+1,
							FALSE)==0)
					{
					vara=TRUE; /* 1 enzyme found */
					break;
					}
			}
		if(vara==FALSE)
			continue;
		/****************************************************************************/
		/* define strategy and display it */
		Strategy.Couple[SIDE_5].NumSite[VECTOR]		= NumSiteI_5;
		Strategy.Couple[SIDE_5].Treatment[VECTOR]	= NO_TREATMENT;
		Strategy.Couple[SIDE_5].Partial[VECTOR]		= 0;
		Strategy.Couple[SIDE_3].NumSite[VECTOR]		= NumSiteI_5;
		Strategy.Couple[SIDE_3].Treatment[VECTOR]	= NO_TREATMENT;
		Strategy.Couple[SIDE_3].Partial[VECTOR]		= 0;
		Strategy.Couple[SIDE_5].NumSite[INSERT]		= NumSiteI_5;
		Strategy.Couple[SIDE_5].Treatment[INSERT]	= NO_TREATMENT;
		Strategy.Couple[SIDE_5].Partial[INSERT]		= partialI_5;
		Strategy.Couple[SIDE_3].NumSite[INSERT]		= NumSiteI_5;
		Strategy.Couple[SIDE_3].Treatment[INSERT]	= NO_TREATMENT;
		Strategy.Couple[SIDE_3].Partial[INSERT]		= partialI_5;
		Strategy.type=DELETION_CARBOXY;
		

		if(NewStrategy(&Strategy)==NULL)
				{
				err_printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
				INKEY;
				DeleteAllStgys();
				return(FALSE);
				}
		}
	}
/*********************************************************************************/
Preference.side_5=stock5; /* set those valors to original */
Preference.side_3=stock3;



if(TotalSolutions==FALSE)
	{
	printf("<H2>No solution was found</H2><UL>");
	if(Preference.allenz_classic==TRUE)
		printf("<LI>May be could you use a complete REBASE file ?");
	if(Preference.partial==0)
		printf("<LI>May be could you allow partial digestions ?");
	if(Preference.allow_T4==FALSE)
		printf("<LI>May be could you allow usage of modifying polymerases ?");
	if(Preference.memory==TRUE)
		printf("<LI>May be could you avoid to discard short sites ?");
	if(Preference.allow_part_overhang==FALSE)
		printf("<LI>May be could you allow partial overhangs ligation ?");
	printf("</UL><BR>");
	return(FALSE);
	}
else
	{
	SaveDeltaSolutions(FORMAT_HTML);
	return(TRUE);
	}
}





  /*********************************/
 /* Display a C terminal deletion */
/*********************************/
int Solution_C_Terminal(short mode,FILE *out,int NumSiteI_5,STRUCT_STRATEGY *Strategy)
	{
	/*******************************/
	int i,vara;
	double varr,vard,varb,varc;
	char s[40];
	/*******************************/
	CountStgy(&i);
	fprintf(out,"<A NAME=\"SOL_%d\">",i);
	sprintf(s,"Solution %d/%d",i,TotalSolutions);
	printTitle(s,2,TRUE);
	fprintf(out,"Those sites that do not necesseraly create in frame deletion but they can be used to make Carboxy-terminal deletions.\n");
							
	/* init strategy in order to use the function Solution_2 */
	Strategy->Couple[SIDE_5].Test_Trans=Strategy->Couple[SIDE_3].Test_Trans=FALSE;
	Strategy->Couple[SIDE_5].Treatment[VECTOR]	= DELTA_TREATMENT;
	Strategy->Couple[SIDE_3].Treatment[VECTOR]	= DELTA_TREATMENT;
	Strategy->Couple[SIDE_5].Treatment[INSERT]	= DELTA_TREATMENT;
	Strategy->Couple[SIDE_3].Treatment[INSERT]	= DELTA_TREATMENT;
	SemiSolution(out,INSERT,Strategy,mode);
	/* display figure */
	varr = (double)(seq[INSERT].var_max-seq[INSERT].var_min)/(double)GRAPHIC;
	if(varr==0.0) {printf("vara=0\n");BEEP;BEEP;BEEP;INKEY;  ERROR_USER;}
	vard = (double)(Sites[NumSiteI_5].Loc-seq[INSERT].var_min)/(varr==0?1:varr);
	varb = ((double)(seq[INSERT].var_max-seq[INSERT].var_min)/(varr==0?1:varr))-(double)vard;
	varc = (double)GRAPHIC-vard-varb;
	fprintf(out,"BEWARE: Check if there is a STOP CODON after ligation.\n");
	fprintf(out,"Cloning box boundaries :[%d-%d] [%d-%d].\n",seq[INSERT].var_min,seq[INSERT].var_min_int,seq[INSERT].var_max_int,seq[INSERT].var_max);
	fprintf(out,"\n\tOriginal: 5' ");
	for(i=0;i<=((int)vard+(int)varb+(int)varc);i++)
				fprintf(out,"=");
	fprintf(out," 3'\n\tDeletion: 5' ");
	for(i=0;i<(int)vard;i++) fprintf(out,"=");
	for(i=0;i<(int)varb;i++) fprintf(out,".");
	for(i=0;i<=(int)varc;i++) fprintf(out,".");
	fprintf(out," 3'\n");
	vara=seq[INSERT].var_max-Sites[NumSiteI_5].Loc;	
	fprintf(out,"\nEquivalent to a  deletion of about %d amino acids [%d %%]\n",
		(int)((double)vara/3.0),
		(int)(((double)(vara)/(double)(seq[INSERT].var_max-seq[INSERT].var_min))*100.0));
	/******************************************************************************/
	/* look at a second enzyme that could be used without T4*/
	fprintf(out," Second Enzyme(s) that do not need(s) polymerase modification:\n  ");
	vara=FALSE; /* no enzyme found */
	for(i=0;i<nbr_sites;i++)
		{
		/*discard if i is not compatible with NumSiteI_5 */
		if(Fct_Test(NumSiteI_5,i,NO_TREATMENT,TRUE,&Strategy->Couple[SIDE_3])==FALSE)
				continue;
		/* i is in the cloning box on the left of NumSiteI_5 , nowhere else */
		if(Fct_N_sites(Sites[i].NumEnz,INSERT,
						Sites[NumSiteI_5].Loc+Enzymes[Sites[NumSiteI_5].NumEnz].taille_site,
						seq[INSERT].var_max,
						TRUE)>0 &&
		  Fct_N_sites(Sites[i].NumEnz,INSERT,
						Sites[NumSiteI_5].Loc-1,
						seq[INSERT].var_max+1,
						FALSE)==0)
				{
				vara=TRUE;/* found an enzyme */
				fprintf(out,"%s (%d) ",Enzymes[Sites[i].NumEnz].nom,Sites[i].Loc);
				break;
				}
		}	
	fprintf(out,"%s",((vara==FALSE) ? "no enzyme found.\n":".\n"));

	/******************************************************************************/
	/* look at a second enzyme that could be used WITH T4*/
	if(Preference.allow_T4!=FALSE)
		{
		fprintf(out," Second Enzyme(s) that need(s) polymerase modification:\n  ");
		vara=FALSE; /* no enzyme found */
		for(i=0;i<nbr_sites;i++)
			{
			/* i is in the cloning box on the left of NumSiteI_5 , nowhere else */
			if(Fct_N_sites(Sites[i].NumEnz,INSERT,
							Sites[NumSiteI_5].Loc+Enzymes[Sites[NumSiteI_5].NumEnz].taille_site,
							seq[INSERT].var_max,
							TRUE)>0 &&
			  Fct_N_sites(Sites[i].NumEnz,INSERT,
							Sites[NumSiteI_5].Loc-1,
							seq[INSERT].var_max+1,
							FALSE)==0)
					{
					vara=TRUE; /* found an enzyme */
					fprintf(out,"%s (%d) ",Enzymes[Sites[i].NumEnz].nom,Sites[i].Loc);
					break;
					}
			}	
		fprintf(out,"%s",((vara==FALSE) ? "no enzyme found.\n":".\n"));
		}
	Save_alignment(NumSiteI_5,NumSiteI_5,NumSiteI_5,NumSiteI_5,DELTA_TREATMENT,TRUE,out,mode);
		
	printTitle(s,2,FALSE);
	return(TRUE);
	}
/********************************************************************************/
  /*********************************/
 /* Display an in frame  deletion */
/*********************************/

int Solution_Frame( short mode,FILE *out,
					int NumSiteI_5, int partialI_5,
					int NumSiteI_3,int partialI_3, 
					int _TREATMENT)
	{
	int i,vara;
	double varr,vard,varb,varc;
	char s[40];
	STRUCT_STRATEGY *Strategy=NULL;
	
	if((Strategy=(STRUCT_STRATEGY*)malloc(sizeof(STRUCT_STRATEGY)))==NULL) return(FALSE);
	
	/* init strategy in order to use the function Solution_2 */
	Strategy->Couple[SIDE_5].Test_Trans=Strategy->Couple[SIDE_3].Test_Trans=FALSE;
	Strategy->Couple[SIDE_5].NumSite[VECTOR]		= NumSiteI_3;
	Strategy->Couple[SIDE_5].Treatment[VECTOR]		= _TREATMENT;
	Strategy->Couple[SIDE_5].Partial[VECTOR]		= 0;
	Strategy->Couple[SIDE_3].NumSite[VECTOR]		= NumSiteI_5;
	Strategy->Couple[SIDE_3].Treatment[VECTOR]		= _TREATMENT;
	Strategy->Couple[SIDE_3].Partial[VECTOR]		= 0;
	
	Strategy->Couple[SIDE_5].NumSite[INSERT]		= NumSiteI_5;
	Strategy->Couple[SIDE_5].Treatment[INSERT]		= _TREATMENT;
	Strategy->Couple[SIDE_5].Partial[INSERT]		= partialI_5;
	Strategy->Couple[SIDE_3].NumSite[INSERT]		= NumSiteI_3;
	Strategy->Couple[SIDE_3].Treatment[INSERT]		= _TREATMENT;
	Strategy->Couple[SIDE_3].Partial[INSERT]		= partialI_3;
	if(_TREATMENT==NO_TREATMENT)
		{
		Strategy->Couple[SIDE_5].Pos5_3[VECTOR]= Sites[NumSiteI_5].Loc + Enzymes[Sites[NumSiteI_5].NumEnz].pos5_3;
		Strategy->Couple[SIDE_5].Pos3_5[VECTOR]= Sites[NumSiteI_5].Loc + Enzymes[Sites[NumSiteI_5].NumEnz].pos3_5;
		Strategy->Couple[SIDE_5].Pos5_3[INSERT]= Sites[NumSiteI_5].Loc + Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3;
		Strategy->Couple[SIDE_5].Pos3_5[INSERT]= Sites[NumSiteI_5].Loc + Enzymes[Sites[NumSiteI_3].NumEnz].pos3_5;

		Strategy->Couple[SIDE_3].Pos5_3[INSERT]= Sites[NumSiteI_3].Loc + Enzymes[Sites[NumSiteI_5].NumEnz].pos5_3;
		Strategy->Couple[SIDE_3].Pos3_5[INSERT]= Sites[NumSiteI_3].Loc + Enzymes[Sites[NumSiteI_5].NumEnz].pos3_5;
		Strategy->Couple[SIDE_3].Pos5_3[VECTOR]= Sites[NumSiteI_3].Loc + Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3;
		Strategy->Couple[SIDE_3].Pos3_5[VECTOR]= Sites[NumSiteI_3].Loc + Enzymes[Sites[NumSiteI_3].NumEnz].pos3_5;

		}
	else if(_TREATMENT==T4_TREATMENT)
		{
		Strategy->Couple[SIDE_5].Pos5_3[VECTOR]= Sites[NumSiteI_5].Loc + Enzymes[Sites[NumSiteI_5].NumEnz].pos3_5;
		Strategy->Couple[SIDE_5].Pos3_5[VECTOR]= Sites[NumSiteI_5].Loc + Enzymes[Sites[NumSiteI_5].NumEnz].pos3_5;
		Strategy->Couple[SIDE_5].Pos5_3[INSERT]= Sites[NumSiteI_3].Loc + Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3;
		Strategy->Couple[SIDE_5].Pos3_5[INSERT]= Sites[NumSiteI_3].Loc + Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3;

		Strategy->Couple[SIDE_3].Pos5_3[INSERT]= Sites[NumSiteI_5].Loc + Enzymes[Sites[NumSiteI_5].NumEnz].pos3_5;
		Strategy->Couple[SIDE_3].Pos3_5[INSERT]= Sites[NumSiteI_5].Loc + Enzymes[Sites[NumSiteI_5].NumEnz].pos3_5;
		Strategy->Couple[SIDE_3].Pos5_3[VECTOR]= Sites[NumSiteI_3].Loc + Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3;
		Strategy->Couple[SIDE_3].Pos3_5[VECTOR]= Sites[NumSiteI_3].Loc + Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3;
		}
	/* print header */
	CountStgy(&i);
	fprintf(out,"<A NAME=\"SOL_%d\">",i);
	sprintf(s,"Solution %d/%d",i,TotalSolutions);
	printTitle(s,2,TRUE);
	/* print semi solution */
	SemiSolution(out,INSERT,Strategy,mode);
	/* display figure */
	varr = (double)(seq[INSERT].var_max-seq[INSERT].var_min)/(double)GRAPHIC;
	if(varr<=0.0) ERROR_USER;
	vard = (double)(Sites[NumSiteI_5].Loc-seq[INSERT].var_min)/varr;
	varb = ((double)(Sites[NumSiteI_3].Loc-seq[INSERT].var_min)/varr)-(double)vard;
	varc = (double)GRAPHIC-vard-varb;
	fprintf(out,"  Cloning box boundaries :[%d-%d] [%d-%d].\n",
					seq[INSERT].var_min,
					seq[INSERT].var_min_int,
					seq[INSERT].var_max_int,
					seq[INSERT].var_max);
	fprintf(out,"\n\tOriginal: 5' ");
	for(i=1;i<=((int)vard+(int)varb+(int)varc);i++)
				fprintf(out,"=");
	fprintf(out," 3'\n\tDeletion: 5' ");
	for(i=1;i<=(int)vard;i++) fprintf(out,"=");
	for(i=1;i<=(int)varb;i++) fprintf(out,".");
	for(i=1;i<=(int)varc;i++) fprintf(out,"=");
	fprintf(out," 3'\n");
	if(_TREATMENT==NO_TREATMENT)
		vara=Sites[NumSiteI_3].Loc+Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3-Sites[NumSiteI_5].Loc-Enzymes[Sites[NumSiteI_5].NumEnz].pos3_5;
	else
		vara=Sites[NumSiteI_3].Loc+Enzymes[Sites[NumSiteI_3].NumEnz].pos3_5-Sites[NumSiteI_5].Loc-Enzymes[Sites[NumSiteI_5].NumEnz].pos5_3;	
	fprintf(out,"\n\tEquivalent to a  deletion of %d amino acids [%d %%]\n",
		(int)((double)vara/3.0),
		(int)(((double)(vara)/(double)(seq[INSERT].var_max-seq[INSERT].var_min))*100.0));
	/* detect if ligation creates a stop codon */
	Detect_Stop(out,mode,NumSiteI_5,NumSiteI_3,&Strategy->Couple[SIDE_5]);
	/**********************************************************/	
	/* search for a post ligation digestion localised between I_5 and I_3*/
	vara=FALSE;
	fprintf(out,"Discard non recombinant molecules by a disgestion post-ligation with:<BR><UL>");
	  for(i=0;i<nbr_enzyme;i++)
		{
		if( Fct_N_sites(i,INSERT,
						Sites[NumSiteI_5].Loc+1,
						Sites[NumSiteI_3].Loc-1,
						TRUE) != 0 &&
			Fct_N_sites(i,INSERT,
						Sites[NumSiteI_5].Loc,
						Sites[NumSiteI_3].Loc,
						FALSE) == 0)
			{vara=TRUE;fprintf(out,"<LI>%s",Enzymes[i].nom);}
		}
	fprintf(out,"%s",(vara==TRUE ? ".\n" : "<LI>no enzyme was found.\n"));
	/**********************************************************/
	/* where are the stop codons ? */
	Find_stop_codon(mode,out,NumSiteI_3,seq[INSERT].npb,SIDE_3);
	fprintf(out,"</UL><P>\n");
	/**********************************************************/
	Save_alignment(NumSiteI_5,partialI_5,NumSiteI_3,partialI_3,_TREATMENT,FALSE,out,mode);
	/**********************************************************/
	if(Preference. allow_T4==FALSE)
		fprintf(out,"Beware: you have not allowed the use of modifying polymerase.\n");
	printTitle(s,2,FALSE);
	free(Strategy);
	return(TRUE);
	}

/*******************************************************************************************/

void Show_All_Shift(FILE *out,short mode)
	{
	char c1;
	int i=0,pos,r,NumSite; 
	STRUCT_STRATEGY *save=NULL,*Stgy=NULL;
	bool_t  found_stop=FALSE ;
	
	save=GetCurrStgy();
	
	fprintf(out,"<UL>");
	
	
		/* Scan all the fragment of INSERT */
	for(r=seq[INSERT].var_min;r<=seq[INSERT].var_max;r=r+(LARGEUR_ECRAN+LARGEUR_ECRAN))
		{
		/*fprintf(out,"r=%d<->%d\n",r,r+(LARGEUR_ECRAN+LARGEUR_ECRAN));*/
		fprintf(out,"<LI>");
		SetFirstStgy();
		while(GetCurrStgy()!=NULL)
			{
			Stgy=CurrStgy;
			NumSite=Stgy->Couple[SIDE_5].NumSite[VECTOR];			
			CountStgy(&pos);
			fprintf(out,"<LI><A HREF=\"#SOL_%d\">",pos);
			fprintf(out,"%03d:",pos);
			fprintf(out,"</A>");
			
			found_stop=FALSE;
			for(i=seq[INSERT].var_min;i<=seq[INSERT].var_max;i++)
				{
				if(found_stop==TRUE) break;
				if(Fct_Frame(INSERT,i)==IS_IN_FRAME)
					{
					c1=Fct_Traduction(	seq[INSERT].sequence[Fct_Pos(INSERT,Fct_BaseShift(NumSite,i  ))],
										seq[INSERT].sequence[Fct_Pos(INSERT,Fct_BaseShift(NumSite,i+1))],
										seq[INSERT].sequence[Fct_Pos(INSERT,Fct_BaseShift(NumSite,i+2))]);
					if(i>= Sites[NumSite].Loc && i<=(Sites[NumSite].Loc+Enzymes[Sites[NumSite].NumEnz].taille_site))
						c1=LOWER(c1);
					if(Is_Codon_Stop(c1)==TRUE && i>=Sites[NumSite].Loc) 
						{found_stop=TRUE;}
					/* if there is a stop after the site, break the loop */
					if(i<=r) continue;
					if(i>=r+(LARGEUR_ECRAN+LARGEUR_ECRAN)) break;
					fprintf(out,"%c",c1);
					if(found_stop==TRUE) break;
					}
				}
			if(i>=seq[INSERT].var_max &&  found_stop==FALSE)
				fprintf(out,"[stop]...");
			NextStgy();
			fprintf(out,"\n");
			}
		}
	fprintf(out,"</UL>");
	SetStgy(save);
	fprintf(out,"\n");
	}



/*******************************************************************************************/
bool_t SaveShiftSolutions(short mode)
	{
	STRUCT_STRATEGY *Stgy=NULL;
	FILE *out=NULL;
	int i;
	
	out=(FILE*)stdout;
	
	printf("<HR><CENTER><H1>Frameshift Solutions</H1></CENTER><BR>");
	/*write header */
		printMainHeader(out,mode,TRUE);
		fprintf(out,"</PRE><UL>\n");
		fprintf(out,"<LI><A HREF=\"#SHIFT_SOLUTIONS\">Solutions</A>\n");
		fprintf(out,"<LI><A HREF=\"#ALIGNMENT\">Alignments</A>\n");
		fprintf(out,"<LI><A HREF=\"#ENZYMES_LIST\">Enzymes</A>\n");
		fprintf(out,"<LI><A HREF=\"#MISC\">Miscellanous</A>\n");
		fprintf(out,"</UL><PRE>\n");
			
	
	/* write title 1: solutions */
		fprintf(out,"<A NAME=\"SHIFT_SOLUTIONS\">\n");
		printTitle("FrameShifts Solutions",1,TRUE);
		SetFirstStgy();
		fprintf(out,"<BR>");
		while((Stgy=GetCurrStgy())!=NULL)
			{
			SolutionShift(mode,out,Stgy);
			Digest(mode,out,INSERT,Stgy->Couple[SIDE_5].NumSite[INSERT],Stgy->Couple[SIDE_5].NumSite[INSERT],Stgy->Couple[SIDE_5].Partial[INSERT]+Stgy->Couple[SIDE_5].Partial[INSERT]);
			NextStgy();
			}
		printTitle("",1,FALSE);
	/* write title 1: alignment */
	fprintf(out,"<A NAME=\"ALIGNMENT\">\n");
		printTitle("Alignments",1,TRUE);
		Show_All_Shift(out,mode);
		printTitle("",1,FALSE);
	/* write title 1:enzymes */
		fprintf(out,"<A NAME=\"ENZYMES_LIST\">\n");
		printTitle("Enzymes Used",1,TRUE);
		fprintf(out,"<UL>");
		for(i=0;i<nbr_enzyme;i++) Enzymes[i].select[VECTOR]=FALSE;
		SetFirstStgy();
		while((Stgy=GetCurrStgy())!=NULL)
			{
			Enzymes[Sites[Stgy->Couple[SIDE_5].NumSite[VECTOR]].NumEnz].select[VECTOR]=TRUE;
			Enzymes[Sites[Stgy->Couple[SIDE_3].NumSite[VECTOR]].NumEnz].select[VECTOR]=TRUE;
			Enzymes[Sites[Stgy->Couple[SIDE_5].NumSite[INSERT]].NumEnz].select[VECTOR]=TRUE;
			Enzymes[Sites[Stgy->Couple[SIDE_3].NumSite[INSERT]].NumEnz].select[VECTOR]=TRUE;
			NextStgy();
			}
		SetFirstStgy();	
		for(i=0;i<nbr_sites;i++)
			{
			if(Enzymes[Sites[i].NumEnz].select[VECTOR]==TRUE)
				{Cmd_Get_Info(mode,out,i);}
			Enzymes[Sites[i].NumEnz].select[VECTOR]=FALSE;
			}
		fprintf(out,"</UL>");
		printTitle("",1,FALSE);
	/* write Footer */
	fprintf(out,"<A NAME=\"MISC\">\n");
	printMainHeader(out,mode,FALSE);
	return(TRUE);
	}





/*******************************************************************************************/
bool_t SaveDeltaSolutions(short mode)
	{
	STRUCT_STRATEGY *Stgy=NULL;
	FILE *out=NULL;
	int i;
	
	out=(FILE*)stdout;
	printf("<HR><CENTER><H1>In-Frame Deletions Solutions</H1></CENTER><BR>");
	/*write header */
		printMainHeader(out,mode,TRUE);
		
		fprintf(out,"</PRE><UL>\n");
		fprintf(out,"<LI><A HREF=\"#FRAME_SOLUTIONS\">Solutions</A>\n");
		fprintf(out,"<LI><A HREF=\"#ALIGNMENT\">Alignments</A>\n");
		fprintf(out,"<LI><A HREF=\"#ENZYMES_LIST\">Enzymes</A>\n");
		fprintf(out,"<LI><A HREF=\"#MISC\">Miscellanous</A>\n");
		fprintf(out,"</UL><PRE>\n");
			
	/* write title 1: solutions */
		fprintf(out,"<A NAME=\"FRAME_SOLUTIONS\">\n");
		printTitle("In-Frame Deletions Solutions",1,TRUE);
		SetFirstStgy();
		fprintf(out,"\n");
		while((Stgy=GetCurrStgy())!=NULL)
			{
			switch(Stgy->type)
				{
				case(DELETION_FRAME_NO):
				case(DELETION_FRAME_T4):
					{
					Solution_Frame( mode,out,Stgy->Couple[SIDE_5].NumSite[VECTOR],
						Stgy->Couple[SIDE_5].Partial[VECTOR],
						Stgy->Couple[SIDE_5].NumSite[INSERT],
						Stgy->Couple[SIDE_5].Partial[INSERT],
						(Stgy->type==DELETION_FRAME_NO?NO_TREATMENT:T4_TREATMENT));
					}
					break;
				case(DELETION_CARBOXY):
					{
					Solution_C_Terminal(mode,out,Stgy->Couple[SIDE_5].NumSite[VECTOR],Stgy);
					}
					break;
				default:printf("Error line %d\n",__LINE__);exit(0);break;
				}

Digest(mode,out,INSERT,Stgy->Couple[SIDE_5].NumSite[VECTOR],Stgy->Couple[SIDE_3].NumSite[VECTOR],Stgy->Couple[SIDE_5].Partial[VECTOR]+Stgy->Couple[SIDE_3].Partial[VECTOR]);
if(Stgy->Couple[SIDE_5].NumSite[INSERT]!=Stgy->Couple[SIDE_5].NumSite[VECTOR])
Digest(mode,out,INSERT,Stgy->Couple[SIDE_5].NumSite[INSERT],Stgy->Couple[SIDE_3].NumSite[INSERT],Stgy->Couple[SIDE_5].Partial[INSERT]+Stgy->Couple[SIDE_3].Partial[INSERT]);

			NextStgy();
			}
		printTitle("",1,FALSE);
	/* write title 1: alignment */
	fprintf(out,"<A NAME=\"ALIGNMENT\">\n");
		printTitle("Alignments",1,TRUE);
		Show_All_Alignment(out,mode);
		printTitle("",1,FALSE);
	/* write title 1:enzymes */
		fprintf(out,"<A NAME=\"ENZYMES_LIST\">\n");
		printTitle("Enzymes Used",1,TRUE);
		fprintf(out,"<UL>");
		for(i=0;i<nbr_enzyme;i++) Enzymes[i].select[VECTOR]=FALSE;
		SetFirstStgy();
		while((Stgy=GetCurrStgy())!=NULL)
			{
			Enzymes[Sites[Stgy->Couple[SIDE_5].NumSite[VECTOR]].NumEnz].select[VECTOR]=TRUE;
			Enzymes[Sites[Stgy->Couple[SIDE_3].NumSite[VECTOR]].NumEnz].select[VECTOR]=TRUE;
			Enzymes[Sites[Stgy->Couple[SIDE_5].NumSite[INSERT]].NumEnz].select[VECTOR]=TRUE;
			Enzymes[Sites[Stgy->Couple[SIDE_3].NumSite[INSERT]].NumEnz].select[VECTOR]=TRUE;
			NextStgy();
			}
		SetFirstStgy();	
		for(i=0;i<nbr_sites;i++)
			{
			if(Enzymes[Sites[i].NumEnz].select[VECTOR]==TRUE)
				{Cmd_Get_Info(mode,out,i);}
			Enzymes[Sites[i].NumEnz].select[VECTOR]=FALSE;
			}
		fprintf(out,"</UL>");
		printTitle("",1,FALSE);
	/* write Footer */
	fprintf(out,"<A NAME=\"MISC\">\n");
	printMainHeader(out,mode,FALSE);
	return(TRUE);
	}











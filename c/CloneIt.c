
/* File 'ADN.h' */
/*
____________________________________________________________________________________
Title:		CloneIt (trade mark)
____________________________________________________________________________________
Version		2.0
Date		1998
Langage:	ANSI-C
Author: 	Pierre LINDENBAUM
Adress:		Domaine de Vilvert INRA CRJ VIM
			Laboratoire de Biologie MolŽculaire des rotavirus.
			78352 JOUY EN JOSAS FRANCE
E-Mail:		linden@biotec.jouy.inra.fr
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
			To compile it, copy this file ( that should be named "CloneIt.c") on
			your UNIX directory.
			
			
			Please before compiling, look at this
	
	

			
			To compile ype:
			
			gcc -oCloneIt CloneIt.c
			
	
			this creates a program called "CloneIt".
			Type "CloneIt" to launch the program
			
	The ANSI-C program use another classic applications usually found on UNIX workstation:
	(pico and lynx...) if you do not have them, supress the line you shuold compile using this
	syntax

			gcc -DNOT_UNIX -oCloneIt CloneIt.c
			
	I used to compile the ANSI-C program with the GNU C compiler (gcc). This compiler
	define the macro __GNUC__, if you compile this program with another compiler
	such as 'cc' use this syntax
	
			cc -D__GNUC__ -oCloneIt CloneIt.c
	
____________________________________________________________________________________
Trade Mark:	National Number 97/704582
			Institut National de la PropriŽtŽ Industrielle
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
____________________________________________________________________________________
I would be very glad to receive suggestions and criticisms about this program from the
users...Pierre LINDENBAUM August 1998
____________________________________________________________________________________*/


/*____________________________________________________________________________________*/



#ifdef __GNUC__
	#ifndef NOT_UNIX
		#define VAR_UNIX
	#endif
	
	#include<stdio.h>
	#include<stdlib.h>
	#include<string.h>
	#include<math.h>
	#include<ctype.h>
	/*#include<assert.h>*/
	#include <errno.h>
	#include <time.h>
	
	#ifdef __MAC__
		#include <console.h>
		#include <sioux.h>
		#include <Dialogs.h>
		#include <Types.h>
		#include <Events.h>
		#include <ToolUtils.h>
		#include <Files.h>
		#define CLS					printf("%c",12)		/** CLear Screen**/
		#define DIRECTORY		;
		
	#else
		#define CLS					system("clear")		/** CLear Screen**/
		#define DIRECTORY	{\
							printf("\n\t-----------------\n");\
							printf("\tCURRENT DIRECTORY\n");\
							printf("\t-----------------\n");\
							system("ls");\
							printf("\n");}
	#endif
#else
	#define CLS					printf("%c",12)		/** CLear Screen**/
	#define DIRECTORY		;
		/*if you compile this program with a macintosh and you DO NOT want to have the macintosh 
	capabilities (files & progress bar),then supress the three line */
	#ifdef macintosh
		/*#define __MAC__*/
	#endif
#endif



#pragma mark var_define
#define MAX_NPB					10000				/** Number max of bp within the sequence    **/
#define NOM_MAX_ENZYME			15					/** Length max of a string such as "EcoR1"  **/
#define PALINDROME_MAX_ENZYME	16					/** Length max of a string such as "GAATTC" **/
#define SITE_MAX_ENZYME			20					/** Length max of a string such as "G^AATTC" **/
#define MAX_SITE				100					/** Number max of sites for an enzyme       **/
#define MAX_LENGHT_POLY			30
#define LARGEUR_ECRAN			80					/** Length of the screen                    **/
#define HAUTEUR_ECRAN			50
#define BEEP 					printf("%c",7)		/** Play a system sound   **/
#define MAX_NOM_FICHIER			FILENAME_MAX					/** Length max of a string such as "MyFile" **/
#define VAR_VERSION				"CloneItª V2.0" 	/** Name of the application Please don't change this**/
#define BACK					printf("%c",8);		/** Go Back */
#define POLYLINKER_FILE			"Polylinkers.Set"	/** Polylinkers File **/
#define INKEY               	while (getchar()!='\n') /*** Wait for return ***/
#define VAR_CITATION			"LINDENBAUM Pierre (1998) CloneIt (tm):Finding Cloning strategies. BioInformatics Vol 14, 5 pp 465-466." 
#define VAR_URL					"http://locus.jouy.inra.fr/soft/cloneit/cloneit.html"
#define SETENV					"CLONEHOME"
#define MAX(a,b)				((a)>(b)?(a):(b))
#define MIN(a,b)				((a)<(b)?(a):(b))
#define LOWER(c)				(((c)<=90 && (c)>=65)? (c)+32:(c))
#define UPPER(c)				(((c)>=97 && (c)<=122)?(c)-32:(c))
#define FORMAT_SCREEN			0
#define FORMAT_TEXT				1
#define FORMAT_RTF				2
#define FORMAT_HTML				3
#define EOS						'\0'
#ifndef TRUE
	#define TRUE					1
	#define FALSE					0
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
#define REBASE_FILE				"allenz"
#define DEFAULT_REBASE			"RebaseSmall"
#ifdef VAR_UNIX
	#define PREF_FILE				".CloneItPref"
#else
	#define PREF_FILE				"CloneItPref"
#endif
#define DEFAULT_PROJECT			"DefaultProject"
#define _ABORT					2
#define GRAPHIC					50	/* length of graphical representation of insert */
#define SMALL_SITE				4 	/* length recognized by an enzyme that will be discarded (see optimizing memory in Main menu)*/
#define VEC_5	0
#define INS_5	1
#define VEC_3	0
#define INS_3	1
#define SIDE_5	0
#define SIDE_3	1
#define DRAW_LINE	DRAW_HR(stdout,'-')
#define MAX_SAVE	10
#define MAX_LINE	200 /* max lenght line from file input */
#define RETURN_NULL				{printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);\
								BEEP;INKEY;\
								exit(0);}
#define ERROR_USER				{printf("Error from the user (file %s line %d)\n",__FILE__,__LINE__);\
								INKEY;\
								exit(0);}
#define MAX_CORP	30 /* max number of society */
#define SUB_CLONING 		0
#define DELETION_FRAME_NO	1
#define DELETION_FRAME_T4	2
#define DELETION_CARBOXY	3
#define FRAME_SHIFT			4
#define VAR_DEBUG			2
#define BLUE_COLOR			"0000FF"
#define RED_COLOR			"FF0000"
#define YELLOW_COLOR		"00FF00"
#define NBR_BUFFER			4
#define REBASE_TEMPERATURE	'T'
#define REBASE_BUFFER		'U'
#define REBASE_BUFFERS		'B'

#ifdef VAR_UNIX
	#define GCG_TAG "GCGCOREROOT"
	#define VECTOR_TEMP "__vector_temp0"
	#define INSERT_TEMP "__insert_temp0"
#endif
	#define OUT_TEMP	"__out_temp0"



#define DEBUG(a) if(VAR_DEBUG==(a)) printf("Debugging line %d file %s\n",__LINE__,__FILE__);

/*
	define MENU : used in main(),
	it will be easer to modify the menu,
	to insert a new item,
	in the future with this system
*/

#pragma mark menu
#define MENU_PROJECT	1
#define MENU_VECTOR_OPEN	MENU_PROJECT+1
#ifdef VAR_UNIX
	#define MENU_EDIT_VECTOR 	MENU_VECTOR_OPEN+1
	#define MENU_VECTOR_BOX		MENU_EDIT_VECTOR+1
#else
	#define MENU_VECTOR_BOX		MENU_VECTOR_OPEN+1
#endif
#define MENU_VECTOR_ATG		MENU_VECTOR_BOX+1
#define MENU_VECTOR_ANTI	MENU_VECTOR_ATG+1
#define MENU_VECTOR_MAP		MENU_VECTOR_ANTI+1
#define MENU_INSERT_OPEN	MENU_VECTOR_MAP+1
#ifdef VAR_UNIX
	#define MENU_EDIT_INSERT 	MENU_INSERT_OPEN+1
	#define MENU_INSERT_BOX		MENU_EDIT_INSERT+1
#else
	#define MENU_INSERT_BOX		MENU_INSERT_OPEN+1
#endif
#define MENU_INSERT_ATG		MENU_INSERT_BOX+1
#define MENU_INSERT_ANTI	MENU_INSERT_ATG+1
#define MENU_INSERT_MAP		MENU_INSERT_ANTI+1
#define MENU_DELETION		MENU_INSERT_MAP+1
#define MENU_CARBOXY		MENU_DELETION+1
#define MENU_FRAMESHIFT		MENU_CARBOXY+1
#define MENU_MIN			MENU_FRAMESHIFT+1
#define MENU_MAX			MENU_MIN+1
#define MENU_CHOOSE_AA		MENU_MAX+1
#define MENU_MUTAGENESIS	MENU_CHOOSE_AA+1
#define MENU_SYSTEMATIC		MENU_MUTAGENESIS+1
#define MENU_REBASE			MENU_SYSTEMATIC+1

#ifdef VAR_UNIX
	#define MENU_UPDATE			MENU_REBASE+1
	#define MENU_SHOW_REBASE	MENU_UPDATE+1
	#define MENU_INFO			MENU_SHOW_REBASE+1
#else
	#define MENU_INFO			MENU_REBASE+1
#endif

#define MENU_CONVERT		MENU_INFO+1

#ifdef VAR_UNIX
	
	#define MENU_GCG_BESTFIT_N	MENU_CONVERT+1
	#define MENU_GCG_GAP_N		MENU_GCG_BESTFIT_N+1
	#define MENU_GCG_BESTFIT_P	MENU_GCG_GAP_N+1
	#define MENU_GCG_GAP_P		MENU_GCG_BESTFIT_P+1
	#define MENU_DISPLAY		MENU_GCG_GAP_P+1
	
	/*#define MENU_GCG_GAP_N		MENU_CONVERT+1
	#define MENU_GCG_GAP_P		MENU_GCG_GAP_N+1
	#define MENU_DISPLAY		MENU_GCG_GAP_P+1*/
	
	
#else
	#define MENU_DISPLAY		MENU_CONVERT+1
#endif


#define MENU_MEMORY			MENU_DISPLAY+1
#define MENU_T4				MENU_MEMORY+1
#define MENU_PARTIAL		MENU_T4+1
#define MENU_BLUNT_PARTIAL	MENU_PARTIAL+1
#define MENU_CONT			MENU_BLUNT_PARTIAL+1
#define MENU_OVERLAPP		MENU_CONT+1

#define MENU_TEMPERATURE	MENU_OVERLAPP+1
#define MENU_BUFFER			MENU_TEMPERATURE+1

#define MENU_CIP			MENU_BUFFER+1
#define MENU_NH2			MENU_CIP+1
#define MENU_COOH			MENU_NH2+1
#define MENU_INTERSECTION	MENU_COOH+1
#define MENU_CLONEIT		MENU_INTERSECTION+1
#define MENU_QUIT			MENU_CLONEIT+1
#ifdef VAR_UNIX
	#define MENU_NEWSEQ			MENU_QUIT+1
	#define MENU_LYNX			MENU_NEWSEQ+1
#endif

/*__________________________________________________________________*/

/************************************************************
 ************************************************************
                   global graphical description
 ************************************************************
 ************************************************************

				SIDE_5                    SIDE_3
				EnzymeV5                  EnzymeV3
				NumSiteV5                 NumSiteV3  
				V5                        V3
				|                         |
			var_min                       |      var_max
		     VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV 
		   V                                      V
		  V                                        V
		 V                                          V
		 V                                          V
		V                                            V
		V                 VECTOR                     V
		V                                            V
		V                                            V
		V                                            V
		 V                                          V
		 V                                          V
		  V                                        V
		   V                                      V
		     V                                   V
		       V                               V
		          V                          V
		             V                    V
		                 V            V
		                   V V V V V
		
				SIDE_5                    SIDE_3
				EnzymeI5                  EnzymeI3
				NumSiteI5                 NumSiteI3  
				I5                        I3
				|                         |
			var_min                       |     var_max
		      IIIIIIIII---------------IIIIIIIIIII
		   I          |               |           I
		  I           sigma5          sigma3       I
		 I                                          I
		 I                                          I
		I                                            I
		I                                            I
		I                  INSERT                    I
		I                                            I
		I                                            I
		 I                                          I
		 I                                          I
		  I                                        I
		   I                                      I
		     I                                   I
		       I                               I
		          I                          I
		             I                    I
		                 I            I
		                   I I I I I

 ************************************************************
 ************************************************************/
/*_______ structures def ___________________________________________________________*/

#pragma mark Typedef

typedef char boolean;

/* this structure describes an enzyme */
typedef struct StructEnz
	{
	char 	nom[NOM_MAX_ENZYME];				/* name */
	char 	site[PALINDROME_MAX_ENZYME];		/* site GAATTC */
	char 	site_complet[SITE_MAX_ENZYME];		/* raw site G^AATTC */
	int 	taille_site;						/* lentgh of site: 6 */
	boolean	select[2];		/* is this enzyme is selected for 0:VECTOR or 1:INSERT */
	int 	pos5_3;			/* distance on strand 5'->3' to the cuting site */
	int 	pos3_5;			/* distance on strand 3'->5' to the cuting site */
	boolean	palindromic;	/* Is this Enzyme palindromic or not ? */
	int		Nbr_N5;			/* number of aNy base on side 5' of site ex: NNATG (used in function GetSite)*/
	int		Nbr_N3;			/* number of aNy base on side 3' of site ex: GTANN */
	/* new ! in this version */
	char	temperature;/* optimal temperature */
	char	NEBBuffer;
	char	NBuffer[NBR_BUFFER]; /*percentage of activity in the 4 standard NEB's buffer */
	short	corp;/* number of customers */
	} STRUCT_ENZYME;

/* this structure describes the preferences */
typedef struct Prefs
	{
	int			partial;		/* max number of partial digestions */
	int			var_pct_min;	/* percentage of insert that will be considered in the cloning box by default  */
	int			DeltaMax;		/* max percentage of insert  that will be deleted see:finding in frame deletion */
	int 		DeltaMin;		/* min percentage of insert  that will be deleted see:finding in frame deletion */
	boolean 	allow_T4;		/* allow modifying polymerase */
	boolean 	allow_CIP;		/* allow phosphatase */
	boolean		allow_all_sol;	/* show all solutions for a couple of enzymes */
	boolean		allow_part_overhang; /* allow ligation between 2 sites that are compatibles but don't have the same overhang lenght */
	boolean		display_messages; 
	boolean 	side_5;	/*clone INSERT in frame in 5' */
	boolean		side_3; /*clone INSERT in frame in 3' */
	boolean		memory; /* discard short ,too degenerate, or ubiquitous enzymes */
	boolean		partial_only_blunt; /* use partial digestion only if the enzyme cuts blunt */
	boolean		search_C_term; /* search for C terminal deletion */
	boolean		mini_display; /* save a few result (developped for cloning projects) */
	boolean		temperature; /* allow non compatible t¡C*/
	boolean		buffer; /* allow non compatible buffer*/
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
	boolean			Test_Trans;	/* TRUE: VECTOR and INSERT are considered FALSE: only INSERT */
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
	boolean				use_CIP;	/* is phosphatase used */
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
#pragma mark declarations
  /*****************************************/
 /* Functions and procedures declarations */
/*****************************************/
void	Handle_error(long errnumber);
void	ChooseATG(int NumSeq);
void	DRAW_HR(FILE *out,char car);
void	EchangeStgy(STRUCT_STRATEGY *A, STRUCT_STRATEGY *B);
void	SortStgys(void);
int 	Sign(int vara);
int 	Fct_Frame(int NumSeq,int nombre);
int		Fct_Pos(int NumSeq,int position);
FILE 	*Fct_Write_File(char var_file_name[MAX_NOM_FICHIER]);
void	Cmd_lire(char *word);
int 	Get_Number(int _vara);
int 	Fct_type(STRUCT_ENZYME Enzyme);
int 	Fct_est_ADN(char lettre);
char	Fct_Complementaire(char lettre);
int 	Fct_Identique(char lettre1,char lettre2);
char	Fct_Reverse(char *mot);
int 	Is_Codon_Stop(char AminoAcid);
char	Fct_Traduction(char c1,char c2,char c3);
char	Translation_at(short NumSeq,int pos);
void	ShowSeq(short mode,FILE *out,STRUCT_SITE site_5, STRUCT_SITE site_3);
void	Display(char *s,int mode);
void	Menu(int vara);
void	Rebase2Strider(void);
void	Strider2Rebase(void);
int 	ResMap(FILE *out,short mode,int var_seq);
boolean	Intersections(FILE *out, short mode);
boolean Are_Same_Asymetric(int NumSite_3,int NumSite_5);
int		Cmd_Get_Info_Choice(int NumEnz1,int NumEnz2,int NumEnz3,int NumEnz4);
int 	Cmd_Get_Info2(short mode, FILE *out,int NumEnz);
void	Cmd_Get_Info(short mode, FILE *out,int NumEnz);
void	search_info(void);
void	Discard_Small(void);
void	Cmd_Get_Rebase(void);
int 	Add_Site(int k,int _NumSeq,int _Loc);
void	Cmd_Get_Site(void);
int 	Get_ATG(int NumSeq);
int 	Cmd_POLYLINKER2(int NumSeq,char mot1[],char mot2[],char mot3[]);
void	Cmd_POLYLINKER1(int NumSeq);
boolean AntiParallele(int NumSeq);
boolean	Cmd_LOAD_SEQUENCE(int NumSeq);
boolean Are_in_Frame(STRUCT_STRATEGY_2 *Stg_2);
int 	Fct_Test(int NumSiteA,int NumSiteB,int Treatment ,int _Frame, STRUCT_STRATEGY_2 *Stg_2);
void	Find_stop_codon(short mode,FILE *out,int NumSite,int var_limit,int var_side);
int 	Fct_N_sites(int Numenz,int NumSeq,int _min,int _max,int Boolean);
void	Init_Enz(void);
void	Digest(short mode,FILE *out, int NumSeq,int NumEnz1, int NumEnz2, int AlertPartial);
boolean Fct_Compatible2( int SeqA, int varA_5_3, int varA_3_5, int SeqB, int varB_5_3, int varB_3_5);
int 	Fct_Compatible_No_Treatment (int NumSiteA,int NumSiteB, STRUCT_STRATEGY_2 *Stg_2);
int 	Fct_Compatible(int NumSiteA,int NumSiteB,int Treatment,STRUCT_STRATEGY_2 *Stg_2);
void	save_it(STRUCT_STRATEGY *Strategy, int CloningType);
void	Orient_Insert_After_CIP(STRUCT_STRATEGY *Strategy,short mode,FILE *out);
int 	Solution_1(STRUCT_STRATEGY Strategy);
boolean	Detect_Stop(FILE *out,short mode,int NumSiteA, int NumSiteB,STRUCT_STRATEGY_2 *Stg_2);
int 	Test_Use_CIP(int NumSeq,STRUCT_STRATEGY Stg);
void	Treatments(int *Treat_5, int *Treat_3,int Type5,int Type3);
int		Sub_Cloning(void);
int 	Fct_BaseShift(int NumSite, int Position);
int 	FrameShift(void);
void	Save_alignment( int NumSiteI_5, int partialI_5,int NumSiteI_3,int partialI_3,boolean _TREATMENT,boolean it_is_a_C_term_del ,FILE *out, short mode);
int 	Solution_C_Terminal(short mode, FILE *out,int NumSiteI_5,STRUCT_STRATEGY *Strategy);
int 	Solution_Frame(short mode, FILE *out,int NumSiteI_5, int partialI_5,int NumSiteI_3,int partialI_3, int _TREATMENT);
int		DeltaFrame(void);
char 	*search_arg(char *word,int argc, char *argv[]);
int 	val(char *word, int var_defaut);
void 	Cmd_Line(int argc, char *argv[]);
void	Cmd_init_preferences(STRUCT_PREFS *Pref);
void 	Save_Pref(STRUCT_PREFS *Pref);
void	Cloning_Project1(void);
void	Cloning_Project2(char *MyChoice);
int		META_CloneIt4(boolean arg_cip, boolean arg_t4,boolean arg_part);
int		META_FrameShift(boolean arg_part);
int		META_DeltaFrame(boolean arg_t4,boolean arg_part);
/* 21/06/98 */
boolean TestStgy(void);
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
int Make_New_enzyme(char *EnzymeName, char *EnzymeSite,int Enzyme_Available, char temp, char buffer, char NEBbuffer[NBR_BUFFER]);
int Discard_iso(void);
void InitStgy(void);
void DeleteAllStgys(void);
void Show_All_Alignment(FILE *out,short mode);
void Show_All_Shift(FILE *out,short mode);
void Show_cloning_solution(FILE *out,STRUCT_STRATEGY *Strategy,short mode);
void Cloning_solution_loop(short mode);
void SemiSolution(FILE *out,int NumSeq,STRUCT_STRATEGY *Stg,short mode);
void printHeader(FILE *out,STRUCT_STRATEGY *Strategy,short mode);
void printEnzymeName(int NumEnz,short mode, FILE *out);
void FrameLoop(short mode);
void ShiftLoop(short mode);
void SolutionShift(short mode, FILE *out,STRUCT_STRATEGY *Strategy);
void printWebEnzyme(char *name,FILE *out);
void printTitle(FILE *out,short mode, char *title,short decal, boolean begin);
void write_date(FILE *out);
void printMainHeader(FILE *out,short mode,boolean begin);
void printEntreprise(short mode, FILE *out,char c);
boolean SaveCloningSolutions(short mode);
boolean SaveDeltaSolutions(short mode);
boolean SaveShiftSolutions(short mode);
boolean freadRebase(FILE *Buffer, char *mot1, char *mot2, char *mot3, char *mot4, char *mot5, char mot6[NBR_BUFFER]);
boolean ReadStriderFormat(int NumSeq);
boolean Temperature_Compatible(int NumSiteA,int NumSiteB);
boolean	Display_Temp_Buffer(FILE *out,short mode,int NumEnzA,int NumEnzB);
boolean	Test_Temp_and_Buffer(int NumSiteA,int NumSiteB);
short	Buffer_Compatible2(int NumSiteA,int NumSiteB);
short 	Buffer_Compatible1(int NumEnzA,int NumEnzB);
boolean Fct_coded(char triplet1[4], char triplet2[4], boolean *x);
boolean IsAACodedBy(char TheAA,char c1,char c2,char c3);
short	Fct_compare_AA(char aa1, char aa2);
void	draw_seq_aound_aa(FILE *out,const short mode,const int pos);
boolean ChooseTheAA(void);
boolean DirectMutagenesis(FILE *out, short mode,int Pos);
boolean SystematicMutageneis(FILE *out, short mode);
boolean PrintStratregy(short type,STRUCT_STRATEGY *Begin,STRUCT_STRATEGY *End);

/* the next functions use the unix programs PICO and LYNX */
#ifdef VAR_UNIX
	/* GCG routines */
	void	GCGReformat(char *TheSeq);
	boolean exportTheSequence(char *SeqName, int NumSeq, boolean SeqIsDNA);
	boolean GCGBestFit(boolean SeqIsDNA);
	boolean GCGGap(boolean SeqIsDNA);
	/* UNIX system routine */
	boolean EditNewSequence(void);
	boolean EditSequence(int NumSeq);
	boolean ShowRebase(void);
	boolean GotoLynx(char *URL);
	boolean  Check_Mail(void);
#endif

/* structure used in with DNA strider */

#ifdef	__GNUC__
	typedef int StriderInt;
#else
	typedef long StriderInt;
#endif

typedef struct strider_header
	{
	char		versionNb;
	char		type;
	char		topology;
	char		reserved1;
	StriderInt	reserved2;
	StriderInt	reserved3;
	StriderInt	reserved4;
	char		reserved5;
	char		filler1;
	short		filler2;
	StriderInt	filler3;
	StriderInt	reserved6;
	StriderInt	nLength;
	StriderInt	nMinus;
	StriderInt	reserved7;
	StriderInt	reserved8;
	StriderInt	reserved9;
	StriderInt	reserved10;
	StriderInt	reserved11;
	char		reserved12[32];
	short		reserved13;
	short		filler4;
	char		reserved14;
	char		reserved15;
	char		reserved16;
	char		filler5;
	StriderInt	com_length;
	StriderInt	reserved17;
	StriderInt	filler6;
	StriderInt	filler7;
	}STRIDER_HEADER;







#ifndef __MAC__
	#define get_char() getchar()
#endif

#pragma mark fonction_macintosh
#ifdef __MAC__
	char Macfgetc(short sRefFichier);
	Boolean	fGetFileName	(Str255 str255Nom, short *psRefVolume, char *name);
	Boolean OpenSelectedFile(Str255 str255NomFichier, short sReferenceVolume,short NumSeq);
	Boolean	fOuvreFichier	(Str255 str255Nom, short sRefVolume, short *psRefFichier);
	void	_FermeFichier	(short sRefFichier);
	Boolean _OuvreFichierSequence(int NumSeq,short sRefFichier);
	Boolean	fLitFichier	(short sRefFichier, TEHandle hteTexte);
	Boolean	fEcritFichier	(short sRefFichier, char *pcTexte, long lTaille);
	Boolean	fEnregistreFichier	(Str255 str255Nom, short sRefVolume, TEHandle hteTexte);
	Boolean	fEnregistreSousFichier (Str255 str255Nom, short *psRefVolume, TEHandle hteTexte);
	Boolean	fSelectionneNouveauFichier	(Str255 str255Nom, short *psRefVolume);
	Boolean	fCreeFichier	(Str255 str255Nom, short sRefVolume, short * psRefFichier);
	long	lCalculeTailleFichier	(Str255 str255Nom, short sRefVolume);
	void  MacAlert(char *word);
	void ShowProgressBar(short state,int count, int the_max,char *TheMessage);
	boolean MacReadStriderFormat(int NumSeq,short sRefFichier);
	boolean MacDialogPref(void);
	boolean ExecuterDialogPref(void);
	void	MacAbout(void);
	char 	MacBrowser(short type, int index);
	FILE 	*MacWriteTextFile(char *Nom_du_Fichier,Str255 originalName);
	#define get_char() getch()
	boolean DoneDrawingPagesOfATextFile(FILE *someFile);
#endif
/*____________________________________________________________________________________*/
#pragma mark variable_extern
#ifndef __GNUC__
	#ifndef _MAIN_C_
		extern int	npb[2],var_max[2],var_min[2],degenerate[2],sigma_5,sigma_3;
		extern char	sequence[2][MAX_NPB]; 
		extern char	FICHIER_ADN[2][MAX_NOM_FICHIER];
		extern char	FICHIER_ENZYME[MAX_NOM_FICHIER];
		extern int	pos_ATG[2];
		extern boolean IsStriderSeq[2];
		extern STRUCT_PREFS	Preference;
		extern int 	nbr_enzyme;
		extern STRUCT_ENZYME   *Enzymes;
		extern int	nbr_sites;
		extern STRUCT_SITE 	*Sites;
		extern char  *string_seq[2];
		extern STRUCT_ENZXY *EnzXY;
		extern int var_nbr_xy;
		extern int NumFragment;
		extern STRUCT_FRAGMENT *TheFragment;
		extern STRUCT_STRATEGY *FirstStgy;
		extern STRUCT_STRATEGY *LastStgy;
		extern STRUCT_STRATEGY *CurrStgy;
		extern STRUCT_STRATEGY *NewStgy;
		extern boolean Search_Done;
		extern int true_nbr_enzyme;
		extern int TotalSolutions;
		extern char TheAminoAcid;
		extern int pos_of_the_AA;
	#endif
#endif




/* File 'Main.c' */
/*ansi_prefix.mac.h*/
#define _MAIN_C_
#ifndef __GNUC__
	#include <stdio.h>
	#include <errno.h>
	#include <stdlib.h>
	#include <string.h>
	#include <stddef.h>
	#include "ADN.h"
	#ifdef __MAC__
		#include <sioux.h>
		#include <Dialogs.h>
		#include <console.h>
	#endif
#endif


int				npb[2],var_max[2],var_min[2],degenerate[2],sigma_5,sigma_3;
/*
	npb 				number of bases in each sequence
	var_min/var_max		cloning boxes boundaries
	degenerate			if sequence's got degenerate bases (should not !)
	sigma_5/sigma_3		internal VECTOR cloning boxes boundaries
*/
char 			sequence[2][MAX_NPB]; /* bases of two sequences */
char			FICHIER_ADN[2][MAX_NOM_FICHIER]; /*names of two sequences */
char			FICHIER_ENZYME[MAX_NOM_FICHIER]; /* name of rebase file */
int				pos_ATG[2]; /* position of a base that should be in frame with the translated sequence */
boolean			IsStriderSeq[2];
STRUCT_PREFS	Preference; /* user's preferences */
int 			nbr_enzyme; /* number of enzymes in the  list of Enzymes*/
int 			true_nbr_enzyme; /* the true number of enzymes ( non palindromic count twice)*/
STRUCT_ENZYME   *Enzymes; /* a pointer to the list of Enzymes */
int 			nbr_sites; /* number of sites in memory */
STRUCT_SITE 	*Sites; /* a pointer to the list of sites */
char 			*string_seq[2]={"VECTOR","INSERT"};
STRUCT_ENZXY	*EnzXY; /* used in function Resmap */
int				var_nbr_xy=0; 	/* used in function Resmap */
int				NumFragment;/* used in Digestion */
STRUCT_FRAGMENT *TheFragment;/* a pointer used in Digestion */
boolean			Search_Done=FALSE;
int				TotalSolutions=0;
char		TheAminoAcid=EOS;
int			pos_of_the_AA=0;
/* 21/06/98 */
STRUCT_STRATEGY *FirstStgy=NULL;
STRUCT_STRATEGY *LastStgy=NULL;
STRUCT_STRATEGY *CurrStgy=NULL;
STRUCT_STRATEGY *NewStgy=NULL;


int main (int argc, char *argv[])
	{
	/***********************************************/
	int		i,var_quit=FALSE,vara,varb,varc;
	char 	selection[2]={' ','X'};
	int 	Logic_choice=MENU_PROJECT;
	/***********************************************/
		
	#ifdef __SIOUX__
		#ifdef __MAC__
	 		/*argc = ccommand(&argv);*/
	 		InitDialogs(0L);
	 		InitGraf(&qd.thePort);
			InitWindows();
			MaxApplZone();
			FlushEvents(everyEvent, 0);
		#endif
 		/*char initializeTB Whether to initialize the Macintosh toolbox.*/
		/*char standalone Whether to use your own event loop or SIOUXÕs.*/
		SIOUXSettings.setupmenus=TRUE;
		SIOUXSettings.autocloseonquit = TRUE;
		SIOUXSettings.asktosaveonclose = FALSE;
		SIOUXSettings.showstatusline = FALSE;
		SIOUXSettings.toppixel = 50;
		SIOUXSettings.leftpixel = 50;
		/*SIOUXSettings.tabspaces = 0;*/
		SIOUXSettings.fontid=4;/*monaco*/
		SIOUXSettings.fontface=normal;
		SIOUXSettings.fontsize=9;
		SIOUXSettings.columns = LARGEUR_ECRAN;
		SIOUXSettings.rows = HAUTEUR_ECRAN;
		SIOUXSetTitle("\pCloneIt");	
	#endif
	
	errno=FALSE;
	nbr_enzyme=0;
	nbr_sites=0;
	Search_Done=FALSE;
	strcpy(FICHIER_ENZYME,"");
	Cmd_init_preferences(&Preference);
	
	for(i=VECTOR;i<=INSERT;i++)
		{
		degenerate[i]=FALSE;
		npb[i] = 0;
		var_max[i] = 0;
		var_min[i] = 0;
		pos_ATG[i] = FALSE;
		strcpy(FICHIER_ADN[i],"");
		}
	InitStgy();
	
	/***********************************************/
	CLS;
	Menu(0);
	printf("Please, report bugs to lindenb@biotec.jouy.inra.fr.\n\n");
	#ifdef __MAC__
		MacAlert("La copie de ce logiciel est la propriete unique de Pierre LINDENBAUM et ne doit pas etre utilisŽe par quelqu'un d'autre sous peine de poursuitesÉ");
		printf("It exists a UNIX oriented version of this software. For more information look at the source code.\n");

	#else
		#ifdef VAR_UNIX
		if(getenv(SETENV)==NULL)
			printf("Did you know it ? You can declare the variable %s in your shell to specify the path to a commun REBASE file.\n",SETENV);
        if(getenv(GCG_TAG)==NULL)
			{
			printf("Did you know it ? %s can use the GCG (c) package to align your sequences.\n",VAR_VERSION,SETENV);
			if(getenv("HOST")!=NULL)
				{
				if(strcmp(getenv("HOST"),"biotec")==ARE_IDENTIC || strcmp(getenv("HOST"),"diamant")==ARE_IDENTIC)
					printf("Petite note pour les utilisateurs INRA du batiment des biotechnologies a Jouy: En vous conectant sur TOPAZE, vous benefierez des fonctionnalites de Gap et BestFit de GCG(c). CIAO ! Pierrot.\n");
				}
			}
		#endif
      
	#endif

	printf("Last compilation: %s (%s)\n",__DATE__,__TIME__);
	/* Command line *********************************/
	if(argc>1)
		{
		Cmd_Line(argc,&argv[0]);
		}
	else
		{
		#ifdef __GNUC__
		
			printf("A command-line is available, type \"%s -HELP\" for more informations\n",argv[0]); 
		#endif
		}
	
	/**********************************************/
	while(var_quit==FALSE)
		{
		PARAGRAF;
		PARAGRAF;
		/***********************************************/
		printf(" PROJECT:\n");
			printf("\t%d  Run a project.\n",MENU_PROJECT);
		/***********************************************/
		printf(" %s:\n",string_seq[VECTOR]);
		if(npb[VECTOR]==0)
			printf("\t%d  Open a DNA sequence for %s.\n",MENU_VECTOR_OPEN,string_seq[VECTOR]);
		else
			{
			printf("\t%d  Open a DNA sequence for %s [%s].\n",MENU_VECTOR_OPEN,string_seq[VECTOR],FICHIER_ADN[VECTOR]);
			#ifdef VAR_UNIX
				if(IsStriderSeq[VECTOR]!=TRUE)
					printf("\t%d  Edit %s.\n",MENU_EDIT_VECTOR,FICHIER_ADN[VECTOR]);
			#endif
			printf("\t%d  Define %s cloning box boundaries [%d-%d].\n",MENU_VECTOR_BOX,string_seq[VECTOR],var_min[VECTOR],var_max[VECTOR]);
			if(pos_ATG[VECTOR]==FALSE)
				printf("\t%d  Define ATG position.\n",MENU_VECTOR_ATG);
			else
				printf("\t%d  Define ATG position [%d].\n",MENU_VECTOR_ATG,pos_ATG[VECTOR]);
			printf("\t%d  Set %s to anti-parallele.\n",MENU_VECTOR_ANTI,FICHIER_ADN[VECTOR]);
			printf("\t%d  %s restriction map.\n",MENU_VECTOR_MAP,FICHIER_ADN[VECTOR]);
			}
		/***********************************************/	
		printf(" %s:\n",string_seq[INSERT]);
		if(npb[INSERT]==0)
			printf("\t%d  Open a DNA sequence for %s.\n",MENU_INSERT_OPEN,string_seq[INSERT]);
		else
			{
			printf("\t%d  Open a DNA sequence for %s [%s].\n",MENU_INSERT_OPEN,string_seq[INSERT],FICHIER_ADN[INSERT]);
			#ifdef VAR_UNIX
				if(IsStriderSeq[VECTOR]!=TRUE)
					printf("\t%d  Edit %s.\n",MENU_EDIT_INSERT,FICHIER_ADN[INSERT]);
			#endif
			printf("\t%d  Define %s cloning boxes boundaries [%d-%d] [%d-%d].\n",MENU_INSERT_BOX,string_seq[INSERT],
					var_min[INSERT],
					var_min[INSERT]+sigma_5,
					var_max[INSERT]-sigma_3,
					var_max[INSERT]);
			if(pos_ATG[INSERT]==FALSE)
				printf("\t%d  Define ATG position.\n",MENU_INSERT_ATG);
			else
				printf("\t%d  Define ATG position [%d].\n",MENU_INSERT_ATG,pos_ATG[INSERT]);
			printf("\t%d  Set %s to anti-parallele.\n",MENU_INSERT_ANTI,FICHIER_ADN[INSERT]);
			printf("\t%d  %s  restriction map.\n",MENU_INSERT_MAP,FICHIER_ADN[INSERT]);
				if(nbr_enzyme>0)
				{
				printf(" %s DELETIONS AND FRAMESHIFTS:\n",FICHIER_ADN[INSERT]);
				printf("\t%d Find in frame deletions in INSERT.\n",MENU_DELETION);
				printf("\t%d\tLook for Carboxy terminal deletions [%c].\n",MENU_CARBOXY,selection[Preference.search_C_term]);
				printf("\t%d Find frameshifts in INSERT.\n",MENU_FRAMESHIFT);
				printf("\t%d\tMinimum percentage of insert [%d %%] (%d bp)\n",MENU_MIN,Preference.DeltaMin,(Preference.DeltaMin==0?var_max[INSERT]-var_min[INSERT]:(int)((double)(var_max[INSERT]-var_min[INSERT])*((double)Preference.DeltaMin/100.0))));
				printf("\t%d\tMaximum percentage of insert [%d %%] (%d bp)\n",MENU_MAX,Preference.DeltaMax,(Preference.DeltaMax==0?var_max[INSERT]-var_min[INSERT]:(int)((double)(var_max[INSERT]-var_min[INSERT])*((double)Preference.DeltaMax/100.0))));
				printf(" %s MUTAGENESIS:\n",FICHIER_ADN[INSERT]);
				if(pos_of_the_AA!=0)
					{
					printf("\t%d Choose an AA for directed mutagenesis. [ '%c' at %d]\n",MENU_CHOOSE_AA,TheAminoAcid,pos_of_the_AA);
					printf("\t%d Find a silent mutation at %d.\n",MENU_MUTAGENESIS,pos_of_the_AA);
					printf("\t%d Find all silent mutations on insert.\n",MENU_SYSTEMATIC);
					}
				else
					printf("\t%d Choose an AA for directed mutagenesis.\n",MENU_CHOOSE_AA);
				
				}
			}
		/***********************************************/	
		printf(" REBASE FILE:\n");
		if(nbr_enzyme==0)
			{
			printf("\t%d Open a REBASE File.\n",MENU_REBASE);
			#ifdef VAR_UNIX
				printf("\t%d Update the REBASE File.\n",MENU_UPDATE);
			#endif
			}
		else
			{
			printf("\t%d Open a REBASE File. [%s] (%d enzymes).\n",MENU_REBASE,FICHIER_ENZYME,true_nbr_enzyme);
			#ifdef VAR_UNIX
				printf("\t%d Update the REBASE File.\n",MENU_UPDATE);
				printf("\t%d Show %s.\n",MENU_SHOW_REBASE,FICHIER_ENZYME);
			#endif
			printf("\t%d Get informations about an enzyme...\n",MENU_INFO);
			}
		printf("\t%d Enzymes Database Convertion...\n",MENU_CONVERT);
		/***********************************************/
		#ifdef VAR_UNIX
		if(getenv(GCG_TAG)!=NULL && npb[INSERT]>0 && npb[VECTOR]>0)
			{
			printf(" GCG Wisconsin Package (c) alignments:\n");
			printf("\t%d BestFit (DNA).\n",MENU_GCG_BESTFIT_N);
			printf("\t%d Gap (DNA).\n",MENU_GCG_GAP_N);
			printf("\t%d BestFit (translated sequence).\n",MENU_GCG_BESTFIT_P);
			printf("\t%d Gap (translated sequence).\n",MENU_GCG_GAP_P);
			}
		#endif
		
		
		
		printf(" PARAMETERS:\n");
		printf("\t%d Display messages [%c].\n",MENU_DISPLAY,selection[Preference.display_messages]);
		printf("\t%d Speed optimization (Discard short sites) [%c].\n",MENU_MEMORY,selection[Preference.memory]);
		if(npb[INSERT]>0 && nbr_enzyme>0)
			{	
			printf("\t%d Use Klenow or T4 DNA Pol.[%c].\n",MENU_T4,selection[Preference.allow_T4]);
			printf("\t%d Number of partial digestion allowed = %d.\n",MENU_PARTIAL,Preference.partial);
			if(Preference.partial>0)
				printf("\t%d \tAllow partial digestion only if enzyme blunt [%c].\n",MENU_BLUNT_PARTIAL,selection[Preference.partial_only_blunt]);
			printf("\t%d Continue to search a solution within a couple [%c].\n",MENU_CONT,selection[Preference.allow_all_sol]);
			printf("\t%d Allow non-overlapping overhangs [%c].\n",MENU_OVERLAPP,selection[Preference.allow_part_overhang]);			
			printf("\t%d Allow non-compatible temperatures. [%c].\n",MENU_TEMPERATURE,selection[Preference.temperature]);
			printf("\t%d Allow non-compatible buffers. [%c].\n",MENU_BUFFER,selection[Preference.buffer]);

		if(npb[VECTOR]>0)
				{
				printf("\t%d Use C.I.P. [%c].\n",MENU_CIP,selection[Preference.allow_CIP]);
				printf("\t%d Clone in frame at 5' of insert (NH2) [%c].\n",MENU_NH2,selection[Preference.side_5]);
				printf("\t%d Clone in frame at 3' of insert (COOH)[%c].\n",MENU_COOH,selection[Preference.side_3]);
				}
			
			}
		/***********************************************/
		if((npb[INSERT]>0 || npb[VECTOR]>0) && nbr_enzyme>0)
			{
			printf(" %s:\n",VAR_VERSION);
			printf("\t%d Intersections.\n",MENU_INTERSECTION);
			}
		if(npb[INSERT]>0 && npb[VECTOR]>0 && nbr_enzyme>0)
			{
			printf("\t%d Find strategies to sub-clone insert into vector.\n",MENU_CLONEIT);
			}				
			
		/***********************************************/
#ifdef VAR_UNIX
		printf("\n\t%d Quit.\n",MENU_QUIT);
		printf("\n\t%d Edit a new sequence.\n\n\n\n",MENU_NEWSEQ);
		printf("\n\t%d Get online help.\n\n\n\n",MENU_LYNX);
#else
		printf("\n\t%d Quit.\n\n\n\n",MENU_QUIT);
#endif
		
		/***********************************************/
			 if(npb[INSERT]!=0 && npb[VECTOR]==0 && nbr_enzyme!=0) Logic_choice=MENU_DELETION;
		else if(npb[VECTOR]==0) Logic_choice=MENU_VECTOR_OPEN;
		else if(npb[INSERT]==0) Logic_choice=MENU_INSERT_OPEN;
		else if(var_min[VECTOR]==0 || var_max[VECTOR]==0) Logic_choice=MENU_VECTOR_BOX;
		else if(var_min[INSERT]==0 || var_max[INSERT]==0) Logic_choice=MENU_INSERT_BOX;
		else if(nbr_enzyme==0) Logic_choice=MENU_REBASE;
		else if(Preference.side_5!=FALSE || Preference.side_3!=FALSE)
				{
				if(pos_ATG[INSERT]==FALSE) Logic_choice=MENU_INSERT_ATG;
				else if(pos_ATG[VECTOR]==FALSE) Logic_choice=MENU_VECTOR_ATG;
				else Logic_choice=MENU_CLONEIT;
				}
		else Logic_choice=MENU_CLONEIT;
		/***********************************************/
#ifdef VAR_UNIX
	Check_Mail();
#endif 


		printf("\t\tYour choice ?:/* %d */",Logic_choice);
		varb = Get_Number(Logic_choice);
		
		     if(npb[VECTOR]==0 && (varb>MENU_VECTOR_OPEN && varb <=MENU_VECTOR_MAP)) varb=0;
		else if(npb[INSERT]==0 && (varb>MENU_INSERT_OPEN && varb <=MENU_MAX)) varb=0;
		else if(npb[INSERT]!=0 && (nbr_enzyme==0) && (varb>=MENU_DELETION && varb <=MENU_MAX)) varb=0;
		else if((npb[INSERT]==0 && nbr_enzyme==0) &&  (varb>=MENU_T4 && varb <=MENU_CONT)) varb=0;
		else if((npb[INSERT]!=0 || npb[VECTOR]!=0) && nbr_enzyme==0 && varb==MENU_INTERSECTION) varb=0;
		else if(npb[INSERT]==0 && npb[VECTOR]==0 && nbr_enzyme==0 && (varb>=MENU_CIP && varb <=MENU_CLONEIT)) varb=0;

		switch(varb)
			{
			case(MENU_PROJECT):Cloning_Project1();break;
			case(MENU_INSERT_OPEN):
			case(MENU_VECTOR_OPEN):/************************************/
				{
				varb=(varb==MENU_VECTOR_OPEN?VECTOR:INSERT);
				DIRECTORY;
				DRAW_LINE;
				if(Preference.display_messages==TRUE)
					Menu(1);
				
				#ifndef __MAC__
				printf("\nInput %s sequence name :",string_seq[varb]);
					Cmd_lire(FICHIER_ADN[varb]);
				
					/*if(strcmp(FICHIER_ADN[INSERT],"")==ARE_IDENTIC)
						strcpy(FICHIER_ADN[INSERT],"p5_1");
					if(strcmp(FICHIER_ADN[VECTOR],"")==ARE_IDENTIC)
						strcpy(FICHIER_ADN[VECTOR],"pET25B_plus");*/
				if(strcmp(FICHIER_ADN[varb],"")==ARE_IDENTIC)
					{BEEP;break;}
				printf(" Opening %s.\n",string_seq[varb]);
				#endif
				
				
				Cmd_LOAD_SEQUENCE(varb);
				if(Preference.display_messages==TRUE)
					{
					if(npb[varb]>0)
						{
						printf("\nThe cloning box(es) boundaries and the ATG position have been searched, please check the displayed datas.\n");
						}
					#ifndef __MAC__
						Menu(12);
					#endif
					}
				CLS;
				}break;
			#ifdef VAR_UNIX
				case(MENU_EDIT_VECTOR):
				case(MENU_EDIT_INSERT):
					{
					varb=(varb==MENU_EDIT_VECTOR?VECTOR:INSERT);
					if(IsStriderSeq[varb]!=TRUE)
						EditSequence(varb);
					else
						{
						BEEP;
						printf("Can't Edit a DNA Strider Sequence !\n");
						Menu(12);
						}
					}break;
			#endif
			case(MENU_VECTOR_ANTI):
			case(MENU_INSERT_ANTI):
				{
				varb=(varb==MENU_VECTOR_ANTI?VECTOR:INSERT);
				AntiParallele(varb);
				}break;
			case(MENU_INSERT_MAP):
			case(MENU_VECTOR_MAP):/************************************/
				{
				varb=(varb==MENU_VECTOR_MAP?VECTOR:INSERT);
				if(Search_Done==FALSE) Cmd_Get_Site();
				vara=ResMap(stdout,FORMAT_SCREEN,varb);
				}break;
			case(MENU_INSERT_BOX):	
			case(MENU_VECTOR_BOX):/************************************/
				{
				varb=(varb==MENU_VECTOR_BOX?VECTOR:INSERT);
				if(npb[varb]==0)
					{BEEP;break;}
				Cmd_POLYLINKER1(varb);
				DRAW_LINE;
				if(Preference.display_messages==TRUE)
					Menu(2);
				vara=FALSE;
				while(vara==FALSE)
					{
					printf("Input localisation of 5' boundary for %s /* 1<%d<%d */:",
						string_seq[varb],var_min[varb],npb[varb]);
					varc = Get_Number(var_min[varb]);
					if(varc==0 || varc<1 || varc>npb[varb])
						{BEEP;continue;}
					else
						{var_min[varb]=varc;vara=TRUE;}
					}
				vara=FALSE;
				while(vara==FALSE)
					{
					printf("Input localisation of 3' boundary for %s /* %d<%d<%d */:",
						string_seq[varb],
						var_min[varb]+1,
						var_max[varb],
						npb[varb]);
					varc = Get_Number(var_max[varb]);
					if(varc==0 || varc<var_min[varb]+1 || varc>npb[varb])
						{BEEP;continue;}
					else
						{var_max[varb]=varc;vara=TRUE;}
					}
					vara=FALSE;
					
				if(varb==INSERT)
					{
					while(vara==FALSE)
						{
						printf("\nInput localisation of cloning box 5' INTERNAL limit: /* %d<%d<%d */:"
							,var_min[INSERT]
							,var_min[INSERT]+sigma_5
							,var_max[INSERT]-1);
						varc = Get_Number(var_min[INSERT]+sigma_5);
						if(varc==0 || varc<var_min[INSERT] || varc>var_max[INSERT]-1)
							{BEEP;continue;}
						else
							{sigma_5=varc-var_min[INSERT];vara=TRUE;}
						}
					vara=FALSE;
					while(vara==FALSE)
						{
						printf("Input localisation of cloning box 3' INTERNAL limit: /* %d<%d<%d */:"
							,var_min[INSERT]+sigma_5+1
							,var_max[INSERT]-sigma_3
							,var_max[INSERT]);
						varc = Get_Number(var_max[INSERT]-sigma_3);
						if(varc==0 || varc<var_min[INSERT]+sigma_5+1 || varc>var_max[INSERT])
							{BEEP;continue;}
						else
							{sigma_3=var_max[INSERT]-varc;vara=TRUE;}
						}
					Search_Done=FALSE;
					}
				vara=FALSE;
				}break;
			case(MENU_INSERT_ATG):
			case(MENU_VECTOR_ATG):/************************************/
				{
				if (varb==MENU_VECTOR_ATG)
					ChooseATG(VECTOR);
				else
					ChooseATG(INSERT);
				}break;
			case(MENU_DELETION):/************************************/
				{
				DeltaFrame();
				}break;
			case(MENU_CARBOXY):/************************************/
				{
				if(Preference.display_messages==TRUE)
					Menu(16);
				Preference.search_C_term=((Preference.search_C_term==TRUE) ? FALSE:TRUE);
				CLS;
				}break;
			case(MENU_FRAMESHIFT):/************************************/
				{
				FrameShift();
				}break;
			case(MENU_MIN):/************************************/
				{
				vara=FALSE;
				while(vara==FALSE)
					{
					DRAW_LINE;
					printf("The Minimum percentage of deletion is the smaller percentage of total insert (100%) ");
					printf("that you want to allow.\n");
					DRAW_LINE;
					printf("\nMinimum percentage of deletion :\n");
					printf("\t    Your choice /* 0<%d<%d */:",Preference.DeltaMin,Preference.DeltaMax);
					varb = Get_Number(Preference.DeltaMax);
					if(varb>-1 && varb<Preference.DeltaMax)
						{
						vara=TRUE;
						Preference.DeltaMin=varb;
						}
					}
				CLS;
				}break;
			case(MENU_MAX):/************************************/
				{
				vara=FALSE;
				while(vara==FALSE)
					{
					printf("The Maximum percentage of deletion is the highest percentage of total insert (100%) ");
					printf("that you want to allow.\n");
					printf("\t    Your choice /* %d<%d<100 */:",Preference.DeltaMin,Preference.DeltaMax);
					varb = Get_Number(Preference.DeltaMax);
					if(varb>Preference.DeltaMin && varb<=100)
						{
						vara=TRUE;
						Preference.DeltaMax=varb;
						}
					}
				CLS;
				}break;
			case(MENU_CHOOSE_AA):/************************************/
				{ChooseTheAA();}break;
			case(MENU_MUTAGENESIS):/************************************/
				{DirectMutagenesis(stdout, FORMAT_SCREEN,pos_of_the_AA);Menu(12);}break;
			case(MENU_SYSTEMATIC):/************************************/
				{SystematicMutageneis(stdout, FORMAT_SCREEN);Menu(12);}break;
			case(MENU_REBASE):/************************************/
				{
				DIRECTORY;
				DRAW_LINE;
				if(Preference.display_messages==TRUE)
					Menu(3);
				strcpy(FICHIER_ENZYME,Preference.RebasePref);
				printf("Input REBASE File name :/* %s */",FICHIER_ENZYME);
				Cmd_lire(FICHIER_ENZYME);
				
				if(strcmp(FICHIER_ENZYME,"")==ARE_IDENTIC)
					strcpy(FICHIER_ENZYME,Preference.RebasePref);
				Cmd_Get_Rebase();
				if(nbr_enzyme==0)
					{BEEP;printf("No enzyme was found !.\n");strcpy(FICHIER_ENZYME,"");}
				else
					strcpy(Preference.RebasePref,FICHIER_ENZYME);
				Search_Done=FALSE;
				#ifndef __MAC__
					if(Preference.display_messages==TRUE) Menu(12);
				#endif
				CLS;
				}break;
			#ifdef VAR_UNIX
				case(MENU_SHOW_REBASE):
					{
					ShowRebase();
					}break;
				case(MENU_UPDATE):/*************************************/
					{
					CLS;
 					if(strcmp(getenv("HOST"),"topaze")==ARE_IDENTIC)
						{
						BEEP;printf("Cette operation est impossible sur TOPAZE\n");Menu(12);
						}
					else
					{					
					Display("Update the REBASE File",2);
					printf("\nYour lynx WWW browser will open this URL:\n\n\t http://www.neb.com/rebase/link_ftpdir/ \n\n");
					printf("Please, do download the 'allenz.*' file and then, quit lynx.\n");
					Menu(12);
					GotoLynx("http://www.neb.com/rebase/link_ftpdir/");
					CLS;
					}
					}break;
			#endif
			case(MENU_INFO):/*************************************/
				{
				if(nbr_enzyme>0)
					search_info();
				else BEEP;
				CLS;
				}break;
			case(MENU_CONVERT):/************************************/
				{
				if(nbr_enzyme==0)
						Strider2Rebase();
				else
					{
					DRAW_LINE;
					printf("1 Convert a DNA Striderª RELibrary to a Rebase FILE...\n");
					printf("  This allow you to use your own old Enzymes database.\n");
					printf("2 Export the memorized Enzymes from \"%s\" to a DNA Striderª RELibrary ...\n",FICHIER_ENZYME);
					printf("  This allow you to update your DNA Sriderª RELibrary.\n");
					DRAW_LINE;
					printf("\n  Your Choice ? /* 1 */:");
					vara=Get_Number(1);
					switch(vara)
						{
						case(1):Strider2Rebase();break;
						case(2):Rebase2Strider();break;
						default:BEEP;CLS;break;
						}
					}
				}break;
			#ifdef VAR_UNIX
				case(MENU_GCG_BESTFIT_N):GCGBestFit(TRUE);break;
				case(MENU_GCG_GAP_N):GCGGap(TRUE);break;
				case(MENU_GCG_BESTFIT_P):GCGBestFit(FALSE);break;
				case(MENU_GCG_GAP_P):GCGGap(FALSE);break;
			#endif
			case(MENU_DISPLAY):/************************************/
				{
				Preference.display_messages=((Preference.display_messages==TRUE) ? FALSE:TRUE);
				/* i'll use this later... */
				/*
				#ifndef __MAC__
					Preference.display_messages=((Preference.display_messages==TRUE) ? FALSE:TRUE);
				#else
					MacDialogPref();
				#endif*/
				CLS;
				}break;
			case(MENU_MEMORY):/************************************/
				{
				if(Preference.memory==FALSE)
					{
					/* Now, optimize memory */
					Preference.memory=TRUE;
					if(Preference.display_messages==TRUE) Menu(5);
					/* if there was any Enzymes in memory, discard short */
					if(nbr_enzyme!=0) {Discard_Small();}
					}
				else 
					{
					/* now don't optimize memory */
					Preference.memory=FALSE;
					/* if there was enzymes in memory, reload the rebase file */ 
					if(strcmp(FICHIER_ENZYME,"")!=ARE_IDENTIC)
						{
						printf("You will have to reload the REBASE File.\n");
						Cmd_Get_Rebase();
						if(nbr_enzyme==0)
							{BEEP;printf("No enzyme was found !.\n");
							Menu(12);
							}
						}
					}
				CLS;
				}break;
			case(MENU_T4):/************************************/
				{
				Preference.allow_T4=((Preference.allow_T4==TRUE) ? FALSE:TRUE);
				DRAW_LINE;
				if(Preference.display_messages==TRUE && Preference.allow_T4==TRUE)
					Menu(4);
				CLS;
				}break;
			case(MENU_PARTIAL):/************************************/
				{
				vara=FALSE;
				while(vara==FALSE)
					{
					printf("\nNumber of Partial sites allowed ?:\n");
					printf("\t    Your choice /* 0<0<%d */:",MAX_VAR_PARTIAL);
					varb = Get_Number(0);
					if(varb>=0 && varb<MAX_VAR_PARTIAL)
						{
						vara=TRUE;
						Preference.partial=varb;
						}
					}
				CLS;
				}break;
			case(MENU_BLUNT_PARTIAL):/************************************/
				{
				Preference.partial_only_blunt=((Preference.partial_only_blunt==TRUE) ? FALSE:TRUE);
				if(Preference.display_messages==TRUE)
					Menu(15);
				CLS;
				}break;
			case(MENU_CONT):/************************************/
				{
				Preference.allow_all_sol=((Preference.allow_all_sol==TRUE) ? FALSE:TRUE);
				if(Preference.display_messages==TRUE)
					Menu(6);
				CLS;
				}break;
			case(MENU_OVERLAPP):/************************************/
				{
				DRAW_LINE;
				if(Preference.display_messages==TRUE)
					Menu(7);
				Preference.allow_part_overhang=((Preference.allow_part_overhang==TRUE) ? FALSE:TRUE);
				CLS;
				}break;
			case(MENU_INTERSECTION):/************************************/
				{
				 Intersections(stdout, FORMAT_SCREEN);
				}break;
			case(MENU_CIP):/************************************/
				{
				Preference.allow_CIP=((Preference.allow_CIP==TRUE) ? FALSE:TRUE);
				if(Preference.display_messages==TRUE && Preference.allow_CIP==TRUE)
					Menu(8);
				CLS;
				}break;
			case(MENU_BUFFER):/************************************/
				{
				if(Preference.display_messages==TRUE) Menu(21);
				Preference.buffer=((Preference.buffer==TRUE) ? FALSE:TRUE);
				CLS;
				}break;
			case(MENU_TEMPERATURE):/************************************/
				{
				if(Preference.display_messages==TRUE) Menu(21);
				Preference.temperature=((Preference.temperature==TRUE) ? FALSE:TRUE);
				CLS;
				}break;
			case(MENU_NH2):/************************************/
				{
				if(Preference.side_5==TRUE)
					Preference.side_5=FALSE;
				else
					{
					Preference.side_5=TRUE;
					if(Preference.display_messages==TRUE && Preference.side_5==TRUE)
						Menu(9);
					if(pos_ATG[INSERT]==FALSE) Get_ATG(INSERT);
					if(pos_ATG[VECTOR]==FALSE) Get_ATG(VECTOR);
					}
				CLS;
				}break;
			case(MENU_COOH):/************************************/
				{
				if(Preference.side_3==TRUE)
					Preference.side_3=FALSE;
				else
					{
					Preference.side_3=TRUE;
					if(Preference.display_messages==TRUE && Preference.side_3==TRUE)
						Menu(10);
					if(pos_ATG[INSERT]==FALSE) Get_ATG(INSERT);
					if(pos_ATG[VECTOR]==FALSE) Get_ATG(VECTOR);
					}
				CLS;
				}break;			
			case(MENU_CLONEIT):/************************************/
				{
				Sub_Cloning();
				}break;
			#ifdef VAR_UNIX
				case(MENU_NEWSEQ):/************************************/
					{
					EditNewSequence();
					}break;
				case(MENU_LYNX):/************************************/
					{
					GotoLynx(VAR_URL);
					}break;
			#endif
			case(MENU_QUIT):/************************************/
				{
				var_quit=TRUE;
				}break;

			default:CLS;BEEP;break;
			}
		PARAGRAF;
		fflush(stdin);
		}
		Save_Pref(&Preference);
		CLS;
		Menu(19);

		return(TRUE);
	}



/* File 'Preferences.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "ADN.h"
#endif

void Cmd_init_preferences(STRUCT_PREFS *Pref)
	{
	FILE *in;
	int i;
	char ASK[FILENAME_MAX];
	char ThePrefFile[FILENAME_MAX];
	
	#ifdef VAR_UNIX
		if(getenv("HOME")!=NULL)
			{
			sprintf(ThePrefFile,"%s/%s",getenv("HOME"),PREF_FILE);
			}
		else
			sprintf(ThePrefFile,"%s",PREF_FILE);
	#else
			sprintf(ThePrefFile,"%s",PREF_FILE);
	#endif
	
	
	if((in=fopen(ThePrefFile,"rb"))==NULL)
		{
		Menu(22);
		Pref->side_5=FALSE;
		Pref->side_3=FALSE;
		Pref->partial=FALSE;
		Pref->partial_only_blunt=FALSE;
		Pref->allow_CIP=FALSE;
		Pref->var_pct_min= 10;
		Pref->DeltaMax=80;
		Pref->DeltaMin=20;
		Pref->allow_T4=FALSE;
		Pref->memory=TRUE;
		Pref->allow_all_sol=FALSE;
		Pref->allow_part_overhang=FALSE;
		Pref->display_messages=TRUE;
		Pref->search_C_term=TRUE;
		Pref->mini_display=FALSE;
		Pref->buffer=TRUE;
		Pref->temperature=TRUE;
		strcpy(Pref->RebasePref,REBASE_FILE);
		}
	else
		{
		
		if(fread(ASK,strlen(VAR_VERSION)+1,1,in)!=1 || strcmp(ASK,VAR_VERSION)!=ARE_IDENTIC)
			{
			fclose(in);
			fprintf(stderr,"There was a minor problem while initialisation (version number). Please re-start.\n");
			if(remove(ThePrefFile)!=0) fprintf(stderr,"Can't delete %s !\n",ThePrefFile);
			INKEY;
			exit(0);
			}
		if(fread(Pref,sizeof(STRUCT_PREFS),1,in)!=1)
			{
			fclose(in);
			fprintf(stderr,"There was a minor problem while initialisation. Please re-start.\n");
			if(remove(ThePrefFile)!=0) fprintf(stderr,"Can't delete %s !\n",ThePrefFile);
			INKEY;
			exit(0);
			}
		fclose(in);
		}
	
	i=Pref->display_messages;
	Pref->display_messages=FALSE;

	
	#ifdef VAR_UNIX
	/* try to open the rebase file of the shell */
	  if(getenv(SETENV)!=NULL)
		{
		strcpy(FICHIER_ENZYME,getenv(SETENV));
		if((in=fopen(FICHIER_ENZYME,"r"))!=NULL)
			{
			fclose(in);
			Cmd_Get_Rebase();
			}
		}
	#endif
	
	if(nbr_enzyme<=0)
		{
		/* try to open the latest rebase file used */
		strcpy(FICHIER_ENZYME,Pref->RebasePref);
		if(strcmp(FICHIER_ENZYME,"")!=ARE_IDENTIC)
			{
			if((in=fopen(FICHIER_ENZYME,"r"))!=NULL)
				{
				fclose(in);
				Cmd_Get_Rebase();
				}
			else
				{
				/* if it can't, try to open the default one */
				strcpy(FICHIER_ENZYME,DEFAULT_REBASE);
				if((in=fopen(FICHIER_ENZYME,"r"))!=NULL)
					{
					fclose(in);
					Cmd_Get_Rebase();
					}
				else
					{
					/* if it can't, try to create this default file */
					printf("Creating the rebase default file '%s'",DEFAULT_REBASE);
					if((in=fopen(FICHIER_ENZYME,"w"))!=NULL)
						{
						fprintf(in,"%s 1998 Pierre LINDENBAUM", VAR_VERSION);
fprintf(in," This is the default Rebase file it contains only a few enzymes,");
fprintf(in," you should download the \"allenz.*\" file from http://www.neb.com/rebase/link_ftpdir/\n\n");
fprintf(in,"<1>EcoRI\n<5>g^aattc\n<7>?\n<1>BamH1\n<5>g^gatcc\n<7>?\n<1>");
fprintf(in,"HindIII\n<5>a^agctt\n<7>?\n<1>NotI\n<5>gc^ggccgc\n<7>?\n<1>XhoI\n");
fprintf(in,"<5>c^tcgag\n<7>?\n<1>PstI\n<5>CTGCA^G\n<7>?\n<1>NotI\n<5>GC^GGCCGC\n");
fprintf(in,"<7>?\n<1>SmaI\n<5>CCC^GGG\n<7>?\n<1>XmaI\n<5>C^CCGGG\n<7>");
fprintf(in,"?\n<1>SalI\n<5>G^TCGAC\n<7>?\n<1>AvrII\n<5>C^CTAGG\n<7>?\n<1>");
fprintf(in,"BglII\n<5>A^GATCT\n<7>?\n<1>NheI\n<5>G^CTAGC\n<7>?");
						fclose(in);
						Cmd_Get_Rebase();
						}
					}
				}
		}
	
	}
	Pref->display_messages=i;
	}

void Save_Pref(STRUCT_PREFS *Pref)
	{
	FILE *out;
	char ThePrefFile[FILENAME_MAX];
	
	#ifdef VAR_UNIX
		if(getenv("HOME")!=NULL)
			{
			sprintf(ThePrefFile,"%s/%s",getenv("HOME"),PREF_FILE);
			}
		else
			sprintf(ThePrefFile,"%s",PREF_FILE);
	#else
			sprintf(ThePrefFile,"%s",PREF_FILE);
	#endif

	if((out=fopen(ThePrefFile,"wb"))!=NULL)
		{
		fwrite(VAR_VERSION,strlen(VAR_VERSION)+1,1,out);
		if(fwrite(Pref,sizeof(STRUCT_PREFS),1,out)!=1)
			printf("Can't save Preferences !\n");
		fclose(out);
		}
	else
		{Handle_error(0);printf("Can't save Preferences !!\n");}
	}



/* File 'Error.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <errno.h>
	#include "ADN.h"
#endif


  /*********************/
 /*Handle files errors*/
/*********************/
void Handle_error(long errnumber)
	{
	#ifndef __MAC__
	if(errno!=FALSE)
		{
		BEEP;
		#ifdef __ANSI__
		switch(errno)
			{
			case(EACCES):fprintf(stderr,"You don't have the autorization to access.\n");break;
			case(EBUSY):fprintf(stderr,"File is busy, can't delete it.\n");break;
			case(EEXIST):fprintf(stderr,"Duplicate Filename.\n");break;
			case(EIO):fprintf(stderr,"Input/Ouput Error.\n");break;
			case(ENFILE):fprintf(stderr,"Can't open any more file.\n");break;
			case(ENOENT):fprintf(stderr,"File not found.\n");break;
			case(ENOSPC):fprintf(stderr,"Disk is full.\n");break;
			case(ENOTDIR):fprintf(stderr,"Directory not found.\n");break;
			case(EROFS):fprintf(stderr,"Disc is write protected.\n");break;
			default:perror("UnDefined Error :");break;
			}
		errno=FALSE;
		#endif
		}
	#else
	if(errnumber !=0)
		{
		BEEP;
		switch(errnumber)
			{
			case(-33):Display("Directory full",3);break;
			case(-34):Display("disk full",3);break;
			case(-35):Display("no such volume",3);break;
			case(-36):Display("I/O error (bummers)",3);break;
			case(-37):Display("there may be no bad names in the final system!",3);break;
			case(-38):Display("File not open",3);break;
			case(-39):Display("End of file",3);break;
			case(-40):Display("tried to position to before start of file (r/w)",3);break;
			case(-41):Display("memory full (open) or file won't fit (load)",3);break;
			case(-42):Display("too many files open",3);break;
			case(-43):Display("File not found",3);break;
			case(-44):Display("diskette is write protected",3);break;
			case(-45):Display("file is locked",3);break;
			case(-46):Display("volume is locked",3);break;
			case(-47):Display("File is busy (delete)",3);break;
			case(-48):Display("duplicate filename (rename)",3);break;
			case(-49):Display("file already open with with write permission",3);break;
			case(-50):Display("error in user parameter list",3);break;
			case(-51):Display("refnum error",3);break;
			case(-52):Display("get file position error",3);break;
			case(-53):Display("volume not on line error (was Ejected)",3);break;
			case(-54):Display("permissions error (on file open)",3);break;
			case(-55):Display("drive volume already on-line at MountVol",3);break;
			case(-56):Display("no such drive (tried to mount a bad drive num)",3);break;
			case(-57):Display("not a mac diskette (sig bytes are wrong)",3);break;
			case(-58):Display("volume in question belongs to an external fs",3);break;
			case(-59):Display("file system internal error: during rename the old entry was deleted but could not be restored&#133;",3);break;
			case(-60):Display("bad master directory block",3);break;
			case(-61):Display("write permissions error",3);break;
			default:Display("unknow error",3);break;
			}
		
		}
	#endif
	}



/* File 'Extern.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include <stdlib.h>
	#include "ADN.h"
#endif


boolean PrintStratregy(short type,STRUCT_STRATEGY *Begin,STRUCT_STRATEGY *End)
	{
	#ifndef __MAC__
		#ifndef VAR_UNIX
		printf("Not available on this platform..\n");
			return(FALSE);
		#else
	STRUCT_STRATEGY *Stgy=NULL,*Save_Stgy=NULL;
		FILE *out=NULL;
		#ifdef VAR_UNIX
			char line[MAX_NOM_FICHIER];
			int i;
		#endif
		#ifdef __MAC__
			GrafPtr	savedPort;
			TPrStatus	prStatus;
			TPPrPort	printPort;
			THPrint	hPrint;
			boolean	ok;
		#endif
		
		fflush(stdin);
		printf("Printing...\n");
		Save_Stgy=GetCurrStgy();

		if((out=fopen(OUT_TEMP,"w"))==NULL)
			{
			BEEP;
			fprintf(stderr,"Error: Can't print !\n");
			return(FALSE);
			}
		printMainHeader(out,FORMAT_TEXT,TRUE);
		SetStgy(Begin);
		while((Stgy=GetCurrStgy())!=NULL)
			{
			switch(type)
				{
				case(1):Show_cloning_solution(out,Stgy,FORMAT_TEXT);break;
				case(2):
					switch(Stgy->type)
						{
						case(DELETION_FRAME_NO):
						case(DELETION_FRAME_T4):
							{
							Solution_Frame( FORMAT_TEXT,out,
								Stgy->Couple[SIDE_5].NumSite[VECTOR],
								Stgy->Couple[SIDE_5].Partial[VECTOR],
								Stgy->Couple[SIDE_5].NumSite[INSERT],
								Stgy->Couple[SIDE_5].Partial[INSERT],
								(Stgy->type==DELETION_FRAME_NO?NO_TREATMENT:T4_TREATMENT));
							}
							break;
						case(DELETION_CARBOXY):
							{
							Solution_C_Terminal(FORMAT_TEXT,out,Stgy->Couple[SIDE_5].NumSite[VECTOR],Stgy);
							}
							break;
						default:printf("Error line %d\n",__LINE__);exit(0);break;
						}
						break;
				case(3):SolutionShift(FORMAT_TEXT,out,Stgy);break;
				default:break;
				}
			if(Stgy==End) break;
			NextStgy();
			}
		fclose(out);
		#ifdef VAR_UNIX
			sprintf(line,"lp -w %s",OUT_TEMP);
			printf("Send %s\n",line);
			system(line);
		#else
			if((out=fopen(OUT_TEMP,"r"))==NULL) return(FALSE);
			GetPort(&savedPort);
			PrOpen();
			hPrint = (THPrint)NewHandle(sizeof(TPrint));/* *not* sizeof(THPrint) */
			PrintDefault(hPrint);
			ok = PrValidate(hPrint);
			ok = PrJobDialog(hPrint);
			if (ok)
				{
				printPort = PrOpenDoc(hPrint, NULL, NULL);
				SetPort(&printPort->gPort);
				PrOpenPage(printPort, NULL);
				while (!DoneDrawingPagesOfATextFile(out))
					{
					PrClosePage(printPort);	/* Close the currently open page. */
					PrOpenPage(printPort, nil);	/* and open a new one. */
					}
				PrClosePage(printPort);
				PrCloseDoc(printPort);
				/* Print spooled document, if any. */
				if ((**hPrint).prJob.bJDocLoop == bSpoolLoop && PrError() ==noErr)
						PrPicFile(hPrint, nil, nil, nil, &prStatus);
				}
			else
				{
				/* You will want to add error handling here... */
				BEEP;
				}
			PrClose();
			SetPort(savedPort);
			fclose(out);
		#endif
		if(remove(OUT_TEMP)!=0) fprintf(stderr,"Can't delete %s !\n",OUT_TEMP);
		SetStgy(Save_Stgy);
		return(TRUE);
	#endif
	#endif
	}
/******************************************************************************/
void printEnzymeName(int NumEnz,short mode, FILE *out)
	{
	switch(mode)
			{
			case(FORMAT_SCREEN):
			case(FORMAT_TEXT  ):fprintf(out,"%s",Enzymes[NumEnz].nom);break;
			case(FORMAT_HTML  ):printWebEnzyme(Enzymes[NumEnz].nom,out);break;
			}
	}

/******************************************************************************/
void printTitle(FILE *out,short mode, char *title,short decal, boolean begin)
	{
	short i;
	if(begin==TRUE)
		{
		switch(mode)
			{
			case(FORMAT_TEXT):
				{
				for(i=0;i<=decal;i++) fprintf(out,"\t\t");
				fprintf(out,"%s\n",title);
				for(i=0;i<=decal;i++) fprintf(out,"\t\t");
				for(i=0;i<strlen(title);i++) fprintf(out,"%c",(decal==1?'=':'-'));
				fprintf(out,"\n\n");
				}break;
			case(FORMAT_HTML):
				{
				fprintf(out,"<UL><LI><H%d>%s</H%d><P>",decal,title,decal);
				}break;
			default:break;
			}
		}
	else
		{
		switch(mode)
			{
			case(FORMAT_TEXT):
				{
				for(i=decal;i<=3;i++) fprintf(out,"\n");
				}break;
			case(FORMAT_HTML):
				{
				fprintf(out,"</UL><P>",decal,title,decal);
				}break;
			default:break;
			}
		}

	}

/******************************************************************************/
void printMainHeader(FILE *out,short mode,boolean begin)
	{
	register short i;
	if(begin==TRUE)
	{
	switch(mode)
		{
		case(FORMAT_TEXT):
			{
			fprintf(out,"%s (c)1998 Pierre LINDENBAUM\n",VAR_VERSION);
			fprintf(out,"\tBiologie Moleculaire des Rotavirus.\n\tVIM INRA.\n\t78350 JOUY-EN-JOSAS.FRANCE. E.Mail:lindenb@biotec.jouy.inra.fr." );
			fprintf(out,"\n\nDate:");
			write_date(out);
			/* sequences */
			DRAW_HR(out,'.');
			fprintf(out,"\n\nSequence file(s) used:\n");			
			for(i=0;i<=1;i++)
				{
				if(npb[i]>0)
					fprintf(out,"\t%s %s (%d pb) [%d-%d] ATG at:%d\n",(i==0?"Vector":"Insert"),FICHIER_ADN[i],npb[i],var_min[i],var_max[i],pos_ATG[VECTOR]);
				}
			/* rebase */
			fprintf(out,"\nREBASE file used:\n\t%s (%d enzymes)\n",FICHIER_ENZYME,true_nbr_enzyme);
			DRAW_HR(out,'.');
			/* parameters */
			fputs("Parameters:\n",out);
			
			if(Preference.memory!=FALSE) fputs("\tSpeed optimization (Discard short sites).\n",out);
			if(Preference.allow_part_overhang!=FALSE) fputs("\tAllow non-overlapping overhangs.\n",out);
			if(Preference.side_5!=FALSE) fputs("\tClone in frame at 5' of insert (NH2).\n",out);
			if(Preference.side_3!=FALSE) fputs("\tClone in frame at 3' of insert (COOH).\n",out);
			fprintf(out,"\tNumber of partial digestion allowed = %d.\n",Preference.partial);
			if(Preference.partial_only_blunt!=FALSE) fputs("\t\tAllow partial digestion only if enzyme blunt.\n",out);
			fprintf(out,"\t%sllow usage of C.I.A.P.\n",(Preference.allow_CIP==TRUE?"A":"Do NOT a"));
			fprintf(out,"\t%sllow incompatible temperatures.\n",(Preference.temperature==TRUE?"A":"Do NOT a"));
			fprintf(out,"\t%sllow incompatible buffers.\n",(Preference.buffer==TRUE?"A":"Do NOT a"));
			fprintf(out,"\t%sllow usage of modifying polymerase.\n",(Preference.allow_T4==TRUE?"A":"Do NOT a"));
			fprintf(out,"\tDeletion length comprised betwween %d%% and %d%% of the insert.\n",Preference.DeltaMin,Preference.DeltaMax);
			fprintf(out,"\t%sook for Carboxy terminal deletions.\n",(Preference.search_C_term==TRUE?"L":"Do NOT l"));
			fputs("\n",out);
			DRAW_HR(out,'.');
			
			}break;
		case(FORMAT_HTML):
			{

			fputs("<HTML><HEAD><TITLE>\n",out);
			fprintf(out,"CloneIt");
			fputs("</TITLE></HEAD><BODY BACKGROUND=\"http://www.jouy.inra.fr/icons/logo-med2.jpg\">\n",out);
			
			fprintf(out,"<CENTER><H1><A HREF=\"%s\">CloneIt</A></H1><P><H3>by <AUTHOR>Pierre LINDENBAUM</Author></H3>",VAR_URL);
			fprintf(out,"<H6>Biologie Moleculaire des Rotavirus.<BR>VIM <A HREF=\"http://www.inra.jouy.inra.fr\">INRA</A>.<BR>");
			fprintf(out,"78350 JOUY-EN-JOSAS.<BR>FRANCE.<P>E.Mail:<A HREF=\"mailto:lindenb@biotec.jouy.inra.fr\">lindenb@biotec.jouy.inra.fr</A>.</H6></CENTER>" );
			fprintf(out,"Date:");
			write_date(out);
			/* sequences */
			
			fprintf(out,"<HR>Sequence file(s) used:<UL>");			
			for(i=0;i<=1;i++)
				{
				if(npb[i]>0)
					fprintf(out,"<LI>%s %s (%d pb) [%d-%d] ATG at:%d\n",(i==0?"Vector":"Insert"),FICHIER_ADN[i],npb[i],var_min[i],var_max[i],pos_ATG[VECTOR]);
				}
			/* rebase */
			fprintf(out,"</UL>REBASE file used:<UL><LI>%s (%d enzymes)</UL><HR>",FICHIER_ENZYME,true_nbr_enzyme);
			/* parameters */
			fputs("Parameters:<UL>",out);
			
			if(Preference.memory!=FALSE) fputs("<LI>Speed optimization (Discard short sites).\n",out);
			if(Preference.allow_part_overhang!=FALSE) fputs("<LI>Allow non-overlapping overhangs.\n",out);
			if(Preference.side_5!=FALSE) fputs("<LI>Clone in frame at 5' of insert (NH2).\n",out);
			if(Preference.side_3!=FALSE) fputs("<LI>Clone in frame at 3' of insert (COOH).\n",out);
			fprintf(out,"<LI>Number of partial digestion allowed = %d.\n",Preference.partial);
			if(Preference.partial_only_blunt!=FALSE) fputs("<LI>Allow partial digestion only if enzyme blunt.\n",out);
			fprintf(out,"<LI>%sllow incompatible temperatures.\n",(Preference.temperature==TRUE?"A":"Do NOT a"));
			fprintf(out,"<LI>%sllow incompatible buffers.\n",(Preference.buffer==TRUE?"A":"Do NOT a"));
			fprintf(out,"<LI>%sllow usage of C.I.A.P.\n",(Preference.allow_CIP==TRUE?"A":"Do NOT a"));
			fprintf(out,"<LI>%sllow usage of modifying polymerase.\n",(Preference.allow_T4==TRUE?"A":"Do NOT a"));
			fprintf(out,"<LI>Deletion length comprised betwween %d%% and %d%% of the insert.\n",Preference.DeltaMin,Preference.DeltaMax);
			fprintf(out,"<LI>%sook for Carboxy terminal deletions.\n",(Preference.search_C_term==TRUE?"L":"Do NOT l"));
			fputs("</UL><HR><PRE>",out);
			}break;
		default:break;
		}
	}
	else
	{
	switch(mode)
		{
		case(FORMAT_TEXT):
			{
			printTitle(out,mode,"Miscellanous",1,TRUE);
			fprintf(out,"\tBibliography\n");
			fprintf(out,"\t\tMarck, C. (1988). 'DNA Strider': a 'C' program for the fast analysis of DNA and sequences on the Apple Macintosh family of computers. Nucleic Acids Res 16, 1829-36.\n");
			fprintf(out,"\t\tRoberts, R. J. & Macelis, D. (1998). REBASE--restriction enzymes and methylases. Nucl. Ac. Res. 26, 338-350.\n");
			fprintf(out,"\t\t%s\n",VAR_CITATION);
			fprintf(out,"\tAcknowledgments\n");
			fprintf(out,"\t\tI would like to thank Audrey Nepveu-de-Villemarceau, Janine L., Dr S. Hazout and his team, Philippe Bessieres, Christine Young, and Maria Piron for their help.\n");
			printTitle(out,mode,"",1,FALSE);
			fprintf(out,"%s by Pierre LINDENBAUM\n",VAR_VERSION);
			fprintf(out,"\tBiologie Moleculaire des Rotavirus.\n\tVIM INRA.\n\t78350 JOUY-EN-JOSAS.FRANCE. E.Mail:lindenb@biotec.jouy.inra.fr." );
			}break;
		case(FORMAT_HTML):
			{
			printTitle(out,mode,"Miscellanous",1,TRUE);
			fprintf(out,"</PRE><UL>");
			fprintf(out,"<LI>Related Links<UL>");
			fprintf(out,"<LI><A HREF=\"http://www.jouy.inra.fr/index.html\">INRA</A> home page");
			fprintf(out,"<LI><A HREF=\"%s\">%s</A> home page",VAR_URL,VAR_VERSION);
			fprintf(out,"<LI>Dr. Richard Robert's<A HREF=\"http://www.neb.com/rebase\"> REBASE</A> at N.E.B.",VAR_URL,VAR_VERSION);
			fprintf(out,"</UL><LI>Bibliography<UL>");
			fprintf(out,"<LI>Marck, C. (1988). 'DNA Strider': a 'C' program for the fast analysis of DNA and sequences on the Apple Macintosh family of computers. Nucleic Acids Res 16, 1829-36.");
			fprintf(out,"<LI>Roberts, R. J. & Macelis, D. (1998). <A HREF=\"http://www.neb.com/rebase/\">REBASE</A>--restriction enzymes and methylases. Nucl. Ac. Res. 26, 338-350.");
			fprintf(out,"<LI><A HREF=\"http://locus.jouy.inra.fr/soft/cloneit/cloneit.html\">%s</A>",VAR_CITATION);
			fprintf(out,"</UL><LI>Acknowledgments");
			fprintf(out,"<P>I would like to thank Audrey Nepveu-de-Villemarceau, Janine L., <A HREF=\"http://duc.urbb.jussieu.fr/equipe.html\">Dr S. Hazout and his team</A>, Philippe Bessieres, Christine Young, and Maria Piron for their help.");
			fprintf(out,"</UL>");
			printTitle(out,mode,"",1,FALSE);
			fprintf(out,"<HR>%s by <AU><A HREF=\"mailto:lindenb@biotec.jouy.inra.fr\">Pierre LINDENBAUM</A></AU>.<BR>",VAR_VERSION);
			fputs("<ADRESS>Biologie Moleculaire des Rotavirus.<BR>\n",out);
			fputs("VIM <A HREF=\"http://jouy.inra.fr\">INRA</A>.<BR>\n",out);
			fputs("78350 JOUY-EN-JOSAS FRANCE.><BR><A HREF=\"http://www.jouy.inra.fr/index.html\"><IMG SRC=\"http://www.jouy.inra.fr/icons/accueil.gif\" ALT=\"Page d'accueil\" BORDER=0></A></ADRESS>\n",out);
			fputs("</BODY></HTML>\n",out);
			}break;
		default:break;
		}
	}
	}

void save_it(STRUCT_STRATEGY *Strategy, int CloningType)
	{
	/*****************************************************************************
	save_it save the sequence from the strategy found by the program. It saves a 
	sequence in FASTA text (ASCII) format.
		- Cloned insert into the vector
		- insert with in frame deletion
		- insert with frameshift
	******************************************************************************/
	FILE *tampon;
	char Nom_du_Fichier[MAX_NOM_FICHIER];
	int i,j,SeqA,SeqB,posV5,posV5I5,posI3V3,posV3,overhang_A,overhang_B;
	/*
	sequence will be write from 1 to posV5 (sequence n¡A)
	then from posV5I5 to posI3V3 (Sequence n¡B it should always be the INSERT)
	then from posV3 to npb[SeqA] (sequence n¡A)
	*/
	char *string_pol[4]={"?","No","modifying polymerase","modifying polymerase"};
	/**************************************************************************/
	DRAW_LINE;
	printf("Saving sequence is experimental in this release, please look carefully at your results.\nInput file's name to save:");
	#ifdef __MAC__
		tampon = MacWriteTextFile(Nom_du_Fichier,"\pCloneIt Sequence");
	#else
		fflush(stdin);
		Cmd_lire(Nom_du_Fichier);
		tampon = Fct_Write_File(Nom_du_Fichier);
	#endif
	if (tampon==NULL || strcmp(Nom_du_Fichier,"")==ARE_IDENTIC) 
		{
		BEEP;
		}
	else
		{
		printf("\n\t----> Saving File \"%s\":",Nom_du_Fichier);
		fprintf(tampon,">\n> File created with %s.\n> Pierre LINDENBAUM 1997\n> lindenb@biotec.jouy.inra.fr\n>\n",VAR_VERSION);
		/* if save_it is called by the function Sub_Cloning.c CloningType==1 **/
		if(CloningType==1)
			{
			/* write caracteristic of VECTOR and of the enzymes used */
			fprintf(tampon,">%s\n",FICHIER_ADN[VECTOR]);
			fprintf(tampon,">  %s [%s] (%d) %s treatment.\n",
										Enzymes[Sites[Strategy->Couple[SIDE_5].NumSite[VECTOR]].NumEnz].nom,
										Enzymes[Sites[Strategy->Couple[SIDE_5].NumSite[VECTOR]].NumEnz].site_complet,
										Sites[Strategy->Couple[SIDE_5].NumSite[VECTOR]].Loc,
										string_pol[Strategy->Couple[SIDE_5].Treatment[VECTOR]]);
			if(Strategy->Couple[SIDE_5].NumSite[VECTOR] != Strategy->Couple[SIDE_3].NumSite[VECTOR])
			fprintf(tampon,">  %s [%s] (%d) %s treatment.\n", 
										Enzymes[Sites[Strategy->Couple[SIDE_3].NumSite[VECTOR]].NumEnz].nom,
										Enzymes[Sites[Strategy->Couple[SIDE_3].NumSite[VECTOR]].NumEnz].site_complet,
										Sites[Strategy->Couple[SIDE_3].NumSite[VECTOR]].Loc,
										string_pol[Strategy->Couple[SIDE_3].Treatment[VECTOR]]);
			SeqA=VECTOR;
			}
		else /* if save_it is called by a deletion in frame resarch */
			SeqA=INSERT;
		/*
		if save_it is called by the function 
			Sub_Cloning.c CloningType==1 or DeltaFrame CloningType==2
		*/
		if(CloningType==1 || CloningType==2)
			{
			SeqB=INSERT;
			/* write caracteristic of INSERT and of the enzymes used */
			fprintf(tampon,">%s\n",FICHIER_ADN[INSERT]);
			
			fprintf(tampon,">  %s [%s] (%d) %s treatment.\n",
										Enzymes[Sites[Strategy->Couple[SIDE_5].NumSite[VECTOR]].NumEnz].nom,
										Enzymes[Sites[Strategy->Couple[SIDE_5].NumSite[VECTOR]].NumEnz].site_complet,
										Sites[Strategy->Couple[SIDE_5].NumSite[VECTOR]].Loc,
										string_pol[Strategy->Couple[SIDE_5].Treatment[VECTOR]]);
			
			if(Strategy->Couple[SIDE_5].NumSite[VECTOR] != Strategy->Couple[SIDE_3].NumSite[INSERT])
			fprintf(tampon,">  %s [%s] (%d) %s treatment.\n", 
										Enzymes[Sites[Strategy->Couple[SIDE_3].NumSite[INSERT]].NumEnz].nom,
										Enzymes[Sites[Strategy->Couple[SIDE_3].NumSite[INSERT]].NumEnz].site_complet,
										Sites[Strategy->Couple[SIDE_3].NumSite[INSERT]].Loc,
										string_pol[Strategy->Couple[SIDE_3].Treatment[INSERT]]);
			}
		/* if save_it is called by the function FrameShift CloningType==3 */
		if(CloningType ==3)
			{
			/* write caracteristic of INSERT and of the unique enzyme used */
			SeqB=INSERT;
			fprintf(tampon,">%s\n",FICHIER_ADN[INSERT]);
			fprintf(tampon,">  %s [%s] (%d) %s treatment (fill-in and ligation).\n",
										Enzymes[Sites[Strategy->Couple[SIDE_5].NumSite[INSERT]].NumEnz].nom,
										Enzymes[Sites[Strategy->Couple[SIDE_5].NumSite[INSERT]].NumEnz].site_complet,
										Sites[Strategy->Couple[SIDE_5].NumSite[INSERT]].Loc,
										string_pol[Strategy->Couple[SIDE_5].Treatment[INSERT]]);
			}
		/****************************************************************************************/
		/* if save_it is called by the function Sub_Cloning.c CloningType==1 **/
		if(CloningType==1)
			{
			overhang_A= Strategy->Couple[SIDE_5].Pos5_3[SeqA]-Strategy->Couple[SIDE_5].Pos3_5[SeqA];
			overhang_B= Strategy->Couple[SIDE_5].Pos5_3[SeqB]-Strategy->Couple[SIDE_5].Pos3_5[SeqB];
		
		/**********************************
		5------    -----------------------3
		3--------    ---------------------5
		or
		5--------    ---------------------3
		3------    -----------------------5
		**********************************/
			if( overhang_A == overhang_B && Sign(overhang_A)==Sign(overhang_B))
				{
				posV5  =Strategy->Couple[SIDE_5].Pos5_3[SeqA];
				posV5I5=Strategy->Couple[SIDE_5].Pos5_3[SeqB];
				}
			else if((overhang_A < overhang_B) && Sign(overhang_A)==Sign(overhang_B))
				{		
				posV5   = MIN(Strategy->Couple[SIDE_5].Pos5_3[SeqA],Strategy->Couple[SIDE_5].Pos3_5[SeqA]);
				posV5I5 = MIN(Strategy->Couple[SIDE_5].Pos5_3[SeqB],Strategy->Couple[SIDE_5].Pos3_5[SeqB]);
				}
			else if((overhang_A > overhang_B) && Sign(overhang_A)==Sign(overhang_B))
				{
				posV5   = MAX(Strategy->Couple[SIDE_5].Pos5_3[SeqA],Strategy->Couple[SIDE_5].Pos3_5[SeqA]);
				posV5I5 = MAX(Strategy->Couple[SIDE_5].Pos5_3[SeqB],Strategy->Couple[SIDE_5].Pos3_5[SeqB]);
				}
			/***************/
			overhang_A= Strategy->Couple[SIDE_3].Pos5_3[SeqB]-Strategy->Couple[SIDE_3].Pos3_5[SeqB];
			overhang_B= Strategy->Couple[SIDE_3].Pos5_3[SeqA]-Strategy->Couple[SIDE_3].Pos3_5[SeqA];
		/**********************************
		5------    -----------------------3
		3--------    ---------------------5
		or
		5--------    ---------------------3
		3------    -----------------------5
		**********************************/
			if( overhang_A == overhang_B && Sign(overhang_A)==Sign(overhang_B))
				{
				posI3V3 =Strategy->Couple[SIDE_3].Pos5_3[SeqB];
				posV3	=Strategy->Couple[SIDE_3].Pos5_3[SeqA];
				}
		/**********************************
		5------        -------------------3
		3--------                       --5
		or
		5--------                       --3
		3------        -------------------5
		**********************************/
			else if((overhang_A < overhang_B) && Sign(overhang_A)==Sign(overhang_B))
				{		
				posI3V3  = MIN(Strategy->Couple[SIDE_3].Pos5_3[SeqB],Strategy->Couple[SIDE_3].Pos3_5[SeqB]);
				posV3	 = MIN(Strategy->Couple[SIDE_3].Pos5_3[SeqA],Strategy->Couple[SIDE_3].Pos3_5[SeqA]);
				}
		/**********************************
		5-------------                   -3
		3-                             ---5
		or
		5-                             ---3
		3-------------                   -5
		**********************************/

			else if((overhang_A > overhang_B) && Sign(overhang_A)==Sign(overhang_B))
				{
				posI3V3 = MAX(Strategy->Couple[SIDE_3].Pos5_3[SeqB],Strategy->Couple[SIDE_3].Pos3_5[SeqB]);
				posV3	= MAX(Strategy->Couple[SIDE_3].Pos5_3[SeqA],Strategy->Couple[SIDE_3].Pos3_5[SeqA]);
				}
			}
		else if(CloningType==2)
			{
			/* if save_it is called by function DeltaFrame CloningType==2 */
			SeqA=VECTOR;
			SeqB=INSERT;
			overhang_A= Strategy->Couple[SIDE_5].Pos5_3[SeqA]-Strategy->Couple[SIDE_5].Pos3_5[SeqA];
			overhang_B= Strategy->Couple[SIDE_5].Pos5_3[SeqB]-Strategy->Couple[SIDE_5].Pos3_5[SeqB];
			
			posV5=0;
			posV5I5=1;
		/**********************************
		5------    -----------------------3
		3--------    ---------------------5
		or
		5--------    ---------------------3
		3------    -----------------------5
		**********************************/
			if( overhang_A == overhang_B && Sign(overhang_A)==Sign(overhang_B))
				{
				posI3V3 =Strategy->Couple[SIDE_5].Pos5_3[SeqA];
				posV3	=Strategy->Couple[SIDE_5].Pos5_3[SeqB];
				}
		/**********************************
		5------        -------------------3
		3--------                       --5
		or
		5--------                       --3
		3------        -------------------5
		**********************************/

			else if((overhang_A < overhang_B) && Sign(overhang_A)==Sign(overhang_B))
				{		
				posI3V3  = MIN(Strategy->Couple[SIDE_5].Pos5_3[SeqA],Strategy->Couple[SIDE_5].Pos3_5[SeqA]);
				posV3	 = MIN(Strategy->Couple[SIDE_5].Pos5_3[SeqB],Strategy->Couple[SIDE_5].Pos3_5[SeqB]);
				}
		/**********************************
		5-------------                   -3
		3-                             ---5
		or
		5-                             ---3
		3-------------                   -5
		**********************************/
			else if((overhang_A > overhang_B) && Sign(overhang_A)==Sign(overhang_B))
				{
				posI3V3 = MAX(Strategy->Couple[SIDE_5].Pos5_3[SeqA],Strategy->Couple[SIDE_5].Pos3_5[SeqA]);
				posV3	= MAX(Strategy->Couple[SIDE_5].Pos5_3[SeqB],Strategy->Couple[SIDE_5].Pos3_5[SeqB]);
				}
			SeqA=INSERT;
			}/* end cloning type 2*/
		else if(CloningType==3)
			{
			/* if save_it is called by the function FrameShift CloningType==3 */
			/*  5----...   -------3          5---    ..--3
			    3-------   ...--- 5          3---..    --5  */
			SeqA=INSERT;
			posV5=0;
			posV5I5=1;
			posI3V3= Sites[Strategy->Couple[SIDE_5].NumSite[INSERT]].Loc+
						Enzymes[Sites[Strategy->Couple[SIDE_5].NumSite[INSERT]].NumEnz].pos3_5;
			posV3=  Sites[Strategy->Couple[SIDE_5].NumSite[INSERT]].Loc+
						Enzymes[Sites[Strategy->Couple[SIDE_5].NumSite[INSERT]].NumEnz].pos5_3;
			}

		/*
		sequence will be write from 1 to posV5 (sequence n¡A)
		then from posV5I5 to posI3V3 (Sequence n¡B it should always be the INSERT)
		then from posV3 to npb[SeqA] (sequence n¡A)
		*/
		fputs(">\n>..\n",tampon);
		j=0;
	
	/*
	printf("je sauve de %d -%d -%d-%d\n",posV5,posV5I5,posI3V3,posV3);INKEY;
	printf("line %d:=%d\n",__LINE__,Strategy->Couple[SIDE_5].Pos5_3[SeqA]);
	printf("line %d:=%d\n",__LINE__,Strategy->Couple[SIDE_5].Pos5_3[SeqB]);
	printf("line %d:=%d\n",__LINE__,Strategy->Couple[SIDE_3].Pos5_3[SeqA]);
	printf("line %d:=%d\n",__LINE__,Strategy->Couple[SIDE_3].Pos5_3[SeqB]);
	printf("line %d:=%d\n",__LINE__,Strategy->Couple[SIDE_5].Pos3_5[SeqA]);
	printf("line %d:=%d\n",__LINE__,Strategy->Couple[SIDE_5].Pos3_5[SeqB]);
	printf("line %d:=%d\n",__LINE__,Strategy->Couple[SIDE_3].Pos3_5[SeqA]);
	printf("line %d:=%d\n",__LINE__,Strategy->Couple[SIDE_3].Pos3_5[SeqB]);
	*/

	
	
	
	/*	 if(posV5<0 || posV5>npb[SeqA]) printf("Error PosV5=%d\n",posV5);
	else if(posV5I5<0 || posV5I5>npb[INSERT]) printf("Error posV5I5=%d\n",posV5I5);
	else if(posI3V3<0 || posI3V3>npb[INSERT]) printf("Error posI3V3=%d\n",posI3V3);
	else if(posV3<0 || posV3>npb[SeqA]) printf("Error posV3=%d\n",posV3);
	else*/
			
			for(i=1;i<posV5;i++)
				{
				if(Fct_Frame(SeqA,i)==IS_IN_FRAME) fputs(" ",tampon);
				fprintf(tampon , "%c",UPPER(sequence[SeqA][i]));
				if(++j>41) {j=0;fputs("\n",tampon);}
				}
			for(i=posV5I5;i<posI3V3;i++)
				{
				if(Fct_Frame(INSERT,i)==IS_IN_FRAME) fputs(" ",tampon);
				fprintf(tampon , "%c",LOWER(sequence[INSERT][i]));
				if(++j>41) {j=0;fputs("\n",tampon);}
				}
			for(i=posV3;i<=npb[SeqA];i++)
				{
				if(Fct_Frame(SeqA,i)==IS_IN_FRAME) fputs(" ",tampon);
				fprintf(tampon , "%c",UPPER(sequence[SeqA][i]));
				if(++j>41) {j=0;fputs("\n",tampon);}
				}
			fputs("\n",tampon);
			
		fclose(tampon);
		printf(" Done.\n");
		INKEY;
		}
	
	}



/* File 'Sequence.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include <stdlib.h>
	#include <errno.h>
	#include"ADN.h"
#endif

void ChooseATG(int NumSeq)
	{
	boolean vara=FALSE;
	int varc,i,j;
	if(npb[NumSeq]!=0)
		{
		Get_ATG(NumSeq);
		for(i=0;i<=2;i++)
			{
			printf("FRAME =%4d:\n",pos_ATG[NumSeq]+i);
			printf("------------\n",pos_ATG[NumSeq]+i);
			for(j=pos_ATG[NumSeq];j<=var_max[NumSeq];j+=3)
				{
				printf("%c",Translation_at(NumSeq,j+i));
				}
			printf("\n\n");
			}
		
		
		while(vara==FALSE)
			{
			printf("\nGive me a position that is in frame\nin %s sequence \"%s\" /* 1<%d<%d */:"
				,(NumSeq==VECTOR?"VECTOR":"INSERT"),FICHIER_ADN[NumSeq],pos_ATG[NumSeq],npb[NumSeq]);
			i = Get_Number(pos_ATG[NumSeq]);
			if(i==0 || i<1 || i>npb[NumSeq])
					{BEEP;continue;}
				else
					{pos_ATG[NumSeq]=i;vara=TRUE;}
			}
		}
	else
		BEEP;
		
	}



/*******************************************************************************
 get the most probable atg frame (0,1 or 2) in the sequence
		that is to say the length max of sequence where there is no stop codon
********************************************************************************/
int Get_ATG(int NumSeq)
	{
	int i,j,vara,_max=0;
	
	/* no boundaries are defined */
	if(var_min[NumSeq]==var_max[NumSeq]) {BEEP;return(0);}
	pos_ATG[NumSeq]=FALSE;

		if(NumSeq==INSERT)
			{
			/* if the sequence is insert, the program is checking the frame, where there is
				the longest ORF without stop codon IN the insert */
			for(i=0;i<=2;i++)
				{
				vara=0;
				for(j=var_min[INSERT]+i;j<=var_max[INSERT];j=j+3)
					{
					if(Is_Codon_Stop(Fct_Traduction(sequence[INSERT][Fct_Pos(INSERT,j)],
													sequence[INSERT][Fct_Pos(INSERT,j+1)],
													sequence[INSERT][Fct_Pos(INSERT,j+2)]))!=TRUE)
						{
						vara++;
					
						if(vara>_max)
							{
							_max=vara;
							pos_ATG[INSERT]=var_min[INSERT]+i;
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
				for(j=var_min[VECTOR]+i;j>=1;j=j-3)
					{
					if(Is_Codon_Stop(Fct_Traduction(sequence[VECTOR][Fct_Pos(VECTOR,j)],
													sequence[VECTOR][Fct_Pos(VECTOR,j+1)],
													sequence[VECTOR][Fct_Pos(VECTOR,j+2)]))==TRUE)
						break;
					else
						vara++;
					}
				if(vara>_max)
					{
					
					_max=vara;
					pos_ATG[VECTOR]=var_min[VECTOR]+i;
					}
				/* search on the right of cloning box */
				vara=0;
				for(j=var_max[VECTOR]+i;j<=npb[VECTOR];j=j+3)
					{
					if(Is_Codon_Stop(Fct_Traduction(sequence[VECTOR][Fct_Pos(VECTOR,j)],
													sequence[VECTOR][Fct_Pos(VECTOR,j+1)],
													sequence[VECTOR][Fct_Pos(VECTOR,j+2)]))==TRUE)
						break;
					else
						vara++;
					}
				if(vara>_max)
					{
					_max=vara;
					pos_ATG[VECTOR]=var_max[VECTOR]+i;
					}
				/* search IN the cloning box */
				vara=0;
				for(j=var_min[VECTOR]+i;j<=var_max[VECTOR];j=j+3)
					{
					if(Is_Codon_Stop(Fct_Traduction(sequence[VECTOR][Fct_Pos(VECTOR,j)],
													sequence[VECTOR][Fct_Pos(VECTOR,j+1)],
													sequence[VECTOR][Fct_Pos(VECTOR,j+2)]))==TRUE)
						break;
					else
						vara++;
					}
				if(vara>_max)
					{
					
					_max=vara;
					pos_ATG[VECTOR]=var_min[VECTOR]+i;
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
		
	p_seq=&sequence[NumSeq][1];
	result=strstr(p_seq, mot2);
	if(result!=NULL)
		{
		vara=FALSE;
		poly_left=npb[NumSeq]-strlen(result)+1+strlen(mot2);
		}
	if(vara==FALSE)
		{
		p_seq=&sequence[NumSeq][poly_left];
		result=strstr(p_seq, mot3);
		if(result!=NULL)
			{
			poly_right=npb[NumSeq]-strlen(result)+1;
			}
		}
	if(poly_left<poly_right && poly_left!=0 && poly_right!=0)
		{
		var_return=TRUE;
		var_min[NumSeq]=poly_left;
		var_max[NumSeq]=poly_right;
		if(Preference.display_messages==TRUE)
			{
			printf("\t  Your sequence seems to be a \"%s\" type plasmid.\n",mot1);
			printf("\t  The cloning boxes may be localised between %d and %d.\n",var_min[NumSeq],var_max[NumSeq]);
			Menu(12);
			}
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
	
	#ifdef VAR_UNIX
	if(getenv("HOME")!=NULL)
		{
		sprintf(ThePolyfFile,"%s/%s",getenv("HOME"),POLYLINKER_FILE);
		}
	else
		sprintf(ThePolyfFile,"%s",POLYLINKER_FILE);
	#else
		sprintf(ThePolyfFile,"%s",POLYLINKER_FILE);
	#endif
		
	
	
	Search_Done=FALSE;
	/* open POLYLINKER_FILE */
	tampon = fopen(ThePolyfFile, "r" );
	if (tampon==NULL) 
		{
		printf("Creating default Polylinkers File \"%s\"...\n",POLYLINKER_FILE);
		if((tampon = fopen(ThePolyfFile, "w" ))!=NULL)
			{
			fprintf(tampon,";\n;%s Pierre LINDENBAUM 1998.\n",VAR_VERSION);
			fprintf(tampon,"; List of Polylinkers (default).\n");
			fprintf(tampon,";\n;Use the oligo that you would have used to amplify your gene:\n");
			fprintf(tampon,";Syntax: Name,Direct,Reverse,\n;\n");
			fprintf(tampon,"pFASTBAC,TATTCCGGATTATTCATACC,GATTATGATCCTCTAGTACTTCTCGAC,\n");
			fprintf(tampon,"Baculo,ttttactgttttcgtaacagtttt,cggatttccttgaagagagta,\n");
			fprintf(tampon,"pBK,ggtctatataagcagagctggt,acaggaaacagctatgaccttg,\n");
			fprintf(tampon,"GEX,atcgaaggtcg,tcagtcagtcacgatg,\n");
			fprintf(tampon,"T7/T3,aattaaccctcactaaaggg,taatacgactcactataggg,\n");
			fprintf(tampon,"pET25B,TAATACGACTCACTATA,CCCGTTTAGAGGCCCCAAGGGGTTA,\n");
			fprintf(tampon,"pIIIMS2-1,ttccggctagaactagtggatcc,tcgactctagaggatcg,\n");
			fprintf(tampon,"pIIIMS2-2,agagtcgacctgcaggcatgcaagctg,gctagaactagtggatcc,\n");
			fprintf(tampon,"pGBT9,CAGTTGACTGTATCGCCG,GCCCGGAATTAGCTTGG,\n");
			fprintf(tampon,"pGAD424,CCAAAAAAAGAGATC,TTCAGTATCTACGATTCAT,\n");
			fprintf(tampon,"pGADGL,CCAAAAAAAGAGATC,ACTATAGGGCGAATTGG,\n");
			fprintf(tampon,"pcDNAFLAG,atggactacaaggacgacgatgacaa,cttggtaccgagctcggatcc,\n");
			fprintf(tampon,"pcDNA3,CACTATAGGGAGACCC,AGGTGACACTATAGAATA,\n");
			fprintf(tampon,"pYX213,TAACGTCAAGGAGAAAAAACCCCGGAT,GAAAAACGTTCATTGTTCCTTAT,\n");
			fclose(tampon);
			if((tampon = fopen(ThePolyfFile, "r" ))==NULL) ERROR_USER;
			}
		}
	if (tampon==NULL) 
		{
		printf("\n---> File error :Can't find '%s' file.\n",POLYLINKER_FILE);
		Handle_error(TRUE);INKEY;Menu(11);Menu(0);Menu(19);
		}
	else
		{
		while(c != EOF )
			{
			/* read file */
			c=fgetc(tampon);
			if(c==';')
				{while((c=fgetc(tampon))!='\n' && c!='\r' && c!=21 && c!=8);vara=0;varb=0;} /* discard commentary */
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
						printf("Beware: polylinker %s should not contain degenerate bases.",mot1);
						INKEY;break;
						}
					if (vara==0) {mot1[varb++]=c;mot1[varb]='\0';} /* complete name of plasmid */
					if (vara==1 && Fct_est_ADN(c)>0 && varb<MAX_LENGHT_POLY)
						{mot2[varb++]=UPPER(c);mot2[varb]='\0';} /* complete oligo1 */
					if (vara==2 && Fct_est_ADN(c)>0 && varb<MAX_LENGHT_POLY)
						{mot3[varb++]=UPPER(c);mot3[varb]='\0';} /* complete oligo1 */
					}
				else if(varb>=MAX_LENGHT_POLY)
					{
					BEEP;
					printf("Reduce your oligonucleotide length (%s).\n",mot1);
					}
				}
			}
		fclose(tampon);
		}
	if(var_min[NumSeq]==var_max[NumSeq]) /* if polylinker is not defined set boundaries to all the sequence*/
			{
			var_min[NumSeq]=1;
			var_max[NumSeq]=npb[NumSeq];
			BEEP;
			BEEP;
fprintf(stderr,"%s can't find the cloning box(es). Please define the ",VAR_VERSION); 
fprintf(stderr,"boundaries. You should also refine your \"%s\" polylinkers file.\n",POLYLINKER_FILE);
Menu(12);fflush(stdin);
			}
	if(NumSeq==INSERT) /* if INSERT, define the internal cloning boxes by default */
		{
		sigma_5=(int)((var_max[INSERT]-var_min[INSERT])/Preference.var_pct_min);
		sigma_3=(int)((var_max[INSERT]-var_min[INSERT])/Preference.var_pct_min);
		printf("The INSERT internal cloning boxes were randomly defined (%d percent of the insert).\n",Preference.var_pct_min);
		printf("You should check the displayed datas...\n");
		}
	}

boolean AntiParallele(int NumSeq)
	{
	int i;
	char c;
	Search_Done=FALSE;
	
	for(i=1;i<20;i++) if(i<=npb[NumSeq]) printf("%c",sequence[NumSeq][i]);
	if(npb[NumSeq]>20)
		{
		printf("...");
		for(i=npb[NumSeq]-20;i<=npb[NumSeq];i++) printf("%c",sequence[NumSeq][i]);
		}
	printf("\n\n");
	printf("processing...\n\n");
	for(i=1;i<=1+(int)((double)npb[NumSeq]/2.0);i++)
		{
		c=sequence[NumSeq][i]; /* i=1->c=npb */ /* i=npb c=0 */
		sequence[NumSeq][i]=Fct_Complementaire(sequence[NumSeq][npb[NumSeq]-i+1]);
		sequence[NumSeq][npb[NumSeq]-i+1]=Fct_Complementaire(c);
		}
	i=var_min[NumSeq];
	var_min[NumSeq]=npb[NumSeq]-var_max[NumSeq]+1;
	var_max[NumSeq]=npb[NumSeq]-i+1;
	if(NumSeq==INSERT)
		{
		i=sigma_5;
		sigma_5=sigma_3;
		sigma_3=i;
		}
	for(i=1;i<20;i++) if(i<=npb[NumSeq]) printf("%c",sequence[NumSeq][i]);
	if(npb[NumSeq]>20)
		{
		printf("...");
		for(i=npb[NumSeq]-20;i<=npb[NumSeq];i++) printf("%c",sequence[NumSeq][i]);
		}
	printf("\n");
	Get_ATG(NumSeq);
	return(TRUE);
	}



#ifndef __MAC__
/********************************  Cmd_LOAD_SEQUENCE  ***/
boolean Cmd_LOAD_SEQUENCE(int NumSeq)
	{
	char c;
	FILE *tampon;
	int degenerate=FALSE;
	/* intit sequence parameters */
	npb[NumSeq]=0;
	pos_ATG[NumSeq]=FALSE;
	var_min[NumSeq]=0;
	var_max[NumSeq]=0;
	IsStriderSeq[NumSeq]=FALSE;
	
	if(ReadStriderFormat(NumSeq)==FALSE)
		{
		if((tampon=fopen(FICHIER_ADN[NumSeq],"r"))==NULL) 
			{
			printf("/* File error : %s */\n",FICHIER_ADN[NumSeq]);
			Handle_error(TRUE);
			strcpy(FICHIER_ADN[NumSeq],"");
			Menu(12);
			return(FALSE);
			}
		else
			{
			rewind(tampon);
			while( (c=fgetc(tampon)) != EOF )
				{
				if (c=='>'|| c==';')
					{
					while((c=fgetc(tampon)!='\n') && c!='\r' && c!=21 && c!=8)
						npb[NumSeq]=0;
					if(c==EOF) {printf("/* File error : %s */\n",FICHIER_ADN[NumSeq]);
								BEEP;
								Menu(19);
								break;}
					}
				if (c!='\n' && c!='/' && c!=' ' && c!='\t')
				if(c<48 || c>57)
					if(Fct_est_ADN(c)>=1)
						{
						npb[NumSeq]++;
						if(npb[NumSeq]<MAX_NPB)
							{
							sequence[NumSeq][npb[NumSeq]]=UPPER(c);
							if(Fct_est_ADN(c)>=AMBIGOUS) degenerate=TRUE;
							}
						else
							{
							printf("Sequence too large! [n>%d]\n",MAX_NPB);Display("Out of memory !",3);
							BEEP;Menu(19);
							}	
						}
				sequence[NumSeq][npb[NumSeq]+1]='\0';
				}
			fclose(tampon);
			printf( "\n\t\t\t%d bp in '%s'.\n\n",npb[NumSeq],FICHIER_ADN[NumSeq]);
			if(degenerate==TRUE)
				{
				BEEP;
				Display("SORRY: This sequence contains degenerate bases",3);
				Menu(12);
				strcpy(FICHIER_ADN[NumSeq],"");
				npb[NumSeq]=0;
				}
			}
		}/* non strider format  */
		
	if(npb[NumSeq]>0)
		{
		Cmd_POLYLINKER1(NumSeq);
		Get_ATG(NumSeq);
		}
	return(npb[NumSeq]>0);
	}

#endif



/* File 'Compte_Enzyme.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <ctype.h>
	#include <errno.h>
	#include <ctype.h>
	#include "ADN.h"
	#ifdef __MAC__
		#include <console.h>
	#endif
#endif

void printEntreprise(short mode, FILE *out,char c)
	{
	if(mode!=FORMAT_HTML)
		{
		switch(c)
			{
			case('A'):fprintf(out,"\tAmersham Life Sciences-USB\n");break;
			case('B'):fprintf(out,"\tLife Technologies Inc, Gibco-BRL\n");break;
			case('C'):fprintf(out,"\tMinotech Molecular Biology Products\n");break;
			case('D'):fprintf(out,"\tAngewandte Gentechnologie Systeme\n");break;
			case('E'):fprintf(out,"\tStratagene\n");break;
			case('F'):fprintf(out,"\tFermentas AB\n");break;
			case('G'):fprintf(out,"\tAppligene Oncor\n");break;
			case('H'):fprintf(out,"\tAmerican Allied Biochemical, Inc.\n");break;
			case('I'):fprintf(out,"\tSibEnzyme Ltd.\n");break;
			case('J'):fprintf(out,"\tNippon Gene Co., Ltd.\n");break;
			case('K'):fprintf(out,"\tTakara Shuzo Co. Ltd.\n");break;
			case('L'):fprintf(out,"\tKramel Biotech\n");break;
			case('M'):fprintf(out,"\tBoehringer-Mannheim\n");break;
			case('N'):fprintf(out,"\tNew England BioLabs\n");break;
			case('O'):fprintf(out,"\tToyobo Biochemicals\n");break;
			case('P'):fprintf(out,"\tPharmacia Biotech Inc.\n");break;
			case('Q'):fprintf(out,"\tCHIMERx\n");break;
			case('R'):fprintf(out,"\tPromega Corporation\n");break;
			case('S'):fprintf(out,"\tSigma\n");break;
			case('T'):fprintf(out,"\tAdvanced Biotechnologies Ltd.\n");break;
			default:fprintf(out,"\t? Unknow '%c'.\n",c);break;
			}
		}
	else
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
	}

int Cmd_Get_Info_Choice(int NumEnz1,int NumEnz2,int NumEnz3,int NumEnz4)
	{
	int MyChoice[4]={0,0,0,0};
	int NumChoice=0,varb=0;
	
	
	fflush(stdin);
	if( strcmp(Enzymes[Sites[NumEnz1].NumEnz].nom,Enzymes[Sites[NumEnz2].NumEnz].nom)==ARE_IDENTIC &&
		strcmp(Enzymes[Sites[NumEnz2].NumEnz].nom,Enzymes[Sites[NumEnz3].NumEnz].nom)==ARE_IDENTIC &&
		strcmp(Enzymes[Sites[NumEnz3].NumEnz].nom,Enzymes[Sites[NumEnz4].NumEnz].nom)==ARE_IDENTIC)
			{Cmd_Get_Info(FORMAT_SCREEN,stdout,NumEnz1);return(TRUE);}
	
	 /*********************************/	
	/* Display Buffer and T¡C        */
   /*********************************/
 	/*Display_Temp_Buffer(stdout,FORMAT_SCREEN,Sites[NumEnz1].NumEnz,Sites[NumEnz2].NumEnz);
	Display_Temp_Buffer(stdout,FORMAT_SCREEN,Sites[NumEnz3].NumEnz,Sites[NumEnz4].NumEnz);*/

	
	
	
	Display("Get informations about:",1);
	
	printf("\t1 %s [%s]\n",Enzymes[Sites[NumEnz1].NumEnz].nom,Enzymes[Sites[NumEnz1].NumEnz].site_complet);
	MyChoice[0]=NumEnz1;
	NumChoice++;
	
	if( strcmp(Enzymes[Sites[NumEnz1].NumEnz].nom,Enzymes[Sites[NumEnz2].NumEnz].nom)!=ARE_IDENTIC)
		{
		printf("\t%d %s [%s]\n",NumChoice+1,Enzymes[Sites[NumEnz2].NumEnz].nom,Enzymes[Sites[NumEnz2].NumEnz].site_complet);
		MyChoice[NumChoice++]=NumEnz2;
		}
	if( strcmp(Enzymes[Sites[NumEnz1].NumEnz].nom,Enzymes[Sites[NumEnz3].NumEnz].nom)!=ARE_IDENTIC &&
		strcmp(Enzymes[Sites[NumEnz2].NumEnz].nom,Enzymes[Sites[NumEnz3].NumEnz].nom)!=ARE_IDENTIC)
		{
		printf("\t%d %s [%s]\n",NumChoice+1,Enzymes[Sites[NumEnz3].NumEnz].nom,Enzymes[Sites[NumEnz3].NumEnz].site_complet);
		MyChoice[NumChoice++]=NumEnz3;
		}
	if( strcmp(Enzymes[Sites[NumEnz1].NumEnz].nom,Enzymes[Sites[NumEnz4].NumEnz].nom)!=ARE_IDENTIC &&
		strcmp(Enzymes[Sites[NumEnz2].NumEnz].nom,Enzymes[Sites[NumEnz4].NumEnz].nom)!=ARE_IDENTIC &&
		strcmp(Enzymes[Sites[NumEnz3].NumEnz].nom,Enzymes[Sites[NumEnz4].NumEnz].nom)!=ARE_IDENTIC)
		{
		printf("\t%d %s [%s]\n",NumChoice+1,Enzymes[Sites[NumEnz4].NumEnz].nom,Enzymes[Sites[NumEnz4].NumEnz].site_complet);
		MyChoice[NumChoice++]=NumEnz4;
		}
	
	while(varb<1 || varb>NumChoice)
		{
		#ifndef __MAC__
			printf("\n  Your choice ? /*1 < 1 < %d*/:",NumChoice);
			varb = Get_Number(1);
		#else
			while((varb=((int)get_char()-'0'))<1 || varb>NumChoice) /* do nothing */;
			INKEY;
		#endif
		}
	Cmd_Get_Info(FORMAT_SCREEN,stdout,MyChoice[varb-1]);
	return(TRUE);
	}
/************************************************************************************/
void printWebEnzyme(char *name,FILE *out)
	{
	int i=0;
	fprintf(out,"<A HREF=\"http://www.neb.com/rebase/enz/");
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
	if(mode!=FORMAT_HTML)
		{
		DRAW_HR(out,'.');
		fprintf(out,"Informations about %s [%s].\n",Enzymes[Sites[NumEnz].NumEnz].nom,Enzymes[Sites[NumEnz].NumEnz].site_complet);
		}
	else
		fprintf(out,"<LI><H2><A NAME=\"TAG%d\"> Informations about %s [%s].</H2>\n",
			Sites[NumEnz].NumEnz,
			Enzymes[Sites[NumEnz].NumEnz].nom,Enzymes[Sites[NumEnz].NumEnz].site_complet);
	
	Display_Temp_Buffer(out,mode,Sites[NumEnz].NumEnz,Sites[NumEnz].NumEnz);
	
	/**Create a link to the rebase www site***/
	fprintf(out,"Internet WWW link:");
	 printWebEnzyme(Enzymes[Sites[NumEnz].NumEnz].nom,out);
	if(mode!=FORMAT_HTML)  DRAW_HR(out,'.');
	if(npb[VECTOR]>0 && Fct_N_sites(Sites[NumEnz].NumEnz,VECTOR,1,npb[VECTOR],TRUE)) {fprintf(out,"VECTOR Pattern ");Digest(mode,out,VECTOR,NumEnz,NumEnz,FALSE);}
	if(npb[INSERT]>0 && Fct_N_sites(Sites[NumEnz].NumEnz,INSERT,1,npb[INSERT],TRUE)) {fprintf(out,"INSERT Pattern ");Digest(mode,out,INSERT,NumEnz,NumEnz,FALSE);}
	if(mode!=FORMAT_HTML) 
		DRAW_HR(out,'.');
	
	Cmd_Get_Info2(mode,out,Sites[NumEnz].NumEnz);
	if(mode==FORMAT_SCREEN)
		Menu(12);
	}


int Cmd_Get_Info2(short mode, FILE *out,int NumEnz)
{
/************************************************************************************
This procedure gets informations about an Enzyme from the Rebase File...
____________________________________________________________________________________*/
/******************************************************/
FILE	*Buffer;
char	c2,mot1[NOM_MAX_ENZYME],mot2[NOM_MAX_ENZYME],mot3[MAX_CORP]="",mot4,mot5,mot6[NBR_BUFFER];
int		i=0,var_find=0,fric=0,find_iso=FALSE,k;
char	line[MAX_LINE];
/******************************************************/

  /*******************************************/
 /* first find information about the enzyme */
/*******************************************/

if ((Buffer = fopen(FICHIER_ENZYME, "r" ))==NULL)
	{
	printf("\n\t----> File error  : '%s'.\n",FICHIER_ENZYME);
	Handle_error(TRUE);
	INKEY;
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
				case(REBASE_TEMPERATURE):fprintf(out," Temperature              :%s\n",&line[3]);break;
				case(REBASE_BUFFER):fprintf(out," Buffer                   :%s\n",&line[3]);break;
				case(REBASE_BUFFERS):fprintf(out," Activity in four buffers :%s\n",&line[3]);break;
				default:break;
				}
			if(c2=='7')
				{
				if(mode==FORMAT_HTML)
					fprintf(out,"<UL>");
				for(k=3;k<strlen(line);k++)
					printEntreprise(mode,out,line[k]);
				if(mode==FORMAT_HTML)
					fprintf(out,"</UL>");
				fprintf(out,"\n");
				}
			if(c2=='8') break;
			}/* end var_find */
		}/*end while*/
	/**************************************
	Now search isoschyzomers in Rebase File
	***************************************/
	if(mode!=FORMAT_HTML) DRAW_HR(out,'.');
	fprintf(out,"Looking for Isoschizomers.\n");
	if(mode==FORMAT_HTML) fprintf(out,"<UL>");
	find_iso=FALSE;
	fseek(Buffer,0L,SEEK_SET);/* redo from start */
	while(freadRebase(Buffer,mot1,mot2,mot3,&mot4,&mot5,&mot6[0])==TRUE)
		{
		if(strcmp(mot1,Enzymes[NumEnz].nom)==ARE_IDENTIC) continue;
		if(strcmp(mot2,Enzymes[NumEnz].site_complet)!=ARE_IDENTIC)  continue;
		if(strlen(mot3)<=0) continue;
		find_iso=TRUE;
		fprintf(out,"%s%s  %s.\n",(mode!=FORMAT_HTML?"\t":"<LI>"),mot1,mot2);

		if(strlen(mot3)>0)
			{
			fprintf(out,"%salso available at:\n",(mode!=FORMAT_HTML?"\t     ":"<LI>"));
			if(mode==FORMAT_HTML) fprintf(out,"<UL>");

			for(fric=0;fric<strlen(mot3);fric++)
				{
				if(mode!=FORMAT_HTML)
					fprintf(out,"\t\t");
				printEntreprise(mode,out,mot3[fric]);
				}
			if(mode==FORMAT_HTML) fprintf(out,"</UL>");
			}
		}
	if(find_iso==FALSE)
		fprintf(out,"\tNo other enzyme was found.\n");
	if(mode==FORMAT_HTML) fprintf(out,"</UL>");
	fclose(Buffer);
	}/* end if file */
if(mode!=FORMAT_HTML)  DRAW_HR(out,'.');
return(TRUE);
}

/************************************************************************************/
/* search an information (user's query ) about an enzyme */
void search_info(void)
	{
	/*===========================================================*/
	char	Word[MAX_NOM_FICHIER];
	int		var_find=FALSE,var_find_in_word,i,j,k,var_len_overhang=0,var_len=0,NumElem=FALSE;
	boolean	var_research_overhang,var_search_error;
	/*===========================================================*/
	DRAW_LINE;
	if(Preference.memory==TRUE)
		printf("Beware, short sites have been discarded.\n");
	Menu(17);
	printf("\nINPUT YOUR QUERY: \tSearch in \"%s\" for :",FICHIER_ENZYME);
	Cmd_lire(Word);
	PARAGRAF;
	if(strcmp(Word,"")==ARE_IDENTIC)
		{BEEP;CLS;}
	else
		{
		var_research_overhang=FALSE;
		var_search_error=FALSE;
		var_find_in_word=TRUE;
		
		/* init the field select */
		for(i=0;i<nbr_enzyme;i++)
			Enzymes[i].select[VECTOR]=FALSE;
		/*****************************************/
		/*is there a research about the length ? */
		/*****************************************/
		if(Word[0]=='$')
			{
			var_find_in_word=FALSE;
			var_research_overhang=TRUE;
			var_len=0;
			j=1;
			while(Word[j]!='\0' && Word[j]!='|' && Word[j]!='\n' && Word[j]!='_')
				{
				if(isdigit(Word[j])!=FALSE)
					var_len=var_len*10+Word[j]-'0';
				else {var_search_error=TRUE;break;}
				j++;
				}
			if(var_search_error==FALSE)
				for(i=0;i<nbr_enzyme;i++)
					if(Enzymes[i].taille_site==var_len)
						{Enzymes[i].select[VECTOR]=TRUE;NumElem=TRUE;}
						
			printf("Looking for an enzyme with a length equal to %d.\n",var_len);
			if(Word[j]=='\0') strcpy(Word,"!´#±");
			
			/* after, is there a research about the overhang ? */
			if(Word[j]=='_') 
				{
				k=j;
				var_find_in_word=TRUE;
				while(Word[j]!='\0' && Word[j]!='\n')
					{
					Word[j-k]=Word[j];Word[j-k+1]='\0';
					j++;
					}
				}
			/* after, is there a research about the name ? */
			else if(Word[j]=='|') 
			  {
			  k=j+1;
			  var_find_in_word=TRUE;
				while(Word[++j]!='\0' && Word[j]!='\n')
					{
					Word[j-k]=Word[j];Word[j-k+1]='\0';
					}
			  }
			if(Word[0]=='|') var_find_in_word=TRUE;
			}
		else
			for(i=0;i<nbr_enzyme;i++)
					Enzymes[i].select[VECTOR]=TRUE;

		/********************************************/
		/* is there a research about the overhang ? */
		/********************************************/
		if(Word[0]=='_' && var_search_error==FALSE)
			{
			var_find_in_word=FALSE;
			var_research_overhang=TRUE;
			switch(Word[1])
				{
				/* Blunt */
				case('B'):
						printf("Looking for BLUNT enzymes ");
						for(i=0;i<nbr_enzyme;i++)
							if(Fct_type(Enzymes[i])!=TYPE_BLUNT)
								Enzymes[i].select[VECTOR]=FALSE;
						break;
				/* 3' over */
				case('+'):
						printf("Looking for 3' overhanged enzymes ");
						for(i=0;i<nbr_enzyme;i++)
							if(Fct_type(Enzymes[i])!=TYPE_3_OVER)
								Enzymes[i].select[VECTOR]=FALSE;
						break;
				/* 5' over */
				case('-'):
						printf("Looking for 5' overhanged enzymes ");
						for(i=0;i<nbr_enzyme;i++)
							if(Fct_type(Enzymes[i])!=TYPE_5_OVER)
								Enzymes[i].select[VECTOR]=FALSE;
						break;
				default:{
						/* the user just want to look after length ? */
						if(isdigit(Word[1])!=FALSE)
							{
							printf("Looking for enzymes ");
							var_len_overhang=Word[1]-'0';
							for(i=0;i<nbr_enzyme;i++)
								Enzymes[i].select[VECTOR]=(Enzymes[i].select[VECTOR]==TRUE?TRUE:FALSE);
							}
						else
							var_search_error=TRUE;
						}break;
				}

			j=2;
			if(var_search_error==FALSE && Word[2]!='\0' && Word[2]!='|' && Word[2]!='\n')
				{
				/* what is the length of the overhang searched ? */ 
				while(Word[j]!='\0' && Word[j]!='|' && Word[j]!='\n')
					{
					if(isdigit(Word[j])!=FALSE)
						var_len_overhang=var_len_overhang*10+Word[j]-'0';
					else {var_search_error=TRUE;break;}
					j++;
					}
				}
								
			if(var_search_error==FALSE && var_len_overhang>0)
				{
				printf(" with an overhang length equal to %d ",var_len_overhang);
				/* discard all enzymes that have not this overhang length */
				for(i=0;i<nbr_enzyme;i++)
					if(abs(Enzymes[i].pos5_3-Enzymes[i].pos3_5)!=var_len_overhang)
					   		Enzymes[i].select[VECTOR]=FALSE;
				}
				
				
			/* does the research is looking for a substring ? ex: _B5|gaattc */
			if(var_search_error==FALSE)
				{
				k=j;
				if(Word[j]=='\0') strcpy(Word,"!´#±");
				if(Word[j]=='|') 
					{
					while(Word[++j]!='\0' && Word[j]!='\n')
						{
						Word[j-k-1]=Word[j];Word[j-k]='\0';
						}
					printf(" and containing the word \"%s\"",Word);
					var_find_in_word=TRUE;
					}
				}
			printf(".\n\n");
			for(i=0;i<nbr_enzyme;i++)
				if(Enzymes[i].select[VECTOR]==TRUE)	NumElem=TRUE;
			}/* end if Word[0] */
		if(var_search_error==TRUE)
			{BEEP;Display("Syntax error",2);strcpy(Word,"!´#±");NumElem=FALSE;}
		else
		if(var_find_in_word==TRUE)
		for(i=0;i<nbr_enzyme;i++)
			{
			var_find_in_word=FALSE;
			  /********************/
			 /* look in the name */
			/********************/
			if(strlen(Enzymes[i].nom) >= strlen(Word))
				{
				for(k=0;k<=strlen(Enzymes[i].nom)-strlen(Word);k++)
					{
					var_find=TRUE;
					for(j=0;j<strlen(Word);j++)
						{
						if(UPPER(Enzymes[i].nom[k+j])!=UPPER(Word[j]))
							{var_find=FALSE;break;}
						}
					if(var_find==TRUE)
						{var_find_in_word=TRUE;break;}
					}
				}
	  		  /*****************************/
			 /* look in the .site_complet */
			/*****************************/
			if(strlen(Enzymes[i].site_complet) >= strlen(Word))
				{
				for(k=0;k<=strlen(Enzymes[i].site_complet)-strlen(Word);k++)
					{
					var_find=TRUE;
					for(j=0;j<strlen(Word);j++)
						{
						if(UPPER(Enzymes[i].site_complet[k+j])!=UPPER(Word[j]))
							{var_find=FALSE;break;}
						}
					if(var_find==TRUE)
						{var_find_in_word=TRUE;break;}
					}
				}
			  /*********************/
			 /* look in the .site */
			/*********************/
			if(strlen(Enzymes[i].site) >= strlen(Word))
				{
				for(k=0;k<=strlen(Enzymes[i].site)-strlen(Word);k++)
					{
					var_find=TRUE;
					for(j=0;j<strlen(Word);j++)
						{
						if(UPPER(Enzymes[i].site[k+j])!=UPPER(Word[j]))
							{var_find=FALSE;break;}
						}
					if(var_find==TRUE)
						{var_find_in_word=TRUE;break;}
					}
				}
			if(var_find_in_word==TRUE)
				{Enzymes[i].select[VECTOR]=(var_research_overhang==TRUE?Enzymes[i].select[VECTOR]:TRUE);NumElem=TRUE;}
			else
				{Enzymes[i].select[VECTOR]=FALSE;}
			}/* end for */
		/***************************************************************/
		for(i=0;i<nbr_enzyme;i++)
				if(Enzymes[i].select[VECTOR]==TRUE)	NumElem=TRUE;
				
		if(NumElem==TRUE)
			{
			for(i=0;i<nbr_enzyme;i++)
			if(Enzymes[i].palindromic!=TRUE && Enzymes[i].select[VECTOR]==TRUE)
				{
				/* is 'i' is the same  enzyme than previous i-1 */
				if(i+1<nbr_enzyme)
					{
					if(strcmp(Enzymes[i].site_complet,Enzymes[i+1].site_complet)==ARE_IDENTIC)
						Enzymes[i].select[VECTOR]=FALSE;
					}
				if(i-1>=0) /* or is 'i' is the same  enzyme than next i+1 */
					{
					if(strcmp(Enzymes[i].site_complet,Enzymes[i-1].site_complet)==ARE_IDENTIC)
						Enzymes[i-1].select[VECTOR]=FALSE;
					}
				}
				
			var_find=TRUE;
			while(var_find!=0)
				{
				for(i=0;i<nbr_enzyme;i++)
					if(Enzymes[i].select[VECTOR]==TRUE)
						printf("\tCode: %d\t%s %s\n",i+1,Enzymes[i].nom,Enzymes[i].site_complet);
				printf("\nInput the code number (type '0' to quit):/* 0 */");
				var_find=Get_Number(0);
				if(var_find>0 && var_find<=nbr_enzyme)
					{
					DRAW_HR(stdout,'.');
					printf("\nInformations about %s [%s].\n",Enzymes[var_find-1].nom,Enzymes[var_find-1].site_complet);
					/**Create a link to the rebase www site***/
					printf("Internet WWW link: <A HREF=\"http://www.neb.com/rebase/enz/");
					i=0;
					while(Enzymes[var_find-1].nom[i]!='\0')
						{
						if(Enzymes[var_find-1].nom[i]!=' ')
							printf("%c",Enzymes[var_find-1].nom[i]);
						i++;
						}
					printf(".html\">%s</A>\n",Enzymes[var_find-1].nom);
					DRAW_HR(stdout,'.');
					Cmd_Get_Info2(FORMAT_SCREEN,stdout,var_find-1);
					}
				}
			}
		else
			{
			if(var_search_error==FALSE)
				printf("No Enzyme was found with this query :[%s]\n",Word);
			Menu(12);
			}
		/***************************************************************/
		}/* end else */
	}



/************************************************************************************/
boolean Are_Same_Asymetric(int NumSite_3,int NumSite_5)
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

if(Preference.display_messages==TRUE)
	{
	Menu(13);
	DRAW_LINE;
	}
/* scan all enzymes */
for(i=0;i<nbr_enzyme;i++)
	{
	vara=(double)(Enzymes[i].taille_site);
	if(Preference.memory==TRUE)
		{
		if(vara<=SMALL_SITE)
			{
			if(Preference.display_messages==TRUE)
				printf("Discard %s because it is too short [%s].\n",Enzymes[i].nom,Enzymes[i].site_complet);
			}
		else
			{
			/* take in charge in _num the degenerate bases in the site */
			for(j=0;j<Enzymes[i].taille_site;j++)
				if( (_num=(double)Fct_est_ADN(Enzymes[i].site[j] )) >= AMBIGOUS)
					vara=vara-(_num/4.0);
			/* discard if vara too short */
			if(vara<=SMALL_SITE && Preference.display_messages==TRUE)
				printf("Discard %-8.8s because it is too short and degenerate [%s].\n",Enzymes[i].nom,Enzymes[i].site_complet);
			}
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
if(Preference.display_messages==TRUE)
	{
	DRAW_LINE;
	printf("\t%d enzymes were found in \"%s\".\n", true_nbr_enzyme, FICHIER_ENZYME);
	DRAW_LINE;
	}
}
/************************************************************************************/
boolean freadRebase(FILE *Buffer, char *mot1, char *mot2, char *mot3, char *mot4, char *mot5, char mot6[4])
	{
	/*int Enzyme_Available=0;*/
	char line[MAX_LINE];
	int var_find=0;
	int i,nbr_slash=0;
	
	*mot4=37;
	*mot5='?';
	memset(&mot6[0],0,4); /* 0 means it is not defined */
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
		else if (line[1]==REBASE_TEMPERATURE) /** find the temperature **/
			{
			*mot4=37;
			if(strlen(line)>3 && strlen(line)<6) /* if T¡C is between <T>0 and <T>99 */
				{
				if(isdigit(line[3])!=FALSE) 
					{
					*mot4=0+line[3]-'0';
					if(isdigit(line[4])!=FALSE)
						{
						*mot4=*mot4*10+line[4]-'0';
						}
					}
				else
					*mot4=37;
				}
			}
		else if (line[1]==REBASE_BUFFER) /** find the buffer **/
			{
			*mot5='?';
			if(strlen(line)>3 && strlen(line)<5) /* just coded by ONE char */
				{
				*mot5=line[3];
				}
			}
		else if (line[1]==REBASE_BUFFERS) /** find the buffer **/
			{
			memset(&mot6[0],0,NBR_BUFFER);
			nbr_slash=0;
			for(i=3;i<strlen(line);i++) if(line[i]=='/') nbr_slash++;
			if(nbr_slash==3)
				{
				nbr_slash=0;
				for(i=3;i<strlen(line);i++)
					{
					if(line[i]=='/')
						{nbr_slash++;continue;}
					if(isdigit(line[i])==FALSE)
						{memset(&mot6[0],0,NBR_BUFFER);break;}
					mot6[nbr_slash]=mot6[nbr_slash]*10+line[i]-'0';
					if(mot6[nbr_slash]>100)
						{memset(&mot6[0],0,NBR_BUFFER);break;}
					}
				}
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
int Make_New_enzyme(char *EnzymeName, char *EnzymeSite,int Enzyme_Available, char temp, char buffer, char NEBbuffer[NBR_BUFFER])
	{
	int i;
	int j,vara,varb;
	int var_parenthese_droite,var_parenthese_gauche,var_slash,nbr_parenthese;
	int pos5_3,pos3_5,taille_site;
	char site[PALINDROME_MAX_ENZYME];
	boolean var_non_palindromic=0;
	/* realloc memory for the new enzyme */
	Enzymes=(STRUCT_ENZYME*)realloc(Enzymes,(nbr_enzyme+1)*sizeof(STRUCT_ENZYME));
	if(Enzymes==NULL)
		{
		BEEP;
		printf("**** OUT OF MEMORY !  (too much enzymes ! %d )****\n",nbr_enzyme);
		INKEY;
		Menu(19);
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
				printf("Alerte: l %s [%s]do not contains '^'\n",EnzymeName,EnzymeSite);
				RETURN_NULL;
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
						{if(Preference.display_messages==TRUE) printf("\n %s can't handle %s  (%s) .\n",VAR_VERSION,EnzymeName,EnzymeSite);return(FALSE);}	
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
				Enzymes[nbr_enzyme].temperature=temp;
				Enzymes[nbr_enzyme].NEBBuffer=buffer;
				memcpy(Enzymes[nbr_enzyme].NBuffer,NEBbuffer,NBR_BUFFER*sizeof(char));
				
				
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
				{if(Preference.display_messages==TRUE) printf("\n Alert: %s can't handle %s (%s).\n",VAR_VERSION,EnzymeName,EnzymeSite);return(FALSE);}
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
			Enzymes[nbr_enzyme].temperature=temp;
			Enzymes[nbr_enzyme].NEBBuffer=buffer;
			memcpy(Enzymes[nbr_enzyme].NBuffer,NEBbuffer,NBR_BUFFER*sizeof(char));

			/* I decided to use this field in order to use Enzyme_Available in Discard_Small procedure */		

			nbr_enzyme++;
			true_nbr_enzyme++;
			/* add an enzyme with the same name BUT with the anti-parallele site */
			Enzymes=(STRUCT_ENZYME*)realloc(Enzymes,(nbr_enzyme+1)*sizeof(STRUCT_ENZYME));
			if(Enzymes==NULL)
				{
				BEEP;
				printf("**** OUT OF MEMORY !  (too much enzymes ! %d )****\n",nbr_enzyme);
				INKEY;
				Menu(19);
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
			Enzymes[nbr_enzyme].temperature=temp;
			Enzymes[nbr_enzyme].NEBBuffer=buffer;
			memcpy(Enzymes[nbr_enzyme].NBuffer,NEBbuffer,NBR_BUFFER*sizeof(char));


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
char	mot4,mot5,mot6[NBR_BUFFER];
/*int Enzyme_Available=0;*/
int Nbr_Asym=0;


true_nbr_enzyme=0;
nbr_enzyme=0;


	if ((Buffer = fopen(FICHIER_ENZYME, "r" ))==NULL)
		{
		printf("\n\t----> File error  : '%s'.\n",FICHIER_ENZYME);
		Handle_error(TRUE);
		INKEY;
		}
	else
		{
		while(freadRebase(Buffer,mot1,mot2,mot3,&mot4,&mot5,&mot6[0])==TRUE)
			{
			if(strlen(mot3)>0)
				{
				Make_New_enzyme(mot1, mot2,(int)strlen(mot3),mot4,mot5,&mot6[0]);
				}
			}/*end while*/
		fclose(Buffer);
		if(Preference.display_messages==TRUE)
			printf("\n\t%d enzymes were found in \"%s\".\n", true_nbr_enzyme, FICHIER_ENZYME);				

		/* now, discard isoschizomeres of enzyme i, note that procedure
			Cmd_Get_Rebase stored the number of companies selling the
			enzymes in enzyme field Enzymes[i].select[VECTOR] */
		Discard_iso();
			
		if(Preference.memory==TRUE) Discard_Small();
			
		if(Preference.display_messages==TRUE)
			 {
			 printf("\n\t%d enzymes were found in \"%s\".\n", true_nbr_enzyme, FICHIER_ENZYME);				
			 Menu(12);
			 }
			
		}/* end if file */
Search_Done=FALSE;
}/* end void */



int Discard_iso(void)
	{
	int i,j,vara;
	
	if(Preference.display_messages==TRUE)
		Menu(14);
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
					if(Preference.display_messages==TRUE)
						printf("%s and %s are identic [%s]: Keep %s and discard the other.\n",
							Enzymes[i].nom,Enzymes[j].nom,Enzymes[j].site_complet,Enzymes[j].nom);
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



/* File 'Convert.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include "ADN.h"
#endif

void Rebase2Strider(void)
	{
	/*************************************************************** 
			This function convert a  Rebase Restriction 
			Enzyme list to a DNA Strider database
	REL File example:
		;
		; Restriction Endonuclease file for DNA Striderª.	
		Aat II,	gacgt/c,
		Acc I,	gt/mkac,	
		Bbe I,	ggcgc/c,
		#Bbs I,	gaagac,2,6,        
		#Bbv I,	gcagc,8,12,
		#BceF I,	acggc,12,13,
		Bcl I,	t/gatca,
		Bcn I,	ccs/gg,
		Bgl I,	gccnnnn/nggc,	
		%%%%%%
		You may store comments here.		
	
	****************************************************************/
	/*_______________________________________________*/
	FILE *Buffer_out;
	char Fileout[MAX_NOM_FICHIER];
	char c;
	int i=0,j=0,NumRebase=0;
	/*_______________________________________________*/
	DRAW_LINE;
	Display("Convert The memorized Rebase File to a DNA Striderª REL File:",2);
	printf("Convert the REBASE  File \"%s\".\n",FICHIER_ENZYME);
	if(Preference.memory==TRUE)
		{
		BEEP;
		printf("BEWARE: Memory has been optimized (short sites have been discarded)\n");
		DRAW_LINE;
		}
	printf("BEWARE: Isoschizomers and non-commercialy available enzymes have been discarded\n");
	DRAW_LINE;
	DIRECTORY;
	
	printf("\tInput the name of the DNA Striderª File to save       :/* RELibrary2 */");
	#ifdef __MAC__
		Buffer_out=MacWriteTextFile(Fileout,"\pRELibrary2");
	#else
		Cmd_lire(Fileout);
		if(strcmp(Fileout,"")==ARE_IDENTIC)
			strcpy(Fileout,"RELibrary2");
		Buffer_out = Fct_Write_File(Fileout);
	#endif
	if (Buffer_out==NULL) 
		{
		printf("\t\t---> File error : %s\n",Fileout);
		Handle_error(TRUE);
		BEEP;
		}
	else
		{
		fprintf(Buffer_out,";%s 1997 Pierre LINDENBAUM\n",VAR_VERSION);
		if(Preference.memory==TRUE)
			fprintf(Buffer_out,";BEWARE: Memory has been optimized (short sites have been discarded)\n");
		fprintf(Buffer_out,";BEWARE: Isoschizomers and non-commercialy available enzymes have been discarded\n");
		fprintf(Buffer_out,";Convert the REBASE  File :%s\n",FICHIER_ENZYME);
		fprintf(Buffer_out,";     to the Strider File :%s\n;\n",Fileout);
		for(i=0;i<nbr_enzyme;i++)
			{
			if(i>0)
				if(strcmp(Enzymes[i].nom,Enzymes[i-1].nom)==ARE_IDENTIC)
					continue;
			if(Enzymes[i].palindromic==TRUE)
				{
				j=0;
				fprintf(Buffer_out,"%s,",Enzymes[i].nom);
				while((c=Enzymes[i].site_complet[j++])!='\0') fprintf(Buffer_out,"%c",(c=='^'?'/':LOWER(c)));
				fprintf(Buffer_out,",\n");
				}
			else
				{
				fprintf(Buffer_out,"#%s,",Enzymes[i].nom);
				j=0;
				while((c=Enzymes[i].site_complet[j++])!='(') fprintf(Buffer_out,"%c",LOWER(c));
				fputs(",",Buffer_out);
				while((c=Enzymes[i].site_complet[j++])!='/') fprintf(Buffer_out,"%c",c);
				fputs(",",Buffer_out);
				while((c=Enzymes[i].site_complet[j++])!=')') fprintf(Buffer_out,"%c",c);
				fputs(",\n",Buffer_out);
				}
			NumRebase++;
			}
		fprintf(Buffer_out,"%%%%%%\nYou may store comments here.\n");
		fprintf(Buffer_out,"%d enzyme%c converted.\n",NumRebase,(NumRebase>1?'s':'\0'));
		if(NumRebase>254)
			{
			fprintf(Buffer_out,"BEWARE Remember that DNA Striderª cannot use more than 255 enzymes !\n");
			BEEP;
			Display("BEWARE !!! Too Many Enzymes !",2);
			printf("Remember that DNA Striderª cannot use more than 255 enzymes !\n");
			}
		fclose(Buffer_out);
		DRAW_LINE;
		printf("\t\t--->%d enzyme%c converted.\n",NumRebase,(NumRebase>1?'s':'\0'));
		DRAW_LINE;	
		}
	Menu(12);
	CLS;
	}

void Strider2Rebase(void)
{
/*************************************************************** 
		This function convert a DNA Striderª Restriction 
		Enzyme list (such as the REL library) to a Rebase Format

	REL File example:
		;
		; Restriction Endonuclease file for DNA Striderª.	
		Aat II,	gacgt/c,
		Acc I,	gt/mkac,	
		Bbe I,	ggcgc/c,
		#Bbs I,	gaagac,2,6,        
		#Bbv I,	gcagc,8,12,
		#BceF I,	acggc,12,13,
		Bcl I,	t/gatca,
		Bcn I,	ccs/gg,
		Bgl I,	gccnnnn/nggc,	
		%%%%%%
		You may store comments here.
		
****************************************************************/
FILE *tampon10,*tampon11;
char Filein[MAX_NOM_FICHIER],Fileout[MAX_NOM_FICHIER];
char c,word[255];
int i=0,vara,EnzStrider;
/*______________________________*/
DRAW_LINE;
printf("%s offers the possibilty to convert a DNA Striderª Restriction Enzyme list (such as the REL library) to a Rebase Format.\n",VAR_VERSION);
printf("\nThe RELibrary must be present in the same folder than %s.\n",VAR_VERSION);
DRAW_LINE;

Display("Convert a DNA Striderª REL File to a REBASE File:",2);

DIRECTORY;
printf("\tInput the name of the DNA Strider File (REL) to open /* RELibrary */:");
Cmd_lire(Filein);
if(strcmp(Filein,"")==ARE_IDENTIC)
	strcpy(Filein,"RELibrary");
#ifndef __MAC__
	printf("\tInput the name of the Rebase File to save       :/* RebaseREL */");
	Cmd_lire(Fileout);
	if(strcmp(Fileout,"")==ARE_IDENTIC)
		strcpy(Fileout,"RebaseREL");
#endif
printf( "\n\tOpening file '%s'É\n",Filein);
tampon10 = fopen(Filein, "r" );
if (tampon10==NULL) 
	{
	printf("\t\t---> File error : %s\n",Filein);
	Handle_error(TRUE);
	BEEP;
	}
else
	{
	printf( "\n\tSaving file '%s'É\n",Fileout);
	#ifdef __MAC__
		printf("\tInput the name of the Rebase File to save\n");
		tampon11 = MacWriteTextFile(Fileout,"\pRebaseREL");
	#else
		tampon11 = Fct_Write_File(Fileout);
	#endif
	if (tampon11==NULL) 
		{
		printf("\t\t---> File error : %s\n",Fileout);
		Handle_error(TRUE);
		BEEP;
		}
	else
		{
		EnzStrider=0; /* no enzyme has been found */
		fprintf(tampon11,"Convert the Strider File :%s\n",Filein);
		fprintf(tampon11,"     to the REBASE  File :%s\n\n",Fileout);

		while( (c=fgetc(tampon10)) != EOF )
			{
			if(c==';') /* if ';' -> it is a commentary , go to next line */
				while(fgetc(tampon10)!='\n');
			else
				{
				i=1; /* it is not just a carriage return */
				vara=0;
				if (c!='\n')
					word[0]=c; /* get the first letter of the line */
				else
					i=0; /* it IS a carriage return */
				c=fgetc(tampon10);
				while( c!= '\n' && c!=EOF) /*read to the end of line and put it in word[]*/
					{
					word[i]=c;
					if (c==',') vara++; /* one ',' is palindromic two ',' is asymetric */
					word[i+1]='\0';
					if(i<254) i++;
					c=fgetc(tampon10);
					}
				/** If it is a palindromic enzyme ***/
				if(vara==2)
					{
					EnzStrider++; /* add an enzyme found */
					if(Preference.display_messages==TRUE)
						printf("\t[%d] %s\n",EnzStrider,word);
					fputs("\n<1>",tampon11);
					i=-1;
					while(word[++i]!= ',')
						fprintf(tampon11,"%c",word[i]); /* get the enzyme name and put it in file */
					fputs("\n<5>",tampon11);
					while(word[++i]!= ',')
						{
						if(word[i]=='/')
							fputs("^",tampon11); /* get the enzyme localisation of cleavage and put it in file */
						else if(word[i]!='\t')
							fprintf(tampon11,"%c",word[i]); /* get the enzyme site of cleavage and put it in file */
						}
					fputs("\n<7>?\n",tampon11); /* Enzyme commercialy available by default */
					}
				/** If it is NOT a palindromic enzyme ***/
				if(vara==4)
					{
					EnzStrider++; /* add an enzyme found */
					if(Preference.display_messages==TRUE)
						printf("\t[%d] %s\n",EnzStrider,word);
					fputs("\n<1>",tampon11);
					i=0;
					while(word[++i]!= ',')
						fprintf(tampon11,"%c",word[i]); /* get the enzyme localisation of cleavage and put it in file */
					fputs("\n<5>",tampon11);
					while(word[++i]!= ',') /* get the enzyme site of cleavage and put it in file */
						if(word[i]!='\t') fprintf(tampon11,"%c",word[i]);
					fputs("(",tampon11);
					while(word[++i]!= ',')
						fprintf(tampon11,"%c",word[i]);
					fputs("/",tampon11);
					while(word[++i]!= ',')
						fprintf(tampon11,"%c",word[i]); /* get the enzyme localisation of cleavage and put it in file */
					fputs(")",tampon11);	
					fputs("\n<7>?\n",tampon11); /* Enzyme commercialy available by default */
					}
				}
			}
		fclose(tampon11);
		DRAW_LINE;
		printf("\t\t--->%d enzyme%c found and converted.\n",EnzStrider,(EnzStrider>1?'s':'\0'));
		DRAW_LINE;		
		}
	fclose(tampon10);
	}
Menu(12);
INKEY;
CLS;
}			



/* File 'Cmd_Line.c' */
#ifndef __GNUC__
	#include<stdio.h>
	#include<stdlib.h>
	#include<string.h>
	#include"ADN.h"
#endif

char *search_arg(char *word,int argc, char *argv[])
	{
	static char *result=NULL;
	int i,j,n=0;
	
	result=NULL;
	
	for(i=0;i<argc;i++)
		if((int)strncmp(word,argv[i],strlen(word))==ARE_IDENTIC)
			{
			if(n++>0) {BEEP;printf("Syntax error in parameter \"%s\" !\n",word);Menu(19);}
			result=(char*)realloc(result,(strlen(argv[i])+1)*sizeof(char));
			if(result==NULL) {BEEP;printf("Memory error (search_arg)\n");Menu(19);}
			if(strlen(word)<=strlen(argv[i]))
				for(j=(int)strlen(word);j<=(int)strlen(argv[i]);j++)
					result[(int)(j-(int)strlen(word))]=argv[i][j];
			}
	return(result);
	}

int val(char *word, int var_defaut)	
	{
	int i,var_return=0;
	if(word==NULL) return(var_defaut);
	if(strspn( word,"0123456789") != strlen(word))
		{BEEP;printf("Text-Error in parameter \"%s\" !\n",word);Menu(19);}
	for(i=0;i<(int)strlen(word);i++)
		var_return = var_return*10+word[i]-'0';
	return(var_return);
	}

void Cmd_Line(int argc, char *argv[])
	{
	char	MyChoice[MAX_NOM_FICHIER];
	boolean arg_cip=TRUE,arg_t4=TRUE,arg_part=TRUE;

	
	Preference.mini_display=TRUE;
		
	if(search_arg("-HELP",argc,&argv[0])!=NULL)
		{
		Menu(18);
		#ifndef __GNUC__
			INKEY;
		#endif
		exit(0);
		}
	
	if(search_arg("-U",argc,&argv[0])!=NULL)
		{
		strcpy(MyChoice,search_arg("-U",argc,&argv[0]));
		Cloning_Project2(MyChoice);
		Menu(19);
		}
	
	  /**************/
	 /* init prefs */
	/**************/
	
	Preference.side_5				=(search_arg("-H",argc,&argv[0])==NULL?Preference.side_5:TRUE);
	Preference.side_3				=(search_arg("-O",argc,&argv[0])==NULL?Preference.side_3:TRUE);
	Preference.memory				=(search_arg("-Y",argc,&argv[0])==NULL?TRUE:FALSE);
	Preference.allow_part_overhang	=(search_arg("-E",argc,&argv[0])==NULL?FALSE:TRUE);
	Preference.temperature			=(search_arg("-TP",argc,&argv[0])==NULL?FALSE:TRUE);
	Preference.buffer				=(search_arg("-BU",argc,&argv[0])==NULL?FALSE:TRUE);
	
	
	arg_part	=(search_arg("-T",argc,&argv[0])==NULL?TRUE:FALSE);
	arg_cip		=(search_arg("-Z",argc,&argv[0])==NULL?TRUE:FALSE);
	arg_t4		=(search_arg("-P",argc,&argv[0])==NULL?TRUE:FALSE);

	
	Preference.DeltaMax				=val(search_arg("-M",argc,&argv[0]),Preference.DeltaMax);
	Preference.DeltaMin				=val(search_arg("-m",argc,&argv[0]),Preference.DeltaMin);
	if(Preference.DeltaMax<=Preference.DeltaMin+1) {BEEP;printf("Cloning box error.\n");Menu(19);}
	
	Preference.search_C_term		=TRUE;
	Preference.display_messages		=FALSE;

	
	  /***************/
	 /* open vector */
	/***************/
	if(search_arg("-VE",argc,&argv[0])!=NULL)
		{
		strcpy(FICHIER_ADN[VECTOR],search_arg("-VE",argc,&argv[0]));
		Cmd_LOAD_SEQUENCE(VECTOR);
		if(npb[VECTOR]==0) Menu(19);
		var_min[VECTOR]	=	val(search_arg("-VL",argc,&argv[0]),var_min[VECTOR]);
		var_max[VECTOR]	=	val(search_arg("-VR",argc,&argv[0]),var_max[VECTOR]);
		pos_ATG[VECTOR] =	val(search_arg("-VA",argc,&argv[0]),pos_ATG[VECTOR]);
		
		if(var_min[VECTOR]>=var_max[VECTOR] || var_max[VECTOR]>npb[VECTOR] || pos_ATG[VECTOR]>npb[VECTOR])
			{BEEP;printf("Error in vector parameters\n");Menu(19);}
		}
	
	  /***************/
	 /* open INSERT */
	/***************/
	if(search_arg("-IN",argc,&argv[0])!=NULL)
		{
		strcpy(FICHIER_ADN[INSERT],search_arg("-IN",argc,&argv[0]));
		Cmd_LOAD_SEQUENCE(INSERT);
		var_min[INSERT]	=	val(search_arg("-IL",argc,&argv[0]),var_min[INSERT]);
		var_max[INSERT]	=	val(search_arg("-IR",argc,&argv[0]),var_max[INSERT]);
		sigma_5			=   val(search_arg("-ILi",argc,&argv[0]),var_min[INSERT]+sigma_5)-var_min[INSERT];
		sigma_3			=   var_max[INSERT]-val(search_arg("-IRi",argc,&argv[0]),var_max[INSERT]-sigma_3);
		pos_ATG[INSERT] =	val(search_arg("-IA",argc,&argv[0]),pos_ATG[INSERT]);
		if(var_min[INSERT]>=var_max[INSERT] || var_max[INSERT]>npb[INSERT] || pos_ATG[INSERT]>npb[INSERT])
			{BEEP;printf("Error in insert parameters\n");Menu(19);}

		}
	
	  /***************/
	 /* open rebase */
	/***************/
	if(search_arg("-R",argc,&argv[0])!=NULL)
		{
		strcpy(FICHIER_ENZYME,search_arg("-R",argc,&argv[0]));
		printf("\t\tOpening \"%s\"\n\n",FICHIER_ENZYME);
		Cmd_Get_Rebase();
		}
	
	
	#ifdef VAR_UNIX
	/* GCG interface */
	if(getenv(GCG_TAG)!=NULL)
		{
		/*if(search_arg("-BN",argc,&argv[0])!=NULL) GCGBestFit(TRUE);*/
		/*if(search_arg("-BP",argc,&argv[0])!=NULL) GCGBestFit(FALSE);*/
		if(search_arg("-GN",argc,&argv[0])!=NULL) GCGGap(TRUE);
		if(search_arg("-GP",argc,&argv[0])!=NULL) GCGGap(FALSE);
		}
	#endif
	  /*************************/
	 /* find cloning strategy */
	/*************************/
	if(search_arg("-A",argc,&argv[0])!=NULL)
		{
		if( npb[VECTOR]==0 || pos_ATG[VECTOR]==FALSE || var_min[VECTOR]==FALSE ||
			var_max[VECTOR]==FALSE || var_max[VECTOR]==FALSE || npb[INSERT]==0 || pos_ATG[INSERT]==FALSE ||
			var_min[INSERT]==FALSE || var_max[INSERT]==FALSE ||
			var_min[INSERT]>=var_min[INSERT]+sigma_5+1 || var_max[VECTOR]<=var_min[VECTOR] || var_max[INSERT]<=var_min[INSERT] || var_max[INSERT]<=var_max[INSERT]-sigma_3)
			{Display("Sequence Error.",1);BEEP;Menu(19);}
		if(nbr_enzyme<=0)
			{Display("No enzyme was found !!!.",1);BEEP;Menu(19);}
		Cmd_Get_Site();
		if(nbr_sites>0)
			{
			Init_Enz();
			if(META_CloneIt4(arg_cip,arg_t4,arg_part)!=TRUE) printf("No solution was found.\n");
			Menu(19);
			}
		else
			{Display("No site was found !!!.",1);BEEP;Menu(19);}
		}
		
	
	  /***************************/
	 /* find in-frame deletions */
	/***************************/
	else if(search_arg("-D",argc,&argv[0])!=NULL)
		{
		if(var_max[INSERT]<=var_min[INSERT])
			{Display("Box Error.",1);BEEP;Menu(19);}
		if(npb[INSERT]<=0)
			{Display("Sequence Error.",1);BEEP;Menu(19);}
		if(nbr_enzyme<=0)
			{Display("No enzyme was found !!!.",1);BEEP;Menu(19);}
		Cmd_Get_Site();
		if(nbr_sites>0)
			{
			/*Init_Enz();*/
			if(META_DeltaFrame(arg_t4,arg_part)!=TRUE) printf("No solution was found.\n");
			Menu(19);
			}
		else
			{Display("No site was found !!!.",1);BEEP;Menu(19);}
		}
	  /*******************/
	 /* find frameshift */
	/*******************/
	else if(search_arg("-F",argc,&argv[0])!=NULL)
		{
		if(var_max[INSERT]<=var_min[INSERT])
			{Display("Box Error.",1);BEEP;Menu(19);}
		if(npb[INSERT]<=0)
			{Display("Sequence Error.",1);BEEP;Menu(19);}
		if(nbr_enzyme<=0)
			{Display("No enzyme was found !!!.",1);BEEP;Menu(19);}
		Cmd_Get_Site();
		if(nbr_sites>0)
			{
			/*Init_Enz();*/
			if(META_FrameShift(arg_part)!=TRUE) printf("No solution was found.\n");
			Menu(19);
			}
		else
			{Display("No site was found !!!.",1);BEEP;Menu(19);}
		}
	  /*******************/
	 /* restriction map */
	/*******************/
	else
		{
		if(search_arg("-N",argc,&argv[0])!=NULL && (npb[VECTOR]>0 || npb[INSERT]>0))
			{
			Cmd_Get_Site();
			Intersections(stdout,FORMAT_SCREEN);
			Menu(19);
			}
		if(search_arg("-VM",argc,&argv[0])!=NULL && npb[VECTOR]>0)
			{
			Cmd_Get_Site();
			ResMap(stdout,FORMAT_SCREEN,VECTOR);
			Menu(19);
			}
		if(search_arg("-IM",argc,&argv[0])!=NULL && npb[INSERT]>0)
			{
			Cmd_Get_Site();
			ResMap(stdout,FORMAT_SCREEN,INSERT);
			Menu(19);
			}
		}
	Preference.mini_display=FALSE;
	}



/* File 'Menu.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "ADN.h"
	#ifdef __MAC__
		#include <console.h>
	#endif
#endif



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
		
		if(mode==FORMAT_HTML) fprintf(out,"<FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));
		
		for(i= site_5.Loc - 4; i<= site_5.Loc +PALINDROME_MAX_ENZYME;i++)
			{
			if (Fct_Frame(N_Seq,i)==IS_IN_FRAME) fprintf(out,".");
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos5_3) {fprintf(out,"/");if(mode==FORMAT_HTML) fprintf(out,"</FONT><FONT COLOR=\"#%s\">",(N_Seq==VECTOR?RED_COLOR:BLUE_COLOR));}
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos3_5) fprintf(out," ");
			fprintf(out,"%c",sequence[N_Seq][Fct_Pos(N_Seq,i)]);
			}
		/* write second part */
		fprintf(out,"--  --");
		for(i= site_3.Loc - 4; i<= site_3.Loc + PALINDROME_MAX_ENZYME;i++)
			{
			if (Fct_Frame(N_Seq,i)==IS_IN_FRAME) fprintf(out,".");
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos5_3) {fprintf(out,"/");if(mode==FORMAT_HTML) fprintf(out,"</FONT><FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));}
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos3_5) fprintf(out," ");
			fprintf(out,"%c",sequence[N_Seq][Fct_Pos(N_Seq,i)]);
			}
		
		if(mode==FORMAT_HTML) fprintf(out,"</FONT>");
		
		/* write anti-first part */
		fprintf(out,"--  3'\n  3'  --");
		
		if(mode==FORMAT_HTML) fprintf(out,"<FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));
		
		for(i= site_5.Loc - 4; i<= site_5.Loc +PALINDROME_MAX_ENZYME;i++)
			{
			if (Fct_Frame(N_Seq,i)==IS_IN_FRAME) fprintf(out,".");
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos3_5) {fprintf(out,"/");if(mode==FORMAT_HTML) fprintf(out,"</FONT><FONT COLOR=\"#%s\">",(N_Seq==VECTOR?RED_COLOR:BLUE_COLOR));}
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos5_3) fprintf(out," ");
			fprintf(out,"%c",Fct_Complementaire(sequence[N_Seq][Fct_Pos(N_Seq,i)]));
			}
		/* write anti-second part */
		fprintf(out,"--  --");
		for(i= site_3.Loc - 4; i<= site_3.Loc + PALINDROME_MAX_ENZYME;i++)
			{
			if (Fct_Frame(N_Seq,i)==IS_IN_FRAME) fprintf(out,".");
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos3_5) {fprintf(out,"/");if(mode==FORMAT_HTML) fprintf(out,"</FONT><FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));}
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos5_3) fprintf(out," ");
			fprintf(out,"%c",Fct_Complementaire(sequence[N_Seq][Fct_Pos(N_Seq,i)]));
			}
			
		if(mode==FORMAT_HTML) fprintf(out,"</FONT>");
			
		fprintf(out,"--  5'\n");
		/* write the translation */
		
		

		
		if(pos_ATG[INSERT]!=FALSE || pos_ATG[VECTOR]!=FALSE)
			{
			/* write first part */
			fprintf(out,"  NH2    ");
			if(mode==FORMAT_HTML) fprintf(out,"<FONT COLOR=\"#%s\">",YELLOW_COLOR);
			for(i= site_5.Loc - 4; i<= site_5.Loc +PALINDROME_MAX_ENZYME;i++)
				{
				if(Fct_Frame(N_Seq,i)==IS_IN_FRAME)
					fprintf(out,"%c ",Fct_Traduction(sequence[N_Seq][Fct_Pos(N_Seq,i)]
												,sequence[N_Seq][Fct_Pos(N_Seq,i+1)]
												,sequence[N_Seq][Fct_Pos(N_Seq,i+2)]));
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
					fprintf(out,"%c ",Fct_Traduction(sequence[N_Seq][Fct_Pos(N_Seq,i)]
											,sequence[N_Seq][Fct_Pos(N_Seq,i+1)]
											,sequence[N_Seq][Fct_Pos(N_Seq,i+2)]));
				else
					fprintf(out," ");
				if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos5_3) fprintf(out," ");
				if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos3_5) fprintf(out," ");
				}
			if(mode==FORMAT_HTML) fprintf(out,"</FONT>");
			fprintf(out,".COOH\n");
			}
		}
	else /* if the two sites are close from each other, write the sequence in one time */
		{
		/* write direct strand */
		fprintf(out,"\n  5'  --");
		
		if(mode==FORMAT_HTML) fprintf(out,"<FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));

		
		for(i= site_5.Loc - 4; i<= site_3.Loc +PALINDROME_MAX_ENZYME;i++)
			{
			if (Fct_Frame(N_Seq,i)==IS_IN_FRAME) fprintf(out,".");
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos5_3) fprintf(out,"/");
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos3_5) fprintf(out," ");
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos5_3) fprintf(out,"/");
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos3_5) fprintf(out," ");
			fprintf(out,"%c",sequence[N_Seq][Fct_Pos(N_Seq,i)]);
			}
		/* write anti-strand */
		if(mode==FORMAT_HTML) fprintf(out,"</FONT>");
		fprintf(out,"--  3'\n  3'  --");
		if(mode==FORMAT_HTML) fprintf(out,"<FONT COLOR=\"#%s\">",(N_Seq==VECTOR?BLUE_COLOR:RED_COLOR));
		for(i= site_5.Loc - 4; i<= site_3.Loc +PALINDROME_MAX_ENZYME;i++)
			{
			if (Fct_Frame(N_Seq,i)==IS_IN_FRAME) fprintf(out,".");
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos5_3) fprintf(out," ");
			if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos3_5) fprintf(out,"/");
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos5_3) fprintf(out," ");
			if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos3_5) fprintf(out,"/");
			fprintf(out,"%c",Fct_Complementaire(sequence[N_Seq][Fct_Pos(N_Seq,i)]));
			}
		if(mode==FORMAT_HTML) fprintf(out,"</FONT>");
		fprintf(out,"--  5'\n");
		
		/* write the translation */
		if(pos_ATG[INSERT]!=FALSE || pos_ATG[VECTOR]!=FALSE)
			{
			fprintf(out,"  NH2    ");
			if(mode==FORMAT_HTML) fprintf(out,"<FONT COLOR=\"#%s\">",YELLOW_COLOR);
			for(i= site_5.Loc - 4; i<= site_3.Loc +PALINDROME_MAX_ENZYME;i++)
				{
				if(Fct_Frame(N_Seq,i)==IS_IN_FRAME)
					fprintf(out,"%c ",Fct_Traduction(sequence[N_Seq][Fct_Pos(N_Seq,i)]
												,sequence[N_Seq][Fct_Pos(N_Seq,i+1)]
												,sequence[N_Seq][Fct_Pos(N_Seq,i+2)]));
				else
					fprintf(out," ");
				if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos5_3) fprintf(out," ");
				if (i==site_5.Loc+Enzymes[site_5.NumEnz].pos3_5) fprintf(out," ");
				if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos5_3) fprintf(out," ");
				if (i==site_3.Loc+Enzymes[site_3.NumEnz].pos3_5) fprintf(out," ");
					
				}
			if(mode==FORMAT_HTML) fprintf(out,"</FONT>");
			fprintf(out,".COOH\n");
			}
		}
	}


/************************************************************************************/
void Display(char *s,int mode)
	{
	/* display the word 's' underligned or in a box */
	int i=0;
	switch(mode)
		{
		case(1):
			printf("%s\n",s);
			while(s[i++]!=0) printf("-");PARAGRAF;
			break;
		case(2):
			for(i=0;i<strlen(s)+4;i++) printf("*");
			printf("\n* %s *\n",s);
			for(i=0;i<strlen(s)+4;i++) printf("*");
			PARAGRAF;
			break;
		case(3):
		default:
		#ifdef __MAC__
			MacAlert(s);
		#else
		printf("***** %s ******\n",s);
		#endif
		break;
		}
	}
/************************************************************************************/
void Menu(int vara)
	{
	/* display various messages */
	switch(vara)
		{
		case(0):
				#ifdef __MAC__
					MacAbout();
				#else
					{
					printf( "\n\n ___________________________________________________\n" );
					printf( "  %s  Pierre LINDENBAUM\n",VAR_VERSION );
					printf( "                 Biologie Moleculaire des Rotavirus.\n" );
					printf( "                 VIM INRA.\n" );
					printf( "                 78350 JOUY-EN-JOSAS.\n" );
					printf( "                 FRANCE.\n" );
					printf( "                 E.Mail:lindenb@biotec.jouy.inra.fr.\n" );
					printf( " ___________________________________________________\n\n" );
					printf( "This program and its source code are protected by a trade mark:\n\t(National Number 97/704582)\n\n" );
					}
				#endif
				break;
		case(1):	
					{
					#ifndef __MAC__
					printf("The sequence in FASTA format or DNA Strider (text ASCII or binary) must be localized in the same folder than %s.",VAR_VERSION);
					#else
					printf("The sequence are in FASTA format or DNA Strider (text ASCII or binary). ");
					#endif
					printf(" Degenerate template are not allowed. Numbers will be discarded. ");
					printf("The DNA sequence is a circular plasmid sequence. Sequence length is limited to %d bp but this can be changed on the level of the source code. ",MAX_NPB);
					printf("A short database of different sequence of classic oligonucleotides is used by the program to try to localize the insert bounds. This database is written as a text ");
					printf("in the \"%s\" file. If %s can't recognize the plasmid, the user must ",POLYLINKER_FILE,VAR_VERSION);
					printf("input the positions of those extremities himself so he will have to know them before running the program. ");
					printf("Translated sequences are supposed to be oriented from NH2 to COOH. Input is case sensitive.\n");
					DRAW_LINE;
					}
					break;
		case(2):
					printf("VECTOR\n");
					printf("\t            cloning box\n");
					printf("\t-----------[===========]------------\n");
					printf("\t           |           |\n");
					printf("\t           |           3' boundary\n");
					printf("\t           5' boundary\n\n");
					printf("INSERT\n");
					printf("\t                5' INTERNAL limit\n");
					printf("\t                | cloning 3' INTERNAL limit\n");
					printf("\t                |   box   |\n");
					printf("\t-----------[====]---------[===]-----\n");
					printf("\t           |                  |\n");
					printf("\t           |                  3' boundary\n");
					printf("\t           5' boundary\n\n");
					DRAW_LINE;
					break;
		case(3):	printf("REBASE File ((c)New England BioLabs) is a database of all known Restriction");
					printf(" Enzymes. It has been created by Dr. Richard J.Roberts and is freely ");
					printf("available on ftp servor: ftp://www.neb.com/pub/rebase/.  Do use an 'allenz' file");
					printf(" to get all known enzymes.A REBASE file is a text file (ASCII)that must be present");
					printf(" in the same folder than %s.Input is case sensitive.",VAR_VERSION);
					printf(" %s use the standard symbolic abreviations. ",VAR_VERSION);
					printf("%s don't handle enzymes cuting on both sides of the site  ex:  Bsp24I (8/13)GACNNNNNNTGG(12/7).\n\n",VAR_VERSION);
					DRAW_LINE;
					break;
		case(4):	printf("Polymerases proprieties.\n");
					printf("\tKlenow Polymerase: 5'->3' pol.\n");
					printf("\t\t..G       -----\\    ..GAATT\n");
					printf("\t\t..CTTAA   -----/    ..CTTAA\n\n");
					printf("\tT4 DNA Polymerase: 5'->3' pol and 3'->5' exo.\n");
					printf("\t\t..CTGCA   -----\\    ..C\n");
					printf("\t\t..C       -----/    ..C\n\n");
					Menu(12);
					break;
		case(5):	DRAW_LINE;
					printf("If two Enzymes with different names cut the same site, one will be discarded.\n");
					printf("Enzymes with short recognicion site (<= %d bases or equivalent) will be discarded.\n",SMALL_SITE);
					Menu(12);
					break;
		case(6):	DRAW_LINE;
					printf("For a couple of Enzyme, if %s find a solution to a cloning problem ",VAR_VERSION);
					printf(" there can be another solutions (Using CIP, modifying Polymerase...) that user could prefer. Nevertheless, the default solution ");
					printf("should be the easier one.\n");
					Menu(12);
					break;
		case(7):	DRAW_LINE;
					printf("If partial overhang ligation is allowed, then such a ligation is possible:\n");
					printf("\t  ..GGGT     CGGA..  --\\  ..GGGTCGGA  --\\  ..GGGTCGGA\n");
					printf("\t  ..CCCAG       T..  --/  ..CCCAG  T  --/  ..CCCAGggT\n");
					Menu(12);
					break;
		case(8):	DRAW_LINE;
					printf("Calf Intestinal Phosphatase is used to dephosphorylate VECTOR. This avoid the self ligation of VECTOR.\n");
					printf("..G      p-AATTC..  --\\  ..G      AATTC..\n");
					printf("..CAATT-p      G..  --/  ..CAATT      G..\n");
					printf("\t---> Self ligation impossible.\n");
					Menu(12);
					break;
		case(9):	DRAW_LINE;
					printf("Ligation will be in frame directed on 5' side of INSERT.\n\n");
					printf(" Vector: 123 123                           123 123 123\n\n");
					printf(" Insert:         123 123 123 123 123 123 1\n\n");
					printf("               _____\n");
					printf("                 In frame ligation (NH2).\n");
					Menu(12);
					break;
		case(10):	DRAW_LINE;
					printf("Ligation will be in frame directed on 3' side of INSERT.\n\n");
					printf(" Vector: 123 1                          23 123 123\n\n");
					printf(" Insert:       123 123 123 123 123 123 1\n\n");
					printf("                                     _____\n");
					printf("                              In frame ligation (COOH).\n");
					Menu(12);
					break;
		case(11):BEEP;
				Display("ABORTED",2);
				break;
		case(12):DRAW_LINE;
				fflush(stdin);
				#ifndef __MAC__
					printf("Press 'Return' to continue.\n");
					INKEY;
				#else
					printf("Press a key to continue.\n");
					while(kbhit()!=1);
				#endif
				break;
		case(13):DRAW_LINE;
				printf("Discard non-useful Enzymes...\n");
				printf("As many short sites have a very low probabilty to be used in a clonage, because \
						they cut a sequence  everywhere, such cuting enzymes will be discarded \
						(Ex: Taq I [t^cga] or NlaIV[ggn^ncc])\n");
				break;
		case(14):DRAW_LINE;
				printf("Discard Isoschyzomers...\n");
				printf("Discard Enzymes with duplicate sites. %s will keep the most often commercialy available.\n",VAR_VERSION);
				DRAW_LINE;
				break;
		case(15):DRAW_LINE;
				printf("You can choose to allow partial digestions only if the used enzyme cus BLUNT. Such sites won't modify the sequence if ");
				printf("a modifying polymerase is used in a cloning strategy.\n");
				DRAW_LINE;
				Menu(12);
				break;
		case(16):
				DRAW_LINE;
				printf("%s can have a look after sites creating Carboxy terminal deletion. Those sites that do not necesseraly create in frame deletion.\n",VAR_VERSION);
				Menu(12);
				break;
		case(17):DRAW_LINE;
				Display("Syntax:",1);
				printf("  Eco\n");printf("  gaattc\n");printf("  g^aattc\n");
				printf("  co\n");printf("  G^A\n");printf("  (1/\n");printf("  (\n");
				printf("  $6        (site length = 6) \n");
				printf("  _4        (all the overhangs of 4 bases)\n");
				printf("  _B        (all the blunt sites)\n");
				printf("  $5_B      (site length = 6 and blunt)\n");
				printf("  _+        (all the 3'overhanged sites)\n");
				printf("  _-        (all the 5'overhanged sites)\n");
				printf("  _-4       (all the 3'overhangs 4 bases)\n");
				printf("  _-4|gaa   (all the 3'overhangs of 4 bases and containing 'gaa')\n");
				printf("  $6_-4|gaa (site length=6, 3'overhangs of 4 bases and containing 'gaa')\n");
				printf("Please, note that if you are searching an enzyme by its NAME, this should be one of those which have been SELECTED from the rebase file, so if you are searching an isoschisomer you'd better ask for its site rather than its name.\n");
				DRAW_LINE;
				PARAGRAF;
				break;
		case(18):
				Display("Command Line Memento.",1);
				printf("-HELP  help (this screen)\n");
				printf("   An help file is available at %s.\n",VAR_URL);
				#ifdef VAR_UNIX
					printf("   This help file is accessible from %s (use menu n¡%d)\n",VAR_VERSION,MENU_LYNX);
				#endif
				PARAGRAF;
				printf("  -U   run a project -UMy_project\n\n");
				
				printf("  -VE  load vector sequence -VEpGBT9\n");
				printf("  -VL  vector cloning box left boundary -VL855\n");
				printf("  -VR  vector cloning box right boundary -VR900\n");
				printf("  -VA  vector ATG position -VA644\n");
				printf("  -VM  vector restriction map\n\n");
				printf("  -IN  load insert sequence -INpbs_myGene\n");
				printf("  -IL  insert cloning box left boundary -IL154\n");
				printf("  -IR  insert cloning box right boundary -IR2001\n");
				printf("  -ILi insert cloning box left internal boundary -ILi160\n");
				printf("  -IRi insert cloning box right internalboundary -IRi1980\n");
				printf("  -IA  insert ATG position -IA155\n");
				printf("  -IM  insert restriction map\n\n");
				
				printf("  -R   open a Rebase file -Rallenz\n");
				printf("  -Y   don't use memory optimization default:FALSE\n");
				printf("  -E   allow non-overlapping overhangs default:FALSE\n");
				printf("  -P   don't allow use of Klenow or T4 DNA Pol default:allow them\n");
				printf("  -T   don't allow partial digestions default:allow it\n\n");
				printf("  -Z   don't allow C.I.P. default:allow C.I.P.\n");
				printf("  -BU  don't allow incompatible buffers. default:allow incompatible buffers\n");
				printf("  -TP  don't allow incompatible T¡C. default:allow incompatible T¡C\n");
				#ifdef VAR_UNIX
				if(getenv(GCG_TAG)!=NULL)
					{
					printf("  -BN  BestFit((c) GCG Wisconsin package) with DNA sequences \n");
					printf("  -BP  BestFit((c) GCG Wisconsin package) with translated sequences \n");
					printf("  -GN  Gap((c) GCG Wisconsin package) with DNA sequences \n");
					printf("  -GP  Gap((c) GCG Wisconsin package) with translated sequences \n");
					}
				else
					printf("  (c)GCG Wisconsin package) with translated sequences \n");
				#endif
				printf("  -N   intersections\n");
				printf("\n Cloning an INSERT into a VECTOR.\n");
				printf("  %s will try to find the best cloning strategy. It will stop when cloning conditions will be found.\n",VAR_VERSION);
				printf("  -A   find Cloning strategies\n");
				printf("  -H   clone in frame at 5' of insert (NH2)  default:FALSE\n");
				printf("  -O   clone in frame at 3' of insert (COOH) default:FALSE\n");
				printf("\n Finding in-frame deletions and frameshifts.\n");
				printf("  %s will try to find a maximum of cloning strategies.\n",VAR_VERSION);
				printf("  -F   find frameshifts\n");
				printf("  -D   find in-frame deletions\n");
				printf("  -M   Maximum percentage of insert -M80\n");
				printf("  -m   minimum percentage of insert -m20\n");
				printf("  -P   don't allow use of Klenow or T4 DNA Pol default:allow them\n");
				printf("  -T   don't allow partial digestions default:allow it\n\n");
				printf("\n Examples:\n");
				printf("  I want to clone the insert of pbs_RF2 into pGBT9 in frame in 5' (NH2)\n  using the Rebase file 'allenz'.\n");
				printf("        ---> arguments are: -VEpGBT9 -INpbs_RF2 -Rallenz -A -H\n");
				printf("  Now, the frame of interest is defined on INSERT at position 430. The\n  cloning box of pGBT9 is defined between 500 and 523. ");
				printf("I don't want to use\n  memory optimization, C.I.P.or partial digestion.\n");
				printf("        ---> arguments are: -VEpGBT9 -INpbs_RF2 -Rallenz -A -H -IA430 -VL500 -VR523 -Y -Z -T\n"); 
				break;
		case(19):
				#ifdef __MAC__
				MacAbout();
				#else
				Menu(0);
				printf("\nIf you use this program, please, cite it in your papers:\n%s\n",VAR_CITATION);
				printf("\n%s WWW Home Page:\n%s\n\n",VAR_VERSION,VAR_URL);
				#ifndef __GNUC__
					INKEY;
				#endif
					exit(0);
				#endif
				break;
		#ifdef VAR_UNIX
		case(20):
				printf(" (c)Genetics Computer Group, Inc. a wholly owned subsidiary of Oxford Molecular Group. All rights reserved.\n");
				DRAW_LINE;
				break;
		#endif
		case(21):printf("Important note: In order to manage incompatible temperatures and buffers, I decided to had 3 new non-official fields in the REBASE database.\n");
		 		 printf("You will have to modify your REBASE by yourself..\n");
		 		 printf("\n\t<%c>65\t(temperature , default is 37)\n",REBASE_TEMPERATURE);
		 		 printf("\t<%c>2\t(NEB specific buffer, one character only)\n",REBASE_BUFFER);
		 		 printf("\t<%c>0/25/100/75\t(percentage of activity in four classical buffers)\n",REBASE_BUFFERS);
		 		 printf("Thoose new fields should be localized between tag <1> and <8>..\n");
		 		 DRAW_LINE;
				 INKEY;
				 break;
		case(22):printf("Copyright Notice: Permission to use, copy, modify, and distribute this software and\n");
				printf("its documentation is hereby granted, subject to the following restrictions\n");
				printf("and understandings:\n");
				printf("1) Any copy of this software or any copy of software derived\n");
				printf("from it must include this copyright notice in full.\n");
				printf("2) All materials or software developed as a consequence of the\n");
				printf("use of this software or software derived requires the express,\n");
				printf("written author permission.\n");
				printf("3) The software may be used by anyone for any purpose, except\n");
				printf("that its redistribution for profit or inclusion in other\n");
				printf("software sold for profit requires the express, written\n");
				printf("permission of the author.\n");
				printf("4) This software is provided AS IS with no warranties of any\n");
				printf("kind. In no event will the author be liable for any lost revenue \n");
				printf("or profits or other special, indirect and consequential damages. \n");
				DRAW_LINE;
				Menu(12);
				break;
		default:PARAGRAF;break;
		}
	}





/* File 'GetSite.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "ADN.h"
#endif
/*** Add a Site in memory see : Cmd_Get_Site  **/
int Add_Site(int k,int _NumSeq,int _Loc)
	{
	
	Sites=(STRUCT_SITE*)realloc(Sites,(nbr_sites+1)*sizeof(STRUCT_SITE));
	if(Sites==NULL)
		{
		BEEP;
		printf("**** OUT OF MEMORY (TOO MUCH SITES %d) !  ****\nPress return to continue.\n",nbr_sites);
		BEEP;
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
#ifndef __MAC__
	DRAW_LINE;
	printf("\t\t\tLooking for sites :");
	if(Preference.display_messages!=TRUE) PARAGRAF;
#else
	ShowProgressBar(1,FALSE,FALSE,"");
	ShowProgressBar(2,FALSE,FALSE,"Looking for sites...");
#endif
nbr_sites=0;
/* scan VECTOR and INSERT */
for(NumSeq=VECTOR;NumSeq<=INSERT;NumSeq++)
	{
	if(npb[NumSeq]>0)
	/* scan all enzymes */
	for(k=0;k<nbr_enzyme;k++)
		{
		#ifndef __MAC__
			if(Preference.display_messages==TRUE)
				{
				printf("%03d %%",abs((int)(((double)NumSeq)*50.0+50.0*(double)k/(double)(nbr_enzyme-1))));
				fflush(stdout);
				printf("%c%c%c%c%c",8,8,8,8,8);
				}
		#else
			ShowProgressBar(3,(int)((NumSeq+1)*(float)k),(int)(2*nbr_enzyme-1),"");
		#endif
		/* if the enzyme site does not contain degenerate base, strpbrk function can be used */
		if(strpbrk(Enzymes[k].site,"YRMKSWBDHVN")==NULL) /* Does not contain degenerate bases */
			{
			p_seq=&sequence[NumSeq][1];
			while((result=strstr(p_seq, Enzymes[k].site ))!=NULL)
				{
				/* if the site has been found, add a new site */
				i=npb[NumSeq]-strlen(result)+1;
				if(Add_Site(k,NumSeq,npb[NumSeq]-strlen(result)+1)==FALSE)
						{NumSeq=INSERT+1;nbr_enzyme=0;break;} /* abort if problem */
				p_seq=result+1;
				}
			start=npb[NumSeq]-Enzymes[k].taille_site+2;
			}
		else
			{
			start=1; /* else scan from the begin of the sequence */
			}

		Nbr_N5 = Enzymes[k].Nbr_N5; /* there is no use to scan the 'N' at the begining of the site see: Get_Rebase */
		for(i=start;i<=npb[NumSeq];i++) /*** Scanning sequence *******/ 
			{
			/* site not found by default */
			vara=FALSE;
			for(j=0;j<Enzymes[k].taille_site;j++)
				if(Fct_Identique(sequence[NumSeq][Fct_Pos(NumSeq,i+j)],Enzymes[k].site[j])!=1)
					{vara=TRUE;break;}
			/* if site was found, add a site in memory */
			if (vara==FALSE)
				{
				if(Add_Site(k,NumSeq,i-Nbr_N5)==FALSE)
					{NumSeq=INSERT+1;i=npb[NumSeq]+1;k=nbr_enzyme;NumSeq=3;nbr_enzyme=0;break;}/* abort if problem */
				}
			}
		}
	}
#ifndef __MAC__
	if(Preference.display_messages==TRUE)
		{
		printf("%03d %%",100);
		fflush(stdout);
		}
	PARAGRAF;
	DRAW_LINE;
	printf("     %d sites found.\n",nbr_sites);
	DRAW_LINE;
#else
	ShowProgressBar(4,FALSE,FALSE,"");
#endif
Search_Done=TRUE;
}



/* File 'FullResMap.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	/*#include <assert.h>*/
	#include "ADN.h"
#endif
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
	register int i,n;
	int var_begin=0,var_end=0; /* begin and end of the sequence description */
	div_t	r;
	/**********************************************************************************/
	var_len_screen=LARGEUR_ECRAN-NOM_MAX_ENZYME; /* length of a on the screen */
	line[0]=EOS;linelen=0;
	if(npb[var_seq]==0) return(FALSE);
	
	if(mode==FORMAT_SCREEN)
		{
		CLS;
		Display("Restriction map.",2);
		printf("Unique sites will be displayed with an UPPER CASE font.\n");
		printf("The first translated line is in the user's defined frame.\n");
		
		if(Preference.mini_display==FALSE)
			{
			while(var_begin<1 || var_begin>npb[var_seq]-1)
				{
				printf(" Begin map from:1< /*  %d  */ < %d:",var_min[var_seq],npb[var_seq]-1); 
				var_begin = Get_Number(var_min[var_seq]);
				}
			while(var_end<=var_begin || var_end>npb[var_seq])
				{
				printf(" end map at :%d< /*  %d  */ < %d:",var_begin,var_max[var_seq],npb[var_seq]); 
				var_end = Get_Number(var_max[var_seq]);
				}
			CLS;
			}
		else
			{
			var_begin = var_min[var_seq];
			var_end   = var_max[var_seq];
			}
		
		/**********************************************************************************/
		if(var_seq==INSERT && npb[INSERT]!=0)
			{
			printf("\nINSERT: %s\n",FICHIER_ADN[INSERT]);
			printf("\tCloning box boundaries [%d-%d] [%d-%d].\n\n",var_min[INSERT],var_min[INSERT]+sigma_5,var_max[INSERT]-sigma_3,var_max[INSERT]);
			}
		else if(var_seq==VECTOR && npb[VECTOR]!=0)
			{
			printf("\nVECTOR: %s\n",FICHIER_ADN[VECTOR]);
			printf("\tCloning box boundaries [%d-%d].\n\n",var_min[VECTOR],var_max[VECTOR]);
			}
		/**********************************************************************************/
		printf("Initialisation... \n");
		}
	else
		{
		var_begin = var_min[var_seq];
		var_end   = var_max[var_seq];
		}
	var_nbr_xy=0; /*no site defined on fragment */

	/* here, I use the field: 'select' to tag the UNIQUE sites in var_seq*/
	for(i=0;i<nbr_enzyme;i++)
		{
		Enzymes[i].select[var_seq]=FALSE; /* this enzyme is not unique by default */
		if(Enzymes[i].palindromic==TRUE)
			{
			if(Fct_N_sites(i,var_seq,1,npb[var_seq],TRUE)==1)
				Enzymes[i].select[var_seq]=TRUE; /* this enzyme is unique */
			}
		else if(i<nbr_enzyme-1)
			if(strcmp(Enzymes[i].site_complet,Enzymes[i+1].site_complet)==ARE_IDENTIC)
			{
			Enzymes[i+1].select[var_seq]=FALSE;
			if(	Fct_N_sites(i,var_seq,1,npb[var_seq],TRUE) + Fct_N_sites(i+1,var_seq,1,npb[var_seq],TRUE)==1 )
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
	/* if  a site displayed on the line n¡:LocY will be write OVER a second site, then add 1 to Loc Y , add the number of line*/
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
	 	/*fprintf(out,"\n");*/
	 	line[linelen]=EOS;
	 	fprintf(out,"%s\n",line);
	 	}
	 linelen=0;
	 /* display one line */
	 for(n=0;n<var_len_screen;n++)
		 	/*fprintf(out,"%c",((var_h[n]==TRUE)?':':' '));*/
		 	line[linelen++]=((var_h[n]==TRUE)?':':' ');
	/*fprintf(out,"\n");*/
	line[linelen]=EOS;
	 	fprintf(out,"%s\n",line);
	 linelen=0;
	 /**********************************************************************/
	 /* display one strand */
	for(j=i;(j< i+var_len_screen && j<= npb[var_seq]);j++)
	 	/*fprintf(out,"%c",sequence[var_seq][Fct_Pos(var_seq,j)]);*/
	 	line[linelen++]=sequence[var_seq][Fct_Pos(var_seq,j)];
	line[linelen]=EOS;
	 	fprintf(out,"%s",line);
	fprintf(out," \\ %d\n",Fct_Pos(var_seq,i));
	/* display ¥ each 10 bases */
	linelen=0;
	for(j=i;(j< i+var_len_screen && j<= npb[var_seq]);j++)
		{
		r = div(j,10);
		/*fprintf(out,"%c",(((int)r.rem==0)?'¥':' '));*/
		line[linelen++]=(((int)r.rem==0)?'¥':' ');
		}
	line[linelen]=EOS;
	 fprintf(out,"%s  \\\n",line);
	/* display anti strand */
	linelen=0;
	for(j=i;(j< i+var_len_screen && j<= npb[var_seq]);j++)
		/*fprintf(out,"%c",Fct_Complementaire(sequence[var_seq][Fct_Pos(var_seq,j)]));*/
		line[linelen++]=Fct_Complementaire(sequence[var_seq][Fct_Pos(var_seq,j)]);
	line[linelen]=EOS;
	fprintf(out,"%s   \\ %d\n",line,Fct_Pos(var_seq,i+var_len_screen-1));
	/* for each 3 frames 5'->3' show the resulting sequence ******************/
	
	for(k=0;k<=2;k++)
	 	{
	 	linelen=0;
	 	for(j=i+k;(j< i+var_len_screen+k && j<= npb[var_seq]);j++)
	 		{
	 		if(Fct_Frame(var_seq,Fct_Pos(var_seq,k+j))==0)
	 			line[linelen++]=Translation_at(var_seq,j);
	 			/*fprintf(out,"%c",Fct_Traduction( sequence[var_seq][Fct_Pos(var_seq,j)],
	 										sequence[var_seq][Fct_Pos(var_seq,j+1)],
	 										sequence[var_seq][Fct_Pos(var_seq,j+2)]));*/
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
	 	for(j=i+k;(j< i+var_len_screen+k && j<= npb[var_seq]);j++)
	 		{
	 		if(Fct_Frame(var_seq,Fct_Pos(var_seq,k+j))==0)
	 			line[linelen++]=Fct_Traduction( sequence[var_seq][Fct_Pos(var_seq,j)],
	 										sequence[var_seq][Fct_Pos(var_seq,j-1)],
	 										sequence[var_seq][Fct_Pos(var_seq,j-2)]);
	 			/*fprintf(out,"%c",Fct_Traduction( sequence[var_seq][Fct_Pos(var_seq,j)],
	 										sequence[var_seq][Fct_Pos(var_seq,j-1)],
	 										sequence[var_seq][Fct_Pos(var_seq,j-2)]));*/
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
	if(Preference.mini_display==FALSE && mode==FORMAT_SCREEN) {Menu(12);CLS;}
	return(TRUE);
	}
/************************************************************************/

  /************************************************************************/
 /* compare the number of site in INSERT/VECTOR and in the cloning boxes */
/************************************************************************/
boolean Intersections(FILE *out, short mode)
	{
	int IN_I,OUT_I,IN_V,OUT_V;
	register int i;
	char c,c1;
	
	if(Search_Done==FALSE) Cmd_Get_Site();
	
	if(mode==FORMAT_SCREEN)
		{
	
		if(nbr_sites==0)
			{
			Display("No site was found !!!.",3);
			BEEP;
			INKEY;
			return(FALSE);
			}


		if(npb[VECTOR]!=0)
			{
			printf("\nVECTOR: %s\n",FICHIER_ADN[VECTOR]);
			printf("\tCloning box boundaries [%d-%d].\n",var_min[VECTOR],var_max[VECTOR]);
			}
				
		if(npb[INSERT]!=0)
			{
			printf("\nINSERT: %s\n",FICHIER_ADN[INSERT]);
			printf("\tCloning box boundaries [%d-%d] [%d-%d].\n",
					var_min[INSERT],
					var_min[INSERT]+sigma_5,
					var_max[INSERT]-sigma_3,
					var_max[INSERT]);
			}
		PARAGRAF;
		if(nbr_enzyme!=0)
		printf("\nREBASE File used :\"%s\".\n",FICHIER_ENZYME);
		}


		
	fprintf(out,"X: Enzyme useful for digestion post-ligation.\n");
	fprintf(out,"+: Enzyme useful to know the INSERT direction.\n");
	if(mode!=FORMAT_HTML) 
		{
		DRAW_HR(out,'.');
		}
		
	fprintf(out,"                            _______________________________\n");
	fprintf(out,"                            |    VECTOR    |    INSERT    |\n");
	fprintf(out,"                            _______________________________\n");
	fprintf(out,"                            |  IN  |  OUT  |  IN  |  OUT  |\n");
	fprintf(out,"___________________________________________________________\n");
	
	for(i=0;i<nbr_enzyme;i++)
		{
		/* get the number of site in INSERT, IN the cloning box */
		IN_I  = Fct_N_sites(i,INSERT,var_min[INSERT],var_max[INSERT],TRUE);
		/* get the number of site in INSERT, out of the cloning box */
		OUT_I = Fct_N_sites(i,INSERT,var_min[INSERT]-1,var_max[INSERT]+1,FALSE);
		/* get the number of site in VECTOR, IN the cloning box */
		IN_V  = Fct_N_sites(i,VECTOR,var_min[VECTOR],var_max[VECTOR],TRUE);
		/* get the number of site in VECTOR, out of the cloning box */
		OUT_V = Fct_N_sites(i,VECTOR,var_min[VECTOR]-1,var_max[VECTOR]+1,FALSE);
		
		if(IN_I==0 && OUT_I>0 && IN_V>0 && OUT_V==0)
			c='X';
		else
			c=' ';
		if(IN_I==1 && OUT_I==1)
			c1='+';
		else
			c1=' ';	
		fprintf(out,"  %-8.8s [%-14.14s] |  %3d |  %3d  |  %3d |  %3d  |%c %c\n",
				Enzymes[i].nom,Enzymes[i].site_complet,
				IN_V,OUT_V,IN_I,OUT_I,c,c1);
		if(i<nbr_enzyme-1)
			{
			if(Enzymes[i].palindromic!=TRUE) i++;
			}
		}
	fprintf(out,"___________________________________________________________\n");
	fprintf(out,"                            |  IN  |  OUT  |  IN  |  OUT  |\n");
	fprintf(out,"                            _______________________________\n");
	fprintf(out,"                            |    VECTOR    |    INSERT    |\n");
	fprintf(out,"                            _______________________________\n");
	fprintf(out,"\n");
	fprintf(out,"X: Enzyme useful for digestion post-ligation.\n");
	fprintf(out,"+: Enzyme useful to know the INSERT direction.\n");
	if(mode==FORMAT_SCREEN) {Menu(12);CLS;}
	return(TRUE);
	}




/* File 'Solutions.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include "ADN.h"
#endif


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
		fprintf(out,"Enzyme(s) that could be used to verify your construction:");
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
				fprintf(out," ");
				printEnzymeName(i,mode,out);
				}
		/* if enzyme is not palindromic do not use the next one that will be the same */
		if(Enzymes[i].palindromic!=TRUE) i++;
		}
	if(vara==FALSE)
		{
		if(mode==FORMAT_HTML)
		fprintf(out,"<I>no enzyme was found.</I>");
			else
		fprintf(out,"no enzyme was found.");
		}
	fprintf(out,"\n");
	}



/* File 'Project.c' */
#ifndef __GNUC__
	#include<stdio.h>
	#include<string.h>
	#include"ADN.h"
#endif

void Cloning_Project1(void)
	{
	char	MyChoice[MAX_NOM_FICHIER];
	/***************************************/
	CLS;
	printf("  R u n   C l o n i n g   P r o j e c t  \n");
	DIRECTORY;
	printf("Cloning projects are used to handle a large number of");
	printf(" sequences with %s (for example: you just have received a new",VAR_VERSION);
	printf(" insert and you want to sub-clone it in your 153 expression ");
	printf("vectors...), for a simple manipulation, you should use the main editor.\n");
	printf("Input name of your project ? :/* %s */",DEFAULT_PROJECT);
	Cmd_lire(MyChoice);
	if(strcmp(MyChoice,"")==ARE_IDENTIC) strcpy(MyChoice,DEFAULT_PROJECT);
	Cloning_Project2(MyChoice);
	}
	
	
void Cloning_Project2(char *MyChoice)
	{
	
	FILE	*file_in;
	char	word[LARGEUR_ECRAN+LARGEUR_ECRAN+LARGEUR_ECRAN];
	int		value=0, line=0;
	STRUCT_PREFS PreferenceCpy;
	boolean arg_cip=TRUE,arg_t4=TRUE,arg_part=TRUE;
	/***************************************/
	memcpy(&PreferenceCpy,&Preference,sizeof(STRUCT_PREFS));/* copy current pref */
	Cmd_init_preferences(&Preference);
	Preference.mini_display=TRUE;
	Preference.display_messages		=FALSE;
	Preference.search_C_term		=TRUE;



	printf("\n Opening cloning project with file named \"%s\".\n",MyChoice);
	if((file_in=fopen(MyChoice,"r"))!=NULL)
		{
		printf("Please wait...\n");
		while(fgets(word,(LARGEUR_ECRAN+LARGEUR_ECRAN+LARGEUR_ECRAN),file_in)!=NULL)
			{
			if(strlen(word)<2) continue;
			if(word[strlen(word)-1]=='\n' || word[strlen(word)-1]=='\r')
				word[strlen(word)-1]='\0';
			/* add a line */
			line++;
			/* comment */
			if(strcmp(word,"")==ARE_IDENTIC)
				continue;
			if(word[0]==';')
				continue;
			  /**********/
			 /* vector */
			/**********/
			if(strcmp(word,"open.vector")==ARE_IDENTIC)
				{
				#ifdef __MAC__
					BEEP;
					printf("%s: Not available for macintosh\n",word);
				#else
					line++;
					if(fscanf(file_in,"%s",word)!=EOF)
						{strcpy(FICHIER_ADN[VECTOR],word);printf("#line %d:Opening vector \"%s\".\n",line,FICHIER_ADN[VECTOR]);Cmd_LOAD_SEQUENCE(VECTOR);}
					else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				#endif
				}
			else if(strcmp(word,"ask.vector")==ARE_IDENTIC)
				{
				#ifndef __MAC__
				DIRECTORY;
				printf("# Interaction line %d:Input VECTOR sequence name :",line);Cmd_lire(FICHIER_ADN[VECTOR]);
				#else
				printf("# Interaction line %d:Input VECTOR sequence name :\n",line);
				#endif
				Cmd_LOAD_SEQUENCE(VECTOR);
				}
			else if(strcmp(word,"set.vector.atg")==ARE_IDENTIC)
				{
				line++;
				if(fscanf(file_in,"%d",&value)!=EOF) {pos_ATG[VECTOR]=value;printf("#line %d: %s=%d\n",line,word,value);}
				else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
			else if(strcmp(word,"set.vector.min")==ARE_IDENTIC)
				{
				line++;
				if(fscanf(file_in,"%d",&value)!=EOF) {var_min[VECTOR]=value;printf("#line %d: %s=%d\n",line,word,value);}
				else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
			else if(strcmp(word,"set.vector.max")==ARE_IDENTIC)
				{
				line++;
				if(fscanf(file_in,"%d",&value)!=EOF) {var_max[VECTOR]=value;printf("#line %d: %s=%d\n",line,word,value);}
				else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
			  /**********/
			 /* insert */
			/**********/
			else if(strcmp(word,"open.insert")==ARE_IDENTIC)
				{
				#ifdef __MAC__
					BEEP;
					printf("%s: Not available for macintosh\n",word);
				#else
					line++;
					if(fscanf(file_in,"%s",word)!=EOF)
						{strcpy(FICHIER_ADN[INSERT],word);printf("#line %d:Opening insert \"%s\".\n",line,FICHIER_ADN[INSERT]);Cmd_LOAD_SEQUENCE(INSERT);}
					else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				#endif
				}
			else if(strcmp(word,"ask.insert")==ARE_IDENTIC)
				{
				#ifndef __MAC__
				DIRECTORY;
				printf("# Interaction line %d:Input INSERT sequence name :",line);Cmd_lire(FICHIER_ADN[INSERT]);
				#else
				printf("# Interaction line %d:Input INSERT sequence name :\n",line);
				#endif
				Cmd_LOAD_SEQUENCE(INSERT);
				}
			else if(strcmp(word,"set.insert.atg")==ARE_IDENTIC)
				{
				line++;
				if(fscanf(file_in,"%d",&value)!=EOF) {pos_ATG[INSERT]=value;printf("#line %d: %s=%d\n",line,word,value);}
				else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
			else if(strcmp(word,"set.insert.min")==ARE_IDENTIC)
				{
				line++;
				if(fscanf(file_in,"%d",&value)!=EOF) {var_min[INSERT]=value;printf("#line %d: %s=%d\n",line,word,value);}
				else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
			else if(strcmp(word,"set.insert.max")==ARE_IDENTIC)
				{
				line++;
				if(fscanf(file_in,"%d",&value)!=EOF) {var_max[INSERT]=value;printf("#line %d: %s=%d\n",line,word,value);}
				else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
			else if(strcmp(word,"set.insert.min.int")==ARE_IDENTIC)
				{
				line++;
				if(fscanf(file_in,"%d",&value)!=EOF) {sigma_5=value-var_min[INSERT];printf("#line %d: %s=%d\n",line,word,value);}
				else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
			else if(strcmp(word,"set.insert.max.int")==ARE_IDENTIC)
				{
				line++;
				if(fscanf(file_in,"%d",&value)!=EOF) {sigma_3=var_max[INSERT]-value;printf("#line %d: %s=%d\n",line,word,value);}
				else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
		#ifdef VAR_UNIX
			else if(strcmp(word,"align.bestfit.dna")==ARE_IDENTIC)
				{
				if(GCGBestFit(TRUE)==FALSE) {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
			else if(strcmp(word,"align.bestfit.protein")==ARE_IDENTIC)
				{
				if(GCGBestFit(FALSE)==FALSE) {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
			else if(strcmp(word,"align.gap.dna")==ARE_IDENTIC)
				{
				if(GCGGap(TRUE)==FALSE) {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
			else if(strcmp(word,"align.gap.protein")==ARE_IDENTIC)
				{
				if(GCGGap(FALSE)==FALSE) {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
		#endif
			  /**********/
			 /* rebase */
			/**********/
			else if(strcmp(word,"open.rebase")==ARE_IDENTIC)
				{
				line++;
				if(fscanf(file_in,"%s",word)!=EOF)
					{strcpy(FICHIER_ENZYME,word);printf("#line %d:Opening rebase file \"%s\".\n",line,FICHIER_ENZYME);Cmd_Get_Rebase();}
				else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
			else if(strcmp(word,"ask.rebase")==ARE_IDENTIC)
				{
				DIRECTORY;
				printf("# Interaction line %d:Input rebase file name :",line,MyChoice);Cmd_lire(FICHIER_ADN[INSERT]);
				Cmd_Get_Rebase();
				}
			  /**************/
			 /* parameters */
			/**************/
			
			else if(strcmp(word,"set.cip.true")==ARE_IDENTIC) {arg_cip=TRUE;printf("#line %d: %s\n",++line,word);}
			else if(strcmp(word,"set.cip.false")==ARE_IDENTIC) {arg_cip=FALSE;printf("#line %d: %s\n",++line,word);}
			
			else if(strcmp(word,"set.cooh.true")==ARE_IDENTIC) {Preference.side_3=TRUE;printf("#line %d: %s\n",++line,word);}
			else if(strcmp(word,"set.cooh.false")==ARE_IDENTIC) {Preference.side_3=FALSE;printf("#line %d: %s\n",++line,word);}
			
			else if(strcmp(word,"set.nh2.true")==ARE_IDENTIC) {Preference.side_5=TRUE;printf("#line %d: %s\n",++line,word);}
			else if(strcmp(word,"set.nh2.false")==ARE_IDENTIC) {Preference.side_5=FALSE;printf("#line %d: %s\n",++line,word);}
			
			else if(strcmp(word,"set.partial.true")==ARE_IDENTIC) {arg_part=TRUE;printf("#line %d: %s\n",++line,word);}
			else if(strcmp(word,"set.partial.false")==ARE_IDENTIC) {arg_part=FALSE;printf("#line %d: %s\n",++line,word);}
			
			else if(strcmp(word,"set.pol.true")==ARE_IDENTIC) {arg_t4=TRUE;printf("#line %d: %s\n",++line,word);}
			else if(strcmp(word,"set.pol.false")==ARE_IDENTIC) {arg_t4=FALSE;printf("#line %d: %s\n",++line,word);}
			
			else if(strcmp(word,"set.memory.true")==ARE_IDENTIC) {Preference.memory=TRUE;printf("#line %d: %s\n",++line,word);}
			else if(strcmp(word,"set.memory.false")==ARE_IDENTIC) {Preference.memory=FALSE;printf("#line %d: %s\n",++line,word);}

			else if(strcmp(word,"set.part.over.true")==ARE_IDENTIC) {Preference.allow_part_overhang=TRUE;printf("#line %d: %s\n",++line,word);}
			else if(strcmp(word,"set.part.over.false")==ARE_IDENTIC) {Preference.allow_part_overhang=FALSE;printf("#line %d: %s\n",++line,word);}

			else if(strcmp(word,"set.carboxy.true")==ARE_IDENTIC) {Preference.search_C_term=TRUE;printf("#line %d: %s\n",++line,word);}
			else if(strcmp(word,"set.carboxy.false")==ARE_IDENTIC) {Preference.search_C_term=FALSE;printf("#line %d: %s\n",++line,word);}


			else if(strcmp(word,"set.temperature.true")==ARE_IDENTIC) {Preference.temperature=TRUE;printf("#line %d: %s\n",++line,word);}
			else if(strcmp(word,"set.temperature.false")==ARE_IDENTIC) {Preference.temperature=FALSE;printf("#line %d: %s\n",++line,word);}

			else if(strcmp(word,"set.buffer.true")==ARE_IDENTIC) {Preference.buffer=TRUE;printf("#line %d: %s\n",++line,word);}
			else if(strcmp(word,"set.buffer.false")==ARE_IDENTIC) {Preference.buffer=FALSE;printf("#line %d: %s\n",++line,word);}


			else if(strcmp(word,"set.pctmin")==ARE_IDENTIC)
				{
				line++;
				if(fscanf(file_in,"%d",&value)!=EOF) {Preference.DeltaMin=value;printf("#line %d: %s=%d\n",line,word,value);}
				else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}
			else if(strcmp(word,"set.pctmax")==ARE_IDENTIC)
				{
				line++;
				if(fscanf(file_in,"%d",&value)!=EOF) {Preference.DeltaMax=value;printf("#line %d: %s=%d\n",line,word,value);}
				else {printf("#Project error line %d:\"%s\"\n",line,word);Menu(19);}
				}			
			  /***********/
			 /* cloneit */
			/***********/
			else  if(strcmp(word,"cloneit")==ARE_IDENTIC)
				{
				if( npb[VECTOR]==0 || pos_ATG[VECTOR]==FALSE || var_min[VECTOR]==FALSE ||
					var_max[VECTOR]==FALSE || var_max[VECTOR]==FALSE || npb[INSERT]==0 || pos_ATG[INSERT]==FALSE ||
					var_min[INSERT]==FALSE  || var_max[VECTOR]<=var_min[VECTOR] || var_max[INSERT]==FALSE || nbr_enzyme<=0 ||
					var_min[INSERT]>=var_min[INSERT]+sigma_5+1 || var_max[INSERT]<=var_min[INSERT] || var_max[INSERT]<=var_max[INSERT]-sigma_3)
					{BEEP;printf("#Syntax error.\n",word);Menu(19);}
				
				if(npb[INSERT]<=0 || npb[VECTOR]<=0)
					{Display("Sequence Error.",1);BEEP;Menu(19);}
				if(nbr_enzyme<=0)
					{Display("No enzyme was found !!!.",1);BEEP;Menu(19);}
				Cmd_Get_Site();
				if(nbr_sites>0)
					{
					Init_Enz();
					if(META_CloneIt4(arg_cip,arg_t4,arg_part)!=TRUE) printf("No solution was found.\n");
					}
				else
					{Display("No site was found !!!.",1);BEEP;Menu(19);}
				}
			  /***********/
			 /* Delta   */
			/***********/
			else  if(strcmp(word,"deltaframe")==ARE_IDENTIC)
				{
				if(var_max[INSERT]<=var_min[INSERT])
					{BEEP;printf("#Syntax error (min=%d>max=%d).\n",var_min[INSERT],var_max[INSERT]);Menu(19);}
				if(npb[INSERT]<=0)
					{Display("Sequence Error.",1);BEEP;Menu(19);}
				if(nbr_enzyme<=0)
					{Display("No enzyme was found !!!.",1);BEEP;Menu(19);}
				Cmd_Get_Site();
				if(nbr_sites>0)
					{
					/*Init_Enz();*/
					if(META_DeltaFrame(arg_t4,arg_part)!=TRUE) printf("No solution was found.\n");
					}
				else
					{Display("No site was found !!!.",1);BEEP;Menu(19);}
				}
			  /***********/
			 /* Frame   */
			/***********/
			else  if(strcmp(word,"frameshift")==ARE_IDENTIC)
				{
				if(var_max[INSERT]<=var_min[INSERT])
					{BEEP;printf("#Syntax error (min=%d>max=%d).\n",var_min[INSERT],var_max[INSERT]);Menu(19);}
				if(npb[INSERT]<=0)
					{Display("Sequence Error.",1);BEEP;Menu(19);}
				if(nbr_enzyme<=0)
					{Display("No enzyme was found !!!.",1);BEEP;Menu(19);}
				Cmd_Get_Site();
				if(nbr_sites>0)
					{
					/*Init_Enz();*/
					if(META_FrameShift(arg_part)!=TRUE) printf("No solution was found.\n");
					}
				else
					{Display("No site was found !!!.",1);BEEP;Menu(19);}
				}
			else
				printf("#Project error line %d:\"%s\"\n",line,word);
			}
		fclose(file_in);
		}
	else printf("#Project error.\n");
	printf("Closing cloning project.\n");
	memcpy(&Preference,&PreferenceCpy,sizeof(STRUCT_PREFS));
	Menu(12);
	}



/* File 'Fonctions.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <time.h>
	#include <stdlib.h>
	#include <string.h>
	#include <errno.h>
	#include "ADN.h"
	#ifdef __MAC__
		#include <console.h>
	#endif
#endif


/* draw a line **********************************/
void DRAW_HR(FILE *out,char car)
	{
	int i;
	for(i=0;i<LARGEUR_ECRAN;i++) fprintf(out,"%c",car);
	fprintf(out,"\n");
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
	if(pos_ATG[NumSeq]==FALSE)
		{
		return(r.rem);
		}
	else
		{
		s = div(pos_ATG[NumSeq],3);
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
	while( position>npb[NumSeq])	position=  position-npb[NumSeq];
	while( position<1)				position=1+position+npb[NumSeq];
	return( position);
	}

/********************************  open a file to write  ***/
FILE *Fct_Write_File(char   var_file_name[MAX_NOM_FICHIER])
	{
	char c;
	FILE *buffer;
	
	fflush(stdin);
	buffer = fopen(var_file_name, "r" );
	if (buffer==NULL)
		{
		buffer = fopen(var_file_name, "w" );
        if(buffer!=NULL)
			return(buffer);
		else
			{
			Handle_error(TRUE);
			return(NULL);
			}
		}
	else /* if the file already exist */
		{
		fclose(buffer);
		BEEP;
		printf("\t '%s' already exists !\n",var_file_name);
		printf("\t\tDo you wish to replace it ? [Yes/No]:");
		c=get_char();fflush(stdin);
		if(LOWER(c)=='y')
			{
			printf("\n\t\t\tFILE DELETED.\n");
			buffer = fopen(var_file_name, "w" );
			if(buffer!=NULL)
				return(buffer);
			else
				{
				Handle_error(TRUE);
				return(NULL);
				}
			}
		else
			{printf("\n\t\t\tCANCELED.\n");return(NULL);}
		}
	}
/*************** Cmd_lire(to read a user's word)***/
void	Cmd_lire(char *word)
	{
	char c;
	int i=0;
	fflush(stdin);
	while ((c = getchar()) != '\n' )
		{
		*word = c;
		word++;
		if(i++==(MAX_NOM_FICHIER-1))
			{
			*word = '\n';
			break;
			}
		}
	*word = 0;
	}
/** return a number , if 'carriage return' it returns _vara by default ********/
int Get_Number(int _vara)
	{
	char ASK[MAX_NOM_FICHIER];
	int var_return=0;
	int i=0;
	char *Num = "0123456789";

	Cmd_lire(ASK);
	if(strcmp(ASK,"")==ARE_IDENTIC)
		return(_vara); /* default */
	if(strspn( ASK, Num ) != strlen(ASK))
		return(FALSE); /* error */
	for(i=0;i<strlen(ASK);i++)
		var_return = var_return*10+ASK[i]-'0'; /* return ASK valor */
	return(var_return);
	}
	
void write_date(FILE *out)
	{
	time_t now;
	char *s;
	now = time( NULL );
	s = ctime( &now );
	fprintf(out,"%s", s );
	}



/* File 'Fonction_ADN.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <ctype.h>
	#include <string.h>
	#include "ADN.h"
#endif

/***returns enzyme type *****/
int Fct_type(Enzyme)
STRUCT_ENZYME Enzyme;
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
	boolean var_font=FALSE;
	
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
char Fct_Reverse(char *mot)
	{
	int l=0,i;
	char c;
	while(mot[l]!='\0') l++;
	for(i=0;i<=(int)((l-1)/2);i++)
		{
		c=mot[i];
		mot[i]=Fct_Complementaire(mot[l-i-1]);
		mot[l-i-1]=Fct_Complementaire(c);
		}
		mot[l]='\0';
	return(*mot);
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
		sequence[NumSeq][Fct_Pos(NumSeq,pos  )],
		sequence[NumSeq][Fct_Pos(NumSeq,pos+1)],
		sequence[NumSeq][Fct_Pos(NumSeq,pos+2)]));
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


short Fct_compare_AA(char aa1, char aa2)
	{
	short var_return=FALSE;
	
	aa2=LOWER(aa2);
	aa1=LOWER(aa1);
	if(aa1==aa2) return(TRUE);

		  if(strchr("GAVFPMIL",aa1)!=NULL)	return(strchr("GAVFPMIL",aa2)==NULL?FALSE:AMBIGOUS);
	else  if(strchr("DE",aa1)!=NULL)		return(strchr("DE",aa2)==NULL?FALSE:AMBIGOUS);
	else  if(strchr("KRH",aa1)!=NULL)		return(strchr("KRH",aa2)==NULL?FALSE:AMBIGOUS);
	else  if(strchr("STYCNQW",aa1)!=NULL)	return(strchr("STYCNQW",aa2)==NULL?FALSE:AMBIGOUS);
	
	BEEP;
	printf("Error Line %d file %s aa1=%c aa2=%c\n",__LINE__,__FILE__,aa1,aa2);
	return(FALSE);
	}



/* File 'Ptr_Strategy.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include <stdlib.h>
	#include "ADN.h"
#endif
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
			if(Stgy->prev==NULL) {printf("Ah ben tu vois 1  !!!!\n");}
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
	if(Stgy==NULL) {printf("Ah ben tu vois  2  !!!!\n");}
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
boolean TestStgy(void)
	{
	if(FirstStgy==NULL || LastStgy==NULL  || CurrStgy==NULL)
		return(FALSE);
	else
		return(TRUE);
	}
/***********************************************************************/
STRUCT_STRATEGY *AddStgy(STRUCT_STRATEGY *Stgy)
	{
	char line[100];
	
	if(Stgy==NULL) return(NULL);
	
	if((LastStgy==NULL || FirstStgy==NULL) && TotalSolutions !=0)
		{
		printf("ARGH ! Internal error %s %d !!!\n",__FILE__,__LINE__);
		BEEP;
		BEEP;
		INKEY;
		}
	if(LastStgy==NULL && FirstStgy==NULL)
		{
		TotalSolutions=1;
		Stgy->prev=NULL;
		Stgy->next=NULL;
		FirstStgy=Stgy;
		LastStgy=Stgy;
		CurrStgy=Stgy;
		#ifndef __MAC__
			printf("%s has found at least one solution.\nPlease wait...\n",VAR_VERSION);
		#else
			ShowProgressBar(2,FALSE,FALSE,"At least one solution was found...");
		#endif
		}
	else
		{
		LastStgy->next=Stgy;
		Stgy->prev=LastStgy;
		Stgy->next=NULL;
		LastStgy=Stgy;
		CurrStgy=Stgy;
		TotalSolutions++;
		
		sprintf(line,"%s found  %d solutions.",VAR_VERSION,TotalSolutions);
		#ifndef __MAC__
			printf("%s\n",line);
		#else
			ShowProgressBar(2,FALSE,FALSE,line);
		#endif
		if(CurrStgy->prev->next==NULL) {printf("Total =%d.next est null\n");}
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
	if(LastStgy==NULL) printf("Error (file %s line %d)\n",__FILE__,__LINE__);
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
		printf("y'a comme un couac %s %d i=%d et TotalSolutions=%d et pos=%d\n",__FILE__,__LINE__,i,TotalSolutions,*pos);
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
	boolean change=TRUE;
	

if(TotalSolutions>0)
	{

	printf("Sorting strategies...\n");

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

/* File 'CloneIt0.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "ADN.h"
#endif

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
	fprintf(out,"%s digestion.\n", FICHIER_ADN[NumSeq]);
	fprintf(out,"\n");
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
			printf("Memory Full ! NumFragment=%d\n",NumFragment);
			printf("*** ERROR ***\n");
			INKEY;
			BEEP;
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
						printf("ERROR !\n");
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
			L1=(L1>0?L1:npb[NumSeq]+L1);
			L2=(Sites[TheFragment[i+1].End].Loc- Sites[TheFragment[i+1].Start].Loc);
			L2=(L2>0?L2:npb[NumSeq]+L2);
			if( L1 < L2 )
				{
				var_change=FALSE;
				memcpy(&Topaze,&TheFragment[i],sizeof(STRUCT_FRAGMENT));
				memcpy(&TheFragment[i],&TheFragment[i+1],sizeof(STRUCT_FRAGMENT));
				memcpy(&TheFragment[i+1],&Topaze,sizeof(STRUCT_FRAGMENT));
				}
			}
		}
	
	if(mode==FORMAT_HTML) fprintf(out,"<TABLE BGCOLOR=\"#888888\">");
	for(j=0;j<NumFragment;j++)
		{
		if(mode==FORMAT_HTML) fprintf(out,"<TR><TD>");
		fprintf(out," %d%s",j+1,(mode==FORMAT_HTML?"<TD>":" "));
		var_change=(Sites[TheFragment[j].End].Loc- Sites[TheFragment[j].Start].Loc);
		fprintf(out,"%d pb%s",(var_change>0?var_change:npb[NumSeq]+var_change),(mode==FORMAT_HTML?"<TD>":"\t"));
		fprintf(out,"%s  %d - %s  %d%s",
			Enzymes[Sites[TheFragment[j].Start].NumEnz].nom,
			Sites[TheFragment[j].Start].Loc,
			Enzymes[Sites[TheFragment[j].End].NumEnz].nom,
			Sites[TheFragment[j].End].Loc,
			(mode==FORMAT_HTML?"":"\n"));
		}
	if(mode==FORMAT_HTML) fprintf(out,"</TABLE>");
	
fprintf(out,"\n");
if(AlertPartial>0)
	fprintf(out,"%sBeware this strategy needs partial digestion(s)%s",(mode==FORMAT_HTML?"<BLINK>":""),(mode==FORMAT_HTML?"</BLINK>":""));
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
	if(mode==FORMAT_HTML) fprintf(out,"<FONT COLOR=\"#888888\">");
	if(Preference.side_3==TRUE || Preference.side_5==TRUE)
		{
		var_Seq=Sites[NumSite].NumSeq;
		/* seach stop codon on the left */
		if(var_side==SIDE_5)
			{
			/* download sequence */
			for(i=Sites[NumSite].Loc-3;i>=var_limit;i--)
			   if(Fct_Frame(var_Seq,i)==IS_IN_FRAME) /* if this position is in frame */
				if(Is_Codon_Stop(Fct_Traduction(sequence[var_Seq][i],
												sequence[var_Seq][i+1],
												sequence[var_Seq][i+2]))==TRUE)
					{
					fprintf(out,"The first stop codon detected BEFORE the %s site (%d) is localized at position %d on %s. ",
						Enzymes[Sites[NumSite].NumEnz].nom,Sites[NumSite].Loc,i,((var_Seq==0)?"vector":"insert"));
					break;
					}
			}
		else if(var_side==SIDE_3) /* seach stop codon on the right */
			{
			/* upload sequence */
			for(i=Sites[NumSite].Loc+Enzymes[Sites[NumSite].NumEnz].taille_site-1;i<=var_limit;i++)
			   if(Fct_Frame(var_Seq,i)==IS_IN_FRAME) /* if this position is in frame */
				if(Is_Codon_Stop(Fct_Traduction(sequence[var_Seq][i],
												sequence[var_Seq][i+1],
												sequence[var_Seq][i+2]))==TRUE)
					{
					fprintf(out,"The first stop codon detected AFTER the %s site (%d) is localized at position %d on %s. ",
						Enzymes[Sites[NumSite].NumEnz].nom,Sites[NumSite].Loc,i,((var_Seq==VECTOR)?"vector":"insert"));
					break;
					}
			}
		}
	if(mode==FORMAT_HTML) fprintf(out,"</FONT>");
	}

			
int Fct_N_sites(int Numenz,int NumSeq,int _min,int _max,int Boolean)
	{
/************************************************************************
 this is one of the most used function in this program
	if Boolean==TRUE
		it returns the number of sites cutes by Enzyme n¡Numenz in the 
			sequence fragment between _min et _max.
	if Boolean==FALSE
		it returns the number of sites cutes by Enzyme n¡Numenz OUT OF the 
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
	printf("Initialisation... \n");
	/* scan all enzymes */
	for(i=0;i<nbr_enzyme;i++)
		{
		Enzymes[i].select[VECTOR]=TRUE; /* TRUE by default */
		Enzymes[i].select[INSERT]=TRUE; /* TRUE by default */
		/* discard if no enz in vector */
		if(Enzymes[i].palindromic==TRUE)
			{
			/* discard if no Enz in vector cloning box*/
			if(Fct_N_sites(i,VECTOR,var_min[VECTOR],var_max[VECTOR],TRUE)==0)
				Enzymes[i].select[VECTOR]=FALSE;
			/* discard if no Enz in insert cloning box */
			if( Fct_N_sites(i,INSERT,var_min[INSERT],var_min[INSERT]+sigma_5,TRUE)==0	&&
				Fct_N_sites(i,INSERT,var_max[INSERT]-sigma_3,var_max[INSERT],TRUE)==0)
				Enzymes[i].select[INSERT]=FALSE;

			}
	/* then there is the problem of the non palindromic enzymes that are coded
		twice in the array *Enzymes (See: Get_Rebase) */

		else if(strcmp(Enzymes[i].site_complet,Enzymes[i+1].site_complet)==ARE_IDENTIC)
			{
			Enzymes[i+1].select[VECTOR]=TRUE; /* TRUE by default */
			Enzymes[i+1].select[INSERT]=TRUE; /* TRUE by default */
			/* discard if no Enz in vector cloning box*/
			if(	Fct_N_sites(i,VECTOR,var_min[VECTOR],var_max[VECTOR],TRUE)==0 &&
				Fct_N_sites(i+1,VECTOR,var_min[VECTOR],var_max[VECTOR],TRUE)==0)
					{
					Enzymes[i].select[VECTOR]=FALSE;
					Enzymes[i+1].select[VECTOR]=FALSE;
					}
			/* discard if no Enz in insert cloning box*/
			if( Fct_N_sites(i,INSERT,var_min[INSERT],var_min[INSERT]+sigma_5,TRUE)==0	&&
				Fct_N_sites(i,INSERT,var_max[INSERT]-sigma_3,var_max[INSERT],TRUE)==0	&&
				Fct_N_sites(i+1,INSERT,var_min[INSERT],var_min[INSERT]+sigma_5,TRUE)==0	&&
				Fct_N_sites(i+1,INSERT,var_max[INSERT]-sigma_3,var_max[INSERT],TRUE)==0)
					{
					Enzymes[i].select[INSERT]=FALSE;
					Enzymes[i+1].select[INSERT]=FALSE;
					}
			i++; /* next enzyme have just been done */
			}
		}
	}
/*********************************************************************************/



/* File 'CloneIt3.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include "ADN.h"
#endif

/*******************************************************************/
boolean Detect_Stop(FILE *out,short mode,int NumSiteA, int NumSiteB,STRUCT_STRATEGY_2 *Stg_2)
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

		c1= sequence[SeqA][Loc5_5_3-2];
		c2= sequence[SeqA][Loc5_5_3-1];
		c3= sequence[SeqB][Loc3_5_3];
		c4= sequence[SeqB][Loc3_5_3+1];
		
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
		
		c1= sequence[SeqA][MIN(Loc5_5_3,Loc5_3_5)-2];
		c2= sequence[SeqA][MIN(Loc5_5_3,Loc5_3_5)-1];
		c3= sequence[SeqB][MIN(Loc3_5_3,Loc3_3_5)  ];
		c4= sequence[SeqB][MIN(Loc3_5_3,Loc3_3_5)+1];
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
		c1= sequence[SeqA][MAX(Loc5_5_3,Loc5_3_5)-2];
		c2= sequence[SeqA][MAX(Loc5_5_3,Loc5_3_5)-1];
		c3= sequence[SeqB][MAX(Loc3_5_3,Loc3_3_5)  ];
		c4= sequence[SeqB][MAX(Loc3_5_3,Loc3_3_5)+1];
		
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
		BEEP;
		fprintf(out,"%sBEWARE: A stop codon is created after the ligation [%c%c%c]!..%s",(mode==FORMAT_HTML?"<BLINK>":""),chara,charb,charc,(mode==FORMAT_HTML?"</BLINK>":""));
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



/* File 'Compatible.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "ADN.h"
#endif

boolean Test_Temp_and_Buffer(int NumSiteA,int NumSiteB)
	{
	if(Preference.buffer==FALSE)
		if(Buffer_Compatible2(NumSiteA,NumSiteB)==FALSE) return(FALSE);
	if(Preference.temperature==FALSE)
		if(Temperature_Compatible(NumSiteA,NumSiteB)==FALSE) return(FALSE);
	return(TRUE);
	}

boolean Temperature_Compatible(int NumSiteA,int NumSiteB)
	{
	if(Enzymes[Sites[NumSiteA].NumEnz].temperature != Enzymes[Sites[NumSiteB].NumEnz].temperature)
		return(FALSE);
	return(TRUE);
	}

short Buffer_Compatible2(int NumSiteA,int NumSiteB)
	{
	return(Buffer_Compatible1(Sites[NumSiteA].NumEnz,Sites[NumSiteB].NumEnz));
	}

short Buffer_Compatible1(int NumEnzA,int NumEnzB)
	{
	register short i;
	boolean r=FALSE,v=0;
	char buf[4]={0,0,0,0};
	/* check if information is available */
	if(Enzymes[NumEnzA].NEBBuffer=='?') return(TRUE);
	if(Enzymes[NumEnzA].NEBBuffer=='?') return(TRUE);
	if(memcmp(Enzymes[NumEnzA].NBuffer,buf,4*sizeof(char))==ARE_IDENTIC) return(TRUE);
	if(memcmp(Enzymes[NumEnzB].NBuffer,buf,4*sizeof(char))==ARE_IDENTIC) return(TRUE);
	
	for(i=0;i<NBR_BUFFER;i++)
		{
		if(Preference.buffer==FALSE)
			{
			if(Enzymes[NumEnzA].NBuffer[i] !=100) continue;
			if(Enzymes[NumEnzB].NBuffer[i] !=100) continue;
			return(i);
			}
		else
			{
			if(abs(Enzymes[NumEnzA].NBuffer[i]-Enzymes[NumEnzB].NBuffer[i])<=abs(Enzymes[NumEnzA].NBuffer[i]-Enzymes[NumEnzB].NBuffer[i]))
				{
				if(Enzymes[NumEnzA].NBuffer[i]>=Enzymes[NumEnzA].NBuffer[r] && Enzymes[NumEnzB].NBuffer[i]>=Enzymes[NumEnzB].NBuffer[r])
					r=i;
				}
			}
		}
	if(Preference.buffer==FALSE)
		return(NBR_BUFFER+1);
	else
		return(r);
	}

boolean Display_Temp_Buffer(FILE *out,short mode,int NumEnzA,int NumEnzB)
	{
	char buf[4]={0,0,0,0};
	
	/* check if information is available */
	if(Enzymes[NumEnzA].NEBBuffer=='?' || Enzymes[NumEnzB].NEBBuffer=='?') return(FALSE);
	
	if(Enzymes[NumEnzA].temperature !=37)
		{
		if(mode==FORMAT_HTML) fprintf(out,"<BLINK>");
		fprintf(out,"\nBeware: %s digests at %d T¡C.\n",Enzymes[NumEnzA].nom,Enzymes[NumEnzA].temperature);
		if(mode==FORMAT_HTML) fprintf(out,"</BLINK>");
		}
	if(Enzymes[NumEnzA].temperature !=37 && NumEnzA!=NumEnzB)
		{
		if(mode==FORMAT_HTML) fprintf(out,"<BLINK>");
		fprintf(out,"Beware: %s digests at %d T¡C.\n",Enzymes[NumEnzB].nom,Enzymes[NumEnzB].temperature);
		if(mode==FORMAT_HTML) fprintf(out,"</BLINK>");
		}
	
	if(Enzymes[NumEnzA].NEBBuffer !='?' && Enzymes[NumEnzA].NEBBuffer !=EOS)
		fprintf(out,"\n%s specific buffer is BUFFER '%c'.\n",Enzymes[NumEnzA].nom,Enzymes[NumEnzA].NEBBuffer);
	if(NumEnzA != NumEnzB  && Enzymes[NumEnzB].NEBBuffer !='?' && Enzymes[NumEnzB].NEBBuffer !=EOS)
		fprintf(out,"%s specific buffer is BUFFER '%c'.\n",Enzymes[NumEnzB].nom,Enzymes[NumEnzB].NEBBuffer);
	
	if(memcmp(Enzymes[NumEnzA].NBuffer,buf,4*sizeof(char))!=ARE_IDENTIC || memcmp(Enzymes[NumEnzB].NBuffer,buf,4*sizeof(char))!=ARE_IDENTIC)
		{
		fprintf(out,"\nPercentage acitivities in  four classical buffers.\n");
		fprintf(out,"               ------------------------\n");
		fprintf(out,"               |  1  |  2  |  3  |  4  |\n");
		fprintf(out,"---------------------------------------\n");
		if(memcmp(Enzymes[NumEnzA].NBuffer,buf,4*sizeof(char))!=ARE_IDENTIC)
		fprintf(out,"  %11s  | %3d | %3d | %3d | %3d |\n",
			Enzymes[NumEnzA].nom,
			(int)Enzymes[NumEnzA].NBuffer[0],(int)Enzymes[NumEnzA].NBuffer[1],
			(int)Enzymes[NumEnzA].NBuffer[2],(int)Enzymes[NumEnzA].NBuffer[3]);
		if(NumEnzA != NumEnzB  && memcmp(Enzymes[NumEnzB].NBuffer,buf,4*sizeof(char))!=ARE_IDENTIC)
			{
			fprintf(out,"  %11s  | %3d | %3d | %3d | %3d |\n",
			Enzymes[NumEnzB].nom,
			(int)Enzymes[NumEnzB].NBuffer[0],(int)Enzymes[NumEnzB].NBuffer[1],
			(int)Enzymes[NumEnzB].NBuffer[2],(int)Enzymes[NumEnzB].NBuffer[3]);
			}
		fprintf(out,"---------------------------------------\n");
		}
	if(NumEnzA != NumEnzB )
		{
		if(Enzymes[NumEnzA].NEBBuffer==Enzymes[NumEnzB].NEBBuffer)
			fprintf(out,"I suggest you use buffer n¡'%c'.\n",Enzymes[NumEnzA].NEBBuffer);
		else if((short)Buffer_Compatible1(NumEnzA,NumEnzB) == NBR_BUFFER+1)
			fprintf(out,"Can't find any compatible buffer.\n");
		}
	return(TRUE);
	}





int Fct_Compatible_No_Treatment (int NumSiteA,int NumSiteB, STRUCT_STRATEGY_2 *Stg_2)
	{
	/***************************************************************************/
	/* this function returns TRUE if two sites are compatible without treatment*/
	/***************************************************************************/
	int var_return=FALSE,typ1,typ2,NumEnzA,NumEnzB,SeqA,SeqB;
	boolean vara=FALSE;
	
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



/* File 'Are_in_frame.c' */

#ifndef __GNUC__
	#include <stdio.h>
	#include "ADN.h"
#endif
/********************************************************/
boolean Are_in_Frame(STRUCT_STRATEGY_2 *Stg_2)
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




/* File 'Compatible_2.c' */

#ifndef __GNUC__
	#include <stdio.h>
	#include "ADN.h"
#endif

boolean Fct_Compatible2( int SeqA, int varA_5_3, int varA_3_5,
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
				if(	Fct_Identique(	sequence[SeqA][Fct_Pos(SeqA,varA_3_5+i)],
									sequence[SeqB][Fct_Pos(SeqB,varB_3_5+i)])!=TRUE)
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
				if(	Fct_Identique(	sequence[SeqA][Fct_Pos(SeqA,varA_5_3-overhang+i)],
									sequence[SeqB][Fct_Pos(SeqB,varB_3_5+i)])!=TRUE)
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
				if(	Fct_Identique(	sequence[SeqA][Fct_Pos(SeqA,varA_5_3+i)],
									sequence[SeqB][Fct_Pos(SeqB,varB_5_3+i)])!=TRUE)
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
				if(	Fct_Identique(	sequence[SeqA][Fct_Pos(SeqA,varA_3_5-overhang+i)],
									sequence[SeqB][Fct_Pos(SeqB,varB_5_3+i)])!=TRUE)
					{vara=FALSE;break;}
				}
			return(vara);	
			}
		else return(FALSE);
		}
	else return(FALSE);
	}





/* File 'SubCloning.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include "ADN.h"
#endif


	/*************************************************************************************/                 
	/*     Sub_Cloning find all the solution to clone INSERT (delimited in the sequence
			by the sites NumSiteI5 and NumSiteI3) into VECTOR (delimited in the sequence
		by the sites NumSiteV5 and NumSiteV3)											 */
	/*************************************************************************************/
	
int Sub_Cloning(void)
	{
	/*****local variables ************************************/
	boolean IsFind=FALSE;
	int				partialV_5=0,partialV_3=0,partialI_5=0,partialI_3=0;	/* number of partial digestions for each site */
	register int	NumSiteV_5,NumSiteV_3,NumSiteI_5,NumSiteI_3;			/* sites used during the search */
	int				TreatV_5,TreatV_3,TreatI_5,TreatI_3;					/* treatment used during the search */
	STRUCT_STRATEGY Strategy;												/* cloning strategy used */
	/**********************************************************/
	Strategy.type=SUB_CLONING;/* sub cloning strategie */
	if(npb[VECTOR]==0)
		{printf("No VECTOR sequence defined !.\n");BEEP;return(FALSE);}
	if(npb[INSERT]==0)
		{printf("No INSERT sequence defined !.\n");BEEP;return(FALSE);}
	if(Preference.side_5==TRUE || Preference.side_3==TRUE)
		{
		if(pos_ATG[VECTOR]==FALSE)
			{printf("No ATG pos. defined for VECTOR !.\n");BEEP;return(FALSE);}
		if(pos_ATG[INSERT]==FALSE)
			{printf("No ATG pos. defined for INSERT!.\n");BEEP;return(FALSE);}
		}
	if(Search_Done==FALSE) Cmd_Get_Site();
	if(nbr_sites<=0)
		{
		Display("No site was found !!!.",3);
		return(FALSE);
		}
	Init_Enz();
	/**********************************************************/
	#ifndef __MAC__
		printf("Searching sub-cloning strategies,please wait...\n");
	#else
		ShowProgressBar(1,FALSE,FALSE,"");
		ShowProgressBar(2,FALSE,FALSE,"Searching sub-cloning strategies...");
	#endif
	
	/* 21/07/98 */
	DeleteAllStgys();
	/* This statement is useful for other functions to know if there is one or two sequence to handle */
	Strategy.Couple[SIDE_5].Test_Trans=TRUE;
	Strategy.Couple[SIDE_3].Test_Trans=TRUE;
	
	
	#pragma mark loop1
	/*************************************************/
	for(NumSiteV_5=0;NumSiteV_5<nbr_sites;NumSiteV_5++)
	/*************************************************/
		{
		#ifdef __MAC__
			ShowProgressBar(3,NumSiteV_5,nbr_sites-1,"");
		#endif
		
			 
		/* Is Enzyme V5 selected ? (cf InitEnz) */
		if(Enzymes[Sites[NumSiteV_5].NumEnz].select[VECTOR]==FALSE)
			continue;
		/** site V5 is on VECTOR ? **/
		if(Sites[NumSiteV_5].NumSeq != VECTOR)
			continue;
		/** site in box **/
		if(Sites[NumSiteV_5].Loc<var_min[VECTOR] || Sites[NumSiteV_5].Loc>var_max[VECTOR])
			continue;
		/** site on VECTOR everywhere in box between V5 and var_max, test partial else **/
		partialV_5=Fct_N_sites(Sites[NumSiteV_5].NumEnz,VECTOR,Sites[NumSiteV_5].Loc-1,var_max[VECTOR]+1,FALSE);
		if(partialV_5 > Preference.partial)
			continue;
		#pragma mark loop2
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
			if(Sites[NumSiteI_5].Loc<var_min[INSERT] || Sites[NumSiteI_5].Loc>var_max[INSERT])
				continue;
			/*** Site on INSERT nowhere in [var_min-sigma5][sigma3-var_max] everywhere else **/
			if( Fct_N_sites(Sites[NumSiteI_5].NumEnz,INSERT,var_min[INSERT],var_max[INSERT],TRUE)>0 	
			&&  (partialI_5=Fct_N_sites(Sites[NumSiteI_5].NumEnz,INSERT,var_min[INSERT]+sigma_5,var_max[INSERT]-sigma_3,TRUE))>Preference.partial)
				continue;			
			/** Must be before sigma5 % */
			if(Sites[NumSiteI_5].Loc>var_min[INSERT]+sigma_5)
					continue;
			#pragma mark loop3
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
				/** compatible temperature and buffer */
				if(Test_Temp_and_Buffer(NumSiteV_5,NumSiteV_3)==FALSE)
					continue;
				/** Must be on the right of V5 ONLY if V3 != V5**/
				if(NumSiteV_5!=NumSiteV_3)
					if((Sites[NumSiteV_5].Loc+Enzymes[Sites[NumSiteV_5].NumEnz].taille_site) > Sites[NumSiteV_3].Loc)
						continue;
				/** site on VECTOR everywhere in box from V5 **/
				if(Sites[NumSiteV_3].Loc<Sites[NumSiteV_5].Loc || Sites[NumSiteV_3].Loc>var_max[VECTOR])
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
				#pragma mark loop4
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
					/** compatible temperature and buffer */
					if(Test_Temp_and_Buffer(NumSiteI_5,NumSiteI_3)==FALSE)
						continue;
					/** site on INSERT **/
					if(Sites[NumSiteI_3].NumSeq != INSERT)
						continue;
					/** this site must be in the box and in sigma_3**/
					if(Sites[NumSiteI_3].Loc<(var_max[INSERT]-sigma_3) || Sites[NumSiteI_3].Loc>var_max[INSERT])
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
#pragma mark decision1
/**************************************/
/*Ligation 5' is OK without treatment */
/**************************************/

if(Fct_Test(NumSiteV_5,NumSiteI_5,NO_TREATMENT,Preference.side_5,&Strategy.Couple[SIDE_5])==TRUE)
	{
	#pragma mark decision1_1
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
		/*Pouet ajoutŽ le 21/06/98 */
		
		if(NewStrategy(&Strategy)==NULL)
			{
			printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
			INKEY;
			DeleteAllStgys();
			return(FALSE);
			}
		IsFind=TRUE;
		}
	#pragma mark decision1_2
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

		/*Pouet ajoutŽ le 21/06/98 */
		
		if(NewStrategy(&Strategy)==NULL)
			{printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
			INKEY;
			DeleteAllStgys();
			return(FALSE);
			}
		IsFind=TRUE;
		}
	}/* end ligation 5' */
#pragma mark decision2
/*****************************************/
/* Else Ligation 5' is OK with treatment */
/*****************************************/
  if(Preference.allow_all_sol==TRUE || IsFind==FALSE)

  if(Fct_Test(NumSiteV_5,NumSiteI_5,T4_TREATMENT,Preference.side_5,&Strategy.Couple[SIDE_5])==TRUE)
	{
	if(Preference.allow_T4==FALSE) /* if usage of modifying polymerases is not allowed */
		continue;
	#pragma mark decision2_1
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
		
		
		/*Pouet ajoutŽ le 21/06/98 */
		if(NewStrategy(&Strategy)==NULL)
			{printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
			INKEY;
			DeleteAllStgys();
			return(FALSE);}
		IsFind=TRUE;
		}
	#pragma mark decision2_2
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
		
		/*Pouet ajoutŽ le 21/06/98 */
		
		if(NewStrategy(&Strategy)==NULL)
			{printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
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
	
	#ifdef __MAC__
		ShowProgressBar(4,FALSE,FALSE,"");
	#else
		printf("End of Research.\n");
	#endif
	
	if(Preference.mini_display==TRUE && TotalSolutions==0) return(FALSE);
		
	if(TotalSolutions==0)
		{
		printf("No solution was found\n");
		if(Preference.display_messages==TRUE)
			{
			if(Preference.partial==0)
				printf("\t¥ May be could you allow partial digestions ?\n");
			if(Preference.allow_CIP==FALSE)
				printf("\t¥ May be could you allow usage of C.I.P. ?\n");
			if(Preference.allow_T4==FALSE)
				printf("\t¥ May be could you allow usage of modifying polymerases ?\n");
			if(Preference.memory==TRUE)
				printf("\t¥ May be could you avoid to discard short sites ?\n");
			if(Preference.allow_part_overhang==FALSE)
				printf("\t¥ May be could you allow partial overhangs ligation ?\n");
			printf("\t¥ May be could you update your Rebase file ? (see: ftp://www.neb.com/pub/rebase/ )\n");
			if(Preference.temperature==FALSE)
				printf("\t¥ May be could you allow non-compatible Temperatures ?\n");
			if(Preference.buffer==FALSE)
				printf("\t¥ May be could you allow non-compatible Buffers ?\n");

			Menu(12);
			}
		return(FALSE);
		}
	else
		{
		Cloning_solution_loop(FORMAT_SCREEN);
		return(TRUE);
		}
	}

/****************************************************************
		search what is the most simple cloning strategy
		it  starts from No C.I.A.P, No Polymerase, No partial digestion
			end with use of all those components
		return TRUE if a solution was found
		
		arg_cip==TRUE means that user allow use of T4
		arg_cip==FALSE means that user don't want to use
 ****************************************************************/
int META_CloneIt4(boolean arg_cip, boolean arg_t4,boolean arg_part)
{

/* CIP:FALSE  T4:FALSE   PARTIAL:FALSE    PART_BLUNT=TRUE*/
	Preference.allow_CIP			=FALSE;
	Preference.allow_T4				=FALSE;
	Preference.partial				=0;
	Preference.partial_only_blunt	=TRUE;
	if(Sub_Cloning()!=FALSE) return(TRUE);
	
/* CIP:FALSE  T4:TRUE   PARTIAL:FALSE    PART_BLUNT=TRUE*/
	Preference.allow_CIP			=FALSE;
	Preference.allow_T4				=TRUE;
	Preference.partial				=0;
	Preference.partial_only_blunt	=TRUE;
	if(arg_t4==TRUE)
		if(Sub_Cloning()!=FALSE) return(TRUE);
	
	
/* CIP:TRUE  T4:FALSE   PARTIAL:FALSE    PART_BLUNT=TRUE*/
	Preference.allow_CIP			=TRUE;
	Preference.allow_T4				=FALSE;
	Preference.partial				=0;
	Preference.partial_only_blunt	=TRUE;
	if(arg_cip==TRUE)
		if(Sub_Cloning()!=FALSE) return(TRUE);

/* CIP:FALSE  T4:FALSE   PARTIAL:1    PART_BLUNT=TRUE*/
	Preference.allow_CIP			=FALSE;
	Preference.allow_T4				=FALSE;
	Preference.partial				=1;
	Preference.partial_only_blunt	=TRUE;
	if(arg_part==TRUE)
		if(Sub_Cloning()!=FALSE) return(TRUE);
	
/* CIP:FALSE  T4:FALSE   PARTIAL:1    PART_BLUNT=FALSE*/
	Preference.allow_CIP			=FALSE;
	Preference.allow_T4				=FALSE;
	Preference.partial				=1;
	Preference.partial_only_blunt	=FALSE;
	if(arg_part==TRUE)
		if(Sub_Cloning()!=FALSE) return(TRUE);

/* CIP:TRUE  T4:TRUE   PARTIAL:FALSE    PART_BLUNT=TRUE*/
	Preference.allow_CIP			=TRUE;
	Preference.allow_T4				=TRUE;
	Preference.partial				=0;
	Preference.partial_only_blunt	=TRUE;
	if(arg_cip==TRUE && arg_t4==TRUE)
		if(Sub_Cloning()!=FALSE) return(TRUE);

/* CIP:TRUE  T4:TRUE   PARTIAL:FALSE    PART_BLUNT=TRUE*/
	Preference.allow_CIP			=TRUE;
	Preference.allow_T4				=TRUE;
	Preference.partial				=1;
	Preference.partial_only_blunt	=TRUE;
	if(arg_cip==TRUE && arg_t4==TRUE && arg_part==TRUE)
		if(Sub_Cloning()!=FALSE) return(TRUE);

return(TRUE);
}



/* File 'CloningSolution.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include <stdlib.h>
	#include "ADN.h"
	#ifdef __MAC__
		#include <console.h>
	#endif
#endif

void Cloning_solution_loop(short mode)
	{
	int i=0;
	STRUCT_STRATEGY *Stgy=NULL;
	char charkey;
	
	SortStgys();
	SetFirstStgy();
	while(GetCurrStgy()!=NULL)
		{
		Stgy=GetCurrStgy();
		Show_cloning_solution(stdout,Stgy,FORMAT_SCREEN);
		if(mode==FORMAT_SCREEN)
			{
			if(Preference.mini_display==TRUE)
				charkey='N';
			else
				{
				#ifdef VAR_UNIX
				        Check_Mail();
				#endif
				charkey=get_char();
				}
			switch(UPPER(charkey))
				{
				case('X'):SaveCloningSolutions(FORMAT_TEXT);break;
				case('H'):SaveCloningSolutions(FORMAT_HTML);break;
				case('N'):
				case('\n'):
				case('\r'):
				case('\0'):{
							DRAW_HR(stdout,'.');
							if(NextStgy()==NULL)
								{
								if(Preference.mini_display==TRUE)
									DeleteAllStgys();
								else
									SetFirstStgy();
								}
							}break;
				case('S'):save_it(Stgy,1);break;
				case('Q'):
				case('A'):DeleteAllStgys();SetStgy(NULL);break;
				case('D'):{DeleteStgy(GetCurrStgy());if(TestStgy()==FALSE) printf("No more Strategies !\n");}break;
				case('P'):DRAW_HR(stdout,'.');if(PrevStgy()==NULL) SetLastStgy();break;
				case('T'):	PARAGRAF;
							Display("VECTOR Pattern",1);
							Digest(FORMAT_SCREEN,stdout,VECTOR,Stgy->Couple[SIDE_5].NumSite[VECTOR],Stgy->Couple[SIDE_3].NumSite[VECTOR],Stgy->Couple[SIDE_5].Partial[VECTOR]+Stgy->Couple[SIDE_3].Partial[VECTOR]);
							PARAGRAF;Display("INSERT Pattern",1);
							Digest(FORMAT_SCREEN,stdout,INSERT,Stgy->Couple[SIDE_5].NumSite[INSERT],Stgy->Couple[SIDE_3].NumSite[INSERT],Stgy->Couple[SIDE_5].Partial[INSERT]+Stgy->Couple[SIDE_3].Partial[INSERT]);
							Menu(12);
							break;
				case('I'):Cmd_Get_Info_Choice(Stgy->Couple[SIDE_5].NumSite[VECTOR],
											Stgy->Couple[SIDE_3].NumSite[VECTOR],
											Stgy->Couple[SIDE_5].NumSite[INSERT],
											Stgy->Couple[SIDE_3].NumSite[INSERT]);break;
				case('Y'):PrintStratregy(1,FirstStgy,LastStgy);break;
				case('W'):PrintStratregy(1,Stgy,Stgy);break;
				default:DRAW_HR(stdout,'.');BEEP;break;
				}/*end switch */
			fflush(stdin);
			}/* end if mode */
		}/* end while */
	DeleteAllStgys();
	}

/*******************************************************************************************/
boolean SaveCloningSolutions(short mode)
	{
	STRUCT_STRATEGY *Stgy=NULL;
	FILE *out=NULL;
	char FileName[MAX_NOM_FICHIER];
	int i;
	
	fflush(stdin);
	DRAW_LINE;
	printf("Save thread.\n");
	switch(mode){
		case(FORMAT_TEXT):printf("(Text format)\n");break;
		case(FORMAT_HTML):printf("(HTML format)\n");break;
		default:break;
		}
		printf("\tPlease input the file name :");
	#ifdef __MAC__
		if((out=MacWriteTextFile(FileName,"\pSubCloning Strategy"))==NULL) return(FALSE);
	#else
		Cmd_lire(FileName);
		if(strcmp(FileName,"")==ARE_IDENTIC) {printf("Canceled\n");return(FALSE);}
		if((out=Fct_Write_File(FileName))==NULL) return(FALSE);
	#endif
	/*write header */
		printMainHeader(out,mode,TRUE);
	printf("\tSaving solutions...\n");
	if(mode==FORMAT_HTML)
			{
			fprintf(out,"</PRE><UL>\n");
			fprintf(out,"<LI><A HREF=\"#CLONING_SOLUTIONS\">Solutions</A>\n");
			fprintf(out,"<LI><A HREF=\"#ENZYMES_LIST\">Enzymes</A>\n");
			fprintf(out,"<LI><A HREF=\"#RESTRICTION_MAP\">Restriction Map</A>\n");
			fprintf(out,"<LI><A HREF=\"#INTERSECTIONS\">Intersections</A>\n");
			fprintf(out,"<LI><A HREF=\"#MISC\">Miscellanous</A>\n");
			fprintf(out,"</UL><PRE>\n");
			}
	/* write title 1: solutions */
	if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"CLONING_SOLUTIONS\">\n");

		printTitle(out,mode,"Cloning Solutions",1,TRUE);
				SetFirstStgy();
		while((Stgy=GetCurrStgy())!=NULL)
			{
			Show_cloning_solution(out,Stgy,mode);
			Digest(mode,out,VECTOR,Stgy->Couple[SIDE_5].NumSite[VECTOR],Stgy->Couple[SIDE_3].NumSite[VECTOR],Stgy->Couple[SIDE_5].Partial[VECTOR]+Stgy->Couple[SIDE_3].Partial[VECTOR]);
			Digest(mode,out,INSERT,Stgy->Couple[SIDE_5].NumSite[INSERT],Stgy->Couple[SIDE_3].NumSite[INSERT],Stgy->Couple[SIDE_5].Partial[INSERT]+Stgy->Couple[SIDE_3].Partial[INSERT]);
			NextStgy();
			}
		printTitle(out,mode,"",1,FALSE);
	/* write title 1:enzymes */
	printf("\tSaving Enzymes...\n");
		printTitle(out,mode,"Enzymes Used",1,TRUE);
		if(mode==FORMAT_HTML) fprintf(out,"<UL>");
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
		if(mode==FORMAT_HTML) fprintf(out,"/<UL>");
		printTitle(out,mode,"",1,FALSE);
	/* write title 1: Restriction Map */
		printf("\tSaving Restriction Map...\n");
		if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"RESTRICTION_MAP\">\n");

		printTitle(out,mode,"Restriction Maps",1,TRUE);
		/* write title 	2: VECTOR  */
			printTitle(out,mode,"VECTOR",2,TRUE);
			ResMap(out,mode,VECTOR);
			printTitle(out,mode,"",2,FALSE);
		/* write title 	2: INSERT  */
			printTitle(out,mode,"INSERT",2,TRUE);
			ResMap(out,mode,INSERT);
			printTitle(out,mode,"",2,FALSE);
		printTitle(out,mode,"",1,FALSE);
	/* write title 1: Intersection */
	printf("\tSaving Intersections...\n");
		printTitle(out,mode,"Intersections",1,TRUE);
		Intersections(out,mode);
		printTitle(out,mode,"",1,FALSE);
	
	/* write Footer */
	printMainHeader(out,mode,FALSE);
	printf("Done.\n");
	fclose(out);
	DRAW_LINE;
	return(TRUE);
	}


/***********************************************************************/
void printHeader(FILE *out,STRUCT_STRATEGY *Strategy,short mode)
	{
	char s[40];
	int i;
	CountStgy(&i);
	sprintf(s,"Solution %d/%d",i,TotalSolutions);
		switch(mode)
			{
			case(FORMAT_SCREEN):
				fprintf(out,"\n%s: %s :\n\n",VAR_VERSION,s);break;
			case(FORMAT_TEXT  ):
			case(FORMAT_HTML  ):
					if (Strategy->type==DELETION_CARBOXY && Preference.display_messages==TRUE)
						fprintf(out,"Those sites that do not necesseraly create in frame deletion but they can be used to make Carboxy-terminal deletions.\n");
					printTitle(out,mode,s,2,TRUE);
			}
		}
/***********************************************************************/

void Show_cloning_solution(FILE *out,STRUCT_STRATEGY *Strategy,short mode)
	{
	
	int vara,i;
	char s[40];
	CountStgy(&i);
	sprintf(s,"Solution %d/%d",i,TotalSolutions);
	switch(mode)
		{
		case(FORMAT_SCREEN):fprintf(out,"\n%s: %s :\n\n",VAR_VERSION,s);break;
		case(FORMAT_TEXT  ):
		case(FORMAT_HTML  ):printTitle(out,mode,s,2,TRUE);break;
		}
	
	
	/**********************************************************/
	/* display strategy for vector */
	SemiSolution(out,VECTOR,Strategy,mode);
	fprintf(out,"\n\n");
	/* display strategy for insert */
	SemiSolution(out,INSERT,Strategy,mode);
	/**********************************************************/
	if(Preference.side_5==TRUE && Preference.side_3==FALSE)
		fprintf(out,"Sites wil be in frame ligated in 5'.\n\n");
	if(Preference.side_5==FALSE && Preference.side_3==TRUE)
		fprintf(out,"Sites wil be in frame ligated in 3'.\n\n");
	if(Preference.side_5==TRUE && Preference.side_3==TRUE)
		fprintf(out,"Sites wil be in frame on both sides.\n\n");
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
			Find_stop_codon(mode,out,Strategy->Couple[SIDE_3].NumSite[VECTOR],npb[VECTOR],SIDE_3);
			}
		fprintf(out,"\n");
		}
	/**********************************************************/
	/* find sites that are localised between the 2 vector site and that could be
	used for post-ligation digestion */
	if(Strategy->Couple[SIDE_5].NumSite[VECTOR] != Strategy->Couple[SIDE_3].NumSite[VECTOR])
		{
		vara=FALSE;
		fprintf(out,"Digestion post-ligation: ");
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
				{vara=TRUE;fprintf(out,"%s ",Enzymes[i].nom);}
			}
		if(vara==TRUE)
			fprintf(out,".\n");
		else
			fprintf(out," no enzyme was found.\n");
		}
	if(mode==FORMAT_SCREEN)
		printf("\n<< (N)ext solution (P)revious (D)iscard (S)ave sequence (A)bort pa(T)tern (I)nformations Te(X)t (H)TML (Y)print thread (W)print this strategy >>\n");
	printTitle(out,mode,"",2,FALSE);
	}



/* File 'SemiCloningSol.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include "ADN.h"
#endif

void SemiSolution(FILE *out,int NumSeq,STRUCT_STRATEGY *Stg,short mode)
	{
	/*************************************************************/
	/* This procedure displays a strategy for vector OR insert   */
	/*************************************************************/
	char *string_pol[4]={"?","?","T4 DNA polymerase","Klenow DNA polymerase"};
	/***********************************************************/

	fprintf(out,"  Digest %s with ",string_seq[NumSeq]);
	if(mode==FORMAT_HTML)
		fprintf(out,"<A HREF=\"#TAG%d\">%s</A>",
			Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz,
			Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ].nom);
	else
		fprintf(out,"%s",Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ].nom);
	
	fprintf(out,"[%s] (%d",
		Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ].site_complet,
		Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].Loc);
	
	/* if site SIDE_3 isn't the same enzyme at SIDE5 (take care of asymetric enz)*/
	if(Are_Same_Asymetric(Stg->Couple[SIDE_5].NumSite[NumSeq],Stg->Couple[SIDE_3].NumSite[NumSeq])==FALSE)
		{
		fprintf(out,") and ");
		
		if(mode==FORMAT_HTML)
			fprintf(out,"<A HREF=\"#TAG%d\">%s</A>",
				Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ,
				Enzymes[ Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ].nom);
		else
			fprintf(out,"%s",Enzymes[ Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ].nom);

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
				fprintf(out,"%sYou will have to dephosphorylate your vector%s.",(mode==FORMAT_HTML?"<BLINK>":""),(mode==FORMAT_HTML?"</BLINK>":""));
	  /**********************************************/
	 /* Temperature                                */
	/**********************************************/
	

	if(Preference.buffer==FALSE)
		{
		Display_Temp_Buffer(out,mode,Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz,Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz);
		}
	else
		{
		if(Enzymes[Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz].temperature !=37)
			{
			if(mode==FORMAT_HTML) fprintf(out,"<BLINK>");
			fprintf(out,"Beware: %s digests at %d T¡C.\n",Enzymes[Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ].nom,Enzymes[Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ].temperature);
			if(mode==FORMAT_HTML) fprintf(out,"</BLINK>");
			}
		if(Enzymes[Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ].temperature !=37 && Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz !=Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz )
			{
			if(mode==FORMAT_HTML) fprintf(out,"<BLINK>");
			fprintf(out,"Beware: %s digests at %d T¡C.\n",Enzymes[Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ].nom,Enzymes[Sites[Stg->Couple[SIDE_3].NumSite[NumSeq]].NumEnz ].temperature);
			if(mode==FORMAT_HTML) fprintf(out,"</BLINK>");
			}
		}
	 /*********************************/	
	/* Test for partials digestions  */
   /*********************************/
	if(Stg->Couple[SIDE_5].Partial[NumSeq] > 0 ||  Stg->Couple[SIDE_3].Partial[NumSeq]>0 )
		{
		/* if it is the same side on side 5 and on side 3 */
		if(Are_Same_Asymetric(Stg->Couple[SIDE_5].NumSite[NumSeq],Stg->Couple[SIDE_3].NumSite[NumSeq])==TRUE)
 				fprintf(out,"%sBeware : %s [%d partial site%c]%s .\n",
 					(mode==FORMAT_HTML?"<BLINK>":""),
					Enzymes[ Sites[Stg->Couple[SIDE_5].NumSite[NumSeq]].NumEnz ].nom,
					MAX(Stg->Couple[SIDE_3].Partial[NumSeq],Stg->Couple[SIDE_5].Partial[NumSeq]),
					(MAX(Stg->Couple[SIDE_3].Partial[NumSeq],Stg->Couple[SIDE_5].Partial[NumSeq])>1?'s':'\0'),
					(mode==FORMAT_HTML?"</BLINK>":""));
		else
			{
			fprintf(out,"%sBeware :",(mode==FORMAT_HTML?"<BLINK>":""));
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
			fprintf(out,"%s.\n",(mode==FORMAT_HTML?"</BLINK>":""));
			}
		}
	fprintf(out,"\n");		
	}



/* File 'FrameShift.c' */
#ifndef __GNUC__
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include "ADN.h"
#endif

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

	int 			varb=0,val_shift=0;
	int				partial=0;
	register int	NumSite;
	double 			pctage=0.0;
	boolean 		stock5,stock3,Mem_T4;
	STRUCT_STRATEGY Strategy;
/*******************************************/
Strategy.type=FRAME_SHIFT;
if(npb[INSERT]==0)
	{printf("No INSERT sequence defined !.\n");BEEP;return(FALSE);}
if(pos_ATG[INSERT]==FALSE)
	{
	Get_ATG(INSERT);
	if(pos_ATG[INSERT]==FALSE)
		printf("No frame can be defined !.\n");BEEP;return(FALSE);
	}
if(Search_Done==FALSE) Cmd_Get_Site();
if(nbr_sites<=0)
	{
	Display("No site was found !!!.",3);
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



#ifndef __MAC__
	printf("Searching Enzyme that could induce frameshift, please wait...\n");
#else
	ShowProgressBar(1,FALSE,FALSE,"");
	ShowProgressBar(2,FALSE,FALSE,"Searching frameshifts...");
#endif


DeleteAllStgys();
#pragma mark loop1
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
	if(	Sites[NumSite].Loc<var_min[INSERT] || 
		Sites[NumSite].Loc>var_max[INSERT])
		continue;
	/** Site must be in the right tolerance percentage **/
	if( ((int)(100.0-100.0*((double)(var_max[INSERT]-Sites[NumSite].Loc)/(double)(var_max[INSERT]-var_min[INSERT])))) < Preference.DeltaMin)
		continue;
	if( ((int)(100.0-100.0*((double)(var_max[INSERT]-Sites[NumSite].Loc)/(double)(var_max[INSERT]-var_min[INSERT])))) > Preference.DeltaMax)
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
		printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
		INKEY;
		DeleteAllStgys();
		return(FALSE);
		}


	}
/***************************************************************************************************/


Preference.side_5=stock5;/* set those valors to original */
Preference.side_3=stock3;
Preference.allow_T4=Mem_T4;

if(Preference.mini_display==TRUE && TotalSolutions==0) return(FALSE);

#ifdef __MAC__
	ShowProgressBar(4,FALSE,FALSE,"");
#else
	printf("End of Research.\n");
#endif
if(TotalSolutions==FALSE)
	{
	printf("No solution was found\n");
	if(Preference.display_messages==TRUE)
		{
		if(Preference.partial==0)
			printf("\t¥ May be could you allow partial digestions ?\n");
		if(Preference.memory==TRUE)
			printf("\t¥ May be could you avoid to discard short sites ?\n");

		printf("\t¥ May be could you update your Rebase file ? (see: ftp://www.neb.com/pub/rebase/ )\n");
		Menu(12);
		}
	return(FALSE);
	}
else 
	{
	ShiftLoop(FORMAT_SCREEN);
	return(TRUE);
	}
}


int META_FrameShift(boolean arg_part)
{
/* CIP:FALSE  T4:FALSE   PARTIAL:0   PART_BLUNT=FALSE*/
	Preference.allow_T4				=FALSE;
	Preference.partial				=0;
	Preference.partial_only_blunt	=FALSE;
	if(arg_part==TRUE)
		if(FrameShift()!=FALSE) return(TRUE);


/* CIP:FALSE  T4:FALSE   PARTIAL:1    PART_BLUNT=FALSE*/
	Preference.allow_T4				=FALSE;
	Preference.partial				=1;
	Preference.partial_only_blunt	=TRUE;
	if(arg_part==TRUE)
		if(FrameShift()!=FALSE) return(TRUE);

/* CIP:FALSE  T4:FALSE   PARTIAL:TRUE    PART_BLUNT=TRUE*/
	Preference.allow_T4				=TRUE;
	Preference.partial				=1;
	Preference.partial_only_blunt	=FALSE;
	if(FrameShift()!=FALSE) return(TRUE);


return(TRUE);
}



/* File 'ShiftSolution.c' */
#ifndef __GNUC__
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include "ADN.h"
	#ifdef __MAC__
		#include <console.h>
	#endif
#endif


void ShiftLoop(short mode)
	{
	int i=0;
	STRUCT_STRATEGY *Stgy=NULL;
	char charkey;
	
	SortStgys();
	SetFirstStgy();
	while(GetCurrStgy()!=NULL)
		{
		if(TestStgy()==FALSE || TotalSolutions<=0)
			break;
		Stgy=GetCurrStgy();
		
		/* show solution */
		SolutionShift(FORMAT_SCREEN,stdout,Stgy);
		if(mode==FORMAT_SCREEN)
			{
			if(Preference.mini_display==TRUE)
				charkey='N';
			else
				{
				#ifdef VAR_UNIX
					Check_Mail();
				#endif
				charkey=get_char();
				}
							
			switch(UPPER(charkey))
				{
				case('X'): SaveShiftSolutions(FORMAT_TEXT);break;
				case('H'): SaveShiftSolutions(FORMAT_HTML);break;
				case('N'):
				case('\n'):
				case('\r'):
				case('\0'):{
							if(NextStgy()==NULL)
								{
								if(Preference.mini_display==TRUE)
									DeleteAllStgys();
								else
									SetFirstStgy();
								}
							}break;
				case('S'):save_it(Stgy,1);break;
				case('L'):Show_All_Shift(stdout,mode);break;
				case('Q'):
				case('A'):DeleteAllStgys();break;
				case('D'):{DeleteStgy(GetCurrStgy());if(TestStgy()==FALSE) printf("No more Strategies !\n");}break;
				case('P'):DRAW_HR(stdout,'.');if(PrevStgy()==NULL) SetLastStgy();break;
				case('T'):	{
							PARAGRAF;
							PARAGRAF;Display("INSERT Pattern",1);
							Digest(FORMAT_SCREEN,stdout,INSERT,Stgy->Couple[SIDE_5].NumSite[INSERT],Stgy->Couple[SIDE_5].NumSite[INSERT],Stgy->Couple[SIDE_5].Partial[INSERT]+Stgy->Couple[SIDE_5].Partial[INSERT]);
							Menu(12);
							}
							break;
				case('I'):Cmd_Get_Info_Choice(Stgy->Couple[SIDE_5].NumSite[VECTOR],
											Stgy->Couple[SIDE_3].NumSite[VECTOR],
											Stgy->Couple[SIDE_5].NumSite[INSERT],
											Stgy->Couple[SIDE_3].NumSite[INSERT]);break;
				case('Y'):PrintStratregy(3,FirstStgy,LastStgy);break;
				case('W'):PrintStratregy(3,Stgy,Stgy);break;
				default:BEEP;break;
				}/*end switch */
			fflush(stdin);
			}/* end if mode */
		}/* end while */
	DeleteAllStgys();
	}
/*******************************************************************/
void SolutionShift(short mode, FILE *out,STRUCT_STRATEGY *Strategy)
	{
	int 			i,j,varb=0,val_shift=0,overhang;
	int				partial=0,site_reconstitued;
	int				NumSite;
	double 			pctage=0.0;
	char c1;
	div_t	r;
	char s[40];
	/********************************************************************/
	CountStgy(&i);
	if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"SOL_%d\">",i);
	sprintf(s,"Solution %d/%d",i,TotalSolutions);
	switch(mode)
		{
		case(FORMAT_SCREEN): fprintf(out,"\n%s: %s :\n\n",VAR_VERSION,s);
		case(FORMAT_TEXT  ):
		case(FORMAT_HTML  ):printTitle(out,mode,s,2,TRUE);break;
		default:break;
		}
	
	
	NumSite=Strategy->Couple[SIDE_5].NumSite[VECTOR];
	partial=Strategy->Couple[SIDE_5].Partial[VECTOR];

	varb =(int)(((double)(Sites[NumSite].Loc-var_min[INSERT])/(double)(var_max[INSERT]-var_min[INSERT]))*(double)GRAPHIC);
	r = div(abs(Enzymes[Sites[NumSite].NumEnz].pos5_3-Enzymes[Sites[NumSite].NumEnz].pos3_5),3);
	val_shift=r.rem;
	pctage=((100.0*((double)(var_max[INSERT]-Sites[NumSite].Loc)/(double)(var_max[INSERT]-var_min[INSERT]))));
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
	if(mode==FORMAT_HTML) fprintf(out,"<FONT COLOR=\"#%s\">",RED_COLOR);

	fprintf(out,"  5'  ");
	for(i= Sites[NumSite].Loc - 4; i<= Sites[NumSite].Loc +PALINDROME_MAX_ENZYME;i++)
		{
		if (Fct_Frame(INSERT,i)==IS_IN_FRAME) fprintf(out,".");
		fprintf(out,"%c",sequence[INSERT][Fct_BaseShift(NumSite,Fct_Pos(INSERT,i))]);
		}
	/* write the anti strand */
	if(mode==FORMAT_HTML) 
		fprintf(out," 3'</FONT>\n  3'  <FONT COLOR=\"#%s\">",RED_COLOR);
	else
		fprintf(out," 3'\n  3'  ");
	for(i= Sites[NumSite].Loc - 4; i<= Sites[NumSite].Loc +PALINDROME_MAX_ENZYME;i++)
		{
		if (Fct_Frame(INSERT,i)==IS_IN_FRAME) fprintf(out,".");
		fprintf(out,"%c",Fct_Complementaire(sequence[INSERT][Fct_BaseShift(NumSite,Fct_Pos(INSERT,i))]));
		}
	if(mode==FORMAT_HTML) fprintf(out,"</FONT>");
	fprintf(out," 5'\n");
	/* write the proteic sequence */
	if(pos_ATG[INSERT]!=FALSE)
		{			
		fprintf(out,"  NH2 ");
		if(mode==FORMAT_HTML) fprintf(out,"<FONT COLOR=\"#%s\">",YELLOW_COLOR);
		for(i= Sites[NumSite].Loc - 4; i<= Sites[NumSite].Loc +PALINDROME_MAX_ENZYME;i++)
			if(Fct_Frame(INSERT,i)==IS_IN_FRAME)
				{
				fprintf(out," %c",(Fct_Traduction(	sequence[INSERT][Fct_BaseShift(NumSite,i)],
									sequence[INSERT][Fct_BaseShift(NumSite,i+1)],
									sequence[INSERT][Fct_BaseShift(NumSite,i+2)])));
				}
			else
				fprintf(out," ");
		fprintf(out,"  \tCOOH\n");
		if(mode==FORMAT_HTML) fprintf(out,"</FONT>");
		}
	/********************************************************************/
	if(mode==FORMAT_HTML) fprintf(out,"<UL>");
	/* write parameters */
	if(Enzymes[Sites[NumSite].NumEnz].pos3_5-Enzymes[Sites[NumSite].NumEnz].pos5_3>=0)
		{
		fprintf(out,"%s FrameShift (+ %d)\n",(mode==FORMAT_HTML?"<LI>":"\n\t¥"),val_shift);
		fprintf(out,"%s %d base%c Added.\n",(mode==FORMAT_HTML?"<LI>":"\t¥"),
		(Enzymes[Sites[NumSite].NumEnz].pos3_5-Enzymes[Sites[NumSite].NumEnz].pos5_3),
		((Enzymes[Sites[NumSite].NumEnz].pos3_5-Enzymes[Sites[NumSite].NumEnz].pos5_3)>1?'s':'\0'));
		}
	else
		{
		fprintf(out,"\n%s FrameShift (- %d)\n",(mode==FORMAT_HTML?"<LI>":"\t¥"),val_shift);
		fprintf(out,"%s %d base%c Deleted.\n",(mode==FORMAT_HTML?"<LI>":"\t¥"),
		(Enzymes[Sites[NumSite].NumEnz].taille_site-(Enzymes[Sites[NumSite].NumEnz].pos5_3-Enzymes[Sites[NumSite].NumEnz].pos3_5)),
		((Enzymes[Sites[NumSite].NumEnz].taille_site-(Enzymes[Sites[NumSite].NumEnz].pos5_3-Enzymes[Sites[NumSite].NumEnz].pos3_5))>1?'s':'\0'));
		}
	fprintf(out,"%s Site is %sreconstitued after ligation.\n",(mode==FORMAT_HTML?"<LI>":"\t¥"),((site_reconstitued==TRUE)?"":"NOT "));
	fprintf(out,"%s [%d %%] percentage of Insert.\n\n",(mode==FORMAT_HTML?"<LI>":"\t¥"),(int)pctage);
	if(mode==FORMAT_HTML) fprintf(out,"</UL>");
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
	for(i=var_min[INSERT];i<=var_max[INSERT];i++)
		{
		if(Fct_Frame(INSERT,i)==IS_IN_FRAME)
			{
			c1=Fct_Traduction(	sequence[INSERT][Fct_Pos(INSERT,Fct_BaseShift(NumSite,i  ))],
								sequence[INSERT][Fct_Pos(INSERT,Fct_BaseShift(NumSite,i+1))],
								sequence[INSERT][Fct_Pos(INSERT,Fct_BaseShift(NumSite,i+2))]);
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
	switch(mode)
		{
		case(FORMAT_SCREEN):
			printf("\n<< (N)ext solution (P)revious (D)iscard (S)ave sequence (A)bort pa(T)tern (I)nformations Te(X)t (H)TML A(L)ignments (Y)print thread (W)print this strategy>>\n");
		case(FORMAT_TEXT  ):
		case(FORMAT_HTML  ):printTitle(out,mode,s,2,FALSE);break;
		default:break;
		}
	}



/* File 'FrameLoop.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include <stdlib.h>
	#include "ADN.h"
	#ifdef __MAC__
		#include <console.h>
	#endif
#endif


void FrameLoop(short mode)
	{
	int i=0;
	STRUCT_STRATEGY *Stgy=NULL;
	char charkey;
	
	SortStgys();
	SetFirstStgy();
	#ifdef __MAC__
	/*MacBrowser(0,0);*/
	#endif
	while(GetCurrStgy()!=NULL)
		{
				if(TestStgy()==FALSE || TotalSolutions<=0)
			break;
				Stgy=GetCurrStgy();
				/* show in-frame solution */
		switch(Stgy->type)
			{
			case(DELETION_FRAME_NO):
			case(DELETION_FRAME_T4):
				{
				Solution_Frame( mode,stdout,
					Stgy->Couple[SIDE_5].NumSite[VECTOR],
					Stgy->Couple[SIDE_5].Partial[VECTOR],
					Stgy->Couple[SIDE_5].NumSite[INSERT],
					Stgy->Couple[SIDE_5].Partial[INSERT],
					(Stgy->type==DELETION_FRAME_NO?NO_TREATMENT:T4_TREATMENT));
				}
				break;
			case(DELETION_CARBOXY):
				{
				Solution_C_Terminal(mode,stdout,Stgy->Couple[SIDE_5].NumSite[VECTOR],Stgy);
				}
				break;
			default:printf("Error line %d\n",__LINE__);exit(0);break;
			}
		if(mode==FORMAT_SCREEN)
			{
			if(Preference.mini_display==TRUE)
				charkey='N';
			else
				{
				#ifdef VAR_UNIX
				        Check_Mail();
				#endif
				charkey=get_char();
				}
			switch(UPPER(charkey))
				{
				case('X'):SaveDeltaSolutions(FORMAT_TEXT);break;
				case('H'):SaveDeltaSolutions(FORMAT_HTML);break;
				case('N'):
				case('\n'):
				case('\r'):
				case('\0'):	{
							DRAW_HR(stdout,'.');
							if(NextStgy()==NULL)
								{
								if(Preference.mini_display==TRUE)
									DeleteAllStgys();
								else
									SetFirstStgy();
								}
							}break;
				case('S'):{if(Stgy->type==DELETION_FRAME_NO || Stgy->type==DELETION_FRAME_T4) save_it(Stgy,2);else BEEP;}break;
				case('Q'):
				case('A'): DeleteAllStgys();break;
				case('D'):{DeleteStgy(GetCurrStgy());if(TestStgy()==FALSE) printf("No more Strategies !\n");}break;
				case('P'):{DRAW_HR(stdout,'.');if(PrevStgy()==NULL)SetLastStgy();}break;
				case('T'):	{
							PARAGRAF;
							Display("INSERT Pattern",1);
							Digest(FORMAT_SCREEN,stdout,INSERT,Stgy->Couple[SIDE_5].NumSite[VECTOR],Stgy->Couple[SIDE_3].NumSite[INSERT],Stgy->Couple[SIDE_5].Partial[VECTOR]+Stgy->Couple[SIDE_3].Partial[INSERT]);
							Menu(12);
							}
							break;
				case('I'):Cmd_Get_Info_Choice(Stgy->Couple[SIDE_5].NumSite[VECTOR],
											Stgy->Couple[SIDE_3].NumSite[VECTOR],
											Stgy->Couple[SIDE_5].NumSite[INSERT],
											Stgy->Couple[SIDE_3].NumSite[INSERT]);break;
				case('L'):Show_All_Alignment(stdout,FORMAT_SCREEN);Menu(12);break;
				case('Y'):PrintStratregy(2,FirstStgy,LastStgy);break;
				case('W'):PrintStratregy(2,Stgy,Stgy);break;
				default:DRAW_HR(stdout,'.');BEEP;break;
				}/*end switch */
			fflush(stdin);
			}/* end if mode */
		}/* end while */
		#ifdef __MAC__
		/*MacBrowser(0,2);*/
		#endif
		DeleteAllStgys();
		}






/* File 'DeltaFrame.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include "ADN.h"
#endif


  /*=================================================================*/
 /*========= show all alignments of truncated protein           ====*/
/*=================================================================*/
	
	
	
void Show_All_Alignment(FILE *out,short mode)
	{
	int m=0,s=0,stop_codon=FALSE,var5_5_3,var3_5_3;
	char c1;
	int i=0,pos,r, partialI_5,partialI_3; 
	register int NumSiteI_5,NumSiteI_3;
	STRUCT_STRATEGY *save=NULL,*Stgy=NULL;
	boolean _TREATMENT, it_is_a_C_term_del ;
	
	save=GetCurrStgy();
	
	if(mode==FORMAT_HTML) fprintf(out,"<UL>");
	
	
		/* Scan all the fragment of INSERT */
	for(r=var_min[INSERT];r<=var_max[INSERT];r=r+(LARGEUR_ECRAN+LARGEUR_ECRAN))
		{
		/*fprintf(out,"r=%d<->%d\n",r,r+(LARGEUR_ECRAN+LARGEUR_ECRAN));*/
		fprintf(out,"%s",(mode==FORMAT_HTML?"<LI>":"\n"));
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
			if(mode==FORMAT_HTML) fprintf(out,"<LI><A HREF=\"#SOL_%d\">",pos);
			fprintf(out,"%03d:",pos);
			if(mode==FORMAT_HTML) fprintf(out,"</A>");
			stop_codon=FALSE;
			for(m=var_min[INSERT];m<=var_max[INSERT];m++)
				{
				if(Fct_Frame(INSERT,m)==IS_IN_FRAME)
					{
					/* c1 will be the amino acid displayed */
					/* if pos m is just at the level of the junction of the ligation */
					if(m==(Sites[NumSiteI_5].Loc+var5_5_3-1))
						c1=Fct_Traduction(sequence[INSERT][Sites[NumSiteI_5].Loc+var5_5_3-1],
										  sequence[INSERT][Sites[NumSiteI_3].Loc+var3_5_3],
										  sequence[INSERT][Sites[NumSiteI_3].Loc+var3_5_3+1]);
					else if(m==(Sites[NumSiteI_5].Loc+var5_5_3-2))
						c1=Fct_Traduction(sequence[INSERT][Sites[NumSiteI_5].Loc+var5_5_3-2],
										  sequence[INSERT][Sites[NumSiteI_5].Loc+var5_5_3-1],
										  sequence[INSERT][Sites[NumSiteI_3].Loc+var3_5_3]);
					/* else if m is BEFORE the junction */
					else if (m <= (Sites[NumSiteI_5].Loc+var5_5_3-1) )
						c1=Fct_Traduction(sequence[INSERT][m],sequence[INSERT][m+1],sequence[INSERT][m+2]);
					/* else if m is AFTER the junction */
					else if (m >= (Sites[NumSiteI_3].Loc+var3_5_3)   )
						c1=Fct_Traduction(sequence[INSERT][m],sequence[INSERT][m+1],sequence[INSERT][m+2]);
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
			if(m>=var_max[INSERT] &&  stop_codon==FALSE)
				fprintf(out,"[stop]...");
			NextStgy();
			fprintf(out,"\n");
			}
		}
	if(mode==FORMAT_HTML) fprintf(out,"</UL>");
	
	SetStgy(save);
	fprintf(out,"\n");
	
	}
	


  /*=================================================================*/
 /*========= save alignment of truncated protein           =========*/
/*=================================================================*/
	
void Save_alignment( int NumSiteI_5, int partialI_5,
					int NumSiteI_3,int partialI_3, 
					boolean _TREATMENT,boolean it_is_a_C_term_del ,FILE *out, short mode)
	{

	int m=0,stop_codon=FALSE,var5_5_3,var3_5_3;
	char c1;
	int i=0;
	
	/* adjust length extremities in function of treatment */
	var5_5_3=(_TREATMENT==NO_TREATMENT ? Enzymes[Sites[NumSiteI_5].NumEnz].pos5_3: Enzymes[Sites[NumSiteI_5].NumEnz].pos3_5);
	var3_5_3=(_TREATMENT==NO_TREATMENT ? Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3: Enzymes[Sites[NumSiteI_3].NumEnz].pos5_3);
	
	switch(mode)
		{
		case(FORMAT_SCREEN): fprintf(out,"Translated truncated sequence:\n");break;
		case(FORMAT_TEXT): fprintf(out,"\nTranslated truncated sequence:\n  ");break;
		case(FORMAT_HTML): fprintf(out,"<DL><DT>Translated truncated sequence:<DD>");break;
		}

	/* Scan all the fragment of INSERT */
	for(m=var_min[INSERT];m<=var_max[INSERT];m++)
		{
		if(m>=var_max[INSERT] &&  stop_codon==FALSE)
			fprintf(out,"[stop]...\n");
		else if(Fct_Frame(INSERT,m)==IS_IN_FRAME)
			{
			/* c1 will be the amino acid displayed */
			/* if pos m is just at the level of the junction of the ligation */
			if(m==(Sites[NumSiteI_5].Loc+var5_5_3-1))
				c1=Fct_Traduction(sequence[INSERT][Sites[NumSiteI_5].Loc+var5_5_3-1],
								  sequence[INSERT][Sites[NumSiteI_3].Loc+var3_5_3],
								  sequence[INSERT][Sites[NumSiteI_3].Loc+var3_5_3+1]);
			else if(m==(Sites[NumSiteI_5].Loc+var5_5_3-2))
				c1=Fct_Traduction(sequence[INSERT][Sites[NumSiteI_5].Loc+var5_5_3-2],
								  sequence[INSERT][Sites[NumSiteI_5].Loc+var5_5_3-1],
								  sequence[INSERT][Sites[NumSiteI_3].Loc+var3_5_3]);
			/* else if m is BEFORE the junction */
			else if (m <= (Sites[NumSiteI_5].Loc+var5_5_3-1) )
				c1=Fct_Traduction(sequence[INSERT][m],sequence[INSERT][m+1],sequence[INSERT][m+2]);
			/* else if m is AFTER the junction */
			else if (m >= (Sites[NumSiteI_3].Loc+var3_5_3)   )
				c1=Fct_Traduction(sequence[INSERT][m],sequence[INSERT][m+1],sequence[INSERT][m+2]);
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
	if(mode==FORMAT_HTML) fprintf(out,"</DD></DL>");
	else fprintf(out,"\n");
	}



/********************************************************************************/
	/*____________________________________________________________
	    This procedure look at in frame deletion and
		carboxy terminal deletion among all the sites in INSERT  
	____________________________________________________________*/
	
int DeltaFrame(void)
	{
	int 			i,vara,var_save=FALSE;
	double			pctage=0.0;
	int				partialI_5=0,partialI_3=0;
	register int	NumSiteI_5,NumSiteI_3;
	register int	Loop_I_5,Loop_I_3;
	boolean 		stock5,stock3,IsFind=FALSE;

	STRUCT_STRATEGY Strategy;
	/****************************************************/
	Strategy.type=DELETION_FRAME_NO;
	if(npb[INSERT]==0)
		{Display("No INSERT sequence defined !.",3);BEEP;return(FALSE);}
	if(pos_ATG[INSERT]==FALSE)
		{
		Get_ATG(INSERT);
		if(pos_ATG[INSERT]==FALSE)
			Display("No frame can be defined !.",3);BEEP;return(FALSE);}
	if(Search_Done==FALSE) Cmd_Get_Site();
	if(nbr_sites<=0)
		{
		Display("No site was found !!!.",3);
		return(FALSE);
		}
	
	
	
	/****************************************************/
	#ifndef __MAC__
		printf("Searching in frame deletions, please wait...\n");
	#else
		ShowProgressBar(1,FALSE,FALSE,"");
		ShowProgressBar(2,FALSE,FALSE,"Searching in frame deletions...");
	#endif

	DeleteAllStgys();
	
	
	stock5=Preference.side_5;   /* Frame will be important in this function so let's */
	stock3=Preference.side_3;	/* memorise the old valor and set them to TRUE */
	Preference.side_5=Preference.side_3=TRUE;
	Strategy.Couple[SIDE_5].Test_Trans=Strategy.Couple[SIDE_3].Test_Trans=FALSE;
#pragma mark loop1
/****************************************************/
for(Loop_I_5=0;Loop_I_5<nbr_sites;Loop_I_5++)
/****************************************************/
	{
	#ifdef __MAC__
			ShowProgressBar(3,Loop_I_5,nbr_sites-1,"");
	#endif
	
	/** site on INSERT **/
	if(Sites[Loop_I_5].NumSeq != INSERT)
		continue;
	/** Site on INSERT must be localized  in the box **/
	if(	Sites[Loop_I_5].Loc<var_min[INSERT] || 
		Sites[Loop_I_5].Loc>var_max[INSERT])
		continue;
	#pragma mark loop1_1
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
		if(	Sites[Loop_I_3].Loc<var_min[INSERT] || 
			Sites[Loop_I_3].Loc>var_max[INSERT])
				continue;
		/* compatible temperature and buffer */
		if(Test_Temp_and_Buffer(Loop_I_5,Loop_I_3)==FALSE)
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
				
		pctage=(double)((double)(Sites[NumSiteI_3].Loc-Sites[NumSiteI_5].Loc)/(double)(var_max[INSERT]-var_min[INSERT]))*100.0;
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
		#pragma mark decision1
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

		#pragma mark decision2
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
	
#pragma mark decision_C_term
/**** DELETIONS C TERMINAL **************************************************/
	if(Preference.search_C_term==TRUE)
		{
		/* check frameshift length */
		NumSiteI_5=Loop_I_5;
		pctage=(double)((double)(var_max[INSERT]-Sites[NumSiteI_5].Loc)/(double)(var_max[INSERT]-var_min[INSERT]))*100.0;
		if(pctage>(double)Preference.DeltaMax || pctage<(double)Preference.DeltaMin)
				continue;
		/** check partial site**/
		partialI_5=Fct_N_sites(Sites[NumSiteI_5].NumEnz,INSERT,Sites[NumSiteI_5].Loc-1,var_max[INSERT]+1,FALSE);
		if(partialI_5>Preference.partial)
			 continue;
		/* if there is a partial site, is it blunt ? */
		if( partialI_5>0 && Preference.partial_only_blunt==TRUE && Fct_type(Enzymes[Sites[NumSiteI_5].NumEnz])!=TYPE_BLUNT)
				continue;
		/* is there a site between I5 and var_max  that could be used as a second enzyme ? **********/
		vara=FALSE; /* no enzyme found */
		for(i=0;i<nbr_sites;i++)
			{
			/* compatible temperature and buffer */
			if(Test_Temp_and_Buffer(Loop_I_5,i)==FALSE)
				continue;

			if(Preference.allow_T4==FALSE)
				if(Fct_Test(NumSiteI_5,i,NO_TREATMENT,TRUE,&Strategy.Couple[SIDE_5])==FALSE)
					continue;
			if(Fct_N_sites(Sites[i].NumEnz,INSERT,
							Sites[NumSiteI_5].Loc+Enzymes[Sites[NumSiteI_5].NumEnz].taille_site,
							var_max[INSERT],
							TRUE)>0 &&
			  Fct_N_sites(Sites[i].NumEnz,INSERT,
							Sites[NumSiteI_5].Loc-1,
							var_max[INSERT]+1,
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
				printf("Out of memory (file %s line %d)\n",__FILE__,__LINE__);
				INKEY;
				DeleteAllStgys();
				return(FALSE);
				}
		}
	}
/*********************************************************************************/
Preference.side_5=stock5; /* set those valors to original */
Preference.side_3=stock3;

#ifdef __MAC__
	ShowProgressBar(4,FALSE,FALSE,"");
#else
	printf("End of Research.\n");
#endif

if(Preference.mini_display==TRUE && TotalSolutions==0) return(FALSE);

if(TotalSolutions==FALSE)
	{
	printf("No solution was found\n");
	if(Preference.display_messages==TRUE)
		{
		if(Preference.partial==0)
			printf("\t¥ May be could you allow partial digestions ?\n");
		if(Preference.allow_T4==FALSE)
			printf("\t¥ May be could you allow usage of modifying polymerases ?\n");
		if(Preference.memory==TRUE)
			printf("\t¥ May be could you avoid to discard short sites ?\n");
		if(Preference.allow_part_overhang==FALSE)
			printf("\t¥ May be could you allow partial overhangs ligation ?\n");
		if(Preference.temperature==FALSE)
			printf("\t¥ May be could you allow non-compatible Temperatures ?\n");
		if(Preference.buffer==FALSE)
			printf("\t¥ May be could you allow non-compatible Buffers ?\n");

		printf("\t¥ May be could you update your Rebase file ? (see: ftp://www.neb.com/pub/rebase/ )\n");
		Menu(12);
		}
	return(FALSE);
	}
else
	{
		FrameLoop(FORMAT_SCREEN);
		return(TRUE);
	}
}

int META_DeltaFrame(boolean arg_t4,boolean arg_part)
{

/* CIP:TRUE  T4:TRUE   PARTIAL:FALSE    PART_BLUNT=TRUE*/
	Preference.allow_T4				=TRUE;
	Preference.partial				=1;
	Preference.partial_only_blunt	=TRUE;
	if(arg_t4==TRUE && arg_part==TRUE)
		if(DeltaFrame()!=FALSE) return(TRUE);

/* CIP:FALSE  T4:FALSE   PARTIAL:1    PART_BLUNT=FALSE*/
	Preference.allow_T4				=FALSE;
	Preference.partial				=1;
	Preference.partial_only_blunt	=FALSE;
	if(arg_part==TRUE)
		if(DeltaFrame()!=FALSE) return(TRUE);


/* CIP:FALSE  T4:TRUE   PARTIAL:FALSE    PART_BLUNT=TRUE*/
	Preference.allow_T4				=TRUE;
	Preference.partial				=0;
	Preference.partial_only_blunt	=TRUE;
	if(arg_t4==TRUE)
		if(DeltaFrame()!=FALSE) return(TRUE);


/* CIP:FALSE  T4:FALSE   PARTIAL:FALSE    PART_BLUNT=TRUE*/
	Preference.allow_T4				=FALSE;
	Preference.partial				=0;
	Preference.partial_only_blunt	=TRUE;
	if(DeltaFrame()!=FALSE) return(TRUE);
	
return(TRUE);
}


/* File 'DeltaSolutions.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "ADN.h"
#endif


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
	if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"SOL_%d\">",i);
	sprintf(s,"Solution %d/%d",i,TotalSolutions);
	switch(mode)
		{
		case(FORMAT_SCREEN): fprintf(out,"\n%s: %s :\n\n",VAR_VERSION,s);break;
		case(FORMAT_TEXT  ):
		case(FORMAT_HTML  ):printTitle(out,mode,s,2,TRUE);
							fprintf(out,"Those sites that do not necesseraly create in frame deletion but they can be used to make Carboxy-terminal deletions.\n");
							break;
		default:break;
		}
	/* init strategy in order to use the function Solution_2 */
	Strategy->Couple[SIDE_5].Test_Trans=Strategy->Couple[SIDE_3].Test_Trans=FALSE;
	Strategy->Couple[SIDE_5].Treatment[VECTOR]	= DELTA_TREATMENT;
	Strategy->Couple[SIDE_3].Treatment[VECTOR]	= DELTA_TREATMENT;
	Strategy->Couple[SIDE_5].Treatment[INSERT]	= DELTA_TREATMENT;
	Strategy->Couple[SIDE_3].Treatment[INSERT]	= DELTA_TREATMENT;
	SemiSolution(out,INSERT,Strategy,mode);
	/* display figure */
	varr = (double)(var_max[INSERT]-var_min[INSERT])/(double)GRAPHIC;
	if(varr==0.0) {printf("vara=0\n");BEEP;BEEP;BEEP;INKEY;  ERROR_USER;}
	vard = (double)(Sites[NumSiteI_5].Loc-var_min[INSERT])/(varr==0?1:varr);
	varb = ((double)(var_max[INSERT]-var_min[INSERT])/(varr==0?1:varr))-(double)vard;
	varc = (double)GRAPHIC-vard-varb;
	fprintf(out,"BEWARE: Check if there is a STOP CODON after ligation.\n");
	fprintf(out,"Cloning box boundaries :[%d-%d] [%d-%d].\n",var_min[INSERT],var_min[INSERT]+sigma_5,var_max[INSERT]-sigma_3,var_max[INSERT]);
	fprintf(out,"\n\tOriginal: 5' ");
	for(i=0;i<=((int)vard+(int)varb+(int)varc);i++)
				fprintf(out,"=");
	fprintf(out," 3'\n\tDeletion: 5' ");
	for(i=0;i<(int)vard;i++) fprintf(out,"=");
	for(i=0;i<(int)varb;i++) fprintf(out,".");
	for(i=0;i<=(int)varc;i++) fprintf(out,".");
	fprintf(out," 3'\n");
	vara=var_max[INSERT]-Sites[NumSiteI_5].Loc;	
	fprintf(out,"\nEquivalent to a  deletion of about %d amino acids [%d %%]\n",
		(int)((double)vara/3.0),
		(int)(((double)(vara)/(double)(var_max[INSERT]-var_min[INSERT]))*100.0));
	/******************************************************************************/
	/* look at a second enzyme that could be used without T4*/
	fprintf(out,"¥ Second Enzyme(s) that do not need(s) polymerase modification:\n  ");
	vara=FALSE; /* no enzyme found */
	for(i=0;i<nbr_sites;i++)
		{
		/*discard if i is not compatible with NumSiteI_5 */
		if(Fct_Test(NumSiteI_5,i,NO_TREATMENT,TRUE,&Strategy->Couple[SIDE_3])==FALSE)
				continue;
		/* i is in the cloning box on the left of NumSiteI_5 , nowhere else */
		if(Fct_N_sites(Sites[i].NumEnz,INSERT,
						Sites[NumSiteI_5].Loc+Enzymes[Sites[NumSiteI_5].NumEnz].taille_site,
						var_max[INSERT],
						TRUE)>0 &&
		  Fct_N_sites(Sites[i].NumEnz,INSERT,
						Sites[NumSiteI_5].Loc-1,
						var_max[INSERT]+1,
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
		fprintf(out,"¥ Second Enzyme(s) that need(s) polymerase modification:\n  ");
		vara=FALSE; /* no enzyme found */
		for(i=0;i<nbr_sites;i++)
			{
			/* i is in the cloning box on the left of NumSiteI_5 , nowhere else */
			if(Fct_N_sites(Sites[i].NumEnz,INSERT,
							Sites[NumSiteI_5].Loc+Enzymes[Sites[NumSiteI_5].NumEnz].taille_site,
							var_max[INSERT],
							TRUE)>0 &&
			  Fct_N_sites(Sites[i].NumEnz,INSERT,
							Sites[NumSiteI_5].Loc-1,
							var_max[INSERT]+1,
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
		
	switch(mode)
		{
		case(FORMAT_SCREEN):
			printf("\n<< (N)ext solution (P)revious (D)iscard (A)bort pa(T)tern (I)nformations Te(X)t (H)TML A(L)ignments (Y)print thread (W)print this strategy >>\n");
		case(FORMAT_TEXT  ):
		case(FORMAT_HTML  ):printTitle(out,mode,s,2,FALSE);break;
		default:break;
		}
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
	if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"SOL_%d\">",i);
	sprintf(s,"Solution %d/%d",i,TotalSolutions);
	switch(mode)
		{
		case(FORMAT_SCREEN):fprintf(out,"\n%s: %s :\n\n",VAR_VERSION,s);
		case(FORMAT_TEXT  ):
		case(FORMAT_HTML  ):printTitle(out,mode,s,2,TRUE);break;
		default:break;
		}
	/* print semi solution */
	SemiSolution(out,INSERT,Strategy,mode);
	/* display figure */
	varr = (double)(var_max[INSERT]-var_min[INSERT])/(double)GRAPHIC;
	if(varr<=0.0) ERROR_USER;
	vard = (double)(Sites[NumSiteI_5].Loc-var_min[INSERT])/varr;
	varb = ((double)(Sites[NumSiteI_3].Loc-var_min[INSERT])/varr)-(double)vard;
	varc = (double)GRAPHIC-vard-varb;
	fprintf(out,"  Cloning box boundaries :[%d-%d] [%d-%d].\n",
					var_min[INSERT],
					var_min[INSERT]+sigma_5,
					var_max[INSERT]-sigma_3,
					var_max[INSERT]);
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
		(int)(((double)(vara)/(double)(var_max[INSERT]-var_min[INSERT]))*100.0));
	/* detect if ligation creates a stop codon */
	Detect_Stop(out,mode,NumSiteI_5,NumSiteI_3,&Strategy->Couple[SIDE_5]);
	/**********************************************************/	
	/* search for a post ligation digestion localised between I_5 and I_3*/
	vara=FALSE;
	fprintf(out,"Digestion post-ligation: ");
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
			{vara=TRUE;fprintf(out,"%s ",Enzymes[i].nom);}
		}
	fprintf(out,"%s",(vara==TRUE ? ".\n" : " no enzyme was found.\n"));
	/**********************************************************/
	/* where are the stop codons ? */
	Find_stop_codon(mode,out,NumSiteI_3,npb[INSERT],SIDE_3);
	fprintf(out,"\n");
	/**********************************************************/
	Save_alignment(NumSiteI_5,partialI_5,NumSiteI_3,partialI_3,_TREATMENT,FALSE,out,mode);
	/**********************************************************/
	if(Preference.display_messages==TRUE && Preference.allow_T4==FALSE)
		fprintf(out,"Beware: you have not allowed the use of modifying polymerase.\n");
	
	switch(mode)
		{
		case(FORMAT_SCREEN):
			printf("\n<< (N)ext solution (P)revious (D)iscard (S)ave sequence (A)bort pa(T)tern (I)nformations Te(X)t (H)TML A(L)ignments (Y)print thread (W)print this strategy >>\n");
		case(FORMAT_TEXT  ):
		case(FORMAT_HTML  ):printTitle(out,mode,s,2,FALSE);break;
		default:break;
		}
	free(Strategy);
	return(TRUE);
	}




/* File 'SaveShiftSolution.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include <stdlib.h>
	#include "ADN.h"
	#ifdef __MAC__
		#include <console.h>
	#endif
#endif

/*******************************************************************************************/

void Show_All_Shift(FILE *out,short mode)
	{
	int m=0,s=0;
	char c1;
	int i=0,pos,r,NumSite; 
	STRUCT_STRATEGY *save=NULL,*Stgy=NULL;
	boolean  found_stop=FALSE ;
	
	save=GetCurrStgy();
	
	if(mode==FORMAT_HTML) fprintf(out,"<UL>");
	
	
		/* Scan all the fragment of INSERT */
	for(r=var_min[INSERT];r<=var_max[INSERT];r=r+(LARGEUR_ECRAN+LARGEUR_ECRAN))
		{
		/*fprintf(out,"r=%d<->%d\n",r,r+(LARGEUR_ECRAN+LARGEUR_ECRAN));*/
		fprintf(out,"%s",(mode==FORMAT_HTML?"<LI>":"\n"));
		SetFirstStgy();
		while(GetCurrStgy()!=NULL)
			{
			Stgy=CurrStgy;
			NumSite=Stgy->Couple[SIDE_5].NumSite[VECTOR];			
			CountStgy(&pos);
			if(mode==FORMAT_HTML) fprintf(out,"<LI><A HREF=\"#SOL_%d\">",pos);
			fprintf(out,"%03d:",pos);
			if(mode==FORMAT_HTML) fprintf(out,"</A>");
			
			found_stop=FALSE;
			for(i=var_min[INSERT];i<=var_max[INSERT];i++)
				{
				if(found_stop==TRUE) break;
				if(Fct_Frame(INSERT,i)==IS_IN_FRAME)
					{
					c1=Fct_Traduction(	sequence[INSERT][Fct_Pos(INSERT,Fct_BaseShift(NumSite,i  ))],
										sequence[INSERT][Fct_Pos(INSERT,Fct_BaseShift(NumSite,i+1))],
										sequence[INSERT][Fct_Pos(INSERT,Fct_BaseShift(NumSite,i+2))]);
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
			if(i>=var_max[INSERT] &&  found_stop==FALSE)
				fprintf(out,"[stop]...");
			NextStgy();
			fprintf(out,"\n");
			}
		}
	if(mode==FORMAT_HTML) fprintf(out,"</UL>");
	SetStgy(save);
	fprintf(out,"\n");
	}



/*******************************************************************************************/
boolean SaveShiftSolutions(short mode)
	{
	STRUCT_STRATEGY *Stgy=NULL;
	FILE *out=NULL;
	char FileName[MAX_NOM_FICHIER];
	int i;
	
	fflush(stdin);
	DRAW_LINE;
	printf("Save thread.\n");
	switch(mode)
		{
		case(FORMAT_TEXT):printf("(Text format)\n");break;
		case(FORMAT_HTML):printf("(HTML format)\n");break;
		default:break;
		}
	printf("\tPlease input the file name :");
	#ifdef __MAC__
		if((out=MacWriteTextFile(FileName,"\pFrameShift Strategy"))==NULL) return(FALSE);
	#else
		Cmd_lire(FileName);
		if(strcmp(FileName,"")==ARE_IDENTIC) {printf("Canceled\n");return(FALSE);}
		if((out=Fct_Write_File(FileName))==NULL) return(FALSE);
	#endif
	/*write header */
		printMainHeader(out,mode,TRUE);
		if(mode==FORMAT_HTML)
			{
			fprintf(out,"</PRE><UL>\n");
			fprintf(out,"<LI><A HREF=\"#SHIFT_SOLUTIONS\">Solutions</A>\n");
			fprintf(out,"<LI><A HREF=\"#ALIGNMENT\">Alignments</A>\n");
			fprintf(out,"<LI><A HREF=\"#ENZYMES_LIST\">Enzymes</A>\n");
			fprintf(out,"<LI><A HREF=\"#RESTRICTION_MAP\">Restriction Map</A>\n");
			fprintf(out,"<LI><A HREF=\"#INTERSECTIONS\">Intersections</A>\n");
			fprintf(out,"<LI><A HREF=\"#MISC\">Miscellanous</A>\n");
			fprintf(out,"</UL><PRE>\n");
			}
	printf("\tSaving solutions...\n");
	/* write title 1: solutions */
		if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"SHIFT_SOLUTIONS\">\n");
		printTitle(out,mode,"FrameShifts Solutions",1,TRUE);
		SetFirstStgy();
		fprintf(out,"\n");
		while((Stgy=GetCurrStgy())!=NULL)
			{
			SolutionShift(mode,out,Stgy);
			Digest(mode,out,INSERT,Stgy->Couple[SIDE_5].NumSite[INSERT],Stgy->Couple[SIDE_5].NumSite[INSERT],Stgy->Couple[SIDE_5].Partial[INSERT]+Stgy->Couple[SIDE_5].Partial[INSERT]);
			NextStgy();
			}
		printTitle(out,mode,"",1,FALSE);
	/* write title 1: alignment */
	printf("\tSaving Alignments...\n");
	if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"ALIGNMENT\">\n");
		printTitle(out,mode,"Alignments",1,TRUE);
		Show_All_Shift(out,mode);
		printTitle(out,mode,"",1,FALSE);
	/* write title 1:enzymes */
	printf("\tSaving Enzymes...\n");
		if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"ENZYMES_LIST\">\n");
		printTitle(out,mode,"Enzymes Used",1,TRUE);
		if(mode==FORMAT_HTML) fprintf(out,"<UL>");
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
		if(mode==FORMAT_HTML) fprintf(out,"</UL>");
		printTitle(out,mode,"",1,FALSE);
	/* write title 1: Restriction Map */
	printf("\tSaving Restriction Map...\n");
		if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"RESTRICTION_MAP\">\n");
		printTitle(out,mode,"Restriction Maps",1,TRUE);
		ResMap(out,mode,INSERT);
		printTitle(out,mode,"",1,FALSE);
	/* write title 1: Intersection */
	printf("\tSaving Intersections...\n");
		if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"INTERSECTIONS\">\n");
		printTitle(out,mode,"Intersections",1,TRUE);
		Intersections(out,mode);
		printTitle(out,mode,"",1,FALSE);
	/* write Footer */
	if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"MISC\">\n");
	printMainHeader(out,mode,FALSE);
	printf("Done.\n");
	fclose(out);
	DRAW_LINE;
	return(TRUE);
	}



/* File 'ShowDeltaSolution.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include <stdlib.h>
	#include "ADN.h"
	#ifdef __MAC__
		#include <console.h>
	#endif
#endif


/*******************************************************************************************/
boolean SaveDeltaSolutions(short mode)
	{
	STRUCT_STRATEGY *Stgy=NULL;
	FILE *out=NULL;
	char FileName[MAX_NOM_FICHIER];
	int i;
	fflush(stdin);
	DRAW_LINE;
	printf("Save thread.\n");
	switch(mode)
		{
		case(FORMAT_TEXT):printf("(Text format)\n");break;
		case(FORMAT_HTML):printf("(HTML format)\n");break;
		default:break;
		}
	
	printf("\tPlease input the file name :");
	#ifdef __MAC__
		if((out=MacWriteTextFile(FileName,"\pDeletion Strategy"))==NULL) return(FALSE);
	#else
		Cmd_lire(FileName);
		if(strcmp(FileName,"")==ARE_IDENTIC) {printf("Canceled\n");return(FALSE);}
		if((out=Fct_Write_File(FileName))==NULL) return(FALSE);
	#endif
	/*write header */
		printMainHeader(out,mode,TRUE);
		if(mode==FORMAT_HTML)
			{
			fprintf(out,"</PRE><UL>\n");
			fprintf(out,"<LI><A HREF=\"#FRAME_SOLUTIONS\">Solutions</A>\n");
			fprintf(out,"<LI><A HREF=\"#ALIGNMENT\">Alignments</A>\n");
			fprintf(out,"<LI><A HREF=\"#ENZYMES_LIST\">Enzymes</A>\n");
			fprintf(out,"<LI><A HREF=\"#RESTRICTION_MAP\">Restriction Map</A>\n");
			fprintf(out,"<LI><A HREF=\"#INTERSECTIONS\">Intersections</A>\n");
			fprintf(out,"<LI><A HREF=\"#MISC\">Miscellanous</A>\n");
			fprintf(out,"</UL><PRE>\n");
			}
	printf("\tSaving solutions...\n");
	/* write title 1: solutions */
		if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"FRAME_SOLUTIONS\">\n");
		printTitle(out,mode,"In-Frame Deletions Solutions",1,TRUE);
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
		printTitle(out,mode,"",1,FALSE);
	/* write title 1: alignment */
	printf("\tSaving Alignments...\n");
	if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"ALIGNMENT\">\n");
		printTitle(out,mode,"Alignments",1,TRUE);
		Show_All_Alignment(out,mode);
		printTitle(out,mode,"",1,FALSE);
	/* write title 1:enzymes */
	printf("\tSaving Enzymes...\n");
		if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"ENZYMES_LIST\">\n");
		printTitle(out,mode,"Enzymes Used",1,TRUE);
		if(mode==FORMAT_HTML) fprintf(out,"<UL>");
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
		if(mode==FORMAT_HTML) fprintf(out,"</UL>");
		printTitle(out,mode,"",1,FALSE);
	/* write title 1: Restriction Map */
	printf("\tSaving Restriction Map...\n");
		if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"RESTRICTION_MAP\">\n");
		printTitle(out,mode,"Restriction Maps",1,TRUE);
		ResMap(out,mode,INSERT);
		printTitle(out,mode,"",1,FALSE);
	/* write title 1: Intersection */
	printf("\tSaving Intersections...\n");
		if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"INTERSECTIONS\">\n");
		printTitle(out,mode,"Intersections",1,TRUE);
		Intersections(out,mode);
		printTitle(out,mode,"",1,FALSE);
	/* write Footer */
	if(mode==FORMAT_HTML) fprintf(out,"<A NAME=\"MISC\">\n");
	printMainHeader(out,mode,FALSE);
	printf("Done.\n");
	fclose(out);
	DRAW_LINE;
	return(TRUE);
	}



/* File 'UNIX.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include "ADN.h"
#endif

#ifdef VAR_UNIX

boolean  Check_Mail(void)
	{
	boolean maillen=FALSE;
	FILE *in;
	if((in=fopen(getenv("MAIL"),"r"))==NULL) return(FALSE);
	if(fgetc(in)!=EOF)
		{
		maillen=TRUE;
		}
	if(maillen==TRUE) printf("You have mail.\n");
	fclose(in);
	return(maillen);
	}

boolean EditNewSequence(void)
	{
	char FileName[FILENAME_MAX];
	char line[100];
	FILE *out=NULL;
	printf("The pico text editor will be used to edit your new sequence./n");
	printf("Please, type your new sequence name:");
	Cmd_lire(FileName);
	if(strcmp(FileName,"")==ARE_IDENTIC)
		{printf("Canceled\n");return(FALSE);}
	if((out=Fct_Write_File(FileName))==NULL)
		return(FALSE);
	fprintf(out,">%s Pierre LINDENBAUM 1998\n",VAR_VERSION);
	fprintf(out,">Name:%s\n",FileName);
	fprintf(out,">Date:");
	write_date(out);
	fprintf(out,">Description:?\n");
	fprintf(out,">..\n");
	fclose(out);
	sprintf(line,"pico +6 -n60 -w %s",FileName);
	printf("Send -%s- to the system.\n",line);
	system(line);
	return(TRUE);
	}

boolean EditSequence(int NumSeq)
	{
	char line[100];
	if(IsStriderSeq[NumSeq]==TRUE) return(FALSE);
	sprintf(line,"pico -n60 -w -t %s",FICHIER_ADN[NumSeq]);
	printf("Send -%s- to the system.\n",line);
	system(line);
	printf("Now, re-loading \"%s\"\n",FICHIER_ADN[NumSeq]);
	Cmd_LOAD_SEQUENCE(NumSeq);
	if(Preference.display_messages==TRUE)
		{
		if(npb[NumSeq]>0)
			{
			printf("\nThe cloning box(es) boundaries and the ATG position have been searched, please check the displayed datas.\n");
			}
		#ifndef __MAC__
			Menu(12);
		#endif
		}
	return(TRUE);
	}

boolean ShowRebase(void)
	{
	char line[100];
	sprintf(line,"pico -n60 -v -w %s",FICHIER_ENZYME);
	printf("Send -%s- to the system.\n",line);
	system(line);
	return(TRUE);
	}

boolean GotoLynx(char *URL)
	{
	char line[1024];
	sprintf(line,"lynx %s",URL);
	printf("Send -%s- to the system.\n",line);
	system(line);
	return(TRUE);
	}
#endif



/* File 'GCG.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include <stdlib.h>
	#ifdef VAR_UNIX
		#include<unistd.h>
	#endif
	#include "ADN.h"
#endif

#ifdef VAR_UNIX

void GCGReformat(char *SeqName)
	{
	char line[1024];
	sprintf(line,"%s/reformat %s -NOCOM",getenv("GCGUTILDIR"),SeqName);
	system(line);
	}

boolean exportTheSequence(char *SeqName, int NumSeq, boolean SeqIsDNA)
	{
	FILE *out;
	int i,j=0;
	char c;
	if((out=fopen(SeqName,"w"))==NULL) return(FALSE);
	fprintf(out,"\n..\n");
	for(i=var_min[NumSeq];i<=var_max[NumSeq];i++)
		{
		if(j++==LARGEUR_ECRAN) {j=0;fprintf(out,"\n");}
		if(SeqIsDNA==TRUE)
			{
			fprintf(out,"%c",sequence[NumSeq][Fct_Pos(NumSeq,i)]);
			}
		else if(Fct_Frame(NumSeq,Fct_Pos(NumSeq,i))==IS_IN_FRAME)
			{	
			c=Translation_at(NumSeq,i);
			c=(Is_Codon_Stop(c)?'*':c);
			fprintf(out,"%c",c);
			}
		}
	fprintf(out,"\n");
	fclose(out);
	GCGReformat(SeqName);
	return(TRUE);
	}


boolean GCGBestFit(boolean SeqIsDNA)
	{
	char line[1024];
	
	if(getenv("GCGUTILDIR")==NULL) return(FALSE);
	if(npb[VECTOR]<=0 || npb[INSERT]<=0) return(FALSE);
	if(Preference.display_messages==TRUE)
		{
		printf("BestFit makes an optimal alignment of the best segment of similarity between the defined two cloning boxes. Optimal alignments are found by inserting gaps to maximize the number of matches.");
		printf("BestFit inserts gaps to obtain the optimal alignment of the best region of similarity between two sequences. BestFit is the most powerful method for identifying the best similarity between two sequences whose relationship is unknown.\n");
		}
	printf("BestFit");Menu(20);
	
	
	
	if(exportTheSequence(VECTOR_TEMP,VECTOR, SeqIsDNA)==FALSE) {printf("Error %d %s \n",__LINE__,__FILE__);return(FALSE);}
	if(exportTheSequence(INSERT_TEMP,INSERT, SeqIsDNA)==FALSE) {printf("Error %d %s \n",__LINE__,__FILE__);return(FALSE);}
	
	sprintf(line,"%s/gap -ProgramName=bestfit %s %s -BEGIN1=1 -END1=%d -BEGIN2=1 -END2=%d -NOREV1 -NOREV2 -OUTFile1=%s",
		getenv("GCGUTILDIR"),VECTOR_TEMP,INSERT_TEMP,var_max[VECTOR]-var_min[VECTOR],var_max[INSERT]-var_min[INSERT],OUT_TEMP);
	printf("--%s--\n",line);
	system(line);
	sprintf(line,"more %s",OUT_TEMP);
	system(line);
	if(remove(VECTOR_TEMP)!=0) fprintf(stderr,"Alert: can't remove %s\n",VECTOR_TEMP);
	if(remove(INSERT_TEMP)!=0) fprintf(stderr,"Alert: can't remove %s\n",INSERT_TEMP);
	remove(OUT_TEMP)!=0;

	return(TRUE);
	}

boolean GCGGap(boolean SeqIsDNA)
	{
	short i;
	char line[1024];

	if(getenv("GCGUTILDIR")==NULL) return(FALSE);
	if(npb[VECTOR]<=0 || npb[INSERT]<=0) return(FALSE);
	if(Preference.display_messages==TRUE)
		{
		printf("Gap finds the alignment of two complete sequences that maximizes the number of matches and minimizes the number of gaps. ");
		printf("Gap considers all possible alignments and gap positions and creates the alignment with the largest number of matched bases and the fewest gaps.\n");
		}
	printf("Gap");Menu(20);

	if(exportTheSequence(VECTOR_TEMP,VECTOR, SeqIsDNA)==FALSE) {printf("Error %d %s \n",__LINE__,__FILE__);return(FALSE);}
	if(exportTheSequence(INSERT_TEMP,INSERT, SeqIsDNA)==FALSE) {printf("Error %d %s \n",__LINE__,__FILE__);return(FALSE);}

	sprintf(line,"%s/gap -INfile1=%s -INfile2=%s -BEGIN1=1 -END1=%d  -BEGIN2=1 -END2=%d -NOREV1 -NOREV2 -OUTfile1=%s -NOSUM",
		getenv("GCGUTILDIR"),VECTOR_TEMP,INSERT_TEMP,var_max[VECTOR]-var_min[VECTOR],var_max[INSERT]-var_min[INSERT],OUT_TEMP);
	printf("%s\n",line);
	system(line);


	sprintf(line,"more %s",OUT_TEMP);
	system(line);

	if(remove(VECTOR_TEMP)!=0) fprintf(stderr,"Alert: can't remove %s\n",VECTOR_TEMP);
	if(remove(INSERT_TEMP)!=0) fprintf(stderr,"Alert: can't remove %s\n",INSERT_TEMP);
	remove(OUT_TEMP);
	return(TRUE);
	}
#endif



/* File 'Mutagenesis.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "ADN.h"
	#ifdef __MAC__
		#include <console.h>
	#endif
#endif

boolean Fct_coded(char triplet1[4], char triplet2[4], boolean *x)
	{
	short i;
	boolean y=TRUE;
	/*printf("entree x=%d codon1=%s codon2=%s\n",*x,triplet1,triplet2);*/
	if(*x==TRUE) return(TRUE);
	y=TRUE;
	for(i=0;i<3;i++)
		{
		/*printf("%c-%c =%d\n",triplet1[i],triplet2[i],Fct_Identique(triplet1[i],triplet2[i]));*/
		if(Fct_Identique(triplet1[i],triplet2[i])!=TRUE) y=FALSE;
		}
	if(y==FALSE)
		*x=FALSE;
	else
		*x=TRUE;
	/*if(strcmp(triplet1,triplet1)==ARE_IDENTIC) INKEY;*/
	/*printf("sortie x=%d\n",*x);*/
	return(*x);
	}

boolean IsAACodedBy(char TheAA,char c1,char c2,char c3)
	{
	char codon[4];
	boolean x=FALSE;
	
	codon[0]=c1;
	codon[1]=c2;
	codon[2]=c3;
	codon[3]=EOS;
	/*printf("entree x=%d codon=%s\n",x,codon);*/
	
	switch(UPPER(TheAA))
		{
		case('!'):	Fct_coded("TGA",&codon[0],&x);break;
		case('$'):	Fct_coded("TAG",&codon[0],&x);break;
		case('*'):	Fct_coded("TAA",&codon[0],&x);break;
		case('A'):	Fct_coded("GCA",&codon[0],&x);
					Fct_coded("GCC",&codon[0],&x);
					Fct_coded("GCG",&codon[0],&x);
					Fct_coded("GCT",&codon[0],&x);break;
		case('C'):	Fct_coded("TGC",&codon[0],&x);
					Fct_coded("TGT",&codon[0],&x);break;
		case('D'):	Fct_coded("GAC",&codon[0],&x);
					Fct_coded("GAT",&codon[0],&x);break;
		case('E'):	Fct_coded("GAA",&codon[0],&x);
					Fct_coded("GAG",&codon[0],&x);break;
		case('F'):	Fct_coded("TTC",&codon[0],&x);
					Fct_coded("TTT",&codon[0],&x);break;
		case('G'):	Fct_coded("GGA",&codon[0],&x);
					Fct_coded("GGC",&codon[0],&x);
					Fct_coded("GGG",&codon[0],&x);
					Fct_coded("GGT",&codon[0],&x);break;
		case('H'):	Fct_coded("CAC",&codon[0],&x);
					Fct_coded("CAT",&codon[0],&x);break;
		case('I'):	Fct_coded("ATA",&codon[0],&x);
					Fct_coded("ATC",&codon[0],&x);
					Fct_coded("ATT",&codon[0],&x);break;
		case('K'):	Fct_coded("AAA",&codon[0],&x);
					Fct_coded("AAG",&codon[0],&x);break;
		case('L'):	Fct_coded("CTA",&codon[0],&x);
					Fct_coded("CTC",&codon[0],&x);
					Fct_coded("CTG",&codon[0],&x);
					Fct_coded("CTT",&codon[0],&x);
					Fct_coded("TTA",&codon[0],&x);
					Fct_coded("TTG",&codon[0],&x);break;
		case('M'):	Fct_coded("ATG",&codon[0],&x);break;
		case('N'):	Fct_coded("AAC",&codon[0],&x);
					Fct_coded("AAT",&codon[0],&x);break;
		case('P'):	Fct_coded("CCA",&codon[0],&x);
					Fct_coded("CCC",&codon[0],&x);
					Fct_coded("CCG",&codon[0],&x);
					Fct_coded("CCT",&codon[0],&x);break;
		case('Q'):	Fct_coded("CAA",&codon[0],&x);
					Fct_coded("CAG",&codon[0],&x);break;
		case('R'):	Fct_coded("AGA",&codon[0],&x);
					Fct_coded("AGG",&codon[0],&x);
					Fct_coded("CGA",&codon[0],&x);
					Fct_coded("CGC",&codon[0],&x);
					Fct_coded("CGG",&codon[0],&x);
					Fct_coded("CGT",&codon[0],&x);break;
		case('S'):	Fct_coded("AGC",&codon[0],&x);
					Fct_coded("AGT",&codon[0],&x);
					Fct_coded("TCA",&codon[0],&x);
					Fct_coded("TCC",&codon[0],&x);
					Fct_coded("TCG",&codon[0],&x);
					Fct_coded("TCT",&codon[0],&x);break;
		case('T'):	Fct_coded("ACA",&codon[0],&x);
					Fct_coded("ACC",&codon[0],&x);
					Fct_coded("ACG",&codon[0],&x);
					Fct_coded("ACT",&codon[0],&x);
					Fct_coded("TGG",&codon[0],&x);break;
		case('V'):	Fct_coded("GTA",&codon[0],&x);
					Fct_coded("GTC",&codon[0],&x);
					Fct_coded("GTG",&codon[0],&x);
					Fct_coded("GTT",&codon[0],&x);break;
		case('Y'):	Fct_coded("TAC",&codon[0],&x);
					Fct_coded("TAT",&codon[0],&x);break;
		
		default:ERROR_USER;break;
		}
	/*if(x==TRUE)
	{	printf("sortie switch x=%d\n",x);INKEY;}*/
	return(x);
	}




void draw_seq_aound_aa(FILE *out,const short mode,const int pos)
	{
	int i;
	char c;
	
	fprintf(out,"\n    ");
	for(i=pos-20;i<=pos+20;i++)
		if(Fct_Frame(INSERT,Fct_Pos(INSERT,i))==IS_IN_FRAME)
		{
		c=Translation_at(INSERT,i);
		if(i!=pos && i!=pos+1 && i!=pos+2) c=LOWER(c);
		fprintf(out,"  %c ",c);
		}
	
	fprintf(out,"\n5' ");
	for(i=pos-20;i<=pos+20;i++)
		{
		if(Fct_Frame(INSERT,Fct_Pos(INSERT,i))==IS_IN_FRAME) fprintf(out," ");
		c=sequence[INSERT][Fct_Pos(INSERT,i)];
		if(i!=pos && i!=pos+1 && i!=pos+2) c=LOWER(c);
		fprintf(out,"%c",c);
		}
	fprintf(out," 3'\n");
	
	fprintf(out,"3' ");
	for(i=pos-20;i<=pos+20;i++)
		{
		if(Fct_Frame(INSERT,Fct_Pos(INSERT,i))==IS_IN_FRAME) fprintf(out," ");
		c=Fct_Complementaire(sequence[INSERT][Fct_Pos(INSERT,i)]);
		if(i!=pos && i!=pos+1 && i!=pos+2) c=LOWER(c);
		fprintf(out,"%c",c);
		}
	fprintf(out," 5'\n");
	}
boolean ChooseTheAA(void)
	{
	boolean done=FALSE;
	char charkey;
	
	if(npb[INSERT]<=0) {BEEP;return(FALSE);}
	CLS;
	
	if(pos_of_the_AA==0)
		{
		pos_of_the_AA=(int)(((float)(float)var_min[INSERT]+(float)var_max[INSERT])/2.0);
		while(Fct_Frame(INSERT,pos_of_the_AA)!=IS_IN_FRAME)
			pos_of_the_AA = Fct_Pos(INSERT,pos_of_the_AA+1);
		}
	
	while(done==FALSE)
		{
		TheAminoAcid=Translation_at(INSERT,pos_of_the_AA);
		printf("The AminoAcid, you want to change is currently '%c' at position %d.\n",TheAminoAcid,pos_of_the_AA);
		draw_seq_aound_aa(stdout,FORMAT_SCREEN,pos_of_the_AA);
		printf("\n<< (L)eft (R)ight (S)elect (C)ancel (J)ump>>\n");
		charkey=get_char();
		fflush(stdin);
		switch(UPPER(charkey))
			{
			case('L'):pos_of_the_AA = Fct_Pos(INSERT,pos_of_the_AA-3);break;
			case('R'):pos_of_the_AA = Fct_Pos(INSERT,pos_of_the_AA+3);break;
			case('C'):pos_of_the_AA=0;TheAminoAcid=EOS;return(FALSE);break;
			case('S'):done=TRUE;INKEY;break;
			case('J'):
					{
					printf("Input the new position please /* %d */:",pos_of_the_AA);
					pos_of_the_AA=Get_Number(pos_of_the_AA);
					while(Fct_Frame(INSERT,pos_of_the_AA)!=IS_IN_FRAME)
						pos_of_the_AA = Fct_Pos(INSERT,pos_of_the_AA+1);
					TheAminoAcid=Translation_at(INSERT,pos_of_the_AA);
					}break;
			default:BEEP;break;
			}
		}
	
	printf("The AminoAcid you want to change is '%c' at position %d.\n",TheAminoAcid,pos_of_the_AA);

	return(TRUE);
	}

boolean DirectMutagenesis(FILE *out, short mode,int Pos)
	{
	char line[2][PALINDROME_MAX_ENZYME+PALINDROME_MAX_ENZYME+PALINDROME_MAX_ENZYME],c;
	int i,j,k,debut,find;
	int partial=0;
	/****************************************/
	if(npb[INSERT]<=0 || pos_of_the_AA<=0 || TheAminoAcid==EOS) {BEEP;return(FALSE);}
	if(mode!=FORMAT_SCREEN) printMainHeader(out,mode,TRUE);
	/********/
	if(Search_Done==FALSE) Cmd_Get_Site();
	if(nbr_sites<=0) fprintf(out,"No site was found !!!.\n");
	printf("Looking for a site at %d in order to change %c...",pos_of_the_AA,TheAminoAcid);
	draw_seq_aound_aa(out,mode,pos_of_the_AA);
	TotalSolutions=0;
	printf("This part of the program was a suggestion from Dr. S.LOPEZ (my manager !) from mexico\nSearching silent mutagenesis,please wait...\n");
	for(i=0;i<nbr_enzyme;i++)
		{
		
		k=0;
		if(Fct_N_sites(i,INSERT,pos_of_the_AA,pos_of_the_AA,TRUE)!=0)
			{
			TotalSolutions++;
			DRAW_HR(out,'-');
				fprintf(out,"\n%s has found a solution to mute %c at %d.\n",VAR_VERSION,TheAminoAcid,pos_of_the_AA);
			fprintf(out,"\nThere is already a  %s [%s] (%s) site at %d !.\n",Enzymes[i].nom,Enzymes[i].site_complet,Enzymes[i].site,pos_of_the_AA);
			partial=Fct_N_sites(i,INSERT,Fct_Pos(INSERT,pos_of_the_AA-1),Fct_Pos(INSERT,pos_of_the_AA+1),FALSE);
			if(partial !=0)
				fprintf(out,"BEWARE: there are also %d sites on INSERT.\n",partial);
			draw_seq_aound_aa(out,mode,pos_of_the_AA);
			printf("Searching silent mutagenesis, please wait...\n");
			continue;
			}
		else
		debut=pos_of_the_AA-Enzymes[i].taille_site;
		for(j=debut;j<pos_of_the_AA+Enzymes[i].taille_site+3;j++)
			{
			c=sequence[INSERT][Fct_Pos(INSERT,j)];
			line[0][k++]=c;line[0][k]=EOS;
			}
		for(j=debut+1;j<=pos_of_the_AA+2;j++)
			{
			strcpy(line[1],line[0]);
			

			memcpy( &line[1][j-(pos_of_the_AA-Enzymes[i].taille_site)],
					&Enzymes[i].site,
					Enzymes[i].taille_site*sizeof(char));
			if(strlen(line[0])!=strlen(line[1])) ERROR_USER;			
			find=TRUE;

			for(k=0;k<strlen(line[0]);k++)
				{
				if(Fct_Frame(INSERT,k+debut)!=IS_IN_FRAME) continue;
				if(k>=(strlen(line[0])-2)) break;
				if(IsAACodedBy(Fct_Traduction(line[0][k],line[0][k+1],line[0][k+2]),line[1][k],line[1][k+1],line[1][k+2])!=TRUE)
						{find=FALSE;break;}
				
				}
			if(find==TRUE)
				{
				TotalSolutions++;
				DRAW_HR(out,'-');
				fprintf(out,"\n%s has found a solution to mute %c at %d.\n",VAR_VERSION,TheAminoAcid,pos_of_the_AA);
				fprintf(out,"Use Enzyme %s [%s] (%s)\n",Enzymes[i].nom,Enzymes[i].site_complet,Enzymes[i].site);
				partial=Fct_N_sites(i,INSERT,1,npb[INSERT],TRUE);
				if(partial !=0)
					fprintf(out,"BEWARE: There is (are) already %d site(s) on INSERT.\n",partial);
				fprintf(out,"\nTranslated  ");
				for(k=0;k<strlen(line[0]);k++)
					{
					if(Fct_Frame(INSERT,k+debut)!=IS_IN_FRAME)
						fprintf(out," ");
					else
						fprintf(out," %c",Translation_at(INSERT,k+debut));
					}
				fprintf(out,"\nOriginal 5' ");
				for(k=0;k<strlen(line[0]);k++)
					{
					if(Fct_Frame(INSERT,k+debut)==IS_IN_FRAME) fprintf(out," ");
					fprintf(out,"%c",line[0][k]);
					}
				fprintf(out," 3'\nMutation 5' ");
				for(k=0;k<strlen(line[1]);k++)
					{
					if(Fct_Frame(INSERT,k+debut)==IS_IN_FRAME) fprintf(out," ");
					fprintf(out,"%c",line[1][k]);
					}
				fprintf(out," 3'\n");

				fprintf(out,"\n");
				printf("Searching silent mutagenesis...\n");
				}
			}
		}
	/********/
	if(TotalSolutions==0) printf("No solution was found !\n");
	if(mode!=FORMAT_SCREEN) printMainHeader(out,mode,FALSE);
	TotalSolutions=0;
	return(TRUE);
	}

boolean SystematicMutageneis(FILE *out, short mode)
	{
	int i;
	for(i=var_min[INSERT];i<=var_max[INSERT];i++)
		{
		if(Fct_Frame(INSERT,i)!=IS_IN_FRAME) continue;
		DRAW_HR(out,'-');
		fprintf(out,"Position = %d\n",i);
		pos_of_the_AA=i;
		TheAminoAcid=Translation_at(INSERT,pos_of_the_AA);
		draw_seq_aound_aa(out,mode,pos_of_the_AA);
		DirectMutagenesis(out,mode,pos_of_the_AA);
		}
	pos_of_the_AA=0;TheAminoAcid=EOS;
	return(TRUE);
	}





/* File 'MacFiles.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <string.h>
	#include "ADN.h"
#endif

#ifdef __MAC__
	#ifdef __MULTIFILE__
		#include <Files.h>
		#include <Dialogs.h>
	#endif

#ifdef __MAC__
boolean DoneDrawingPagesOfATextFile(FILE *someFile)
	{
	short	i;
	Str255	aStringOfText;

	TextSize(10);
	TextFont(4);
	for (i = 1; i <= 50; ++i)
		{
		if (feof(someFile))
				return TRUE;
		fgets((char *)aStringOfText, 255, someFile);
		CtoPstr((char *)aStringOfText);
		
		if (aStringOfText[aStringOfText[0]] == '\n')
			aStringOfText[aStringOfText[0]] = ' ';
		MoveTo(10, 14 * i);
		DrawString(aStringOfText);
		}
	return FALSE;
	}


/********************************  open a file to write  ***/
FILE *MacWriteTextFile(char *Nom_du_Fichier,Str255 originalName)
	{
	SFReply		theReply;
	SFTypeList	theTypeList = {'TEXT'};
	Point		where = {40, 60};
	OSErr		PotentialErr;
	FILE		*in;
	SFPutFile(where,"\pSave file as...",originalName,NULL,&theReply);
	if (theReply.good)
		{
		PotentialErr  = SetVol(NULL, theReply.vRefNum);
		p2cstr(theReply.fName);
		strcpy(Nom_du_Fichier,(char*)theReply.fName);
		in = fopen((char*)theReply.fName, "w");
		c2pstr((char *)theReply.fName);
		return(in);
		}
	return (NULL);
	}
#endif


char Macfgetc(short sRefFichier)
	{
	OSErr		oserror;
	long count=1;
	char c;
	char var_return=EOF;
	
	oserror = FSRead (sRefFichier,&count,&c);
	
	if((long)oserror==-39) return(EOF);/* end of file */
	 	
	if( count != 1 ||  c==EOF  || oserror!=0)
		{
		if(oserror!=0 && c!=EOF ) Handle_error((long)oserror);
		var_return=EOF;
		}
	else
		{
		if(c=='\r' || c==21 || c==8) c='\n';
		var_return=c;
		}
		return(var_return);
	}
/***********************************************************************/
boolean Cmd_LOAD_SEQUENCE(int NumSeq)
	{
	Str255	str255NomFichier;
	short	 sReferenceVolume;
	Search_Done=FALSE;
	if (fGetFileName(str255NomFichier, &sReferenceVolume,FICHIER_ADN[NumSeq])==TRUE)
		{
		OpenSelectedFile(str255NomFichier, sReferenceVolume,NumSeq);
		}
	return(npb[NumSeq]>0);
	} /* GetSequence() */
	
/***********************************************************************/
Boolean	fGetFileName	(Str255 str255Nom, short *psRefVolume, char *name)
	{
	SFTypeList	listeTypes;
	Point		pointPosition;
	SFReply		infoFichier;
	
	listeTypes[0]= 'xDNA';
	listeTypes[1]= 'TEXT';
	
	pointPosition.h = 100;
	pointPosition.v = 100;
	
	SFGetFile(pointPosition, "\p", 0L, 2, listeTypes, 0L, &infoFichier);
	if (infoFichier.good==FALSE)
		return(FALSE);
	p2cstr(infoFichier.fName);
	strcpy((char *)str255Nom, (const char *)infoFichier.fName);
	strcpy(name, (const char *)infoFichier.fName);
	c2pstr((char *)infoFichier.fName);
	*psRefVolume = infoFichier.vRefNum;
	
	return(TRUE);
	} /* fGetFileName() */

/***********************************************************************/
Boolean OpenSelectedFile(Str255 str255NomFichier, short sReferenceVolume,short NumSeq)
	{
	short				sReferenceFichier;
	Boolean var_return=FALSE;
	

	if ((var_return =fOuvreFichier(str255NomFichier, sReferenceVolume, &sReferenceFichier ))==TRUE)
		{
		_OuvreFichierSequence(NumSeq,sReferenceFichier);
		_FermeFichier(sReferenceFichier);
		}
	return(var_return);
	} /* OpenSelectedFile() */
/***********************************************************************/
Boolean	fOuvreFichier	(Str255 str255Nom, short sRefVolume, short *psRefFichier)
	{
	OSErr	oserror;

	c2pstr((char *)str255Nom);
	if((oserror = FSOpen(str255Nom, sRefVolume, psRefFichier))!=noErr)
		printf("impossible d'ouvrir ce fichier !\n");
	p2cstr(str255Nom);
	return (noErr == oserror);
	}
/***********************************************************************/
void	_FermeFichier	(short sRefFichier)
	{
		FSClose(sRefFichier);
	}
/***********************************************************************/
Boolean _OuvreFichierSequence(int NumSeq,short sRefFichier)
	{
	SFTypeList	listeTypes;
	Point		pointPosition;
	OSErr		oserror;
	char		c='\0';
	long		fileSize=0;
	Boolean		var_return=TRUE;
	
	listeTypes[0]   = 'TEXT';
	pointPosition.h = 100;
	pointPosition.v = 100;
	npb[NumSeq]=0;
	pos_ATG[NumSeq]=FALSE;
	var_min[NumSeq]=0;
	var_max[NumSeq]=0;
	degenerate[NumSeq]=FALSE;
	IsStriderSeq[NumSeq]=FALSE;
	
	if((oserror = GetEOF (sRefFichier, &fileSize))!=noErr)
		{
		Handle_error((long)oserror);
		return(FALSE);
		}
	if ((oserror =SetFPos(sRefFichier, fsFromStart, 0))!=noErr)
		{
		Handle_error((long)oserror);
		return(FALSE);
		}
	if (fileSize == 0)
		{
		printf("Le fichier est vide !\n");
		return(FALSE);
		}
	
	if(MacReadStriderFormat(NumSeq,sRefFichier)==FALSE)
		{
		if(SetFPos(sRefFichier, fsFromStart, 0)) return(FALSE);
		while((c=Macfgetc(sRefFichier))!=EOF)
			{
			if (c=='>'|| c==';')
					{
					while(Macfgetc(sRefFichier)!='\n')
						npb[NumSeq]=0;
					if(c==EOF) {printf("/* File error : %s */\n",FICHIER_ADN[NumSeq]);
								BEEP;
								Menu(19);
								break;}
					}
				
				if (c!='\n' && c!='/' && c!=' ' && c!='\t')
				if(c<48 || c>57)
					if(Fct_est_ADN(c)>=1)
						{
						npb[NumSeq]++;
						if(npb[NumSeq]<=MAX_NPB)
							{
							sequence[NumSeq][npb[NumSeq]]=UPPER(c);
							if(Fct_est_ADN(c)>=AMBIGOUS) degenerate[NumSeq]=TRUE;
							}
						else
							{
							printf("Sequence too large! [n>%d]\n\tOut of memory !\n",MAX_NPB);
							BEEP;Menu(19);}	
						}
				sequence[NumSeq][npb[NumSeq]+1]='\0';
				}

		printf( "\n\t\t\t%d bp in '%s'.\n\n",npb[NumSeq],FICHIER_ADN[NumSeq]);
			if(degenerate[NumSeq]==TRUE)
				{
				BEEP;
				Display("SORRY: This sequence contains degenerate bases",1);
				Menu(12);
				strcpy(FICHIER_ADN[NumSeq],"");
				npb[NumSeq]=0;
				}
				}
	if(npb[NumSeq]>0)
			{
			Cmd_POLYLINKER1(NumSeq);
			Get_ATG(NumSeq);
			}
	return(var_return);
	}

/***********************************************************************/
Boolean	fLitFichier	(short sRefFichier, TEHandle hteTexte)
	{
	long	lTailleTexte;
	OSErr   oserrResultat;
	
	TESetSelect(0, (*hteTexte)->teLength, hteTexte);
	TEDelete(hteTexte);

	if (noErr == GetEOF(sRefFichier, &lTailleTexte))
		{
		if (noErr == SetFPos(sRefFichier, fsFromStart, 0))
			{
			SetHandleSize((*hteTexte)->hText, lTailleTexte);
			(*hteTexte)->teLength = lTailleTexte;
			HLock((Handle)(*hteTexte)->hText);
			oserrResultat = FSRead(sRefFichier, &lTailleTexte, *(*hteTexte)->hText);
			HUnlock((Handle)(*hteTexte)->hText);
			}
		}
	return (noErr == oserrResultat);
	}/* fLitFichier() */


/***********************************************************************/
Boolean	fEcritFichier	(short sRefFichier, char *pcTexte, long lTaille)
{
	OSErr	oserrResultat;

	oserrResultat = SetFPos(sRefFichier, fsFromStart, 0);
	if (noErr == oserrResultat)
		oserrResultat = FSWrite(sRefFichier, &lTaille, pcTexte);
	if (noErr == oserrResultat)
		oserrResultat = SetEOF(sRefFichier, lTaille );
	return (noErr == oserrResultat);
} /* fEcritFichier() */
/***********************************************************************/
Boolean	fEnregistreFichier	(Str255 str255Nom, short sRefVolume, TEHandle hteTexte)
{
	short	  sRefFichier;
	Boolean	fResultat = FALSE;
	
	if (TRUE == fOuvreFichier(str255Nom, sRefVolume, &sRefFichier)) {
		HLock((Handle)(*hteTexte)->hText);
		fResultat = fEcritFichier(sRefFichier, *((*hteTexte)->hText),
								  (long)(*hteTexte)->teLength);
		HUnlock((Handle)(*hteTexte)->hText);
		_FermeFichier(sRefFichier);
	}
	return fResultat;
} /* fEnregistreFichier() */
/***********************************************************************/
Boolean	fEnregistreSousFichier (Str255 str255Nom, short *psRefVolume, TEHandle hteTexte)
{
	short	  sRefFichier;
	Boolean	fResultat = FALSE;

	if (fSelectionneNouveauFichier(str255Nom, psRefVolume)) {
		if (fCreeFichier(str255Nom, *psRefVolume, &sRefFichier)) {
			HLock((Handle)(*hteTexte)->hText);
			fResultat = fEcritFichier(sRefFichier, *((*hteTexte)->hText),
									  (long)(*hteTexte)->teLength);
			HUnlock((Handle)(*hteTexte)->hText);
			_FermeFichier(sRefFichier);
		}
	}
	return fResultat;
} /* fEnregistreSousFichier() */
/***********************************************************************/
Boolean	fSelectionneNouveauFichier	(Str255 str255Nom, short *psRefVolume)
{
	Str255	 str255Fichier;
	Point	  pointPosition;
	SFReply	infoFichier;

	strcpy((char *)str255Fichier, (const char *)str255Nom);
	pointPosition.h = 100;
	pointPosition.v = 100;

	c2pstr((char *)str255Fichier);
	SFPutFile(pointPosition, NULL, str255Fichier, 0L, &infoFichier);
	p2cstr(str255Fichier);
	if (!infoFichier.good) return FALSE;

	p2cstr(infoFichier.fName);
	strcpy((char *)str255Nom, (const char *)infoFichier.fName);
	c2pstr((char *)infoFichier.fName);
	*psRefVolume = infoFichier.vRefNum;
	
	return TRUE;
} /* fSelectionneNouveauFichier() */
/***********************************************************************/
Boolean	fCreeFichier	(Str255 str255Nom, short sRefVolume, short * psRefFichier)
{
	OSErr	oserrResultat;
	
	c2pstr((char *)str255Nom);
	oserrResultat = Create(str255Nom, sRefVolume, 'EDIT', 'TEXT');
	p2cstr(str255Nom);
	if ((noErr == oserrResultat) || (dupFNErr == oserrResultat))
		return fOuvreFichier(str255Nom, sRefVolume, psRefFichier);
	return FALSE;
} /* fCreeFichier() */
/***********************************************************************/
long	lCalculeTailleFichier	(Str255 str255Nom, short sRefVolume)
{
	long	lTaille = -1;
	short   sRefFichier;
	
	if (TRUE == fOuvreFichier(str255Nom, sRefVolume, &sRefFichier)) {
	   if (noErr != GetEOF(sRefFichier, &lTaille))
		  lTaille = -1;
		_FermeFichier(sRefFichier);
	}
	return lTaille;
} /* lCalculeTailleFichier() */
/***********************************************************************/

void  MacAlert(char *word)
	{
	Str255 word2;
	p2cstr(word2);
	strcpy((char *)word2,word);
	c2pstr((char *)word2);
	ParamText(word2,"\p","\p","\p");
	StopAlert(128, NULL);
	
	}
/***********************************************************************/
boolean MacReadStriderFormat(int NumSeq,short sRefFichier)
	{
	boolean r=TRUE,isStriderType=FALSE;
	STRIDER_HEADER signature;
	int i;
	char c;
	long count=1;
	long curseur;
	
	npb[NumSeq]=0;
	pos_ATG[NumSeq]=FALSE;
	var_min[NumSeq]=0;
	var_max[NumSeq]=0;
	IsStriderSeq[NumSeq]=FALSE;
	
	
	if((count=(long)sizeof(STRIDER_HEADER))!=112)
		{
		
		BEEP;
		fprintf(stderr,"SORRY, This version of %s is not able to analyse DNA Strider files:\n",VAR_VERSION);
		fprintf(stderr,"int=%d bytes should be 4 bytes.\n",(int)sizeof(StriderInt));
		fprintf(stderr,"short=%d bytes should be 2 bytes.\n",(int)sizeof(short));
		fprintf(stderr,"char=%d bytes should be 1 byte.\n",(int)sizeof(char));
		return(FALSE);
		}
	
	
	if(FSRead(sRefFichier,&count,&signature)==0)
		{
		/*
		printf("version=%d\n",(int)signature.versionNb);
		printf("type=%d\n",(int)signature.type);		
		printf("topology=%d\n",(int)signature.topology);
		printf("nLength=%d\n",(int)signature.nLength);
		printf("nMinus=%d\n",(int)signature.nMinus);
		printf("com_length=%d\n",(int)signature.com_length);
		printf("cOUNT=%d\n",(int)count);
		/*if(count!=(long)sizeof(STRIDER_HEADER)) return(FALSE);*/
		
		
		if(signature.versionNb!=0)
			r=FALSE;
		}
	else r=FALSE;
	
	
	if(r==TRUE)
		{
		for(i=0;i<(int)signature.nLength;i++)
			{
			if((c=Macfgetc(sRefFichier))==EOF)
				{r=FALSE;break;}
			if(Fct_est_ADN(c)!=TRUE)
				{r=FALSE;break;}
			npb[NumSeq]++;
			if(npb[NumSeq]<MAX_NPB)
				sequence[NumSeq][npb[NumSeq]]=UPPER(c);
			else
				{r=FALSE;break;}
			}
		}
	
	
	/* memorise comment length position */
	if(r==TRUE)
		{
		if(GetFPos(sRefFichier,&curseur)!=0)
			r=FALSE;
		}


	/* check comment length */
	if(r==TRUE)
		{
		for(i=0;i<signature.com_length;i++)
			{
			if((c=Macfgetc(sRefFichier))==EOF)
				{r=FALSE;break;}
			}
		}


	/* check it is file EOF */
	if(r==TRUE)
		{
		if((c=Macfgetc(sRefFichier))!=EOF)
			r=FALSE;
		}


	if(r==TRUE)
		{
		isStriderType=TRUE;
		printf("\nIdentified as a DNA Strider Sequence.\n");
		if(Preference.display_messages==TRUE && SetFPos(sRefFichier,fsFromStart,curseur)==0 && signature.com_length>0)
			{
			printf("= Included comment =============\n");
			while((c=Macfgetc(sRefFichier))!=EOF)
				{
				printf("%c",c);
				}
			printf("\n");
			Menu(12);
			}
		}
	
	
	if(r==TRUE)
		{
		if((int)signature.nLength>=MAX_NPB)
			{
			printf("Sorry this sequence is too large (%d>%d)..!\n",(int)signature.nLength,(int)MAX_NPB);
			r=FALSE;
			}
		if(signature.type!=1)
			{
			printf("Sorry this sequence is not a non-degenerate DNA sequence !\n");
			r=FALSE;
			}
		if(signature.topology!=1)
			{
			printf("Sorry this sequence is not a circular sequence !\n");
			r=FALSE;
			}
		/*if(signature.nMinus!=1)
			{
			printf("Sorry this sequence contains non authorized base (nMinus)\n");
			r=FALSE;
			}*/
		if(r==FALSE)
			{
			npb[NumSeq]=0;
			INKEY;
			}
		}
	IsStriderSeq[NumSeq]=isStriderType;

	return(isStriderType);
	}

#endif




/* File 'MacProgressBar.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "ADN.h"
#endif
#ifdef __MAC__

#ifdef __MULTIFILE__
	#include <console.h>
	#include <Dialogs.h>
	#include <QuickDraw.h>
	#include <Types.h>
	#include <Events.h>
	#include <ToolUtils.h>
	
#endif

#define DI_ABOUT	131

void MacAbout(void)
	{
	DialogPtr	pDialogue;
    GrafPtr		pPortCourant;
	short 		sTypeControl,sControlSelectionne=FALSE;
	Rect		rectControl;
	Handle		hText,hImage;
	Str255		Str255Text;
	
	
 
 	pDialogue = GetNewDialog(DI_ABOUT, NULL, (WindowPtr)-1);
 	
 	
	MoveWindow((WindowPtr)pDialogue,
		(short)((qd.screenBits.bounds.right-(pDialogue->portRect.right-pDialogue->portRect.left))/2.0),
		(short)((qd.screenBits.bounds.bottom-(pDialogue->portRect.bottom-pDialogue->portRect.top))/2.0),
		TRUE);

	
	GetPort(&pPortCourant);
    SetPort((GrafPtr) pDialogue);
    GetDialogItem(pDialogue,4, &sTypeControl,&hText, &rectControl);
    GetDialogItem(pDialogue,3, &sTypeControl,&hImage, &rectControl);
    /*BeginUpdate((WindowPtr)pDialogue);*/
    /*ShowDialogItem(pDialogue,3);*/
    p2cstr(Str255Text);
    sprintf((char*)Str255Text,"%s (compiled %s at %s)",VAR_VERSION,__DATE__,__TIME__);	
    c2pstr((char *)Str255Text);
    SetDialogItemText(hText,Str255Text);
    /*EndUpdate((WindowPtr)pDialogue);*/
	ShowWindow((WindowPtr)pDialogue);
	FlushEvents(everyEvent, 0);
	do {
		ModalDialog(NULL,&sControlSelectionne);
		
        
        
        /*SetDialogItemText(hText,Str255Text);
       
    	}while (sControlSelectionne!=3);*/
    	}while(!Button());
    DisposeDialog(pDialogue);
    SetPort(pPortCourant);
	}
	
void ShowProgressBar(short state,int count, int the_max,char *TheMessage)
	{
	static DialogPtr	barDlogPtr;
	static GrafPtr		pPortCourant;
	static Handle		iHndl;
	static Rect		box;
	static short	i, old_right,itemType;
	static float	Lepas;
	/*unsigned long		counter;*/
	static	int		MemoryState=0;
	static Pattern ThePat;
	
	Str255 word2;
	
	
	
	
	switch(state)
	{
	case(1):
		{
		if(MemoryState!=0) RETURN_NULL;
		MemoryState=1;
		barDlogPtr = GetNewDialog(128, nil, (WindowPtr)-1L);
		GetPort(&pPortCourant);
		SetPort(barDlogPtr);
		GetDialogItem(barDlogPtr,3,&itemType,&iHndl,&box);/* get window caracterisitic */
		FrameRect(&box);/* draw a box */
		old_right = box.right;/* droite */
		GetIndPattern(&ThePat,0,4);
		}break;
	case(2):
		{
		if(MemoryState!=1) RETURN_NULL;
		p2cstr(word2);
		strcpy((char *)word2,TheMessage);
		c2pstr((char *)word2);
		GetDialogItem(barDlogPtr,2,&itemType,&iHndl,&box);
		SetDialogItemText(iHndl, word2);
		GetDialogItem(barDlogPtr,3,&itemType,&iHndl,&box);
		}break;
	case(3):
		{
		if(MemoryState!=1) RETURN_NULL;
		Lepas = (old_right - box.left)/*/100.0*/; /* use a float to get good*/
		box.right = box.left + (short)(((float)count/(float)the_max)*Lepas);
		FillRect(&box,&ThePat);
		}break;
	case(4):
		{
		if(MemoryState!=1) RETURN_NULL;
		MemoryState=0;
		box.right = old_right;
		GetIndPattern(&ThePat,0,4);
		FillRect(&box,&ThePat);		/* when we get close, */
		/*Delay(20,&counter);*/	/* wait a 1/3 of a second*/
		GetIndPattern(&ThePat,0,20);
		FillRect(&box,&ThePat);	/* and get back to the original*/
		FrameRect(&box);
		DisposeDialog(barDlogPtr);
		SetPort(pPortCourant);
		}break;
	default:RETURN_NULL;break;
	}
}
#endif



/* File 'MacPrefs.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "ADN.h"
#endif
#ifdef __MAC__

#ifdef __MULTIFILE__
	#include <Dialogs.h>
	#include <QuickDraw.h>
	#include <Types.h>
	#include <ToolUtils.h>
#endif

#define ID_DIALOG_WINDOW	130
#define DI_DIPLAY_OK		1
#define DI_DIPLAY_CANCEL	2
#define DI_DIPLAY_MESSAGE   3
#define DI_MEMORY	  		4
#define DI_KLENOW  			5
#define DI_EDIT_PARTIAL		6
#define DI_TEXT_PARTIAL		7
#define DI_CONT				8
#define DI_OVERHANG			9
#define DI_BUFFER			10
#define DI_TEMP				11
#define DI_CIP				12
#define DI_5_PRIME          13
#define DI_3_PRIME          14
#define DI_ONLY_BLUNT		15

boolean MacDialogPref(void)
	{
	return(ExecuterDialogPref());
	}

boolean ExecuterDialogPref(void)
	{
	DialogPtr	pDialogue;
	GrafPtr		pPortCourant;
	short		sTypeControl,sControlSelectionne;
	Handle		hItemType[DI_ONLY_BLUNT+1];
    Rect		rectControl;
    Str255 		word2;
	STRUCT_PREFS savepref;
	memcpy(&savepref,&Preference,sizeof(STRUCT_PREFS));
	
    pDialogue = GetNewDialog(ID_DIALOG_WINDOW, NULL, (WindowPtr)-1);
	if(pDialogue==NULL) return(FALSE);

	GetPort(&pPortCourant);
    SetPort((GrafPtr) pDialogue);
    
    GetDialogItem(pDialogue, DI_DIPLAY_MESSAGE	, &sTypeControl,&hItemType[DI_DIPLAY_MESSAGE],	&rectControl);
    GetDialogItem(pDialogue, DI_MEMORY			, &sTypeControl,&hItemType[DI_MEMORY],			&rectControl);
    GetDialogItem(pDialogue, DI_KLENOW			, &sTypeControl,&hItemType[DI_KLENOW],			&rectControl);
    GetDialogItem(pDialogue, DI_TEXT_PARTIAL		, &sTypeControl,&hItemType[DI_TEXT_PARTIAL],	&rectControl);
    GetDialogItem(pDialogue, DI_EDIT_PARTIAL		, &sTypeControl,&hItemType[DI_EDIT_PARTIAL],	&rectControl);
    GetDialogItem(pDialogue, DI_ONLY_BLUNT			, &sTypeControl,&hItemType[DI_ONLY_BLUNT],	&rectControl);
    GetDialogItem(pDialogue, DI_CONT				, &sTypeControl,&hItemType[DI_CONT],			&rectControl);
    GetDialogItem(pDialogue, DI_OVERHANG			, &sTypeControl,&hItemType[DI_OVERHANG],		&rectControl);
    GetDialogItem(pDialogue, DI_BUFFER			, &sTypeControl,&hItemType[DI_BUFFER],			&rectControl);
    GetDialogItem(pDialogue, DI_TEMP				, &sTypeControl,&hItemType[DI_TEMP],			&rectControl);
    GetDialogItem(pDialogue, DI_CIP				, &sTypeControl,&hItemType[DI_CIP],				&rectControl);
    GetDialogItem(pDialogue, DI_5_PRIME			, &sTypeControl,&hItemType[DI_5_PRIME],			&rectControl);
    GetDialogItem(pDialogue, DI_3_PRIME			, &sTypeControl,&hItemType[DI_3_PRIME],			&rectControl);
    
    p2cstr(word2);
	sprintf((char*)word2,"%d",Preference.partial);
	c2pstr((char *)word2);
	SetDialogItemText(hItemType[DI_EDIT_PARTIAL],word2);
      do
       {
	    SetControlValue((ControlHandle)hItemType[DI_DIPLAY_MESSAGE],Preference.display_messages);
	    SetControlValue((ControlHandle)hItemType[DI_MEMORY],Preference.memory);
	    SetControlValue((ControlHandle)hItemType[DI_KLENOW],Preference.allow_T4);
		
		HiliteControl((ControlHandle)hItemType[DI_ONLY_BLUNT],(Preference.partial==0?255:0));
	    SetControlValue((ControlHandle)hItemType[DI_ONLY_BLUNT],Preference.partial_only_blunt);
	    SetControlValue((ControlHandle)hItemType[DI_CONT],Preference.allow_all_sol);
	    SetControlValue((ControlHandle)hItemType[DI_OVERHANG],Preference.allow_part_overhang);
	    SetControlValue((ControlHandle)hItemType[DI_BUFFER],Preference.buffer);
	    SetControlValue((ControlHandle)hItemType[DI_TEMP],Preference.temperature);
	    SetControlValue((ControlHandle)hItemType[DI_CIP],Preference.allow_CIP);
	    SetControlValue((ControlHandle)hItemType[DI_5_PRIME],Preference.side_5);
	    SetControlValue((ControlHandle)hItemType[DI_3_PRIME],Preference.side_3);
		
		ModalDialog(NULL, &sControlSelectionne);
		
		switch(sControlSelectionne)
			{
			case(DI_DIPLAY_MESSAGE):
				Preference.display_messages=(Preference.display_messages==TRUE?FALSE:TRUE);
				/*SetControlValue((ControlHandle)hItemType[DI_DIPLAY_MESSAGE],!(GetControlValue((ControlHandle)hItemType[DI_DIPLAY_MESSAGE])));*/
				break;
			case(DI_MEMORY):	Preference.memory=(Preference.memory==TRUE?FALSE:TRUE);break;
			case(DI_KLENOW):	Preference.allow_T4=(Preference.allow_T4==TRUE?FALSE:TRUE);break;
			case(DI_ONLY_BLUNT):Preference.partial_only_blunt=(Preference.partial_only_blunt==TRUE?FALSE:TRUE);break;
			case(DI_CONT):		Preference.allow_all_sol=(Preference.allow_all_sol==TRUE?FALSE:TRUE);break;
			case(DI_OVERHANG):	Preference.allow_part_overhang=(Preference.allow_part_overhang==TRUE?FALSE:TRUE);break;
			case(DI_TEMP):		Preference.temperature=(Preference.temperature==TRUE?FALSE:TRUE);break;
			case(DI_BUFFER):	Preference.buffer=(Preference.buffer==TRUE?FALSE:TRUE);break;
			case(DI_CIP):		Preference.allow_CIP=(Preference.allow_CIP==TRUE?FALSE:TRUE);break;
			case(DI_5_PRIME):	Preference.side_5=(Preference.side_5==TRUE?FALSE:TRUE);break;
			case(DI_3_PRIME):	Preference.side_3=(Preference.side_3==TRUE?FALSE:TRUE);break;
			case(DI_EDIT_PARTIAL):
				{
				GetDialogItemText(hItemType[DI_EDIT_PARTIAL],word2);
				p2cstr(word2);
				printf("¥%s¥\n",word2);
				if(strlen((char*)word2)>1 || (word2[0]!='1' && word2[0]!='2' && word2[0]!='3'))
					{
					sprintf((char*)word2,"");
					c2pstr((char *)word2);
					SetDialogItemText(hItemType[DI_EDIT_PARTIAL],word2);
					}
				else
					{
					Preference.partial=word2[0]-'0';
					c2pstr((char *)word2);
					}
				}
				break;
						
			default:break;
			}
		}while (sControlSelectionne!=1 && sControlSelectionne!=2);

		DisposeDialog(pDialogue);
		SetPort(pPortCourant);
    
    if(sControlSelectionne==1)
    	{
        return(TRUE);
    	}
    else
       {
		memcpy(&Preference,&savepref,sizeof(STRUCT_PREFS));
		return FALSE;
       }
      
 }
#endif




/* File 'MacBrowser.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "ADN.h"
#endif
#ifdef __MAC__

#ifdef __MULTIFILE__
	#include <Dialogs.h>
	#include <QuickDraw.h>
	#include <Types.h>
	#include <Events.h>
	#include <ToolUtils.h>
	
#endif

#define ID_BROWSER_WINDOW 132

char MacBrowser(short type, int index)
	{
	static DialogPtr	pDialogue;
	static GrafPtr		pPortCourant;
	short		sTypeControl,sControlSelectionne;
	static Handle hItemType[14];
    Rect		rectControl;
	short i=0, j;

	switch(i)
		{
		case(0):pDialogue = GetNewDialog(ID_BROWSER_WINDOW, NULL, (WindowPtr)-1);
				if(pDialogue==NULL) return('A');
				
				GetPort(&pPortCourant);
			    SetPort((GrafPtr) pDialogue);
			    for(j=0;j<14;j++) GetDialogItem(pDialogue,128+j, &sTypeControl,&hItemType[j],&rectControl);
	    		i=1;
	    		ShowWindow((WindowPtr)pDialogue);
	    		return('&');
	    		break;
	    case(1):if(i!=1) return('A');
				ModalDialog(NULL, &sControlSelectionne);
				switch(sControlSelectionne)
					{
					case(1):return('P');break;
					case(7):return('N');break;
					case(3):return('D');break;
					case(2):return('A');break;
					case(8):return('X');break;
					case(4):return('H');break;
					case(5):return('S');break;
					case(6):return('L');break;
					case(9):return('T');break;
					case(10):return('I');break;
					case(11):MacAbout();return('&');break;
					default:break;
					}
				return('&');
				break;
		case(2):if(i!=1) return('&');
				DisposeDialog(pDialogue);
				SetPort(pPortCourant);
				return('&');
				break;
		}
	}

#endif



/* File 'Strider.c' */
#ifndef __GNUC__
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "ADN.h"
#endif

#ifndef __MAC__

boolean ReadStriderFormat(int NumSeq)
	{
	FILE *in;
	boolean r=TRUE,isStriderType=FALSE;
	STRIDER_HEADER signature;
	fpos_t	curseur;
	int i;
	char c;
	
	npb[NumSeq]=0;
	pos_ATG[NumSeq]=FALSE;
	var_min[NumSeq]=0;
	var_max[NumSeq]=0;
	IsStriderSeq[NumSeq]=FALSE;
	
	if(sizeof(STRIDER_HEADER)!=112)
		{
		BEEP;
		fprintf(stderr,"SORRY, This version of %s is not able to analyse DNA Strider files:\n",VAR_VERSION);
		fprintf(stderr,"int (StriderInt) =%d bytes should be 4 bytes.\n",(int)sizeof(StriderInt));
		fprintf(stderr,"short=%d bytes should be 2 bytes.\n",(int)sizeof(short));
		fprintf(stderr,"char=%d bytes should be 1 byte.\n",(int)sizeof(char));
		return(FALSE);
		}
	if((in=fopen(FICHIER_ADN[NumSeq],"rb"))==NULL)
		return(FALSE);
	
	if(fread(&signature,(size_t)sizeof(STRIDER_HEADER),1,in)==1)
		{
		/*
		printf("version=%d\n",(int)signature.versionNb);
		printf("type=%d\n",(int)signature.type);		
		printf("topology=%d\n",(int)signature.topology);
		printf("nLength=%d\n",(int)signature.nLength);
		printf("nMinus=%d\n",(int)signature.nMinus);
		printf("com_length=%d\n",(int)signature.com_length);
		*/
		if(signature.versionNb!=0)
			r=FALSE;
		}
	else r=FALSE;
	
	if(r==TRUE)
		{
		for(i=0;i<(int)signature.nLength;i++)
			{
			if((c=fgetc(in))==EOF)
				{r=FALSE;break;}
			if(Fct_est_ADN(c)!=TRUE)
				{r=FALSE;break;}
			npb[NumSeq]++;
			if(npb[NumSeq]<MAX_NPB)
				sequence[NumSeq][npb[NumSeq]]=UPPER(c);
			else
				{r=FALSE;break;}
			}
		}
	/* memorise comment length position */
	if(r==TRUE)
		{
		if(fgetpos(in,&curseur)!=0) r=FALSE;
		}

	/* check comment length */
	if(r==TRUE)
		{
		for(i=0;i<signature.com_length;i++)
			{
			if((c=fgetc(in))==EOF)
				{r=FALSE;break;}
			}
		}

	/* check it is file EOF */
	if(r==TRUE)
		{
		if((c=fgetc(in))!=EOF)
			r=FALSE;
		}

	if(r==TRUE)
		{
		isStriderType=TRUE;
		printf("\nIdentified as a DNA Strider Sequence.\n");
		if(Preference.display_messages==TRUE && fsetpos(in,&curseur)==0 && signature.com_length>0)
			{
			printf("= Included comment =============\n");
			while((c=fgetc(in))!=EOF)
				{
				printf("%c",c);
				}
			printf("\n");
			Menu(12);
			}
		}
	
	if(r==TRUE)
		{
		if((int)signature.nLength>=MAX_NPB)
			{
			printf("Sorry this sequence is too large (%d>%d)..!\n",(int)signature.nLength,(int)MAX_NPB);
			r=FALSE;
			}
		if(signature.type!=1)
			{
			printf("Sorry this sequence is not a non-degenerate DNA sequence !\n");
			r=FALSE;
			}
		if(signature.topology!=1)
			{
			printf("Sorry this sequence is not a circular sequence !\n");
			r=FALSE;
			}
		/*if(signature.nMinus!=1)
			{
			printf("Sorry this sequence contains non authorized base (nMinus)\n");
			r=FALSE;
			}*/
		if(r==FALSE)
			{
			npb[NumSeq]=0;
			INKEY;
			}
		}
	IsStriderSeq[NumSeq]=isStriderType;
	fclose(in);
	return(isStriderType);
	}
#endif

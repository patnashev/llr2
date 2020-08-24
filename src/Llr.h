#ifndef LLRDOTH
#define LLRDOTH

/* Constants */

#define LLR_VERSION		"1.0.2"
#define LLR_VERSIONC	1,0,2

/* Global variables */

extern char INI_FILE[80];		/* Name of the prime INI file */

extern int volatile ERRCHK;			/* 1 to turn on error checking */
extern unsigned int PRIORITY;		/* Desired priority level */
extern unsigned int CPU_AFFINITY;	/* NT Processor affinity */

/* Factoring limits based on complex formulas given the speed of the */
/* factoring code vs. the speed of the Lucas-Lehmer code */
/* As an example, examine factoring to 2^68 (finding all 68-bit factors). */
/* First benchmark a machine to get LL iteration times and trial factoring */
/* times for a (16KB sieve of p=35000011). */
/*	We want to find when time spend eliminating an exponent with */
/* trial factoring equals time saved running 2 LL tests. */

/*	runs to find a factor (68) *
	#16KB sections (2^68-2^67)/p/(120/16)/(16*1024*8) *
	factoring_benchmark = 2.0 * LL test time (p * ll_benchmark)

	simplifying:

	68 * (2^68-2^67)/p/(120/16)/(16*1024*8) * facbench = 2 * p * llbench
	68 * 2^67 / p / (120/16) / 2^17 * facbench = 2 * p * lltime
	68 * 2^49 / p / (120/16) * facbench = p * lltime
	68 * 2^49 / (120/16) * facbench = p^2 * lltime
	68 * 2^53 / 120 * facbench = p^2 * lltime
	68 * 2^53 / 120 * facbench / lltime = p^2
	sqrt (68 * 2^53 / 120 * facbench / lltime) = p
*/

/* Now lets assume 30% of these factors would have been found by P-1.  So
   we only save a relatively quick P-1 test instead 2 LL tests.  Thus:
	sqrt (68 / 0.7 * 2^53 / 120 * facbench / lltime) = p
*/

/* Now factor in that 35000000 does 19 squarings, but 70000000 requires 20.
   Thus, if maxp is the maximum exponent that can be handled by an FFT size:
	sqrt (68 / 0.7 * 2^53 / 120 *
	      facbench * (1 + LOG2 (maxp/35000000) / 19) / lltime) = p
*/

/* Now factor in that errors sometimes force us to run more than 2 LL tests.
   Assume, 2.04 on average:
	sqrt (68 / 0.7 * 2^53 / 120 *
	      facbench * (1 + LOG2 (maxp/35000000) / 19) / lltime / 1.02) = p
*/

/* These breakeven points we're calculated on a 2.0 GHz P4 Northwood: */

#define FAC80	516000000L
#define FAC79	420400000L
#define FAC78	337400000L
#define FAC77	264600000L
#define FAC76	227300000L
#define FAC75	186400000L
#define FAC74	147500000L
#define FAC73	115300000L
#define FAC72	96830000L
#define FAC71	75670000L
#define FAC70	58520000L
#define FAC69	47450000L
#define FAC68	37800000L
#define FAC67	29690000L
#define FAC66	23390000L

/* These breakevens we're calculated a long time ago on unknown hardware: */

#define FAC65	13380000L
#define FAC64	8250000L
#define FAC63	6515000L
#define FAC62	5160000L
#define FAC61	3960000L
#define FAC60	2950000L
#define FAC59	2360000L
#define FAC58	1930000L
#define FAC57	1480000L
#define FAC56	1000000L

short default_work_type (void);
#define WORK_FACTOR		0
#define WORK_TEST		1
#define WORK_ADVANCEDTEST	2
#define WORK_DBLCHK		3
#define WORK_ECM		4
#define WORK_PMINUS1		5
#define WORK_PFACTOR		6
#define WORK_ADVANCEDFACTOR	7


extern unsigned long volatile ITER_OUTPUT;/* Iterations between outputs */
extern unsigned long volatile ITER_OUTPUT_RES;/* Iterations between results */
					/* file outputs */
extern unsigned long volatile DISK_WRITE_TIME;
					/* Number of minutes between writing */
					/* intermediate results to disk */
extern int TWO_BACKUP_FILES;		/* TRUE for 2 backup files(qXXXXXXX) */
extern int RUN_ON_BATTERY;		/* Run program even on battery power */
extern int TRAY_ICON;			/* Display tiny tray icon */
extern int HIDE_ICON;			/* Display no icon */
extern unsigned int PRECISION;		/* Number of decimal places to output*/
					/* in percent complete lines */
extern int CUMULATIVE_TIMING;		/* True if outputting cumulative time*/

/* Common routines */


void getCpuInfo (); 
 
int isPrime (unsigned long p);


#include "gwnum/gwini.h"

struct IniCache {						// defined in gwini.c but not in gwini.h
	char	*filename;
	int	immediate_writes;
	int	dirty;
	unsigned int num_lines;
	unsigned int array_size;
	struct IniLine **lines;
};

#ifdef __cplusplus
extern "C" {
#endif

struct IniCache *openIniFile (			// defined in gwini.c but not in gwini.h
	char	*,
	int);
void writeIniFile (						// defined in gwini.c but not in gwini.h
	struct IniCache *);
void IniFileOpen (char *, int);
void nameIniFiles (int);
void readIniFiles ();

#ifdef __cplusplus
}
#endif

void OutputBoth (char *);
void OutputSomewhere (char *);
void LogMsg (char *);
void ReplaceableLine (int);

unsigned long pick_fft_length (unsigned long);

void tempFileName (char	*, char, giant);
int fileExists (char *);
int readFileHeader (char *, int *, short *, unsigned long *);
int writeResults (char	*);

int communicateWithServer ();

#ifndef X86_64

EXTERNC void setupf();
EXTERNC int factor64();
EXTERNC void psetupf();
EXTERNC int pfactor64();

#endif

EXTERNC void* aligned_malloc (unsigned long, unsigned long);
EXTERNC void  aligned_free (void *);
EXTERNC void  truncated_strcpy (char*, unsigned long, const char*);


/* Routines called by common routines */

void OutputStr (char *);
int isHighResTimerAvailable (void); 
double getHighResTimer (void); 
double getHighResTimerFrequency (void); 
void guessCpuType ();
void getCpuDescription (char *, int); 
unsigned long num_cpus ();
#define stopCheck escapeCheck
int escapeCheck ();
#define	WORKING_ICON	0
#define	IDLE_ICON	1
void ChangeIcon (int);
void BlinkIcon (int);


/* Common routines */

int primeContinue ();

/* Routines used to time code chunks */

extern double __timers[10];		/* 10 timers are available */

/* Routines called by common routines */

void title (char *);
void flashWindowAndBeep ();
void SetPriority ();

#endif

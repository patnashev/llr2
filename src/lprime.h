/* Handy definitions */

#define FALSE	0
#define TRUE	1

/* This controls whether we want to pause computation if the load average */
/* becomes too great.  This does not apply to OS/2. */

#if defined (__linux__) || defined (__FreeBSD__) || defined (__APPLE__)
#define MPRIME_LOADAVG

#if !defined (__APPLE__)

/* Handle the difference between the naming conventions in the two */
/* C compilers.  We only need to to this for global variables that */
/* referenced by the assembly routines */

//#define CPU_FLAGS	_CPU_FLAGS 
#define FACHSW _FACHSW
#define FACMSW _FACMSW
#define FACLSW _FACLSW
#define FACPASS _FACPASS
#define FACTEST _FACTEST
#define SRCARG	_SRCARG
#define CPUID_EAX _CPUID_EAX
#ifndef X86_64
#define setupf _setupf
#define factor64 _factor64
#define psetupf _psetupf
#define pfactor64 _pfactor64
#define erdtsc _erdtsc
#endif

#endif

#define max(x, y) ((x) > (y))? x : y
#define min(x, y) ((x) < (y))? x : y
void Sleep (long);

/* Handle differences between Windows and Linux runtime libraries */

#define stricmp(x,y)	strcasecmp(x,y)
#define _commit(f)	fsync(f)
#define _open		open
#define _close		close
#define _read		read
#define _write		write
#define _lseek		lseek
#define _unlink		unlink
#define _creat		creat
#define _chdir		chdir
#define _ftime		ftime
#define _timeb		timeb
#define IsCharAlphaNumeric(c) isalnum(c)
#define _O_APPEND	O_APPEND
#define _O_RDONLY	O_RDONLY
#define _O_WRONLY	O_WRONLY
#define _O_RDWR		O_RDWR
#define _O_CREAT	O_CREAT
#define _O_TRUNC	O_TRUNC
#define _O_BINARY 	0
#define _O_TEXT		0

#endif

#define EXTERNC

/* The common include files */

#include <time.h>
#include <assert.h>
extern int NO_GUI;
/* #if defined (__linux__) || defined (__FreeBSD__) || defined (__APPLE__)
#include "gwrename.h"
#endif */
#include "./gwnum/cpuid.h"
#include "./gwnum/giants.h"
#include "./gwnum/gwnum.h"
#include "./gwnum/gwcommon.h"
#include "Llr.h"

/* Global variables */

extern int volatile THREAD_STOP;	/* TRUE if thread should stop */
extern int volatile THREAD_KILL;	/* TRUE if program should terminate */
extern int MENUING;			/* TRUE when main menu active */

/* Internal routines */

void main_menu ();
void linuxContinue (char *);


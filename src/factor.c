/**************************************************************
 *
 *	factor.c
 *
 *	General purpose factoring program
 *
 *	Updates:
 *		18 May 97   REC - invoked new, fast divide functions
 *		26 Apr 97	RDW - fixed tabs and unix EOL
 *		20 Apr 97	RDW - conversion to TC4.5
 *
 *	c. 1997 Perfectly Scientific, Inc.
 *	All Rights Reserved.
 *
 *
 *************************************************************/

/* include files */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#ifdef _WIN32 
#include <process.h>
#endif
#include <string.h>
#include "gwnum/giants.h"


/* definitions */

#define DD 100
#define NUM_PRIMES 6542 /* PrimePi[2^16]. */


/* compiler options */

#ifdef _WIN32
#pragma warning( disable : 4127 4706 ) /* disable conditional is constant warning */
#endif


/* global variables */

extern int volatile THREAD_STOP;

int		display = 0;
int 		pr[NUM_PRIMES];
giant 	xr = NULL, xs = NULL, zr = NULL, zs = NULL, xorg = NULL,
		zorg = NULL, t1 = NULL, t2 = NULL, t3 = NULL,
		gg = NULL, An = NULL, Ad = NULL;
giant 	xb[DD+2] = {NULL}, zb[DD+2] = {NULL}, xzb[DD+2] = {NULL};
int 	modcount = 0;
char	prline[80];


/* function prototypes */

void		ell_even(giant x1, giant z1, giant x2, giant z2, giant An,
						giant Ad, giant N);
void		ell_odd(giant x1, giant z1, giant x2, giant z2, giant xor,
						giant zor, giant N);
void		ell_mul(giant xx, giant zz, int n, giant An, giant Ad, giant N);
int 		least_2(int n);
void		dot(void);
int		psi_rand();


/**************************************************************
 *
 *	Functions
 *
 **************************************************************/


int			
psi_rand(
	void
)
{
	unsigned short	hi;
	unsigned short	low;
	time_t			tp;
	int				result;

	time(&tp);
	low = (unsigned short)rand();
	hi = (unsigned short)rand();
	result = ((hi << 16) | low) ^ ((int)tp);

	return (result & 0x7fffffff);
}


void
set_random_seed(
	void
)
{
	/* Start the random number generator at a new position. */
	time_t		tp;

	time(&tp);
	srand((int)tp + (int)getpid());
}


int
isprime(
	int 	odd
)
{
	int 	j;
	int 	p;

	for (j=1; ; j++)
	{
		p = pr[j];
		if (p*p > odd)
			return(1);
		if (odd % p == 0)
			return(0);
	}
}


int
primeq(
	int				p
)
{
	register int	j=3;

	if ((p&1)==0)
		return ((p==2)?1:0);
	if (j>=p)
		return (1);
	while ((p%j)!=0)
	{
		j+=2;
		if (j*j>p)
			return(1);
	}
	return(0);
}


void
ell_even(
	giant 	x1,
	giant 	z1,
	giant 	x2,
	giant 	z2,
	giant 	An,
	giant 	Ad,
	giant 	N
)
{
	gtog(x1, t1);
	addg(z1, t1);
	squareg(t1);
	modg(N, t1);
	gtog(x1, t2);
	subg(z1, t2);
	squareg(t2);
	modg(N, t2);
	gtog(t1, t3);
	subg(t2, t3);
	gtog(t2, x2);
	mulg(t1, x2);
	gshiftleft(2, x2);
	modg(N, x2);
	mulg(Ad, x2);
	modg(N, x2);
	mulg(Ad, t2);
	gshiftleft(2, t2);
	modg(N, t2);
	gtog(t3, t1);
	mulg(An, t1);
	modg(N, t1);
	addg(t1, t2);
	mulg(t3, t2);
	modg(N, t2);
	gtog(t2,z2);
}


void
ell_odd(
	giant 	x1,
	giant 	z1,
	giant 	x2,
	giant 	z2,
	giant 	xor,
	giant 	zor,
	giant 	N
)
{
	gtog(x1, t1);
	subg(z1, t1);
	gtog(x2, t2);
	addg(z2, t2);
	mulg(t1, t2);
	modg(N, t2);
	gtog(x1, t1);
	addg(z1, t1);
	gtog(x2, t3);
	subg(z2, t3);
	mulg(t3, t1);
	modg(N, t1);
	gtog(t2, x2);
	addg(t1, x2);
	squareg(x2);
	modg(N, x2);
	gtog(t2, z2);
	subg(t1, z2);
	squareg(z2);
	modg(N, z2);
	mulg(zor, x2);
	modg(N, x2);
	mulg(xor, z2);
	modg(N, z2);
}


void
ell_mul(
	giant 			xx,
	giant 			zz,
	int 			n,
	giant 			An,
	giant 			Ad,
	giant 			N
)
{
	unsigned int 	c = (unsigned int)0x80000000;

	if (n==1)
		return;
	if (n==2)
	{
		ell_even(xx, zz, xx, zz, An, Ad, N);
		return;
	}
	gtog(xx, xorg);
	gtog(zz, zorg);
	ell_even(xx, zz, xs, zs, An, Ad, N);

	while((c&n) == 0)
	{
		c >>= 1;
	}

	c>>=1;
	do
	{
		if (c&n)
		{
			ell_odd(xs, zs, xx, zz, xorg, zorg, N);
			ell_even(xs, zs, xs, zs, An, Ad, N);
		}
		else
		{
			ell_odd(xx, zz, xs, zs, xorg, zorg, N);
			ell_even(xx, zz, xx, zz, An, Ad, N);
		}
		c >>= 1;
	} while(c);
}



/* From R. P. Brent, priv. comm. 1996:
Let s > 5 be a pseudo-random seed (called $\sigma$ in the Tech. Report),

	u/v = (s^2 - 5)/(4s)

Then starting point is (x_1, y_1) where

	x_1 = (u/v)^3
and
	a = (v-u)^3(3u+v)/(4u^3 v) - 2
*/

void
choose12(
	giant 	x,
	giant 	z,
	int 	k,
	giant 	An,
	giant 	Ad,
	giant 	N
)
{
	itog(k, zs);
	gtog(zs, xs);
	squareg(xs);
	itog(5, t2);
	subg(t2, xs);
	modg(N, xs);
	addg(zs, zs);
	addg(zs, zs);
	modg(N, zs);
	gtog(xs, x);
	squareg(x);
	modg(N, x);
	mulg(xs, x);
	modg(N, x);
	gtog(zs, z);
	squareg(z);
	modg(N, z);
	mulg(zs, z);
	modg(N, z);

	/* Now for A. */
	gtog(zs, t2);
	subg(xs, t2);
	gtog(t2, t3);
	squareg(t2);
	modg(N, t2);
	mulg(t3, t2);
	modg(N, t2);  /* (v-u)^3. */
	gtog(xs, t3);
	addg(t3, t3);
	addg(xs, t3);
	addg(zs, t3);
	modg(N, t3);
	mulg(t3, t2);
	modg(N, t2);  /* (v-u)^3 (3u+v). */
	gtog(zs, t3);
	mulg(xs, t3);
	modg(N, t3);
	squareg(xs);
	modg(N, xs);
	mulg(xs, t3);
	modg(N, t3);
	addg(t3, t3);
	addg(t3, t3);
	modg(N, t3);
	gtog(t3, Ad);
	gtog(t2, An);  /* An/Ad is now A + 2. */
}


void
ensure(
	int	q
)
{
	int 	nsh, j;

/*	N = newgiant(INFINITY);
	if(!q)
	{
		gread(N,stdin);
		q = bitlen(N) + 1;
	}
*/
	nsh = q/4; /* Allowing (easily) enough space per giant,
					since N is about 2^q, which is q bits, or
				   q/16 shorts.  But squaring, etc. is allowed,
					so we need at least q/8, and we choose q/4
					to be conservative. */
	if (!xr)
		xr = newgiant(nsh);
	if (!zr)
		zr = newgiant(nsh);
	if (!xs)
		xs = newgiant(nsh);
	if (!zs)
		zs = newgiant(nsh);
	if (!xorg)
		xorg = newgiant(nsh);
	if (!zorg)
		zorg = newgiant(nsh);
	if (!t1)
		t1 = newgiant(nsh);
	if (!t2)
		t2 = newgiant(nsh);
	if (!t3)
		t3 = newgiant(nsh);
	if (!gg)
		gg = newgiant(nsh);
	if (!An)
		An = newgiant(nsh);
	if (!Ad)
		Ad = newgiant(nsh);
	for (j=0;j<DD+2;j++)
	{
		xb[j] = newgiant(nsh);
		zb[j] = newgiant(nsh);
		xzb[j] = newgiant(nsh);
	}
}

void	wclear () {
	int j;
	if (xr) {
		free (xr);
		xr = NULL;
	}
	if (zr) {
		free (zr);
		zr = NULL;
	}
	if (xs) {
		free (xs);
		xs = NULL;
	}
	if (zs) {
		free (zs);
		zs = NULL;
	}
	if (xorg) {
		free (xorg);
		xorg = NULL;
	}
	if (zorg) {
		free (zorg);
		zorg = NULL;
	}
	if (t1) {
		free (t1);
		t1 = NULL;
	}
	if (t2) {
		free (t2);
		t2 = NULL;
	}
	if (t3) {
		free (t3);
		t3 = NULL;
	}
	if (gg) {
		free (gg);
		gg = NULL;
	}
	if (An) {
		free (An);
		An = NULL;
	}
	if (Ad) {
		free (Ad);
		Ad = NULL;
	}
	for (j=0;j<DD+2;j++) {
		if (xb[j]) {
			free (xb[j]);
			xb[j] = NULL;
		}
		if (zb[j]) {
			free (zb[j]);
			zb[j] = NULL;
		}
		if (xzb[j]) {
			free (xzb[j]);
			xzb[j] = NULL;
		}
	}
}

int
bigprimeq(
	giant 	x
)
{
	itog(1, t1);
	gtog(x, t2);
	subg(t1, t2);
	itog(5, t1);
	powermodg(t1, t2, x);
	if (isone(t1))
		return(1);
	return(0);
}


void
dot(
	void
)
{
	if (display)
		OutputStr (".");
}

/**************************************************************
 *
 *	Main Function
 *
 **************************************************************/

void gfact (
	giant	N,
	giant DV,
	int randmode,
	int limitbits,
	int disp
)
{
	int 	j, k, C, nshorts, cnt, count,
			pass, npr;
	long	B;
	display = disp;
	pr[0] = 2;
	for (k=0, npr=1;; k++)
	{
		if (primeq(3+2*k))
		{
			pr[npr++] = 3+2*k;
			if (npr >= NUM_PRIMES)
				break;
		}
	}

	ensure (bitlen (N));		// Allocate memory

	if (randmode == 0)
	{
		if (display) {
			OutputStr("\n");
			OutputStr("Commencing Pollard rho...\n");
		}
		itog(1, gg);
		itog(3, t1); itog(3, t2);

		for (j=0; j < 15000; j++)
		{
			if (stopCheck ())
				goto STOPFACT;
			if((j%100) == 0)
			{
				dot();
				gcdg(N, gg);
				if (!isone(gg))
					break;
			}
			squareg(t1);
			iaddg(2, t1);
			modg(N, t1);
			squareg(t2);
			iaddg(2, t2);
			modg(N, t2);
			squareg(t2);
			iaddg(2, t2);
			modg(N, t2);
			gtog(t2, t3);
			subg(t1, t3);
			t3->sign = abs(t3->sign);
			mulg(t3, gg);
			modg(N, gg);
		}
		gcdg(N, gg);

		if ((gcompg(N,gg) != 0) && (!isone(gg)))
		{				// gg != N and gg != 1
			gtog(gg,DV);	// divisor candidate
			divg(gg, N);	// quotient candidate
			if (isone(N))	// factorization is complete (?)
			{
				wclear();
				if (display)
					OutputStr ("\n");
				return;
			}
			if (bigprimeq(N)) // accept only a PRP divisor.
			{
				if (gcompg (N, DV) < 0) {
					gtog (N, DV);
					gtog (gg, N);	// assure DV <= cofactor
				}
				wclear();
				if (display)
					OutputStr ("\n");
				return;
			}
		}

		if (display) {
			OutputStr ("\n");
			OutputStr ("Commencing Pollard (p-1)...\n");
		}
		itog(1, gg);
		itog(3, t1);
		for (j=0; j< NUM_PRIMES; j++)
		{
			if (stopCheck ())
				goto STOPFACT;
			cnt = (int)(8*log(2.0)/log(1.0*pr[j]));
			if (cnt < 2)
				cnt = 1;
			for (k=0; k< cnt; k++)
			{
				powermod(t1, pr[j], N);
			}
			itog(1, t2);
			subg(t1, t2);
			mulg(t2, gg);
			modg(N, gg);

			if (j % 100 == 0)
			{
				dot();
				gcdg(N, gg);
				if	(!isone(gg))
					break;
			}
		}
		gcdg(N, gg);
		if ((gcompg(N,gg) != 0) && (!isone(gg)))
		{				// gg != N and gg != 1
			gtog (gg, DV);	// divisor candidate
			divg(gg, N);	// quotient candidate
			if (isone(N))	// factorization is complete (?)
			{
				wclear();
				if (display)
					OutputStr ("\n");
				return;
			}
			if (bigprimeq(N)) // accept only a PRP divisor.
			{
				if (gcompg (N, DV) < 0) {	// assure DV <= cofactor
					gtog (N, DV);
					gtog (gg, N);
				}
				wclear();
				if (display)
					OutputStr ("\n");
				return;
			}
		}
	} /* This is the end of (randmode == 0) */
	if (display) {
		OutputStr ("\n");
		OutputStr ("Commencing ECM...\n");
	}

	if (randmode)
		set_random_seed();
	pass = 0;
	while (++pass)
	{
		if (randmode == 0)
		{
			if (pass <= 3)
			{
				B = 1000;
			}
			else if (pass <= 10)
			{
				B = 10000;
			}
			else if (pass <= 100)
			{
				B = 100000L;
			} else
			{
				B = 1000000L;
			}
		}
		else
		{
			B = 2000000L;
		}
		C = 50*((int)B);

		/* Next, choose curve with order divisible by 16 and choose
		 *	a point (xr/zr) on said curve.
		 */

		/* Order-div-12 case. 
		 * cnt = 8020345;   Brent's parameter for stage one discovery
		 * of 27-digit factor of F_13.
		 */

		cnt = psi_rand(); /* cnt = 8020345; */
		choose12(xr, zr, cnt, An, Ad, N);
		if (display) {
			sprintf(prline, "Choosing curve %d, with s = %d, B = %ld, C = %d:\n", pass, cnt, B, C);
			OutputStr(prline);
		}
		cnt = 0;
		nshorts = 1;
		count = 0;
		for (j=0;j<nshorts;j++)
		{
			ell_mul(xr, zr, 1<<16, An, Ad, N);
			ell_mul(xr, zr, 3*3*3*3*3*3*3*3*3*3*3, An, Ad, N);
			ell_mul(xr, zr, 5*5*5*5*5*5*5, An, Ad, N);
			ell_mul(xr, zr, 7*7*7*7*7*7, An, Ad, N);
			ell_mul(xr, zr, 11*11*11*11, An, Ad, N);
			ell_mul(xr, zr, 13*13*13*13, An, Ad, N);
			ell_mul(xr, zr, 17*17*17, An, Ad, N);
		}
		k = 19;
		while (k<B)
		{
			if (stopCheck ())
				goto STOPFACT;
			if (isprime(k))
			{
				ell_mul(xr, zr, k, An, Ad, N);
				if (k<100)
					ell_mul(xr, zr, k, An, Ad, N);
				if (cnt++ %100==0)
					dot();
			}
			k += 2;
		}
		count = 0;

		gtog(zr, gg);
		if (!isZero(gg))
			gcdg(N, gg);
		else
			gtog(N, gg);
		if ((!isone(gg))&&(bitlen(gg)>limitbits))
		{				// gg != 1 and not too small
			gtog(gg, DV);	// divisor candidate
			divg(gg, N);	// quotient candidate
			if (isone(N))	// factorization is complete (?)
			{
				wclear();
				if (display)
					OutputStr ("\n");
				return;
			}
			if (bigprimeq(N)) // accept only a PRP divisor.
			{
				if (gcompg (N, DV) < 0) {	// assure DV <= cofactor
					gtog (N, DV);
					gtog (gg, N);
				}
				wclear();
				if (display)
					OutputStr ("\n");
				return;
			}
			continue;		// Continue Stage 1
		}
		else
		{
			if (display)
				OutputStr("\n");	// End Stage 1
		}

		/* Continue;  Invoke, to test Stage 1 only. */
		k = ((int)B)/DD;
		gtog(xr, xb[0]);
		gtog(zr, zb[0]);
		ell_mul(xb[0], zb[0], k*DD+1 , An, Ad, N);
		gtog(xr, xb[DD+1]);
		gtog(zr, zb[DD+1]);
		ell_mul(xb[DD+1], zb[DD+1], (k+2)*DD+1 , An, Ad, N);

		for (j=1; j <= DD; j++)
		{
			gtog(xr, xb[j]);
			gtog(zr, zb[j]);
			ell_mul(xb[j], zb[j], 2*j , An, Ad, N);
			gtog(zb[j], xzb[j]);
			mulg(xb[j], xzb[j]);
			modg(N, xzb[j]);
		}
		modcount = 0;
		if (display) {
			sprintf(prline, "\nCommencing second stage, curve %d...\n",pass);
			OutputStr(prline);
		}
		count = 0;
		itog(1, gg);

		while (1)
		{
			gtog(zb[0], xzb[0]);
			mulg(xb[0], xzb[0]);
			modg(N, xzb[0]);
			mulg(zb[0], gg);
			modg(N,gg); /* Accumulate. */
			for (j = 1; j < DD; j++)
			{
				if (stopCheck ())
					goto STOPFACT;
				if (!isprime(k*DD+1+ 2*j))
					continue;

				/* Next, accumulate (xa - xb)(za + zb) - xa za + xb zb. */
				gtog(xb[0], t1);
				subg(xb[j], t1);
				gtog(zb[0], t2);
				addg(zb[j], t2);
				mulg(t1, t2);
				modg(N, t1);
				subg(xzb[0], t2);
				addg(xzb[j], t2);
				modg(N, t2);
				--modcount;
				mulg(t2, gg);
				modg(N, gg);
				if((++count)%1000==0)
					dot();
			}

			k += 2;
			if(k*DD > C)
				break;
			gtog(xb[DD+1], xs);
			gtog(zb[DD+1], zs);
			ell_odd(xb[DD], zb[DD], xb[DD+1], zb[DD+1], xb[0], zb[0], N);
			gtog(xs, xb[0]);
			gtog(zs, zb[0]);
		}

		if (!isZero(gg))
			gcdg(N, gg);
		else
			gtog(N, gg);
		if((!isone(gg))&&(bitlen(gg)>limitbits))
		{				// gg != 1 and not too small
			gtog(gg, DV);	// divisor candidate
			divg(gg, N);	// quotient candidate
			if (isone(N))	// factorization is complete (?)
			{
				wclear();
				if (display)
					OutputStr ("\n");
				return;
			}
			if (bigprimeq(N)) // accept only a PRP divisor.
			{
				if (gcompg (N, DV) < 0) {	// assure DV <= cofactor
					gtog (N, DV);
					gtog (gg, N);
				}
				wclear();
				if (display)
					OutputStr ("\n");
				return;
			}
			continue;
		}
		if (display)
			OutputStr("\n");
	}
STOPFACT:
	THREAD_STOP = 0;
	itog (1, DV);
	wclear();
	return;
}


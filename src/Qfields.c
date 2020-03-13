#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "./gwnum/giants.h"

#define  ULONG unsigned long
#define MAXULONG 0xFFFFFF00

void trace(int);

EXTERNC ULONG gcd (ULONG, ULONG);
extern unsigned long globalb;
extern double globalk;

/************************** Giants.c extensions *************************************/

#define BITSINULONG 32

char pbuf[256];

int gmodi (	/* Returns g%i as an integer */ 
	uint32_t den, giant g) { 
//	uint32_t wordweight, j, k, size, value;
	uint32_t size, value;
	uint32_t denval = den;
	giant gmod;
	giantstruct gdenstruct = {1, &denval};
	giant gden = &gdenstruct;
	int sign;
	if (den==1 || g->sign==0) return 0;
	if (g->sign < 0) {
	    sign = -1;
	    size = -g->sign;
	}
	else {
	    sign = 1;
	    size = g->sign;
	}
	gmod = newgiant (size*sizeof(uint32_t)+8);
	gtog (g, gmod);
	modg (gden, gmod);
	value = gmod->n[0];
/*	wordweight = 1;
	value = 0;
	for (j=0; j<size; j++) {
	    value += (uint32_t)((__int64)(g->n[j]%i)*wordweight)%i;
	    if (value >= i) value -= i;
	    for (k=1; k<=BITSINULONG; k++) {
			wordweight <<=1;
			if (wordweight >= i) wordweight -= i;
	    }
	} */
	free (gmod);
	return (sign*value);
}

void uldivg (uint32_t den, giant num) {
	uint32_t denval = den;
	giantstruct gdenstruct = {1, &denval};
	giant gden = &gdenstruct;
	divg (gden, num);
}

/************************************************************************************/

int ispower (unsigned long x, unsigned long b) {
//	Test if non negative integer x is a power of non negative integer b.
	if (b == 0)
		return ((x == 0)? TRUE : FALSE);	// Only zero can be a power of zero.
	if (x == 0)
		return (FALSE);						// No power of positive b can be zero.
	if (b == 1)
		return ((x == 1)? TRUE : FALSE);	// All positive powers of one are one.
	while (x%b == 0)						// General case.
		x /= b;
	return ((x == 1)? TRUE : FALSE);
}

void Reduce (uint32_t x, uint32_t *d, uint32_t *b) {

// Reduce a Discriminant to a square free integer.
// Given x, compute d, whithout square factor, and b, such as x = d*b^2

	uint32_t div, sq;

	*d = x;
	*b = 1;

	if (x<4)
		return;

	while (!((*d)%4)) {			// Divide by even power of two.
		*d /= 4;
		*b *= 2;
	}

	for (div = 3; (sq = div*div) <= *d; div += 2)
		while (!((*d)%sq)) {	// Divide by even powers of odd factors.
			*d /= sq;
			*b *= div;
		}
}

uint32_t issquare (uint32_t n) {
// This function returns the square root of an integer square, or zero.
	uint32_t s;
	s = (uint32_t)floor(sqrt((double) n));
	if (s*s == n)
		return s;
	else
		return 0;
}

uint32_t twopownmodm (uint32_t n, uint32_t m, uint32_t *order, uint32_t *nmodorder) {
	uint32_t tpnmodm, tp, i, work, mask = 1<<31;
	unsigned __int64 ltp;
	tpnmodm = 0;				// This function computes 2^n modulo m
	if (!(m&1)||(m==1)) {		// It returns this value, and also :
		*order = 0;				// The order of 2 modulo m, and the remainder of n modulo this order.
		*nmodorder = n;
	}
	if (m==1)
		return (tpnmodm);
	tp = 1;
	for (i=1; i<m; i++) {
	    tp <<= 1;
	    if (tp >= m) tp -= m;	// Modular reduction
	    if (i==n) tpnmodm = tp;	// If m is even or n < order of 2 modulo m the calculus is completed here.
	    if (tp==1) {			// But continue to compute the order of 2
			*order = i;
			break;
	    }
	}
	if (*order)
		*nmodorder = n%*order;	// Compute the remainder
	if (!tpnmodm) {				// If n >= order, continue
		work = *nmodorder;
		if (!work)				// n is a multiple of the order...
			return (1);
		while (!(work&mask))
			mask >>= 1;
		ltp = 1;				// init the result
		while (mask) {
			ltp *= ltp;			// square the result
			ltp %= m;			// modular reduction
			if (work&mask) {	// test the current bit of the exponent
				ltp <<= 1;		// multiply the result by the base
				if (ltp >= m) ltp -= m;
			}
			mask >>= 1;			// shift the mask
		}

/*	    for (i=1, tp=1; i<=*nmodorder; i++) {
			tp <<= 1;
			if (tp >= m) tp -= m;
	    } */
	    tpnmodm = (uint32_t)ltp;	// 2^n modulo m == 2^(n%order) modulo m
	}
	return tpnmodm;
}

uint32_t Bachet(uint32_t u, uint32_t v, long *a, long *b) {

//  Computes a and b such as a*u+b*v = gcd(u,v),
//  returns gcd(u,v).

	uint32_t n=0, m11=1, m12=0, m21=0, m22=1;
	uint32_t q, newm11, newm21, newu;

	while (v!=0) {
		q = u/v;
		newm11 = m11*q+m12;
		m12 = m11;
		m11 = newm11;
		newm21 = m21*q+m22;
		m22 = m21;
		m21 = newm21;
		newu = v;
		v = u-q*v;
		u = newu;
		n++;
	}
	if (n&1) {
		*a = -(int)m22;
		*b = (int)m12;
	}
	else {
		*a = (int)m22;
		*b = -(int)m12;
	}
	return u;
}

int gen_v1(giant k, uint32_t n, int general, int eps2, int debug) {
	long sign, jNd, jNa, v;
	uint32_t kmod3, rawd, d, dred, kmodd, tpnmd, i, orderd, Nmodd;
	uint32_t X, Y, aplus, aminus, b, rplus, rminus;
	uint32_t nmodorderd, ared, kmoda, tpnma, ordera, nmodordera, Nmoda;


	kmod3 = gmodi (3, k);		// kmod3 == k modulo 3

	if (kmod3 && !general) {	// Consider only the simple case !
		if ((kmod3 == 1 && !(n&1)) || (kmod3 == 2 && (n&1))) {
			if (debug) {		// 1*2^(2m) = 2*2^(2m+1) = 1 (modulo 3), so, N = 0 (modulo 3)
				sprintf (pbuf,"d = 3 divides N !\n");
				OutputBoth (pbuf);
			}
			return (-3);
		}
		else {					// 1*2^(2m+1) = 2*2^(2m) = 2 (modulo 3), so, N = 1 (modulo 3)
			if (debug) {		// Display the result if required.
				sprintf (pbuf,"epsilon = 2+sqrt(3)\n");
				OutputBoth (pbuf);
				sprintf (pbuf, "k = %d (mod 3), n = %d (mod 2)\n", kmod3, n&1);
				OutputBoth (pbuf);
				sprintf (pbuf, "v1 = 4, d = 3, a = 6, b = 2, r = 24, +1\n");
				OutputBoth (pbuf);
				sprintf (pbuf, "v1 = 4, d = 3, a = 2, b = 2, r = 8, -1\n");
				OutputBoth (pbuf);
			}					// Jacobi(3,N) = -Jacobi(N,3) = -1, if n>=3, and Jacobi(2,N) = 1 with minus sign.
			return (4);			// The conditions for Riesel theorem 5 are satisfied!
		}
	}											// End simple case

	if (eps2) {
		for (v=1; (rawd = v*v+4)<MAXULONG; v++) {	// General case ; searching units with norm == -1
			Reduce (rawd, &d, &b);			// v^2+4 == d*b^2 with d square free.
			dred = (d&1)? d : d>>1;			// dred == odd part of d.
			sign = ((n>2)||(d&1))? 1 : -1;	// Jacobi(d,N) == Jacobi(dred,N) if n>2 or d is odd, else they are opposite.
			kmodd = gmodi (dred, k);		// Now, all modular reductions are done modulo odd part of d.
			if (!kmodd) continue;			// N == -1 (mod dred) ==> Jacobi(dred,N) == 1 ==> d is not valid for n>2...
			if (n>1 && (((dred-1)/2) & 1))	// Jacobi(dred,N) == Jacobi(N,dred)*(-1)^((N-1)/2)*((dred-1)/2)
				sign = - sign;
			tpnmd = twopownmodm (n, dred, &orderd, &nmodorderd);// tpnmd == 2^n modulo dred.
			Nmodd = (kmodd*tpnmd-1+dred)%dred;	// Nmodd = N modulo dred ; be careful to avoid unsigned overflow...
			if (!Nmodd && (dred != 1)) {		// dred divides N!
				if (debug) {
					if (d&1) {
						sprintf (pbuf, "d = %d divides N !\n", dred);
						OutputBoth (pbuf);
					}
					else {
						sprintf (pbuf, "d/2 = %d divides N !\n", dred);
						OutputBoth (pbuf);
					}
				}
				return (-(int)dred);		// Return the divisor with minus sign.
			}
			if ((jNd = jacobi(Nmodd,dred)) > 1) {
				if (debug) {
					sprintf (pbuf, "%ld divides d = %u and N !\n", jNd, d);
					OutputBoth (pbuf);
				}
				return (-jNd);				// Return the divisor with minus sign.
			}
			if ((sign*jNd) != -1) continue;	// This value of d cannot be used.
			if (debug) {					// OK, we have found the fundamental unit. 
				if (v&1 || b&1) {			// Display this unit if required.
					if (b != 1) {
						sprintf (pbuf, "epsilon = [%ld+%u*sqrt(%u)]/2\n", v, b, d);
						OutputBoth (pbuf);
					}
					else {
						sprintf (pbuf, "epsilon = [%ld+sqrt(%u)]/2\n", v, d);
						OutputBoth (pbuf);
					}
				}
				else {
					if (b/2 != 1) {
						sprintf (pbuf, "epsilon = %ld+%u*sqrt(%u)\n", v/2, b/2, d);
						OutputBoth (pbuf);
					}
					else {
						sprintf (pbuf, "epsilon = %ld+sqrt(%u)\n", v/2, d);
						OutputBoth (pbuf);
					}
				}
			}
			b = b*v;					// Compute the square of the unit.
			v = v*v+2;					// This is the solution!
			aplus = v+2;
			aminus = v-2;				// aminus is a square, so, it is valid.
			rplus = 4*aplus;
			rminus = 4*aminus;
			if (debug) {				// Display the result if required.
				sprintf (pbuf, "k = %u (mod %u), n = %u (mod %u)\n", kmodd, dred, nmodorderd, orderd);
				OutputBoth (pbuf);
				sprintf (pbuf, "v1 = %ld, d = %u, a = %u, b = %u, r = %u, +1,eps2\n", v, d, aplus, b, rplus);
				OutputBoth (pbuf);
				sprintf (pbuf, "v1 = %ld, d = %u, a = %u, b = %u, r = %u, -1,eps2\n", v, d, aminus, b, rminus);
				OutputBoth (pbuf);
			}
			return v;
		}								// End for (v=1; (rawd = v*v+4)<MAXULONG; v++)
    return -1;							// Unable to find a value for v...
	}									// End if (eps2).

	for (v=3; (rawd = v*v-4)<MAXULONG; v++) {	// General case
		Reduce (rawd, &d, &b);			// v^2-4 == d*b^2 with d square free.
		dred = (d&1)? d : d>>1;			// dred == odd part of d.
		sign = ((n>2)||(d&1))? 1 : -1;	// Jacobi(d,N) == Jacobi(dred,N) if n>2 or d is odd, else they are opposite.
		kmodd = gmodi (dred, k);		// Now, all modular reductions are done modulo odd part of d.
		if (!kmodd) continue;			// N == -1 (mod dred) ==> Jacobi(dred,N) == 1 ==> d is not valid for n>2...
		if (n>1 && (((dred-1)/2) & 1))	// Jacobi(dred,N) == Jacobi(N,dred)*(-1)^((N-1)/2)*((dred-1)/2)
			sign = - sign;
		tpnmd = twopownmodm (n, dred, &orderd, &nmodorderd);// tpnmd == 2^n modulo dred.
		Nmodd = (kmodd*tpnmd-1+dred)%dred;	// Nmodd = N modulo dred ; be careful to avoid unsigned overflow...
		if (!Nmodd && (dred != 1)) {		// dred divides N!
			if (debug) {
				if (d&1) {
					sprintf (pbuf, "d = %d divides N !\n", dred);
					OutputBoth (pbuf);
				}
				else {
					sprintf (pbuf, "d/2 = %d divides N !\n", dred);
					OutputBoth (pbuf);
				}
			}
			return (-(int)dred);		// Return the divisor with minus sign.
		}
		if ((jNd = jacobi(Nmodd,dred)) > 1) {
			if (debug) {
				sprintf (pbuf, "%ld divides d = %u and N !\n", jNd, d);
				OutputBoth (pbuf);
			}
			return (-jNd);				// Return the divisor with minus sign.
		}
		if ((sign*jNd) != -1) continue;	// This value of d cannot be used.
		aplus = v+2;
		aminus = v-2;
		rplus = 4*aplus;
		rminus = 4*aminus;

//		Search if the quadratic unit candidate is the square of the fundamental one.

		if ((X=issquare(v-2))) {
			Y=issquare((v+2)/d);		// Yes,the fundamental unit has norm == -1
			if (debug) {				// And then, the candidate is already valid,
				if (X&1 || Y&1) {		// because Jacobi(v-2,N) is +1 and sign is -1
					if (Y != 1) {
						sprintf (pbuf, "epsilon = [%u+%u*sqrt(%u)]/2\n", X, Y, d);
						OutputBoth (pbuf);
					}
					else {
						sprintf (pbuf, "epsilon = [%u+sqrt(%u)]/2\n", X, d);
						OutputBoth (pbuf);
					}
				}
				else {
					if (Y/2 != 1) {
						sprintf (pbuf, "epsilon = %u+%u*sqrt(%d)\n", X/2, Y/2, d);
						OutputBoth (pbuf);
					}
					else {
						sprintf (pbuf, "epsilon = %u+sqrt(%u)\n", X/2, d);
						OutputBoth (pbuf);
					}
				}
				sprintf (pbuf, "k = %u (mod %u), n = %u (mod %u)\n", kmodd, dred, nmodorderd, orderd);
				OutputBoth (pbuf);
				sprintf (pbuf, "v1 = %ld, d = %u, a = %u, b = %u, r = %u, +1,eps2\n", v, d, aplus, b, rplus);
				OutputBoth (pbuf);
				sprintf (pbuf, "v1 = %ld, d = %u, a = %u, b = %u, r = %u, -1,eps2\n", v, d, aminus, b, rminus);
				OutputBoth (pbuf);
			}
			return v;
		}									// End v-2 is a square.

		else {								// No,the candidate is a fundamental unit of norm == +1
			for (ared=aminus, i=0; !(ared & 1); ared >>= 1) i++;// ared is the odd part of a.
			sign = -1;						// aminus^2-d*b^2 == (v-2)^2-d*b^2 == -4*aminus is negative.
			sign = ((n>2)||(aminus&1)||!(i&1))? sign : -sign;	// Jacobi(a,N) == Jacobi(ared,N) if n>2,
			kmoda = gmodi (ared, k);							// or a is odd, or i is even,else they are opposite.
			if (n>1 && (((ared-1)/2) & 1))	// Jacobi(ared,N) == Jacobi(N,ared)*(-1)^((N-1)/2)*((ared-1)/2)
				sign = - sign;
			tpnma = twopownmodm (n, ared, &ordera, &nmodordera);// tpnma == 2^n modulo ared.
			Nmoda = (kmoda*tpnma-1+ared)%ared;	// Nmoda == N modulo ared ; avoid unsigned overflow...
			if (!Nmoda && (ared != 1)) {		// ared divides N!
				if (debug) {
					if (aminus != ared)
						sprintf (pbuf, "a/%u = %u divides N !\n", aminus/ared, ared);
					else
						sprintf (pbuf, "a = %u divides N !\n", ared);
					OutputBoth (pbuf);
				}
				return (-(int)ared);		// Return the divisor with minus sign.
			}
			if ((jNa = jacobi(Nmoda,ared)) > 1) {
				if (debug) {
					sprintf (pbuf, "%ld divides a = %d and N !\n", jNa, aminus);
					OutputBoth (pbuf);
				}
				return (-jNa);				// Return the divisor with minus sign.
			}
			if ((sign*jNa) != -1) continue;	// Candidate not valid...
			if (debug) {					// OK, display the result if required.
				if (v&1 || b&1) {
					if (b != 1) {
						sprintf (pbuf, "epsilon = [%ld+%u*sqrt(%u)]/2\n", v, b, d);
						OutputBoth (pbuf);
					}
					else {
						sprintf (pbuf, "epsilon = [%ld+sqrt(%u)]/2\n", v, d);
						OutputBoth (pbuf);
					}
				}
				else {
					if (b/2 != 1) {
						sprintf (pbuf, "epsilon = %ld+%u*sqrt(%u)\n", v/2, b/2, d);
						OutputBoth (pbuf);
					}
					else {
						sprintf (pbuf, "epsilon = %ld+sqrt(%u)\n", v/2, d);
						OutputBoth (pbuf);
					}
				}
				sprintf (pbuf, "k = %u (mod %u), n = %u (mod %u) and n = %u (mod %u)\n", kmodd, dred, nmodorderd, orderd, nmodordera, ordera);
				OutputBoth (pbuf);
				sprintf (pbuf, "v1 = %ld, d = %u, a = %u, b = %u, r = %u, +1\n", v, d, aplus, b, rplus);
				OutputBoth (pbuf);
				sprintf (pbuf, "v1 = %ld, d = %u, a = %u, b = %u, r = %u, -1\n",v, d, aminus, b, rminus);
				OutputBoth (pbuf);
			}
			return v;
		}			// End v-2 is not a square.
    }				// End for (v=3; (rawd = v*v-4)<MAXULONG; v++)
    return -1;		// Unable to find a value for v...
}

int genProthBase(giant k, uint32_t n) {
	uint32_t Nmodp, kmodp, p, tpnmp, orderp, nmodorderp, kw;
	int jNp;

//	Return the least prime p such Jacobi (N, p) = -1

	if (k->sign == 1 && n < 3) {	//	Eliminate some trivial cases
		kw = k->n[0];
		if (n == 1 && kw == 1)
			return (2);
		else if (n == 2)
			return (2);
		else
			return (-1);
	}
	else {							// General case
		for (p = 3; p<=2147483647; p+=2) {
			if (!isPrime(p))
				continue;
			kmodp = gmodi (p, k);
			if (!kmodp)
				continue;
			tpnmp = twopownmodm (n, p, &orderp, &nmodorderp);
			Nmodp = (kmodp*tpnmp+1)%p;
			if (!Nmodp) {
				return (-(int)p);
			}
			if ((jNp = jacobi(Nmodp, p)) > 1) {
				return (-jNp);
			}
			if (jNp != -1)
				continue;
			return (p);
		}
		return (-1);
	}
}

int genProthBase1(giant N)
{
	uint32_t NmodD, D, dred, Nmod8;
	int jNp, chgsign;

//	Return the least D such as kronecker (D, N) = -1

//	Rem : for n>1, (N-1)/2 = k*2^(n-1) is even, so (D / N) = (N / D)

//  Then, if D = 2^s*t, (N / D) = (N / 2)^s * (N / t) = (N / 2)^s * (Nmodt / t)

	Nmod8 = N->n[0] & 7;

	for (D = 2; D<=2147483647; D++) {
		dred = D;
		chgsign = 1;
		while (!(dred&1)) {
			dred >>= 1;					// Compute the odd part of D
			if (Nmod8 == 3 || Nmod8 == 5)
				chgsign = -chgsign;
		}
		if (dred == 1)
			jNp = 1;
		else {
			NmodD = gmodi (dred, N);
			if (!NmodD)
				return (-(int)dred);
			if ((jNp = kronecker(NmodD, dred)) > 1)
				return (-jNp);
		}
		if ((jNp*chgsign) != -1)
			continue;
		return ((int)D);
		}
		return (-1);
}

int genLucasBaseQ(giant N, uint32_t D0) {
	uint32_t D, NmodD, dred, Nmod8;
	int jNp, chgsign;

//	Return the least D = 1+4*Q such as kronecker (D, N) = -1

//  if D = 2^s*t, (N / D) = (N / 2)^s * (N / t) = (N / 2)^s * (Nmodt / t)

	Nmod8 = N->n[0] & 7;

	for (D=D0; D<=2147483647; D+=4) {
		dred = D;
		chgsign = 1;
		while (!(dred&1)) {
			dred >>= 1;					// Compute the odd part of D
			if (Nmod8 == 3 || Nmod8 == 5)
				chgsign = -chgsign;
		}
		if (dred == 1)
			jNp = 1;
		else {
			NmodD = gmodi (dred, N);
			if (!NmodD)
				return (-(int)dred);
			if ((jNp = kronecker(NmodD, dred)) > 1)
				return (-jNp);
		}
		ulsubg (1, N);					// Compute N-1
		if (((dred-1) & 2) && (N->n[0] & 2))	// Quadratic reciprocity
			chgsign = -chgsign;
		iaddg (1, N);					// Restore N
		if ((jNp*chgsign) != -1)
			continue;
		return ((int)D);
	}
	return (-1);
}


int isLucasBaseQ(giant N, uint32_t D, int sign) {
	uint32_t NmodD, dred, Nmod8;
	int jNp, chgsign;

//	Return TRUE if D = 1+4*Q is such as kronecker (D, N) = sign

//  if D = 2^s*t, (N / D) = (N / 2)^s * (N / t) = (N / 2)^s * (Nmodt / t)

	Nmod8 = N->n[0] & 7;

//	for (D; D<=2147483647; D+=4) {
		dred = D;
		chgsign = 1;
		while (!(dred&1)) {
			dred >>= 1;					// Compute the odd part of D
			if (Nmod8 == 3 || Nmod8 == 5)
				chgsign = -chgsign;
		}
		if (dred == 1)
			jNp = 1;
		else {
			NmodD = gmodi (dred, N);
			if (!NmodD)
				return (-(int)dred);
			if ((jNp = kronecker(NmodD, dred)) > 1)
				return (-jNp);
		}
		iaddg (-1, N);					// Compute N-1
		if (((dred-1) & 2) && (N->n[0] & 2))	// Quadratic reciprocity
			chgsign = -chgsign;
		iaddg (1, N);					// Restore N
		if ((jNp*chgsign) == sign)
			return (TRUE);
		else
			return (FALSE);
//	}
}

int genLucasBaseP(giant N, uint32_t P0) {
	uint32_t P, NmodD, D, dred, Nmod8;
	int jNp, chgsign;

//	Return the least P such as D = P^2-4 and kronecker (D, N) = -1

//  if D = 2^s*t, (N / D) = (N / 2)^s * (N / t) = (N / 2)^s * (Nmodt / t)

	Nmod8 = N->n[0] & 7;

	for (P=P0; P*P<=2147483647; P++) {
		D = P*P-4;
		dred = D;
		chgsign = 1;
		while (!(dred&1)) {
			dred >>= 1;					// Compute the odd part of D
			if (Nmod8 == 3 || Nmod8 == 5)
				chgsign = -chgsign;
		}
		if (dred == 1)
			jNp = 1;
		else {
			NmodD = gmodi (dred, N);
			if (!NmodD)
				return (-(int)dred);
			if ((jNp = kronecker(NmodD, dred)) > 1)
				return (-jNp);
		}
		ulsubg (1, N);					// Compute N-1
		if (((dred-1) & 2) && (N->n[0] & 2))	// Quadratic reciprocity
			chgsign = -chgsign;
		iaddg (1, N);					// Restore N
		if ((jNp*chgsign) != -1)
			continue;
		return ((int)P);
	}
	return (-1);
}

long generalLucasBase(giant N, uint32_t *P0, uint32_t *Q) {
	uint32_t *P, NmodD, D, dred, Nmod8, NmodPQD, gcdNPQD;
	int jNp, chgsign;

//	Return the least D = P^2-4*Q such as kronecker (D, N) = -1

//  if D = 2^s*t, (N / D) = (N / 2)^s * (N / t) = (N / 2)^s * (Nmodt / t)

	Nmod8 = N->n[0] & 7;

	for (P=P0; (*P)*(*P)<=2147483647; (*P)++) {
		for ((*Q)=2; 4*(*Q)<(*P)*(*P); (*Q)++) {
			D = (*P)*(*P)-4*(*Q);
			if ((uint32_t)(floor(sqrt ((double)D)) * floor(sqrt ((double)D))) == D) {
				continue;
			}
			dred = D;
			chgsign = 1;
			while (!(dred&1)) {
				dred >>= 1;					// Compute the odd part of D
				if (Nmod8 == 3 || Nmod8 == 5)
					chgsign = -chgsign;
			}
			if (dred == 1)
				jNp = 1;
			else {
				NmodD = gmodi (dred, N);
				if (!NmodD)
					return (-(long)dred);
				if ((jNp = kronecker(NmodD, dred)) > 1)
					return (-jNp);
			}
			ulsubg (1, N);					// Compute N-1
			if (((dred-1) & 2) && (N->n[0] & 2))	// Quadratic reciprocity
				chgsign = -chgsign;
			iaddg (1, N);					// Restore N
			if (((jNp*chgsign) != -1) || ((globalk == 1.0) && ispower((*Q), globalb)))
				continue;
			NmodPQD = gmodi ((*P)*(*Q)*D, N);
			gcdNPQD = gcd (NmodPQD, (*P)*(*Q)*D);
			if (gcdNPQD != 1)
				return (-(long)gcdNPQD);
			return (D);
		}
	}
	return (-1);
}

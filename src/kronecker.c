#include "gwnum/gwcommon.h"

int kronecker(uint32_t a, uint32_t b) {
	/* Computes Kronecker (a, b) = (a/b) */
	/* Note : returns gcd (a, b) if a and b are not coprime. */

uint32_t  jdvs, jdvd, jq, jr, jmul;
int resul; uint32_t s,t,u,v,amod8;
jdvs = a;
jdvd = b;
jmul = 1;
amod8 = a & 7;
resul=1;
while (!(jdvd&1)) {			/* if b is even... */
	jdvd >>= 1;
	if ((amod8 == 3) || (amod8 == 5))
		resul = -resul;
	else if (!(jdvs & 1)) {
		jdvs >>= 1;
		jmul <<= 1;
	}
}		/* Now, jdvd is odd and jmul is the largest power of 2 dividing both a and b. */
if (jdvd == 1 || jdvs == 1)
	return((jmul == 1) ? resul : jmul);	/* Finished ! */
while(jdvs)
	if(jdvs==1)				/* Finished ! */
		return((jmul == 1) ? resul : jmul);
	else {
		v=jdvd;
		s=(v-1)>>1;  		/* (dvd-1)/2 */
		t=(v+1)>>1;  		/* (dvd+1)/2 */
		while(!(jdvs & 1)) {	/* While dvs is even */
						/* resul *= (-1)**(dvd**2-1)/8; */
			if (t & 1) {
				if ((s>>1) & 1)
					resul = -resul;
			}
                  else
				if ((t>>1) & 1)
					resul = -resul;
			jdvs >>= 1;		/* dvs /= 2; */
		}
		if(jdvs==1)			/* Finished ! */
			return((jmul == 1) ? resul : jmul);
		else {
			u=(jdvs-1)>>1; 	/* (dvs-1)/2 */
			if (s & u & 1)	/* resul *= (-1)**(dvd-1)*(dvs-1)/4; */
				resul = -resul;
			jq = jdvd/jdvs;
			jr = jdvd%jdvs;
			jdvd = jdvs;
			jdvs = jr;
		}      			/* dvs != 1 */
	}					/* dvs != 1 */
return(jdvd*jmul);			/* a and b are not coprime, */
						/* so, return their gcd. */
}

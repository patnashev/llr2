int jacobi(unsigned long a, unsigned long b) {
     /* Computes Jacobi (a, b) */

unsigned long  jdvs, jdvd, jq, jr;
int resul; unsigned long s,t,u,v;
jdvs = a;
jdvd = b;
resul=1;
while(jdvs)
	if(jdvs==1)				/* Finished ! */
		return(resul);
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
			return(resul);
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
return(jdvd);				/* a and b are not coprime, */
						/* so, return their gcd. */
}

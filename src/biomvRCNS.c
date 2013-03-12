#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Constants.h>
#include <stdlib.h>

#define m_e(x, i, j, nrow) x[(long)(i)+(long)(j)*(long)(nrow)]

int min(int a,int b) {
  if(a<b) return(a);
  else return(b);	
}

void print_matrix_double(double* x, int nrow, int ncol) {
    int i, j;
    for(i=0; i<nrow; i++) {
	Rprintf("%2d: ", i);
	for(j=0; j<ncol; j++) 
	    Rprintf("%7g ", m_e(x, i, j, nrow));
        Rprintf("\n");
    }
}

void print_matrix_int(int* x, int nrow, int ncol) {
    int i, j;
    for(i=0; i<nrow; i++) {
	Rprintf("%2d: ", i);
	for(j=0; j<ncol; j++) 
	    Rprintf("%d ", m_e(x, i, j, nrow));
        Rprintf("\n");
    }
}

/*dynamic programming algorithm of Picard et al. */
/*some codes borrowed from tilingArray*/

void univaRseg(double *C, int *maxcpr, int *maxkr, int *nr, double *mincost, int *minpos, double *mc, int *ss){
	/*  C variables */
	int maxcp = maxcpr[0];
	int maxk = maxkr[0];
	int n = nr[0];	
	double z, zmin;  
	int i, imin, cp, j, k, k0;  

	// start of dp 
	for(k=0; k<maxk; k++)	mincost[k] = C[k];
	for(k=maxk; k<n; k++)	mincost[k] = R_PosInf;

	for (cp=1; cp<maxcp; cp++) {	
		R_CheckUserInterrupt();
		for (j=0; j<n; j++) {   
			zmin = R_PosInf;
			imin = j;
			k0 = (j<maxk) ? j : maxk;
			for (k=0; k<k0; k++) { 
				z = m_e(mincost, j-k-1, cp-1, n);
				if (R_FINITE(z))	z += m_e(C, k, j-k, maxk);
				if(z<zmin) {
					zmin = z;
					imin = j-k;
				} 
			}  
			m_e(mincost, j, cp, n) = zmin; 
			m_e(minpos, j, cp-1, n) = imin; 
		} 
	} 

	for(cp=0;  cp<maxcp; cp++) {
		/* the minimal cost to be returned */
		mc[cp] = m_e(mincost, n-1, cp, n);
		/* Backtrack to get segs, i is always the changepoint to the right */
		for(j=cp+1; j<maxcp; j++)	m_e(ss, cp, j, maxcp) = -1; 
		m_e(ss, cp, cp, maxcp) = i = n; 
		for(j=cp-1; j>=0; j--) {
			m_e(ss, cp, j, maxcp) = i = m_e(minpos, i-1, j, n);
		}
	}
	/* add 1 to indices for R*/
	for(cp=0; cp<maxcp*maxcp; cp++){
		ss[cp] += 1; 
	}
	// end of dp
}




/* viterbi, forward and backward algorithm described in Guedon 2003 */
void forward(double *a, double *pi, double *b, double *d, double *D, int T, int DL, int J, int maxk, double *F, double *N, double *si) {

/* setup local var*/
	double obs;
	int u,i,j,t;
	
/*start*/
	for(t=0;t<T;t++) {
		R_CheckUserInterrupt();
		N[t]=0;
		for(j=0;j<J;j++) {
			m_e(F, t, j, T)=0; 
			obs=m_e(b, t, j, T);
			if(t<T-1) {
				for(u=1;u<=min(t+1,maxk);u++) {
					if(u<t+1) {
						if(DL>maxk){
							m_e(F, t, j, T)+=obs*m_e(d, t*maxk+u-1, j, DL)*m_e(si, t-u+1, j, T); 
							N[t]+=obs*m_e(D, t*maxk+u-1, j, DL)*m_e(si, t-u+1, j, T);
						} else {
							m_e(F, t, j, T)+=obs*m_e(d, u-1, j, maxk)*m_e(si, t-u+1, j, T); 
							N[t]+=obs*m_e(D, u-1, j, maxk)*m_e(si, t-u+1, j, T);
						}
/*						Rprintf("j= %d,	t= %d, u= %d, obs= %.3g, Nt-u= %.3g, m_e(si, t-u+1, j, T)= %.3g,\n",j, t, u, obs,N[t-u], m_e(si, t-u+1, j, T));*/
						obs*=m_e(b, t-u, j, T)/N[t-u];
/*						Rprintf("j= %d,	t= %d, u= %d, obs= %.3g, Nt= %.3g,\n",j, t, u, obs,N[t]);*/
					} else {
						if(DL>maxk){
							m_e(F, t, j, T)+=obs*m_e(d, t*maxk+t, j, DL)*pi[j];	
							N[t]+=obs*m_e(D, t*maxk+t, j, DL)*pi[j]; 
						} else {
							m_e(F, t, j, T)+=obs*m_e(d, t, j, maxk)*pi[j];	
							N[t]+=obs*m_e(D, t, j, maxk)*pi[j]; 
						}
/*						Rprintf("j= %d,	t= %d, u= %d, obs= %.3g, Nt= %.3g,\n",j, t, u, obs,N[t]);*/
					}
				}
			} else {
				for(u=1;u<=min(t+1,maxk);u++) {
					if(u<T) {
						if (DL>maxk){
							m_e(F, T-1, j, T)+=obs*m_e(D, t*maxk+u-1, j, DL)*m_e(si, T-u, j, T); 
						} else {
							m_e(F, T-1, j, T)+=obs*m_e(D, u-1, j, maxk)*m_e(si, T-u, j, T);  
						}
/*						Rprintf("j= %d,	t= %d, u= %d, obs= %.3g, NT-1-u= %.3g, m_e(si, T-u, j, T)\n",j, t, u, obs,N[T-1-u], m_e(si, T-u, j, T));*/
						obs*=m_e(b, T-1-u, j, T)/N[T-1-u];
/*						Rprintf("j= %d,	t= %d, u= %d, obs= %.3g, Nt= %.3g,\n",j, t, u, obs,N[t]);*/
					} else{
						if(DL>maxk){
							m_e(F, T-1, j, T)+=obs*m_e(D, t*maxk+T-1, j, DL)*pi[j];					 
						} else {
							m_e(F, T-1, j, T)+=obs*m_e(D, T-1, j, maxk)*pi[j];							
						}
/*						Rprintf("j= %d,	t= %d, u= %d, obs= %.3g, Nt= %.3g,\n",j, t, u, obs,N[t]);*/
					}
				}
				N[T-1]+=m_e(F, T-1, j, T); 
			}
		}

		for(j=0;j<J;j++){
			m_e(F, t, j, T) /= N[t];
			m_e(F, t, j, T)+= DOUBLE_XMIN;
/*			Rprintf("j= %d,	t= %d, u= %d, m_e(F, t, j, T)= %.3g,\n",j, t, u, m_e(F, t, j, T));*/
		}
		
		if(t<T-1) {
			for(j=0;j<J;j++){
				m_e(si, t+1, j, T) =0; 
				for(i=0;i<J;i++) {
					m_e(si, t+1, j, T)+=m_e(a, i, j, J)*m_e(F, t, i, T);
				}
/*				Rprintf("j= %d,	t= %d, u= %d, m_e(si, t+1, j, T)= %.3g,\n",j, t, u, m_e(si, t+1, j, T)); */
			}
		}
		
	}
	
/*	return nothing*/
}

/*the R .C interface for backwrd-forward algorithm*/
void backward(double *a, double *pi, double *b, double *d, double *D, int *maxk_r, int *DL_r, int *T_r, int *J_r, 
			double *eta, double *L, double *N,double *ahat, double *pihat, double *F, double *G, double *L1,double *si){

/*	checking inputs*/
	int maxk = maxk_r[0];	
	int J = J_r[0];
	int T = T_r[0];
	int DL = DL_r[0];

/*setup local variables */
	double obs, *den, *num;
	int i, j, t, u;
	den  = (double *)malloc(sizeof(double)*J);
	num  = (double *)malloc(sizeof(double)*J*J); 
	
/* all pointer used in forward should have been setup already*/
	forward(a, pi, b, d, D, T, DL, J, maxk, F, N, si);

/*start*/
	for(j=0;j<J;j++) {
		m_e(L, T-1, j, T)=m_e(F, T-1, j, T) ;
		for(u=0;u<DL;u++){
			m_e(eta, u, j, DL)=0;
		}
	}
	
	for(t=T-2;t>=0;t--) {
		 R_CheckUserInterrupt();
		for(j=0;j<J;j++) {
			m_e(G, t+1, j, T) =0;
			obs=1;
			for(u=1;u<=min(T-1-t, maxk);u++) {		
				obs*=m_e(b, t+u, j, T)/N[t+u]; //fixme, there is some larger variation in N-> G-> L1 -> L when bad soj is supplied
				
				if(u<T-1-t) {
					if(DL>maxk){
/*						Rprintf("j=%d, t=%d, u=%d, m_e(G, t+1, j, T)=%.3g, m_e(L1, t+u, j, T)=%.3g, m_e(F, t+u, j, T)=%.3g, obs=%.3g, m_e(d, t*maxk+u-1, j, DL)=%.3g\n", j, t, u, m_e(G, t+1, j, T), m_e(L1, t+u, j, T), m_e(F, t+u, j, T), obs, m_e(d, t*maxk+u-1, j, DL));*/
						m_e(G, t+1, j, T)+=m_e(L1, t+u, j, T)/m_e(F, t+u, j, T)*obs*m_e(d, t*maxk+u-1, j, DL);
						m_e(eta, t*maxk+u-1, j, DL)+=m_e(L1, t+u, j, T)/m_e(F, t+u, j, T)*obs*m_e(d, t*maxk+u-1, j, DL)*m_e(si, t+1, j, T);
					} else {
/*						Rprintf("j=%d, t=%d, u=%d, m_e(G, t+1, j, T)=%.3g, m_e(L1, t+u, j, T)=%.3g, m_e(F, t+u, j, T)=%.3g, obs=%.3g, m_e(d, u-1, j, maxk)=%.3g\n", j, t, u, m_e(G, t+1, j, T), m_e(L1, t+u, j, T), m_e(F, t+u, j, T), obs, m_e(d, u-1, j, maxk));*/
						m_e(G, t+1, j, T)+=m_e(L1, t+u, j, T)/m_e(F, t+u, j, T)*obs*m_e(d, u-1, j, maxk);
						m_e(eta, u-1, j, maxk)+=m_e(L1, t+u, j, T)/m_e(F, t+u, j, T)*obs*m_e(d, u-1, j, maxk)*m_e(si, t+1, j, T);
					}
				} else {
					if(DL>maxk){
						m_e(G, t+1, j, T)+=obs*m_e(D, t*maxk+T-1-t-1, j, DL);
						m_e(eta, t*maxk+u-1, j, DL)+=obs*m_e(D, t*maxk+u-1, j, DL)*m_e(si, t+1, j, T);
					} else {
						m_e(G, t+1, j, T) += obs*m_e(D, T-1-t-1, j, maxk);
						m_e(eta, u-1, j, maxk) += obs*m_e(d, u-1, j, maxk)*m_e(si, t+1, j, T);
					}
				}
				if(t==0){
					if(u>T-1) {
						if(DL>maxk){
							m_e(eta, t*maxk+u-1, j, DL)+=m_e(L1, t+u, j, T)/m_e(F, t+u, j, T)*obs*m_e(d, t*maxk+u-1, j, DL)*pi[j];
						} else {
							m_e(eta, u-1, j, maxk)+=m_e(L1, t+u, j, T)/m_e(F, t+u, j, T)*obs*m_e(d, u-1, j, maxk)*pi[j];
						}
					} else {
						if(DL>maxk){
							m_e(eta, t*maxk+u-1, j, DL)+= obs*m_e(d, t*maxk+u-1, j, DL)*pi[j];
						} else {
							m_e(eta, u-1, j, maxk)+= obs*m_e(d, u-1, j, maxk)*pi[j];
						}
					}
				}	
			} // end for u
		} // end for j	
		for(j=0;j<J;j++){
			m_e(L1, t, j, T) =0;
			for(i=0;i<J;i++) {
				m_e(L1, t, j, T) +=m_e(G, t+1, i, T)*m_e(a, j, i, J) ; //3.15
			}
			m_e(L1, t, j, T) *=m_e(F, t, j, T);//3.15
			m_e(L, t, j, T)=m_e(L1, t, j, T)+m_e(L, t+1, j, T)-m_e(G, t+1, j, T)*m_e(si, t+1, j, T); //3.12
		} // end for j
	} // end for t
	
	
/*	reset initial pi*/
	for(i=0;i<J;i++) {
		pihat[i]=0;
		den[i]=0; 
		for(j=0;j<J;j++) {
			m_e(num, j, i, J)=0;
			m_e(ahat, j, i, J)=0; 
		}
	}
/*	new estimates for a and pi*/
	for(i=0;i<J;i++){
		pihat[i]+=m_e(L, 0, i, T);					
		for(t=0;t<T-2;t++) {
			den[i]+=m_e(L1, t, i, T); 
			for(j=0;j<J;j++){
				m_e(num, j, i, J)+=m_e(G, t+1, j, T)*m_e(a, i, j, J)*m_e(F, t, i, T) ;			
			}		
		}	
	}

	for(i=0;i<J;i++) {
		for(j=0;j<J;j++) {
			m_e(ahat, i, j, J) =m_e(num, j, i, J)/den[i];			
		}
	}
	
	free(den);
	free(num);
}

void viterbi(double *a, double *pi, double *b, double *d, double *D, int *maxk_r, int *DL_r, int *T_r, int *J_r, 
	double *alpha, int *shat, double *si, int*opt, int* ops) {	 

/*	checking inputs*/
	int maxk = maxk_r[0];	
	int J = J_r[0];
	int T = T_r[0];
	int DL = DL_r[0];
	
/*setup local variables */	     
	double obs, prod, maxprod=-DOUBLE_XMAX;
	int i, j, t, u;//, *ops, *opt;
	
	for(t=0;t<T;t++) {
		R_CheckUserInterrupt();
		for(j=0;j<J;j++) {								
			obs=1; 
			if(t<T-1) {					
				for(u=1;u<=min(t+1,maxk);u++) {
					if(u<t+1) {
						if(DL>maxk){ 
							prod = obs*m_e(d, t*maxk+u-1, j, DL)*m_e(si, t-u+1, j, T);
/*							N[t]+=obs*m_e(D, t*maxk+u-1, j, DL)*m_e(si, t-u+1, j, T);*/
						} else {
							prod = obs*m_e(d, u-1, j, maxk)*m_e(si, t-u+1, j, T);
/*							N[t]+=obs*m_e(D, u-1, j, maxk)*m_e(si, t-u+1, j, T);*/
						}
						if(u==1 || maxprod < prod) {
							maxprod = prod;
							m_e(opt, t, j, T)=u;
						}	
						obs*=m_e(b, t-u, j, T);///N[t-u];
					}	else {
						if(DL>maxk){ 
							prod = obs*m_e(d, t*maxk+t, j, DL)*pi[j];
/*							N[t]+=obs*m_e(D, t*maxk+t, j, DL)*pi[j]; */
						} else {
							prod = obs*m_e(d, t, j, maxk)*pi[j];
/*							N[t]+=obs*m_e(D, t, j, maxk)*pi[j]; */
						}
						if(u==1 || maxprod < prod) {
							maxprod = prod;
							m_e(opt, t, j, T)=u;
						}	
					}

				}			
			} else {
				for(u=1;u<=min(t+1,maxk);u++) {
					if(u<T) {
						if (DL>maxk){
							prod =obs*m_e(D, t*maxk+u-1, j, DL)*m_e(si, T-u, j, T); 
						} else {
							prod =obs*m_e(D, u-1, j, maxk)*m_e(si, T-u, j, T);  
						}					
						obs *= m_e(b, T-1-u, j, T);///N[T-1-u];
					} else {
						if(DL>maxk){
							prod=obs*m_e(D, t*maxk+T-1, j, DL)*pi[j];					 
						} else {
							prod=obs*m_e(D, T-1, j, maxk)*pi[j];							
						}
					}
					if(u==1 || maxprod < prod) {
						maxprod = prod;
						m_e(opt, t, j, T)=u;
					}						
				}
			}
			m_e(alpha, t, j, T) = maxprod*m_e(b, t, j, T);	

		}

		if(t<T-1) {
			for(j=0;j<J;j++){
				i=0;
				m_e(si, t+1, j, T)=m_e(a, i, j, J)*m_e(alpha, t, i, T);
				m_e(ops, t+1, j, T)=0;
				for(i=1;i<J;i++){
					if(i!=j) {
						prod = m_e(a, i, j, J)*m_e(alpha, t, i, T);
						if(m_e(si, t+1, j, T) <= prod) {
							m_e(si, t+1, j, T) = prod;
							m_e(ops, t+1, j, T)=i;
						}
					}
				} 
			}
		}
	}

	shat[T-1] = 0;
	for(j=1;j<J;j++) {
		if(m_e(alpha, T-1, shat[T-1] , T) < m_e(alpha, T-1, j, T)) shat[T-1] = j;
	}
	u=1;
	for(t=T-2;t>=0;t--) {
		if(u < m_e(opt, t+u, shat[t+u], T)) {
			shat[t] = shat[t+u];
			u++;
		} else {
			shat[t] = m_e(ops, t+u, shat[t+u], T);
			u=1;
		}
	}

}
void logviterbi(double *a, double *pi, double *b, double *d, double *D, int *maxk_r, int *DL_r, int *T_r, int *J_r, 
	double *alpha, int *shat, double *si, int*opt, int* ops) {	 

/*	checking inputs*/
	int maxk = maxk_r[0];	
	int J = J_r[0];
	int T = T_r[0];
	int DL = DL_r[0];
	
/*setup local variables */	     
	double obs, prod, maxprod=-DOUBLE_XMAX;
	int i, j, t, u;//, *ops, *opt;

	for(t=0;t<T;t++) {
		R_CheckUserInterrupt();
		for(j=0;j<J;j++) {								
			obs=0; 
			if(t<T-1) {					
				for(u=1;u<=min(t+1,maxk);u++) {
					if(u<t+1) {
						if(DL>maxk){ 
							prod = obs+m_e(d, t*maxk+u-1, j, DL)+m_e(si, t-u+1, j, T);
						} else {
							prod = obs+m_e(d, u-1, j, maxk)+m_e(si, t-u+1, j, T);
						}
						obs+=m_e(b, t-u, j, T);
					}	else {	
						if(DL>maxk){ 
							prod = obs+m_e(d, t*maxk+t, j, DL)+pi[j];
						} else {
							prod = obs+m_e(d, t, j, maxk)+pi[j];
						}
					}
					if(u==1 || maxprod < prod) {
						maxprod = prod;
						m_e(opt, t, j, T)=u;
					}	
				}			
			} else {
				for(u=1;u<=min(t+1,maxk);u++) {
					if(u<T) {
						if (DL>maxk){
							prod =obs+m_e(D, t*maxk+u-1, j, DL)+m_e(si, T-u, j, T); 
						} else {
							prod =obs+m_e(D, u-1, j, maxk)+m_e(si, T-u, j, T);  
						}					
						obs += m_e(b, T-1-u, j, T);
					} else {	
						if(DL>maxk){
							prod=obs+m_e(D, t*maxk+T-1, j, DL)+pi[j];					 
						} else {
							prod=obs+m_e(D, T-1, j, maxk)+pi[j];							
						}
					}
					if(u==1 || maxprod < prod) {
						maxprod = prod;
						m_e(opt, t, j, T)=u;
					}						
				}
			}
			m_e(alpha, t, j, T) = maxprod+m_e(b, t, j, T);	
		}

		if(t<T-1) {
			for(j=0;j<J;j++){
				i=0;
				m_e(si, t+1, j, T)=m_e(a, i, j, J)+m_e(alpha, t, i, T);
				m_e(ops, t+1, j, T)=0;
				for(i=1;i<J;i++){
					if(i!=j) {
						prod = m_e(a, i, j, J)+m_e(alpha, t, i, T);
						if(m_e(si, t+1, j, T) <= prod) {
							m_e(si, t+1, j, T) = prod;
							m_e(ops, t+1, j, T)=i;
						}
					}
				} 
			}
		}
	}

	shat[T-1] = 0;
	for(j=1;j<J;j++) {
		if(m_e(alpha, T-1, shat[T-1] , T) < m_e(alpha, T-1, j, T)) shat[T-1] = j;
	}
	u=1;
	for(t=T-2;t>=0;t--) {
		if(u < m_e(opt, t+u, shat[t+u], T)) {
			shat[t] = shat[t+u];
			u++;
		} else {
			shat[t] = m_e(ops, t+u, shat[t+u], T);
			u=1;
		}
	}
}

static R_NativePrimitiveArgType univaRseg_t[] = {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType backward_t[] = {REALSXP, REALSXP, REALSXP, REALSXP,REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType viterbi_t[] = {REALSXP, REALSXP, REALSXP, REALSXP,REALSXP,INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP};

static R_NativePrimitiveArgType logviterbi_t[] = {REALSXP, REALSXP, REALSXP, REALSXP,REALSXP,INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP};

static const R_CMethodDef cMethods[] = {
  {"univaRseg", (DL_FUNC) &univaRseg, 8, univaRseg_t},
  {"backward", (DL_FUNC) &backward, 18, backward_t},
  {"viterbi", (DL_FUNC) &viterbi, 14, viterbi_t},
  {"logviterbi", (DL_FUNC) &logviterbi, 14, logviterbi_t},
  {NULL, NULL, 0}
};

void R_init_biomvRCNS(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

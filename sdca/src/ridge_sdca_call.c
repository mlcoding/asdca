#include "mymath.h"

SEXP ridge_sdca_call(SEXP Y_, SEXP X_, SEXP X_maxrn_, SEXP n_, SEXP d_, SEXP lambda_, SEXP nnlambda_, SEXP mmax_ite_, SEXP pprec_, SEXP idx_){
//void ridge_sdca_call(double *Y, double * X, double * X_maxrn, double * beta, double * intcpt, int * nn, int * dd, int * ite_lamb, int * ite_in, double *runt, double *lambda, int *nnlambda, int *mmax_ite, double *pprec, int *idx){
    
    int i, j, max_ite2, ite1, ite2, c_idx;
    double prec2, ilambda, dif1, dif2, dbn;
	double R, kappa, mu, rho, eta, bt, lambdak, one, delta_alp;
    clock_t start, stop;
    
    double X_maxrn = REAL(X_maxrn_)[0];
    int n = INTEGER(n_)[0];
    int d = INTEGER(d_)[0];
    int nlambda = INTEGER(nnlambda_)[0];
    int max_ite1 = INTEGER(mmax_ite_)[0];
    max_ite2 = n*3;
    double prec1 = REAL(pprec_)[0];
    prec2 = 1e-3;
    dbn = (double)n;
    one = 1;
    
    double *X = REAL(X_);
    double *Y = REAL(Y_);
    double *lambda = REAL(lambda_);
    int *idx = INTEGER(idx_);
    
    SEXP beta0, intcpt0, ite_lamb0, ite_in0, runt0;
    PROTECT(beta0 = allocVector(REALSXP, d*nlambda));
    double *beta = REAL(beta0);
    PROTECT(intcpt0 = allocVector(REALSXP, nlambda));
    double *intcpt = REAL(intcpt0);
    PROTECT(ite_lamb0 = allocVector(INTSXP, nlambda));
    int *ite_lamb = INTEGER(ite_lamb0);
    PROTECT(ite_in0 = allocVector(INTSXP, nlambda));
    int *ite_in = INTEGER(ite_in0);
    PROTECT(runt0 = allocVector(REALSXP, nlambda));
    double *runt = REAL(runt0);
    
    double *beta2 = (double *) Calloc(d, double);
    double *beta1 = (double *) Calloc(d, double);
    double *betay = (double *) Calloc(d, double);
    double *v = (double *) Calloc(d, double);
    double *z = (double *) Calloc(d, double);
    double *alp = (double *) Calloc(n, double);
    double *y_tild = (double *) Calloc(n, double);
    //int *idx = (int *) Calloc(max_ite2, int);
    for(i=0;i<d;i++){
        beta2[i] = 0;
        beta1[i] = 0;
        betay[i] = 0;
    }
    for(i=0;i<n;i++){
        alp[i] = 0;
    }
    for(i=0;i<d*nlambda;i++){
        beta[i] = 0;
    }
    for(i=0;i<nlambda;i++){
        intcpt[i] = 0;
        ite_lamb[i] = 0;
        ite_in[i] = 0;
        runt[i] = 0;
    }
    //for(i=0;i<max_ite2;i++){
        //idx[i] = i%n;
    //}
	//printf("idx: 0=%d,5=%d,203=%d,416=%d \n",idx[0],idx[5],idx[203],idx[416]);
    start = clock();
    
    for (i=0; i<nlambda; i++) {
		intcpt[i] = 0;
        ilambda = lambda[i];
		R = max(X_maxrn, sqrt(11*dbn*ilambda)); // parameter set-up
		kappa = R*R/dbn - ilambda;
		mu = ilambda/2;
		rho = mu + kappa;
		eta = sqrt(mu/rho);
		bt = (1-eta)/(1+eta);
		lambdak = ilambda+kappa;
        ite1 = 0;
        dif1 = 1;
        while (dif1>prec1 && ite1<max_ite1) {
			prod_vec_const(z, betay, kappa/lambdak, d); // z=betay*kappa/(lambda+kappa)
			get_residual_dense(y_tild, Y, X, z, n, d); // y_tild = Y - X * z;
			vec_mat_prod_coef(v, alp, X, 1/(lambdak*dbn), n, d); // v = X^T alp / ((lambda+kappa)*n);
			//shuffle(idx, max_ite2);
			ite2 = 0;
			dif2 = 1;
			while (dif2>prec2 && ite2<max_ite2) {
                c_idx = idx[ite2];
                //if(c_idx<0 || c_idx>n-1) printf("c_idx=%d \n",c_idx);
				delta_alp = -(alp[c_idx]+vec_matrow_inprod(X+c_idx, v, one, n, d)-y_tild[c_idx])/(1+norm2sq_matrow_coef(X+c_idx, 1/(lambdak*dbn), n, d));
				//if(i==0&&ite1==0&&ite2==2)printf("y_tild=%f,Xv=%f,X2=%f,c_idx=%d,delta_alp=%f,v=%f,beta2=%f \n",y_tild[0],vec_matrow_inprod(X+c_idx, v, one, n, d),norm2sq_matrow_coef(X+c_idx, 1/(lambdak*dbn), n, d),c_idx,delta_alp,v[0],beta2[0]);
				alp[c_idx] = alp[c_idx]+delta_alp;
				sum_vec_matrow(v, X+c_idx, delta_alp/(lambdak*dbn), n, d); // v = v+X[c_idx,]*delta_alp/(lambdak*dbn)
				sum_vec(beta2,v,z,d); //beta2 = v+z
				dif2 = (residual_l2sq_dense(Y,X,beta2,n,d)+vec_sum_l2sq_dense(alp,Y,n)-vec_2normsq(Y,n))/(2*dbn)+lambdak*vec_inprod(beta2,v,d);
				ite2++;
				//if(ite2%100==0) printf("i=%d,ite1=%d,ite2=%d,dif2=%f \n",i,ite1,ite2,dif2);
			}
			ite_in[i] += ite2;
			dif1 = dif_2norm_dense(beta1, beta2, d); // ||beta1-beta2||_2
			//printf("dif1=%f \n",dif1);
			for(j=0;j<d;j++){
				betay[j] = beta2[j]+bt*(beta2[j]-beta1[j]);
				beta1[j] = beta2[j];
			}
			ite1++;
		}
		ite_lamb[i] = ite1;
		for(j=0;j<d;j++){
			beta[i*d+j] = beta2[j];
        }
        printf("i=%d,ite1=%d,ite2=%d,b1=%f,b2=%f,b3=%f,b4=%f,b5=%f \n",i,ite1,ite_in[i],beta[i*d+0],beta[i*d+1],beta[i*d+2],beta[i*d+3],beta[i*d+4]);
        stop = clock();
        runt[i] = (double)(stop - start)/CLOCKS_PER_SEC;
    }
    Free(beta2);
    Free(beta1);
    Free(betay);
	Free(v);
    Free(z);
	Free(alp);
    Free(y_tild);
    //Free(idx);
    SEXP res;
    PROTECT(res = allocVector(VECSXP, 5));
    SET_VECTOR_ELT(res, 0, beta0);
    SET_VECTOR_ELT(res, 1, intcpt0);
    SET_VECTOR_ELT(res, 2, ite_lamb0);
    SET_VECTOR_ELT(res, 3, ite_in0);
    SET_VECTOR_ELT(res, 4, runt0);
    UNPROTECT(6);
    return(res);
}

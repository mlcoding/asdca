#include "mymath.h"

void stocml_svm_sdca_idx(double *Y, double * X, double * X_maxrn, double * beta, double * intcpt, int * nn, int * dd, int * ite_lamb, int * ite_in, double *runt, double *lambda, int *nnlambda, int *mmax_ite, double *pprec, double * X_rn, double * dual){
  
  int i, j, n, d, max_ite1, max_ite2, nlambda, ite1, ite2, c_idx;
  double prec1, prec2, ilambda, dif1, dif2, dbn, gamma, epst;
  double R, kappa, mu, rho, eta, bt, lambdak, one, delta_alp;
  clock_t start, stop;
  
  n = *nn;
  d = *dd;
  max_ite1 = *mmax_ite;
  max_ite2 = n*5;
  prec1 = *pprec;
  prec2 = 1e-3;
  nlambda = *nnlambda;
  dbn = (double)n;
  one = 1;
  
  double *beta2 = (double *) Calloc(d, double);
  double *beta1 = (double *) Calloc(d, double);
  double *betay = (double *) Calloc(d, double);
  double *Xb = (double *) Calloc(n, double);
  double *z = (double *) Calloc(d, double);
  double *alp = (double *) Calloc(n, double);
  double *y_tild = (double *) Calloc(n, double);
  int *idx = (int *) Calloc(max_ite2, int);
  for(i=0;i<n;i++){
    idx[i] = i;
    /*     printf("%d ",i);
    if(i%50 ==0) printf("\n"); */
  }
  //printf("idx: 0=%d,5=%d,203=%d,416=%d \n",idx[0],idx[5],idx[203],idx[416]);
  start = clock();
  
  for(j=0;j<d;j++){
    beta2[j] = 0;
    beta1[j] = 0;
    betay[j] = 0;
  }
  for (i=0; i<nlambda; i++) {
    for(j=0;j<n;j++){
      alp[j] = 0;
    }
    intcpt[i] = 0;
    ilambda = lambda[i];
    gamma = prec2;
    R = max(*X_maxrn, sqrt(11*dbn*ilambda*1)); // parameter set-up
    kappa = R*R/(dbn*1)- ilambda;
    mu = ilambda/2;
    rho = mu + kappa;
    eta = sqrt(mu/rho);
    bt = (1-eta)/(1+eta);
    lambdak = ilambda+kappa;
    epst = pd_gap_svm1(X,beta2,Xb,alp,gamma,ilambda,n,d)*(1+pow(eta,-2));
    ite1 = 0;
    dif1 = 1;
    max_ite1 = (int)(1+2*log(epst/prec1)/eta)+1;
	//printf("ilambda=%f,gamma=%f,R=%f,kappa=%f,mu=%f,rho=%f,eta=%f,bt=%f,epst=%f,max_ite1=%d \n",ilambda,gamma,R,kappa,mu,rho,eta,bt,epst,max_ite1);
	//break;
    while (dif1>prec1 && ite1<max_ite1) {
      prod_vec_const(z, betay, kappa/lambdak, d); // z=betay*kappa/(lambda+kappa)
      //shuffle(idx, n);
      prec2 = epst*eta/(2*(1+pow(eta,-2)));
      ite2 = 0;
      dif2 = 1;
	  //printf("z1==%f,z1==%f,z1==%f,z1==%f,z1==%f,prec2=%f",z[0],z[1],z[2],z[3],z[4],prec2);
	  //break;
      while (dif2>prec2 && ite2<max_ite2) {
        c_idx = idx[ite2%n];
        //if(c_idx<0 || c_idx>n-1) printf("c_idx=%d \n",c_idx);
        delta_alp = max(-alp[c_idx],min(1-alp[c_idx],(1-vec_matrow_inprod(X+c_idx, beta2, one, n, d)-gamma*alp[c_idx])/(X_rn[c_idx]/(ilambda*dbn)+gamma)));
        //if(i==0&&ite1==0&&ite2==2)printf("y_tild=%f,Xv=%f,X2=%f,c_idx=%d,delta_alp=%f,v=%f,beta2=%f \n",y_tild[0],vec_matrow_inprod(X+c_idx, v, one, n, d),norm2sq_matrow_coef(X+c_idx, 1/(lambdak*dbn), n, d),c_idx,delta_alp,v[0],beta2[0]);
        alp[c_idx] = alp[c_idx]+delta_alp;
        sum_vec_matrow(beta2, X+c_idx, delta_alp/(ilambda*dbn), n, d); // beta2 = beta2+X[c_idx,]*delta_alp/(lambdak*dbn)
		//printf("beta2 %f,%f,%f,%f \n",beta2[0],beta2[1],beta2[2],beta2[3]);
        dif2 = pd_gap_svm2(X,beta2,Xb,alp,z,gamma,ilambda,n,d);
        ite2++;
        //printf("c_idx=%d,delta_alp=%f,dif2=%f,gamma=%f \n",c_idx,delta_alp,dif2,gamma);
		//break;
      }
      ite_in[i] += ite2;
      epst = pow(1-eta/2,ite1)*epst;
      //dif1 = dif_2norm_dense(beta1, beta2, d); // ||beta1-beta2||_2
      dif1 = (1+rho/mu)*epst+(rho*kappa*dif_2norm_dense(beta2,betay,d))/(2*mu);
      //printf("dif1=%f \n",dif1);
      for(j=0;j<d;j++){
        betay[j] = beta2[j]+bt*(beta2[j]-beta1[j]);
        beta1[j] = beta2[j];
      }
      ite1++;
	  //break;
    }
    ite_lamb[i] = ite1;
    for(j=0;j<d;j++){
      beta[i*d+j] = beta2[j];
    }
    for(j=0;j<n;j++){
      dual[i*n+j] = alp[j];
    }
    //printf("i=%d,ite1=%d,b1=%f,b2=%f,b3=%f,b4=%f,b5=%f \n",i,ite1,beta[i*d+0],beta[i*d+1],beta[i*d+2],beta[i*d+3],beta[i*d+4]);
    stop = clock();
    runt[i] = (double)(stop - start)/CLOCKS_PER_SEC;
	//break;
  }
  Free(beta2);
  Free(beta1);
  Free(betay);
  Free(Xb);
  Free(z);
  Free(alp);
  Free(y_tild);
  Free(idx);
}
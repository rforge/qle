/**
 * @file        optimize.cc
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       R (Interface) Functions for optimization
 *
 * Explanation: Quasi-Fisher scoring iteration
 *
 */

#include "qsoptim.h"

#include <R_ext/Applic.h>
#include <R_ext/Constants.h>

#define STPMAX 100
#define TOLSTEP 1e-11

typedef void (*fnCompare) (double *x, void *data, double *f, int *info);

qfs_result qfscoring(double *x, int n, double *f,
			 qfs_options cond, int *info);

void backtr(int n, double *xold, double fold, double *d, double *g, double *x,
			double *f, int *check, fnCompare *monitor, fnCompare *fnold, int *fntype,
			 double stepmax, double *slope, double *delta, void *data, int *info);

// projection into box constraints
double med(double x, double y, double z, int *info);
void projmid(double *x, int nx, double *lb, double *ub, int *info);

// final status of QS iteration
SEXP getStatus( qfs_result status );

/** \brief Compute (projected) norm of quasi-score
 *
 * @param x 	variable
 * @param data  data pointer
 * @param f 	pointer to function value
 * @param info  integer, >0 if it reaches any bound constraints
 */
void fnQSNorm(double *x, void *data, double *f, int *info) {
  ql_model qlm = (ql_model) data;
  projmid(x,qlm->glkm->dx,qlm->lower,qlm->upper, info);
  *f = qlm->intern_qfScoreNorm(x);
  if(qlm->info > 0L)
	WRR("Could not compute monitor function value (norm of quasi-score) for line search in quasi-scoring iteration.")
}

void fnQDev(double *x, void *data, double *f, int *info) {
  ql_model qlm = (ql_model) data;
  projmid(x,qlm->glkm->dx,qlm->lower,qlm->upper, info);
  *f = qlm->intern_qfScoreStat(x);
  if(qlm->info > 0L)
  	WRR("Could not compute monitor function value (quasi-deviance) for line search in quasi-scoring iteration.")
}


/** \brief C Interface:
 *      Local root finding of score equations
 *
 * @param R_start start vector
 * @param R_args argument list for kriging
 * @param R_opt options for qfs iteration
 *
 * @return result object
 */

SEXP QSopt(SEXP R_start, SEXP R_qsd, SEXP R_qlopts, SEXP R_X, SEXP R_Vmat, SEXP R_cm, SEXP R_opt) {

    int nProtected = 0, info = 0,
		xdim = LENGTH(R_start);

    /*! set start vector */
    SEXP R_sol = R_NilValue;
    PROTECT(R_sol = allocVector(REALSXP, xdim));
    ++nProtected;

    double *xsol = REAL(R_sol);
    double *start = REAL(AS_NUMERIC(R_start));
    MEMCPY(xsol,start,xdim);

    ql_model_t qlm(R_qsd, R_qlopts, R_X, R_Vmat, R_cm, COPY_ZERO);

    /*! Set optimization options*/
    qfs_options_t qfs(&qlm,R_opt);

    /* scoring iteration */
    double fval=HUGE_VAL, qd=HUGE_VAL;
    qfs_result status = qfscoring(xsol,xdim,&fval,&qfs,&info);

    /* return objects */
    SEXP R_S, R_jac, R_I, R_Iobs=R_NilValue;
    PROTECT(R_S = allocVector(REALSXP,xdim));
    ++nProtected;
    PROTECT(R_I = allocMatrix(REALSXP,xdim,xdim));
    ++nProtected;
    PROTECT(R_jac = allocMatrix(REALSXP,xdim,qlm.nCov));
    ++nProtected;

    /* copy results  */
    MEMCPY(REAL(R_S),qlm.score,xdim);
    MEMCPY(REAL(R_I),qlm.qimat,xdim*xdim);
    MEMCPY(REAL(R_jac),qlm.jac,xdim*qlm.nCov);
    if(qfs.doIobs){
      	PROTECT(R_Iobs = allocMatrix(REALSXP,xdim,xdim));
       	++nProtected;
       	info = qlm.intern_quasiObs(xsol, qlm.score, REAL(R_Iobs));
       	if( info != NO_ERROR )
       	  XWRR( info, "inter_quasiObs")
    }

    /* add prediction variances to return list */
    SEXP R_sig2 = R_NilValue, R_varS = R_NilValue;

    // names variance matrix
    SEXP R_VmatNames;
    PROTECT(R_VmatNames = allocVector(VECSXP,2));
    ++nProtected;
    SEXP R_obs = getListElement( R_qsd, "obs" );
    SET_DIMNAMES_MATRIX(R_VmatNames,R_obs)

    if(qlm.glkm->krigType) {
      	PROTECT(R_sig2 = allocVector(REALSXP,qlm.nCov));
      	++nProtected;
      	MEMCPY(REAL(R_sig2),qlm.glkm->krigr[0]->sig2,qlm.nCov);
       	setAttrib(R_sig2,R_NamesSymbol,getAttrib(R_obs,R_NamesSymbol));

       	/* variance quasi-score */
       	PROTECT(R_varS = allocMatrix(REALSXP,qlm.dx,qlm.dx));
       	++nProtected;
       	qlm.intern_varScore(REAL(R_varS));
    }

    // compute quasi-deviance
    qd = qlm.qfValue(qlm.score,qlm.qimat);

#ifdef DEBUG
    Rprintf("value: %f \n", qd);
    Rprintf("Qnorm: %f \n", fval);
    printMatrix("vmat",qlm.qld->vmat, &xdim,&xdim);
    printMatrix("jac",REAL(R_jac), &xdim,&qlm.nCov);
    printMatrix("I",REAL(R_I), &xdim,&xdim);
    printVector("start:", xsol, &xdim);
    printVector("score:", REAL(R_S), &xdim);
#endif

    SEXP R_dimnames = R_NilValue;
    PROTECT(R_dimnames = allocVector(VECSXP,2));
    ++nProtected;
    SET_DIMNAMES_MATRIX(R_dimnames,R_start)
    setAttrib(R_I, R_DimNamesSymbol, R_dimnames);
    setAttrib(R_sol, R_NamesSymbol, getAttrib(R_start,R_NamesSymbol));

    static const char *nms[] =
     {"convergence", "message", "iter", "value", "par",
     "score", "sig2", "I", "Iobs", "varS", "start", "Qnorm",
	 "method", "criterion", ""};

    SEXP R_ret = R_NilValue;
    PROTECT(R_ret = mkNamed(VECSXP, nms));
    ++nProtected;

    SET_VECTOR_ELT(R_ret, 0, ScalarInteger((int)status));
    SET_VECTOR_ELT(R_ret, 1, getStatus(status));
    SET_VECTOR_ELT(R_ret, 2, ScalarInteger(qfs.num_iter));
    SET_VECTOR_ELT(R_ret, 3, ScalarReal(qd));
    SET_VECTOR_ELT(R_ret, 4, R_sol);
    SET_VECTOR_ELT(R_ret, 5, R_S);
    SET_VECTOR_ELT(R_ret, 6, R_sig2);
    SET_VECTOR_ELT(R_ret, 7, R_I);
    SET_VECTOR_ELT(R_ret, 8, R_Iobs);
    SET_VECTOR_ELT(R_ret, 9, R_varS);
    SET_VECTOR_ELT(R_ret, 10, R_start);
    SET_VECTOR_ELT(R_ret, 11, ScalarReal(fval));
    SET_VECTOR_ELT(R_ret, 12, mkString("qscoring"));
    SET_VECTOR_ELT(R_ret, 13, mkString("qle"));
    setVmatAttrib(&qlm, R_VmatNames, R_ret);
    SET_CLASS_NAME(R_ret,"QSResult")

    UNPROTECT(nProtected);
    return R_ret;
}


#define FREE_WORK { \
   FREE(g)     	    \
   FREE(d)			\
   FREE(xold)  	    \
}

/** \brief  Quasi-Fisher scoring iteration
 *        Comment: Either use step length equal to 1
 *                 or a line search based on the Fisher-Score statistic
 *
 * @param x start vector
 * @param n dimension of start vector
 * @param f monitor function value at solution
 * @param fnMonitor quasi-ddeviance function
 * @param qfs options
 * @param info information code
 *
 * @return result flag
 */

qfs_result qfscoring(double *x,			 	/* start */
					 int n,      		 	/* parameter length */
					 double *f,  		 	/* objective value */
					 qfs_options qfs,    	/* options for scoring */
					 int *info)
{
	// temp pointers
   ql_model qlm = qfs->qlm;
   ql_storage_t qlsolve = qlm->qlsolve;
   // default return status
   qfs_result status = QFS_NO_CONVERGENCE;
   /*! start with quasi-score norm as a monitor function*/
   fnCompare fnMonitor = &fnQSNorm,
		    fnMonitor2 = &fnQDev;
   // calculate value, score, quasi-info, etc.
   fnMonitor(x, qlm, f, info);

   /*! test for f first */
   if(*f < qfs->ftol_stop ){
     *info=QFS_STOPVAL_REACHED;
     return QFS_CONVERGENCE;
   }

   int i=0, niter=0, check=0, fntype=0,
       Nmax=qfs->max_iter, pl=qfs->pl;

   double fold=0, test=0, tmp=0,
		  delta=1.0, slope=0, den=0;

   double *d=0, *xold=0, *g=0;
   CALLOCX(g,n,double);
   CALLOCX(d,n,double);
   CALLOCX(xold,n,double);

   // score vector from `fnMonitor`
   double *qimat = qlm->qimat;
   double *score = qlm->score;

   test=0.0;
   for (i=0;i<n; ++i) {
      if (std::fabs(score[i]) > test)
        test=std::fabs(score[i]);
   }
   if (test < qfs->score_tol) {
      FREE_WORK
      *info=QFS_SCORETOL_REACHED;
      return QFS_CONVERGENCE;
    }

   for (test=0.0,i=0;i<n;++i)
	 test += SQR(x[i]);
   double stepmax = STPMAX*MAX(std::sqrt(test),(double)n);

   /*! optimization loop */
   for(niter=0; niter < Nmax; ++niter)
   {
	     fold = *f;
	     MEMCPY(xold,x,n);
	     MEMCPY(d,score,n);
	     // solve for direction d = I^{-1}Q < 0
	     gsiSolve(qimat,n,d,1,qlsolve.qimat,info,Chol);
         if(*info != 0){
        	 status=QFS_BAD_DIRECTION;
        	 PRINT_MSG("Cannot compute quasi-score correction vector.")
        	 LOG_ERROR(LAPACK_SOLVE_ERROR, "gsiSolve");
        	 break;
         }
         // compute (approximate) gradient of 1/2||Q||
         // using QI as the Jacobian of QS
         for(i=0;i<n;++i)
          qlsolve.score[i] = -score[i];
         matmult(qimat,n,n,qlsolve.score,n,1,g,info);
         if(*info > 0){
          LOG_WARNING(*info,"matmult")
		  WRR("Could not compute gradient of monitor function.")
		  // just use quasi-score as a gradient if it fails
		  MEMCPY(g,score,n);
         }
         // line search: dynamically switch between both monitor functions
         backtr(n,xold,fold,d,g,x,f,&check,&fnMonitor,&fnMonitor2,
        		 &fntype,stepmax,&slope,&delta,(void*)qlm,info);

         /*! display information */
         if(pl >= 10) {
           Rprintf("Quasi-scoring:\n\n");
           Rprintf("iteration......... %d \n", niter);
           Rprintf("at bounds......... %d \n", *info);
           Rprintf("objective......... %3.12f \n", *f);
           Rprintf("step size......... %3.12f (check=%d) \n", delta, check);
           Rprintf("slope............ %3.12f (monitor=%d) \n\n", slope, fntype);
           printVector("current", x, &n);
           Rprintf("\n");
           printVector("score", score, &n);
           Rprintf("\n");
           printVector("direction (d)", d, &n);
           Rprintf(" length: %3.12f\n\n", denorm(d,n));
           printVector("gradient (g)", g, &n);
           Rprintf("\n");
           Rprintf("--------------------------------------------------------------\n\n");
         }

         /*! test for score being zero */
         test=0.0;
         for (i=0;i<n;++i) {
        	tmp=std::fabs(score[i]);
            if(tmp > test)
              test=tmp;
         }
         if(test < qfs->score_tol) {
           FREE_WORK
		   qfs->num_iter = niter;
           status=QFS_SCORETOL_REACHED;
           return status;
         }
         /*! test for f */
         if(*f < qfs->ftol_stop) {
             FREE_WORK
			 qfs->num_iter=niter;
             return QFS_STOPVAL_REACHED;
         }

         /*! test for local min */
         if(check > 0) {
        	 test=0.0;
        	 den=MAX(*f,0.5*n);
        	 for (i=0;i<n;++i) {
        	   tmp=std::fabs(g[i])*MAX(std::fabs(x[i]),1.0)/den;
        	   if(tmp > test)
        		 test=tmp;
        	 }
        	 /* use gradient for testing local convergence */
        	 if(test < qfs->grad_tol) {
        	    status=(*f < qfs->ftol_abs ? QFS_CONVERGENCE : QFS_GRADTOL_REACHED);
			 } else {
				status=QFS_LINESEARCH_FAILURE;
			 }
			 FREE_WORK
			 qfs->num_iter=niter;
			 return status;

         } else {  /* line search success */

			 test=(std::fabs(*f-fold))/MAX(std::fabs(*f),1.0);
			 if (test < qfs->ftol_rel) {
				 FREE_WORK
				 qfs->num_iter=niter;
				 return QFS_FTOLREL_REACHED;
			 }
			 /*! test for relative change in x */
			 test=0.0;
			 for (i=0;i<n; ++i) {
				   tmp=(std::fabs(x[i]-xold[i]))/MAX(std::fabs(x[i]),1.0);
				   if(tmp > test)
				    test=tmp;
			 }
			 if(test < qfs->xtol_rel && test > 0) {
				FREE_WORK
				qfs->num_iter=niter;
				return QFS_XTOL_REACHED;
			 }
			 /*! test for zero slope if monitor is norm of QS */
			 if(!fntype && std::fabs(slope) < qfs->slope_tol) {
				FREE_WORK
				qfs->num_iter=niter;
				return QFS_SLOPETOL_REACHED;
			 }

         }
   } /*! end for */

   if(niter == Nmax)
	status=QFS_MAXITER_REACHED;

   FREE_WORK
   qfs->num_iter=niter;
   return status;
}

#undef FREE_WORK

void backtr(int n, double *xold, double fold,  double *d, double *g, double *x,
		double *f, int *check, fnCompare *fn, fnCompare *fn2, int *fntype, double stepmax,
		 double *slope, double *delta, void *data, int *info) {

  int i=0, flag=0, type=*fntype;
  double s=1.0, tmp=0.0, alpha=1e-3, slope2=1.0;

  fnCompare monitor = *fn;							// temp. monitor
  double sum = denorm(d,n);
  if (sum > stepmax)								// not too big steps
   for (i=0;i<=n;i++)
	d[i] *= stepmax/sum;

  slope2=( type ? -fold : innerProduct(g,d,n)); 	// equals monitor if
  *slope=slope2;       						     	// quasi-deviance is used
  if(!R_FINITE(slope2) || slope2 >= 0.0){
  	   WRR("Roundoff problems in line search.")
  	   *check=1;
  	   return;
  }

  *check=0;
  for(;;) {
	  *delta=s;
	  for (i=0;i<n;++i)
	    x[i]=xold[i]+s*d[i];
	  (*fn)(x, data, f, info);
	  tmp=alpha*s*slope2;							// < 0 always
	  if(*f <= fold + tmp) return;					// found valid step
	  if(std::fabs(MIN(s,tmp)) < TOLSTEP){
		  if(flag){
			  *check=1;
			  for (i=0;i<n; ++i)
			    x[i]=xold[i];						// reset current iterate
			  (*fn)(x, data, f, info);				// and compute things a last time
			  return;
		  } else {
			 *fn=*fn2;								// switch monitor
			 *fn2=monitor;
			 // use new monitor for old value
			 (*fn)(xold, data, &fold, info);		// recompute monitor at fold
			 if(type){								// set new slope parameter
				type=0;								// use norm of quasi-score next iteration of LS
			 	slope2=innerProduct(g,d,n);			// would be the slope of norm QS if we had the Jacobian of QS (use QI instead as approximation)
			 } else {
				type=1;								// use quasi-deviance next iteration of LS
				slope2=-fold;						// heuristic: some required decrease in monitor function QD
			 }
			 *slope=slope2;
			 if(!R_FINITE(slope2) || slope2 >= 0.0){
			   WRR("Roundoff problems in line search.")
			   *check=1;
			   return;
			 }
			 s=1.0;
			 flag=1;								// changed the monitor function
			 *fntype=type;							// store current type of monitor
		  }
	  } else s=0.5*s;								// decrease step
  }
}

double med(double x, double y, double z, int *info) {
   if ( (x - y) * (z - x) >= 0 ) {
      if((x - y) * (z - x) == 0)
       *info=1;
      return x;
   } else if ( (y - x) * (z - y) >= 0 ) {
       *info=1;
       return y;
   } else {
       *info=1;
       return z;
  }
}


void projmid(double *xproj, int nx, double *lb, double *ub, int *info) {
  *info=0;
  for(int i=0;i < nx; ++i)
    xproj[i] = med(xproj[i],lb[i],ub[i], info);
}


/** \brief Convert status message
 *
 * @param status
 * @return message
 */
SEXP getStatus( qfs_result status ) {
   // convert message to an R object
   SEXP R_message;
   PROTECT(R_message = allocVector(STRSXP, 1));

   switch ( status ) {
       // (= +0)
       case QFS_CONVERGENCE:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_CONVERGENCE: Optimization stopped because of approximate root."));
           break;
       // (= +1)
       case QFS_SCORETOL_REACHED:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_SCORETOL_REACHED: Optimization stopped because score_tol was reached."));
           break;
       // (= +2)
       case QFS_FTOLREL_REACHED:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_FTOLREL_REACHED: Optimization stopped because ftol_rel was reached."));
           break;
       // (= +3)
       case QFS_STOPVAL_REACHED:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_STOPVAL_REACHED: Optimization stopped because ftol_stop or ftol_abs was reached."));
           break;
       // (= +4)
       case QFS_XTOL_REACHED:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_XTOL_REACHED: Optimization stopped because xtol_rel was reached."));
           break;
       // (= +5)
       case QFS_GRADTOL_REACHED:
            SET_STRING_ELT(R_message, 0, mkChar("QFS_GRAD_TOL_REACHED: Optimization stopped because grad_tol was reached."));
            break;
	   // (= +6)
	   case QFS_SLOPETOL_REACHED:
			SET_STRING_ELT(R_message, 0, mkChar("QFS_SLOPETOL_REACHED: Optimization stopped because slope_tol was reached."));
			break;

       // (= +10)
       case QFS_LOCAL_CONVERGENCE:
            SET_STRING_ELT(R_message, 0, mkChar("QFS_LOCAL_CONVERGENCE: Optimization stopped because of local convergence."));
            break;

       /* Error codes (negative return values): */

       // (= -1)
       case QFS_NO_CONVERGENCE:
		   SET_STRING_ELT(R_message, 0, mkChar("QFS_NO_CONVERGENCE: Optimization stopped because no convergence could be detected."));
		   break;
       // (= -2)
       case QFS_BAD_DIRECTION:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_FAILURE: Could not calculate search direction."));
           break;
       case QFS_LINESEARCH_FAILURE:
            SET_STRING_ELT(R_message, 0, mkChar("QFS_LINESEARCH_FAILURE: Optimization stopped because of line search failure."));
            break;
       case QFS_LINESEARCH_ZEROSLOPE:
            SET_STRING_ELT(R_message, 0, mkChar("QFS_LINESEARCH_ZEROSLOPE: Optimization stopped because of nearly zero slope during line search."));
            break;
       case QFS_ERROR:
		    SET_STRING_ELT(R_message, 0, mkChar("QFS_FAILURE: Generic failure code."));
		    break;
       case QFS_MAXITER_REACHED:
            SET_STRING_ELT(R_message, 0, mkChar("QFS_MAXTIME_REACHED: Optimization stopped because maxiter (above) was reached."));
            break;
       default:
           SET_STRING_ELT(R_message, 0, mkChar("Unknown return status."));
           break;
       }

   UNPROTECT( 1 );
   return R_message;
}

/**
 *  @brief:
 *   testing internal function `invMatrix`,
 *   do not forget to reset dynamic symbols
 *
 *  */
//SEXP invertMatrix(SEXP R_A,SEXP R_type)
//{
//	int err = 0;
//	int n = GET_DIMS(R_A)[0];
//	inversion_type type = (inversion_type) asInteger(AS_INTEGER(R_type));
//
//	SEXP R_B = PROTECT(allocMatrix(REALSXP, n, n));
//	SEXP R_AP = PROTECT(allocMatrix(REALSXP, n, n));
//
//	MEMCPY(REAL(R_AP),REAL(R_A),SQR(n));
//	invMatrix(REAL(R_AP), n, REAL(R_B), &err, type);
//
//  UNPROTECT(2);
//	return R_B;
//}

/**
 * @brief
 * 	 testing internal function `gsiSolve`
 *   do not forget to reset dynamic symbols
 */
//SEXP RSolve(SEXP R_A, SEXP R_B, SEXP R_type){
//	int err = 0,
//	 n = GET_DIMS(R_A)[0],
//	 m = GET_DIMS(R_B)[1];
//
//	Rprintf("n: %d, m: %d \n",n,m);
//	Rprintf("A -> [%d,%d] \n",GET_DIMS(R_A)[0],GET_DIMS(R_A)[1]);
//	Rprintf("B -> [%d,%d] \n",GET_DIMS(R_B)[0],GET_DIMS(R_B)[1]);
//
//	inversion_type type = (inversion_type) asInteger(AS_INTEGER(R_type));
//
//	SEXP R_AP = PROTECT(allocMatrix(REALSXP, n, n));
//	MEMCPY(REAL(R_AP),REAL(R_A),n*n);
//
//	SEXP R_Z = PROTECT(allocMatrix(REALSXP, n, m));
//	MEMCPY(REAL(R_Z),REAL(R_B),n*m);
//
//	double *Awork = CALLOC(SQR(n),double);
//	gsiSolve(REAL(R_AP), n, REAL(R_Z), m, Awork, &err, type);
//	FREE(Awork);
//
//	UNPROTECT(2);
//	return R_Z;
//}


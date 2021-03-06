/**
 * @file        kriging.cc
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       R (Interface) Functions for kriging estimation
 *
 * Explanation: Kriging statistics and derivatives
 *
 */

#ifndef KRIGING_CC_
#define KRIGING_CC_

#define FD_EPS 1e-4
#define KRIGE_TOLERANCE 1e-10

#include "kriging.h"

static int ONE_ELEMENT = 1;
static int ZERO_ELEMENT = 0;
static double ZERO_DBL = 0.0;
static ql_model qlm_global = NULL ;

int intern_dualKmat(double *C, int Cdim, double *F, int Fdim, double *K);
//
int intern_dualKrigingWeights(double *Cmat, int nlx, double *Fmat, int fddim, double *data, double *w);

//
//void solveKriging(double *Cinv,double *Fmat, double *X1mat, double *Qmat,
//      int lx, int fddim, double *sigma2, double *w, double *sigma0, double *f0 );


/** \brief Internal function: Dual kriging matrix K,
 *              Chilès, J.-P. and Delfiner, P. (2008) Frontmatter, in Geostatistics:
 *              Modeling Spatial Uncertainty, John Wiley & Sons, Inc., Hoboken, NJ, USA.
 *              Chapter 3.4.8
 *
 * @param C covariance matrix
 * @param Cdim dimension of covariance matrix
 * @param F trend matrix
 * @param Fdim columns of F
 * @param K dual kriging matrix
 */
int intern_dualKmat(double *C, int nc, double *F, int fddim, double *K) {
	int i=0, j=0, k=nc+fddim;

	for (j = 0; j < nc; j++)
	      for (i = 0; i < nc; i++)
			K[i+j*(nc+fddim)] = C[i+j*nc];

	for (j = 0; j < fddim; j++)
	      for (i = 0; i < nc; i++)
			K[nc*k+i+j*k] = F[i+j*nc];

	for (j = 0; j < fddim; j++)
	      for (i = 0; i < fddim; i++)
	          K[nc*k+nc+i+j*k] = 0;

	for (j = 0; j < nc; j++)
		for (i = 0; i < fddim; i++)
        	K[(i+nc)+j*(nc+fddim)] = F[j+i*nc]; 	// F^t

	return NO_ERROR;
}


/**\brief Interface to low level C functions
   *     - Get dual kriging weights for a single statistic
   *     - Comment: This function is independent of the
   *                underlying C data structure for the kriging models, which is not initialized
   *
   * @param R_Cmat Covariance matrix
   * @param R_Fmat Trend matrix
   * @param R_data Data frame of sampled statistic values at the design sites
   *
   * @return  R vector of kriging weights
   */

SEXP getDualKrigingWeights(SEXP R_Cmat, SEXP R_Fmat, SEXP R_data) {
    int nx = GET_DIMS(R_Cmat)[0],
        fddim = GET_DIMS(R_Fmat)[1];
    int info=0, dK = nx+fddim;

    SEXP R_weights;
    PROTECT(R_weights = NEW_NUMERIC(dK));
    if( (info = intern_dualKrigingWeights(REAL(R_Cmat), nx, REAL(R_Fmat), fddim, REAL(R_data), REAL(R_weights)))  != NO_ERROR){
    	PRINT_MSG("Failed to compute dual kriging weights.")
    	XWRR(info,"intern_dualKrigingWeights")
    }
    UNPROTECT(1);
    return R_weights;

}


/**\brief Kriging prediction: calculate dual kriging weights
   *     - Kriging of general data
   *     - Comment: The kriging variance is not calculated
   *
   * @param Cmat Covariance matrix (with estimated parameters)
   * @param nlx Rows/Cols of Cmat
   * @param Fmat trend matrix
   * @param fddim Columns of trend matrix
   * @param data Data vector
   * @param w The weights
   *
   * return (void)
   */

int intern_dualKrigingWeights(double *Cmat, int nlx, double *Fmat, int fddim, double *data, double *w) {
  int info=0,
	  dK=nlx+fddim;

  double *K = NULL;
  CALLOCX(K, dK*dK, double);
  info = intern_dualKmat(Cmat,nlx,Fmat,fddim,K);

  MEMZERO(w,dK);
  MEMCPY(w,data,nlx);

  solveLU(K,dK,w,ONE_ELEMENT,info);
  FREE(K)

  return info;
}

/**\brief Kriging prediction (dual)
 *
 *    Comment: Kriging variances not calculated
 *
 * @param s0           covariance vector between design sites and x0
 * @param f0           trend vector of x0
 * @param fddim        length of f0 (columns of trend matrix)
 * @param w  [IN]      kriging weights
 * @param zx [OUT]     predicted values: z*(x0)
 *
 * return (void)
 */

int intern_dualKrigingPrediction(double *s0, double *f0, int fddim, double *w, int nz, double *mean) {
  int i=0, have_na=0;
  double zx=0;

  for(i = 0; i < nz; i++)
	zx += w[i]*s0[i];
  for(i = 0; i < fddim; i++)
	zx += w[i+nz]*f0[i];
  if(!R_FINITE(zx))
   { have_na=1; }

  *mean = zx;
  return have_na;
}

SEXP estimateJacobian(SEXP R_Xmat, SEXP R_data, SEXP R_points, SEXP R_covT, SEXP R_krigtype) {
   if(!isMatrix(R_Xmat))
	 ERR("Expected matrix of sample points.");

   int *dimX = GET_DIMS(R_Xmat),
		nCov = LENGTH(R_covT);

   int info = 0,
	   dx = dimX[1],
   	   npoints = LENGTH(R_points);

   glkrig_models glkm(R_covT,R_Xmat,R_data,R_krigtype);

   SEXP R_names, R_jacobians, R_jac_element, R_jac;
   PROTECT(R_jacobians = allocVector(VECSXP,npoints));

   // names Jacobian
   SEXP R_dimT = R_NilValue;
   SEXP R_nI = VECTOR_ELT(R_points,0);
   PROTECT(R_dimT = allocVector(VECSXP,2));
   SET_DIMNAMES_MATRIX2(R_dimT,R_nI,R_data)

   double *point=0, *fdwork=0;
   CALLOCX(fdwork,nCov,double);

   //krig_result_s krigr(nCov,lx);
   for(int i=0; i < npoints; ++i)
   {
	  point = REAL(AS_NUMERIC(VECTOR_ELT(R_points,i)));
	  //printVector("point",point,&dx);
	  if ( (info = glkm.intern_kriging(point)) != NO_ERROR){
		XWRR(info,"intern_kriging")
		SET_VECTOR_ELT(R_jacobians,i,R_NilValue);
		continue;
	  }
	  PROTECT(R_jac = allocMatrix(REALSXP, dx, nCov));
	  setAttrib(R_jac, R_DimNamesSymbol, R_dimT);

	  // get jacobian
	  if( (info = glkm.intern_jacobian(point,REAL(R_jac),fdwork)) != NO_ERROR)
		 XWRR(info,"intern_jacobian")

	  PROTECT(R_jac_element = allocVector(VECSXP, 1));
	  PROTECT(R_names = allocVector(STRSXP, 1));

	  SET_STRING_ELT(R_names, 0, mkChar("mean"));
	  setAttrib(R_jac_element, R_NamesSymbol, R_names);

	  SET_VECTOR_ELT(R_jac_element,0,R_jac);
	  SET_VECTOR_ELT(R_jacobians,i,R_jac_element);
	  UNPROTECT(3);
   }
   FREE(fdwork)
   UNPROTECT(2);
   return R_jacobians;
}

#define CHECK_UNPROTECT(fun) {	\
  if( info != NO_ERROR ) {			\
	 XWRR( info, fun)				\
	 UNPROTECT( nprotect );			\
	 continue;						\
  }									\
}


/**\brief C Interface: Kriging of multiple data
   *
   * @param R_Xmat      sample (design) matrix
   * @param R_data      list of observation vectors (as a data.frame)
   * @param R_points    matrix of unobserved points
   * @param R_covSList  list of covariance structures
   *
   * return List of points,
   *        each containes a list of kriging mean and variance
   *        and a list of kriging weights for each point
   */


SEXP kriging(SEXP R_Xmat, SEXP R_data, SEXP R_points, SEXP R_covList, SEXP R_krigType)
{
     int i = 0, k = 0,
    	 npoints = 0, info = 0,
         nCov = length(R_covList),
        *dimX = GET_DIMS(R_Xmat);

     int lx = dimX[0]; 			// number of sample points as rows of matrix

     /* is matrix or vector */
     if(isMatrix(R_points))
       npoints = GET_DIMS(R_points)[0];  // number of points to  predict gradient
     else ERR("Expected sample points as a matrix.")

     /* init all kriging models */
     glkrig_models glkm(R_covList,R_Xmat,R_data,R_krigType);

     SEXP R_mean, R_sigma2, R_weights, R_nT,
	 	  R_weights_tmp, R_tmp, R_retlist;

     PROTECT(R_retlist = allocVector(VECSXP, npoints));
     PROTECT(R_nT = getAttrib(R_data, R_NamesSymbol));

     double *mean = 0,
    		*sigma2 = 0,
    		*points = REAL(R_points);

     if(glkm.krigType) {
    	 const char *nms[] = {"mean", "sigma2", "weights", ""};
    	 for(; i<npoints; i++, points++) {
             PROTECT(R_mean = allocVector(REALSXP, nCov));
        	 PROTECT(R_sigma2 = allocVector(REALSXP, nCov));
             PROTECT(R_weights_tmp = allocVector(VECSXP, nCov));

             mean = REAL(R_mean);
             sigma2 = REAL(R_sigma2);
             for(k=0; k < nCov; ++k) {
                 PROTECT(R_weights = allocVector(REALSXP, lx));
                 glkm.km[k]->univarKriging(points,npoints,mean+k,sigma2+k,REAL(R_weights),VARIANCE,info);
                 if(info != NO_ERROR){
                   XWRR(info,"univarKriging")
                 }
                 SET_VECTOR_ELT(R_weights_tmp, k, R_weights);
                 UNPROTECT(1);
             }
             setAttrib(R_mean,R_NamesSymbol,R_nT);
             setAttrib(R_sigma2,R_NamesSymbol,R_nT);
             setAttrib(R_weights_tmp,R_NamesSymbol,R_nT);

             PROTECT(R_tmp = mkNamed(VECSXP, nms));
             SET_VECTOR_ELT(R_tmp, 0, R_mean);
             SET_VECTOR_ELT(R_tmp, 1, R_sigma2);
             SET_VECTOR_ELT(R_tmp, 2, R_weights_tmp);

             SET_VECTOR_ELT(R_retlist, i, R_tmp);
             UNPROTECT(4);
    	 }
     } else {
    	  const char *nms[] = {"mean", ""};
    	  for(; i<npoints; i++, points++) {
    		 PROTECT(R_mean = allocVector(REALSXP, nCov));

    		 mean = REAL(R_mean);
             for(k=0; k < nCov; k++){
               	 glkm.km[k]->dualKriging(points,npoints,mean+k,info);
               	 if(info != NO_ERROR)
               	   XWRR(info,"dualKriging")
             }
             setAttrib(R_mean,R_NamesSymbol,R_nT);
             PROTECT(R_tmp = mkNamed(VECSXP, nms));
             SET_VECTOR_ELT(R_tmp, 0, R_mean);
             SET_VECTOR_ELT(R_retlist, i, R_tmp);
             UNPROTECT(2);
    	  }
    	  PROTECT(R_weights_tmp = allocVector(VECSXP, nCov));
    	  for(k=0; k < nCov; k++){
    		  PROTECT(R_weights = allocVector(REALSXP, lx));
    		  MEMCPY(REAL(R_weights),glkm.km[k]->dw,lx);
    		  SET_VECTOR_ELT(R_weights_tmp, k, R_weights);
    		  UNPROTECT(1);
    	  }
    	  setAttrib(R_weights_tmp,R_NamesSymbol,R_nT);
    	  setAttrib(R_retlist,install("weights"),R_weights_tmp);
    	  UNPROTECT(1);
     }

     SET_CLASS_NAME(R_retlist,"krigResult");

     UNPROTECT(2);
     return R_retlist;
}


SEXP initQL(SEXP R_qsd, SEXP R_qlopts, SEXP R_cvT)
{
  if(qlm_global) {
	 WRR("QL model is already initialized. Finalize now.");
	 if(!finalizeQL())
	   ERR("Could not free memory of `qlm_model`.")
  }
  if( (qlm_global = new(std::nothrow) ql_model_s(R_qsd, R_qlopts, R_cvT)) == NULL)
  	MEM_ERR(1,ql_model_s);
  return ScalarLogical(TRUE);
}

/**
 * @brief  FREE 'qldata' storage from R level
 * @return NULL
 */
SEXP finalizeQL() {
  if(qlm_global)
	DELETE(qlm_global)
  return ScalarLogical(TRUE);
}


SEXP qDValue(SEXP R_point) {
  if(!qlm_global)
    ERR("Pointer to `qldata` object not set (NULL).");
  return ScalarReal(qlm_global->intern_qfScoreStat(REAL(AS_NUMERIC(R_point))));
}

SEXP internQD(SEXP R_point) {
  if(!qlm_global)
    ERR("Pointer to `qldata` object not set (NULL).");

  int info = 0,
		dx = qlm_global->dx,
	  nCov = qlm_global->nCov;

  SEXP R_sig2, R_Iobs, R_S, R_jac, R_I, R_varS, R_stats, R_Vmat;
  PROTECT(R_S = allocVector(REALSXP,dx));
  PROTECT(R_sig2 = allocVector(REALSXP,nCov));
  PROTECT(R_stats = allocVector(REALSXP,nCov));
  PROTECT(R_I = allocMatrix(REALSXP,dx,dx));
  PROTECT(R_jac = allocMatrix(REALSXP,dx,nCov));
  PROTECT(R_varS = allocMatrix(REALSXP,dx,dx));
  PROTECT(R_Iobs = allocMatrix(REALSXP,dx,dx));
  PROTECT(R_Vmat = allocMatrix(REALSXP,nCov,nCov));

  double *x = REAL(R_point);
  double fval = qlm_global->qfScoreStat(x,REAL(R_jac),REAL(R_S),REAL(R_I),REAL(R_varS));

  if( (info = qlm_global->intern_quasiObs(x,REAL(R_S),REAL(R_Iobs))) != NO_ERROR)
   XWRR(info, "intern_quasiObs")

  /* copy prediction variance */
  MEMCPY(REAL(R_sig2),qlm_global->glkm->krigr[0]->sig2,nCov);
  MEMCPY(REAL(R_stats),qlm_global->glkm->krigr[0]->mean,nCov);

  const char *nms[] = {"value", "par", "I", "score", "sig2", "stats", "jac","varS", "Iobs", ""};
  SEXP R_ans = R_NilValue;
  PROTECT(R_ans = mkNamed(VECSXP, nms));
  SET_VECTOR_ELT(R_ans, 0, ScalarReal(fval));
  SET_VECTOR_ELT(R_ans, 1, R_point);
  SET_VECTOR_ELT(R_ans, 2, R_I);
  SET_VECTOR_ELT(R_ans, 3, R_S);
  SET_VECTOR_ELT(R_ans, 4, R_sig2);
  SET_VECTOR_ELT(R_ans, 5, R_stats);
  SET_VECTOR_ELT(R_ans, 6, R_jac);
  SET_VECTOR_ELT(R_ans, 7, R_varS);
  SET_VECTOR_ELT(R_ans, 8, R_Iobs);

  MEMCPY(REAL(R_Vmat),qlm_global->qld->vmat,qlm_global->nCov2);
  setAttrib(R_ans, install("Sigma"), R_Vmat);

  UNPROTECT(9);
  return R_ans;
}

SEXP internMD(SEXP R_point) {
  if(!qlm_global)
     ERR("Pointer to `qldata` object not set (NULL).");

  int info = 0,
		dx = qlm_global->dx,
	  nCov = qlm_global->nCov;

  SEXP R_sig2, R_Iobs, R_S, R_jac, R_I, R_varS, R_stats, R_Vmat;
  PROTECT(R_S = allocVector(REALSXP,dx));
  PROTECT(R_sig2 = allocVector(REALSXP,nCov));
  PROTECT(R_stats = allocVector(REALSXP,nCov));
  PROTECT(R_I = allocMatrix(REALSXP,dx,dx));
  PROTECT(R_jac = allocMatrix(REALSXP,dx,nCov));
  PROTECT(R_varS = allocMatrix(REALSXP,dx,dx));
  PROTECT(R_Iobs = allocMatrix(REALSXP,dx,dx));
  PROTECT(R_Vmat = allocMatrix(REALSXP,nCov,nCov));

  double *x = REAL(R_point);
  double fval = qlm_global->mahalDist(x,REAL(R_jac),REAL(R_S));

  if( (info = qlm_global->intern_varScore(x,REAL(R_jac),REAL(R_I),REAL(R_varS))) != NO_ERROR)
	XWRR(info,"intern_varScore")

  if( (info = qlm_global->intern_quasiObs(x,REAL(R_S),REAL(R_Iobs))) != NO_ERROR)
 	XWRR(info,"intern_quasiObs")

	/* copy prediction variance */
  MEMCPY(REAL(R_sig2),qlm_global->glkm->krigr[0]->sig2,nCov);
  MEMCPY(REAL(R_stats),qlm_global->glkm->krigr[0]->mean,nCov);

  const char *nms[] = {"value", "par", "I", "score", "sig2", "stats", "jac", "varS", "Iobs", ""};
  SEXP R_ans = R_NilValue;
  PROTECT(R_ans = mkNamed(VECSXP, nms));
  SET_VECTOR_ELT(R_ans, 0, ScalarReal(fval));
  SET_VECTOR_ELT(R_ans, 1, R_point);
  SET_VECTOR_ELT(R_ans, 2, R_I);
  SET_VECTOR_ELT(R_ans, 3, R_S);
  SET_VECTOR_ELT(R_ans, 4, R_sig2);
  SET_VECTOR_ELT(R_ans, 5, R_stats);
  SET_VECTOR_ELT(R_ans, 6, R_jac);
  SET_VECTOR_ELT(R_ans, 7, R_varS);
  SET_VECTOR_ELT(R_ans, 8, R_Iobs);

  MEMCPY(REAL(R_Vmat),qlm_global->qld->vmat,qlm_global->nCov2);
  setAttrib(R_ans, install("Sigma"), R_Vmat);
  UNPROTECT(9);
  return R_ans;
}

/**
 * \brief
 * 		Calculate Mahalanobis distance of the statistics and its gradient.
 *
 * 		If R_Sigma is used as inverse covariance matrix of statistics
 * 		(simulated), then the score vector as the gradient is returned,
 * 		otherwise the covariance matrix is continuously computed at each point
 * 		according to the settings stored in `qlopts`. The quasi-score as a gradient
 * 		makes sense only in case of constant variance matrix 'Sigma^{-1}' or if an
 * 		average approximation without kriging variances is used, that is, for kriging
 * 		type 'krigType'=="dual'.
 */
SEXP mahalValue(SEXP R_point) {
  if(!qlm_global)
     ERR("Pointer to `qldata` object not set (NULL).");
  GLKM glkm = qlm_global->glkm;
  ql_data qld = qlm_global->qld;

  SEXP R_score = R_NilValue;
  PROTECT(R_score = allocVector(REALSXP,glkm->dx));
  double f = qlm_global->mahalDist(REAL(AS_NUMERIC(R_point)),qld->jactmp,REAL(R_score));
  if(qlm_global->info != NO_ERROR)
    XERR(qlm_global->info,"mahalValue")

  SEXP R_ans = R_NilValue;
  PROTECT(R_ans = ScalarReal(f));
  setAttrib(R_ans,install("score"),R_score);
  UNPROTECT(2);
  return R_ans;
}

/* add variance matrix to result list as
 * an attribute if kriging is used for
 * variance matrix approximation
*/
void setVmatAttrib(ql_model qlm, SEXP R_VmatNames, SEXP R_ans) {
    SEXP R_Vmat = R_NilValue;
    PROTECT(R_Vmat = allocMatrix(REALSXP,qlm->nCov,qlm->nCov));
    MEMCPY(REAL(R_Vmat),qlm->qld->vmat,qlm->nCov2);

    setAttrib(R_Vmat, R_DimNamesSymbol, R_VmatNames);
    setAttrib(R_ans, install("Sigma"), R_Vmat);
    UNPROTECT(1);
}


SEXP mahalanobis(SEXP R_points, SEXP R_qsd, SEXP R_qlopts, SEXP R_cvT, SEXP R_qdValue, SEXP R_w)
{
	int i = 0, info = 0, np = LENGTH(R_points), nqsd = LENGTH(R_qsd);
	value_type type = (value_type) asInteger(AS_INTEGER(R_qdValue));

	/* init model storage */
	ql_model_t qlm(VECTOR_ELT(R_qsd,0),R_qlopts,R_cvT);
	/* append CV models as linked list */
	ql_model p = &qlm;
	for(int i=1; i<nqsd; i++){
		if((p->cvnext = new(std::nothrow)
		 ql_model_t(VECTOR_ELT(R_qsd,i),R_qlopts,R_NilValue)) == NULL)
		{
		  MEM_ERR(1,ql_model_s);
		}
		p = p->cvnext;
	}
	qlm.cvlast = p;
	qlm.nCV = nqsd-1;

	GLKM glkm = qlm.glkm;
	int dx = glkm->dx,
	  nCov = glkm->nCov;

	if(type > COPY_ZERO) {														/* only return scalar values */
 		  SEXP Rval = R_NilValue;
 		  PROTECT(Rval= allocVector(REALSXP,np));
 		  double *fx = REAL(Rval);
 		  for(i=0; i < np; ++i)
 			fx[i] = qlm.intern_mahalanobis(REAL(AS_NUMERIC(VECTOR_ELT(R_points,i))));
		  if (qlm.info != NO_ERROR)
		    PRINT_MSG("At least one of Mahalanobis distance computations produced errors.")
		  UNPROTECT(1);
		  return Rval;
	} else {

		  SEXP R_S, R_jac, R_I, R_ans, R_varS, R_sig2, R_stats, R_Iobs;

		  /* results */
		  SEXP R_ret = R_NilValue;
		  PROTECT(R_ret = allocVector(VECSXP,np));

		  /* names QI */
		  SEXP R_nI = VECTOR_ELT(R_points,0);
		  SEXP R_dimnames = R_NilValue;
		  PROTECT(R_dimnames = allocVector(VECSXP,2));
		  SET_DIMNAMES_MATRIX(R_dimnames,R_nI)

		  /* names Jacobian */
		  SEXP R_dimT = R_NilValue;
		  PROTECT(R_dimT = allocVector(VECSXP,2));
		  SEXP R_nT = getListElement( VECTOR_ELT(R_qsd,0), "obs" );
		  SET_DIMNAMES_MATRIX2(R_dimT,R_nI,R_nT)

		  // names variance matrix
 		  SEXP R_VmatNames = R_NilValue;
		  PROTECT(R_VmatNames = allocVector(VECSXP,2));
		  SET_DIMNAMES_MATRIX(R_VmatNames,R_nT)

	 	  int nprotect = 5;
		  int &info = qlm.info;
		  double fval = 0;
		  const char *nms[] = {"value", "par", "I", "score", "sig2", "stats", "jac", "varS", "Iobs", ""};

		  for(i=0; i<np; i++)
		  {
			 PROTECT(R_S = allocVector(REALSXP,dx));
			 PROTECT(R_sig2 = allocVector(REALSXP,nCov));
			 PROTECT(R_stats = allocVector(REALSXP,nCov));
			 PROTECT(R_I = allocMatrix(REALSXP,dx,dx));
			 PROTECT(R_varS = allocMatrix(REALSXP,dx,dx));
			 PROTECT(R_jac = allocMatrix(REALSXP,dx,nCov));
			 PROTECT(R_Iobs = allocMatrix(REALSXP,dx,dx));

			 double *x = REAL(AS_NUMERIC(VECTOR_ELT(R_points,i)));
			 /* mahalanobis distance */
			 fval = qlm.mahalDist(x,REAL(R_jac),REAL(R_S));
			 CHECK_UNPROTECT("mahalDist")

			 info = qlm.intern_varScore(x,REAL(R_jac),REAL(R_I),REAL(R_varS));
			 CHECK_UNPROTECT("intern_varScore")

			 info = qlm.intern_quasiObs(x,REAL(R_S),REAL(R_Iobs));
			 CHECK_UNPROTECT("intern_quasiObs")

			 /* copy prediction variance */
			 MEMCPY(REAL(R_sig2),qlm.glkm->krigr[0]->sig2,nCov);
			 MEMCPY(REAL(R_stats),qlm.glkm->krigr[0]->mean,nCov);

			 /*  set names but not for varS! */
			 setAttrib(R_jac, R_DimNamesSymbol, R_dimT);
			 setAttrib(R_I, R_DimNamesSymbol, R_dimnames);
			 setAttrib(R_varS, R_DimNamesSymbol, R_dimnames);
			 setAttrib(R_Iobs, R_DimNamesSymbol, R_dimnames);
			 setAttrib(R_sig2,R_NamesSymbol,getAttrib(R_nT,R_NamesSymbol));
			 setAttrib(R_stats,R_NamesSymbol,getAttrib(R_nT,R_NamesSymbol));

			 PROTECT(R_ans = mkNamed(VECSXP, nms));
			 SET_VECTOR_ELT(R_ans, 0, ScalarReal(fval));
			 SET_VECTOR_ELT(R_ans, 1, VECTOR_ELT(R_points,i));
			 SET_VECTOR_ELT(R_ans, 2, R_I);
			 SET_VECTOR_ELT(R_ans, 3, R_S);
			 SET_VECTOR_ELT(R_ans, 4, R_sig2);
			 SET_VECTOR_ELT(R_ans, 5, R_stats);
			 SET_VECTOR_ELT(R_ans, 6, R_jac);
			 SET_VECTOR_ELT(R_ans, 7, R_varS);
			 SET_VECTOR_ELT(R_ans, 8, R_Iobs);

			 setVmatAttrib(&qlm, R_VmatNames, R_ans);
			 SET_VECTOR_ELT(R_ret, i, R_ans);
			 UNPROTECT(8);
		  }

		  UNPROTECT(4);
		  return R_ret;
	}
}

double ql_model_s::mahalDist(double *x, double *jac, double *score) {
	/* kriging */
	  if( (info = glkm->intern_kriging(x)) != NO_ERROR){
		LOG_ERROR(info,"intern_kriging")
		return R_NaN;
	  }
	  krig_result krig = glkm->krigr[0];
	  if ( (info = glkm->intern_jacobian(x,jac,qld->fdwork)) != NO_ERROR){
	  	LOG_ERROR(info,"intern_jacobian")
		return R_NaN;
	  }

	  if(qld->qlopts.varType == CONST) {
		  /* constant 'Sigma':without added kriging prediction variances */
		  for(int k = 0; k < nCov; ++k)
		  	qld->tmp[k] = qld->qtheta[k] = qld->obs[k] - krig->mean[k];

		  /* vmat is already inverted at R level */
		  matmult(qld->vmat,nCov,nCov,qld->qtheta,nCov,ONE_ELEMENT,qld->tmp,info);
		  if(info != NO_ERROR ){
		    LOG_WARNING(info,"matmult")
			return R_NaN;
		  }
		  /* score vector as gradient of Mahalanobis distance */
		  mat_trans(qld->jactmp,nCov,jac,dx,dx,nCov,info);
		  if(info > 0){
		    WRR("`NaN` detected in `mat_trans`.")
			return R_NaN;
		  }
		  /* quasi-information */
		  //matmult(qld->vmat,nCov,nCov,qld->jactmp,nCov,dx,qld->Atmp,info);
		  //matmult(jac,dx,nCov,qld->Atmp,nCov,dx,varS,info);

		  /* quasi-score */
		  matmult(jac,dx,nCov,qld->tmp,nCov,ONE_ELEMENT,score,info);
	  	  if(info != NO_ERROR){
	  	    WRR("`NaN` detected in `matmult`.")
 		    return R_NaN;
	  	  }

	  } else {

		  /* variance matrix if kriging the variance matrix switched on */
		  varMatrix(x,qld->vmat,info);
		  if(info != NO_ERROR){
		  	LOG_ERROR(info,"varMatrix")
			return R_NaN;
		  }
		  /* quasi-score */
		  if( (info = intern_quasiScore(jac,score)) != NO_ERROR){
		  	LOG_ERROR(info,"intern_quasi-score")
			return R_NaN;
		  }
		  for(int k = 0; k < nCov; ++k)
	  	  	qld->tmp[k] = qld->qtheta[k];

	  	  /* use restored diagonal terms plus kriging variances
	  	   * for average approximation of variance matrix */
	  	  if(qld->qlopts.varType == MEAN) {
	  	  	 info = addVar(glkm->krigr[0]->sig2,nCov,qld->vmat_work,qld->work);
	  	  } else {
	  		  /* add CV error terms */
	  		  if(qld->qlopts.useCV) {
				 cvmod->cvError(x,glkm->km,glkm->krigr[0]->sig2,info);
				 if(info != NO_ERROR)
				  LOG_ERROR(info,"cvError")
			 }
	  		 MEMCPY(qld->vmat_work,qld->vmat,nCov2);
	  		 info = add2diag(qld->vmat_work,nCov,glkm->krigr[0]->sig2);
	  	  }
	  	  if(info != NO_ERROR){
	  		 WRR("`NaN` detected in `add2diag` adding kriging variances to variance matrix approximation.")
	  		 return R_NaN;
	  	  }
		  /*  Note that 'vmat_work' from 'intern_varScore'
		   *  contains additional kriging variances as diagonal terms */
		  gsiSolve(qld->vmat_work,nCov,qld->tmp,ONE_ELEMENT,qlsolve.vmat,info,Chol);
		  if(info != NO_ERROR){
			LOG_ERROR(info,"gsiSolve")
			return R_NaN;
		  }
	  }

	  double sum = 0;
	  for(int k = 0; k < nCov; ++k)
		 sum += qld->tmp[k] * qld->qtheta[k];

	  if( !R_FINITE(sum))
		 LOG_WARNING(NaN_ERROR,"mahalDist")
	  return sum;
}

/**
 *  \brief A test implementation for
 *  	   Cross-validation based prediction errors
 *
 */
//SEXP cvError(SEXP R_point, SEXP R_Xmat, SEXP R_data, SEXP R_covT, SEXP R_krigType, SEXP R_cm)
//{
//	int err=0;
//	/* CV models (might have different number of left out points! */
//	if(!isMatrix(R_Xmat))
//	 ERR("Expected sample matrix.");
//
//	cv_model_s cvmod(R_Xmat,R_data,R_cm);
//	/* full kriging models */
//	glkrig_models glkm(R_covT,R_Xmat,R_data,R_krigType,FALSE);
//
//	SEXP R_ret = PROTECT(allocVector(REALSXP,glkm.nCov));
//	cvmod.cvError(REAL(AS_NUMERIC(R_point)),glkm.km,REAL(R_ret),err);
//
//	UNPROTECT(1);
//	return R_ret;
//}


SEXP quasiDeviance(SEXP R_points, SEXP R_qsd, SEXP R_qlopts, SEXP R_cvT, SEXP R_qdValue, SEXP R_w)
{
      int i = 0, info = 0, np = LENGTH(R_points), nqsd = LENGTH(R_qsd);
      /* type of return value */
      value_type type = (value_type) asInteger(AS_INTEGER(R_qdValue));
      /* init QL model (including any type of CV models */
      ql_model_t qlm(R_qsd, R_qlopts, R_cvT);

      if(type > COPY_ZERO) {							// 1
    	  SEXP Rval = R_NilValue;
    	  if(type == COPY_ONE) {
    		  PROTECT(Rval = allocVector(REALSXP,np));
    		  double *fx = REAL(Rval);
    		  for(; i < np; ++i)
    		     fx[i] = qlm.intern_qfScoreStat(REAL(AS_NUMERIC(VECTOR_ELT(R_points,i))));
    	  } else if(type == COPY_MOD) {					// 2
    		  if(qlm.cvnext) {
    			  PROTECT(Rval = allocMatrix(REALSXP,np,qlm.nCV+1));
    			  double *fx = REAL(Rval);
    			  qlm.qDPressVec(R_points,fx);
    		  } else {
    			  PRINT_MSG("CV models are not initialized.")
    			  XERR(info,"quasiDeviance")
    		  }
    		  //for(;i < np; ++i)
			  // fx[i] = qlm.intern_wlogdet(REAL(AS_NUMERIC(VECTOR_ELT(R_points,i))),REAL(R_w)[0]);
    	  } else if(type == COPY_HIGHER) {				// 3
    		  PROTECT(Rval = allocVector(REALSXP,np));
    		  double *fx = REAL(Rval);
    		  for(;i < np; ++i)
    		 	fx[i] = qlm.trVarScore(REAL(AS_NUMERIC(VECTOR_ELT(R_points,i))));
    	  } else {
    		  PRINT_MSG("Unknown type of return value.")
    		  XERR(info,"quasiDeviance")
    	  }
          if(qlm.info != NO_ERROR)
        	PRINT_MSG("Quasi-deviance computations produced errors.")
		  UNPROTECT(1);
          return Rval;

      } else {

    	  int dx = qlm.dx,
    		  nCov = qlm.nCov;

    	  double *x = NULL, fval = 0;

    	  SEXP R_ans, R_sig2, R_Iobs, R_S, R_jac, R_I, R_varS, R_stats;
          SEXP R_ret = R_NilValue;
          PROTECT(R_ret = allocVector(VECSXP,np));

          // names QI
          SEXP R_nI = VECTOR_ELT(R_points,0);
          SEXP R_dimnames = R_NilValue;
          PROTECT(R_dimnames = allocVector(VECSXP,2));
          SET_DIMNAMES_MATRIX(R_dimnames,R_nI)

          // names Jacobian
          SEXP R_dimT = R_NilValue;
          SEXP R_nT = getListElement( VECTOR_ELT(R_qsd,0), "obs" );
          PROTECT(R_dimT = allocVector(VECSXP,2));
          SET_DIMNAMES_MATRIX2(R_dimT,R_nI,R_nT)

          // names variance matrix
          SEXP R_VmatNames = R_NilValue;
          PROTECT(R_VmatNames = allocVector(VECSXP,2));
          SET_DIMNAMES_MATRIX(R_VmatNames,R_nT)

          int nprotect = 6;
          const char *nms[] = {"value", "par", "I", "score", "sig2", "stats", "jac","varS", "Iobs", ""};

          for(; i < np; ++i)
          {
			  PROTECT(R_S = allocVector(REALSXP,dx));
			  PROTECT(R_sig2 = allocVector(REALSXP,nCov));
			  PROTECT(R_stats = allocVector(REALSXP,nCov));
			  PROTECT(R_I = allocMatrix(REALSXP,dx,dx));
			  PROTECT(R_jac = allocMatrix(REALSXP,dx,nCov));
			  PROTECT(R_varS = allocMatrix(REALSXP,dx,dx));
			  PROTECT(R_Iobs = allocMatrix(REALSXP,dx,dx));

			  x = REAL(AS_NUMERIC(VECTOR_ELT(R_points,i)));
			  fval = qlm.qfScoreStat(x,REAL(R_jac),REAL(R_S),REAL(R_I),REAL(R_varS));

			  info = qlm.intern_quasiObs(x,REAL(R_S),REAL(R_Iobs));
			  CHECK_UNPROTECT("intern_quasiObs")

			  /* copy kriging mena and iits prediction variance */
			  MEMCPY(REAL(R_sig2),qlm.glkm->krigr[0]->sig2,nCov);
			  MEMCPY(REAL(R_stats),qlm.glkm->krigr[0]->mean,nCov);

			  setAttrib(R_jac, R_DimNamesSymbol, R_dimT);
			  setAttrib(R_I, R_DimNamesSymbol, R_dimnames);
			  setAttrib(R_varS, R_DimNamesSymbol, R_dimnames);
			  setAttrib(R_Iobs, R_DimNamesSymbol, R_dimnames);
			  setAttrib(R_sig2,R_NamesSymbol,getAttrib(R_nT,R_NamesSymbol));
			  setAttrib(R_stats,R_NamesSymbol,getAttrib(R_nT,R_NamesSymbol));

			  PROTECT(R_ans = mkNamed(VECSXP, nms));
			  SET_VECTOR_ELT(R_ans, 0, ScalarReal(fval));
			  SET_VECTOR_ELT(R_ans, 1, VECTOR_ELT(R_points,i));
			  SET_VECTOR_ELT(R_ans, 2, R_I);
			  SET_VECTOR_ELT(R_ans, 3, R_S);
			  SET_VECTOR_ELT(R_ans, 4, R_sig2);
			  SET_VECTOR_ELT(R_ans, 5, R_stats);
			  SET_VECTOR_ELT(R_ans, 6, R_jac);
			  SET_VECTOR_ELT(R_ans, 7, R_varS);
			  SET_VECTOR_ELT(R_ans, 8, R_Iobs);

			  setVmatAttrib(&qlm, R_VmatNames, R_ans);
			  SET_VECTOR_ELT(R_ret, i, R_ans);
			  UNPROTECT(8);
          }

          UNPROTECT(4);
          return R_ret;
      }
}

/* Prepare computation of quasi-score statistic and norm of quasi-score */
int ql_model_s::qfScore(double *x, double *jac, double *score, double *qimat, double *varS) {
	/* kriging statistics */
	if ( (info = glkm->intern_kriging(x)) != NO_ERROR){
		LOG_ERROR(info,"intern_kriging")
		return info;
	}
	//printVector("krig.mean",glkm->krigr[0]->mean,&dx);
	//printVector("krig.var",glkm->krigr[0]->sig2,&dx);
	if ( (info = glkm->intern_jacobian(x,jac,qld->fdwork)) != NO_ERROR){
		LOG_ERROR(info,"intern_jacobian")
		return info;
	}
	/* compute variance matrix approximation */
	varMatrix(x,qld->vmat,info);
	if(info != NO_ERROR) {
		LOG_ERROR(info,"varMatrix")
		return info;
	}
	//printMatrix("vmat",qld->vmat,&nCov,&nCov);
	if( (info = intern_quasiScore(jac,score)) != NO_ERROR){
		 LOG_ERROR(info,"intern_quasi-score")
		 return info;
	}
	// modified quasi-information
	if ( (info = intern_varScore(x,jac,qimat,varS)) != NO_ERROR){
		 LOG_ERROR(info,"intern_varScore")
		 return info;
	}
	//printMatrix("varS",varS,&dx,&dx);
	return NO_ERROR;
}

/**
 *  \brief Modified quasi-information matrix, no quasi-score!
 */
int ql_model_s::varScore(double *x, double *varS) {
	    if((info = glkm->intern_kriging(x)) != NO_ERROR){
			LOG_ERROR(info,"intern_kriging")
			return info;
		}
		if( (info = glkm->intern_jacobian(x,jac,qld->fdwork)) != NO_ERROR){
			LOG_ERROR(info,"intern_jacobian")
			return info;
		}
		varMatrix(x,qld->vmat,info);
		if(info != NO_ERROR) {
			LOG_ERROR(info,"varMatrix")
			return info;
		}
		// compute Atmp
		mat_trans(qld->Atmp,nCov,jac,dx,dx,nCov,info);
		if(info > 0){
		  WRR("`NaN` detected in `mat_trans` computing 'Atmp'.")
		  return info;
		}
		gsiSolve(qld->vmat,nCov,qld->Atmp,dx,qlsolve.vmat,info,Chol);
		if(info != NO_ERROR){
		  LOG_ERROR(info,"gsiSolve");
		  return info;
		}
		// modified quasi-information
		if((info = intern_varScore(x,jac,qimat,varS)) != NO_ERROR){
			 LOG_ERROR(info,"intern_varScore")
			 return info;
		}
		return NO_ERROR;
}

double ql_model_s::trVarScore(double *x) {
	if((info = varScore(x,varS)) != NO_ERROR){
		LOG_ERROR(info,"varScore")
		return R_NaN;
	}
	double sum=0;
	for(int i=0; i < dx; ++i)
	  sum += varS[i*dx+i];
	if(!R_FINITE(sum)){
	  WRR("`NaN` detected in `trVarScore`.")
	  return R_NaN;
	}
	return sum;
}

/**
 * \brief Quasi-Fisher score statistic
 *
 */
double ql_model_s::qfScoreStat(double *x, double *jac, double *score, double *qimat, double *varS) {
 	if(info != (qfScore(x,jac,score,qimat,varS) != NO_ERROR)){
		LOG_ERROR(info,"qfScore");
		return R_NaN;
	}
 	return qfValue(score,qimat);
}

// TODO: does not work yet, produces wrong results!
void ql_model_s::qDPressVec(SEXP R_points, double *mpress) {
	double *x = 0, *pp = 0;
	int info = 0, np=LENGTH(R_points);
	for(int i=0;i < np; ++i, mpress += i) {
		x = REAL(AS_NUMERIC(VECTOR_ELT(R_points,i)));
		/* full model */
		*mpress = intern_qfScoreStat(x);
		/* vector of LOO-CV quasi-deviances */
		pp = mpress+np;
		for(ql_model p=cvnext; p != 0; p=p->cvnext, pp += np){
		 *pp = p->intern_qfScoreStat(x);
		 printVector("score: ",p->score,&p->dx);
		}
	}
}


double ql_model_s::intern_wlogdet(double *x, double w) {
	double val = intern_qfScoreStat(x);
	if(info != NO_ERROR) {
	  WRR("`NaN` detected in `intern_qfScoreStat`.")
	  return R_NaN;
	}
	double tmp = -(1.0-w)*val;
	//double tmp = -(1.0-w)*std::log(val);
	if(w > 0.0) {
		double sum = logdet(varS, dx, 1, info);
		if(info != NO_ERROR) {
			  WRR("`NaN` detected in `logdet`.")
			  LOG_ERROR(info,"intern_wlogdet")
			  return info;
		}
		//tmp += w/dx*sum;
		tmp += w*sum;
	}
	return tmp;
}

/**
 * Computation of modified quasi-deviance.
 * Use variance of quasi-score (varS) as weighting matrix instead of QI.
 */
double ql_model_s::qfValue(double *score, double *varS) {
	MEMCPY(qlsolve.score,score,dx);
	gsiSolve(varS,dx,qlsolve.score,ONE_ELEMENT,qlsolve.varS,info,Chol);
	if(info != NO_ERROR){
	   LOG_ERROR(info,"gsiSolve");
	   return R_NaN;
	}
	double sum=0.;
	for(int i=0; i<dx; ++i)
	  sum += qlsolve.score[i]*score[i];
	return sum;
}

/* IN:
 * 	jac  = dE[T(X)]/dTheta  ( dx (rows) * nCov (cols) )
 * 	Atmp = Sigma_{Theta}^{-1} %*% t(jac) = B
 *
 * OUT:
 * 	qimat = jac %*% B
 *
 */
int ql_model_s::intern_quasiInfo(double *jac, double *qimat) {
	matmult(jac,dx,nCov,qld->Atmp,nCov,dx,qimat,info);
	if(info != NO_ERROR){
	  WRR("`NaN` detected in `matmult`.")
	}
	return info;
}


/* IN:
 * 	x 	 = current point, only for CV error of quasi-score
 * 	jac  = dE[T(X)]/dTheta ( dx (rows) * nCov (cols) )
 * 	Atmp = Sigma_{Theta}^{-1} %*% t(jac) = B
 *
 * OUT:
 * 	qimat = t(B) %*% (Sigma_{Theta} + hat{Sigma}_k) %*% B
 *
 */
int ql_model_s::intern_varScore(double *x, double *jac, double *qimat, double *varS){
	if(qld->qlopts.varType == CONST) {
		matmult(jac,dx,nCov,qld->Atmp,nCov,dx,varS,info);
		MEMCPY(qimat,varS,dxdx);
	} else if(nCV > 0) {
		/* quasi-information without additional CV error terms */
		if( (info = intern_quasiInfo(jac,qimat)) != NO_ERROR){
			 LOG_ERROR(info,"intern_quasiInfo")
		     return info;
		}
		/* CV computation of varS */
		int i = 0, have_na =0;
		double *smean = qld->fdscore, *tmp = qlsolve.qimat, *stmp = qlsolve.score;

		for(ql_model p=cvnext; p != 0; p=p->cvnext) {
			wrap_intern_quasiScore(x, (void*) p, p->score, info);
			for(i=0;i<dx;i++)
			 smean[i] += p->score[i];
		}
		/* mean of quasi-score over CV quasi-scores*/
		for(i=0;i<dx;i++) smean[i] /= (double)nCV;
		/* init variance matrix of quasi-score */
		for(i=0;i<dxdx;i++) varS[i]=.0;

		/* cross-validation MSE of quasi-score */
		for(ql_model p=cvnext; p != 0; p=p->cvnext) {
			for(i=0;i<dx;i++)
			 stmp[i] = p->score[i]-smean[i];
			matmult_trans(stmp,ONE_ELEMENT,dx,stmp,ONE_ELEMENT,dx,tmp,info);
			if(info > 0){
		      WRR("`NaN` detected in `matmult_trans`.")
			  return info;
			}
			for(i=0;i<dxdx;i++) varS[i] += tmp[i];
		}
		/* add to quasi-information (varS[i]/nCV is MSE of quasi-score) */
		for(i=0;i<dxdx;i++) {
			varS[i] = qimat[i] + varS[i]/nCV;
			if(!R_FINITE(varS[i]))
			 { have_na=1; break;}
		}
		if(have_na>0)
		 WRR("`NaN` detected in summation of quasi-information matrix and cross-validation error of quasi-score.")
	} else {
		/* original quasi-information */
		matmult(jac,dx,nCov,qld->Atmp,nCov,dx,qimat,info);
		/* variance of quasi-score */
		if(!glkm->krigType) {
			MEMCPY(varS,qimat,dxdx);
		} else {
			/* CV error of statistics only */
			if(qld->qlopts.useCV) {
			 cvmod->cvError(x,glkm->km,glkm->krigr[0]->sig2,info);
			 if(info != NO_ERROR)
			  LOG_ERROR(info,"cvError")
		    }
			if(qld->qlopts.varType == MEAN) {
			 info = addVar(glkm->krigr[0]->sig2,nCov,qld->vmat_work,qld->work);
			} else {
			 /* and for kriging approximation of variance matrix add prediction variances temporarily */
			 MEMCPY(qld->vmat_work,qld->vmat,nCov2);
			 info = add2diag(qld->vmat_work,nCov,glkm->krigr[0]->sig2);
			}
			if(info != NO_ERROR){
			  WRR("`NaN` detected in `add2diag` adding kriging variances to variance matrix approximation.")
			  return info;
			}
			/* modified quasi-information*/
			matmult_trans(qld->Atmp,nCov,dx,qld->vmat_work,nCov,nCov,qld->jactmp,info);
			matmult(qld->jactmp,dx,nCov,qld->Atmp,nCov,dx,varS,info);
		}
	}
	if(info != NO_ERROR)
	 WRR("`NaN` detected in `matmult`.")
	return info;
}


/* computes Atmp=B (see above) first */
void ql_model_s::quasiScore(double *mean, double *jac, double *vmat, double *score, int &err) {
	for(int i=0; i<nCov; i++)
		qld->qtheta[i] = qld->obs[i]-mean[i];

	if(qld->qlopts.varType == CONST) {
		/* transpose jac */
		mat_trans(qld->jactmp,nCov,jac,dx,dx,nCov,info);
		if(info > 0){
		   WRR("`NaN` detected in `mat_trans`.")
		   return;
		}
		/* vmat already inverted at R level */
		matmult(vmat,nCov,nCov,qld->jactmp,nCov,dx,qld->Atmp,info);
		if(info != NO_ERROR){
			WRR("`NaN` detected in `matmult`.")
			return;
	    }
		/* quasi-score with constant variance matrix */
		matmult(qld->qtheta,ONE_ELEMENT,nCov,qld->Atmp,nCov,dx,score,err);
		if(err != NO_ERROR){
		  WRR("`NaN` detected in `matmult`.")
		  return;
		}

	} else {

		/* transpose jac */
		mat_trans(qld->Atmp,nCov,jac,dx,dx,nCov,info);
		if(info > 0){
		  WRR("`NaN` detected in `mat_trans`.")
		  return;
		}
		/* solve for inverted vmat  */
		gsiSolve(vmat,nCov,qld->Atmp,dx,qlsolve.vmat,info,Chol);
		if(info != NO_ERROR){
		  LOG_ERROR(info,"gsiSolve");
		  return;
		}
		/* quasi-score */
		matmult(qld->qtheta,ONE_ELEMENT,nCov,qld->Atmp,nCov,dx,score,err);
		if(err != NO_ERROR){
		  WRR("`NaN` detected in `matmult`.")
		  return;
		}

	}
}


/**  \brief Wrapper function: Score vector
 *
 * @param x     Parameter vector
 * @param data  Pointer data
 * @param score Score vector allocated
 */
void
wrap_intern_quasiScore(double *x, void *data, double *score, int &err ) {
   ql_model qlm = (ql_model) data;
   GLKM glkm = qlm->glkm;
   ql_data qld = qlm->qld;

   krig_result krig_tmp = glkm->krigr[1];
   glkm->kriging(x,krig_tmp->mean,krig_tmp->sig2,krig_tmp->w,err);
   if(err != NO_ERROR){
	   LOG_ERROR(err,"kriging");
	   return;
   }
   /* compute numerical derivatives (Jacobian of statistics) */
   glkm->jacobian(x,krig_tmp->mean,qld->fdjac,qld->fdwork,err);
   if(err != NO_ERROR){
	   LOG_ERROR(err,"jacobian");
	   return;
   }
   /* prediction error of mean values of statistics */
   if(qld->qlopts.varType != CONST) {
	 //  if(qld->qlopts.useCV){
     //    qlm->cvmod->cvError(x,glkm->km,krig_tmp->sig2,err);
     //    if(err != NO_ERROR)
     //      LOG_ERROR(err," cvError");
     //  }
     qlm->varMatrix(x,qld->vmat_work,err);
     if(err != NO_ERROR){
    	 LOG_ERROR(err,"varMatrix")
    	 return;
     }
   }
   qlm->quasiScore(krig_tmp->mean,qld->fdjac,qld->vmat_work,score,err);
   if(err != NO_ERROR)
	 LOG_ERROR(err,"wrap_intern_quasiScore");
   return;
}

/**  \brief Wrapper function: Quasi-deviance
 *
 * @param x     Parameter vector
 * @param data  Pointer data
 * @param score Score vector allocated
 */
double
wrap_intern_quasiDeviance(double *x, void *data, int &err ) {
   ql_model qlm = (ql_model) data;
   return qlm->intern_qfScoreStat(x);
}


/** \brief Wrapper function: intern_kriging
 *
 *     Comment: Kriging result is constructed and destroyed
 *
 * @param x point
 * @param data data pointer
 * @param mean kriging mean vector
 */
void
wrap_intern_kriging(double *x, void *data, double *mean, int &err) {
  GLKM glkm = (GLKM) data;
  krig_result krig_tmp = glkm->krigr[1];
  glkm->kriging(x,mean,krig_tmp->sig2,krig_tmp->w,err);
  if(err != NO_ERROR)
    LOG_ERROR(err,"kriging");
}

void
ql_model_s::varMatrix(double *x, double *vmat, int &err) {
	if(qld->qlopts.varType > MEAN) {
	   if(varkm == NULL)
		 ERR("Null pointer exception in `varMatrix`. This seems to be a severe bug.");
	   if( (err = varkm->intern_kriging(x)) != NO_ERROR){
		    LOG_ERROR(err,"intern_kriging");
		  	return;
       }
	   if(qld->qlopts.varType == LOGKRIG) {
		   err = exp_mergeMatrix(varkm->krigr[0]->mean,vmat,nCov,qld->workx);
	   } else {
		   err = chol2var(varkm->krigr[0]->mean,vmat,nCov,qld->workx);
	   }
	   if(err != NO_ERROR) {
		  LOG_ERROR(err,"chol2var");
		  return;
	   }
	}
#ifdef DEBUG
	    printVector("varkm->mean",varkm->krigr[0]->mean,&nCov);
	    printMatrix("vmat",vmat,&nCov,&nCov);
#endif
}

int ql_model_s::intern_quasiObs(double *x, double *score, double *qiobs) {
   fdJacobian(x,dx,score,dx,qiobs,qld->fdscore,&wrap_intern_quasiScore,(void*) this,FD_EPS,ONE_ELEMENT,info);
   if(info != NO_ERROR)
	 WRR("`NaN` values detected in `fdJac`.")
   return info;
}

/** @brief Allocate storage for a kriging model
 *
 * @param Xmat
 * @param lx
 * @param dx
 * @param trend
 * @param type
 * @return
 */

void
krig_model_s::alloc() {
  trend = cov.trend;
  switch (trend) {
    case 0: fddim = 1;
            break;
    case 1: fddim = dx+1;
            break;
    case 2: fddim = (dx+1)*(dx+2)/2;
            break;
    default: error(_("Invalid drift term specified.\n"));
             break;
  }
  CALLOCX(s0,lx,double);
  CALLOCX(f0,fddim,double);
  CALLOCX(Fmat,lx * fddim, double);
  CALLOCX(Cmat,lx * lx, double);

  /** allocate for inverse covariance matrix */
  if( (ks = new(std::nothrow) krig_storage_s(lx,fddim))==NULL)
    MEM_ERR(1,krig_storage_s);
  CALLOCX(Cinv,lx*lx, double);
  CALLOCX(X1mat,lx*fddim, double);
  CALLOCX(Qmat,fddim*fddim, double);
  CALLOCX(dw,lx+fddim, double);
  CALLOCX(dw,lx+fddim, double);
}

void
krig_model_s::setup(double *_data)
{
	  int info=0;
	  // allocation of matrices
	  alloc();
	  // copy data or not
	  CALLOCX(data,lx,double);
	  MEMCPY(data,_data,lx);

	  // trend matrix
	  Fmatrix(Xmat,Fmat,lx,dx,trend);
	  // REML covariance matrix
	  if( (info = intern_covMatrix(Xmat,dx,lx,Cmat,&cov)) != NO_ERROR){
	      PRINT_MSG("Failed to setup kriging model (including prediction variances).")
		  XERR(info,"intern_covMatrix")
	  }

	  /*  init matrices */
	  invMatrix(Cmat,lx,Cinv,info,Chol);
	  if(info != NO_ERROR){
	    PRINT_MSG("Setting up kriging model failed due to inversion error of the covariance matrix.")
		XERR(info,"invMatrix")
	  }
 	  /* store matrices for solving kriging equations */
	  matmult(Cinv,lx,lx,Fmat,lx,fddim,X1mat,info);
	  if(info > 0)
	   WRR("`NaN` detected in matrix multiplication.")
	  matmult_trans(Fmat,lx,fddim,X1mat,lx,fddim,Qmat,info);
	  if(info > 0)
	   WRR("`NaN` detected in matrix multiplication.")
	  /* also initialize dual kriging equations */
	  if( (info = intern_dualKrigingWeights(Cmat,lx,Fmat,fddim,data,dw)) != NO_ERROR){
	     PRINT_MSG("Failed to setup (dual) kriging model.")
		 XWRR(info,"intern_dualKrigingWeights")
	  }
}

void cv_model_s::set(SEXP R_Xmat, SEXP R_data, SEXP R_cm) {
	int i=0, j=0, k=0, l=0, m=0,
		lx=0, len=0;

	int *id=0, *dim=GET_DIMS(R_Xmat);

	double tmp=0;
	double *Xmatp=NULL,
		  **datap=NULL;

	np  = dim[0];
	dx  = dim[1];

	nc  = LENGTH(R_cm);
	if(nc == 0)
	  ERR("Length of covariance models for cross-validation is zero. There might be bug somewhere.")

	fnc = 1.0/(nc*(nc-1));
	Xmat = REAL(R_Xmat);
	nCov = LENGTH(VECTOR_ELT(VECTOR_ELT(R_cm,0),0)); // because each is in a list again
    //errType = !std::strcmp("max",translateChar(asChar(getAttrib(R_cm,install("type")))));

	/* container of CV models*/
	if( (cm = new(std::nothrow) GLKM[nc]) == NULL)
		MEM_ERR(1,GLKM);

	CALLOCX(s2,nCov, double);
	CALLOCX(ybar,nCov, double);
	CALLOCX(ytil,nc*nCov, double);
	CALLOCX(y0,nCov,double);

	SEXP R_covList, R_id;
	/* construct new CV model data */
	for(k = 0; k < nc; k++) {
		R_covList = VECTOR_ELT(VECTOR_ELT(R_cm,k),0);
		R_id = AS_INTEGER(getAttrib(R_covList,install("id")));
		len = LENGTH(R_id);
		id = INTEGER(R_id);
		lx = np - len;

		CALLOCX(datap,nCov,double*);
		CALLOCX(Xmatp,lx*dx,double);

		for(j = 0; j < nCov; j++)
		 CALLOCX(datap[j],lx,double);

		// points
		for(m = 0, i = 0; i < np; i++) {
			 Rboolean isin = FALSE;
			 for(l = 0; l < len; l++) {
				if(i == (id[l]-1)) {
				  isin = TRUE;
				  break;
				}
			 }
			 if(isin) continue;
			 for(j = 0; j < nCov; j++) {
			   tmp = REAL(AS_NUMERIC(VECTOR_ELT(AS_LIST(R_data),j+dx)))[i];
			   if (!R_FINITE(tmp) )
				 WRR("`NaN` detected in data vector.")
			   datap[j][m] = tmp;
			 }
			 for(j = 0; j < dx; j++)
			   Xmatp[lx*j+m] = Xmat[np*j+i];
			 ++m;
		}
		if( (cm[k] = new(std::nothrow) glkrig_models_s(R_covList,Xmatp,datap,lx,dx,VARIANCE)) == NULL)
			MEM_ERR(1,glkrig_models_s);

		// free memory
		for(j = 0; j < nCov; j++)
	      FREE(datap[j]);

		FREE(datap);
		FREE(Xmatp);
	}
}

void cv_model_s::cvError(double *x, krig_model *km, double *cv, int &info) {
	int i,k;
	double yhat = 0;

	for(k = 0; k < nCov; ++k) {
		s2[k] = ybar[k] = 0;
		km[k]->dualKriging(x,ONE_ELEMENT,y0+k,info);
		if(info != NO_ERROR)
		 break;
	}
	if(info != NO_ERROR){
	 LOG_ERROR(info,"dualKriging");
	 return;
	}

	GLKM cmp=NULL;
	krig_model *kmp=NULL;

	for(i = 0; i < nc; ++i) {
			cmp = cm[i];
			kmp = cmp->km;
			for(k = 0; k < nCov; ++k) {
				kmp[k]->dualKriging(x,ONE_ELEMENT,&yhat,info);
				if (!R_FINITE(yhat))
				 { info=1; continue; }
				ytil[k*nc+i] = nc*y0[k] - (nc-1)*yhat;
				ybar[k] += ytil[k*nc+i];
			}
	}

	for(i = 0; i < nc; ++i) {
		for(k = 0; k < nCov; ++k){
			if (!R_FINITE(ybar[k]))
			  { info=1; continue; }
			s2[k] += SQR(ytil[k*nc+i]-ybar[k]/nc);
		}
	}
    for(k = 0; k < nCov; ++k) {
	 if (!R_FINITE(s2[k]))
	  { info=1; continue; }
	 cv[k] = fnc*s2[k];
	}
}


inline void
glkrig_models::jacobian(double *x, double *mean, double *jac, double *fdwork, int &info) {
	fdJacobian(x,dx,mean,nCov,jac,fdwork,&wrap_intern_kriging,(void*)this,FD_EPS,ZERO_ELEMENT,info);
	if(info != NO_ERROR)
	  LOG_WARNING(info,"fdJac")
}

/**
 * Kriging a single point,
 * for all nCov statistics
 */
inline void
glkrig_models::kriging(double *x, double *m, double *s, double *w, int &info) {
	if(krigType) {
		for(int k=0; k<nCov; k++){
			km[k]->univarKriging(x,ONE_ELEMENT,m+k,s+k,w+k,krigType,info);
			if(info != NO_ERROR){
			  LOG_ERROR(info,"univarKriging")
			}
		}
	} else {
		for(int k=0; k<nCov; k++){
			km[k]->dualKriging(x,ONE_ELEMENT,m+k,info);
			if(info != NO_ERROR)
			  LOG_ERROR(info,"dualKriging")
		}
	}
}


/**
 * Dual kriging a single point, single statistic
 */
void
krig_model_s::dualKriging(double *x, int nx, double *m, int &info) {
	if( (info = intern_covVector(Xmat,dx,lx,x,nx,s0,&cov)) != NO_ERROR)
	  WRR("`NaN` detected in `intern_covVector`.")

	trendfunc(x,nx,dx,f0,cov.trend);

	/** weights do not change */
	if( (info = intern_dualKrigingPrediction(s0,f0,fddim,dw,lx,m)) != NO_ERROR)
	 WRR("`NaN` detected in `intern_dualKrigingPrediction`.")

}

/** \brief Internal function: Kriging
 * Chilès, J.-P. and Delfiner, P. (2008) Frontmatter, in Geostatistics:
 *              Modeling Spatial Uncertainty, John Wiley & Sons, Inc., Hoboken, NJ, USA.
 *              page: 169
 *
 * @param Xmat sample matrix
 * @param data observed data vector
 * @param Fmat trend matrix
 * @param fddim columns of trend matrix
 * @param Cinv inverse covariance matrix (stored)
 * @param point prediction point x0 as matrix
 * @param Npoints rows of point
 * @param mean predicted (kriging) mean
 * @param sigma2 prediction (kriging) variance
 * @param w kriging weights
 * @param cov covariance model
 */
void
krig_model_s::univarKriging(double *x, int nx, double *mean, double *sigma2, double *lambda, krig_type type, int &err) {

    int j=0, k=0, info=0;

	/**
	 *  solving Kriging equations
	 *    X1 = C^{-1} F
     *     Q = F^{t} X1
	 *
	 *    lambdak = C^-1 s0
	 *    R = F^t lambdak -f0
	 *    mu = Q^-1 R
	 *    lambda = lambdak - X1 mu
	 */

    /* storage */
    double *lambdak = ks->lambdak,
    	   *Qwork = ks->Qwork,
    	   *Rvec = ks->Rvec,
		   *mu = ks->mu;

	/*! calculates the covariance vector: sigma0 */
	if( (info = intern_covVector(Xmat,dx,lx,x,nx,s0,&cov)) != NO_ERROR )
	 WRR("`NaN` detected in `intern_covVector`.");

	/*! calculates the trend vector: f0 */
	trendfunc(x,nx,dx,f0,cov.trend);

	//printVector("inter_covVector",s0,lx);
	//printVector("inter_trendfunc",f0,fddim);
    //printMatrix("Qmat",Qmat,fddim,fddim);

	matmult(Cinv,lx,lx,s0,lx,ONE_ELEMENT,lambdak,info);
	if(info > 0)
	  WRR("`NaN` detected in matrix multiplication")

	matmult_trans(Fmat,lx,fddim,lambdak,lx,ONE_ELEMENT,Rvec,info);
	if(info > 0)
	  WRR("`NaN` detected in matrix multiplication")

	for(k=0; k<fddim; k++)
	  mu[k] = Rvec[k] - f0[k];

	/** solve Q mu = R */
	gsiSolve(Qmat,fddim,mu,ONE_ELEMENT,Qwork,err,Chol);
	if(info != 0){
	    *mean = R_NaN;
	  *sigma2 = R_NaN;
	  *lambda = R_NaN;
	  LOG_ERROR(err,"gsiSolve");
	  return;
	}
	matmult(X1mat, lx, fddim, ks->mu, fddim, ONE_ELEMENT, lambda, info);
	if(info > 0)
	  WRR("`NaN` detected in matrix multiplication")

	for(j=0; j < lx; j++)
	  lambda[j] = lambdak[j] - lambda[j];

	/* calculate kriging variance */
	double sum = cov.cf(&cov,&ZERO_DBL);
	//Rprintf("%s: %u: %f", __FILE__, __LINE__, sum);;

	for (j=0; j < lx; j++)
	  sum -= lambda[j] * s0[j];
	for (j=0; j < fddim; j++)
	  sum -= mu[j] * f0[j];

	/* for numerical stability */
	if (sum < 0) {
		sum = KRIGE_TOLERANCE;
	} else if(sum<KRIGE_TOLERANCE) {
		sum = 0.0;
	}
	*sigma2=sum;

	/* manipulate Kriging MSE */
	//*sigma2 += cov.nugget;
	// Rprintf("%s:%u: %f", __FILE__, __LINE__, *sigma2);

	/** kriging mean */
	if(type == VARIANCE) {
	 for (sum=j=0; j<lx; j++)
	    sum += lambda[j] * data[j];
	 *mean = sum;
	} else {
		ERR("Kriging type 'dual' is deprecated.")
	}
	err = info;
}

#endif /* KRIGING_CC_ */

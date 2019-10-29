/**
 * @file        kriging.h
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       R (Interface) Functions for kriging
 *
 *
 */

#ifndef KRIGING_H_
#define KRIGING_H_

#include "covariance.h"

/** typedefs */
typedef enum { DUAL = 0, VARIANCE = 1} krig_type;
typedef enum { CONST = 0, MEAN = 1, KRIG = 2, LOGKRIG = 3} var_type;
typedef enum { COPY_ZERO = 0, COPY_ONE = 1, COPY_MOD = 2, COPY_HIGHER = 3} value_type;

typedef struct krig_storage_s {
	double *Qwork,
	       *lambdak,
		   *Rvec,
		   *mu;

	krig_storage_s() :
		Qwork(0),lambdak(0),Rvec(0),mu(0)
	{}

	krig_storage_s(int lx, int fddim) :
			Qwork(0),lambdak(0),Rvec(0),mu(0)
	{
		CALLOCX(mu,fddim, double);
		CALLOCX(Rvec,fddim,double);
		CALLOCX(lambdak,lx, double);
		CALLOCX(Qwork,fddim*fddim,double);
	}


	~krig_storage_s(){
		FREE(mu)
		FREE(Rvec)
		FREE(Qwork)
		FREE(lambdak)
	}

} krig_storage_t, *krig_storage;

typedef struct krig_result_s {
  double *mean, *sig2, *w;
  krig_type krigType;

  krig_result_s() : mean(0), sig2(0), w(0), krigType(VARIANCE) {}

  krig_result_s(int nCov, int lx, krig_type type) :  mean(0), sig2(0), w(0), krigType(type)
  { alloc(nCov,lx); }

  ~krig_result_s() {
	 FREE(w)
	 FREE(mean)
	 FREE(sig2)
  }

  void alloc(int nCov, int lx) {
	  if(nCov==0 || lx==0)
	    ERR("number of covariance models or sample points must not be zero!");

	  CALLOCX(w,nCov*lx,double);     // kriging weights as matrix: nCov rows, lx columns
	  CALLOCX(mean,nCov,double);     // kriging means (at single point)
	  CALLOCX(sig2,nCov,double);     // kriging variances (at single point) and nCov covariances
  }

} krig_result_t, *krig_result;

typedef struct ql_options_s {
	Rboolean useCV;
	var_type varType;

  ql_options_s(SEXP R_opts)  {
	  useCV    = (Rboolean) asLogical(getListElement(R_opts, "useCV"));
	  SEXP R_varType = getListElement(R_opts, "varType");
	  if(!isString(R_varType))
	    Rf_error(_("Variance interpolation type should be a string: varType"));
	  const char *var_type = translateChar(asChar(R_varType));
	  if(!strcmp("const",var_type)) {
		  varType = CONST;
	  } else if( !strcmp("kriging",var_type) ) {
		  varType = KRIG;
	  } else  if( !strcmp("logKrig",var_type) ) {
		  varType = LOGKRIG;
	  } else varType = MEAN;
  }
} ql_options_t, *ql_options;


typedef struct krig_model_s {
  int lx,      /* amount of sampled points */
	  dx;      /* dimension */

  krig_type krigType;

  double *Xmat,	 /* design (always only pointer) */
  	  	 *data,  /* stats vector */
         *Fmat,  /* Trend matrix */
         *Cmat,  /* Covariance matrix */
         *Cinv,  /* Inverse Covariance */
         *X1mat, /* X1 matrix, solving kriging equations */
         *Qmat,  /* Q matrix, solving kriging equations */
	     *dw;    /* Dual kriging weights */

  double *f0,    /* trend vector*/
         *s0; 	 /* and covariance vector work */

  int   trend,   /* trend order */
        fddim;   /* columns of trend matrix F*/

  cov_model cov;
  krig_storage ks;

  krig_model_s(SEXP _R_cov, double *_Xmat, double *_data,
		  int _lx, int _dx, krig_type _type ) :
	 lx(_lx), dx(_dx), krigType(_type), Xmat(_Xmat),
	 data(0),Fmat(0),Cmat(0),Cinv(0),X1mat(0),Qmat(0),dw(0),f0(0),s0(0),trend(0),fddim(0),
	 cov(_R_cov), ks(0)
  {
	 setup(_data);
  }

  ~krig_model_s() {
	  FREE(Cinv)
	  FREE(X1mat)
	  FREE(Qmat)
	  FREE(dw)
	  FREE(dw)
	  FREE(s0)
	  FREE(f0)
	  FREE(Fmat)
  	  FREE(Cmat)
  	  FREE(data)
	  DELETE(ks)
  }

  void alloc();

  void setup(double *data);

  void dualKriging(double *x, int nx, double *m, int &err);

  void univarKriging(double *x, int nx, double *m, double *s, double *w, krig_type type, int &err);

} krig_model_t, *krig_model;


typedef struct glkrig_models_s
{
  double *Xmat;	  /* design matrix */

  int lx,         /* rows of sample matrix */
      dx,         /* columns, dimension of sample points */
      nCov,
	  info;

  krig_type krigType;

  krig_model *km;
  krig_result krigr[2]; 		  /* working storage kriging: 0: default, 1: temporary */

  glkrig_models_s(SEXP _R_covList, SEXP _R_Xmat, SEXP _R_data, SEXP R_krigType, int idx = 0) :
	  Xmat(0), lx(0), dx(0), nCov(0), info(0)
  {
	  int *dim = GET_DIMS(_R_Xmat);
	  lx = dim[0]; dx = dim[1];

	  if(isNull(_R_covList))
		ERR("Covariance model is NULL.");
	  nCov = LENGTH(_R_covList);

	  if(!isString(R_krigType))
	   Rf_error(_("Kriging type estimation should be a string: `R_krigType`."));
	  const char *krig_type = translateChar(asChar(R_krigType));
	  if( !strcmp("dual",krig_type) )
		  krigType = DUAL;
	  else if( !strcmp("var",krig_type) )
		  krigType = VARIANCE;
	  else ERR("unknown kriging type: choose `var` or `dual`.");

	  /* kriging results storage */
	  for(int i=0; i<2; i++) {
		if( (krigr[i] = new krig_result_s(nCov,lx,krigType)) == NULL)
		  MEM_ERR(1,krig_result_s)
	  }
	  CALLOCX(Xmat,lx*dx,double);
	  MEMCPY(Xmat,REAL(_R_Xmat),lx*dx);
	  if( (km =  new(std::nothrow) krig_model[nCov]) == NULL)
		  MEM_ERR(1,krig_model)
	  for(int i=0; i<nCov; i++) {
		 if( ( km[i] = new(std::nothrow) krig_model_s(VECTOR_ELT(_R_covList,i),
				 	 	 	 	  Xmat, REAL(AS_NUMERIC(VECTOR_ELT(AS_LIST(_R_data),i+idx))),
								  lx,dx,krigType) ) == NULL )
			 MEM_ERR(1,krig_model_s)
	  }

  }

  glkrig_models_s(SEXP _R_covList, double *_Xmat, double **_data, int _lx, int _dx,
		  krig_type _krigType, int idx = 0) :
		Xmat(0), lx(_lx), dx(_dx), nCov(LENGTH(_R_covList)),
		info(0), krigType(_krigType)
  {
  	  /* kriging results storage */
  	  for(int i=0; i<2; i++) {
  		if( (krigr[i] = new krig_result_s(nCov,lx,krigType)) == NULL)
  	  	  MEM_ERR(1,krig_result_s)
  	  }
	  CALLOCX(Xmat,lx*dx,double);
	  MEMCPY(Xmat,_Xmat,lx*dx);
  	  if( (km =  new(std::nothrow) krig_model[nCov]) == NULL)
  		  MEM_ERR(1,krig_model)

  	  for(int i=0; i<nCov; i++)
  	    if( (km[i] = new(std::nothrow) krig_model_s(VECTOR_ELT(_R_covList,i),Xmat,
  	    		 	 	 _data[i+idx],lx,dx,krigType)) == NULL)
  	    	MEM_ERR(1,krig_model_s)
  }

  ~glkrig_models_s() {
	 int k=0;
	  if(km != NULL) {
	  for(k=0; k<nCov;k++)
	    DELETE(km[k])
	  DELETE_ARR(km)
	 }
	 for(k=0; k<2; k++)
	   DELETE(krigr[k])
	 FREE(Xmat)
   }

  void kriging(double *x, double *m, double *s, double *w, int &err);

  inline int intern_kriging(double *x) {
	  kriging(x, krigr[0]->mean, krigr[0]->sig2, krigr[0]->w, info);
	  return info;
   }

  void jacobian(double *x, double *mean, double *jac, double *fdwork, int &err);

  inline int intern_jacobian(double *x, double *jac, double *fdwork) {
	  jacobian(x, krigr[0]->mean, jac, fdwork, info);
	  return info;
  }

} glkrig_models, *GLKM;

typedef struct cv_model_s {
	int nc,			/* number of CV models, each made up of nCov statistics */
		np,dx,	    /* overall number and dimension of sample point(s) 	*/
		nCov;
	double fnc, tmp;

	GLKM *cm;
	krig_type krigType;
	double *Xmat, *s2, *ybar, *ytil, *y0;

	cv_model_s(SEXP R_Xmat, SEXP R_data, SEXP R_cm) :
	 nc(0), np(0), dx(0), nCov(0), fnc(0), tmp(0),
	 cm(0), krigType(VARIANCE),
	 Xmat(0), s2(0), ybar(0), ytil(0), y0(0)
	{
		set(R_Xmat, R_data, R_cm);
	}

	~cv_model_s() {
	  if(cm != NULL) {
		 for(int i=0; i<nc; i++)
		   DELETE(cm[i])
		 DELETE_ARR(cm)
	  }
	  FREE(y0)
	  FREE(s2)
	  FREE(ybar)
	  FREE(ytil)
	}

	void set(SEXP R_Xmat, SEXP R_data, SEXP R_cm);

	void cvError(double *x, krig_model *km, double *cv, int &info);

} cv_model_t, *CVKM;

typedef struct {
	double *qimat,
		   *vmat,
		   *varS,
		   *score;

} ql_storage_t;


typedef struct ql_data_s {

  // always to be allocated
  double *obs,          /* Statistics of reference parameter */
         *qtheta,       /* working array of length nCov, (T(x)-E[T(X)]) */
		 *vmat,         /* variance matrix of statistics */
		 *vmat_work,    /* variance matrix of statistics */
		 *work,	        /* copy of diagonal terms of variance matrix of statistics */
		 *workx,		/* for kriging the variance matrix */
		 *jactmp,*fdjac,/* jac working array */
		 *fdwork,		/* working array for FD approxmiation ETX/dtheta (length: number of statistics) */
		 *fdscore,		/* working array for FD approxmiation   Q/dtheta (length: number of parameters) */
		 *Atmp, 	    /* temporary working matrix */
  	  	 *tmp;			/* working array with length number of statistics */

  ql_options_t qlopts;

  int dx,lx,nCov; /* dimension parameter, number of samples, number of statistics */

  // alloc and init
  ql_data_s(SEXP R_obs, SEXP R_Vmat, SEXP R_qlopts, int _nCov, int *_dims) :
		  obs(0), qtheta(0), vmat(0), vmat_work(0), work(0),
		  workx(0), jactmp(0), fdjac(0), fdwork(0), fdscore(0), Atmp(0), tmp(0), qlopts(R_qlopts),
		  dx(_dims[1]), lx(_dims[0]), nCov(_nCov)
  {
	 int nCov2 = SQR(nCov);

	 CALLOCX(qtheta,nCov,double);
	 CALLOCX(fdwork,nCov,double);
	 CALLOCX(fdscore,dx,double);
	 CALLOCX(Atmp,dx*nCov,double);
	 CALLOCX(obs,nCov,double);
	 CALLOCX(tmp,nCov,double);
	 CALLOCX(jactmp,dx*nCov,double);
	 CALLOCX(fdjac,dx*nCov,double);
	 CALLOCX(vmat,nCov2,double);
	 CALLOCX(vmat_work,nCov2,double);

	 for(int i=0;i<dx*nCov;i++)
		 fdjac[i]=jactmp[i]=Atmp[i]=0;

	 double *_obs = REAL(AS_NUMERIC(R_obs));
	 for(int i=0;i<nCov;i++) {
		 tmp[i]=obs[i]=_obs[i];
	 }

	 if(qlopts.varType < KRIG) {
		if(isNull(R_Vmat) || !isMatrix(R_Vmat))
		  ERR("Variance `R_Vmat` is either NULL or not a matrix.");
		double *_vmat = REAL(R_Vmat);
		MEMCPY(vmat,_vmat,nCov2);
		MEMCPY(vmat_work,_vmat,nCov2);
		if(qlopts.varType == MEAN) {
		  CALLOCX(work,nCov,double);
		  splitDiag(vmat_work,nCov,work);
		}
	 } else { /* for kriging variance matrix */
		 CALLOCX(workx,nCov2,double);
	 }
  };

  ~ql_data_s() {
	  FREE(jactmp)
	  FREE(fdjac)
	  FREE(fdwork)
	  FREE(fdscore)
	  FREE(tmp)
	  FREE(obs)
	  FREE(Atmp)
	  FREE(qtheta)
	  FREE(vmat)
	  FREE(vmat_work)
	  FREE(work)
	  FREE(workx)
  };

} ql_data_t, *ql_data;

typedef struct ql_model_s {
	//ql_data_t qld;

	GLKM glkm;			  			/* Kriging models for statistics */
	GLKM varkm;  					/* Kriging models for variance matrices */
	CVKM cvmod;						/* CV models (Cross-Validation kriging models) */
	ql_data qld;
	ql_model_s *cvnext,	*cvlast;	/* CV models of QL models */
    ql_storage_t qlsolve;			/* general storage for solving kriging equations */

	double *qimat,  /* expected Quasi-Information */
		   *varS,	/* modified Quasi-Information */
		   *score,  /* quasi score vector */
		   *jac,    /* gradient statistics */
		   *qiobs;  /* observed Quasi-Information = dS/dtheta by FD */

	double *lower, *upper;   /* bounds on the parameter space */

	int dx,dxdx,lx,nCov,nCov2,nCV,info;

	// R_qsd, R_qlopts, R_Vmat, R_cm, (Rboolean) type
	ql_model_s(SEXP R_qlm, SEXP R_qlopts, SEXP R_cvm) :
			glkm(NULL), varkm(NULL), cvmod(NULL), qld(NULL),cvnext(0), cvlast(0),
			qimat(0), score(0), jac(0), qiobs(0), lower(0), upper(0), nCV(0), info(0)
	{
		cvlast=this;
		SEXP R_qsd  = getListElement(R_qlm, "qlm");
		SEXP R_Xmat = getListElement(R_qlm, "X");
		SEXP R_Vmat = getListElement(R_qlm, "VTX");
		/* bounds */
		lower = REAL(getListElement( R_qsd, "lower" ));
		upper = REAL(getListElement( R_qsd, "upper" ));
		/* list of covariance models*/
		SEXP R_covT = getListElement( R_qsd, "covT" );
		/* list of covariance models for variances*/
		SEXP R_covL = getListElement( R_qsd, "covL" );
		/* reference statistic values*/
		SEXP R_obs = getListElement( R_qsd, "obs" );
		/* ql data list */
		SEXP R_qldata = getListElement( R_qsd, "qldata" );
		/* always set to "var" at R level */
		SEXP R_krigType = getListElement( R_qsd, "krig.type" );

		nCov = LENGTH(R_covT);
		if( (qld = new(std::nothrow) ql_data_s(R_obs,R_Vmat,R_qlopts,nCov,GET_DIMS(R_Xmat))) == NULL)
		    MEM_ERR(1,ql_data_s)

		lx = qld->lx;
		dx = qld->dx;
		dxdx = SQR(dx);
		nCov2 = SQR(nCov);

		 if(qld->qlopts.varType > MEAN){
			 if(isNull(R_covL))
			   ERR("Covariance model for kriging variance matrix was not set.");
			 int idx = dx+2*nCov;
			 // always use VARIANCE as krigType for kriging variance matrix of statistics
			 if( (varkm = new(std::nothrow) glkrig_models_s(R_covL,R_Xmat,R_qldata,R_krigType,idx)) == NULL)
			   MEM_ERR(1,glkrig_models_s)
		 }
		 if(!isNull(R_cvm)) {								/* Cross-validation models initialization */
			 if(qld->qlopts.useCV) {			 			/* CV for statistics */
			  if( (cvmod = new(std::nothrow) cv_model_s(R_Xmat,R_qldata,R_cvm)) == NULL)
			    { MEM_ERR(1,cv_model_s) }
			 } else { 										/* CV for quasi-score error */
		  	  nCV = LENGTH(R_cvm);
			  ql_model_s *p = this;
			  for(int i=0; i<nCV; i++){
			  	  if((p->cvnext = new(std::nothrow) ql_model_s(VECTOR_ELT(R_cvm,i),R_qlopts,R_NilValue)) == NULL)
			       { MEM_ERR(1,ql_model_s); }
			   	  p = p->cvnext;
			  }
			  cvlast = p;
			 }
		 }

		 if( (glkm = new(std::nothrow) glkrig_models_s(R_covT,R_Xmat,R_qldata,R_krigType,dx)) == NULL)
			MEM_ERR(1,glkrig_models_s)

		 CALLOCX(score,dx,double);
		 CALLOCX(qiobs,dxdx,double);
		 CALLOCX(qimat,dxdx,double);
		 CALLOCX(varS,dxdx,double);
		 CALLOCX(jac,dx*nCov,double);

		 /* ql storage for solving */
		 CALLOCX(qlsolve.score,dx,double);
		 CALLOCX(qlsolve.qimat,dxdx,double);
		 CALLOCX(qlsolve.varS,dxdx,double);
		 CALLOCX(qlsolve.vmat,nCov2,double);
	}

	~ql_model_s()
	{
	  DELETE(glkm)
	  DELETE(varkm)
	  DELETE(cvmod)
	  DELETE(cvnext);  				/* delete linked list */
	  DELETE(qld)

	  FREE(qiobs)
	  FREE(qimat)
	  FREE(varS)
	  FREE(jac)
	  FREE(score)

	  /* storage */
	  FREE(qlsolve.score)
	  FREE(qlsolve.qimat)
	  FREE(qlsolve.varS)
	  FREE(qlsolve.vmat)
	}

	/* variance of quasi-score within the approximating kriging models */
	int intern_varScore(double *x, double *jac, double *qimat, double *vars);
	/* observed quasi-information matrix (FD approximation of quasiscore vector) */
	int intern_quasiObs(double *x, double *score, double *qiobs);

	inline int intern_quasiScore(double *jac, double *score) {
		 quasiScore(glkm->krigr[0]->mean,jac,qld->vmat,score,info);
		 return info;
	}

	double qfValue(double *score, double *varS);
    int intern_qfScore(double *x){ return qfScore(x,jac,score,qimat,varS); }

    void qDPressVec(SEXP R_points, double *mpress);

    double intern_wlogdet(double *x, double w);
    inline double intern_qfValue() { return qfValue(score,varS); }
    inline double intern_qfScoreStat(double *x) { return qfScoreStat(x,jac,score,qimat,varS); }

	void varMatrix(double *x, double *vmat, int &err);

	int varScore(double *x, double *varS);
	double trVarScore(double *x);

	int qfScore(double *x, double *jac, double *score, double *qimat, double *varS);
	double qfScoreStat(double *x, double *jac, double *score, double *qimat, double *varS);

	int intern_quasiInfo(double *jac, double *qimat);
	void quasiScore(double *mean, double *jac, double *vmat, double *score, int &err);

	double mahalDist(double *x,double *jac,double *score);
	double intern_mahalanobis(double *x) { return mahalDist(x,jac,score); }


} ql_model_t, *ql_model;

int addMatrixDiag(double *vmat, int nx, double *w);
int addVar(double *sig2, int nc, double *vmat, double *work);

//////////////////////////// wrapper ///////////////////////////////////////////////////////

void wrap_intern_kriging(double *x, void *data,  double *mean, int &err);
void wrap_intern_quasiScore(double *x, void *data, double *score, int &err);
double wrap_intern_quasiDeviance(double *x, void *data, int &err );

//////////////////////////// attributes setting ////////////////////////////////////////////

void setVmatAttrib(ql_model qlm, SEXP R_VmatNames, SEXP R_ans);

///////////////////// KRIGING //////////////////////////////////////////////////////////////

SEXP kriging(SEXP R_Xmat, SEXP R_data, SEXP R_points, SEXP R_CovT, SEXP R_krigType);
/* either kriging all statistics in dual form or with calculating variances
 * and gradient estimation type is the same for all covariance structures  */
SEXP estimateJacobian(SEXP R_Xmat, SEXP R_data, SEXP R_points, SEXP R_covList, SEXP R_krigtype);
SEXP getDualKrigingWeights(SEXP R_Cmat, SEXP R_Fmat, SEXP R_data);

///////////////////// QUASI_LIKELIHOOD  ////////////////////////////////////////////////////

SEXP finalizeQL();
SEXP initQL(SEXP R_qsd, SEXP R_qlopts, SEXP R_cvT);

SEXP mahalanobis(SEXP R_points, SEXP R_qsd, SEXP R_qlopts, SEXP R_cm, SEXP R_qdValue, SEXP R_w);
SEXP quasiDeviance(SEXP R_points, SEXP R_qsd, SEXP R_qlopts, SEXP R_cm, SEXP R_qdValue, SEXP R_w);

/* return only quasi-deviance and Mahalanobis distance */
SEXP qDValue(SEXP R_point);
SEXP mahalValue(SEXP R_point);
/* return list of quasi-score, quasi-information and others */
SEXP internQD(SEXP R_point);
SEXP internMD(SEXP R_point);


/* commented out - testing only */
//SEXP cvError(SEXP R_point, SEXP R_Xmat, SEXP R_data, SEXP R_covT, SEXP R_krigType, SEXP R_cm);

#endif /* KRIGING_H_ */

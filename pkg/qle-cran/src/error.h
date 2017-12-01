/**
 * @file        error.h
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       Debugging definitions and error functions
 *
 *
 * Explanation: Trigger errors and warnings either directly by calling
 *              ERR, WRR, XERR, XWRR makros or store them in some
 *              objects to read from R and trigger them as you want
 */

#ifndef DEBUG_H_
#define DEBUG_H_

#include "basics.h"

#include <R.h>
#include <Rdefines.h>

#define NO_ERROR 0
#define NO_WARNING 0

#define PRINT_MSG(s) Rprintf("[%s: %u] %s \n", __FILE__, (unsigned)__LINE__, s)

extern int  PL;
extern char C_MSG_BUFFER[100];

#define LOG_ERROR(X, MSG) { \
  errorMSG(X, C_MSG_BUFFER); \
  Rprintf("%s (code=%d)\n %s", C_MSG_BUFFER, X, MSG); \
}

#define LOG_WARNING(X, MSG) { \
  warningMSG(X, C_MSG_BUFFER); \
  Rprintf("%s (code=%d)\n %s", C_MSG_BUFFER, X, MSG); \
}


/* trigger error and warnings */
#define ERR(MSG) { \
  Rprintf("%s (line=%u)\n", __FILE__, (unsigned)__LINE__); \
  Rf_error(_(MSG)); \
}

#define XERR(X,MSG) { \
     errorMSG(X, C_MSG_BUFFER); \
     Rprintf("%s in %s  (line=%u)\n", C_MSG_BUFFER, __FILE__, (unsigned)__LINE__); \
     Rf_error(_(MSG)); \
  }

// X is message string
#define WRR(MSG) { \
  Rprintf("%s (line=%u)\n", __FILE__, (unsigned)__LINE__); \
  Rf_warning(_(MSG)); \
}

// X is integer error code
#define XWRR(X,MSG) { \
    errorMSG(X, C_MSG_BUFFER); \
    Rprintf("%s \n %s (line=%u)\n", C_MSG_BUFFER, __FILE__, (unsigned)__LINE__); \
    Rf_warning(_(MSG)); \
}

// memory allocation errors
#define MEM_ERR(n, t) { \
  Rprintf("%s (line=%u)\n", __FILE__, (unsigned)__LINE__); \
  Rprintf("(%.0f of %u bytes) \n", (double) (n), (unsigned)sizeof(t)); \
  Rf_error("Could not allocate memory."); \
}

typedef enum  {
        SYSTEM_ERROR = 1,
        ALLOC_ERROR = 2,
        MEMORY_ERROR = 3,
        MATH_ERROR = 4,
        NaN_ERROR = 5,
        LAPACK_ERROR = 100,
        LAPACK_QR_ERROR = 101,
        LAPACK_PMAT_ERROR = 102,
        LAPACK_SOLVE_ERROR = 103,
        LAPACK_INVERSION_ERROR = 104,
        LAPACK_FACTORIZE_ERROR = 105,
        PARAM_ERROR = 500,
        FINAL_ERRROR = 1000      /*not to be changed */

} error_type;

template<class Type>
void printArray(const char fmt[], const Type *v, int *lx) {
  Rprintf("[");
  for(int i=0; i < *lx; i++)
     Rprintf( fmt , v[i]);
  Rprintf("]");
}


/////////////////////// Error /////////////////////////////////////////////////////////////////////////////////////
extern void errorMSG(int, char*);
extern void warningMSG(int , char *, char* );

void printMatrix(const char ch[], const double *mat, int *irows, int *icols);
void printVector(const char ch[], const double *vec, int *lx);

void print_R_matrix( SEXP A, const std::string& matrix_name);
void print_R_vector( SEXP v, const std::string& vector_name);

#define DEBUG_INFO(_x) \
	do { \
		fprintf(stderr, "===== BEGIN: DEBUG block =====\n" \
				"file: %s, line: %u\n", \
				__FILE__, (unsigned)__LINE__); \
				_x; \
				fprintf(stderr, "===== E N D: DEBUG block =====\n"); \
				fflush(stderr); \
	} while (0)

#define DEBUG_PRINT_VECTOR_R( c , cstr ) { print_R_vector( c , cstr ); }
#define DEBUG_DUMP_VAR(x,fmt) { Rprintf("%s:%u: %s=\n" fmt, __FILE__, (unsigned)__LINE__, #x, x); }
#define DEBUG_PRINT_MATRIX_R( C , cstr) { print_R_matrix( C , cstr); }

/////////////////// Some R convenience macros ////////////////////////////////////

#define RMATRIX(m,i,j) (REAL(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
#define GET_DIMS(A) INTEGER(coerceVector (getAttrib( (A), (R_DimSymbol) ) , INTSXP) )

#define SET_CLASS_NAME(RObject,ClassName) {            \
  SEXP RClass;                                         \
  PROTECT(RClass=allocVector(STRSXP,1));               \
  SET_STRING_ELT( (RClass) ,0, mkChar( (ClassName) ));  \
  classgets( (RObject), (RClass) );                    \
  UNPROTECT(1);                                        \
}

// Dimension names of a matrix, UNPROTECT later
#define SET_DIMNAMES_MATRIX(RObject,RNamedObject){ 		\
  SEXP R_names = getAttrib(RNamedObject,R_NamesSymbol); \
  SET_VECTOR_ELT(RObject, 0, R_names);					\
  SET_VECTOR_ELT(RObject, 1, R_names);					\
}

// Dimension names of a matrix, UNPROTECT later
#define SET_DIMNAMES_MATRIX2(RObject,RNamedObject0,RNamedObject1){ 		\
  SEXP R_names0 = getAttrib(RNamedObject0,R_NamesSymbol); 				\
  SEXP R_names1 = getAttrib(RNamedObject1,R_NamesSymbol); 				\
  SET_VECTOR_ELT(RObject, 0, R_names0);									\
  SET_VECTOR_ELT(RObject, 1, R_names1);									\
}


#endif /* DEBUG_H */

/*!
 \file  proto.h
 \brief This file contains function prototypes

 \author David C. Anastasiu
 */
#ifndef _L2AP_PROTO_H_
#define _L2AP_PROTO_H_

#include "includes.h"

#ifdef __cplusplus
extern "C"
{
#endif


/** MACROS  **/

#define F_NEQUALS_P(a,b,p) \
	do { if(fabs((a)-(b) ) > p) return 0; } while(0)


#define FL2INT(a,ll,f) \
	( (uint16_t) round((a - ll) * f) )


/**
 * Debug macros
 */
#define __BACKTRACE() \
  do { \
    void * _buffer[255]; \
    int _size = backtrace(_buffer,255); \
    char ** _strings = backtrace_symbols(_buffer,_size); \
    int _i; \
    fprintf(stderr,"\n"); \
    for (_i=0;_i<_size;++_i) { \
      fprintf(stderr,"%d:[%p] %s\n",_i,_buffer[_i],_strings[_i]); \
    } \
    free(_strings); \
  } while(0)

#define eprintf(fmt, ...) \
  do { \
    fprintf(stderr , "ERROR: "fmt , ##__VA_ARGS__); \
    fflush(stderr); \
  } while (0)

#ifndef NDEBUG
  #define wprintf(fmt, ...) \
    do { \
      fprintf(stderr , "WARN: "fmt , ##__VA_ARGS__); \
      fflush(stderr); \
    } while (0)

  #define dprintf(fmt, ...) \
    do { \
      fprintf(stdout , "DEBUG: "fmt , ##__VA_ARGS__); \
      fflush(stdout); \
    } while (0)

  #define DA_ASSERT(cond,fmt, ...) \
    do { \
      if (!(cond)) { \
        eprintf(fmt, ##__VA_ARGS__); \
        /* print the stack trace */ \
        __BACKTRACE(); \
        assert(0); /* will always fail */ \
      } \
    } while (0)

  #define DA_ASSERT_EQUALS(a,b,fmt) \
    DA_ASSERT(a == b,"("#a" = "fmt") != ("#b" = "fmt")",a,b)
#else
	#define wprintf(fmt, ...)
	#define dprintf(fmt, ...)
	#define DL_ASSERT(cnd,fmt, ...)
	#define DL_ASSERT_EQUALS(a,b,fmt)
#endif


/**
 * Create memory allocation function prototypes for generated types
 */
// base types: ptr, idx, val, etc
GK_MKALLOC_PROTO(da_p,      ptr_t)
GK_MKALLOC_PROTO(da_i,      idx_t)
GK_MKALLOC_PROTO(da_v,      val_t)
GK_MKALLOC_PROTO(da_c,      char)
GK_MKALLOC_PROTO(da_uc,     unsigned char)
GK_MKALLOC_PROTO(da_i8,     int8_t)
GK_MKALLOC_PROTO(da_ui8,    uint8_t)
GK_MKALLOC_PROTO(da_i16,    int16_t)
GK_MKALLOC_PROTO(da_ui16,   uint16_t)
GK_MKALLOC_PROTO(da_i32,    int32_t)
GK_MKALLOC_PROTO(da_ui32,   uint32_t)
GK_MKALLOC_PROTO(da_i64,    int64_t)
GK_MKALLOC_PROTO(da_ui64,   uint64_t)
GK_MKALLOC_PROTO(da_z,      ssize_t)
GK_MKALLOC_PROTO(da_uz,     size_t)
GK_MKALLOC_PROTO(da_f,      float)
GK_MKALLOC_PROTO(da_d,      double)
GK_MKALLOC_PROTO(da_l,      long)
GK_MKALLOC_PROTO(da_ul,     unsigned long)
GK_MKALLOC_PROTO(da_s,      da_sim_t)
// ptr_t key, various vals
GK_MKALLOC_PROTO(da_ppkv,   da_ppkv_t)
GK_MKALLOC_PROTO(da_pikv,   da_pikv_t)
GK_MKALLOC_PROTO(da_pvkv,   da_pvkv_t)
GK_MKALLOC_PROTO(da_pckv,   da_pckv_t)
GK_MKALLOC_PROTO(da_pi32kv, da_pi32kv_t)
GK_MKALLOC_PROTO(da_pi64kv, da_pi64kv_t)
GK_MKALLOC_PROTO(da_pzkv,   da_pzkv_t)
GK_MKALLOC_PROTO(da_pfkv,   da_pfkv_t)
GK_MKALLOC_PROTO(da_pdkv,   da_pdkv_t)
GK_MKALLOC_PROTO(da_pskv,   da_pskv_t)
/* idx_t key, various vals */
GK_MKALLOC_PROTO(da_ipkv,   da_ipkv_t)
GK_MKALLOC_PROTO(da_iikv,   da_iikv_t)
GK_MKALLOC_PROTO(da_ivkv,   da_ivkv_t)
GK_MKALLOC_PROTO(da_ickv,   da_ickv_t)
GK_MKALLOC_PROTO(da_ii32kv, da_ii32kv_t)
GK_MKALLOC_PROTO(da_ii64kv, da_ii64kv_t)
GK_MKALLOC_PROTO(da_izkv,   da_izkv_t)
GK_MKALLOC_PROTO(da_ifkv,   da_ifkv_t)
GK_MKALLOC_PROTO(da_idkv,   da_idkv_t)
GK_MKALLOC_PROTO(da_iskv,   da_iskv_t)
/* size_t key, various vals */
GK_MKALLOC_PROTO(da_upkv,   da_upkv_t)
GK_MKALLOC_PROTO(da_uikv,   da_uikv_t)
GK_MKALLOC_PROTO(da_uvkv,   da_uvkv_t)
GK_MKALLOC_PROTO(da_uckv,   da_uckv_t)
GK_MKALLOC_PROTO(da_ui32kv, da_ui32kv_t)
GK_MKALLOC_PROTO(da_ui64kv, da_ui64kv_t)
GK_MKALLOC_PROTO(da_uzkv,   da_uzkv_t)
GK_MKALLOC_PROTO(da_ufkv,   da_ufkv_t)
GK_MKALLOC_PROTO(da_udkv,   da_udkv_t)
GK_MKALLOC_PROTO(da_uskv,   da_uskv_t)

/**
 * BLAS prototypes
 */

DA_MKBLAS_PROTO(da_p,   ptr_t,    ptr_t)
DA_MKBLAS_PROTO(da_i,   idx_t,    idx_t)
DA_MKBLAS_PROTO(da_v,   val_t,    val_t)
DA_MKBLAS_PROTO(da_c,   char,     int)
DA_MKBLAS_PROTO(da_i32, int32_t,  int32_t)
DA_MKBLAS_PROTO(da_i64, int64_t,  int64_t)
DA_MKBLAS_PROTO(da_z,   ssize_t,  ssize_t)
DA_MKBLAS_PROTO(da_f,   float,    float)
DA_MKBLAS_PROTO(da_d,   double,   double)

/**
 * RAND prototypes
 */
GK_MKRANDOM_PROTO(da_p,   size_t, ptr_t)
GK_MKRANDOM_PROTO(da_i,   size_t, idx_t)
GK_MKRANDOM_PROTO(da_v,   size_t, val_t)
GK_MKRANDOM_PROTO(da_c,   size_t, char)
GK_MKRANDOM_PROTO(da_f,   size_t, float)
GK_MKRANDOM_PROTO(da_d,   size_t, double)
GK_MKRANDOM_PROTO(da_z,   size_t, ssize_t)


/* sort.c */

#define DA_MKSORT_PROTO(PRFX, TYPE) \
void     PRFX ## sorti(size_t n, TYPE *base);\
void     PRFX ## sortd(size_t n, TYPE *base);\

DA_MKSORT_PROTO(da_p, ptr_t)
DA_MKSORT_PROTO(da_i, idx_t)
DA_MKSORT_PROTO(da_v, val_t)
DA_MKSORT_PROTO(da_ppkv, da_ppkv_t)
DA_MKSORT_PROTO(da_pikv, da_pikv_t)
DA_MKSORT_PROTO(da_pvkv, da_pvkv_t)
DA_MKSORT_PROTO(da_pckv, da_pckv_t)
DA_MKSORT_PROTO(da_pi32kv, da_pi32kv_t)
DA_MKSORT_PROTO(da_pi64kv, da_pi64kv_t)
DA_MKSORT_PROTO(da_pzkv, da_pzkv_t)
DA_MKSORT_PROTO(da_pfkv, da_pfkv_t)
DA_MKSORT_PROTO(da_pdkv, da_pdkv_t)
DA_MKSORT_PROTO(da_pskv, da_pskv_t)
DA_MKSORT_PROTO(da_ipkv, da_ipkv_t)
DA_MKSORT_PROTO(da_iikv, da_iikv_t)
DA_MKSORT_PROTO(da_ivkv, da_ivkv_t)
DA_MKSORT_PROTO(da_ickv, da_ickv_t)
DA_MKSORT_PROTO(da_ii32kv, da_ii32kv_t)
DA_MKSORT_PROTO(da_ii64kv, da_ii64kv_t)
DA_MKSORT_PROTO(da_izkv, da_izkv_t)
DA_MKSORT_PROTO(da_ifkv, da_ifkv_t)
DA_MKSORT_PROTO(da_idkv, da_idkv_t)
DA_MKSORT_PROTO(da_iskv, da_iskv_t)
DA_MKSORT_PROTO(da_upkv, da_upkv_t)
DA_MKSORT_PROTO(da_uikv, da_uikv_t)
DA_MKSORT_PROTO(da_uvkv, da_uvkv_t)
DA_MKSORT_PROTO(da_uckv, da_uckv_t)
DA_MKSORT_PROTO(da_ui32kv, da_ui32kv_t)
DA_MKSORT_PROTO(da_ui64kv, da_ui64kv_t)
DA_MKSORT_PROTO(da_uzkv, da_uzkv_t)
DA_MKSORT_PROTO(da_ufkv, da_ufkv_t)
DA_MKSORT_PROTO(da_udkv, da_udkv_t)
DA_MKSORT_PROTO(da_uskv, da_uskv_t)

/* omp.c */
#if !defined(_OPENMP)
void omp_set_num_threads(int num_threads);
int omp_get_num_threads(void);
int omp_get_max_threads(void);
int omp_get_thread_num(void);
int omp_get_num_procs(void);
int omp_in_parallel(void);
void omp_set_dynamic(int num_threads);
int omp_get_dynamic(void);
void omp_set_nested(int nested);
int omp_get_nested(void);
#endif


/* main.c */
//void readInputData(params_t *params);
//void preProcessData(params_t * params);
//void da_testMatricesEqual(params_t *params);
//void da_matrixInfo(params_t *params);
//void da_matrixIo(params_t *params);
//void simSearchSetup(params_t *params);
//void simSearchFinalize(params_t *params);
//void freeParams(params_t** params);
//
///* idxjoin.c */
//void ijFindNeighbors(params_t *params);
//
///* ap.c */
//void apFindNeighbors(params_t *params);
//
///* mmj.c */
//void mmjFindNeighbors(params_t *params);
//
///* l2ap.c */
//void l2apFindNeighbors(params_t *params);
//void l2apReorderDocs(params_t *params, da_csr_t **docs, idx_t *rsizes, idx_t *csizes, val_t *rwgts,
//		val_t *cwgts, idx_t *rperm, idx_t *cperm);
//
///* l2apblsh.c */
//void l2apblshFindNeighbors(params_t *params);


/* csr.c  */
da_csr_t *da_csr_Create();
void da_csr_Init(da_csr_t *mat);
void da_csr_Free(da_csr_t **mat);
void da_csr_FreeAll(da_csr_t **ptr1,...);
void da_csr_FreeBase(da_csr_t *mat, char type);
void da_csr_LoadBases(da_csr_t *csr);
void da_csr_Transfer(da_csr_t *from, da_csr_t *to);
void da_csr_FreeContents(da_csr_t *mat);
da_csr_t *da_csr_Copy(da_csr_t *mat);
da_csr_t *da_csr_ExtractSubmatrix(da_csr_t *mat, idx_t rstart, idx_t nrows);
da_csr_t *da_csr_ExtractRows(da_csr_t *mat, idx_t nrows, idx_t *rind);
void da_csr_ExtractRowsInto(da_csr_t *mat, da_csr_t *nmat, idx_t nrows, idx_t *rind);
da_csr_t *da_csr_ExtractPartition(da_csr_t *mat, idx_t *part, idx_t pid);
da_csr_t **da_csr_Split(da_csr_t *mat, idx_t *color);
da_csr_t *da_csr_Read(char *filename, char format, char readvals, char numbering);
void da_csr_Write(da_csr_t *mat, char *filename, char format, char writevals, char numbering);
void da_csr_Print(da_csr_t *mat);
char da_csr_isClutoOrCsr(char *file);
da_csr_t *da_csr_Prune(da_csr_t *mat, char what, idx_t minf, idx_t maxf);
da_csr_t *da_csr_LowFilter(da_csr_t *mat, char what, char norm, float fraction);
da_csr_t *da_csr_topKPlusFilter(da_csr_t *mat, char what, idx_t topk, val_t keepval);
da_csr_t *da_csr_ZScoreFilter(da_csr_t *mat, char what, float zscore);
void da_csr_CompactColumns(da_csr_t *mat);
void da_csr_CompactRows(da_csr_t *mat);
void da_csr_SortIndices(da_csr_t *mat, char what);
char da_csr_CheckSortedIndex(da_csr_t *mat, char what);
void da_csr_CreateIndex(da_csr_t *mat, char what);
void da_csr_Normalize(da_csr_t *mat, char what, char norm);
void da_csr_Scale(da_csr_t *mat, char type);
void da_csr_ComputeSums(da_csr_t *mat, char what);
void da_csr_ComputeSquaredNorms(da_csr_t *mat, char what);
val_t da_csr_ComputeSimilarity(da_csr_t *mat, idx_t i1, idx_t i2, char what, char simtype);
char da_csr_Compare(da_csr_t *a, da_csr_t *b, double p);
void da_csr_Grow(da_csr_t *mat, ptr_t newNnz);
val_t da_csr_partialDotProduct(const ptr_t *rowptr, const ptr_t *endptr,
		const idx_t *rowind, const val_t *rowval, idx_t a, idx_t b);
val_t da_csr_dotProduct(const ptr_t *rowptr, const idx_t *rowind,
		const val_t *rowval, idx_t a, idx_t b);
void da_getfilestats(char *fname, size_t *r_nlines, size_t *r_ntokens,
		size_t *r_max_nlntokens, size_t *r_nbytes);
idx_t da_csr_GetSimilarSmallerRows(da_csr_t *mat, idx_t rid, char noSelfSim,
		idx_t nqterms, idx_t *qind, val_t *qval, char simtype, idx_t nsim,
        float minsim, da_ivkv_t *hits, idx_t *i_marker, da_ivkv_t *i_cand);


/* select.c */
idx_t da_uvkvkselectd(size_t n, idx_t topk, da_uvkv_t *cand);
idx_t da_uvkvkselecti(size_t n, idx_t topk, da_uvkv_t *cand);
idx_t da_uikvkselectd(size_t n, idx_t topk, da_uikv_t *cand);
idx_t da_uikvkselecti(size_t n, idx_t topk, da_uikv_t *cand);
idx_t da_ufkvkselectd(size_t n, idx_t topk, da_ufkv_t *cand);
idx_t da_ufkvkselecti(size_t n, idx_t topk, da_ufkv_t *cand);
idx_t da_udkvkselectd(size_t n, idx_t topk, da_udkv_t *cand);
idx_t da_udkvkselecti(size_t n, idx_t topk, da_udkv_t *cand);

idx_t da_pvkvkselectd(size_t n, idx_t topk, da_pvkv_t *cand);
idx_t da_pvkvkselecti(size_t n, idx_t topk, da_pvkv_t *cand);
idx_t da_pikvkselectd(size_t n, idx_t topk, da_pikv_t *cand);
idx_t da_pikvkselecti(size_t n, idx_t topk, da_pikv_t *cand);
idx_t da_pfkvkselectd(size_t n, idx_t topk, da_pfkv_t *cand);
idx_t da_pfkvkselecti(size_t n, idx_t topk, da_pfkv_t *cand);
idx_t da_pdkvkselectd(size_t n, idx_t topk, da_pdkv_t *cand);
idx_t da_pdkvkselecti(size_t n, idx_t topk, da_pdkv_t *cand);

idx_t da_ivkvkselectd(size_t n, idx_t topk, da_ivkv_t *cand);
idx_t da_ivkvkselecti(size_t n, idx_t topk, da_ivkv_t *cand);
idx_t da_iikvkselectd(size_t n, idx_t topk, da_iikv_t *cand);
idx_t da_iikvkselecti(size_t n, idx_t topk, da_iikv_t *cand);
idx_t da_ifkvkselectd(size_t n, idx_t topk, da_ifkv_t *cand);
idx_t da_ifkvkselecti(size_t n, idx_t topk, da_ifkv_t *cand);
idx_t da_idkvkselectd(size_t n, idx_t topk, da_idkv_t *cand);
idx_t da_idkvkselecti(size_t n, idx_t topk, da_idkv_t *cand);



/* util.c */
char* da_getStringKey(const gk_StringMap_t *strmap, char id);
int da_getStringID(const gk_StringMap_t *strmap, char *key);
char da_getFileFormat(char *file, const char format);
char* da_getDataset(params_t *params);
void da_vWriteVector(char* filename, val_t* vec, idx_t size, char *separator);
void da_vWriteMatrix(char* filename, val_t** mat, idx_t nrows, idx_t ncols);
void da_dWriteMatrix(char* filename, double** mat, idx_t nrows, idx_t ncols);
void da_iWriteVector(char* filename, idx_t* vec, idx_t size, char *separator);
void da_pWriteVector(char* filename, idx_t* vec, idx_t size, char *separator);
void da_storeSim(params_t *params, idx_t i, idx_t j, val_t v);
char da_isFmtBinary(char fmt);
void da_printTimer(char *name, double time);
void da_printTimerLong(char *name, double time);
void da_addHigherNeighbors(da_csr_t* mat);
int da_log2(idx_t a);
int da_ispow2(idx_t a);
float da_flog2(float a);
double da_dlog2(double a);
val_t da_vlog2(val_t a);
void da_csrCompare(da_csr_t *a, da_csr_t *b, float eps, char compVals);
float da_fsqrt(float number);
void da_inversePermuteMatrix(da_csr_t **matP, idx_t* rowPerm, idx_t* colPerm);
void printCompileChoices();

/* cmdline.c */
void cmdline_parse(params_t *ctrl, int argc, char *argv[]);


#ifdef __cplusplus
}
#endif

#endif 

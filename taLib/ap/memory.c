/*!
 \file  memory.c
 \brief This file contains functions for allocating and freeing memory

 \author David C. Anastasiu
 */

#include "includes.h"

/**
 * Memory allocation functions
 */
// base types: ptr, idx, and val
GK_MKALLOC(da_p,      ptr_t)
GK_MKALLOC(da_i,      idx_t)
GK_MKALLOC(da_v,      val_t)
GK_MKALLOC(da_c,      char)
GK_MKALLOC(da_uc,     unsigned char)
GK_MKALLOC(da_i8,     int8_t)
GK_MKALLOC(da_ui8,    uint8_t)
GK_MKALLOC(da_i16,    int16_t)
GK_MKALLOC(da_ui16,   uint16_t)
GK_MKALLOC(da_i32,    int32_t)
GK_MKALLOC(da_ui32,   uint32_t)
GK_MKALLOC(da_i64,    int64_t)
GK_MKALLOC(da_ui64,   uint64_t)
GK_MKALLOC(da_z,      ssize_t)
GK_MKALLOC(da_uz,     size_t)
GK_MKALLOC(da_f,      float)
GK_MKALLOC(da_d,      double)
GK_MKALLOC(da_l,      long)
GK_MKALLOC(da_ul,     unsigned long)
GK_MKALLOC(da_s,      da_sim_t)
// ptr_t key, various vals
GK_MKALLOC(da_ppkv,   da_ppkv_t)
GK_MKALLOC(da_pikv,   da_pikv_t)
GK_MKALLOC(da_pvkv,   da_pvkv_t)
GK_MKALLOC(da_pckv,   da_pckv_t)
GK_MKALLOC(da_pi32kv, da_pi32kv_t)
GK_MKALLOC(da_pi64kv, da_pi64kv_t)
GK_MKALLOC(da_pzkv,   da_pzkv_t)
GK_MKALLOC(da_pfkv,   da_pfkv_t)
GK_MKALLOC(da_pdkv,   da_pdkv_t)
GK_MKALLOC(da_pskv,   da_pskv_t)
// idx_t key, various vals
GK_MKALLOC(da_ipkv,   da_ipkv_t)
GK_MKALLOC(da_iikv,   da_iikv_t)
GK_MKALLOC(da_ivkv,   da_ivkv_t)
GK_MKALLOC(da_ickv,   da_ickv_t)
GK_MKALLOC(da_ii32kv, da_ii32kv_t)
GK_MKALLOC(da_ii64kv, da_ii64kv_t)
GK_MKALLOC(da_izkv,   da_izkv_t)
GK_MKALLOC(da_ifkv,   da_ifkv_t)
GK_MKALLOC(da_idkv,   da_idkv_t)
GK_MKALLOC(da_iskv,   da_iskv_t)
// size_t key, various vals
GK_MKALLOC(da_upkv,   da_upkv_t)
GK_MKALLOC(da_uikv,   da_uikv_t)
GK_MKALLOC(da_uvkv,   da_uvkv_t)
GK_MKALLOC(da_uckv,   da_uckv_t)
GK_MKALLOC(da_ui32kv, da_ui32kv_t)
GK_MKALLOC(da_ui64kv, da_ui64kv_t)
GK_MKALLOC(da_uzkv,   da_uzkv_t)
GK_MKALLOC(da_ufkv,   da_ufkv_t)
GK_MKALLOC(da_udkv,   da_udkv_t)
GK_MKALLOC(da_uskv,   da_uskv_t)



/**
 * BLAS functions
 */
DA_MKBLAS(da_p,   da_up,    ptr_t,    ptr_t)
DA_MKBLAS(da_i,   da_ui,    idx_t,    idx_t)
DA_MKBLAS(da_v,   da_uv,    val_t,    val_t)
DA_MKBLAS(da_c,   da_uc,    char,     int)
DA_MKBLAS(da_i32, da_ui32,  int32_t,  int32_t)
DA_MKBLAS(da_i64, da_ui64,  int64_t,  int64_t)
DA_MKBLAS(da_z,   da_uz,    ssize_t,  ssize_t)
DA_MKBLAS(da_f,   da_uf,    float,    float)
DA_MKBLAS(da_d,   da_ud,    double,   double)


/**
 * RAND functions
 */
GK_MKRANDOM(da_p,   size_t, ptr_t)
GK_MKRANDOM(da_i,   size_t, idx_t)
GK_MKRANDOM(da_v,   size_t, val_t)
GK_MKRANDOM(da_c,   size_t, char)
GK_MKRANDOM(da_f,   size_t, float)
GK_MKRANDOM(da_d,   size_t, double)
GK_MKRANDOM(da_z,   size_t, ssize_t)


/**
 * Deal with no OMP config
 */
#if !defined(_OPENMP)
void omp_set_num_threads(int num_threads) { return; }
int omp_get_num_threads(void) { return 1; }
int omp_get_max_threads(void) { return 1; }
int omp_get_thread_num(void) { return 0; }
int omp_get_num_procs(void) { return 1; }
int omp_in_parallel(void) { return 0; }
void omp_set_dynamic(int num_threads) { return; }
int omp_get_dynamic(void) { return 0; }
void omp_set_nested(int nested) { return; }
int omp_get_nested(void) { return 0; }
#endif

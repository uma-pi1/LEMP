/*!
\file  sort.c
\brief This file contains various sorting routines

These routines are implemented using the GKSORT macro that is defined
in gk_qsort.h and is based on GNU's GLIBC qsort() implementation.

Additional sorting routines can be created using the same way that
these routines where defined.

\date   Started 2/27/2013
\author David C. Anastasiu
*/

#include "includes.h"



/*************************************************************************/
/*! Sorts an array of ptr_t in increasing order */
/*************************************************************************/
void da_psorti(size_t n, ptr_t *base)
{
#define ptr_lt(a, b) ((*a) < (*b))
  GK_MKQSORT(ptr_t, base, n, ptr_lt);
#undef ptr_lt
}


/*************************************************************************/
/*! Sorts an array of ptr_t in decreasing order */
/*************************************************************************/
void da_psortd(size_t n, ptr_t *base)
{
#define ptr_gt(a, b) ((*a) > (*b))
  GK_MKQSORT(ptr_t, base, n, ptr_gt);
#undef ptr_gt
}


/*************************************************************************/
/*! Sorts an array of idx_t in increasing order */
/*************************************************************************/
void da_isorti(size_t n, idx_t *base)
{
#define idx_lt(a, b) ((*a) < (*b))
  GK_MKQSORT(idx_t, base, n, idx_lt);
#undef idx_lt
}


/*************************************************************************/
/*! Sorts an array of idx_t in decreasing order */
/*************************************************************************/
void da_isortd(size_t n, idx_t *base)
{
#define idx_gt(a, b) ((*a) > (*b))
  GK_MKQSORT(idx_t, base, n, idx_gt);
#undef idx_gt
}


/*************************************************************************/
/*! Sorts an array of val_t in increasing order */
/*************************************************************************/
void da_vsorti(size_t n, val_t *base)
{
#define val_lt(a, b) ((*a) < (*b))
  GK_MKQSORT(val_t, base, n, val_lt);
#undef val_lt
}


/*************************************************************************/
/*! Sorts an array of val_t in decreasing order */
/*************************************************************************/
void da_vsortd(size_t n, val_t *base)
{
#define float_gt(a, b) ((*a) > (*b))
  GK_MKQSORT(val_t, base, n, float_gt);
#undef float_gt
}



/*************************************************************************/
/*! Sorts an array of da_ppkv_t in increasing order */
/*************************************************************************/
void da_ppkvsorti(size_t n, da_ppkv_t *base)
{
#define da_ppkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_ppkv_t, base, n, da_ppkv_lt);
#undef da_ppkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ppkv_t in decreasing order */
/*************************************************************************/
void da_ppkvsortd(size_t n, da_ppkv_t *base)
{
#define da_ppkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_ppkv_t, base, n, da_ppkv_gt);
#undef da_ppkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pikv_t in increasing order */
/*************************************************************************/
void da_pikvsorti(size_t n, da_pikv_t *base)
{
#define da_pikv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_pikv_t, base, n, da_pikv_lt);
#undef da_pikv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pikv_t in decreasing order */
/*************************************************************************/
void da_pikvsortd(size_t n, da_pikv_t *base)
{
#define da_pikv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_pikv_t, base, n, da_pikv_gt);
#undef da_pikv_gt
}



/*************************************************************************/
/*! Sorts an array of da_pvkv_t in increasing order */
/*************************************************************************/
void da_pvkvsorti(size_t n, da_pvkv_t *base)
{
#define da_pvkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_pvkv_t, base, n, da_pvkv_lt);
#undef da_pvkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pvkv_t in decreasing order */
/*************************************************************************/
void da_pvkvsortd(size_t n, da_pvkv_t *base)
{
#define da_pvkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_pvkv_t, base, n, da_pvkv_gt);
#undef da_pvkv_gt
}



/*************************************************************************/
/*! Sorts an array of da_pckv_t in increasing order */
/*************************************************************************/
void da_pckvsorti(size_t n, da_pckv_t *base)
{
#define da_pckv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_pckv_t, base, n, da_pckv_lt);
#undef da_pckv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pckv_t in decreasing order */
/*************************************************************************/
void da_pckvsortd(size_t n, da_pckv_t *base)
{
#define da_pckv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_pckv_t, base, n, da_pckv_gt);
#undef da_pkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pi32kv_t in increasing order */
/*************************************************************************/
void da_pi32kvsorti(size_t n, da_pi32kv_t *base)
{
#define da_pi32kv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_pi32kv_t, base, n, da_pi32kv_lt);
#undef da_pi32kv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pi32kv_t in decreasing order */
/*************************************************************************/
void da_pi32kvsortd(size_t n, da_pi32kv_t *base)
{
#define da_pi32kv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_pi32kv_t, base, n, da_pi32kv_gt);
#undef da_pi32kv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pi64kv_t in increasing order */
/*************************************************************************/
void da_pi64kvsorti(size_t n, da_pi64kv_t *base)
{
#define da_pi64kv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_pi64kv_t, base, n, da_pi64kv_lt);
#undef da_pi64kv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pi64kv_t in decreasing order */
/*************************************************************************/
void da_pi64kvsortd(size_t n, da_pi64kv_t *base)
{
#define da_pi64kv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_pi64kv_t, base, n, da_pi64kv_gt);
#undef da_pi64kv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pzkv_t in increasing order */
/*************************************************************************/
void da_pzkvsorti(size_t n, da_pzkv_t *base)
{
#define da_pzkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_pzkv_t, base, n, da_pzkv_lt);
#undef da_pzkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pzkv_t in decreasing order */
/*************************************************************************/
void da_pzkvsortd(size_t n, da_pzkv_t *base)
{
#define da_pzkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_pzkv_t, base, n, da_pzkv_gt);
#undef da_pzkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pfkv_t in increasing order */
/*************************************************************************/
void da_pfkvsorti(size_t n, da_pfkv_t *base)
{
#define da_pfkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_pfkv_t, base, n, da_pfkv_lt);
#undef da_pfkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pfkv_t in decreasing order */
/*************************************************************************/
void da_pfkvsortd(size_t n, da_pfkv_t *base)
{
#define da_pfkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_pfkv_t, base, n, da_pfkv_gt);
#undef da_pfkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pdkv_t in increasing order */
/*************************************************************************/
void da_pdkvsorti(size_t n, da_pdkv_t *base)
{
#define da_pdkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_pdkv_t, base, n, da_pdkv_lt);
#undef da_pdkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pdkv_t in decreasing order */
/*************************************************************************/
void da_pdkvsortd(size_t n, da_pdkv_t *base)
{
#define da_pdkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_pdkv_t, base, n, da_pdkv_gt);
#undef da_pdkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pskv_t in increasing order */
/*************************************************************************/
void da_pskvsorti(size_t n, da_pskv_t *base)
{
#define da_pskv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_pskv_t, base, n, da_pskv_lt);
#undef da_pskv_lt
}


/*************************************************************************/
/*! Sorts an array of da_pskv_t in decreasing order */
/*************************************************************************/
void da_pskvsortd(size_t n, da_pskv_t *base)
{
#define da_pskv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_pskv_t, base, n, da_pskv_gt);
#undef da_pskv_gt
}


/* da_i ## kv pairs */


/*************************************************************************/
/*! Sorts an array of da_ipkv_t in increasing order */
/*************************************************************************/
void da_ipkvsorti(size_t n, da_ipkv_t *base)
{
#define da_ipkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_ipkv_t, base, n, da_ipkv_lt);
#undef da_ipkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ipkv_t in decreasing order */
/*************************************************************************/
void da_ipkvsortd(size_t n, da_ipkv_t *base)
{
#define da_ipkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_ipkv_t, base, n, da_ipkv_gt);
#undef da_ipkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pikv_t in increasing order */
/*************************************************************************/
void da_iikvsorti(size_t n, da_iikv_t *base)
{
#define da_iikv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_iikv_t, base, n, da_iikv_lt);
#undef da_iikv_lt
}


/*************************************************************************/
/*! Sorts an array of da_iikv_t in decreasing order */
/*************************************************************************/
void da_iikvsortd(size_t n, da_iikv_t *base)
{
#define da_iikv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_iikv_t, base, n, da_iikv_gt);
#undef da_iikv_gt
}



/*************************************************************************/
/*! Sorts an array of da_ivkv_t in increasing order */
/*************************************************************************/
void da_ivkvsorti(size_t n, da_ivkv_t *base)
{
#define da_ivkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_ivkv_t, base, n, da_ivkv_lt);
#undef da_ivkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ivkv_t in decreasing order */
/*************************************************************************/
void da_ivkvsortd(size_t n, da_ivkv_t *base)
{
#define da_ivkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_ivkv_t, base, n, da_ivkv_gt);
#undef da_ivkv_gt
}



/*************************************************************************/
/*! Sorts an array of da_ickv_t in increasing order */
/*************************************************************************/
void da_ickvsorti(size_t n, da_ickv_t *base)
{
#define da_ickv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_ickv_t, base, n, da_ickv_lt);
#undef da_ickv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ickv_t in decreasing order */
/*************************************************************************/
void da_ickvsortd(size_t n, da_ickv_t *base)
{
#define da_ickv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_ickv_t, base, n, da_ickv_gt);
#undef da_pkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_ii32kv_t in increasing order */
/*************************************************************************/
void da_ii32kvsorti(size_t n, da_ii32kv_t *base)
{
#define da_ii32kv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_ii32kv_t, base, n, da_ii32kv_lt);
#undef da_ii32kv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ii32kv_t in decreasing order */
/*************************************************************************/
void da_ii32kvsortd(size_t n, da_ii32kv_t *base)
{
#define da_ii32kv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_ii32kv_t, base, n, da_ii32kv_gt);
#undef da_ii32kv_gt
}




/*************************************************************************/
/*! Sorts an array of da_ii64kv_t in increasing order */
/*************************************************************************/
void da_ii64kvsorti(size_t n, da_ii64kv_t *base)
{
#define da_ii64kv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_ii64kv_t, base, n, da_ii64kv_lt);
#undef da_ii64kv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ii64kv_t in decreasing order */
/*************************************************************************/
void da_ii64kvsortd(size_t n, da_ii64kv_t *base)
{
#define da_ii64kv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_ii64kv_t, base, n, da_ii64kv_gt);
#undef da_ii64kv_gt
}




/*************************************************************************/
/*! Sorts an array of da_izkv_t in increasing order */
/*************************************************************************/
void da_izkvsorti(size_t n, da_izkv_t *base)
{
#define da_izkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_izkv_t, base, n, da_izkv_lt);
#undef da_izkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_izkv_t in decreasing order */
/*************************************************************************/
void da_izkvsortd(size_t n, da_izkv_t *base)
{
#define da_izkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_izkv_t, base, n, da_izkv_gt);
#undef da_izkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_ifkv_t in increasing order */
/*************************************************************************/
void da_ifkvsorti(size_t n, da_ifkv_t *base)
{
#define da_ifkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_ifkv_t, base, n, da_ifkv_lt);
#undef da_ifkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ifkv_t in decreasing order */
/*************************************************************************/
void da_ifkvsortd(size_t n, da_ifkv_t *base)
{
#define da_ifkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_ifkv_t, base, n, da_ifkv_gt);
#undef da_ifkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_idkv_t in increasing order */
/*************************************************************************/
void da_idkvsorti(size_t n, da_idkv_t *base)
{
#define da_idkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_idkv_t, base, n, da_idkv_lt);
#undef da_idkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_idkv_t in decreasing order */
/*************************************************************************/
void da_idkvsortd(size_t n, da_idkv_t *base)
{
#define da_idkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_idkv_t, base, n, da_idkv_gt);
#undef da_idkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_iskv_t in increasing order */
/*************************************************************************/
void da_iskvsorti(size_t n, da_iskv_t *base)
{
#define da_iskv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_iskv_t, base, n, da_iskv_lt);
#undef da_iskv_lt
}


/*************************************************************************/
/*! Sorts an array of da_iskv_t in decreasing order */
/*************************************************************************/
void da_iskvsortd(size_t n, da_iskv_t *base)
{
#define da_iskv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_iskv_t, base, n, da_iskv_gt);
#undef da_iskv_gt
}






/* da_u ## kv pairs */


/*************************************************************************/
/*! Sorts an array of da_upkv_t in increasing order */
/*************************************************************************/
void da_upkvsorti(size_t n, da_upkv_t *base)
{
#define da_upkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_upkv_t, base, n, da_upkv_lt);
#undef da_upkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_upkv_t in decreasing order */
/*************************************************************************/
void da_upkvsortd(size_t n, da_upkv_t *base)
{
#define da_upkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_upkv_t, base, n, da_upkv_gt);
#undef da_upkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_pikv_t in increasing order */
/*************************************************************************/
void da_uikvsorti(size_t n, da_uikv_t *base)
{
#define da_uikv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_uikv_t, base, n, da_uikv_lt);
#undef da_uikv_lt
}


/*************************************************************************/
/*! Sorts an array of da_uikv_t in decreasing order */
/*************************************************************************/
void da_uikvsortd(size_t n, da_uikv_t *base)
{
#define da_uikv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_uikv_t, base, n, da_uikv_gt);
#undef da_uikv_gt
}



/*************************************************************************/
/*! Sorts an array of da_uvkv_t in increasing order */
/*************************************************************************/
void da_uvkvsorti(size_t n, da_uvkv_t *base)
{
#define da_uvkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_uvkv_t, base, n, da_uvkv_lt);
#undef da_uvkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_uvkv_t in decreasing order */
/*************************************************************************/
void da_uvkvsortd(size_t n, da_uvkv_t *base)
{
#define da_uvkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_uvkv_t, base, n, da_uvkv_gt);
#undef da_uvkv_gt
}



/*************************************************************************/
/*! Sorts an array of da_uckv_t in increasing order */
/*************************************************************************/
void da_uckvsorti(size_t n, da_uckv_t *base)
{
#define da_uckv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_uckv_t, base, n, da_uckv_lt);
#undef da_uckv_lt
}


/*************************************************************************/
/*! Sorts an array of da_uckv_t in decreasing order */
/*************************************************************************/
void da_uckvsortd(size_t n, da_uckv_t *base)
{
#define da_uckv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_uckv_t, base, n, da_uckv_gt);
#undef da_pkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_ui32kv_t in increasing order */
/*************************************************************************/
void da_ui32kvsorti(size_t n, da_ui32kv_t *base)
{
#define da_ui32kv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_ui32kv_t, base, n, da_ui32kv_lt);
#undef da_ui32kv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ui32kv_t in decreasing order */
/*************************************************************************/
void da_ui32kvsortd(size_t n, da_ui32kv_t *base)
{
#define da_ui32kv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_ui32kv_t, base, n, da_ui32kv_gt);
#undef da_ui32kv_gt
}




/*************************************************************************/
/*! Sorts an array of da_ui64kv_t in increasing order */
/*************************************************************************/
void da_ui64kvsorti(size_t n, da_ui64kv_t *base)
{
#define da_ui64kv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_ui64kv_t, base, n, da_ui64kv_lt);
#undef da_ui64kv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ui64kv_t in decreasing order */
/*************************************************************************/
void da_ui64kvsortd(size_t n, da_ui64kv_t *base)
{
#define da_ui64kv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_ui64kv_t, base, n, da_ui64kv_gt);
#undef da_ui64kv_gt
}




/*************************************************************************/
/*! Sorts an array of da_uzkv_t in increasing order */
/*************************************************************************/
void da_uzkvsorti(size_t n, da_uzkv_t *base)
{
#define da_uzkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_uzkv_t, base, n, da_uzkv_lt);
#undef da_uzkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_uzkv_t in decreasing order */
/*************************************************************************/
void da_uzkvsortd(size_t n, da_uzkv_t *base)
{
#define da_uzkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_uzkv_t, base, n, da_uzkv_gt);
#undef da_uzkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_ufkv_t in increasing order */
/*************************************************************************/
void da_ufkvsorti(size_t n, da_ufkv_t *base)
{
#define da_ufkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_ufkv_t, base, n, da_ufkv_lt);
#undef da_ufkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_ufkv_t in decreasing order */
/*************************************************************************/
void da_ufkvsortd(size_t n, da_ufkv_t *base)
{
#define da_ufkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_ufkv_t, base, n, da_ufkv_gt);
#undef da_ufkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_udkv_t in increasing order */
/*************************************************************************/
void da_udkvsorti(size_t n, da_udkv_t *base)
{
#define da_udkv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_udkv_t, base, n, da_udkv_lt);
#undef da_udkv_lt
}


/*************************************************************************/
/*! Sorts an array of da_udkv_t in decreasing order */
/*************************************************************************/
void da_udkvsortd(size_t n, da_udkv_t *base)
{
#define da_udkv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_udkv_t, base, n, da_udkv_gt);
#undef da_udkv_gt
}




/*************************************************************************/
/*! Sorts an array of da_uskv_t in increasing order */
/*************************************************************************/
void da_uskvsorti(size_t n, da_uskv_t *base)
{
#define da_uskv_lt(a, b) ((a)->val < (b)->val)
  GK_MKQSORT(da_uskv_t, base, n, da_uskv_lt);
#undef da_uskv_lt
}


/*************************************************************************/
/*! Sorts an array of da_uskv_t in decreasing order */
/*************************************************************************/
void da_uskvsortd(size_t n, da_uskv_t *base)
{
#define da_uskv_gt(a, b) ((a)->val > (b)->val)
  GK_MKQSORT(da_uskv_t, base, n, da_uskv_gt);
#undef da_uskv_gt
}



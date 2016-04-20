/*!
\file  select.c
\brief Sorts only the largest k values

\date   Started 2/28/2013
\author David C. Anastasiu
 */


#include "includes.h"

/* Byte-wise swap two items of size SIZE. */
#define DA_QSSWAP(a, b, stmp) do { stmp = (a); (a) = (b); (b) = stmp; } while (0)

/**** kselect functions for ux type kv arrays  ****/

/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t da_uvkvkselectd(size_t n, idx_t topk, da_uvkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_uvkv_t stmp;
	val_t pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val < cand[mid].val)
			mid = lo;
		if (cand[hi].val > cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val < cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val >= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t da_uvkvkselecti(size_t n, idx_t topk, da_uvkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_uvkv_t stmp;
	val_t pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val > cand[mid].val)
			mid = lo;
		if (cand[hi].val < cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val > cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val <= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}





/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t da_uikvkselectd(size_t n, idx_t topk, da_uikv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_uikv_t stmp;
	idx_t pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val < cand[mid].val)
			mid = lo;
		if (cand[hi].val > cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val < cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val >= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t da_uikvkselecti(size_t n, idx_t topk, da_uikv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_uikv_t stmp;
	idx_t pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val > cand[mid].val)
			mid = lo;
		if (cand[hi].val < cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val > cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val <= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}




/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t da_ufkvkselectd(size_t n, idx_t topk, da_ufkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_ufkv_t stmp;
	float pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val < cand[mid].val)
			mid = lo;
		if (cand[hi].val > cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val < cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val >= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t da_ufkvkselecti(size_t n, idx_t topk, da_ufkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_ufkv_t stmp;
	float pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val > cand[mid].val)
			mid = lo;
		if (cand[hi].val < cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val > cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val <= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}




/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t da_udkvkselectd(size_t n, idx_t topk, da_udkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_udkv_t stmp;
	double pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val < cand[mid].val)
			mid = lo;
		if (cand[hi].val > cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val < cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val >= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t da_udkvkselecti(size_t n, idx_t topk, da_udkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_udkv_t stmp;
	double pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val > cand[mid].val)
			mid = lo;
		if (cand[hi].val < cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val > cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val <= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}




/**** kselect functions for px type kv arrays  ****/

/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t da_pvkvkselectd(size_t n, idx_t topk, da_pvkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_pvkv_t stmp;
	val_t pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val < cand[mid].val)
			mid = lo;
		if (cand[hi].val > cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val < cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val >= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t da_pvkvkselecti(size_t n, idx_t topk, da_pvkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_pvkv_t stmp;
	val_t pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val > cand[mid].val)
			mid = lo;
		if (cand[hi].val < cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val > cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val <= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}





/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t da_pikvkselectd(size_t n, idx_t topk, da_pikv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_pikv_t stmp;
	idx_t pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val < cand[mid].val)
			mid = lo;
		if (cand[hi].val > cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val < cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val >= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t da_pikvkselecti(size_t n, idx_t topk, da_pikv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_pikv_t stmp;
	idx_t pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val > cand[mid].val)
			mid = lo;
		if (cand[hi].val < cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val > cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val <= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}




/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t da_pfkvkselectd(size_t n, idx_t topk, da_pfkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_pfkv_t stmp;
	float pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val < cand[mid].val)
			mid = lo;
		if (cand[hi].val > cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val < cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val >= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t da_pfkvkselecti(size_t n, idx_t topk, da_pfkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_pfkv_t stmp;
	float pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val > cand[mid].val)
			mid = lo;
		if (cand[hi].val < cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val > cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val <= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}




/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t da_pdkvkselectd(size_t n, idx_t topk, da_pdkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_pdkv_t stmp;
	double pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val < cand[mid].val)
			mid = lo;
		if (cand[hi].val > cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val < cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val >= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t da_pdkvkselecti(size_t n, idx_t topk, da_pdkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_pdkv_t stmp;
	double pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val > cand[mid].val)
			mid = lo;
		if (cand[hi].val < cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val > cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val <= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}






/**** kselect functions for ix type kv arrays  ****/

/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t da_ivkvkselectd(size_t n, idx_t topk, da_ivkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_ivkv_t stmp;
	val_t pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val < cand[mid].val)
			mid = lo;
		if (cand[hi].val > cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val < cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val >= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t da_ivkvkselecti(size_t n, idx_t topk, da_ivkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_ivkv_t stmp;
	val_t pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val > cand[mid].val)
			mid = lo;
		if (cand[hi].val < cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val > cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val <= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}





/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t da_iikvkselectd(size_t n, idx_t topk, da_iikv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_iikv_t stmp;
	idx_t pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val < cand[mid].val)
			mid = lo;
		if (cand[hi].val > cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val < cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val >= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t da_iikvkselecti(size_t n, idx_t topk, da_iikv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_iikv_t stmp;
	idx_t pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val > cand[mid].val)
			mid = lo;
		if (cand[hi].val < cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val > cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val <= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}




/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t da_ifkvkselectd(size_t n, idx_t topk, da_ifkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_ifkv_t stmp;
	float pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val < cand[mid].val)
			mid = lo;
		if (cand[hi].val > cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val < cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val >= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t da_ifkvkselecti(size_t n, idx_t topk, da_ifkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_ifkv_t stmp;
	float pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val > cand[mid].val)
			mid = lo;
		if (cand[hi].val < cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val > cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val <= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}




/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
idx_t da_idkvkselectd(size_t n, idx_t topk, da_idkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_idkv_t stmp;
	double pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val < cand[mid].val)
			mid = lo;
		if (cand[hi].val > cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val < cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val >= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}


/******************************************************************************/
/*! This function puts the 'topk' smallest values in the beginning of the array */
/*******************************************************************************/
idx_t da_idkvkselecti(size_t n, idx_t topk, da_idkv_t *cand)
{
	idx_t i, j, lo, hi, mid;
	da_idkv_t stmp;
	double pivot;

	if (n <= topk)
		return n; /* return if the array has fewer elements than we want */

	for (lo=0, hi=n-1; lo < hi;) {
		mid = lo + ((hi-lo) >> 1);

		/* select the median */
		if (cand[lo].val > cand[mid].val)
			mid = lo;
		if (cand[hi].val < cand[mid].val)
			mid = hi;
		else
			goto jump_over;
		if (cand[lo].val > cand[mid].val)
			mid = lo;

		jump_over:
		DA_QSSWAP(cand[mid], cand[hi], stmp);
		pivot = cand[hi].val;

		/* the partitioning algorithm */
		for (i=lo-1, j=lo; j<hi; j++) {
			if (cand[j].val <= pivot) {
				i++;
				DA_QSSWAP(cand[i], cand[j], stmp);
			}
		}
		i++;
		DA_QSSWAP(cand[i], cand[hi], stmp);


		if (i > topk)
			hi = i-1;
		else if (i < topk)
			lo = i+1;
		else
			break;
	}

	return topk;
}








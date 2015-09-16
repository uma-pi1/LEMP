/*!
 \file  csr.c
 \brief Functions for dealing with csr structures

 \author David C. Anastasiu
 */
#include "includes.h"


#define DA_OMPMINOPS       50000

//forward declaration
void da_csr_WriteClean(int32_t **bptr, int32_t **bind, float **bval, char * filename, FILE *fpout);

/*************************************************************************/
/*! Allocate memory for a CSR matrix and initializes it
    \returns the allocated matrix. The various fields are set to NULL.
 */
/**************************************************************************/
da_csr_t *da_csr_Create()
{
	da_csr_t *mat;

	mat = (da_csr_t *)gk_malloc(sizeof(da_csr_t), "da_csr_Create: mat");

	da_csr_Init(mat);

	return mat;
}


/*************************************************************************/
/*! Initializes the matrix
    \param mat is the matrix to be initialized.
 */
/*************************************************************************/
void da_csr_Init(da_csr_t *mat)
{
	memset(mat, 0, sizeof(da_csr_t));
	mat->nrows = mat->ncols = -1;
}


/*************************************************************************/
/*! Frees all the memory allocated for matrix.
    \param mat is the matrix to be freed.
 */
/*************************************************************************/
void da_csr_Free(da_csr_t **mat)
{
	if (*mat == NULL)
		return;
	da_csr_FreeContents(*mat);
	gk_free((void **)mat, LTERM);
}


/*************************************************************************/
/*! Frees a variable sized list of matrix pointers. Last item in the list must be LTERM.
    \param ptr1 is the first matrix to be freed,
    \param ... are additional matrices to be freed.
 */
/*************************************************************************/
void da_csr_FreeAll(da_csr_t **ptr1,...)
{
	va_list plist;
	void **ptr;

	if (*ptr1 != NULL)
		da_csr_Free(ptr1);

	va_start(plist, ptr1);
	while ((ptr = va_arg(plist, void **)) != LTERM) {
		if (*ptr != NULL) {
			da_csr_Free((da_csr_t **)ptr);
		}
	}
	va_end(plist);
}


/*************************************************************************/
/*! Frees the memory allocated for the matrix's fields for the given type and
    sets them to NULL.
    \param mat is the matrix whose partial contents will be freed,
    \param type is the internal representation to be freed: row (DA_ROW),
    		or column (DA_COL).
 */
/*************************************************************************/
void da_csr_FreeBase(da_csr_t *mat, char type)
{
	if(type & DA_ROW)
		//drop the row index
		gk_free((void **)&mat->rowptr, &mat->rowind, &mat->rowval, LTERM);

	if(type & DA_COL)
		//drop the column index
		gk_free((void **)&mat->colptr, &mat->colind, &mat->colval, LTERM);

}


/*************************************************************************/
/*! Create missing indexes for the matrix such that both indexes exist.
    \param mat is the matrix whose bases need loading.
 */
/*************************************************************************/
void da_csr_LoadBases(da_csr_t *csr)
{
	if(csr->rowptr && csr->colptr) return;
	if(!csr->rowptr && !csr->colptr) return;
	if(csr->rowptr){
		da_csr_CreateIndex(csr, DA_COL);
		// sort column indices for input matrices
		da_csr_SortIndices(csr, DA_COL);
	} else {
		da_csr_CreateIndex(csr, DA_ROW);
		// sort column indices for input matrices
		da_csr_SortIndices(csr, DA_ROW);
	}
}



/*************************************************************************/
/*! Frees only the memory allocated for the matrix's different fields and
    sets them to NULL.
    \param mat is the matrix whose contents will be freed.
 */
/*************************************************************************/
void da_csr_FreeContents(da_csr_t *mat)
{
	gk_free((void *)&mat->rowptr, &mat->rowind, &mat->rowval, &mat->rowids,
			&mat->colptr, &mat->colind, &mat->colval, &mat->colids,
			&mat->rnorms, &mat->cnorms, &mat->rsums, &mat->csums,
			&mat->rsizes, &mat->csizes, &mat->rvols, &mat->cvols,
			&mat->rwgts, &mat->cwgts,
			LTERM);
}


/*************************************************************************/
/*! Transfers the memory within from to the to matrix.
    \params from is the matrix from which we transfer memory. Its pointers
    	will be set to NULL.
    \param to is the matrix whose contents will be freed that will receive
    	memory from the from matrix.
 */
/*************************************************************************/
void da_csr_Transfer(da_csr_t *from, da_csr_t *to)
{
	da_csr_FreeContents(to);
	to->rowptr = from->rowptr;
	to->rowind = from->rowind;
	to->rowval = from->rowval;
	to->rowids = from->rowids;
	to->colptr = from->colptr;
	to->colind = from->colind;
	to->colval = from->colval;
	to->colids = from->colids;
	to->rnorms = from->rnorms;
	to->cnorms = from->cnorms;
	to->rsums = from->rsums;
	to->csums = from->csums;
	to->rsizes = from->rsizes;
	to->csizes = from->csizes;
	to->rvols = from->rvols;
	to->cvols = from->cvols;
	to->rwgts = from->rwgts;
	to->cwgts = from->cwgts;

	to->nrows = from->nrows;
	to->ncols = from->ncols;

	from->rowptr = NULL;
	from->rowind = NULL;
	from->rowval = NULL;
	from->rowids = NULL;
	from->colptr = NULL;
	from->colind = NULL;
	from->colval = NULL;
	from->colids = NULL;
	from->rnorms = NULL;
	from->cnorms = NULL;
	from->rsums = NULL;
	from->csums = NULL;
	from->rsizes = NULL;
	from->csizes = NULL;
	from->rvols = NULL;
	from->cvols = NULL;
	from->rwgts = NULL;
	from->cwgts = NULL;
}


/*************************************************************************/
/*! Returns a copy of a matrix.
    \param mat is the matrix to be duplicated.
    \returns the newly created copy of the matrix.
 */
/**************************************************************************/
da_csr_t *da_csr_Copy(da_csr_t *mat)
{
	da_csr_t *nmat;

	nmat = da_csr_Create();

	nmat->nrows  = mat->nrows;
	nmat->ncols  = mat->ncols;

	/* copy the row structure */
	if (mat->rowptr)
		nmat->rowptr = da_pcopy(mat->nrows+1, mat->rowptr,
				da_pmalloc(mat->nrows+1, "da_csr_Dup: rowptr"));
	if (mat->rowids)
		nmat->rowids = da_icopy(mat->nrows, mat->rowids,
				da_imalloc(mat->nrows, "da_csr_Dup: rowids"));
	if (mat->rnorms)
		nmat->rnorms = da_vcopy(mat->nrows, mat->rnorms,
				da_vmalloc(mat->nrows, "da_csr_Dup: rnorms"));
	if (mat->rowind)
		nmat->rowind = da_icopy(mat->rowptr[mat->nrows], mat->rowind,
				da_imalloc(mat->rowptr[mat->nrows], "da_csr_Dup: rowind"));
	if (mat->rowval)
		nmat->rowval = da_vcopy(mat->rowptr[mat->nrows], mat->rowval,
				da_vmalloc(mat->rowptr[mat->nrows], "da_csr_Dup: rowval"));

	/* copy the col structure */
	if (mat->colptr)
		nmat->colptr = da_pcopy(mat->ncols+1, mat->colptr,
				da_pmalloc(mat->ncols+1, "da_csr_Dup: colptr"));
	if (mat->colids)
		nmat->colids = da_icopy(mat->ncols, mat->colids,
				da_imalloc(mat->ncols, "da_csr_Dup: colids"));
	if (mat->cnorms)
		nmat->cnorms = da_vcopy(mat->ncols, mat->cnorms,
				da_vmalloc(mat->ncols, "da_csr_Dup: cnorms"));
	if (mat->colind)
		nmat->colind = da_icopy(mat->colptr[mat->ncols], mat->colind,
				da_imalloc(mat->colptr[mat->ncols], "da_csr_Dup: colind"));
	if (mat->colval)
		nmat->colval = da_vcopy(mat->colptr[mat->ncols], mat->colval,
				da_vmalloc(mat->colptr[mat->ncols], "da_csr_Dup: colval"));

	return nmat;
}


/*************************************************************************/
/*! Returns a submatrix containing a set of consecutive rows.
    \param mat is the original matrix.
    \param rstart is the starting row.
    \param nrows is the number of rows from rstart to extract.
    \returns the row structure of the newly created submatrix.
 */
/**************************************************************************/
da_csr_t *da_csr_ExtractSubmatrix(da_csr_t *mat, idx_t rstart, idx_t nrows)
{
	ssize_t i;
	da_csr_t *nmat;

	if (rstart+nrows > mat->nrows)
		return NULL;

	nmat = da_csr_Create();

	nmat->nrows  = nrows;
	nmat->ncols  = mat->ncols;

	/* copy the row structure */
	if (mat->rowptr)
		nmat->rowptr = da_pcopy(nrows+1, mat->rowptr+rstart,
				da_pmalloc(nrows+1, "da_csr_ExtractSubmatrix: rowptr"));
	for (i=nrows; i>=0; i--)
		nmat->rowptr[i] -= nmat->rowptr[0];
	ASSERT(nmat->rowptr[0] == 0);

	if (mat->rowids)
		nmat->rowids = da_icopy(nrows, mat->rowids+rstart,
				da_imalloc(nrows, "da_csr_ExtractSubmatrix: rowids"));
	if (mat->rnorms)
		nmat->rnorms = da_vcopy(nrows, mat->rnorms+rstart,
				da_vmalloc(nrows, "da_csr_ExtractSubmatrix: rnorms"));

	if (mat->rsums)
		nmat->rsums = da_vcopy(nrows, mat->rsums+rstart,
				da_vmalloc(nrows, "da_csr_ExtractSubmatrix: rsums"));

	ASSERT(nmat->rowptr[nrows] == mat->rowptr[rstart+nrows]-mat->rowptr[rstart]);
	if (mat->rowind)
		nmat->rowind = da_icopy(mat->rowptr[rstart+nrows]-mat->rowptr[rstart],
				mat->rowind+mat->rowptr[rstart],
				da_imalloc(mat->rowptr[rstart+nrows]-mat->rowptr[rstart],
						"da_csr_ExtractSubmatrix: rowind"));
	if (mat->rowval)
		nmat->rowval = da_vcopy(mat->rowptr[rstart+nrows]-mat->rowptr[rstart],
				mat->rowval+mat->rowptr[rstart],
				da_vmalloc(mat->rowptr[rstart+nrows]-mat->rowptr[rstart],
						"da_csr_ExtractSubmatrix: rowval"));

	return nmat;
}


/*************************************************************************/
/*! Returns a submatrix containing a certain set of rows.
    \param mat is the original matrix.
    \param nrows is the number of rows to extract.
    \param rind is the set of row numbers to extract.
    \returns the row structure of the newly created submatrix.
 */
/**************************************************************************/
da_csr_t *da_csr_ExtractRows(da_csr_t *mat, idx_t nrows, idx_t *rind)
{
	ssize_t i, ii, j, nnz;
	da_csr_t *nmat;

	nmat = da_csr_Create();

	nmat->nrows = nrows;
	nmat->ncols = mat->ncols;

	for (nnz=0, i=0; i<nrows; i++)
		nnz += mat->rowptr[rind[i]+1]-mat->rowptr[rind[i]];

	nmat->rowptr = da_pmalloc(nmat->nrows+1, "da_csr_ExtractPartition: rowptr");
	nmat->rowind = da_imalloc(nnz, "da_csr_ExtractPartition: rowind");
	nmat->rowval = da_vmalloc(nnz, "da_csr_ExtractPartition: rowval");

	nmat->rowptr[0] = 0;
	for (nnz=0, j=0, ii=0; ii<nrows; ii++) {
		i = rind[ii];
		da_icopy(mat->rowptr[i+1]-mat->rowptr[i], mat->rowind+mat->rowptr[i], nmat->rowind+nnz);
		da_vcopy(mat->rowptr[i+1]-mat->rowptr[i], mat->rowval+mat->rowptr[i], nmat->rowval+nnz);
		nnz += mat->rowptr[i+1]-mat->rowptr[i];
		nmat->rowptr[++j] = nnz;
	}
	ASSERT(j == nmat->nrows);

	return nmat;
}


/*************************************************************************/
/*! Returns a submatrix containing a certain set of rows.
    \param mat is the original matrix.
    \param nrows is the number of rows to extract.
    \param rind is the set of row numbers to extract.
    \returns the row structure of the newly created submatrix.
 */
/**************************************************************************/
void da_csr_ExtractRowsInto(da_csr_t *mat, da_csr_t *nmat, idx_t nrows, idx_t *rind)
{
	ssize_t i, ii, j, nnz;

	for (nnz=0, i=0; i<nrows; i++)
		nnz += mat->rowptr[rind[i]+1]-mat->rowptr[rind[i]];
	ASSERT(nmat->rowptr[nmat->nrows] >= nnz)

	nmat->ncols = mat->ncols;
	nmat->nrows = nrows;
	nmat->rowptr[0] = 0;
	for (nnz=0, j=0, ii=0; ii<nrows; ii++) {
		i = rind[ii];
		da_icopy(mat->rowptr[i+1]-mat->rowptr[i], mat->rowind+mat->rowptr[i], nmat->rowind+nnz);
		da_vcopy(mat->rowptr[i+1]-mat->rowptr[i], mat->rowval+mat->rowptr[i], nmat->rowval+nnz);
		nnz += mat->rowptr[i+1]-mat->rowptr[i];
		nmat->rowptr[++j] = nnz;
	}
	ASSERT(j == nmat->nrows);

}


/*************************************************************************/
/*! Returns a submatrix corresponding to a specified partitioning of rows.
    \param mat is the original matrix.
    \param part is the partitioning vector of the rows.
    \param pid is the partition ID that will be extracted.
    \returns the row structure of the newly created submatrix.
 */
/**************************************************************************/
da_csr_t *da_csr_ExtractPartition(da_csr_t *mat, idx_t *part, idx_t pid)
{
	ssize_t i, j, nnz;
	da_csr_t *nmat;

	nmat = da_csr_Create();

	nmat->nrows = 0;
	nmat->ncols = mat->ncols;

	for (nnz=0, i=0; i<mat->nrows; i++) {
		if (part[i] == pid) {
			nmat->nrows++;
			nnz += mat->rowptr[i+1]-mat->rowptr[i];
		}
	}

	nmat->rowptr = da_pmalloc(nmat->nrows+1, "da_csr_ExtractPartition: rowptr");
	nmat->rowind = da_imalloc(nnz, "da_csr_ExtractPartition: rowind");
	nmat->rowval = da_vmalloc(nnz, "da_csr_ExtractPartition: rowval");

	nmat->rowptr[0] = 0;
	for (nnz=0, j=0, i=0; i<mat->nrows; i++) {
		if (part[i] == pid) {
			da_icopy(mat->rowptr[i+1]-mat->rowptr[i], mat->rowind+mat->rowptr[i], nmat->rowind+nnz);
			da_vcopy(mat->rowptr[i+1]-mat->rowptr[i], mat->rowval+mat->rowptr[i], nmat->rowval+nnz);
			nnz += mat->rowptr[i+1]-mat->rowptr[i];
			nmat->rowptr[++j] = nnz;
		}
	}
	ASSERT(j == nmat->nrows);

	return nmat;
}


/*************************************************************************/
/*! Splits the matrix into multiple sub-matrices based on the provided
    color array.
    \param mat is the original matrix.
    \param color is an array of size equal to the number of non-zeros
           in the matrix (row-wise structure). The matrix is split into
           as many parts as the number of colors. For meaningfull results,
           the colors should be numbered consecutively starting from 0.
    \returns an array of matrices for each supplied color number.
 */
/**************************************************************************/
da_csr_t **da_csr_Split(da_csr_t *mat, idx_t *color)
{
	ssize_t i, j;
	idx_t nrows, ncolors;
	ptr_t *rowptr;
	idx_t *rowind;
	val_t *rowval;
	da_csr_t **smats;

	nrows  = mat->nrows;
	rowptr = mat->rowptr;
	rowind = mat->rowind;
	rowval = mat->rowval;

	ncolors = da_imax(rowptr[nrows], color)+1;

	smats = (da_csr_t **)gk_malloc(sizeof(da_csr_t *)*ncolors, "da_csr_Split: smats");
	for (i=0; i<ncolors; i++) {
		smats[i] = da_csr_Create();
		smats[i]->nrows  = mat->nrows;
		smats[i]->ncols  = mat->ncols;
		smats[i]->rowptr = da_psmalloc(nrows+1, 0, "da_csr_Split: smats[i]->rowptr");
	}

	for (i=0; i<nrows; i++) {
		for (j=rowptr[i]; j<rowptr[i+1]; j++)
			smats[color[j]]->rowptr[i]++;
	}
	for (i=0; i<ncolors; i++)
		MAKECSR(j, nrows, smats[i]->rowptr);

	for (i=0; i<ncolors; i++) {
		smats[i]->rowind = da_imalloc(smats[i]->rowptr[nrows], "da_csr_Split: smats[i]->rowind");
		smats[i]->rowval = da_vmalloc(smats[i]->rowptr[nrows], "da_csr_Split: smats[i]->rowval");
	}

	for (i=0; i<nrows; i++) {
		for (j=rowptr[i]; j<rowptr[i+1]; j++) {
			smats[color[j]]->rowind[smats[color[j]]->rowptr[i]] = rowind[j];
			smats[color[j]]->rowval[smats[color[j]]->rowptr[i]] = rowval[j];
			smats[color[j]]->rowptr[i]++;
		}
	}

	for (i=0; i<ncolors; i++)
		SHIFTCSR(j, nrows, smats[i]->rowptr);

	return smats;
}


/**************************************************************************/
/*! Reads a CSR matrix from the supplied file and stores it the matrix's
    forward structure.
    \param filename is the file that stores the data.
    \param format is either DA_FMT_METIS, DA_FMT_CLUTO,
           DA_FMT_CSR, DA_FMT_BINROW, DA_FMT_BINCOL
           specifying the type of the input format.
           The DA_FMT_CSR does not contain a header
           line, whereas the DA_FMT_BINROW is a binary format written
           by da_csr_Write() using the same format specifier.
    \param readvals is either 1 or 0, indicating if the CSR file contains
           values or it does not. It only applies when DA_FMT_CSR is
           used.
    \param numbering is either 1 or 0, indicating if the numbering of the
           indices start from 1 or 0, respectively. If they start from 1,
           they are automatically decreamented during input so that they
           will start from 0. It only applies when DA_FMT_CSR is
           used.
    \returns the matrix that was read.
 */
/**************************************************************************/
da_csr_t *da_csr_Read(char *filename, char format, char readvals, char numbering)
{
	ssize_t i, j, k, l;
	size_t nrows, ncols, nnz, nfields, fmt, ncon, lnlen, read, size;
	ptr_t *rowptr;
	idx_t *rowind;
	int32_t *iinds, *jinds, *bptr = NULL, *bind = NULL, rid, len, nnz2, nr, nc, nz, rv;
	val_t *rowval = NULL, *vals, fval;
	float *bval = NULL;
	char readsizes, readwgts;
	char *line = NULL, *head, *tail, fmtstr[256];
	FILE *fpin;
	da_csr_t *mat = NULL;

	if (!gk_fexists(filename))
		gk_errexit(SIGERR, "File %s does not exist!\n", filename);

	if (format == DA_FMT_BINROW) {
		mat = da_csr_Create();

		fpin = gk_fopen(filename, "rb", "da_csr_Read: fpin");
		if (fread(&(nr), sizeof(int32_t), 1, fpin) != 1)
			gk_errexit(SIGERR, "Failed to read the nrows from file %s!\n", filename);
		if (fread(&(nc), sizeof(int32_t), 1, fpin) != 1)
			gk_errexit(SIGERR, "Failed to read the ncols from file %s!\n", filename);
		mat->nrows = nr;
		mat->ncols = nc;
		mat->rowptr = da_pmalloc(nr+1, "da_csr_Read: rowptr");
		if(sizeof(ptr_t) != sizeof(int32_t)){
			bptr = da_i32malloc(nr+1, "da_csr_Read: bptr");
			if (fread(bptr, sizeof(int32_t), nr+1, fpin) != nr+1)
				gk_errexit(SIGERR, "Failed to read the rowptr from file %s!\n", filename);
			for(i=0; i < nr+1; i++)
				mat->rowptr[i] = bptr[i];
			gk_free((void**)&bptr, LTERM);
		} else if (fread(mat->rowptr, sizeof(ptr_t), nr+1, fpin) != nr+1)
			gk_errexit(SIGERR, "Failed to read the rowptr from file %s!\n", filename);

		nnz = mat->rowptr[mat->nrows];
		mat->rowind = da_imalloc(nnz, "da_csr_Read: rowind");
		if(sizeof(idx_t) != sizeof(int32_t)){
			bind = da_i32malloc(nnz, "da_csr_Read: bind");
			if (fread(bind, sizeof(int32_t), nnz, fpin) != nnz)
				gk_errexit(SIGERR, "Failed to read the rowind from file %s!\n", filename);
			for(i=0; i < nnz; i++)
				mat->rowind[i] = bind[i];
			gk_free((void**)&bind, LTERM);
		} else if (fread(mat->rowind, sizeof(idx_t), nnz, fpin) != nnz)
			gk_errexit(SIGERR, "Failed to read the rowind from file %s!\n", filename);

		if (readvals == 1) {
			mat->rowval = da_vmalloc(nnz, "da_csr_Read: rowval");
			if(sizeof(val_t) != sizeof(float)){
				bval = da_fmalloc(nnz, "da_csr_Read: bval");
				if (fread(bind, sizeof(float), nnz, fpin) != nnz)
					gk_errexit(SIGERR, "Failed to read the rowind from file %s!\n", filename);
				for(i=0; i < nnz; i++)
					mat->rowval[i] = bval[i];
				gk_free((void**)&bval, LTERM);
			} else if (fread(mat->rowval, sizeof(val_t), nnz, fpin) != nnz)
				gk_errexit(SIGERR, "Failed to read the rowval from file %s!\n", filename);
		}

		gk_fclose(fpin);
		return mat;
	}

	if (format == DA_FMT_BINAP) {
		mat = da_csr_Create();

		fpin = gk_fopen(filename, "rb", "da_csr_Read: fpin");

		if(fread(&nr, sizeof(int32_t), 1, fpin) != 1)
			gk_errexit(SIGERR, "Could not read the number of rows...\n");
		mat->nrows = nr;

		struct stat st;
		stat(filename, &st);
		long size = st.st_size;
		nnz = (size - sizeof(int32_t) * mat->nrows - 1)
				/ (sizeof(int32_t) + sizeof(float));

		mat->rowptr = da_pmalloc(nr+1, "da_csr_Read: rowptr");
		if(sizeof(ptr_t) != sizeof(int32_t))
			bptr = da_i32malloc(nr+1, "da_csr_Read: bptr");
		mat->rowind = da_imalloc(nnz, "da_csr_Read: rowind");
		if(sizeof(idx_t) != sizeof(int32_t))
			bind = da_i32malloc(nnz, "da_csr_Read: bind");
		if (readvals == 1){
			mat->rowval = da_vmalloc(nnz, "da_csr_Read: rowval");
			if(sizeof(val_t) != sizeof(float))
				bval = da_fmalloc(nnz, "da_csr_Read: bval");
		}

		mat->rowptr[0] = 0;
		for ( i=0; i < nr; i++ )
		{
			if(fread(&len, sizeof(int32_t), 1, fpin) != 1)
				gk_errexit(SIGERR, "Failed to read the number of elements in row %d "
						"from file %s!\n", i+1, filename);

			if(sizeof(idx_t) != sizeof(int32_t)){
				if(fread( bind, sizeof(int32_t), len, fpin) != len)
					gk_errexit(SIGERR, "Failed to read the indicator elements in row %d "
						"from file %s!\n", i+1, filename);
				for(j=0; j < len; j++)
					mat->rowind[mat->rowptr[i]+j] = bind[j];
			} else if(fread( mat->rowind+mat->rowptr[i], sizeof(idx_t), len, fpin) != len)
				gk_errexit(SIGERR, "Failed to read the indicator elements in row %d "
					"from file %s!\n", i+1, filename);

			if(sizeof(val_t) != sizeof(float)){
				if(fread( bval, sizeof(float), len, fpin) != len)
					gk_errexit(SIGERR, "Failed to read the value elements in row %d "
						"from file %s!\n", i+1, filename);
				for(j=0; j < len; j++)
					mat->rowval[mat->rowptr[i]+j] = bval[j];
			} else if(fread( mat->rowval+mat->rowptr[i], sizeof(val_t), len, fpin) != len)
				gk_errexit(SIGERR, "Failed to read the value elements in row %d "
					"from file %s!\n", i+1, filename);

			mat->rowptr[i+1] = mat->rowptr[i] + len;
		}

		if ( mat->rowptr[mat->nrows] != nnz )
		{
			nnz = mat->rowptr[mat->nrows];
			mat->rowind = da_irealloc(mat->rowind, nnz, "da_csr_Read: rowind resize");
			mat->rowval = da_vrealloc(mat->rowval, nnz, "da_csr_Read: rowval resize");
		}

		for ( i=0, ncols=0; i<nnz; i++ ){
			if(mat->rowind[i] > ncols)
				ncols = mat->rowind[i];
			mat->rowind[i]--;
		}
		mat->ncols = ncols;

		da_csr_SortIndices(mat, DA_ROW);

		gk_fclose(fpin);
		gk_free((void**)&bptr, &bind, &bval, LTERM);
		return mat;
	}

	if (format == DA_FMT_BINAPB) {
		mat = da_csr_Create();

		fpin = gk_fopen(filename, "rb", "da_csr_Read: fpin");

		struct stat st;
		stat(filename, &st);
		long size = st.st_size;
		nnz = size/(sizeof(int32_t));

		// Going to make the pessimistic assumption that there are as
		// many records as there are nnz
		mat->rowptr = da_pmalloc(nnz+1, "da_csr_Read: rowptr");
		mat->rowind = da_imalloc(nnz, "da_csr_Read: rowind");
		mat->nrows = 0;
		mat->rowptr[0] = 0;
		mat->rowval = NULL;

		nnz2 = 0;
		while ( fread(&rid, sizeof(int32_t), 1, fpin) == 1 ) {
			if(fread(&len, sizeof(int32_t), 1, fpin) != 1)
				gk_errexit(SIGERR, "Could not read row length in file %s!\n", filename);
			if(sizeof(idx_t) != sizeof(int32_t)){
				if(fread( bind, sizeof(int32_t), len, fpin) != len)
					gk_errexit(SIGERR, "Failed to read the indicator elements in row %d "
						"from file %s!\n", mat->nrows+1, filename);
				for(j=0; j < len; j++)
					mat->rowind[mat->rowptr[mat->nrows]+j] = bind[j];
			} else if(fread( mat->rowind+mat->rowptr[mat->nrows], sizeof(idx_t), len, fpin) != len)
				gk_errexit(SIGERR, "Failed to read the indicator elements in row %d "
					"from file %s!\n", mat->nrows+1, filename);
			nnz2 += len;
			ASSERT(nnz2 <= nnz);
			mat->rowptr[mat->nrows+1] = mat->rowptr[mat->nrows] + len;
			mat->nrows++;
			ASSERT(mat->nrows <= nnz);
		}

		mat->rowptr = da_prealloc(mat->rowptr, (mat->nrows+1), "da_csr_Read: rowptr resize");
		if ( mat->rowptr[mat->nrows] != nnz )
		{
			nnz = mat->rowptr[mat->nrows];
			mat->rowind = da_irealloc(mat->rowind, nnz, "da_csr_Read: rowind resize");
		}

		for ( i=0, ncols=0; i<nnz; i++ ){
			if(mat->rowind[i] > ncols)
				ncols = mat->rowind[i];
			mat->rowind[i]--;
		}
		mat->ncols = ncols;

		da_csr_SortIndices(mat, DA_ROW);

		gk_fclose(fpin);
		gk_free((void**) &bind, LTERM);
		return mat;
	}

	if (format == DA_FMT_BINCOL) {
		mat = da_csr_Create();

		fpin = gk_fopen(filename, "rb", "da_csr_Read: fpin");

		if (fread(&(nr), sizeof(int32_t), 1, fpin) != 1)
			gk_errexit(SIGERR, "Failed to read the nrows from file %s!\n", filename);
		if (fread(&(nc), sizeof(int32_t), 1, fpin) != 1)
			gk_errexit(SIGERR, "Failed to read the ncols from file %s!\n", filename);
		mat->nrows = nr;
		mat->ncols = nc;
		mat->colptr = da_pmalloc(nc+1, "da_csr_Read: colptr");
		if(sizeof(ptr_t) != sizeof(int32_t)){
			bptr = da_i32malloc(nc+1, "da_csr_Read: bptr");
			if (fread(bptr, sizeof(int32_t), nc+1, fpin) != nr+1)
				gk_errexit(SIGERR, "Failed to read the colptr from file %s!\n", filename);
			for(i=0; i < nc+1; i++)
				mat->colptr[i] = bptr[i];
			gk_free((void**)&bptr, LTERM);
		} else if (fread(mat->colptr, sizeof(ptr_t), nc+1, fpin) != nc+1)
			gk_errexit(SIGERR, "Failed to read the colptr from file %s!\n", filename);

		nnz = mat->colptr[mat->ncols];
		mat->colind = da_imalloc(nnz, "da_csr_Read: colind");
		if(sizeof(idx_t) != sizeof(int32_t)){
			bind = da_i32malloc(nnz, "da_csr_Read: bind");
			if (fread(bind, sizeof(int32_t), nnz, fpin) != nnz)
				gk_errexit(SIGERR, "Failed to read the colind from file %s!\n", filename);
			for(i=0; i < nnz; i++)
				mat->colind[i] = bind[i];
			gk_free((void**)&bind, LTERM);
		} else if (fread(mat->colind, sizeof(idx_t), nnz, fpin) != nnz)
			gk_errexit(SIGERR, "Failed to read the colind from file %s!\n", filename);

		if (readvals == 1) {
			mat->colval = da_vmalloc(nnz, "da_csr_Read: colval");
			if(sizeof(val_t) != sizeof(float)){
				bval = da_fmalloc(nnz, "da_csr_Read: bval");
				if (fread(bind, sizeof(float), nnz, fpin) != nnz)
					gk_errexit(SIGERR, "Failed to read the colind from file %s!\n", filename);
				for(i=0; i < nnz; i++)
					mat->colval[i] = bval[i];
				gk_free((void**)&bval, LTERM);
			} else if (fread(mat->colval, sizeof(val_t), nnz, fpin) != nnz)
				gk_errexit(SIGERR, "Failed to read the colval from file %s!\n", filename);
		}

		gk_fclose(fpin);
		return mat;
	}

	if (format == DA_FMT_BIJV) {
		mat = da_csr_Create();

		fpin = gk_fopen(filename, "rb", "da_csr_Read: fpin");
		if (fread(&(nr), sizeof(int32_t), 1, fpin) != 1)
			gk_errexit(SIGERR, "Failed to read the nrows from file %s!\n", filename);
		if (fread(&(nc), sizeof(int32_t), 1, fpin) != 1)
			gk_errexit(SIGERR, "Failed to read the ncols from file %s!\n", filename);
		if (fread(&nnz, sizeof(size_t), 1, fpin) != 1)
			gk_errexit(SIGERR, "Failed to read the nnz from file %s!\n", filename);
		if (fread(&rv, sizeof(int32_t), 1, fpin) != 1)
			gk_errexit(SIGERR, "Failed to read the readvals from file %s!\n", filename);

		if(readvals != rv)
			gk_errexit(SIGERR, "Readvals requested but file %s does not contain values!\n", filename);

		iinds = da_i32malloc(nnz, "iinds");
	    jinds = da_i32malloc(nnz, "jinds");
	    vals  = (readvals ? da_vmalloc(nnz, "vals") : NULL);
	    nrows = nr;
	    ncols = nc;

	    for (i=0; i<nnz; i++) {
			if (fread(&(iinds[i]), sizeof(int32_t), 1, fpin) != 1)
				gk_errexit(SIGERR, "Failed to read iinds[i] from file %s!\n",
						filename);
			if (fread(&(jinds[i]), sizeof(int32_t), 1, fpin) != 1)
				gk_errexit(SIGERR, "Failed to read jinds[i] from file %s!\n",
						filename);
			if (readvals) {
				if (fread(&(vals[i]), sizeof(float), 1, fpin) != 1)
					gk_errexit(SIGERR, "Failed to read vals[i] from file %s!\n",
							filename);
			}
		}
		gk_fclose(fpin);

	}

	if (format == DA_FMT_IJV) {
		gk_getfilestats(filename, &nrows, &nnz, NULL, NULL);

		if (readvals == 1 && 3 * nrows != nnz)
			gk_errexit(SIGERR,
					"Error: The number of numbers (%zd %d) in the file %s is not a multiple of 3.\n",
					nnz, readvals, filename);
		if (readvals == 0 && 2 * nrows != nnz)
			gk_errexit(SIGERR,
					"Error: The number of numbers (%zd %d) in the file %s is not a multiple of 2.\n",
					nnz, readvals, filename);

		mat = da_csr_Create();
		nnz = nrows;
		numbering = (numbering ? -1 : 0);
		iinds = da_i32malloc(nnz, "iinds");
		jinds = da_i32malloc(nnz, "jinds");
		vals  = (readvals ? da_vmalloc(nnz, "vals") : NULL);

		fpin = gk_fopen(filename, "r", "da_csr_Read: fpin");
		for (nrows = 0, ncols = 0, i = 0; i < nnz; i++) {
			if (readvals) {
				if (fscanf(fpin, "%d %d %f", &iinds[i], &jinds[i], &vals[i]) != 3)
					gk_errexit(SIGERR, "Error: Failed to read (i, j, val) for nnz: %zd.\n", i);
			} else {
				if (fscanf(fpin, "%d %d", &iinds[i], &jinds[i]) != 2)
					gk_errexit(SIGERR, "Error: Failed to read (i, j) value for nnz: %zd.\n", i);
			}
			iinds[i] += numbering;
			jinds[i] += numbering;

			if (nrows < iinds[i])
				nrows = iinds[i];
			if (ncols < jinds[i])
				ncols = jinds[i];
		}
		nrows++;
		ncols++;
		gk_fclose(fpin);

	}

	if (format == DA_FMT_IJV || format == DA_FMT_BIJV) {
		/* convert (i, j, v) into a CSR matrix */
		mat->nrows = nrows;
		mat->ncols = ncols;
		rowptr = mat->rowptr = da_psmalloc(mat->nrows + 1, 0, "rowptr");
		rowind = mat->rowind = da_imalloc(nnz, "rowind");
		if (readvals)
			rowval = mat->rowval = da_vmalloc(nnz, "rowval");

		for (i = 0; i < nnz; i++)
			rowptr[iinds[i]]++;
		MAKECSR(i, mat->nrows, rowptr);

		for (i = 0; i < nnz; i++) {
			rowind[rowptr[iinds[i]]] = jinds[i];
			if (readvals)
				rowval[rowptr[iinds[i]]] = vals[i];
			rowptr[iinds[i]]++;
		}
		SHIFTCSR(i, mat->nrows, rowptr);

		gk_free((void **) &iinds, &jinds, &vals, LTERM);

		return mat;
	}

	if (format == DA_FMT_CLUTO) {
		fpin = gk_fopen(filename, "r", "da_csr_Read: fpin");
		do {
			if (gk_getline(&line, &lnlen, fpin) <= 0)
				gk_errexit(SIGERR, "Premature end of input file: file:%s\n", filename);
		} while (line[0] == '%');

		if (sscanf(line, "%zu %zu %zu", &nrows, &ncols, &nnz) != 3)
			gk_errexit(SIGERR, "Header line must contain 3 integers.\n");

		readsizes = 0;
		readwgts  = 0;
		readvals  = 1;
		numbering = 1;
	}
	else if (format == DA_FMT_METIS) {
		fpin = gk_fopen(filename, "r", "da_csr_Read: fpin");
		do {
			if (gk_getline(&line, &lnlen, fpin) <= 0)
				gk_errexit(SIGERR, "Premature end of input file: file:%s\n", filename);
		} while (line[0] == '%');

		fmt = ncon = 0;
		nfields = sscanf(line, "%zu %zu %zu %zu", &nrows, &nnz, &fmt, &ncon);
		if (nfields < 2)
			gk_errexit(SIGERR, "Header line must contain at least 2 integers (#vtxs and #edges).\n");

		ncols = nrows;
		nnz *= 2;

		if (fmt > 111)
			gk_errexit(SIGERR, "Cannot read this type of file format [fmt=%zu]!\n", fmt);

		sprintf(fmtstr, "%03zu", fmt%1000);
		readsizes = (fmtstr[0] == '1');
		readwgts  = (fmtstr[1] == '1');
		readvals  = (fmtstr[2] == '1');
		numbering = 1;
		ncon      = (ncon == 0 ? 1 : ncon);
	}
	else {
		readsizes = 0;
		readwgts  = 0;

		gk_getfilestats(filename, &nrows, &nnz, NULL, NULL);

		if (readvals == 1 && nnz%2 == 1)
			gk_errexit(SIGERR, "Error: The number of numbers (%zd %d) in the file %s is not even.\n",
				nnz, readvals, filename);
		if (readvals == 1)
			nnz = nnz/2;
		fpin = gk_fopen(filename, "r", "da_csr_Read: fpin");
	}

	mat = da_csr_Create();

	mat->nrows = nrows;

	rowptr = mat->rowptr = da_pmalloc(nrows+1, "da_csr_Read: rowptr");
	rowind = mat->rowind = da_imalloc(nnz, "da_csr_Read: rowind");
	if (readvals != 2)
		rowval = mat->rowval = da_vsmalloc(nnz, 1.0, "da_csr_Read: rowval");

	if (readsizes)
		mat->rsizes = da_vsmalloc(nrows, 0.0, "da_csr_Read: rsizes");

	if (readwgts)
		mat->rwgts = da_vsmalloc(nrows*ncon, 0.0, "da_csr_Read: rwgts");

	/*----------------------------------------------------------------------
	 * Read the sparse matrix file
	 *---------------------------------------------------------------------*/
	numbering = (numbering ? - 1 : 0);
	for (ncols=0, rowptr[0]=0, k=0, i=0; i<nrows; i++) {
		do {
			if (gk_getline(&line, &lnlen, fpin) == -1)
				gk_errexit(SIGERR, "Premature end of input file: file while reading row %d\n", i);
		} while (line[0] == '%');

		head = line;
		tail = NULL;

		/* Read vertex sizes */
		if (readsizes) {
#ifdef __MSC__
			mat->rsizes[i] = (float)strtod(head, &tail);
#else
			mat->rsizes[i] = strtof(head, &tail);
#endif
			if (tail == head)
				gk_errexit(SIGERR, "The line for vertex %zd does not have size information\n", i+1);
			if (mat->rsizes[i] < 0)
				errexit("The size for vertex %zd must be >= 0\n", i+1);
			head = tail;
		}

		/* Read vertex weights */
		if (readwgts) {
			for (l=0; l<ncon; l++) {
#ifdef __MSC__
				mat->rwgts[i*ncon+l] = (float)strtod(head, &tail);
#else
				mat->rwgts[i*ncon+l] = strtof(head, &tail);
#endif
				if (tail == head)
					errexit("The line for vertex %zd does not have enough weights "
							"for the %d constraints.\n", i+1, ncon);
				if (mat->rwgts[i*ncon+l] < 0)
					errexit("The weight vertex %zd and constraint %zd must be >= 0\n", i+1, l);
				head = tail;
			}
		}


		/* Read the rest of the row */
		while (1) {
			len = (int)strtol(head, &tail, 0);
			if (tail == head)
				break;
			head = tail;

			if ((rowind[k] = len + numbering) < 0)
				gk_errexit(SIGERR, "Error: Invalid column number %d at row %zd.\n", len, i);

			ncols = gk_max(rowind[k], ncols);

			if (readvals == 1) {
#ifdef __MSC__
				fval = (float)strtod(head, &tail);
#else
				fval = strtof(head, &tail);
#endif
				if (tail == head)
					gk_errexit(SIGERR, "Value could not be found for column! Row:%zd, NNZ:%zd\n", i, k);
				head = tail;

				rowval[k] = fval;
			}
			k++;
		}
		rowptr[i+1] = k;
	}

	if (format == DA_FMT_METIS) {
		ASSERT(ncols+1 == mat->nrows);
		mat->ncols = mat->nrows;
	}
	else {
		mat->ncols = ncols+1;
	}

	if (k != nnz)
		gk_errexit(SIGERR, "da_csr_Read: Something wrong with the number of nonzeros in "
				"the input file. NNZ=%zd, ActualNNZ=%zd.\n", nnz, k);

	gk_fclose(fpin);

	gk_free((void **)&line, LTERM);

	return mat;
}




void inline da_csr_WriteClean(int32_t **bptr, int32_t **bind, float **bval, char * filename, FILE *fpout){
	if (sizeof(ptr_t) != sizeof(int32_t))
		gk_free((void**) bptr, LTERM);

	if (sizeof(idx_t) != sizeof(int32_t))
		gk_free((void**) bind, LTERM);

	if (sizeof(val_t) != sizeof(float))
		gk_free((void**) bval, LTERM);

	if (fpout)
		gk_fclose(fpout);
}

/**************************************************************************/
/*! Writes the row-based structure of a matrix into a file.
    \param mat is the matrix to be written,
    \param filename is the name of the output file.
    \param format is one of: DA_FMT_CLUTO, DA_FMT_CSR,
           DA_FMT_BINROW, DA_FMT_BINCOL.
    \param writevals is either 1 or 0 indicating if the values will be
           written or not. This is only applicable when DA_FMT_CSR
           is used.
    \param numbering is either 1 or 0 indicating if the internal 0-based
           numbering will be shifted by one or not during output. This
           is only applicable when DA_FMT_CSR is used.
 */
/**************************************************************************/
void da_csr_Write(da_csr_t *mat, char *filename, char format, char writevals, char numbering)
{
	ssize_t i, j;
	size_t nnz;
	idx_t len;
	int32_t nr, nc, vId;
	int32_t edge[2], *bind = NULL, *bptr = NULL;
	float *bval = NULL;
	da_csr_t *tmp = NULL;
	FILE *fpout = NULL;

	if (!mat->rowval)
		writevals = 0;

	nr = mat->nrows;
	nc = mat->ncols;
	nnz = (format == DA_FMT_BINCOL) ? mat->colptr[mat->ncols] : mat->rowptr[mat->nrows];
	if (sizeof(ptr_t) != sizeof(int32_t)){
		if(format == DA_FMT_BINCOL){
			if(!mat->colptr)
				gk_errexit(SIGERR, "DA_FMT_BINCOL and no mat->colptr!\n");
			bptr = da_i32malloc(mat->ncols+1, "da_csr_Write: bptr");
			for(i=0; i <= mat->ncols; i++)
				bptr[i] = mat->colptr[i];
		} else {
			bptr = da_i32malloc(mat->nrows+1, "da_csr_Write: bptr");
			for(i=0; i <= mat->nrows; i++)
				bptr[i] = mat->rowptr[i];
		}
	} else
		bptr = (format == DA_FMT_BINCOL) ? (int32_t*)mat->colptr : (int32_t*)mat->rowptr;
	if (sizeof(idx_t) != sizeof(int32_t)){
		bind = da_i32malloc(nnz, "da_csr_Write: bind");
		if(format == DA_FMT_BINCOL){
			for(i=0; i < nnz; i++)
				bind[i] = mat->colind[i];
		} else {
			for(i=0; i < nnz; i++)
				bind[i] = mat->rowind[i];
		}
	} else
		bind = (format == DA_FMT_BINCOL) ? (int32_t*)mat->colind : (int32_t*)mat->rowind;
	if (writevals)
		if (sizeof(val_t) != sizeof(float)){
			bval = da_fmalloc(nnz, "da_csr_Write: bval");
			if(format == DA_FMT_BINCOL)
				for(i=0; i < nnz; i++)
					bval[i] = mat->colval[i];
			else
				for(i=0; i < nnz; i++)
					bval[i] = mat->rowval[i];
		} else
			bval = (format == DA_FMT_BINCOL) ? (float*)mat->colval : (float*)mat->rowval;

	if (format == DA_FMT_BINAP) {
		if (filename == NULL)
			gk_errexit(SIGERR, "The filename parameter cannot be NULL.\n");

		assert(mat->nrows <= INT32_MAX);
		assert(mat->ncols <= INT32_MAX);
		assert(mat->rowptr[mat->nrows] <= SIZE_MAX);

		fpout = gk_fopen(filename, "wb", "da_csr_Write: fpout");

		// first add 1 to all features
		for (i = 0; i < nnz; i++)
			bind[i]++;

		if (!bval)
			gk_errexit(SIGERR, "da_csr_Write: Format requires values but rowval not present.\n");

		// write the number of rows
		fwrite(&nr, sizeof(int32_t), 1, fpout);

		for (i = 0; i < nr; i++) {
			len = bptr[i+1] - bptr[i];
			fwrite(&len, sizeof(int32_t), 1, fpout);
			if (len > 0) {
				fwrite(bind + bptr[i], sizeof(int32_t), len, fpout);
				fwrite(bval + bptr[i], sizeof(float), len, fpout);
			}
		}

		// return rowind to its initial state
		for (i = 0; i < nnz; i++)
			bind[i]--;

		da_csr_WriteClean(&bptr, &bind, &bval, filename, fpout);
		return;
	}

	if (format == DA_FMT_BINAPB) {
		if (filename == NULL)
			gk_errexit(SIGERR, "The filename parameter cannot be NULL.\n");

		assert(mat->nrows <= INT32_MAX);
		assert(mat->ncols <= INT32_MAX);
		assert(mat->rowptr[mat->nrows] <= SIZE_MAX);

		fpout = gk_fopen(filename, "wb", "da_csr_Write: fpout");

		// first add 1 to all features
		for (i = 0; i < nnz; i++)
			bind[i]++;

		for (i = 0, vId = 1; i < nr; i++, vId++) {
			fwrite(&vId, sizeof(int32_t), 1, fpout);

			len = bptr[i+1] - bptr[i];
			fwrite(&len, sizeof(int32_t), 1, fpout);

			if (len > 0)
				fwrite(bind + bptr[i], sizeof(int32_t), len, fpout);
		}

		// return rowind to its initial state
		for (i = 0; i < nnz; i++)
			bind[i]--;

		da_csr_WriteClean(&bptr, &bind, &bval, filename, fpout);
		return;
	}

	if (format == DA_FMT_BINROW) {
		if (filename == NULL)
			gk_errexit(SIGERR, "The filename parameter cannot be NULL.\n");

		assert(mat->nrows <= INT32_MAX);
		assert(mat->ncols <= INT32_MAX);
		assert(mat->rowptr[mat->nrows] <= SIZE_MAX);

		fpout = gk_fopen(filename, "wb", "da_csr_Write: fpout");

		fwrite(&nr, sizeof(int32_t), 1, fpout);
		fwrite(&nc, sizeof(int32_t), 1, fpout);
		fwrite(bptr, sizeof(int32_t), nr+1, fpout);
		fwrite(bind, sizeof(int32_t), nnz, fpout);
		if (writevals)
			fwrite(bval, sizeof(float), nnz, fpout);

		da_csr_WriteClean(&bptr, &bind, &bval, filename, fpout);
		return;
	}

	if (format == DA_FMT_BINCOL) {
		if (filename == NULL)
			gk_errexit(SIGERR, "The filename parameter cannot be NULL.\n");

		assert(mat->nrows <= INT32_MAX);
		assert(mat->ncols <= INT32_MAX);
		assert(mat->rowptr[mat->nrows] <= SIZE_MAX);

		fpout = gk_fopen(filename, "wb", "da_csr_Write: fpout");

		fwrite(&nr, sizeof(int32_t), 1, fpout);
		fwrite(&nc, sizeof(int32_t), 1, fpout);
		fwrite(bptr, sizeof(int32_t), nc + 1, fpout);
		fwrite(bind, sizeof(int32_t), nnz, fpout);
		if (writevals)
			fwrite(bval, sizeof(val_t), nnz, fpout);

		da_csr_WriteClean(&bptr, &bind, &bval, filename, fpout);
		return;
	}

	if (format == DA_FMT_BIJV) {
		if (filename == NULL)
			gk_errexit(SIGERR, "The filename parameter cannot be NULL.\n");

		assert(mat->nrows <= INT32_MAX);
		assert(mat->ncols <= INT32_MAX);
		assert(mat->rowptr[mat->nrows] <= SIZE_MAX);

		fpout = gk_fopen(filename, "wb", "gk_csr_Write: fpout");

		fwrite(&nr, sizeof(int32_t), 1, fpout);
		fwrite(&nc, sizeof(int32_t), 1, fpout);
		fwrite(&nnz, sizeof(size_t), 1, fpout);
		fwrite(&writevals, sizeof(int32_t), 1, fpout);

		for (i = 0; i < nr; i++) {
			edge[0] = i;
			for (j = bptr[i]; j < bptr[i+1]; j++) {
				edge[1] = bind[j];
				fwrite(edge, sizeof(int32_t), 2, fpout);
				if (writevals)
					fwrite(bval+j, sizeof(float), 1, fpout);
			}
		}

		da_csr_WriteClean(&bptr, &bind, &bval, filename, fpout);
		return;
	}

	if (format == DA_FMT_METIS)
		gk_errexit(SIGERR, "METIS output format is not currently supported.\n");

	if (filename)
		fpout = gk_fopen(filename, "w", "da_csr_Write: fpout");
	else
		fpout = stdout;

	if (format == DA_FMT_IJV) {
		numbering = (numbering ? 1 : 0);
		for (i = 0; i < mat->nrows; i++) {
			for (j = bptr[i]; j < bptr[i+1]; j++) {
				if (writevals)
					fprintf(fpout, "%zd %d %.8f\n", i + numbering,
							bind[j] + numbering, bval[j]);
				else
					fprintf(fpout, "%zd %d\n", i + numbering,
							bind[j] + numbering);
			}
		}

		if (fpout)
			gk_fclose(fpout);
		return;
	}

	if (format == DA_FMT_SMAT) {
		// writevals and numbering are ignored for this format
		if (!mat->rowptr)
			gk_errexit(3,
					"da_csr_Write: Row data needed for MV_FMT_SMAT format SMAT file.");

		// Matlab requires columns to be ordered - copy matrix temporarily to order columns
		tmp = da_csr_Copy(mat);
		da_csr_LoadBases(tmp); // also sorts index
		idx_t row, col;

		for (col = row = i = 0; i < tmp->ncols; i++)
			for (j = tmp->colptr[i]; j < tmp->colptr[i + 1]; j++) {
				col = i + 1;
				row = tmp->colind[j] + 1;
				fprintf(fpout, PRNT_IDXTYPE "\t" PRNT_IDXTYPE "\t%.17g\n", row,
						col, tmp->colval[j]);
			}

		if (row != mat->nrows || col != mat->ncols)
			fprintf(fpout, PRNT_IDXTYPE "\t" PRNT_IDXTYPE "\t%.17g\n",
					mat->nrows, mat->ncols, 0.0);

		da_csr_Free(&tmp);

		if (fpout)
			gk_fclose(fpout);
		return;
	}

	assert(mat->nrows <= INT32_MAX);
	assert(mat->ncols <= INT32_MAX);
	assert(mat->rowptr[mat->nrows] <= SIZE_MAX);

	if (format == DA_FMT_CLUTO) {
		fprintf(fpout, "%d %d %zu\n",
				nr, nc, nnz);
		if (mat->rowval)
			writevals = 1;
		numbering = 1;
	}

	for (i = 0; i < nr; i++) {
		for (j = bptr[i]; j < bptr[i+1]; j++) {
			fprintf(fpout, " %d",
					bind[j] + (numbering ? 1 : 0));
			if (writevals)
				fprintf(fpout, " %f", bval[j]);
		}
		fprintf(fpout, "\n");
	}

	da_csr_WriteClean(&bptr, &bind, &bval, filename, fpout);
}


/**************************************************************************/
/*! Prints the row based representation of the matrix to screen.
    \param mat is the matrix to be printed.
*/
/**************************************************************************/
void da_csr_Print(da_csr_t *mat)
{
	da_csr_Write(mat, NULL, DA_FMT_CLUTO, 1, 1);
}


/*************************************************************************/
/*! Check if a text (non-binary) text file containing a csr is in CLUTO or CSR format.
    \param file is the matrix file to be checked.
    \return the CSR format: MV_FMT_CLUTO or MV_FMT_CSR
 */
/*************************************************************************/
char da_csr_isClutoOrCsr(char *file)
{
	size_t nrows, nnz;

	da_getfilestats(file, &nrows, &nnz, NULL, NULL);
	return (nnz%2 == 1) ? DA_FMT_CLUTO : DA_FMT_CSR;
}

/*************************************************************************/
/*! Prunes certain rows/columns of the matrix. The pruning takes place
    by analyzing the row structure of the matrix. The pruning takes place
    by removing rows/columns but it does not affect the numbering of the
    remaining rows/columns.

    \param mat the matrix to be pruned,
    \param what indicates if the rows (DA_ROW) or the columns (DA_COL)
           of the matrix will be pruned,
    \param minf is the minimum number of rows (columns) that a column (row) must
           be present in order to be kept,
    \param maxf is the maximum number of rows (columns) that a column (row) must
          be present at in order to be kept.
    \returns the pruned matrix consisting only of its row-based structure.
          The input matrix is not modified.
 */
/**************************************************************************/
da_csr_t *da_csr_Prune(da_csr_t *mat, char what, idx_t minf, idx_t maxf)
{
	ssize_t i, j, nnz;
	idx_t nrows, ncols;
	ptr_t *rowptr, *nrowptr;
	idx_t *rowind, *nrowind, *collen;
	val_t *rowval, *nrowval;
	da_csr_t *nmat;

	nmat = da_csr_Create();

	nrows = nmat->nrows = mat->nrows;
	ncols = nmat->ncols = mat->ncols;

	rowptr = mat->rowptr;
	rowind = mat->rowind;
	rowval = mat->rowval;

	nrowptr = nmat->rowptr = da_pmalloc(nrows+1, "da_csr_Prune: nrowptr");
	nrowind = nmat->rowind = da_imalloc(rowptr[nrows], "da_csr_Prune: nrowind");
	nrowval = nmat->rowval = da_vmalloc(rowptr[nrows], "da_csr_Prune: nrowval");


	switch (what) {
	case DA_COL:
		collen = da_ismalloc(ncols, 0, "da_csr_Prune: collen");

		for (i=0; i<nrows; i++) {
			for (j=rowptr[i]; j<rowptr[i+1]; j++) {
				ASSERT(rowind[j] < ncols);
				collen[rowind[j]]++;
			}
		}
		for (i=0; i<ncols; i++)
			collen[i] = (collen[i] >= minf && collen[i] <= maxf ? 1 : 0);

		nrowptr[0] = 0;
		for (nnz=0, i=0; i<nrows; i++) {
			for (j=rowptr[i]; j<rowptr[i+1]; j++) {
				if (collen[rowind[j]]) {
					nrowind[nnz] = rowind[j];
					nrowval[nnz] = rowval[j];
					nnz++;
				}
			}
			nrowptr[i+1] = nnz;
		}
		gk_free((void **)&collen, LTERM);
		break;

	case DA_ROW:
		nrowptr[0] = 0;
		for (nnz=0, i=0; i<nrows; i++) {
			if (rowptr[i+1]-rowptr[i] >= minf && rowptr[i+1]-rowptr[i] <= maxf) {
				for (j=rowptr[i]; j<rowptr[i+1]; j++, nnz++) {
					nrowind[nnz] = rowind[j];
					nrowval[nnz] = rowval[j];
				}
			}
			nrowptr[i+1] = nnz;
		}
		break;

	default:
		da_csr_Free(&nmat);
		gk_errexit(SIGERR, "Unknown pruning type of %d\n", what);
		return NULL;
	}

	return nmat;
}


/*************************************************************************/
/*! Eliminates certain entries from the rows/columns of the matrix. The
    filtering takes place by keeping only the highest weight entries whose
    sum accounts for a certain fraction of the overall weight of the
    row/column.

    \param mat the matrix to be pruned,
    \param what indicates if the rows (DA_ROW) or the columns (DA_COL)
           of the matrix will be pruned,
    \param norm indicates the norm that will be used to aggregate the weights
           and possible values are 1 or 2,
    \param fraction is the fraction of the overall norm that will be retained
           by the kept entries.
    \returns the filtered matrix consisting only of its row-based structure.
           The input matrix is not modified.
 */
/**************************************************************************/
da_csr_t *da_csr_LowFilter(da_csr_t *mat, char what, char norm, float fraction)
{
	ssize_t i, j, nnz;
	idx_t nrows, ncols, ncand, maxlen=0;
	ptr_t *rowptr, *colptr, *nrowptr;
	idx_t *rowind, *colind, *nrowind;
	val_t *rowval, *colval, *nrowval, rsum, tsum;
	da_csr_t *nmat;
	da_ivkv_t *cand;

	nmat = da_csr_Create();

	nrows = nmat->nrows = mat->nrows;
	ncols = nmat->ncols = mat->ncols;

	rowptr = mat->rowptr;
	rowind = mat->rowind;
	rowval = mat->rowval;
	colptr = mat->colptr;
	colind = mat->colind;
	colval = mat->colval;

	nrowptr = nmat->rowptr = da_pmalloc(nrows+1, "da_csr_LowFilter: nrowptr");
	nrowind = nmat->rowind = da_imalloc(rowptr[nrows], "da_csr_LowFilter: nrowind");
	nrowval = nmat->rowval = da_vmalloc(rowptr[nrows], "da_csr_LowFilter: nrowval");


	switch (what) {
	case DA_COL:
		if (mat->colptr == NULL)
			gk_errexit(SIGERR, "Cannot filter columns when column-based structure has not been created.\n");

		da_pcopy(nrows+1, rowptr, nrowptr);

		for (i=0; i<ncols; i++)
			maxlen = gk_max(maxlen, colptr[i+1]-colptr[i]);

		#pragma omp parallel private(i, j, ncand, rsum, tsum, cand)
		{
			cand = da_ivkvmalloc(maxlen, "da_csr_LowFilter: cand");

			#pragma omp for schedule(static)
			for (i=0; i<ncols; i++) {
				for (tsum=0.0, ncand=0, j=colptr[i]; j<colptr[i+1]; j++, ncand++) {
					cand[ncand].key = colind[j];
					cand[ncand].val = colval[j];
					tsum += (norm == 1 ? colval[j] : colval[j]*colval[j]);
				}
				da_ivkvsortd(ncand, cand);

				for (rsum=0.0, j=0; j<ncand && rsum<=fraction*tsum; j++) {
					rsum += (norm == 1 ? cand[j].val : cand[j].val*cand[j].val);
					nrowind[nrowptr[cand[j].key]] = i;
					nrowval[nrowptr[cand[j].key]] = cand[j].val;
					nrowptr[cand[j].key]++;
				}
			}

			gk_free((void **)&cand, LTERM);
		}

		/* compact the nrowind/nrowval */
		for (nnz=0, i=0; i<nrows; i++) {
			for (j=rowptr[i]; j<nrowptr[i]; j++, nnz++) {
				nrowind[nnz] = nrowind[j];
				nrowval[nnz] = nrowval[j];
			}
			nrowptr[i] = nnz;
		}
		SHIFTCSR(i, nrows, nrowptr);

		break;

	case DA_ROW:
		if (mat->rowptr == NULL)
			gk_errexit(SIGERR, "Cannot filter rows when row-based structure has not been created.\n");

		for (i=0; i<nrows; i++)
			maxlen = gk_max(maxlen, rowptr[i+1]-rowptr[i]);

#pragma omp parallel private(i, j, ncand, rsum, tsum, cand)
		{
			cand = da_ivkvmalloc(maxlen, "da_csr_LowFilter: cand");

#pragma omp for schedule(static)
			for (i=0; i<nrows; i++) {
				for (tsum=0.0, ncand=0, j=rowptr[i]; j<rowptr[i+1]; j++, ncand++) {
					cand[ncand].key = rowind[j];
					cand[ncand].val = rowval[j];
					tsum += (norm == 1 ? rowval[j] : rowval[j]*rowval[j]);
				}
				da_ivkvsortd(ncand, cand);

				for (rsum=0.0, j=0; j<ncand && rsum<=fraction*tsum; j++) {
					rsum += (norm == 1 ? cand[j].val : cand[j].val*cand[j].val);
					nrowind[rowptr[i]+j] = cand[j].key;
					nrowval[rowptr[i]+j] = cand[j].val;
				}
				nrowptr[i+1] = rowptr[i]+j;
			}

			gk_free((void **)&cand, LTERM);
		}

		/* compact nrowind/nrowval */
		nrowptr[0] = nnz = 0;
		for (i=0; i<nrows; i++) {
			for (j=rowptr[i]; j<nrowptr[i+1]; j++, nnz++) {
				nrowind[nnz] = nrowind[j];
				nrowval[nnz] = nrowval[j];
			}
			nrowptr[i+1] = nnz;
		}

		break;

	default:
		da_csr_Free(&nmat);
		gk_errexit(SIGERR, "Unknown pruning type of %d\n", what);
		return NULL;
	}

	return nmat;
}


/*************************************************************************/
/*! Eliminates certain entries from the rows/columns of the matrix. The
    filtering takes place by keeping only the highest weight top-K entries
    along each row/column and those entries whose weight is greater than
    a specified value.

    \param mat the matrix to be pruned,
    \param what indicates if the rows (DA_ROW) or the columns (DA_COL)
           of the matrix will be pruned,
    \param topk is the number of the highest weight entries to keep.
    \param keepval is the weight of a term above which will be kept. This
           is used to select additional terms past the first topk.
    \returns the filtered matrix consisting only of its row-based structure.
           The input matrix is not modified.
 */
/**************************************************************************/
da_csr_t *da_csr_topKPlusFilter(da_csr_t *mat, char what, idx_t topk, val_t keepval)
{
	ssize_t i, j, k, nnz;
	idx_t nrows, ncols, ncand;
	ptr_t *rowptr, *colptr, *nrowptr;
	idx_t *rowind, *colind, *nrowind;
	val_t *rowval, *colval, *nrowval;
	da_csr_t *nmat;
	da_ivkv_t *cand;

	nmat = da_csr_Create();

	nrows = nmat->nrows = mat->nrows;
	ncols = nmat->ncols = mat->ncols;

	rowptr = mat->rowptr;
	rowind = mat->rowind;
	rowval = mat->rowval;
	colptr = mat->colptr;
	colind = mat->colind;
	colval = mat->colval;

	nrowptr = nmat->rowptr = da_pmalloc(nrows+1, "da_csr_LowFilter: nrowptr");
	nrowind = nmat->rowind = da_imalloc(rowptr[nrows], "da_csr_LowFilter: nrowind");
	nrowval = nmat->rowval = da_vmalloc(rowptr[nrows], "da_csr_LowFilter: nrowval");


	switch (what) {
	case DA_COL:
		if (mat->colptr == NULL)
			gk_errexit(SIGERR, "Cannot filter columns when column-based structure has not been created.\n");

		cand = da_ivkvmalloc(nrows, "da_csr_LowFilter: cand");

		da_pcopy(nrows+1, rowptr, nrowptr);
		for (i=0; i<ncols; i++) {
			for (ncand=0, j=colptr[i]; j<colptr[i+1]; j++, ncand++) {
				cand[ncand].key = colind[j];
				cand[ncand].val = colval[j];
			}
			da_ivkvsortd(ncand, cand);

			k = gk_min(topk, ncand);
			for (j=0; j<k; j++) {
				nrowind[nrowptr[cand[j].key]] = i;
				nrowval[nrowptr[cand[j].key]] = cand[j].val;
				nrowptr[cand[j].key]++;
			}
			for (; j<ncand; j++) {
				if (cand[j].val < keepval)
					break;

				nrowind[nrowptr[cand[j].key]] = i;
				nrowval[nrowptr[cand[j].key]] = cand[j].val;
				nrowptr[cand[j].key]++;
			}
		}

		/* compact the nrowind/nrowval */
		for (nnz=0, i=0; i<nrows; i++) {
			for (j=rowptr[i]; j<nrowptr[i]; j++, nnz++) {
				nrowind[nnz] = nrowind[j];
				nrowval[nnz] = nrowval[j];
			}
			nrowptr[i] = nnz;
		}
		SHIFTCSR(i, nrows, nrowptr);

		gk_free((void **)&cand, LTERM);
		break;

	case DA_ROW:
		if (mat->rowptr == NULL)
			gk_errexit(SIGERR, "Cannot filter rows when row-based structure has not been created.\n");

		cand = da_ivkvmalloc(ncols, "da_csr_LowFilter: cand");

		nrowptr[0] = 0;
		for (nnz=0, i=0; i<nrows; i++) {
			for (ncand=0, j=rowptr[i]; j<rowptr[i+1]; j++, ncand++) {
				cand[ncand].key = rowind[j];
				cand[ncand].val = rowval[j];
			}
			da_ivkvsortd(ncand, cand);

			k = gk_min(topk, ncand);
			for (j=0; j<k; j++, nnz++) {
				nrowind[nnz] = cand[j].key;
				nrowval[nnz] = cand[j].val;
			}
			for (; j<ncand; j++, nnz++) {
				if (cand[j].val < keepval)
					break;

				nrowind[nnz] = cand[j].key;
				nrowval[nnz] = cand[j].val;
			}
			nrowptr[i+1] = nnz;
		}

		gk_free((void **)&cand, LTERM);
		break;

	default:
		da_csr_Free(&nmat);
		gk_errexit(SIGERR, "Unknown pruning type of %d\n", what);
		return NULL;
	}

	return nmat;
}


/*************************************************************************/
/*! Eliminates certain entries from the rows/columns of the matrix. The
    filtering takes place by keeping only the terms whose contribution to
    the total length of the document is greater than a user-supplied multiple
    over the average.

    This routine assumes that the vectors are normalized to be unit length.

    \param mat the matrix to be pruned,
    \param what indicates if the rows (DA_ROW) or the columns (DA_COL)
           of the matrix will be pruned,
    \param zscore is the multiplicative factor over the average contribution
           to the length of the document.
    \returns the filtered matrix consisting only of its row-based structure.
           The input matrix is not modified.
 */
/**************************************************************************/
da_csr_t *da_csr_ZScoreFilter(da_csr_t *mat, char what, float zscore)
{
	ssize_t i, j, nnz;
	idx_t nrows;
	ptr_t *rowptr, *nrowptr;
	idx_t *rowind, *nrowind;
	val_t *rowval, *nrowval, avgwgt;
	da_csr_t *nmat;

	nmat = da_csr_Create();

	nmat->nrows = mat->nrows;
	nmat->ncols = mat->ncols;

	nrows  = mat->nrows;
	rowptr = mat->rowptr;
	rowind = mat->rowind;
	rowval = mat->rowval;

	nrowptr = nmat->rowptr = da_pmalloc(nrows+1, "da_csr_ZScoreFilter: nrowptr");
	nrowind = nmat->rowind = da_imalloc(rowptr[nrows], "da_csr_ZScoreFilter: nrowind");
	nrowval = nmat->rowval = da_vmalloc(rowptr[nrows], "da_csr_ZScoreFilter: nrowval");


	switch (what) {
	case DA_COL:
		gk_errexit(SIGERR, "This has not been implemented yet.\n");
		break;

	case DA_ROW:
		if (mat->rowptr == NULL)
			gk_errexit(SIGERR, "Cannot filter rows when row-based structure has not been created.\n");

		nrowptr[0] = 0;
		for (nnz=0, i=0; i<nrows; i++) {
			avgwgt = zscore/(rowptr[i+1]-rowptr[i]);
			for (j=rowptr[i]; j<rowptr[i+1]; j++) {
				if (rowval[j] > avgwgt) {
					nrowind[nnz] = rowind[j];
					nrowval[nnz] = rowval[j];
					nnz++;
				}
			}
			nrowptr[i+1] = nnz;
		}
		break;

	default:
		da_csr_Free(&nmat);
		gk_errexit(SIGERR, "Unknown pruning type of %d\n", what);
		return NULL;
	}

	return nmat;
}


/*************************************************************************/
/*! Compacts the column-space of the matrix by removing empty columns.
    As a result of the compaction, the column numbers are renumbered.
    The compaction operation is done in place and only affects the row-based
    representation of the matrix.
    The new columns are ordered in decreasing frequency.

    \param mat the matrix whose empty columns will be removed.
 */
/**************************************************************************/
void da_csr_CompactColumns(da_csr_t *mat)
{
	ssize_t i;
	idx_t nrows, ncols, nncols;
	ptr_t *rowptr;
	idx_t *rowind, *colmap;
	da_iikv_t *clens;

	nrows  = mat->nrows;
	ncols  = mat->ncols;
	rowptr = mat->rowptr;
	rowind = mat->rowind;

	colmap = da_imalloc(ncols, "da_csr_CompactColumns: colmap");

	clens = da_iikvmalloc(ncols, "da_csr_CompactColumns: clens");
	for (i=0; i<ncols; i++) {
		clens[i].key = i;
		clens[i].val = 0;
	}

	for (i=0; i<rowptr[nrows]; i++)
		clens[rowind[i]].val++;
	da_iikvsortd(ncols, clens);

	for (nncols=0, i=0; i<ncols; i++) {
		if (clens[i].val > 0)
			colmap[clens[i].key] = nncols++;
		else
			break;
	}

	for (i=0; i<rowptr[nrows]; i++)
		rowind[i] = colmap[rowind[i]];

	mat->ncols = nncols;

	gk_free((void **)&colmap, &clens, LTERM);
}

/*************************************************************************/
/*! Compacts the row-space of the matrix by removing empty rows.
    As a result of the compaction, the row numbers are renumbered.
    The compaction operation is done in place and only affects the row-based
    representation of the matrix.

    \param mat the matrix whose empty rows will be removed.
 */
/**************************************************************************/
void da_csr_CompactRows(da_csr_t *mat)
{
	ssize_t i, j;
	idx_t nrows;
	ptr_t *rowptr;

	nrows  = mat->nrows;
	rowptr = mat->rowptr;

	for(j=0, i=0; i < nrows; i++){
		rowptr[j] = rowptr[i];
		if(rowptr[i+1] > rowptr[i])
			j++;
	}

	rowptr[j+1] = rowptr[nrows];
	mat->nrows = j;
	mat->rowptr = da_prealloc(mat->rowptr, j+1, "da_csr_CompactRows: mat->rowptr realloc");
}


/*************************************************************************/
/*! Sorts the indices in increasing order (whether matrix has values or not)
    \param mat the matrix itself,
    \param what is either DA_ROW or DA_COL indicating which set of
           indices to sort.
*/
/**************************************************************************/
void da_csr_SortIndices(da_csr_t *mat, char what)
{
	ptr_t n, nn=0;
	ptr_t *ptr;
	idx_t *ind;
	val_t *val;

	switch (what) {
	case DA_ROW:
		if (!mat->rowptr)
			gk_errexit(SIGERR, "Row-based view of the matrix does not exists.\n");

		n   = mat->nrows;
		ptr = mat->rowptr;
		ind = mat->rowind;
		val = mat->rowval;
		break;

	case DA_COL:
		if (!mat->colptr)
			gk_errexit(SIGERR, "Column-based view of the matrix does not exists.\n");

		n   = mat->ncols;
		ptr = mat->colptr;
		ind = mat->colind;
		val = mat->colval;
		break;

	default:
		gk_errexit(SIGERR, "Invalid index type of %d.\n", what);
		return;
	}

	#pragma omp parallel if (n > 100)
	{
		ssize_t i, j, k;
		da_pikv_t *cand;
		val_t *tval = NULL;

		#pragma omp single
		for (i=0; i<n; i++)
			nn = gk_max(nn, ptr[i+1]-ptr[i]);

		cand = da_pikvmalloc(nn, "da_csr_SortIndices: cand");
		if(val){
			tval = da_vmalloc(nn, "da_csr_SortIndices: tval");

			#pragma omp for schedule(static)
			for (i=0; i<n; i++) {
				for (k=0, j=ptr[i]; j<ptr[i+1]; j++) {
					if (j > ptr[i] && ind[j] < ind[j-1])
						k = 1; /* an inversion */
					cand[j-ptr[i]].key = j-ptr[i];
					cand[j-ptr[i]].val = ind[j];
					tval[j-ptr[i]]     = val[j];
				}
				if (k) {
					da_pikvsorti(ptr[i+1]-ptr[i], cand);
					for (j=ptr[i]; j<ptr[i+1]; j++) {
						ind[j] = cand[j-ptr[i]].val;
						val[j] = tval[cand[j-ptr[i]].key];
					}
				}
			}

		} else {

			#pragma omp for schedule(static)
			for (i=0; i<n; i++) {
				for (k=0, j=ptr[i]; j<ptr[i+1]; j++) {
					if (j > ptr[i] && ind[j] < ind[j-1])
						k = 1; /* an inversion */
					cand[j-ptr[i]].key = j-ptr[i];
					cand[j-ptr[i]].val = ind[j];
				}
				if (k) {
					da_pikvsorti(ptr[i+1]-ptr[i], cand);
					for (j=ptr[i]; j<ptr[i+1]; j++)
						ind[j] = cand[j-ptr[i]].val;
				}
			}

		}
		gk_free((void **)&cand, &tval, LTERM);
	}

}



/*************************************************************************/
/*! Checks that an index is sorted
    \param mat the matrix itself,
    \param what is either DA_ROW or DA_COL indicating which set of
           indices to check.
*/
/**************************************************************************/
char da_csr_CheckSortedIndex(da_csr_t *mat, char what)
{
	ptr_t n, nn=0;
	ptr_t *ptr;
	idx_t *ind;
	ssize_t i, j;
	char k=1;

	switch (what) {
	case DA_ROW:
		if (!mat->rowptr)
			gk_errexit(SIGERR, "Row-based view of the matrix does not exists.\n");

		n   = mat->nrows;
		ptr = mat->rowptr;
		ind = mat->rowind;
		break;

	case DA_COL:
		if (!mat->colptr)
			gk_errexit(SIGERR, "Column-based view of the matrix does not exists.\n");

		n   = mat->ncols;
		ptr = mat->colptr;
		ind = mat->colind;
		break;

	default:
		gk_errexit(SIGERR, "Invalid index type of %d.\n", what);
		return k;
	}

	for (k=1, i=0; i < n && k == 1; i++) {
		for (j=ptr[i]; j < ptr[i+1]; j++) {
			if (j > ptr[i] && ind[j] < ind[j-1]){
				k = 0; /* an inversion */
				break;
			}
		}
	}

	return k;
}



/*************************************************************************/
/*! Creates a row/column index from the column/row data.
    \param mat the matrix itself,
    \param what is either DA_ROW or DA_COL indicating which index
           will be created.
 */
/**************************************************************************/
void da_csr_CreateIndex(da_csr_t *mat, char what)
{
	/* 'f' stands for forward, 'r' stands for reverse */
	ssize_t i, j, k, nf, nr;
	ptr_t *fptr, *rptr;
	idx_t *find, *rind;
	val_t *fval, *rval;

	switch (what) {
	case DA_COL:
		nf   = mat->nrows;
		fptr = mat->rowptr;
		find = mat->rowind;
		fval = mat->rowval;

		if (mat->colptr) gk_free((void **)&mat->colptr, LTERM);
		if (mat->colind) gk_free((void **)&mat->colind, LTERM);
		if (mat->colval) gk_free((void **)&mat->colval, LTERM);

		nr   = mat->ncols;
		rptr = mat->colptr = da_psmalloc(nr+1, 0, "da_csr_CreateIndex: rptr");
		rind = mat->colind = da_imalloc(fptr[nf], "da_csr_CreateIndex: rind");
		rval = mat->colval = (fval ? da_vmalloc(fptr[nf], "da_csr_CreateIndex: rval") : NULL);
		break;
	case DA_ROW:
		nf   = mat->ncols;
		fptr = mat->colptr;
		find = mat->colind;
		fval = mat->colval;

		if (mat->rowptr) gk_free((void **)&mat->rowptr, LTERM);
		if (mat->rowind) gk_free((void **)&mat->rowind, LTERM);
		if (mat->rowval) gk_free((void **)&mat->rowval, LTERM);

		nr   = mat->nrows;
		rptr = mat->rowptr = da_psmalloc(nr+1, 0, "da_csr_CreateIndex: rptr");
		rind = mat->rowind = da_imalloc(fptr[nf], "da_csr_CreateIndex: rind");
		rval = mat->rowval = (fval ? da_vmalloc(fptr[nf], "da_csr_CreateIndex: rval") : NULL);
		break;
	default:
		gk_errexit(SIGERR, "Invalid index type of %d.\n", what);
		return;
	}


	for (i=0; i<nf; i++) {
		for (j=fptr[i]; j<fptr[i+1]; j++)
			rptr[find[j]]++;
	}
	MAKECSR(i, nr, rptr);

	if (rptr[nr] > 6*nr) {
		for (i=0; i<nf; i++) {
			for (j=fptr[i]; j<fptr[i+1]; j++)
				rind[rptr[find[j]]++] = i;
		}
		SHIFTCSR(i, nr, rptr);

		if (fval) {
			for (i=0; i<nf; i++) {
				for (j=fptr[i]; j<fptr[i+1]; j++)
					rval[rptr[find[j]]++] = fval[j];
			}
			SHIFTCSR(i, nr, rptr);
		}
	}
	else {
		if (fval) {
			for (i=0; i<nf; i++) {
				for (j=fptr[i]; j<fptr[i+1]; j++) {
					k = find[j];
					rind[rptr[k]]   = i;
					rval[rptr[k]++] = fval[j];
				}
			}
		}
		else {
			for (i=0; i<nf; i++) {
				for (j=fptr[i]; j<fptr[i+1]; j++)
					rind[rptr[find[j]]++] = i;
			}
		}
		SHIFTCSR(i, nr, rptr);
	}
}


/*************************************************************************/
/*! Normalizes the rows/columns of the matrix to be unit
    length.
    \param mat the matrix itself,
    \param what indicates what will be normalized and is obtained by
           specifying DA_ROW, DA_COL, DA_ROW|DA_COL.
    \param norm indicates what norm is to normalize to, 1: 1-norm, 2: 2-norm
 */
/**************************************************************************/
void da_csr_Normalize(da_csr_t *mat, char what, char norm)
{
	ssize_t i, j;
	idx_t n;
	ptr_t *ptr;
	val_t *val, sum;

	if ((what & DA_ROW) && mat->rowval) {
		n   = mat->nrows;
		ptr = mat->rowptr;
		val = mat->rowval;

		#pragma omp parallel if (ptr[n] > DA_OMPMINOPS)
		{
			#pragma omp for private(j,sum) schedule(static)
			for (i=0; i<n; i++) {
				for (sum=0.0, j=ptr[i]; j<ptr[i+1]; j++){
					if (norm == 2)
						sum += val[j]*val[j];
					else if (norm == 1)
						sum += val[j]; /* assume val[j] > 0 */
				}
				if (sum > 0) {
					if (norm == 2)
						sum=1.0/sqrt(sum);
					else if (norm == 1)
						sum=1.0/sum;
					for (j=ptr[i]; j<ptr[i+1]; j++)
						val[j] *= sum;

				}
			}
		}
	}

	if ((what & DA_COL) && mat->colval) {
		n   = mat->ncols;
		ptr = mat->colptr;
		val = mat->colval;

		#pragma omp parallel if (ptr[n] > DA_OMPMINOPS)
		{
			#pragma omp for private(j,sum) schedule(static)
			for (i=0; i<n; i++) {
				for (sum=0.0, j=ptr[i]; j<ptr[i+1]; j++)
					if (norm == 2)
						sum += val[j]*val[j];
					else if (norm == 1)
						sum += val[j];
				if (sum > 0) {
					if (norm == 2)
						sum=1.0/sqrt(sum);
					else if (norm == 1)
						sum=1.0/sum;
					for (j=ptr[i]; j<ptr[i+1]; j++)
						val[j] *= sum;
				}
			}
		}
	}
}


/*************************************************************************/
/*! Applies different row scaling methods.
    \param mat the matrix itself,
    \param type indicates the type of row scaling. Possible values are:
           DA_SCALE_MAXTF, DA_SCALE_SQRT, DA_SCALE_LOG, DA_SCALE_IDF, DA_SCALE_MAXTF2.
 */
/**************************************************************************/
void da_csr_Scale(da_csr_t *mat, char type)
{
	ssize_t i, j;
	idx_t nrows, ncols, nnzcols, bgfreq;
	ptr_t *rowptr;
	idx_t *rowind, *collen;
	val_t *rowval, *cscale, maxtf;

	nrows  = mat->nrows;
	rowptr = mat->rowptr;
	rowind = mat->rowind;
	rowval = mat->rowval;

	switch (type) {
	case DA_SCALE_MAXTF: /* TF' = .5 + .5*TF/MAX(TF) */
	#pragma omp parallel if (rowptr[nrows] > DA_OMPMINOPS)
	{
		#pragma omp for private(j, maxtf) schedule(static)
		for (i=0; i<nrows; i++) {
			maxtf = fabs(rowval[rowptr[i]]);
			for (j=rowptr[i]; j<rowptr[i+1]; j++)
				maxtf = (maxtf < fabs(rowval[j]) ? fabs(rowval[j]) : maxtf);

			for (j=rowptr[i]; j<rowptr[i+1]; j++)
				rowval[j] = .5 + .5*rowval[j]/maxtf;
		}
	}
	break;

	case DA_SCALE_MAXTF2: /* TF' = .1 + .9*TF/MAX(TF) */
	#pragma omp parallel if (rowptr[nrows] > DA_OMPMINOPS)
	{
		#pragma omp for private(j, maxtf) schedule(static)
		for (i=0; i<nrows; i++) {
			maxtf = fabs(rowval[rowptr[i]]);
			for (j=rowptr[i]; j<rowptr[i+1]; j++)
				maxtf = (maxtf < fabs(rowval[j]) ? fabs(rowval[j]) : maxtf);

			for (j=rowptr[i]; j<rowptr[i+1]; j++)
				rowval[j] = .1 + .9*rowval[j]/maxtf;
		}
	}
	break;

	case DA_SCALE_SQRT: /* TF' = .1+SQRT(TF) */
	#pragma omp parallel if (rowptr[nrows] > DA_OMPMINOPS)
	{
		#pragma omp for private(j) schedule(static)
		for (i=0; i<nrows; i++) {
			for (j=rowptr[i]; j<rowptr[i+1]; j++) {
				if (rowval[j] != 0.0)
					rowval[j] = .1+sign(rowval[j], sqrt(fabs(rowval[j])));
			}
		}
	}
	break;

	case DA_SCALE_POW25: /* TF' = .1+POW(TF,.25) */
	#pragma omp parallel if (rowptr[nrows] > DA_OMPMINOPS)
	{
		#pragma omp for private(j) schedule(static)
		for (i=0; i<nrows; i++) {
			for (j=rowptr[i]; j<rowptr[i+1]; j++) {
				if (rowval[j] != 0.0)
					rowval[j] = .1+sign(rowval[j], sqrt(sqrt(fabs(rowval[j]))));
			}
		}
	}
	break;

	case DA_SCALE_POW65: /* TF' = .1+POW(TF,.65) */
	#pragma omp parallel if (rowptr[nrows] > DA_OMPMINOPS)
	{
		#pragma omp for private(j) schedule(static)
		for (i=0; i<nrows; i++) {
			for (j=rowptr[i]; j<rowptr[i+1]; j++) {
				if (rowval[j] != 0.0)
					rowval[j] = .1+sign(rowval[j], powf(fabs(rowval[j]), .65));
			}
		}
	}
	break;

	case DA_SCALE_POW75: /* TF' = .1+POW(TF,.75) */
	#pragma omp parallel if (rowptr[nrows] > DA_OMPMINOPS)
	{
		#pragma omp for private(j) schedule(static)
		for (i=0; i<nrows; i++) {
			for (j=rowptr[i]; j<rowptr[i+1]; j++) {
				if (rowval[j] != 0.0)
					rowval[j] = .1+sign(rowval[j], powf(fabs(rowval[j]), .75));
			}
		}
	}
	break;

	case DA_SCALE_POW85: /* TF' = .1+POW(TF,.85) */
	#pragma omp parallel if (rowptr[nrows] > DA_OMPMINOPS)
	{
		#pragma omp for private(j) schedule(static)
		for (i=0; i<nrows; i++) {
			for (j=rowptr[i]; j<rowptr[i+1]; j++) {
				if (rowval[j] != 0.0)
					rowval[j] = .1+sign(rowval[j], powf(fabs(rowval[j]), .85));
			}
		}
	}
	break;

	case DA_SCALE_LOG: /* TF' = 1+log_2(TF) */
	#pragma omp parallel if (rowptr[nrows] > DA_OMPMINOPS)
	{
		double logscale = 1.0/log(2.0);
		#pragma omp for schedule(static,32)
		for (i=0; i<rowptr[nrows]; i++) {
			if (rowval[i] != 0.0)
				rowval[i] = 1+(rowval[i]>0.0 ? log(rowval[i]) : -log(-rowval[i]))*logscale;
		}
	}
	break;

	case DA_SCALE_IDF: /* TF' = TF*IDF */
		ncols  = mat->ncols;
		cscale = da_vmalloc(ncols, "da_csr_Scale: cscale");
		collen = da_ismalloc(ncols, 0, "da_csr_Scale: collen");

		for (i=0; i<nrows; i++) {
			for (j=rowptr[i]; j<rowptr[i+1]; j++)
				collen[rowind[j]]++;
		}

		#pragma omp parallel if (ncols > DA_OMPMINOPS)
		{
			#pragma omp for schedule(static)
			for (i=0; i<ncols; i++)
				cscale[i] = (collen[i] > 0 ? log(1.0*nrows/collen[i]) : 0.0);
		}

		#pragma omp parallel if (rowptr[nrows] > DA_OMPMINOPS)
		{
			#pragma omp for private(j) schedule(static)
			for (i=0; i<nrows; i++) {
				for (j=rowptr[i]; j<rowptr[i+1]; j++)
					rowval[j] *= cscale[rowind[j]];
			}
		}

		gk_free((void **)&cscale, &collen, LTERM);
		break;

	case DA_SCALE_IDF2: /* TF' = TF*IDF */
		ncols  = mat->ncols;
		cscale = da_vmalloc(ncols, "da_csr_Scale: cscale");
		collen = da_ismalloc(ncols, 0, "da_csr_Scale: collen");

		for (i=0; i<nrows; i++) {
			for (j=rowptr[i]; j<rowptr[i+1]; j++)
				collen[rowind[j]]++;
		}

		nnzcols = 0;
		#pragma omp parallel if (ncols > DA_OMPMINOPS)
		{
			#pragma omp for schedule(static) reduction(+:nnzcols)
			for (i=0; i<ncols; i++)
				nnzcols += (collen[i] > 0 ? 1 : 0);

			bgfreq = gk_max(10, (ssize_t)(.5*rowptr[nrows]/nnzcols));
			printf("nnz: " PRNT_PTRTYPE ", nnzcols: " PRNT_IDXTYPE
					", bgfreq: " PRNT_IDXTYPE "\n", rowptr[nrows], nnzcols, bgfreq);

			#pragma omp for schedule(static)
			for (i=0; i<ncols; i++)
				cscale[i] = (collen[i] > 0 ? log(1.0*(nrows+2*bgfreq)/(bgfreq+collen[i])) : 0.0);
		}

		#pragma omp parallel if (rowptr[nrows] > DA_OMPMINOPS)
		{
			#pragma omp for private(j) schedule(static)
			for (i=0; i<nrows; i++) {
				for (j=rowptr[i]; j<rowptr[i+1]; j++)
					rowval[j] *= cscale[rowind[j]];
			}
		}

		gk_free((void **)&cscale, &collen, LTERM);
		break;

	default:
		gk_errexit(SIGERR, "Unknown scaling type of %d\n", type);
		break;
	}

}


/*************************************************************************/
/*! Computes the sums of the rows/columns
    \param mat the matrix itself,
    \param what is either DA_ROW or DA_COL indicating which
           sums to compute.
 */
/**************************************************************************/
void da_csr_ComputeSums(da_csr_t *mat, char what)
{
	ssize_t i;
	idx_t n;
	ptr_t *ptr;
	val_t *val, *sums;

	switch (what) {
	case DA_ROW:
		n   = mat->nrows;
		ptr = mat->rowptr;
		val = mat->rowval;

		if (mat->rsums)
			gk_free((void **)&mat->rsums, LTERM);

		sums = mat->rsums = da_vsmalloc(n, 0, "da_csr_ComputeSums: sums");
		break;
	case DA_COL:
		n   = mat->ncols;
		ptr = mat->colptr;
		val = mat->colval;

		if (mat->csums)
			gk_free((void **)&mat->csums, LTERM);

		sums = mat->csums = da_vsmalloc(n, 0, "da_csr_ComputeSums: sums");
		break;
	default:
		gk_errexit(SIGERR, "Invalid sum type of %d.\n", what);
		return;
	}

	#pragma omp parallel for if (ptr[n] > DA_OMPMINOPS) schedule(static)
	for (i=0; i<n; i++)
		sums[i] = da_vsum(ptr[i+1]-ptr[i], val+ptr[i], 1);
}


/*************************************************************************/
/*! Computes the squared of the norms of the rows/columns
    \param mat the matrix itself,
    \param what is either DA_ROW or DA_COL indicating which
           squared norms to compute.
 */
/**************************************************************************/
void da_csr_ComputeSquaredNorms(da_csr_t *mat, char what)
{
	ssize_t i;
	idx_t n;
	ptr_t *ptr;
	val_t *val, *norms;

	switch (what) {
	case DA_ROW:
		n   = mat->nrows;
		ptr = mat->rowptr;
		val = mat->rowval;

		if (mat->rnorms) gk_free((void **)&mat->rnorms, LTERM);

		norms = mat->rnorms = da_vsmalloc(n, 0, "da_csr_ComputeSums: norms");
		break;
	case DA_COL:
		n   = mat->ncols;
		ptr = mat->colptr;
		val = mat->colval;

		if (mat->cnorms) gk_free((void **)&mat->cnorms, LTERM);

		norms = mat->cnorms = da_vsmalloc(n, 0, "da_csr_ComputeSums: norms");
		break;
	default:
		gk_errexit(SIGERR, "Invalid norm type of %d.\n", what);
		return;
	}

	#pragma omp parallel for if (ptr[n] > DA_OMPMINOPS) schedule(static)
	for (i=0; i<n; i++)
		norms[i] = da_vdot(ptr[i+1]-ptr[i], val+ptr[i], 1, val+ptr[i], 1);
}


/*************************************************************************/
/*! Computes the similarity between two rows/columns

    \param mat the matrix itself. The routine assumes that the indices
           are sorted in increasing order.
    \param i1 is the first row/column,
    \param i2 is the second row/column,
    \param what is either DA_ROW or DA_COL indicating the type of
           objects between the similarity will be computed,
    \param simtype is the type of similarity and is one of DA_SIM_COS,
           DA_SIM_JAC, DA_SIM_MIN, DA_SIM_AMIN
    \returns the similarity between the two rows/columns.
 */
/**************************************************************************/
val_t da_csr_ComputeSimilarity(da_csr_t *mat, idx_t i1, idx_t i2, char what, char simtype)
{
	idx_t nind1, nind2;
	idx_t *ind1, *ind2;
	val_t *val1, *val2, stat1, stat2, sim;

	switch (what) {
	case DA_ROW:
		if (!mat->rowptr)
			gk_errexit(SIGERR, "Row-based view of the matrix does not exists.\n");
		nind1 = mat->rowptr[i1+1]-mat->rowptr[i1];
		nind2 = mat->rowptr[i2+1]-mat->rowptr[i2];
		ind1  = mat->rowind + mat->rowptr[i1];
		ind2  = mat->rowind + mat->rowptr[i2];
		val1  = mat->rowval + mat->rowptr[i1];
		val2  = mat->rowval + mat->rowptr[i2];
		break;

	case DA_COL:
		if (!mat->colptr)
			gk_errexit(SIGERR, "Column-based view of the matrix does not exists.\n");
		nind1 = mat->colptr[i1+1]-mat->colptr[i1];
		nind2 = mat->colptr[i2+1]-mat->colptr[i2];
		ind1  = mat->colind + mat->colptr[i1];
		ind2  = mat->colind + mat->colptr[i2];
		val1  = mat->colval + mat->colptr[i1];
		val2  = mat->colval + mat->colptr[i2];
		break;

	default:
		gk_errexit(SIGERR, "Invalid index type of %d.\n", what);
		return 0.0;
	}


	switch (simtype) {
	case DA_SIM_COS:
	case DA_SIM_JAC:
		sim = stat1 = stat2 = 0.0;
		i1 = i2 = 0;
		while (i1<nind1 && i2<nind2) {
			if (i1 == nind1) {
				stat2 += val2[i2]*val2[i2];
				i2++;
			}
			else if (i2 == nind2) {
				stat1 += val1[i1]*val1[i1];
				i1++;
			}
			else if (ind1[i1] < ind2[i2]) {
				stat1 += val1[i1]*val1[i1];
				i1++;
			}
			else if (ind1[i1] > ind2[i2]) {
				stat2 += val2[i2]*val2[i2];
				i2++;
			}
			else {
				sim   += val1[i1]*val2[i2];
				stat1 += val1[i1]*val1[i1];
				stat2 += val2[i2]*val2[i2];
				i1++;
				i2++;
			}
		}
		if (simtype == DA_SIM_COS)
			sim = (stat1*stat2 > 0.0 ? sim/sqrt(stat1*stat2) : 0.0);
		else
			sim = (stat1+stat2-sim > 0.0 ? sim/(stat1+stat2-sim) : 0.0);
		break;

	case DA_SIM_MIN:
		sim = stat1 = stat2 = 0.0;
		i1 = i2 = 0;
		while (i1<nind1 && i2<nind2) {
			if (i1 == nind1) {
				stat2 += val2[i2];
				i2++;
			}
			else if (i2 == nind2) {
				stat1 += val1[i1];
				i1++;
			}
			else if (ind1[i1] < ind2[i2]) {
				stat1 += val1[i1];
				i1++;
			}
			else if (ind1[i1] > ind2[i2]) {
				stat2 += val2[i2];
				i2++;
			}
			else {
				sim   += gk_min(val1[i1],val2[i2]);
				stat1 += val1[i1];
				stat2 += val2[i2];
				i1++;
				i2++;
			}
		}
		sim = (stat1+stat2-sim > 0.0 ? sim/(stat1+stat2-sim) : 0.0);

		break;

	case DA_SIM_AMIN:
		sim = stat1 = stat2 = 0.0;
		i1 = i2 = 0;
		while (i1<nind1 && i2<nind2) {
			if (i1 == nind1) {
				stat2 += val2[i2];
				i2++;
			}
			else if (i2 == nind2) {
				stat1 += val1[i1];
				i1++;
			}
			else if (ind1[i1] < ind2[i2]) {
				stat1 += val1[i1];
				i1++;
			}
			else if (ind1[i1] > ind2[i2]) {
				stat2 += val2[i2];
				i2++;
			}
			else {
				sim   += gk_min(val1[i1],val2[i2]);
				stat1 += val1[i1];
				stat2 += val2[i2];
				i1++;
				i2++;
			}
		}
		sim = (stat1 > 0.0 ? sim/stat1 : 0.0);

		break;

	default:
		gk_errexit(SIGERR, "Unknown similarity measure %d\n", simtype);
		return -1;
	}

	return sim;

}





/*************************************************************************/
/*! Check of two matrices are equal. Note that one of the matrices may have an
 * 	additional index created to facilitate the comparison
    \param a is the matrix to be compared.
    \param b is the matrix to be compared against.
    \param p is the precision to be used in the comparison.
    \returns 1 if matrices have the same elements, 0 otherwise.
 */
/**************************************************************************/
char da_csr_Compare(da_csr_t *a, da_csr_t *b, double p)
{
	if(!(a && b)) return 0;
	if(a->ncols != b->ncols || a->nrows != b->nrows) return 0;
	ptr_t nnz = 0;

	assert((a->rowptr && b->rowptr) || (a->colptr && b->colptr));

	if(a->rowptr && b->rowptr){
		if(a->rowptr[a->nrows] != b->rowptr[b->nrows])
			return 0;
		nnz = a->rowptr[a->nrows];
		if(!da_parreq(a->nrows+1, a->rowptr, b->rowptr))
			return 0;
		if(!da_iarreq(nnz, a->rowind, b->rowind))
			return 0;
		if(a->rowval && b->rowval && !da_varreq_p(nnz, a->rowval, b->rowval, p))
			return 0;
		else if((a->rowval && !b->rowval) || (!a->rowval && b->rowval))
			return 0;
	} else if(a->colptr && b->colptr){
		if(a->colptr[a->ncols] != b->colptr[b->ncols])
			return 0;
		nnz = a->colptr[a->ncols];
		if(!da_parreq(a->ncols+1, a->colptr, b->colptr))
			return 0;
		if(!da_iarreq(nnz, a->colind, b->colind))
			return 0;
		if(a->rowval && b->rowval && !da_varreq_p(nnz, a->colval, b->colval, p))
			return 0;
		else if((a->rowval && !b->rowval) || (!a->rowval && b->rowval))
			return 0;
	} else {
		return 0;
	}
	return 1;
}



/**
 * Increase the space needed for storing nnzs in a csr matrix
 * 	\param mat The matrix to have its nnz space re-allocated
 * 	\param newNnz The new size of the row/col ind/val arrays
 */
void da_csr_Grow(da_csr_t *mat, ptr_t newNnz)
{
	if(mat->rowind){
		mat->rowind = da_irealloc(mat->rowind, newNnz, "da_csr_matrixNNzRealloc: mat->rowind");
		if(mat->rowind == NULL)
			gk_errexit(SIGERR, "da_csr_matrixNNzRealloc: Could not reallocate mat->rowind size %zu.\n", newNnz);
	}
	if(mat->rowval){
		mat->rowval = da_vrealloc(mat->rowval, newNnz, "da_csr_matrixNNzRealloc: mat->rowval");
		if(mat->rowval == NULL)
			gk_errexit(SIGERR, "da_csr_matrixNNzRealloc: Could not reallocate mat->rowval size %zu.\n", newNnz);
	}
	if(mat->colind){
		mat->colind = da_irealloc(mat->colind, newNnz, "da_csr_matrixNNzRealloc: mat->colind");
		if(mat->colind == NULL)
			gk_errexit(SIGERR, "da_csr_matrixNNzRealloc: Could not reallocate mat->colind size %zu.\n", newNnz);
	}
	if(mat->colval){
		mat->colval = da_vrealloc(mat->colval, newNnz, "da_csr_matrixNNzRealloc: mat->colval");
		if(mat->colval == NULL)
			gk_errexit(SIGERR, "da_csr_matrixNNzRealloc: Could not reallocate mat->colval size %zu.\n", newNnz);
	}
}

/**
 * Compute a partial dot product of two vectors
 */
val_t da_csr_partialDotProduct(const ptr_t *rowptr, const ptr_t *endptr,
		const idx_t *rowind, const val_t *rowval, idx_t a, idx_t b)
{
	val_t sum = 0.0;
	ptr_t aidx, bidx;
	if ( endptr[a] <= rowptr[a] || endptr[b] <= rowptr[b] ) {
		// one of the two is a null vector.
		return sum;
	}
	aidx = rowptr[a];
	bidx = rowptr[b];
	for ( ; aidx < endptr[a] && bidx < endptr[b]; )
	{
		if ( rowind[aidx] < rowind[bidx] )
			aidx++;
		else if ( rowind[bidx] < rowind[aidx] )
			bidx++;
		else {
			sum += rowval[aidx] * rowval[bidx];
			aidx++;
			bidx++;
		}
	}
	return sum;

}

/**
 * Compute the dot product of two vectors
 */
val_t da_csr_dotProduct(const ptr_t *rowptr, const idx_t *rowind,
		const val_t *rowval, idx_t a, idx_t b)
{
	val_t sum = 0.0;
	ptr_t aidx, bidx;

	if ( rowptr[a+1] <= rowptr[a] || rowptr[b+1] <= rowptr[b] ) {
		// one of the two is a null vector.
		return sum;
	}
	for (aidx = rowptr[a], bidx = rowptr[b]; aidx < rowptr[a+1] && bidx < rowptr[b+1]; )
	{
		if ( rowind[aidx] < rowind[bidx] )
			aidx++;
		else if ( rowind[bidx] < rowind[aidx] )
			bidx++;
		else {
			sum += rowval[aidx] * rowval[bidx];
			aidx++;
			bidx++;
		}
	}
	return sum;

}



/**
 * Find similar rows in the matrix
 */
idx_t da_csr_GetSimilarSmallerRows(da_csr_t *mat, idx_t rid, char noSelfSim,
		idx_t nqterms, idx_t *qind, val_t *qval, char simtype, idx_t nsim,
        float minsim, da_ivkv_t *hits, idx_t *i_marker, da_ivkv_t *i_cand)
{
	ssize_t i, ii, j, k;
	idx_t nrows, ncols, ncand;
	ptr_t *colptr;
	idx_t *colind, *marker;
	val_t *colval, *rnorms, mynorm, *rsums, mysum;
	da_ivkv_t *cand;

	if (nqterms == 0)
		return 0;

	nrows  = mat->nrows;
	ncols  = mat->ncols;
	colptr = mat->colptr;
	colind = mat->colind;
	colval = mat->colval;

	marker = (i_marker ? i_marker : da_ismalloc(nrows, -1, "da_csr_GetSimilarSmallerRows: marker"));
	cand   = (i_cand   ? i_cand   : da_ivkvmalloc(nrows, "da_csr_GetSimilarSmallerRows: cand"));

	switch (simtype) {
	case DA_SIM_COS:
		for (ncand=0, ii=0; ii<nqterms; ii++) {
			i = qind[ii];
			if (i < ncols) {
				for (j=colptr[i]; j<colptr[i+1]; j++) {
					k = colind[j];
					if(k > rid || (noSelfSim && k == rid))
						continue;
					if (marker[k] == -1) {
						cand[ncand].key = k;
						cand[ncand].val = 0;
						marker[k]       = ncand++;
					}
					cand[marker[k]].val += colval[j]*qval[ii];
				}
			}
		}
		break;

	case DA_SIM_JAC:
		for (ncand=0, ii=0; ii<nqterms; ii++) {
			i = qind[ii];
			if (i < ncols) {
				for (j=colptr[i]; j<colptr[i+1]; j++) {
					k = colind[j];
					if(k > rid || (noSelfSim && k == rid))
						continue;
					if (marker[k] == -1) {
						cand[ncand].key = k;
						cand[ncand].val = 0;
						marker[k]       = ncand++;
					}
					cand[marker[k]].val += colval[j]*qval[ii];
				}
			}
		}

		rnorms = mat->rnorms;
		mynorm = da_vdot(nqterms, qval, 1, qval, 1);

		for (i=0; i<ncand; i++)
			cand[i].val = cand[i].val/(rnorms[cand[i].key]+mynorm-cand[i].val);
		break;

	case DA_SIM_MIN:
		for (ncand=0, ii=0; ii<nqterms; ii++) {
			i = qind[ii];
			if (i < ncols) {
				for (j=colptr[i]; j<colptr[i+1]; j++) {
					k = colind[j];
					if(k > rid || (noSelfSim && k == rid))
						continue;
					if (marker[k] == -1) {
						cand[ncand].key = k;
						cand[ncand].val = 0;
						marker[k]       = ncand++;
					}
					cand[marker[k]].val += gk_min(colval[j], qval[ii]);
				}
			}
		}

		rsums = mat->rsums;
		mysum = da_vsum(nqterms, qval, 1);

		for (i=0; i<ncand; i++)
			cand[i].val = cand[i].val/(rsums[cand[i].key]+mysum-cand[i].key);
		break;

		/* Assymetric MIN  similarity */
	case DA_SIM_AMIN:
		for (ncand=0, ii=0; ii<nqterms; ii++) {
			i = qind[ii];
			if (i < ncols) {
				for (j=colptr[i]; j<colptr[i+1]; j++) {
					k = colind[j];
					if(k > rid || (noSelfSim && k == rid))
						continue;
					if (marker[k] == -1) {
						cand[ncand].key = k;
						cand[ncand].val = 0;
						marker[k]       = ncand++;
					}
					cand[marker[k]].val += gk_min(colval[j], qval[ii]);
				}
			}
		}

		mysum = da_vsum(nqterms, qval, 1);

		for (i=0; i<ncand; i++)
			cand[i].val = cand[i].val/mysum;
		break;

	default:
		gk_errexit(SIGERR, "Unknown similarity measure %d\n", simtype);
		return -1;
	}

	/* go and prune the hits that are bellow minsim */
	for (j=0, i=0; i<ncand; i++) {
		marker[cand[i].key] = -1;
		if (cand[i].val >= minsim)
			cand[j++] = cand[i];
	}
	ncand = j;

	if (nsim == -1 || nsim >= ncand) {
		nsim = ncand;
	}
	else {
		nsim = gk_min(nsim, ncand);
		da_ivkvkselectd(ncand, nsim, cand);
		da_ivkvsortd(nsim, cand);
	}

	da_ivkvcopy(nsim, cand, hits);

	if (i_marker == NULL)
		gk_free((void **)&marker, LTERM);
	if (i_cand == NULL)
		gk_free((void **)&cand, LTERM);

	return nsim;
}




/*************************************************************************/
/*! This function gets some basic statistics about the file. Same as GK's version
 *  but handles comment lines.
    \param fname is the name of the file
    \param r_nlines is the number of lines in the file. If it is NULL,
           this information is not returned.
    \param r_ntokens is the number of tokens in the file. If it is NULL,
           this information is not returned.
    \param r_max_nlntokens is the maximum number of tokens in any line
           in the file. If it is NULL this information is not returned.
    \param r_nbytes is the number of bytes in the file. If it is NULL,
           this information is not returned.
*/
/*************************************************************************/
void da_getfilestats(char *fname, size_t *r_nlines, size_t *r_ntokens,
		size_t *r_max_nlntokens, size_t *r_nbytes)
{
	size_t nlines=0, ntokens=0, max_nlntokens=0, nbytes=0, oldntokens=0, nread;
	int32_t intoken=0, lineType=0; //lineType 0-not started, 1-started ok line, 2-comment line
	char buffer[2049], *cptr;
	FILE *fpin;

	fpin = gk_fopen(fname, "r", "da_getfilestats");

	while (!feof(fpin)) {
		nread = fread(buffer, sizeof(char), 2048, fpin);
		nbytes += nread;

		buffer[nread] = '\0';  /* There is space for this one */
		for (cptr=buffer; *cptr!='\0'; cptr++) {
			if (*cptr == '%' && lineType == 0){
				lineType = 2;
			}
			else if (*cptr == '\n') {
				if(lineType != 2){
					nlines += lineType;
					ntokens += intoken;
				}
				intoken = 0;
				lineType = 0;
				if (max_nlntokens < ntokens-oldntokens)
					max_nlntokens = ntokens-oldntokens;
				oldntokens = ntokens;
			}
			else if (*cptr == ' ' || *cptr == '\t') {
				ntokens += intoken;
				intoken = 0;
			}
			else if(lineType != 2){
				intoken = 1;
				lineType = 1;
			}
		}
	}
	ntokens += intoken;
	if (max_nlntokens < ntokens-oldntokens)
		max_nlntokens = ntokens-oldntokens;

	gk_fclose(fpin);

	if (r_nlines != NULL)
		*r_nlines  = nlines;
	if (r_ntokens != NULL)
		*r_ntokens = ntokens;
	if (r_max_nlntokens != NULL)
		*r_max_nlntokens = max_nlntokens;
	if (r_nbytes != NULL)
		*r_nbytes  = nbytes;
}

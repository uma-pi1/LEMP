/*!
 \file  basetypes.h
 \brief This file contains base type definitions

 \author David C. Anastasiu
 */

#ifndef _L2AP_BASETYPES_H_
#define _L2AP_BASETYPES_H_



/* Index type should be a signed int type */
#ifdef IDXTYPE_LONG
	#define PRNT_IDXTYPE "%ld"
	typedef int64_t idx_t;
#elif defined IDXTYPE_SHORT
	#define PRNT_IDXTYPE "%hi"
	typedef short idx_t;
#else
	#define IDXTYPE_INT
	#define PRNT_IDXTYPE "%d"
	typedef int32_t idx_t;
#endif

/* Pointer type should be a signed int type, but can also be unsigned */
#ifdef PTRTYPE_SSZIE
	#define PRNT_PTRTYPE "%zu"
	typedef ssize_t ptr_t;
//#elif defined PTRTYPE_SIZE	/* TODO: Bug when using unsigned pointer types */
//	#define PRNT_PTRTYPE "%zd"
//	typedef size_t ptr_t;
//#elif defined PTRTYPE_ULONG
//	#define PRNT_PTRTYPE "%lu"
//	typedef uint64_t ptr_t;
#elif defined PTRTYPE_LONG
	#define PRNT_PTRTYPE "%ld"
	typedef int64_t ptr_t;
//#elif defined PTRTYPE_UINT
//	#define PRNT_PTRTYPE "%u"
//	typedef uint32_t ptr_t;
#else
	#define PTRTYPE_INT
	#define PRNT_PTRTYPE "%d"
	typedef int32_t ptr_t;
#endif

/* Value type should be a float or double */
#ifdef VALTYPE_DOUBLE
	typedef double val_t;
	#define PRNT_VALTYPE "%g"
#else
	#define VALTYPE_FLOAT
	#define PRNT_VALTYPE "%f"
	//typedef float val_t;
        typedef double val_t;
#endif

/* Accumulator type should be a float or double */
#ifdef ACCUMTYPE_DOUBLE
	typedef double accum_t;
	#define PRNT_ACCUMTYPE "%g"
#else
	#define ACCUMTYPE_FLOAT
	#define PRNT_ACCUMTYPE "%f"
	//typedef float accum_t;
        typedef double accum_t;
#endif

#define MAXIDX	(1<<8*sizeof(idx_t)-2)


#endif

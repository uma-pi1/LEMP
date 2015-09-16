/*!
 \file  defs.h
 \brief This file contains various constant and parameter definitions

 \author David C. Anastasiu
 */
#ifndef _L2AP_DEFS_H_
#define _L2AP_DEFS_H_

/* Versions */
#define PROGRAM_NAME        "apss"
#define VER_MAJOR           0
#define VER_MINOR           1
#define VER_SUBMINOR        0
#define VER_COMMENT         "initial public release"

/** General parameter definitions **/
//#define EXTRACOUNTS  # display extra counts at the end of program execution
//#define EXTRATIMES   # display times for sub-sections of the program execution

/** parameter definitions for apss **/
#define L2PS   // index reduction based on l-2 norm
//#define LENPS  // index reduction based on length (Lee et al. - MMJoin)
//#define RS1    // use remaining score bound to reduce candidate pool during C.G.
//#define RS2    // use the MMJoin length bound to reduce candidate pool during C.G.
//#define RS3    // RS1 + keep track of increasing cmax values in reverse/inverted index
#define RS4    // use the l2-norm bound to reduce candidate pool during C.G.
#define L2CG   // check accum + ||x_p|| * ||y_p|| during candidate generation
#define L2CV   // check accum + ||x_p|| * ||y_p|| during candidate verification
//#define LENCG  // check accum + 1/2||x_p||^2 + 1/2||y_p||^2 during candidate generation
//#define LENCV  // check accum + 1/2||x_p||^2 + 1/2||y_p||^2 during candidate verification
//#define SZ3    // reduce index based on minsize bound
//#define SZ1    // reduce index based on Bayardo's minsize bound
#define PSCV   // store pscore - last term and use during candidate pruning
//#define DP1    // Awekar's MMJoin bound: \dotp(\vec x,\vec y) &\le A[y] + \min(\rmax_x \times \rsum_y', \rmax_y' \times \rsum_x )
//#define DP2    // DP1 + prefix max weights: \dotp(\vec x,\vec y) &\le A[y] + \min(\rmax_{x_p}' \times \rsum_{y_p}', \rmax_y' \times \rsum_x)
//#define DP3    // DP1 + prefix sums: \dotp(\vec x,\vec y) &\le A[y] + \min(\rmax_x \times \rsum_{y_p}', \rmax_y' \times \rsum_{x_p}')
//#define DP4    // DP1 + prefix max weights and sums: \dotp(\vec x,\vec y) &\le A[y] + \min(\rmax_{x_p}' \times \rsum_{y_p}', \rmax_y' \times \rsum_{x_p}')
#define DP5    // Bayardo's AP bound: \dotp(\vec x,\vec y) &\le A[y] + \min(|\vec x|,|\vec y'|) \times \rmax_x \times \rmax_y'
#define DP6    // DP5 + prefix max weights: \dotp(\vec x,\vec y) &\le A[y] + \min(|\vec x|,|\vec y_p'|) \times \rmax_{x_p}' \times \rmax_{y_p}'
//#define DP7    // DP5 + prefix sizes: \dotp(\vec x,\vec y) &\le A[y] + \min(|\vec x_p'|,|\vec y_p'|) \times \rmax_x \times \rmax_{y_p}'
//#define DP8    // DP5 + prefix max weights and sizes: \dotp(\vec x,\vec y) &\le A[y] + \min(|\vec x_p'|,|\vec y_p'|) \times \rmax_{x_p}' \times \rmax_{y_p}'

/* Bayes LSH params */
#define PI	3.14159265
#define NBITS_SKETCH         32             /* must be 32 */
#define F2I2F                65535          /* conversion factor for floatToInt and intToFloat */
#define LARGE_PRIME          2147483647
#define PERC_NUM_PRUNED      0.3            /* percent of numPruned of earlier stage to exhibit to go to next stage */
#define LENBNDEPS            1.5e-4         /* epsilon for length bound 2 */

/* initial number of neighbors per row to allocate in neighborhood array */
#ifndef NINITNEIGHBORS
	#define NINITNEIGHBORS 10
#endif
/* size of neighbors cache: NSIMSBUFFER * (2 idx + 1 val) */
#ifndef NSIMSBUFFER
	#define NSIMSBUFFER 16384
#endif

/* housekeeping for pruning params */
#if defined(L2CG) || defined(L2CV) || defined(RS4)
#define IDXL2
#endif
#if defined(LENCG) || defined(LENCV) || defined(RS2)
#define IDXLEN
#endif
#if defined(SZ1) && !defined(SZ3)
#define SZ3
#endif
#if defined(RS3) && !defined(RS1)
#define RS1
#endif

/* Command-line option codes */
#define CMD_SIMTHRESHOLD        21
#define CMD_EPSILON             22
#define CMD_NSKETCHES           23
#define CMD_NIM                 25
#define CMD_FMT_WRITE           31
#define CMD_FMT_WRITE_NUM       32
#define CMD_WRITE_VALS          33
#define CMD_WRITE_NUMBERING     34
#define CMD_FMT_READ            35
#define CMD_FMT_READ_NUM        36
#define CMD_READ_VALS           37
#define CMD_READ_NUMBERING      38
#define CMD_SEED                40
#define CMD_SCALE               51
#define CMD_NORM                52
#define CMD_COMPACT_COLS        53
#define CMD_COMPACT_ROWS        54
#define CMD_PR                  55
#define CMD_PRMINLEN            56
#define CMD_PRMAXLEN            57
#define CMD_PC                  58
#define CMD_PCMINLEN            59
#define CMD_PCMAXLEN            60
#define CMD_IGNORE_SIZE         65
#define CMD_FLDELTA             90
#define CMD_VERBOSITY           105
#define CMD_VERSION             109
#define CMD_HELP                110

#define USE_GKRAND              1

/* Execution modes */
#define MODE_TESTEQUAL          99  /* Test whether two matrices contain the same values */
#define MODE_IO                 98  /* Transform a matrix from some format into another */
#define MODE_INFO               97  /* Find information about a matrix */
#define MODE_IDXJOIN            1   /* IdxJoin */
#define MODE_AP                 2   /* AllPairs */
#define MODE_AP2                3   /* AllPairs new ordering */
#define MODE_L2AP               4   /* L2-norm AllPairs */
#define MODE_L2APBLSH           5   /* L2-norm AllPairs + BayesLSH-Lite pruning */
#define MODE_MMJOIN             6   /* MMJoin algo */

/* Similarity type */
#define DA_SIM_COS        1     /* Cosine similarity */
#define DA_SIM_JAC        2     /* Jaccard index/coefficient */
#define DA_SIM_MIN        3     /* Minimum distance/similarity */
#define DA_SIM_AMIN       4     /* Assymetric MIN similarity */
#define DA_SIM_EUC        55    /* Euclidean distance/similarity */
#define DA_SIM_DICE       56    /* Sorensen-Dice coefficient */
#define DA_SIM_OVER       57    /* Overlap distance */
#define DA_SIM_MAN        58    /* Manhattan distance */
#define DA_SIM_TAN        59    /* Tanimoto coefficient */

/* CSR structure components */
#define DA_ROW                  1   /* row-based structure */
#define DA_COL                  2   /* col-based structure */
#define DA_ROWCOL               3   /* both row and col-based */

/* scaling types */
#define DA_SCALE_MAXTF    1    /* TF' = .5 + .5*TF/MAX(TF) */
#define DA_SCALE_MAXTF2   10   /* TF' = .1 + .9*TF/MAX(TF) */
#define DA_SCALE_SQRT     2    /* TF' = .1+SQRT(TF) */
#define DA_SCALE_POW25    3    /* TF' = .1+POW(TF,.25) */
#define DA_SCALE_POW65    4    /* TF' = .1+POW(TF,.65) */
#define DA_SCALE_POW75    5    /* TF' = .1+POW(TF,.75) */
#define DA_SCALE_POW85    6    /* TF' = .1+POW(TF,.85) */
#define DA_SCALE_LOG      7    /* TF' = 1+log_2(TF) */
#define DA_SCALE_IDF      8    /* TF' = TF*IDF */
#define DA_SCALE_IDF2     9    /* TF' = TF*IDF ?? */

/* CSR input formats */
#define DA_FMT_CSR          2
#define DA_FMT_METIS        3
#define DA_FMT_CLUTO        1
#define DA_FMT_BINROW       4
#define DA_FMT_BINCOL       5
#define DA_FMT_IJV          6
#define DA_FMT_BIJV         7
#define DA_FMT_SMAT         50  /* matlab sparse matrix */
#define DA_FMT_BINAP        51  /* satubin - binary format used by allPairs */
#define DA_FMT_BINAPB       52  /* satubin - binary format used by allPairs and PPJoin programs for binary (non-real) input */

extern const gk_StringMap_t mode_options[];
extern const gk_StringMap_t sim_options[];
extern const gk_StringMap_t fmt_options[];
extern const gk_StringMap_t scale_options[];

#endif
  
 

/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Enable assertion. */
#define EXAFMM_ASSERT 1

/* Count number of M2L and P2P kernel calls. */
/* #undef EXAFMM_COUNT_KERNEL */

/* Count interaction list per cell. */
/* #undef EXAFMM_COUNT_LIST */

/* Define to enable extra debugging options. */
/* #undef EXAFMM_DEBUG */

/* Define to enable AVX optimizations. */
/* #undef EXAFMM_HAVE_AVX */

/* Define to enable MIC optimizations. */
/* #undef EXAFMM_HAVE_MIC */

/* Define if you have the MPI library. */
#define EXAFMM_HAVE_MPI 1

/* Define to enable ARM NEON optimizations. */
/* #undef EXAFMM_HAVE_NEON */

/* Define if OpenMP is enabled */
#define EXAFMM_HAVE_OPENMP 1

/* Define to enable SSE/SSE3 optimizations. */
#define EXAFMM_HAVE_SSE3 1

/* Define to compile in single precision. */
//#define EXAFMM_SINGLE 1

/* Enable DAG recorder. */
/* #undef EXAFMM_USE_DAG */

/* Enable Kahan summation. */
/* #undef EXAFMM_USE_KAHAN */

/* Enable PAPI performance counter. */
/* #undef EXAFMM_USE_PAPI */

/* Disable SIMD optimizations. */
/* #undef EXAFMM_USE_SIMD */

/* Enable thread tracing. */
/* #undef EXAFMM_USE_TRACE */

/* Enable weight for partitioning. */
/* #undef EXAFMM_USE_WEIGHT */

/* Use Intel Cilk */
/* #undef EXAFMM_WITH_CILK */

/* Use MassiveThreads */
/* #undef EXAFMM_WITH_MTHREAD */

/* Use QThreads */
/* #undef EXAFMM_WITH_QTHREAD */

/* Use Strumpack */
/* #undef EXAFMM_WITH_STRUMPACK */

/* Use Intel TBB */
//#define EXAFMM_WITH_TBB 1

/* Name of package */
#define PACKAGE "exafmm"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME "exaFMM"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "exaFMM 1.0.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "exafmm"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0.1"

/* Version number of package */
#define VERSION "1.0.1"

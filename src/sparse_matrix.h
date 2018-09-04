/*
   Basically, this implements a hash and some helper functions for using
   the hash as an n-dimensional sparse matrix. Reading and writing should
   be fast, deleting relatively fast. Getting the keys and values in order
   goes through quicksort. No arithmetic functions are implemented.
*/
#ifndef __SPARSE_MATRIX_HXX
#define __SPARSE_MATRIX_HXX

#ifdef __cplusplus
extern "C" {
#endif

#ifndef HUGE_MATRICES
/* If we don't need too big matrices, use int as array indices (faster,
   more memory efficient). Otherwise use ssize_t. Must be signed type.
*/
typedef int indextype;
#define PRIME_NUMBERS                                                          \
  {                                                                            \
    5, 1009, 10007, 49999, 99991, 500009, 1064281, 10000019, 25000033,         \
        50000021, 100000007, 250000013, 500000003, 1000000007                  \
  }
#define NUM_OF_PRIMES 14
#else
typedef ssize_t indextype;
* /
#define PRIME_NUMBERS                                                          \
  {                                                                            \
    5, 1009, 10007, 49999, 99991, 500009, 1064281, 10000019, 25000033,         \
        50000021, 100000007, 250000013, 500000003, 1000000007, 2500000001      \
  }
#define NUM_OF_PRIMES 15
#endif

#include "matrix_common.h"
#include <stdlib.h>
//#include "unistd.h"
//#define no_inline_funcs // FIXME: This is required for debug compile, figure
//out why?

typedef unsigned char byte;

struct matrix {
  /* General information */
  size_t key_size;
  indextype hashsize;
  indextype datasize;
  void *default_value;
  size_t size_of_entry;
  int allow_shrinking;

  /* The real data storage */
  indextype *hash; /* Mapping from hash value to table indices */
  byte *keys;      /* Table containing all of the indices */
  byte *data;      /* All of the stored data              */
  indextype *next; /* The links for filled hashes         */
  indextype *prev; /* This is for making deletions fast */
                   /* positive values: ref from hash, negative-1
                      ref from next-array. */
  indextype num_entries;

  int type_of_entry; /* not necessary */
};

struct matrix *CreateMatrixI(const int d, const indextype hashs,
                             const int def_value);
struct matrix *CreateMatrixF(const int d, const indextype hashs,
                             const float def_value);
struct matrix *CreateMatrixD(const int d, const indextype hashs,
                             const double def_value);
void DeleteMatrix(struct matrix *);
void ClearMatrix(struct matrix *);

void SetValueI(struct matrix *m, const int *indices, const int value);
void SetValueF(struct matrix *m, const int *indices, const float value);
void SetValueD(struct matrix *m, const int *indices, const double value);

int GetValueI(struct matrix *m, const int *indices);
float GetValueF(struct matrix *m, const int *indices);
double GetValueD(struct matrix *m, const int *indices);

int IncrementI(struct matrix *m, const int *indices, const int value);
float IncrementF(struct matrix *m, const int *indices, const float value);

/* If you need to use you own data structures, you can use these */
void *OrderedStepThrough(struct matrix *m, byte *key, void *value);
void *StepThrough(struct matrix *mat, byte *key, void *data);
/*void Create_Intkeys(struct matrix *m, const int dims, const int *dsizes);*/
struct matrix *CreateMatrix(const int keysize, const indextype hashs,
                            const void *def_value, const int toe,
                            const int soe);

/* private funcs */
indextype FindEntry(struct matrix *m, const byte *indices, const int create);
void RemoveEntryIdx(struct matrix *m, const indextype idx);
void DeleteCurrentST();
void RemoveDefaultValues(struct matrix *m);

#include <assert.h>
#include <stdlib.h>
#include <string.h>
/* Wrappers for the raw functions */
#if !defined no_inline_funcs
// Inline next functions in c++ for extra speed
// Does ansi c define inline ?
#define sm_inline_void static inline void
#define sm_inline_int static inline int
#define sm_inline_float static inline float
#define sm_inline_double static inline double
#define sm_inline_indextype static inline indextype
#include "sparse_matrix_pinline.h"

#else
void SetRawValue(struct matrix *m, const byte *indices, const void *value);
void GetRawValue(struct matrix *m, const byte *indices, void *value);
indextype CalculateHashIndex(const struct matrix *m, const byte *indices);
void *OrderedStepThroughI(struct matrix *m, int *indices, int *value);
void *OrderedStepThroughF(struct matrix *m, int *indices, float *value);
void *OrderedStepThroughD(struct matrix *m, int *indices, double *value);
void *StepThroughI(struct matrix *m, int *indices, int *value);
void *StepThroughF(struct matrix *m, int *indices, float *value);
void *StepThroughD(struct matrix *m, int *indices, double *value);
#endif
/*********************************************************************/
/* These functions are for internal use only. They are included in   */
/* this header only for those, who really need to poke around the    */
/* internals of this library.                                        */
/*********************************************************************/
indextype NewEntry(struct matrix *m);
void CheckConsistency(struct matrix *m);          /* for debugging */
void showvalues(struct matrix *m, indextype idx); /* for debugging */
#ifdef __cplusplus
}
#endif

#endif // __SPARSE_MATRIX_HXX

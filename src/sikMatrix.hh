// C++ interface to the sparce matrix library
#ifndef SMATRIX_HH
#define SMATRIX_HH

#include "sparse_matrix.h"
#include <vector>

template <typename K, typename V> class sikMatrix {
public:
  sikMatrix(const int dims, const indextype hashsize, const V def_value);
  ~sikMatrix();

  void setvalue(const K *indices, const V value);
  void setvalue(const K *indices, const V *value);
  V getvalue(const K *indices);
  void getvalue(const K *indices, V *value);
  V increment(const K *indices, const V &value);
  V increment_wo_del(const K *indices, const V &value);
  V increment(const K *indices, const V *value);
  inline void *stepthrough(const bool init, K *indices, V *value);
  inline void delete_current_st() { DeleteCurrentST(); }
  inline void *ordered_stepthrough(const bool init, K *indices, V *value);
  inline void *stept(const bool init, K *indices, V *value);

  inline K *Idx2Keyp(indextype idx);
  inline V *Idx2Valp(indextype idx);
  inline indextype num_entries() { return m->num_entries; }

  inline void clear() { ClearMatrix(m); }

  /* These should probably be private */
  struct matrix *m;
  int dims;
  bool stept_sortsearch;

  /* This is for debugging */
  inline void printmatrix();
};

template <typename K, typename V>
sikMatrix<K, V>::sikMatrix(const int dims_in, const indextype hashsize,
                           const V def_value) {
  dims = dims_in;

  // fprintf(stderr,"CreateMatrix with %d %d %d\n", sizeof(K), dims, (int)
  // hashsize);
  m = CreateMatrix(sizeof(K) * dims, hashsize, (void *)&def_value, UNKNOWN,
                   sizeof(V));
  stept_sortsearch = false;
}

template <typename K, typename V> sikMatrix<K, V>::~sikMatrix() {
  DeleteMatrix(m);
}

template <typename K, typename V>
void sikMatrix<K, V>::setvalue(const K *indices, const V value) {
  SetRawValue(m, (byte *)indices, (void *)&value);
}

template <typename K, typename V>
void sikMatrix<K, V>::setvalue(const K *indices, const V *value) {
  SetRawValue(m, (byte *)indices, (void *)value);
}

template <typename K, typename V>
V sikMatrix<K, V>::getvalue(const K *indices) {
  V value;
  GetRawValue(m, (byte *)indices, (void *)&value);
  return (value);
}

template <typename K, typename V>
void sikMatrix<K, V>::getvalue(const K *indices, V *value) {
  GetRawValue(m, (byte *)indices, (void *)value);
}

template <typename K, typename V>
V sikMatrix<K, V>::increment(const K *indices, const V &value) {
  indextype idx = FindEntry(m, (byte *)indices, 1);
  ;
  V *v = (V *)&(m->data[idx * m->size_of_entry]);
  *v += value;
  if (memcmp(v, m->default_value, m->size_of_entry))
    return (*v);
  RemoveEntryIdx(m, idx);
  return (*((V *)(m->default_value)));
}

template <typename K, typename V>
V sikMatrix<K, V>::increment_wo_del(const K *indices, const V &value) {
  indextype idx = FindEntry(m, (byte *)indices, 1);
  ;
  V *v = (V *)&(m->data[idx * m->size_of_entry]);
  *v += value;
  return (*v);
}

template <typename K, typename V>
V sikMatrix<K, V>::increment(const K *indices, const V *value) {
  indextype idx = FindEntry(m, (byte *)indices, 1);
  ;
  V *v = (V *)&(m->data[idx * m->size_of_entry]);
  *v += *value;
  if (memcmp(v, m->default_value, m->size_of_entry))
    return (*v);
  RemoveEntryIdx(m, idx);
  return (*((V *)(m->default_value)));
}

template <typename K, typename V>
void *sikMatrix<K, V>::stepthrough(const bool init, K *indices, V *value) {
  if (!init)
    return (StepThrough(NULL, (byte *)indices, (void *)value));
  return (StepThrough(m, (byte *)indices, (void *)value));
}

template <typename K, typename V>
void *sikMatrix<K, V>::ordered_stepthrough(const bool init, K *indices,
                                           V *value) {
  if (!init)
    return (OrderedStepThrough(NULL, (byte *)indices, (void *)value));
  return (OrderedStepThrough(m, (byte *)indices, (void *)value));
}

template <typename K, typename V>
void *sikMatrix<K, V>::stept(const bool init, K *indices, V *value) {
  if (stept_sortsearch)
    return (ordered_stepthrough(init, indices, value));
  return (stepthrough(init, indices, value));
}

template <typename K, typename V> K *sikMatrix<K, V>::Idx2Keyp(indextype idx) {
  return ((K *)&(m->keys[idx * m->key_size]));
}

template <typename K, typename V> V *sikMatrix<K, V>::Idx2Valp(indextype idx) {
  return ((V *)&(m->data[idx * m->size_of_entry]));
}

// The rest is for debugging
#include "def.hh"
template <typename K, typename V> void sikMatrix<K, V>::printmatrix() {
  for (indextype i = 0; i < num_entries(); i++) {
    print_indices(stderr, Idx2Keyp(i), dims);
    fprintf(stderr, "=%d\n", (int)*Idx2Valp(i));
  }
}
#endif

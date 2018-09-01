#ifndef SPARSE_MATRIX_PINLINE_HH
#define SPARSE_MATRIX_PINLINE_HH
sm_inline_void SetRawValue(struct matrix *m, const byte *key,
                           const void *value) {
  indextype idx = FindEntry(m, key, 1);
  if (memcmp(value, m->default_value, m->size_of_entry)) {
    memcpy(&(m->data[idx * m->size_of_entry]), value, m->size_of_entry);
  } else
    RemoveEntryIdx(m, idx);
}

sm_inline_void GetRawValue(struct matrix *m, const byte *key, void *value) {
  indextype idx = FindEntry(m, key, 0);
  if (idx < 0) {
    memcpy(value, m->default_value, m->size_of_entry);
  } else
    memcpy(value, &(m->data[idx * m->size_of_entry]), m->size_of_entry);
}

#if 0
sm_inline_indextype CalculateHashIndex(const struct matrix *m,
                                       const int *indices) {
  int i,multiply=1;
  unsigned long pureindex=0;
  /* Using unsigned int here guaranteers, that we alway get a positive index,
     even if the number would overflow. Not very good, since it introduces
     another modulo to the executable, but good enough for me
     (for the moment) */
  for (i=0;i<m->dims;i++) {
    pureindex+=indices[i]*multiply;
    multiply*=m->dimsize[i];
  }
  return(pureindex%m->hashsize);
}
#else
/* Let's try rotating hash, see http://burtleburtle.net/bob/hash/doobs.html */
sm_inline_indextype CalculateHashIndex(const struct matrix *m,
                                       const byte *key) {
  unsigned long hash, i;
  for (hash = m->key_size, i = 0; i < m->key_size; ++i) {
    hash = (hash << 4) ^ (hash >> 28) ^ key[i];
  }
  return (hash % m->hashsize);
}
#endif

sm_inline_void *OrderedStepThroughI(struct matrix *m, int *indices,
                                    int *value) {
  return (OrderedStepThrough(m, (byte *)indices, (void *)value));
}

sm_inline_void *OrderedStepThroughF(struct matrix *m, int *indices,
                                    float *value) {
  return (OrderedStepThrough(m, (byte *)indices, (void *)value));
}

sm_inline_void *OrderedStepThroughD(struct matrix *m, int *indices,
                                    double *value) {
  return (OrderedStepThrough(m, (byte *)indices, (void *)value));
}

sm_inline_void *StepThroughI(struct matrix *m, int *indices, int *value) {
  return (StepThrough(m, (byte *)indices, (void *)value));
}

sm_inline_void *StepThroughF(struct matrix *m, int *indices, float *value) {
  return (StepThrough(m, (byte *)indices, (void *)value));
}

sm_inline_void *StepThroughD(struct matrix *m, int *indices, double *value) {
  return (StepThrough(m, (byte *)indices, (void *)value));
}

sm_inline_int *Key2Intp(struct matrix *m, indextype idx) {
  return ((int *)&(m->keys[idx * m->key_size]));
}
#endif

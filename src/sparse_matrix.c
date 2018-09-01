#include "sparse_matrix.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void check_malloc(void *p, char *erstring) {
  if (p != NULL)
    return;
  fprintf(stderr, "Memory allocation failed: %s\nBailing out\n", erstring);
  exit(-1);
}

#if 0
void Create_Intkeys(struct matrix *m, const int dims,
                                    const int *dsizes) {
  /* To make the code take any kind of keys, we should get rid
     of dependencies to these variables. Currently, there has been no
     need for anyhting else than integer keys. Only thing we still depend
     on is m->dims */
  int i;
  m->dims=dims;
  m->dimsize=(int *) malloc(m->key_size);
  check_malloc(m->dimsize,"dims");
  for (i=0;i<m->dims;i++)
    m->dimsize[i]=dsizes[i];
}
#endif

struct matrix *CreateMatrix(const int keysize, const indextype hashs,
                            const void *def_value, const int toe,
                            const int soe) {
  struct matrix *m;
  indextype pr[NUM_OF_PRIMES] = PRIME_NUMBERS;
  indextype i;

  m = (struct matrix *)malloc(sizeof(struct matrix));
  check_malloc(m, "structure");
  m->key_size = keysize;
  m->type_of_entry = toe;
  m->allow_shrinking = 1;

  switch (m->type_of_entry) {
  case INT:
    m->size_of_entry = sizeof(int);
    break;
  case FLOAT:
    m->size_of_entry = sizeof(float);
    break;
  case DOUBLE:
    m->size_of_entry = sizeof(double);
    break;
  case UNKNOWN:
    m->size_of_entry = soe;
  }
  if (soe != m->size_of_entry) {
    fprintf(stderr,
            "Conflicting memory requirement, type of entry %d, size of entry "
            "%d\nExit.\n",
            m->type_of_entry, (int)m->size_of_entry);
    exit(-1);
  }
  m->default_value = (void *)malloc(m->size_of_entry);
  check_malloc(m->default_value, "default value");
  memcpy(m->default_value, def_value, m->size_of_entry);

  m->num_entries = 0;
  m->datasize = hashs;
  // fprintf(stderr,"mallocing %d*%d=%d\n", (int) m->key_size, hashs, (int)
  // (m->key_size *hashs));
  m->keys = (byte *)malloc(m->key_size * hashs);
  check_malloc(m->keys, "keys");
  m->data = (byte *)malloc(m->size_of_entry * hashs);
  check_malloc(m->data, "data");
  m->next = (indextype *)malloc(sizeof(indextype) * hashs);
  check_malloc(m->next, "next pointers");
  m->prev = (indextype *)malloc(sizeof(indextype) * hashs);
  check_malloc(m->prev, "prev pointers");

  /* Find a good hash size */
  i = -1;
  while (pr[++i] < hashs) {
    if (i == NUM_OF_PRIMES - 1) {
      i = NUM_OF_PRIMES - 1;
      fprintf(stderr,
              "Trying to create a too big hash.\nLimited the hash size to %d ! "
              "To be able to create bigger hashes, add a sufficiently big "
              "prime to sparse_matrix.h.\n",
              (int)pr[i]);
      // exit(-1);
    }
  }
  m->hashsize = pr[i];
  m->hash = (indextype *)malloc(m->hashsize * sizeof(indextype));
  check_malloc(m->hash, "hash table");
  for (i = 0; i < m->hashsize; i++) {
    m->hash[i] = -1;
  }
  return (m);
}

indextype NewEntry(struct matrix *m) {
  /* WARNING: This function modifies the m->keys, m->data and
     m->next. Please do not use the return value of this function directly
     in an assigment. m->next[a]=NewEntry(m) can fail and cause lots of
     trouble.
  */

  if (m->num_entries >= m->datasize) {
    m->datasize *= 2;
    m->keys = (byte *)realloc(m->keys, m->datasize * m->key_size);
    check_malloc(m->keys, "realloc keys");
    m->data = (byte *)realloc(m->data, m->datasize * m->size_of_entry);
    check_malloc(m->data, "realloc data");
    m->next = (indextype *)realloc(m->next, m->datasize * sizeof(indextype));
    check_malloc(m->next, "realloc next pointers");
    m->prev = (indextype *)realloc(m->prev, m->datasize * sizeof(indextype));
    check_malloc(m->prev, "realloc prev pointers");
  } else if (m->allow_shrinking && (m->num_entries < m->datasize / 3 &&
                                    m->num_entries > m->hashsize)) {
    m->datasize /= 2;
    m->keys = (byte *)realloc(m->keys, m->datasize * m->key_size);
    check_malloc(m->keys, "shrink");
    m->data = (byte *)realloc(m->data, m->datasize * m->size_of_entry);
    check_malloc(m->data, "shrink");
    m->next = (indextype *)realloc(m->next, m->datasize * sizeof(indextype));
    check_malloc(m->next, "shrink");
    m->prev = (indextype *)realloc(m->prev, m->datasize * sizeof(indextype));
    check_malloc(m->prev, "shrink");
  }
  return (m->num_entries++);
}

void RemoveEntryIdx(struct matrix *m, const indextype idx) {
  /* Fix the next and prev fields for the next and prev node
     of the entry to be removed*/
  assert(idx < m->num_entries);
  // fprintf(stderr,"rem%d/%d \n",idx,m->num_entries);

  if (m->prev[idx] >= 0)
    m->hash[m->prev[idx]] = m->next[idx];
  else
    m->next[-m->prev[idx] - 1] = m->next[idx];

  if (m->next[idx] >= 0) {
    m->prev[m->next[idx]] = m->prev[idx];
  }

  m->num_entries--;
  /* Special case, this is already the last */
  if (idx == m->num_entries)
    return;

  /* Fix the next fields for the prev of entry that will be moved */
  if (m->prev[m->num_entries] >= 0)
    m->hash[m->prev[m->num_entries]] = idx;
  else
    m->next[-m->prev[m->num_entries] - 1] = idx;
  if (m->next[m->num_entries] >= 0)
    m->prev[m->next[m->num_entries]] = -idx - 1;

  /* copy the cur fields */
  m->next[idx] = m->next[m->num_entries];
  m->prev[idx] = m->prev[m->num_entries];

  /* copy the data and keys */
  memcpy(&(m->keys[idx * m->key_size]),
         &(m->keys[m->num_entries * m->key_size]), m->key_size);
  memcpy(&(m->data[idx * m->size_of_entry]),
         &(m->data[m->num_entries * m->size_of_entry]), m->size_of_entry);
  return;
}

indextype FindEntry(struct matrix *m, const byte *key, const int create) {
  indextype hashidx, *storeidx, idx;

  hashidx = CalculateHashIndex(m, key);
  storeidx = &(m->hash[hashidx]);

  if (*storeidx < 0) { /* There is no node here */
    if (create <= 0)
      return (-1); /* Delete or read, nothing to do */
    *storeidx = NewEntry(m);
    m->next[*storeidx] = -1;
    m->prev[*storeidx] = hashidx;
    idx = *storeidx;
  } else {
    const int cr =
        memcmp(&(m->keys[*storeidx * m->key_size]), key, m->key_size);
    if (cr > 0) {
      /* The node is to be created as the first of the list */
      if (create <= 0)
        return (-1);
      idx = *storeidx;
      *storeidx = NewEntry(m);
      m->next[*storeidx] = idx;
      m->prev[*storeidx] = hashidx;
      m->prev[idx] = -*storeidx - 1;
      idx = *storeidx;
    } else if (cr == 0) {
      /* The node is the first node */
      if (create >= 0)
        return (*storeidx);
      /* Ok, the node should be deleted then ... */
      RemoveEntryIdx(m, *storeidx);
      return (-1);
    } else {
      /* Node is somewhere in the list ? */
      idx = *storeidx;
      while (m->next[idx] >= 0 &&
             (memcmp(&(m->keys[m->next[idx] * m->key_size]), key,
                     m->key_size) <= 0)) {
        idx = m->next[idx];
      }
      if (!memcmp(&(m->keys[idx * m->key_size]), key, m->key_size)) {
        /* Ok, found and return */
        if (create >= 0)
          return (idx);
        RemoveEntryIdx(m, idx);
        return (-1);
      }
      /* Was not in the list */
      if (create <= 0)
        return (-1);
      {
        const indextype idx3 = idx;
        idx = NewEntry(m);
        m->next[idx] = m->next[idx3];
        m->next[idx3] = idx;
        m->prev[idx] = -idx3 - 1;
        if (m->next[idx] != -1)
          m->prev[m->next[idx]] = -idx - 1;
      }
    }
  }
  // assert(idx>=0);
  memcpy(&(m->keys[idx * m->key_size]), key, m->key_size);
  memcpy(&(m->data[idx * m->size_of_entry]), m->default_value,
         m->size_of_entry);
  return (idx);
}

struct matrix *CreateMatrixF(const int d, const indextype hash,
                             const float def_value) {
  return (CreateMatrix(sizeof(int) * d, hash, (void *)&def_value, FLOAT,
                       sizeof(float)));
}

struct matrix *CreateMatrixI(const int d, const indextype hash,
                             const int def_value) {
  return (CreateMatrix(sizeof(int) * d, hash, (void *)&def_value, INT,
                       sizeof(int)));
}

struct matrix *CreateMatrixD(const int d, const indextype hash,
                             const double def_value) {
  return (CreateMatrix(sizeof(int) * d, hash, (void *)&def_value, DOUBLE,
                       sizeof(double)));
}

void DeleteMatrix(struct matrix *m) {
  if (m == NULL)
    return;
  free(m->default_value);
  free(m->hash);
  free(m->keys);
  free(m->data);
  free(m->next);
  free(m->prev);
  free(m);
}

int qindcmp(const void *a, const void *b) {
  static struct matrix *m = NULL;

  /* Ugliness... If the first pointer is NULL, initialize
     the matrix pointer to second value. This is because
     we can only pass two paramters to this funciton to
     have it work with qsort */
  if (a)
    return (memcmp(&(m->keys[(*(indextype *)a) * m->key_size]),
                   &(m->keys[(*(indextype *)b) * m->key_size]), m->key_size));

  m = (struct matrix *)b;
  return (0);
}

int IncrementI(struct matrix *m, const int *indices, const int value) {
  indextype idx;
  int *v;
  assert(m->type_of_entry == INT);
  idx = FindEntry(m, (byte *)indices, 1);
  v = (int *)&(m->data[idx * m->size_of_entry]);

  if ((*v += value) != *((int *)(m->default_value)))
    return (*v);

  /* Default value, delete node */
  RemoveEntryIdx(m, idx);
  return (*((int *)(m->default_value)));
}

float IncrementF(struct matrix *m, const int *indices, const float value) {
  indextype idx;
  float *v;
  assert(m->type_of_entry == FLOAT);
  idx = FindEntry(m, (byte *)indices, 1);
  v = (float *)&(m->data[idx * m->size_of_entry]);

  if ((*v += value) != *((float *)(m->default_value)))
    return (*v);

  /* Default value, delete node */
  RemoveEntryIdx(m, idx);
  return (*((float *)(m->default_value)));
}

void ClearMatrix(struct matrix *m) {
  int i;
  m->num_entries = 0;
  for (i = 0; i < m->hashsize; i++) {
    m->hash[i] = -1;
  }
}

void SetValueI(struct matrix *m, const int *indices, const int value) {
  assert(m->type_of_entry == INT);
  SetRawValue(m, (byte *)indices, (void *)&value);
}

void SetValueF(struct matrix *m, const int *indices, const float value) {
  assert(m->type_of_entry == FLOAT);
  SetRawValue(m, (byte *)indices, (void *)&value);
}

void SetValueD(struct matrix *m, const int *indices, const double value) {
  assert(m->type_of_entry == DOUBLE);
  SetRawValue(m, (byte *)indices, (void *)&value);
}

int GetValueI(struct matrix *m, const int *indices) {
  int value;
  assert(m->type_of_entry == INT);
  GetRawValue(m, (byte *)indices, (void *)&value);
  return (value);
}

float GetValueF(struct matrix *m, const int *indices) {
  float value;
  assert(m->type_of_entry == FLOAT);
  GetRawValue(m, (byte *)indices, (void *)&value);
  return (value);
}

double GetValueD(struct matrix *m, const int *indices) {
  double value;
  assert(m->type_of_entry == DOUBLE);
  GetRawValue(m, (byte *)indices, (void *)&value);
  return (value);
}

/* The raw prototypes for *xStepThroughs cannot be inlined, since they
   have static variables. Inlined the wrappers instead */

static indextype sm_STidx;
static struct matrix *sm_STm;

void *StepThrough(struct matrix *mat, byte *key, void *data) {
  /* This function steps through given matrix. First invocation
     StepThrough(matrix,key,data) initializes data and returns
     a random value.. After this, StepThrough(NULL,key,data) returns
     values in order. Calling again with matrix pointer
     resets the counters to new data. After returning the last value,
     function returns NULL. It is users responsability to not call this
     function with NULL argument after this */
  // static struct matrix_node *mnp;
  // static indextype idx;
  // static struct matrix *m;

  if (mat != NULL) {
    sm_STm = mat;
    sm_STidx = -1;
    return (NULL);
  }

  if (++sm_STidx >= sm_STm->num_entries)
    return (NULL);

  /* Return indices */
  memcpy(key, &(sm_STm->keys[sm_STidx * sm_STm->key_size]), sm_STm->key_size);

  /* Return value */
  memcpy(data, &(sm_STm->data[sm_STidx * sm_STm->size_of_entry]),
         sm_STm->size_of_entry);
  return (&(sm_STm->data[sm_STidx * sm_STm->size_of_entry]));
}

void DeleteCurrentST() {
  assert(sm_STidx >= 0);
  RemoveEntryIdx(sm_STm, sm_STidx--);
}

void *OrderedStepThrough(struct matrix *m, byte *key, void *value) {
  indextype i;
  static indextype *sarray = NULL;
  static indextype count = 0;
  static struct matrix *sm;

  if (m) {
    sm = m;
    /* Create table with pointers to all of the nodes */
    if (sarray != NULL)
      free(sarray);
    sarray = (indextype *)malloc(m->num_entries * sizeof(indextype));
    for (i = 0; i < sm->num_entries; i++) {
      sarray[i] = i;
    }

    /* Quicksort the array */
    qindcmp(NULL, (void *)sm);
    qsort(sarray, sm->num_entries, sizeof(indextype), qindcmp);
    count = 0;
    return (NULL);
  }
  if (count != sm->num_entries) {
    memcpy(key, &(sm->keys[sarray[count] * sm->key_size]), sm->key_size);
    memcpy(value, (&(sm->data[sarray[count] * sm->size_of_entry])),
           sm->size_of_entry);
    return (&(sm->data[sarray[count++] * sm->size_of_entry]));
  }
  /* last entry, clean up */
  free(sarray);
  sarray = NULL;
  return (NULL);
}

void CheckConsistency(struct matrix *m) {
  /* This is for debugging */
  int i = 0, pn;
  for (i = 0; i < m->num_entries; i++) {
    /* Check the prev of next */
    if (m->next[i] > -1) {
      if (!(m->next[i] < m->num_entries)) {
        fprintf(stderr, "m_next[i] %d >= %d. i=%d. exit\n", (int)m->next[i],
                (int)m->num_entries, i);
        exit(-1);
      }
      pn = m->prev[m->next[i]];
      assert(pn < 0);
      if (-pn - 1 != i) {
        fprintf(stderr, "i%d next%d, pn%d. exit\n", i, (int)m->next[i],
                (int)m->prev[m->next[i]]);
        exit(-1);
      }
    }

    /* Check the next of prev */
    if (m->prev[i] >= 0) {
      // assert(m->hash[m->prev[i]]==i);
      if (!(m->hash[m->prev[i]] == i)) {
        fprintf(stderr, "m->hash[m->prev[i]]: %d != %d (m->prev[i] %d). exit\n",
                (int)m->hash[m->prev[i]], (int)i, (int)m->prev[i]);
        exit(-1);
      }
    } else {
      // assert(m->next[- m->prev[i] -1 ] == i );
      if (!(m->next[-m->prev[i] - 1] == i)) {
        fprintf(stderr, "-mp-1 %d, mn[-mp-1] %d, i %d. exit.\n",
                (int)(-m->prev[i] - 1), (int)m->next[-m->prev[i] - 1], (int)i);
        exit(-1);
      }
    }
  }
  /* Check the consistency of the first nodes */
  for (i = 0; i < m->hashsize; i++) {
    if (m->hash[i] == -1)
      continue;
    assert(m->hash[i] < m->num_entries);
    assert(m->prev[m->hash[i]] == i);
  }
}

void showvalues(struct matrix *m, indextype i) {
  int a = -999;
  fprintf(stderr, "%d: prev%d next%d, max%d, ", (int)i, (int)m->prev[i],
          (int)m->next[i], (int)m->num_entries);
  if (m->prev[i] >= 0)
    a = m->hash[m->prev[i]];
  else
    a = m->next[-m->prev[i] - 1];
  fprintf(stderr, "consistency p%d ", (int)a);
  if (m->next[i] > 0)
    fprintf(stderr, "n%d ", (int)(-m->prev[m->next[i]] - 1));
  fprintf(stderr, "\n");
}

void RemoveDefaultValues(struct matrix *m) {
  indextype i;
  for (i = 0; i < m->num_entries; i++) {
    if (!memcmp(&(m->data[i * m->size_of_entry]), m->default_value,
                m->size_of_entry))
      RemoveEntryIdx(m, i);
  }
}

#if defined no_inline_funcs
#define sm_inline_void void
#define sm_inline_int int
#define sm_inline_float float
#define sm_inline_double double
#define sm_inline_indextype indextype
#include "sparse_matrix_pinline.h"
#endif

/*
 * Local Variables:
 * mode: c
 * make-backup-files: t
 * End:
 */

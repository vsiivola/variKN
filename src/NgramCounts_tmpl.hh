// Helper library for storing and modifying the n-gram counts
#include "NgramCounts.hh"
#include "str.hh"
#include <limits.h>

template <typename K, typename T>
NgramCounts_t<K, T>::NgramCounts_t(const int n, const int max_vocab,
                                   indextype hashsize) {
  if (n <= 0) {
    fprintf(stderr, "Impossible n (%d). Exit.\n", n);
    exit(-1);
  }
  if (max_vocab)
    m_max_vocab = max_vocab;
  else
    m_max_vocab = 5000000;

  for (int i = 0; i < n; i++)
    m_ds.push_back(m_max_vocab);

  if (!hashsize)
    hashsize = 6000000;
  counts = new sikMatrix<K, T>(n, hashsize, 0);
}

template <typename K, typename T> NgramCounts_t<K, T>::~NgramCounts_t() {
  delete counts;
}

template <typename K, typename T>
long NgramCounts_t<K, T>::count(FILE *file, bool grow_vocab) {
  char charbuf[MAX_WLEN + 1];
  long n_entries = 0;
  long num_read = 0;
  K idx;

  while (fscanf(file, MAX_WLEN_FMT_STRING, charbuf) != EOF) {
    if (grow_vocab)
      idx = vocab->add_word(charbuf);
    else
      idx = vocab->word_index(charbuf);

    if (idx >= m_max_vocab - 1) {
      fprintf(stderr,
              "Exceeded maximum vocab size %d.\n"
              "Please increase the max size\n",
              m_max_vocab);
      exit(-1);
    }

    /* Cycle indices */
    for (int i = 0; i < m_ds.size() - 1; i++) {
      m_ds[i] = m_ds[i + 1];
    }
    m_ds[m_ds.size() - 1] = idx;

    if (++n_entries >= m_ds.size()) {
      counts->increment(&m_ds[0], 1);
      num_read++;
    }

    // if (n_entries%10000==0) fprintf(stderr,"\rread %d entries  ",n_entries);
  }
  return (num_read);
}

template <typename K, typename T>
long NgramCounts_t<K, T>::count(Storage_t<K, T> *data) {
  long n_entries = 0;
  long num_read = 0;
  K idx;

  for (size_t di = 0; di < data->size(); di++) {
    idx = data->data(di);

    if (idx >= m_max_vocab - 1) {
      fprintf(stderr,
              "Datastorage has different vocab than ngramcounts. Exit.\n");
      exit(-1);
    }

    /* Cycle indices */
    for (int i = 0; i < m_ds.size() - 1; i++) {
      m_ds[i] = m_ds[i + 1];
    }
    m_ds[m_ds.size() - 1] = idx;

    if (++n_entries >= m_ds.size()) {
      counts->increment(&m_ds[0], 1);
      num_read++;
    }
  }
  return (num_read);
}

template <typename K, typename T>
void NgramCounts_t<K, T>::write(FILE *out, FILE *vocabout, bool sort) {
  T value;

  if (sort)
    counts->stept_sortsearch = true;

  counts->stept(true, &m_ds[0], &value);
  while (counts->stept(false, &m_ds[0], &value)) {
    if (!vocabout) {
      for (int i = 0; i < m_ds.size(); i++) {
        fprintf(out, "%s ", vocab->word(m_ds[i]).c_str());
      }
    } else {
      /* indice 0 is for UNK */
      for (int i = 0; i < m_ds.size(); i++) {
        fprintf(out, "%d ", m_ds[i]);
      }
    }
    write_num(out, value);
    fprintf(out, "\n");
  }
  /* Write vocab if requested */
  if (vocabout)
    vocab->write(vocabout);
}

template <typename K, typename T>
void NgramCounts_t<K, T>::read_vocab(FILE *vocabfile) {
  std::string sbuf;
  while (str::read_line(&sbuf, vocabfile, true)) {
    int idx = vocab->add_word(sbuf);
    // fprintf(stderr,"Word %s, idx %d size %d\n",sbuf,idx,num_words());
    if (idx >= m_max_vocab - 1) {
      fprintf(stderr,
              "Exceeded maximum vocab size %d.\n"
              "Please increase the max size\n",
              m_max_vocab);
      exit(-1);
    }
  }
}

template <typename K, typename T>
void NgramCounts_t<K, T>::read(FILE *countsfile, FILE *vocabfile) {
  /* Read in files. Convert all indices to words and the
     the words to inidices corresponding to the internale
     vocabulary. This makes merging different countfiles simpler, but
     slows everything down a bit. */

  char sbuf[MAX_WLEN + 1], *cptr;

  if (vocabfile)
    read_vocab(vocabfile);

  char tmpchar[MAX_WLEN + 1];
  while (fgets(sbuf, MAX_WLEN, countsfile)) {
    // fprintf(stderr,"Looking at line %s",sbuf);
    cptr = strtok(sbuf, " ");
    for (int i = 0; i < m_ds.size(); i++) {
      if (!cptr) {
        fprintf(stderr, "Problem with input, skip entry %s\n", sbuf);
        break;
      }
      if (vocabfile) {
        // fprintf(stderr,"Adding cptr %s :",cptr);
        // fprintf(stderr,"%s\n", tmp_vocab.word(atoi(cptr)).c_str());
        int tmp = atoi(cptr);
        if (tmp >= vocab->num_words()) {
          fprintf(stderr, "Errors in input, word index %d unknown\n", tmp);
          continue;
        }
        m_ds[i] = tmp;
      } else {
        sscanf(cptr, "%s", tmpchar);
        m_ds[i] = vocab->add_word(tmpchar);
      }
      cptr = strtok(NULL, " ");
    }
    counts->increment(&m_ds[0], ascii2num(cptr, (T *)NULL));
  }
}

/* Don't mix ordered and unordered stepthroughs, you'll have problems */

template <typename K, typename T>
void NgramCounts_t<K, T>::shrink(float ndrop, int nfirst) {
  if (nfirst == 0)
    nfirst = INT_MAX;
  if (vocab->num_words() < nfirst && ndrop <= 0)
    return;

  std::vector<K> v(order());
  T count;

  std::vector<struct sortstruct> wordtable(vocab->num_words());
  Vocabulary newvoc;

  for (int i = 0; i < vocab->num_words(); i++)
    wordtable[i].wstring = vocab->word(i);
  /* Get the word counts */
  long num_grams = 0;
  counts->stepthrough(true, &v[0], &count);
  while (counts->stepthrough(false, &v[0], &count)) {
    for (int i = 0; i < order(); i++) {
      wordtable[v[i]].oldidx = v[i];
      wordtable[v[i]].count += count;
    }
    num_grams++;
  }

  /* sort by counts */
  std::sort(wordtable.begin() + 1, wordtable.end());

  /* create a backwards mapping */
  std::vector<int> backwards(wordtable.size(), 0);
  int i;
  for (i = 1; (i <= nfirst && wordtable[i].count > ndrop); i++)
    backwards[wordtable[i].oldidx] = i;
  nfirst = i - 1;

  /* Redo the vocabulary */
  vocab->set_oov(wordtable[0].wstring); // clears the current vocabulary
  for (i = 1; i <= nfirst; i++) {
    vocab->add_word(wordtable[i].wstring);
  }

  /* Redo the counts */
  for (i = 0; i < order(); i++)
    v[i] = nfirst + 1;
  sikMatrix<K, T> *new_counts =
      new sikMatrix<K, T>(i, std::min(num_grams + 1, (long)10000019), 0);
  counts->stepthrough(true, &v[0], &count);
  while (counts->stepthrough(false, &v[0], &count)) {
    for (int i = 0; i < order(); i++)
      v[i] = backwards[v[i]];
    new_counts->increment(&v[0], count);
  }

  /* Replace the current matrix */
  delete counts;
  counts = new_counts;
}

#if 0
template <typename K, typename T>
void NgramCounts_t<K, T>::convert_clustered(ClusterMap<K> *clmap) {

  std::vector<K> v(order());
  T count;
  sikMatrix<K, T> *new_counts = 
    new sikMatrix<K, T>(counts->dims, counts->num_entries(), 0);
  /* Get the word counts */
  counts->stepthrough(true,&v[0],&count);
  while (counts->stepthrough(false,&v[0],&count)) {
    clmap->wv2cv(v);
    new_counts->increment(&v[0], count);
  }

  delete counts;
  counts=new_counts;
}
#endif

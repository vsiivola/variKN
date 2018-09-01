// Library for storing and modifying the n-gram counts

/**********************************************************************
  These functions are for top level class (BOT abstracted away)
**********************************************************************/

template <typename KT, typename CT>
CT MultiOrderCounts<KT, CT>::GetCount(const std::vector<KT> &v) {
  return (GetCount(v.size(), &v[0]));
}

template <typename KT, typename CT>
CT MultiOrderCounts<KT, CT>::IncrementCount(const std::vector<KT> &v,
                                            const CT value) {
  return (IncrementCount(v.size(), &v[0], value));
}

/***********************************************************************
  Functions for the lower level template implementing most of the class
**********************************************************************/

template <typename KT, typename CT>
long MultiOrderCounts<KT, CT>::InitializeCountsFromText(
    FILE *in, Vocabulary *vocab, const bool grow_vocab, const int read_order,
    const std::string &sent_start_sym) {
  char charbuf[MAX_WLEN + 1];
  long num_read = 0;
  int sent_start_idx;
  KT idx;
  std::vector<KT> v;

  if (grow_vocab) {
    if (sent_start_sym.size())
      sent_start_idx = vocab->add_word(sent_start_sym);
    else
      sent_start_idx = -1;
    vocabsize = 64000;
  } else {
    vocabsize = vocab->num_words();
    if (!sent_start_sym.size())
      sent_start_idx = -1;
    else if (!(sent_start_idx = vocab->word_index(sent_start_sym))) {
      fprintf(stderr, "No sentence start symbol %s in vocabulary, exit.\n",
              sent_start_sym.c_str());
      exit(-1);
    }
  }

  while (fscanf(in, MAX_WLEN_FMT_STRING, charbuf) != EOF) {
    num_read++;
    // if (num_read % 1000000 == 0) fprintf(stderr,"Read %lld
    // words\n",num_read);
    if (grow_vocab)
      idx = vocab->add_word(charbuf);
    else
      idx = vocab->word_index(charbuf);

    // fprintf(stderr,"Word %s = %d\n", charbuf, idx);

    if (idx == sent_start_idx)
      v.clear();

    if (v.size() < read_order)
      v.push_back(idx);
    else
      v.back() = idx;

    IncrementCount(v, 1);
    // fprintf(stderr,"Cadd ");print_indices(stderr,v);fprintf(stderr," %d\n",
    // 1);

    if (v.size() == read_order)
      for (int j = 0; j < read_order - 1; j++)
        v[j] = v[j + 1];
  }
  fprintf(stderr, "Finished reading %ld words.\n", num_read);
  return (num_read);
}

template <typename KT, typename CT>
long MultiOrderCounts<KT, CT>::InitializeCountsFromStorage(
    Storage_t<KT, CT> *data, const int read_order, const int sent_start_idx) {
  long num_read = 0;
  std::vector<KT> v;
  KT idx;
  // fprintf(stderr,"Init from storage\n");

  for (size_t di = 0; di < data->size(); di++) {
    // fprintf(stderr,"adding %d: %d\n", di, data->data(di));
    num_read++;
    idx = data->data(di);
    if (idx == sent_start_idx)
      v.clear();

    if (v.size() < read_order)
      v.push_back(idx);
    else
      v.back() = idx;

    IncrementCount(v, 1);

    if (v.size() == read_order)
      for (int j = 0; j < read_order - 1; j++)
        v[j] = v[j + 1];
  }
  return (num_read);
}

template <typename KT, typename CT>
CT MultiOrderCounts<KT, CT>::Increment_wo_del(struct matrix *m,
                                              const KT *indices,
                                              const CT value) {

  indextype idx = FindEntry(m, (byte *)indices, 1);
  CT *valp = (CT *)&(m->data[idx * m->size_of_entry]);
  *valp += value;
  return (*valp);
}

template <typename KT, typename CT>
void MultiOrderCounts<KT, CT>::SetCount(const int order, const KT *v,
                                        const CT value) {
  allocate_matrices_counts(order);
  m_counts[order]->setvalue(v, value);
}

template <typename KT, typename CT>
CT MultiOrderCounts<KT, CT>::GetCount(const int order, const KT *v) {
  if (order >= m_counts.size())
    return (0);
  return (m_counts[order]->getvalue(v));
}

template <typename KT, typename CT>
CT MultiOrderCounts<KT, CT>::IncrementCount(const int order, const KT *v,
                                            const CT value) {
  allocate_matrices_counts(order);
  return (Increment_wo_del(m_counts[order]->m, v, value));
}

template <typename KT, typename CT>
void *MultiOrderCounts<KT, CT>::StepCountsOrder(const bool init,
                                                const int order, KT *indices,
                                                CT *value) {
  if (order >= m_counts.size())
    return (NULL);
  return (m_counts[order]->stepthrough(init, indices, value));
}

template <typename KT, typename CT>
void *MultiOrderCounts<KT, CT>::OrderedStepCountsOrder(const bool init,
                                                       const int order,
                                                       KT *indices, CT *value) {
  if (order >= m_counts.size())
    return (NULL);
  return (m_counts[order]->ordered_stepthrough(init, indices, value));
}

template <typename KT, typename CT>
CT MultiOrderCounts<KT, CT>::GetBackoffDen(const std::vector<KT> &v) {
  return (GetBackoffDen(v.size(), &v[0]));
}

template <typename KT, typename CT>
CT MultiOrderCounts<KT, CT>::GetBackoffNzer(const std::vector<KT> &v) {
  return (GetBackoffNzer(v.size(), &v[0]));
}

template <typename KT, typename CT>
MultiOrderCounts<KT, CT>::~MultiOrderCounts() {
  for (int i = 1; i < m_counts.size(); i++)
    if (std::find(m_do_not_delete.begin(), m_do_not_delete.end(), i) ==
        m_do_not_delete.end()) // The matrix was internally malloced
      delete m_counts[i];
}

template <typename KT, typename CT>
void MultiOrderCounts<KT, CT>::allocate_matrices_counts(int o) {
  if (o < m_counts.size())
    return;
  if (vocabsize == 0) {
    fprintf(
        stderr,
        "MultiOrderCounts: Please set a reasonable vocabulary size. Exit.\n");
    exit(-1);
  }
  if (hashsize == 0)
    hashsize = 600000;
  indextype real_hashsize;

  int old_size = m_counts.size();
  m_counts.resize(o + 1, NULL);
  for (int i = std::max(1, old_size); i < m_counts.size(); i++) {
    assert(m_counts[i] == NULL);
    std::vector<KT> v(i, (KT)vocabsize);
    // Some heuristics to try to get reasonable hash sizes. Not too well tested
    real_hashsize = std::min(
        std::max(1000, (indextype)(vocabsize * pow((float)i, 3))), hashsize);
    // fprintf(stderr,"min(%d*%d^3=%.0f,", vocabsize,i, vocabsize*pow((float) i,
    // 3)); fprintf(stderr,"%d)\n", (int) hashsize);
    if (i > 4 && order_size(i - 1) > 1)
      real_hashsize = order_size(i - 1) * 2 + 1;
    // if (i>2) fprintf(stderr,"Allocating counts matrices for order %d, size %d
    // (prev size %d, vocabsize %d)\n", i, real_hashsize, order_size(i-1),
    // this->vocabsize);
    m_counts[i] = new sikMatrix<KT, CT>(i, real_hashsize, 0);
  }
}

template <typename KT, typename CT>
bool MultiOrderCounts<KT, CT>::NextVector(std::vector<KT> &v) {
  if (m_cur_ng >= m_counts[m_cur_order]->num_entries()) {
    m_cur_ng = 0;
    m_cur_order++;
    while (m_cur_order < m_counts.size() &&
           m_counts[m_cur_order]->num_entries() == 0)
      m_cur_order++;
    if (m_cur_order >= m_counts.size()) {
      m_cur_order = 1;
      m_cur_ng = 0;
      return (false);
    }
  }
  v.resize(m_cur_order);
  const KT *keys = m_counts[m_cur_order]->Idx2Keyp(m_cur_ng);
  for (int i = 0; i < m_cur_order; i++) {
    v[i] = keys[i];
  }
  m_cur_ng++;
  return (true);
}

template <typename KT, typename CT>
CT MultiOrderCounts<KT, CT>::IncrementCountCache(const int order,
                                                 const KT *indices,
                                                 const CT value) {
  allocate_matrices_counts(order);
  CT *v;
  c_cache.resize(c_cache.size() + 1);
  c_cache_t &c = c_cache.back();

  c.order = order;
  c.val = value;

  struct matrix *m = m_counts[order]->m;
  const indextype idx = FindEntry(m, (byte *)indices, 1);
  c.index = idx;

  v = m_counts[order]->Idx2Valp(idx); //(CT *) &(m->data[idx*m->size_of_entry]);
  *v += value;
  return (*v);
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::WriteCounts(FILE *out) {
  std::vector<KT> v;
  CT val;
  BOT bo_val;

  fprintf(out, "\\vocabsize %d\n", MultiOrderCounts<KT, CT>::vocabsize);
  for (int o = 1; o <= order(); o++) {
    v.resize(o);
    fprintf(out, "\\%d-gram counts\n", o);
    this->StepCountsOrder(true, o, &v[0], &val);
    while (this->StepCountsOrder(false, o, &v[0], &val)) {
      if (val == 0)
        continue;
      print_indices(out, v);
      fprintf(out, " ");
      this->write_num(out, val);
      fprintf(out, "\n");
    }
    v.resize(o - 1);
    fprintf(out, "\n\\%d-gram backoffs\n", o);
    if (o == 1) {
      fprintf(out, "[ ] ");
      WriteCounts_BOhelper(out, &m_uni_bo);
      fprintf(out, "\n\n");
      continue;
    }
    StepBackoffsOrder(true, o, &v[0], &bo_val);
    while (StepBackoffsOrder(false, o, &v[0], &bo_val)) {
      if (bo_val.den == 0)
        continue;
      print_indices(out, v);
      fprintf(out, " ");
      WriteCounts_BOhelper(out, &bo_val);
      fprintf(out, "\n");
    }
    fprintf(out, "\n");
  }
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::ReadCounts(FILE *in) {
  int order = -1;
  std::string line;
  std::vector<KT> v;
  std::vector<std::string> vec;
  CT val;
  BOT bo_val;
  long lineno = 0;
  bool read_backoffs = false;
  bool ok = true;

  if (!str::read_line(&line, in, true) ||
      1 != sscanf(line.c_str(), "\\vocabsize %d", &(this->vocabsize))) {
    fprintf(stderr, "Error reading counts file.\nExit.\n");
    exit(-1);
  }

  while (str::read_line(&line, in, true)) {
    // fprintf(stderr,"READ %s\n", line.c_str());
    lineno++;

    // Skip empty lines
    if (line.find_first_not_of(" \t\n") == line.npos)
      continue;

    if (line[0] == '\\') {
      str::split(&line, "-", true, &vec);
      vec[0][0] = ' ';
      order = str::str2long(&vec[0], &ok);
      // fprintf(stderr,"ordered %s=%d(%d)\n", vec[0].c_str(), order, ok);
      // fprintf(stderr,"comparestring *%s*\n", vec[1].c_str());
      if (!strcmp("gram counts", vec[1].c_str())) {
        read_backoffs = false;
        v.resize(order);
      } else if (!strcmp("gram backoffs", vec[1].c_str())) {
        read_backoffs = true;
        v.resize(order - 1);
      } else
        ok = false;
      if (!ok) {
        fprintf(stderr, "Error reading line %ld.\n%s\nExit.\n", lineno,
                line.c_str());
        exit(-1);
      }
      continue;
    }

    str::split(&line, " \t", true, &vec);
    // fprintf(stderr,"ORDER %d\n", order);
    if (!read_backoffs) {
      for (int i = 1; i <= order; i++) {
        v[i - 1] = str::str2long(&vec[i], &ok);
        // fprintf(stderr,"Putting %s=%d, %d\n", vec[i].c_str(), v[i-1], ok);
      }
      this->read_num(&val, &vec[order + 2], &ok);
      // fprintf(stderr,"Val %s=%d (%d)\n", vec[order+2].c_str(), val, ok);
      if (!ok) {
        fprintf(stderr, "Error reading line %ld:\n%s\nExit2.\n", lineno,
                line.c_str());
        exit(-1);
      }
      // fprintf(stderr,"Setcount ");print_indices(stderr,v);
      // fprintf(stderr,"val %d\n", val);
      this->SetCount(v, val);
      continue;
    }
    for (int i = 1; i < order; i++) {
      v[i - 1] = str::str2long(&vec[i], &ok);
    }
    ReadCounts_BOhelper(&bo_val, &vec[order + 1], &ok);
    if (!ok) {
      fprintf(stderr, "Error reading line %ld. Exit3.\n", lineno);
      exit(-1);
    }
    SetBackoff(order, &v[0], &bo_val);
  }
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::IncrementBackoffCache(
    const int order, const KT *indices, const BOT *value) {
  bo_cache.resize(bo_cache.size() + 1);
  bo_cache_t &b = bo_cache.back();
  b.order = order;
  b.bo = *value;

  if (order == 1) {
    m_uni_bo += *value;
    return;
  }

  allocate_matrices_backoffs(order);
  BOT *v;
  struct matrix *m = m_backoffs[order]->m;

  const indextype idx = FindEntry(m, (byte *)indices, 1);
  b.index = idx;

  v = (BOT *)&(m->data[idx * m->size_of_entry]);
  *v += *value;
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::ResetCaches() {
  this->c_cache.resize(0);
  bo_cache.resize(0);

  this->min_cc_cache.resize(this->m_counts.size() + 1);
  for (int i = 1; i < this->m_counts.size(); i++) {
    this->min_cc_cache[i] = this->m_counts[i]->num_entries();
    // fprintf(stderr,"cc[%d]=%d\n",i,m_counts[i]->num_entries);
  }
  this->min_cc_cache[this->m_counts.size()] = 0;

  min_bo_cache.resize(m_backoffs.size() + 1);
  for (int i = 2; i < m_backoffs.size(); i++) {
    min_bo_cache[i] = m_backoffs[i]->num_entries();
    // fprintf(stderr,"bo[%d]=%d\n",i,m_backoffs[i]->num_entries);
  }
  min_bo_cache[m_backoffs.size()] = 0;
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::UndoCached() {
  /* This could be speeded up by assuming that all cached
     values are of the same order
  */
  for (long i = this->c_cache.size() - 1; i >= 0; i--) {
    struct MultiOrderCounts<KT, CT>::c_cache_t &c = this->c_cache[i];
    struct matrix *m = this->m_counts[c.order]->m;
    *(CT *)(&(m->data[c.index * m->size_of_entry])) -= c.val;
  }

  for (int j = 1; j < this->m_counts.size(); j++) {
    for (long i = this->m_counts[j]->num_entries() - 1;
         i >= this->min_cc_cache[j]; i--) {
      RemoveEntryIdx(MultiOrderCounts<KT, CT>::m_counts[j]->m, i);
    }
  }

  /* Copy of the beginning of the func modified for backoffs.. */
  for (long i = bo_cache.size() - 1; i >= 0; i--) {
    struct bo_cache_t &c = bo_cache[i];
    if (c.order == 1) {
      m_uni_bo -= c.bo;
      continue;
    }
    if (c.index >= min_bo_cache[c.order])
      continue;
    struct matrix *m = m_backoffs[c.order]->m;
    BOT *bo = (BOT *)(&(m->data[c.index * m->size_of_entry]));
    *bo -= c.bo;
  }

  for (int j = 2; j < m_backoffs.size(); j++) {
    for (long i = m_backoffs[j]->num_entries() - 1; i >= min_bo_cache[j]; i--)
      RemoveEntryIdx(m_backoffs[j]->m, i);
  }
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::IncrementBackoff(
    const int order, const KT *v, const BOT *value) {
  if (order == 1) {
    m_uni_bo += *value;
    return;
  }
  allocate_matrices_backoffs(order);
  m_backoffs[order]->increment(v, value);
}

template <typename KT, typename CT, typename BOT>
CT MultiOrderCounts_Generic_BOT<KT, CT, BOT>::GetBackoffDen(const int order,
                                                            const KT *v) {
  if (order == 1)
    return (m_uni_bo.den);
  if (order >= m_backoffs.size())
    return (0);
  return (m_backoffs[order]->getvalue(v).den);
}

template <typename KT, typename CT>
void MultiOrderCounts<KT, CT>::UseAsCounts(sikMatrix<KT, CT> *mat) {
  /* Should get rid of m->dims, so that libsparsematrix can
     be cleaned up */
  allocate_matrices_counts(mat->dims);
  if (std::find(m_do_not_delete.begin(), m_do_not_delete.end(), mat->dims) ==
      m_do_not_delete.end()) // The matrix was internally malloced
    delete m_counts[mat->dims];
  m_counts[mat->dims] = mat;
  m_do_not_delete.push_back(mat->dims);
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::RemoveOrder(int order) {
  fprintf(stderr, "cur_o %d(%ld): ", order,
          (long)MultiOrderCounts<KT, CT>::m_counts.size() - 1);
  assert(order == (MultiOrderCounts<KT, CT>::m_counts.size() - 1));
  delete MultiOrderCounts<KT, CT>::m_counts.back();
  MultiOrderCounts<KT, CT>::m_counts.pop_back();

  assert(order == m_backoffs.size() - 1);
  delete m_backoffs.back();
  m_backoffs.pop_back();

  fprintf(stderr, "cur_o %d(%ld): ", order,
          (long)MultiOrderCounts<KT, CT>::m_counts.size() - 1);
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::IncrementBackoff(
    const std::vector<KT> &v, const BOT *value) {
  IncrementBackoff(v.size(), &v[0], value);
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::SetBackoff(const int order,
                                                           const KT *v,
                                                           const BOT *bo) {
  allocate_matrices_backoffs(order);
  if (order > 1) {
    m_backoffs[order]->setvalue(v, bo);
    return;
  }
  m_uni_bo = *bo;
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::SetBackoff(
    const std::vector<KT> &v, const BOT *bo) {
  SetBackoff(v.size() + 1, &v[0], bo);
}
template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::GetBackoff(const int order,
                                                           const KT *v,
                                                           BOT *value) {
  if (order >= m_backoffs.size()) {
    memcpy(value, &m_bb_init, sizeof(BOT));
    return;
  }
  if (order > 1) {
    m_backoffs[order]->getvalue(v, value);
    return;
  }
  memcpy(value, &m_uni_bo, sizeof(BOT));
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::IncrementBackoffCacheDen(
    const int order, const KT *v, const CT value) {
  BOT bo;
  zero_bo(bo);
  bo.den = value;
  IncrementBackoffCache(order, v, &bo);
}

template <typename KT, typename CT, typename BOT>
MultiOrderCounts_Generic_BOT<KT, CT, BOT>::MultiOrderCounts_Generic_BOT() {
  m_bb_init.den = 0;
  m_bb_init.nzer = 0;
  m_bb_init.prune_den = 0;
  m_bb_fp_init.den = 0;
  m_bb_fp_init.nzer = 0;
  m_bb_fp_init.lost = 0;
  m_bb_fp_init.lost_den = 0;
  m_3bb_init.den = 0;
  m_3bb_init.nzer[0] = 0;
  m_3bb_init.nzer[1] = 0;
  m_3bb_init.nzer[2] = 0;
  m_3bb_init.prune_den = 0;
  m_3bb_fp_init.den = 0;
  m_3bb_fp_init.nzer[0] = 0;
  m_3bb_fp_init.nzer[1] = 0;
  m_3bb_fp_init.nzer[2] = 0;
  m_3bb_fp_init.prune_den = 0;
  m_3bb_fp_init.prune_den_left[0] = 0;
  m_3bb_fp_init.prune_den_left[1] = 0;
  m_3bb_fp_init.prune_den_left[2] = 0;
  m_3bb_fp_init.prune_den_den[0] = 0;
  m_3bb_fp_init.prune_den_den[1] = 0;
  m_3bb_fp_init.prune_den_den[2] = 0;

  zero_bo(m_uni_bo);
  zero_bo(bo_init);
}

template <typename KT, typename CT, typename BOT>
MultiOrderCounts_Generic_BOT<KT, CT, BOT>::~MultiOrderCounts_Generic_BOT() {
  for (int i = 2; i < m_backoffs.size(); i++) {
    delete m_backoffs[i];
  }
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::IncrementBackoffDen(
    const int order, const KT *indices, const CT den) {
  allocate_matrices_backoffs(order);
  if (order > 1) {
    struct matrix *m = m_backoffs[order]->m;
    indextype idx = FindEntry(m, (byte *)indices, 1);
    BOT *bop = (BOT *)(&(m->data[idx * m->size_of_entry]));
    bop->den += den;
    // Delete node, if matches default value
    if (!memcmp(bop, m->default_value, m->size_of_entry))
      RemoveEntryIdx(m, idx);
    return;
  }
  m_uni_bo.den += den;
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::allocate_matrices_backoffs(
    int o) {
  if (o < m_backoffs.size())
    return;
  if (this->vocabsize == 0) {
    fprintf(
        stderr,
        "MultiOrderCounts_t: Please set a reasonable vocabulary size. Exit.\n");
    exit(-1);
  }
  if (this->hashsize == 0)
    this->hashsize = 300000;
  indextype real_hashsize;

  int old_size = m_backoffs.size();
  m_backoffs.resize(o + 1, NULL);
  for (int i = std::max(2, old_size); i < m_backoffs.size(); i++) {
    assert(m_backoffs[i] == NULL);
    // Some heuristics to try to get reasonable hash sizes. Not too well tested
    real_hashsize = std::min(
        std::max(1000, (indextype)(this->vocabsize * pow((float)i, 3))),
        this->hashsize);
    if (i > 4 && bo_order_size(i - 1) > 1)
      real_hashsize = bo_order_size(i - 1) * 2 + 1;
    fprintf(stderr, "Allocating backoff matrices for order %d, size %ld", i,
            (long)real_hashsize);
    if (i > 2)
      fprintf(stderr, "(prev_size %d, vocabsize %d)\n", bo_order_size(i - 1),
              this->vocabsize);
    else
      fprintf(stderr, "\n");
    m_backoffs[i] = new sikMatrix<KT, BOT>(i - 1, real_hashsize, bo_init);
    fprintf(stderr, "allocation succesful\n");
  }
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::clear_derived_counts() {
  zero_bo(m_uni_bo);
  this->m_counts[1]->clear();
  for (int i = 2; i < this->m_counts.size() - 1; i++) {
    this->m_counts[i]->clear();
    m_backoffs[i]->clear();
  }
  m_backoffs.back()->clear();
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::clear_nzer(int o) {
  if (o > 1) {
    for (indextype i = 0; i < m_backoffs[o]->num_entries(); i++)
      zero_nz(m_backoffs[o]->Idx2Valp(i));
    return;
  }
  zero_nz(&m_uni_bo);
}

template <typename KT, typename CT>
void MultiOrderCounts_1nzer_fp<KT, CT>::clear_lden(const int o) {
  if (o > 1) {
    for (indextype i = 0;
         i <
         MultiOrderCounts_Generic_BOT<
             KT, CT, MultiOrderCounts_counter_types::bo_c_fp<CT>>::m_backoffs[o]
             ->num_entries();
         i++) {
      MultiOrderCounts_counter_types::bo_c_fp<CT> *bop =
          MultiOrderCounts_Generic_BOT<
              KT, CT,
              MultiOrderCounts_counter_types::bo_c_fp<CT>>::m_backoffs[o]
              ->Idx2Valp(i);
      bop->lost_den = 0;
    }
    return;
  }
  MultiOrderCounts_Generic_BOT<
      KT, CT, MultiOrderCounts_counter_types::bo_c_fp<CT>>::m_uni_bo.lost_den =
      0;
}

template <typename KT, typename CT, typename BOT>
void *MultiOrderCounts_Generic_BOT<KT, CT, BOT>::StepBackoffsOrder(
    const bool init, const int order, KT *indices, BOT *value) {
  if (order >= MultiOrderCounts<KT, CT>::m_counts.size())
    return (NULL);
  assert(order >= 2);
  return (m_backoffs[order]->stepthrough(init, indices, value));
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::IncrementBackoffNzer_1nzer(
    const int order, const KT *indices, const CT nz) {
  if (order > 1) {
    allocate_matrices_backoffs(order);
    struct matrix *m = m_backoffs[order]->m;
    indextype idx = FindEntry(m, (byte *)indices, 1);
    BOT *bop = (BOT *)(&(m->data[idx * m->size_of_entry]));
    bop->nzer += nz;
    // Delete node, if matches default value
    if (!memcmp(bop, m->default_value, m->size_of_entry))
      RemoveEntryIdx(m, idx);
    return;
  }
  m_uni_bo.nzer += nz;
}

template <typename KT, typename CT, typename BOT>
CT MultiOrderCounts_Generic_BOT<KT, CT, BOT>::GetBackoffNzer_1nzer(
    const int order, const KT *v) {
  if (order == 1)
    return (m_uni_bo.nzer);
  if (order >= m_backoffs.size())
    return (0);
  return (m_backoffs[order]->getvalue(v).nzer);
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::IncrementBackoffCacheNzer_1nzer(
    const int order, const KT *v, const CT value) {
  BOT bo;
  zero_bo(bo);
  bo.nzer = value;
  IncrementBackoffCache(order, v, &bo);
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::GetBackoffNzer_3nzer(
    const int order, const KT *v, CT *res) {
  if (order == 1) {
    for (int i = 0; i < 3; i++)
      res[i] = m_uni_bo.nzer[i];
    return;
  }
  if (order >= m_backoffs.size()) {
    for (int i = 0; i < 3; i++)
      res[i] = 0;
    return;
  }
  BOT bo;
  m_backoffs[order]->getvalue(v, &bo);
  for (int i = 0; i < 3; i++)
    res[i] = bo.nzer[i];
}

template <typename KT, typename CT, typename BOT>
CT MultiOrderCounts_Generic_BOT<KT, CT, BOT>::GetBackoffNzer_3nzer(
    const int order, const KT *v, const int which) {
  CT res[3];
  GetBackoffNzer(order, v, res);
  return (res[which]);
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::IncrementBackoffCacheNzer_3nzer(
    const int order, const KT *v, const CT *value) {
  BOT bo;
  zero_bo(bo);
  bo.nzer[0] = value[0];
  bo.nzer[1] = value[1];
  bo.nzer[2] = value[2];
  IncrementBackoffCache(order, v, &bo);
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::IncrementBackoffCacheNzer_3nzer(
    const int order, const KT *v, const int pos, const CT value) {
  BOT bo;
  zero_bo(bo);
  bo.nzer[pos] = value;
  IncrementBackoffCache(order, v, &bo);
}

template <typename KT, typename CT, typename BOT>
void MultiOrderCounts_Generic_BOT<KT, CT, BOT>::IncrementBackoffNzer_3nzer(
    const int order, const KT *indices, const int pos, const CT nz) {
  if (order > 1) {
    struct matrix *m = m_backoffs[order]->m;
    indextype idx = FindEntry(m, (byte *)indices, 1);
    BOT *bop = (BOT *)(&(m->data[idx * m->size_of_entry]));
    bop->nzer[pos] += nz;

    // Delete node, if matches default value
    if (!memcmp(bop, m->default_value, m->size_of_entry))
      RemoveEntryIdx(m, idx);
    return;
  }
  m_uni_bo.nzer[pos] += nz;
}

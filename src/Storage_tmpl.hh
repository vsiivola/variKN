// Library for storing corpuses in memory
#include "def.hh"
#include "io.hh"
#include "str.hh"

template <typename T> void Storage<T>::read(FILE *in, Vocabulary &voc) {
  char buf[MAX_WLEN + 1];
  datavec.reserve(100000);
  // fprintf(stderr,"storage read:");
  while (fscanf(in, MAX_WLEN_FMT_STRING, buf) == 1) {
    datavec.push_back(voc.word_index(buf));
    // fprintf(stderr," %s(%d)", buf, voc.word_index(buf));
  }
  // fprintf(stderr,"\n");
}

template <typename T, typename ICT>
void Storage_t<T, ICT>::initialize_fast_search_lists(
    const int order, sikMatrix<T, ICT> *refmat, sikMatrix<T, ICT> *curmat) {
  assert(refmat != NULL);
  m_refmat = refmat;
  m_lists.clear();
  this->m_lists2.clear();
  m_last_init_mapped = (m_refmat->dims <= 2);
  if (m_last_init_mapped)
    m_lists.resize(m_refmat->num_entries());
  else
    this->m_lists2.resize(m_refmat->num_entries());
  std::vector<T> v;
  size_t i;

  if (this->clear_lm_history == -1) {
    for (i = 0; i < order - 1; i++)
      v.push_back(this->datavec[i]);
    v.push_back(0); // Make the vector right size
    for (; i < this->size(); i++) {
      const indextype idx = FindEntry(m_refmat->m, (byte *)(&v[0]), 0);
      v[order - 1] = this->datavec[i];
      if (idx >= 0) // The prefix exists ?
        if (!curmat || FindEntry(curmat->m, (byte *)(&v[0]), 0) == -1) {
          // fprintf(stderr,"Adding to list ");print_indices(v);
          // fprintf(stderr,": %ld", curmat); if (curmat) fprintf(stderr," %d",
          // curmat->getvalue(&v[0])); fprintf(stderr,"\n");
          if (m_last_init_mapped)
            m_lists[idx][this->datavec[i]]++;
          else
            this->m_lists2[idx].push_back(this->datavec[i]);
        }
      for (int j = 0; j < order - 1; j++)
        v[j] = v[j + 1];
    }
    return;
  }

  for (size_t i = 0; i < this->size(); i++) {
    if (this->datavec[i] == this->clear_lm_history) {
      v.clear();
    }

    if (v.size() < order)
      v.push_back(this->datavec[i]);
    else
      v[order - 1] = this->datavec[i];

    if (v.size() == order) {
      const indextype idx = FindEntry(m_refmat->m, (byte *)(&v[0]), 0);
      if (idx >= 0) { // The prefix exists ?
        // This is for growing-pruning iterations, slows down otherwise
        if (!curmat || FindEntry(curmat->m, (byte *)(&v[0]), 0) == -1) {
          // fprintf(stderr,"adding %d", datavec[i]);
          // if (curmat) fprintf(stderr,"(%d)", curmat->getvalue(&v[0]));
          // fprintf(stderr,"\n");
          if (m_last_init_mapped)
            m_lists[idx][this->datavec[i]]++;
          else
            this->m_lists2[idx].push_back(this->datavec[i]);
        }
      }
      for (int j = 0; j < order - 1; j++)
        v[j] = v[j + 1];
    }
  }
  return;
}

template <typename T, typename ICT>
void Storage_t<T, ICT>::initialize_fast_search_lists_for_pruning(
    const int order, sikMatrix<T, ICT> *refmat) {
  assert(refmat != NULL);
  m_refmat = refmat;
  m_last_init_mapped = false;
  m_lists.clear();
  this->m_lists2.clear();
  this->prune_lists.clear();
  this->prune_lists.resize(m_refmat->num_entries(), 0);
  std::vector<T> v;

  for (indextype i = 0; i < this->size(); i++) {
    if (this->datavec[i] == this->clear_lm_history) {
      v.clear();
    }

    if (v.size() < order)
      v.push_back(this->datavec[i]);
    else
      v[order - 1] = this->datavec[i];

    if (v.size() == order) {
      const indextype idx = FindEntry(m_refmat->m, (byte *)(&v[0]), 0);
      if (idx >= 0)
        this->prune_lists[idx] += 1;
      for (int j = 0; j < order - 1; j++)
        v[j] = v[j + 1];
    }
  }
  return;
}

template <typename T, typename ICT>
void Storage_t<T, ICT>::init_fsl_file(const int order,
                                      sikMatrix<T, ICT> *refmat,
                                      std::string &fname, Vocabulary *voc) {

  assert(refmat != NULL);
  m_refmat = refmat;
  m_last_init_mapped = (m_refmat->dims <= 2);
  m_lists.clear();
  this->m_lists2.clear();

  if (m_last_init_mapped)
    m_lists.resize(m_refmat->num_entries());
  else
    this->m_lists2.resize(m_refmat->num_entries());

  io::Stream in(fname, "r", io::REOPENABLE);
  std::vector<T> v;
  char sbuf[MAX_WLEN];

  if (this->clear_lm_history == -1) {
    for (int i = 0; i < order - 1; i++) {
      if (fscanf(in.file, "%s", sbuf) != 1)
        return;
      v.push_back(voc->word_index(sbuf));
    }
    v.push_back(0); // vector is of right size

    while (fscanf(in.file, "%s", sbuf) == 1) {
      indextype idx = FindEntry(m_refmat->m, (byte *)(&v[0]), 0);
      // fprintf(stderr,"%d: idx %d(<%d) ",i,idx,m_lists.size());
      T widx = voc->word_index(sbuf);
      v[order - 1] = widx;
      if (idx >= 0) { // The prefix exists ?
        if (m_last_init_mapped)
          m_lists[idx][widx]++;
        else
          this->m_lists2[idx].push_back(widx);
      }
      for (int j = 0; j < order - 1; j++)
        v[j] = v[j + 1];
    }
  } else {
    while (fscanf(in.file, "%s", sbuf) == 1) {
      T widx = voc->word_index(sbuf);
      if (widx == this->clear_lm_history)
        v.clear();
      if (v.size() < order)
        v.push_back(widx);
      else
        v[order - 1] = widx;

      if (v.size() == order) {
        indextype idx = FindEntry(m_refmat->m, (byte *)(&v[0]), 0);
        // fprintf(stderr,"%d: idx %d(<%d) ",i,idx,m_lists.size());
        if (idx >= 0) { // The prefix exists ?
          if (m_last_init_mapped)
            m_lists[idx][widx]++;
          else
            this->m_lists2[idx].push_back(widx);
        }
        for (int j = 0; j < order - 1; j++)
          v[j] = v[j + 1];
      }
    }
  }
}

template <typename T, typename ICT>
void Storage_t<T, ICT>::fast_search_next(std::vector<T> *v, int *ridx,
                                         ICT *rval) {
  *ridx = -1;

  /* if v!=NULL, init new search */
  if (v != NULL) {
    // fprintf(stderr,"a%d",(*v)[0]);
    m_cur_vec = FindEntry(m_refmat->m, (byte *)(&(*v)[0]), 0);
    // fprintf(stderr,"MC %d(%ld)\n",m_cur_vec,m_lists[m_cur_vec].size());
    if (m_last_init_mapped) {
      // FIXME: Why is the next if clause needed ?
      if (m_cur_vec >= m_lists.size()) {
        m_cur_vec = -1;
        return;
      }
      m_cur_vec_idx = m_lists[m_cur_vec].begin();
    } else
      this->m_cur_vec_idx2 = 0;
    return;
  }
  // fprintf(stderr,"b%d_%d",m_cur_vec_idx,m_lists[m_cur_vec].size());

  if (m_cur_vec == -1)
    return;
  if (m_last_init_mapped && (m_cur_vec >= m_lists.size() ||
                             m_cur_vec_idx == m_lists[m_cur_vec].end()))
    return;
  if (!m_last_init_mapped &&
      (m_cur_vec >= this->m_lists2.size() ||
       this->m_cur_vec_idx2 >= this->m_lists2[m_cur_vec].size()))
    return;

  if (m_last_init_mapped) {
    *ridx = m_cur_vec_idx->first;
    *rval = m_cur_vec_idx->second;
    m_cur_vec_idx++;
  } else {
    *ridx = this->m_lists2[m_cur_vec][this->m_cur_vec_idx2++];
    *rval = 1;
  }
  return;
}

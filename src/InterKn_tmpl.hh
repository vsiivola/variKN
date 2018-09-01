// The main library for the n-gram model estimation
#include "GramSorter.hh"
#include <limits>
#ifdef _WIN32
inline double log2(double n) { return log(n) / log(2.0); }
#endif

/******************************************************************/
/* Inlined functions differing between float, int backoff and     */
/* int3 backoff                                                   */
/******************************************************************/

template <typename KT, typename ICT>
void InterKn_int_disc<KT, ICT>::init_disc(float x) {
  std::vector<float> di(this->m_order, x);
  flatv2disc(di);
}

template <typename KT, typename ICT>
void InterKn_int_disc3<KT, ICT>::init_disc(float x) {
  std::vector<float> di(this->m_order * 3, x);
  flatv2disc(di);
}

template <typename KT, typename ICT>
void InterKn_int_disc<KT, ICT>::disc2flatv(std::vector<float> &v) {
  v.resize(this->m_order);
  for (int i = 0; i < this->m_order; i++)
    v[i] = m_discount[i + 1];
}

template <typename KT, typename ICT>
void InterKn_int_disc3<KT, ICT>::disc2flatv(std::vector<float> &v) {
  v.resize(this->m_order * 3);
  for (int i = 1; i < m_discount.size(); i++)
    for (int j = 0; j < 3; j++)
      v[(i - 1) * 3 + j] = m_discount[i][j];
}

template <typename KT, typename ICT>
float InterKn_int_disc<KT, ICT>::flatv2disc(std::vector<float> &v) {
  float debug_overthrow = 0.0; // This is no longer needed, should remove
  fprintf(stderr, "[");
  for (int i = 1; i <= this->m_order; i++) {
    m_discount[i] = v[i - 1];
    fprintf(stderr, " %f", m_discount[i]);
    if (m_discount[i] < 0.0) {
      debug_overthrow -= m_discount[i];
      m_discount[i] = 0.0;
    } else if (m_discount[i] > 1.0) {
      debug_overthrow += m_discount[i] - 1.0;
      m_discount[i] = 1.0;
    }
  }
  fprintf(stderr, "]");
  return (debug_overthrow);
}

template <typename KT, typename ICT>
float InterKn_int_disc3<KT, ICT>::flatv2disc(std::vector<float> &v) {
  float debug_overthrow = 0.0; // This is no longer needed, should remove
  for (int i = 1; i <= this->m_order; i++) {
    fprintf(stderr, "[");
    for (int j = 0; j < 3; j++) {
      m_discount[i][j] = v[(i - 1) * 3 + j];
      fprintf(stderr, " %.2f", m_discount[i][j]);
      if (m_discount[i][j] < 0.0) {
        debug_overthrow -= m_discount[i][j];
        m_discount[i][j] = 0.0;
      } else if (m_discount[i][j] > j + 1) {
        debug_overthrow += m_discount[i][j] - j - 1;
        m_discount[i][j] = j + 1;
      }
    }
    fprintf(stderr, "]");
  }
  return (debug_overthrow);
}

/************************************************************
 Non-inlined functions for the different child classes      *
*************************************************************/

template <typename KT, typename ICT>
InterKn_int_disc<KT, ICT>::InterKn_int_disc(
    const bool absolute, const std::string data, const std::string vocab,
    const std::string optisource, const int read_counts, const int order,
    const int ndrop, const int nfirst, Storage_t<KT, ICT> *datastorage,
    const std::string prunedata_name, const std::string sent_boundary,
    const indextype hashsize)
    : InterKn_t<KT, ICT>(absolute, data, optisource, prunedata_name) {
  this->moc = new MultiOrderCounts_1nzer<KT, ICT>;
  this->constructor_helper(vocab, read_counts, order, ndrop, nfirst,
                           datastorage, hashsize, sent_boundary);
}

template <typename KT, typename ICT>
InterKn_int_disc3<KT, ICT>::InterKn_int_disc3(
    const bool absolute, const std::string data, const std::string vocab,
    const std::string optisource, const int read_counts, const int order,
    const int ndrop, const int nfirst, Storage_t<KT, ICT> *datastorage,
    const std::string prunedata_name, const std::string sent_boundary,
    const indextype hashsize)
    : InterKn_t<KT, ICT>(absolute, data, optisource, prunedata_name) {
  this->moc = new MultiOrderCounts_3nzer<KT, ICT>;
  this->constructor_helper(vocab, read_counts, order, ndrop, nfirst,
                           datastorage, hashsize, sent_boundary);
}

template <typename KT, typename ICT>
void InterKn_int_disc<KT, ICT>::estimate_nzer_counts() {
  std::vector<KT> v(this->m_order);
  ICT value;

  for (int o = 1; o <= this->m_order; o++) {
    this->moc->StepCountsOrder(true, o, &v[0], &value);
    while (this->moc->StepCountsOrder(false, o, &v[0], &value)) {
      if (value == 0)
        continue;
      this->moc->IncrementBackoffNzer(o, &v[0], 1);
    }
  }
}

template <typename KT, typename ICT>
std::vector<float> InterKn_t<KT, ICT>::calculate_leaveoneout_discounts(
    int order, std::vector<float> cur_disc) {
  const int num_d_coeffs(cur_disc.size());
  std::vector<ICT> count_of_counts(num_d_coeffs + 1);
  std::vector<KT> v(order);
  ICT value;
  std::vector<float> result(cur_disc);

  this->moc->StepCountsOrder(true, order, &v[0], &value);
  while (this->moc->StepCountsOrder(false, order, &v[0], &value)) {
    if (value > num_d_coeffs + 1 || value < 1)
      continue;
    count_of_counts[value - 1] += 1;
  }
  // Leave-one-out estimates for discounts
  if (count_of_counts[0] == 0 || count_of_counts[1] == 0) {
    fprintf(stderr,
            "Count of counts zero, skipping leave-one-out estimation.\n");
    return cur_disc;
  }

  float Y =
      count_of_counts[0] / (count_of_counts[0] + 2.0f * count_of_counts[1]);
  // fprintf(stderr, "Y = %d / (%d + 2 * %d) = %.2f\n", count_of_counts[0],
  // count_of_counts[0], count_of_counts[1], Y);
  fprintf(stderr, "set loo_disc order %d -> [", order);
  for (int i = 0; i < num_d_coeffs; ++i) {
    if (count_of_counts[i] > 0) {
      result[i] = std::max(
          0.1f, std::min(i + 1 - 0.2f,
                         i + 1 -
                             float((i + 2) * Y * count_of_counts[i + 1]) /
                                 count_of_counts[i]));
      // fprintf(stderr, " {%d - (%d * %.2f * %d) / %d} = ",
      //        i+1, i+2, Y, count_of_counts[i+1], count_of_counts[i]);
    }
    fprintf(stderr, " %.2f", result[i]);
  }
  fprintf(stderr, " ]\n");
  return result;
}

template <typename KT, typename ICT>
void InterKn_int_disc3<KT, ICT>::estimate_nzer_counts() {
  std::vector<KT> v(this->m_order);
  ICT value;

  for (int o = 1; o <= this->m_order; o++) {
    this->moc->StepCountsOrder(true, o, &v[0], &value);
    while (this->moc->StepCountsOrder(false, o, &v[0], &value)) {
      /*
      assert(o>1 || value>0);
      if (value<=0) {
        fprintf(stderr,"Asserted failed order %d:",o);
        for (int i=0;i<o;i++) {
          fprintf(stderr," %d",v[i]);
        }
        fprintf(stderr," =%d\n",value);
        exit(-1);
      }
*/

      if (value == 0) { // With clear hist this is needed
        // assert(m_sent_boundary!=-1);
        continue;
      }
      this->moc->IncrementBackoffNzer(o, &v[0], std::min((ICT)2, value - 1), 1);
    }
  }
}

template <typename KT, typename ICT>
template <typename BOT>
void InterKn_int_disc<KT, ICT>::remove_sent_start_prob_fbase(BOT *dummy) {
  std::vector<KT> tmp(1, (KT)this->m_sent_boundary);
  const ICT value = this->moc->GetCount(tmp);
  this->moc->IncrementCount(tmp, -value);
  BOT bo;
  this->moc->zero_bo(bo);
  bo.den = -value;
  bo.nzer = -1;
  this->moc->IncrementBackoff(1, NULL, &bo);
}

template <typename KT, typename ICT>
template <typename BOT>
void InterKn_int_disc3<KT, ICT>::remove_sent_start_prob_fbase(BOT *dummy) {
  std::vector<KT> tmp(1, (KT)this->m_sent_boundary);
  const ICT value = this->moc->GetCount(tmp);
  this->moc->IncrementCount(tmp, -value);
  BOT bo;
  this->moc->zero_bo(bo);
  bo.den = -value;
  bo.nzer[std::min((ICT)2, value - 1)] = -1;
  this->moc->IncrementBackoff(1, NULL, &bo);
}

template <typename KT, typename ICT>
float InterKn_int_disc<KT, ICT>::kn_prob(const int order, const KT *i,
                                         const ICT num) {
  double prob = 0;
  if (order == 1) {
    // Unigram smoothing
    prob = (m_discount[1] * this->moc->GetBackoffNzer(1, NULL)) /
           (float)this->moc->GetBackoffDen(1, NULL) /
           ((float)this->vocab.num_words());
    // fprintf(stderr,"Unigram smooth %.4f\n", prob);
  }

  if (num <= 0)
    return (prob);
  const ICT den = this->moc->GetBackoffDen(order, i);
  if (den <= 0)
    return (prob);
  prob += (num - m_discount[order]) / den;
  // fprintf(stderr,"knp: (%ld - %.4f) / %ld = %f -> %f\n", (long) num,
  // m_discount[order], (long) den, (num-m_discount[order])/den, prob);
  return prob;
}

template <typename KT, typename ICT>
float InterKn_int_disc3<KT, ICT>::kn_prob(const int order, const KT *i,
                                          const ICT num) {
  double prob = 0;
  if (order == 1) {
    typename MultiOrderCounts<KT, ICT>::bo_3c bo;
    this->moc->GetBackoff(1, i, &bo);

    prob = (m_discount[1][0] * bo.nzer[0] + m_discount[1][1] * bo.nzer[1] +
            m_discount[1][2] * bo.nzer[2] + bo.prune_den) /
           (float)(bo.den) / ((float)this->vocab.num_words());
  }
  if (num <= 0)
    return (prob);
  const float den = this->moc->GetBackoffDen(order, i);
  if (den <= 0)
    return (prob);
  prob += (num - m_discount[order][std::min((ICT)2, num - 1)]) / den;
  return (prob);
}

template <typename KT, typename ICT>
float InterKn_int_disc<KT, ICT>::kn_coeff(const int order, const KT *i) {
  if (order > this->m_order)
    return (1.0);

  MultiOrderCounts_counter_types::bo_c<ICT> bb;
  this->moc->GetBackoff(order, i, &bb);
  // fprintf(stderr," BO (%.5f +%d)/ %d\n", m_discount[order] *bb.nzer,
  // bb.prune_den ,bb.den);
  if (bb.den)
    return ((bb.nzer * m_discount[order] + bb.prune_den) / (float)bb.den);
  return (1.0);
}

template <typename KT, typename ICT>
template <typename BOT>
float InterKn_int_disc3<KT, ICT>::kn_coeff_3nzer(const int order, const KT *i,
                                                 const BOT *dummy) {
  if (order > this->m_order)
    return (1.0);

  BOT bb;
  this->moc->GetBackoff(order, i, &bb);
  if (!bb.den)
    return (1.0);

  if (bb.den < 0 ||
      (bb.den > 0 && (bb.nzer[0] < 0 || bb.nzer[1] < 0 || bb.nzer[2] < 0))) {
    fprintf(stderr, "Weird values %ld [%ld %ld %ld], bailing out\n",
            (long)bb.den, (long)bb.nzer[0], (long)bb.nzer[1], (long)bb.nzer[2]);
    exit(-1);
  }

  // fprintf(stderr,"%d [%d %d %d]\n", bb.den,
  //	    bb.nzer[0], bb.nzer[1], bb.nzer[2]);
  if (bb.den)
    return ((bb.nzer[0] * m_discount[order][0] +
             bb.nzer[1] * m_discount[order][1] +
             bb.nzer[2] * m_discount[order][2] + bb.prune_den) /
            (float)bb.den);
  return (1.0);
}

template <typename KT, typename CT>
void InterKn_t<KT, CT>::create_model(float prunetreshold) {
  // if (m_sent_boundary>0) clear_lm_sentence_boundaries();
  for (int i = 1; i <= this->m_order; ++i) {
    this->set_leaveoneout_discounts(i);
  }

  find_coeffs(-0.1 * m_ori_treshold, 1e-1, 5e-2);
  if (prunetreshold > 0.0 || this->cutoffs.size() ||
      this->discard_ngrams_with_unk) {
    if (!this->prune_with_real_counts) {
      prune_model(prunetreshold, 1, (Storage_t<KT, CT> *)NULL);
    } else {
      io::Stream::verbose = true;
      io::Stream in(this->m_prunedata_name, "r", io::REOPENABLE);
      Storage_t<KT, CT> real_counts;
      real_counts.read(in.file, this->vocab);
      fprintf(stderr, "real counts size %ld\n", (long)real_counts.size());
      prune_model(prunetreshold, 1, &real_counts);
    }
  }
  find_coeffs(-0.1 * m_ori_treshold, 1e-1, 5e-2);
}

template <typename KT, typename CT>
template <typename BOT>
void InterKn_t<KT, CT>::add_counts_for_backoffs_fbase(BOT *dummy) {
  float coeff;
  for (int o = 1; o < this->m_order; o++) {
    BOT bo;
    std::vector<KT> gram(o);
    moc->StepBackoffsOrder(true, o + 1, &gram[0], &bo);
    while (moc->StepBackoffsOrder(false, o + 1, &gram[0], &bo)) {
      coeff = safelogprob(kn_coeff(o + 1, &gram[0]));
      if (coeff >= -1e-3)
        continue;
      FindEntry(moc->m_counts[o]->m, (byte *)(&gram[0]), 1);
    }
  }
}

template <typename KT, typename CT>
template <typename BOT>
void InterKn_t<KT, CT>::prune_model_fbase(float treshold, bool recorrect_kn,
                                          Storage_t<KT, CT> *real_counts,
                                          BOT *dummy) {
  std::vector<KT> v;
  CT num;
  float logprobdelta, safelogprob_mult;

  treshold = treshold * this->model_cost_scale;

  this->set_order(moc->order());
  if (this->m_absolute_discounting)
    recorrect_kn = false;

  for (int o = this->order(); o >= 2; o--) {
    if (real_counts) {
      fprintf(stderr, "Using real counts\n");
      real_counts->initialize_fast_search_lists_for_pruning(o,
                                                            moc->m_counts[o]);
    }
    fprintf(stderr, "Pruning order %d\n", o);
    v.resize(o);
    moc->StepCountsOrder(true, o, &v[0], &num);
    while (moc->StepCountsOrder(false, 0, &v[0], &num)) {
      if (num == 0) {
        // fprintf(stderr,"SKIPPING "); print_indices(stderr,v);
        // fprintf(stderr,"\n");
        if (!real_counts)
          moc->DeleteCurrentST(o);
        continue;
      }
      assert(num > 0);
      moc->ResetCaches();
      if (num <= this->cutoff(o)) {
        if (this->discard_cutoffs)
          prune_gram(v, num, false, (BOT *)NULL);
        else
          prune_gram(v, num, recorrect_kn, (BOT *)NULL);
        if (!real_counts)
          moc->DeleteCurrentST(o);
        continue;
      }

      if (this->discard_ngrams_with_unk) {
        bool flag = false;
        for (int ngu = 0; ngu < o; ngu++) {
          if (v[ngu] == 0) {
            flag = true;
            break;
          }
        }
        if (flag) {
          prune_gram(v, num, recorrect_kn, (BOT *)NULL);
          if (!real_counts)
            moc->DeleteCurrentST(o);
          continue;
        }
      }
      if (!(num > 0))
        fprintf(stderr, "Weird num %ld\n", (long)num);
      if (real_counts) {
        const indextype idx =
            FindEntry(moc->m_counts[o]->m, (byte *)(&v[0]), 0);
        // print_indices(stderr,v); fprintf(stderr,"IDX %d\n", idx);
        safelogprob_mult = real_counts->prune_lists[idx];
      } else if (!this->m_ehist_estimate) {
        safelogprob_mult = num;
      } else {
        safelogprob_mult = this->m_ehist_estimate;
        for (int o2 = v.size() - 1; o2 >= 1; o2--) {
          std::vector<KT> g2(v);
          g2.resize(o2);
          safelogprob_mult *= tableprob(g2);
        }
      }
      logprobdelta = safelogprob(tableprob(v));

      /* Try removing the gram */
      prune_gram(v, num, recorrect_kn, (BOT *)NULL);

      logprobdelta -= safelogprob(tableprob(v));
      logprobdelta *= safelogprob_mult;

      if (logprobdelta - treshold > 0) {
        moc->UndoCached();

      } else {
        if (!real_counts)
          moc->DeleteCurrentST(o);
      }
    }

    // Clean up for realcounts
    if (real_counts) {
      moc->StepCountsOrder(true, o, &v[0], &num);
      while (moc->StepCountsOrder(false, 0, &v[0], &num)) {
        if (num == 0)
          moc->DeleteCurrentST(o);
      }
    }
  }
  // Clean backoffs
  moc->RemoveDefaultBackoffs();
  if (real_counts)
    real_counts->prune_lists.clear();
}

template <typename KT, typename ICT>
void InterKn_int_disc<KT, ICT>::prune_gram(
    std::vector<KT> &v, ICT num, bool recorrect_kn,
    MultiOrderCounts_counter_types::bo_c<ICT> *dummy) {
  MultiOrderCounts_counter_types::bo_c<ICT> bo;
  const int o = v.size();

  this->moc->IncrementCountCache(o, &v[0], -num);
  this->moc->GetBackoff(o, &v[0], &bo);
  if (bo.den == bo.prune_den + num) {
    bo.den = -bo.prune_den - num;
    bo.prune_den = -bo.prune_den;
  } else {
    bo.den = 0;
    bo.prune_den = num;
  }
  bo.nzer = -1;
  this->moc->IncrementBackoffCache(o, &v[0], &bo);

  if (recorrect_kn && num != 1) {
    ICT orig = this->moc->IncrementCountCache(o - 1, &v[1], num - 1);
    if (orig == num - 1) {
      this->moc->IncrementCountCache(o - 1, &v[1], -num + 1);
    } else {
      bo.den = num - 1;
      bo.nzer = 0;
      bo.prune_den = 0;
      this->moc->IncrementBackoffCache(o - 1, &v[1], &bo);
    }
  }
}

template <typename KT, typename ICT>
void InterKn_int_disc3<KT, ICT>::prune_gram(
    std::vector<KT> &v, ICT num, bool recorrect_kn,
    MultiOrderCounts_counter_types::bo_3c<ICT> *dummy) {

  MultiOrderCounts_counter_types::bo_3c<ICT> bo;
  const int o = v.size();

  this->moc->IncrementCountCache(o, &v[0], -num);
  this->moc->GetBackoff(o, &v[0], &bo);
  if (bo.den == bo.prune_den + num) {
    bo.den = -bo.prune_den - num;
    bo.prune_den = -bo.prune_den;
  } else {
    bo.den = 0;
    bo.prune_den = num;
  }
  bo.nzer[0] = 0;
  bo.nzer[1] = 0;
  bo.nzer[2] = 0;
  bo.nzer[std::min(num - 1, (ICT)2)] = -1;
  this->moc->IncrementBackoffCache(o, &v[0], &bo);

  if (recorrect_kn && num != 1) {
    ICT origl = this->moc->IncrementCountCache(o - 1, &v[1], num - 1);
    if (origl == num - 1) {
      this->moc->IncrementCountCache(o - 1, &v[1], -num + 1);
    } else {
      /* Take care of the possible changing nzer class */
      bo.den = num - 1;
      bo.nzer[0] = 0;
      bo.nzer[1] = 0;
      bo.nzer[2] = 0;
      bo.prune_den = 0;
      bo.nzer[std::min((ICT)2, origl - num)] = -1;
      bo.nzer[std::min((ICT)2, origl - 1)] += 1;
      this->moc->IncrementBackoffCache(o - 1, &v[1], &bo);
    }
  }
}

template <typename KT, typename CT>
template <typename BOT>
void InterKn_t<KT, CT>::add_zeroprob_grams_fbase(BOT *dummy) {
  std::vector<KT> v;
  CT num;
  if (std::numeric_limits<KT>::max() < this->vocab.num_words()) {
    fprintf(stderr, "Too large vocab for the given key type!!! Abort.\n");
    exit(-1);
  }

  this->set_order(moc->order());

  for (int o = this->m_order; o >= 2; o--) {
    v.resize(o);
    this->moc->StepCountsOrder(true, o, &v[0], &num);
    while (this->moc->StepCountsOrder(false, o, &v[0], &num)) {
      this->moc->IncrementCount(o - 1, &v[0], 0);
    }
    BOT b;

    if (o == 2)
      continue;
    this->moc->StepBackoffsOrder(true, o, &v[0], &b);
    while (this->moc->StepBackoffsOrder(false, o, &v[0], &b))
      if (b.den > 0)
        this->moc->IncrementCount(o - 1, &v[0], 0);
  }

  // Add zeroprob unigrams:
  for (KT i = 0; i < this->vocab.num_words(); i++) {
    // fprintf(stderr, "%d/%d %s\n", i, this->vocab.num_words(),
    // this->vocab.word(i).c_str());
    this->moc->IncrementCount(1, &i, 0);
  }
}

template <typename KT, typename CT>
void InterKn_t<KT, CT>::counts2asciilm(FILE *out) {
  fprintf(stderr, "Warning, writing interpolated format to arpa. "
                  "Don't do this unless you know what you are doing. "
                  "You can save this model in binary format and use bin2arpa "
                  "to turn it to arpa format. "
                  "\"arpa2arpa model.fakearpa model.realarpa\" "
                  "should convert the model "
                  "to real arpa.\n");

  /*********************************************/
  /* Put smoothed model into ngram            */
  /*******************************************/
  std::vector<std::string> strbuf;
  strbuf.reserve(100000);

  fprintf(out, "\\interpolated\n");
  fprintf(out, "\\data\\\n");

  remove_zeroprob_grams();
  if (this->zeroprobgrams)
    add_zeroprob_grams();

  // Insert to counts MINLOGPROB values for each required backoff
  add_counts_for_backoffs();

  std::vector<int> num_grams(this->m_order + 1, 0);
  // A dry run to figure out the number of gram that will be printed
  for (int o = 1; o <= this->m_order; o++) {
    std::vector<KT> gram(o);
    float logprob, bo_den;
    CT num;
    moc->StepCountsOrder(true, o, &gram[0], &num);
    while (moc->StepCountsOrder(false, o, &gram[0], &num)) {
      logprob = safelogprob(kn_prob(o, &gram[0], num));
      bo_den = moc->GetBackoffDen(o + 1, &gram[0]);
      if (!this->zeroprobgrams && (logprob <= MINLOGPROB && bo_den <= 0))
        continue;
      num_grams[o]++;
    }
  }

  for (int i = 1; i <= this->m_order; i++) {
    fprintf(out, "ngram %d=%d\n", i, num_grams[i]);
  }

  // The real stuff. Print.
  for (int o = 1; o <= this->m_order; o++) {
    fprintf(out, "\n\\%d-grams:\n", o);
    std::vector<KT> gram(o);
    float logprob, bo_den;
    CT num;
    moc->StepCountsOrder(true, o, &gram[0], &num);
    while (moc->StepCountsOrder(false, o, &gram[0], &num)) {
      logprob = safelogprob(kn_prob(o, &gram[0], num));
      bo_den = moc->GetBackoffDen(o + 1, &gram[0]);
      // fprintf(stderr,"consider %d -> %g", moc->GetBackoffDen(o+1, &gram[0]),
      //      safelogprob(kn_coeff(o+1,&gram[0])));
      if (!this->zeroprobgrams && (logprob <= MINLOGPROB && bo_den <= 0))
        continue;
      fprintf(out, "%.4f", logprob);
      for (int i = 0; i < o; i++)
        fprintf(out, " %s", this->vocab.word(gram[i]).c_str());
      if (bo_den > 0) {
        fprintf(out, " %.4f", safelogprob(kn_coeff(o + 1, &gram[0])));
        // fprintf(stderr," ok");
      }
      // fprintf(stderr,"\n");
      fprintf(out, "\n");
    }
  }
  fprintf(out, "\\end\\\n");
}

/************************************************************
  Inlined functions for the template implementing most of the class
***********************************************************/

template <typename KT, typename CT>
float InterKn_t<KT, CT>::kn_prob(const int order, const KT *i) {
  return (kn_prob(order, i, moc->GetCount(order, i)));
}

/****************************************************************************/
/* Noninlined functions for template InterKn_t                              */
/****************************************************************************/
#ifdef sgi
#define log2 1 / flog10(2) * flog10
#endif

template <typename KT, typename CT>
void InterKn_t<KT, CT>::constructor_helper(
    const std::string &vocabname, const int read_counts, const int order,
    const int ndrop, const int nfirst_i, Storage_t<KT, CT> *datastorage,
    const indextype hashsize, const std::string &sent_start_symbol) {

  int nfirst = nfirst_i;
  if (nfirst == -1)
    nfirst = 9999999;
  moc->hashsize = hashsize;
  io::Stream::verbose = true;
  this->input_data_size = 0;

  /* Should we read the counts from a file? In that case, we also need
     a vocabulary*/
  if (read_counts) {
    if (vocabname.size() == 0) {
      fprintf(stderr, "Can't use counts without a vocabulary. Exit.\n");
      exit(-1);
    }
    io::Stream datain(this->m_data_name, "r");
    io::Stream vocabin(vocabname, "r");
    /* if read_counts==-1 the counts file contains only the highest order
       counts and other counts should be estimated later */
    if (read_counts == -1) {
      if (this->vocab.num_words() > 1)
        fprintf(stderr, "Warning: something is going wrong. The indices must "
                        "match, this will not be checked.\n");
      fprintf(stderr, "Reading counts for the highest order\n");
      NgramCounts_t<KT, CT> *ng = new NgramCounts_t<KT, CT>(order, 0, hashsize);
      ng->read(datain.file, vocabin.file);
      this->moc->UseAsCounts(ng->counts);
    } else {
      /* The counts include all orders and have already been KN-discounted */
      fprintf(stderr, "Reading previously KN-discounted counts\n");
      this->vocab.read(vocabin.file);
      this->moc->ReadCounts(datain.file);
    }
    this->moc->vocabsize = this->vocab.num_words();
  } else {
    fprintf(stderr, "Reading data\n");
    if (vocabname.length()) {
      fprintf(stderr, "Using vocab %s\n", vocabname.c_str());
      if (this->vocab.num_words() > 1)
        fprintf(stderr, "Warning: something is going wrong. The vocabularies "
                        "must be the same (not checked)\n");
      if (this->vocab.num_words() > 65534) {
        fprintf(stderr, "Too big vocabulary for --smallvocab (%d). Exit.\n",
                this->vocab.num_words());
        exit(-1);
      }
      io::Stream vocabin(vocabname, "r");
      this->vocab.read(vocabin.file);
      if (this->vocab.num_words() < 1) {
        fprintf(stderr, "Warning: no words from vocab file? Exit\n");
        exit(-1);
      }
    }
    if (this->vocab.num_words() > 1) {
      fprintf(stderr, "Restricted vocab\n");
      io::Stream datain(this->m_data_name, "r");
      if (!datastorage) {
        // fprintf(stderr,"Init from text\n");
        this->input_data_size = moc->InitializeCountsFromText(
            datain.file, &(this->vocab), false, order, sent_start_symbol);
      } else {
        // fprintf(stderr,"Init from file\n");
        datastorage->read(datain.file, this->vocab);
        // fprintf(stderr,"All read\n");
        int sent_start_idx = this->vocab.word_index(sent_start_symbol);
        if (sent_start_idx == 0 && sent_start_symbol.length()) {
          fprintf(stderr, "No sentence start %s(len %d) in vocab, exit.\n",
                  sent_start_symbol.c_str(), (int)sent_start_symbol.length());
          exit(-1);
        }
        this->input_data_size = moc->InitializeCountsFromStorage(
            datastorage, order, sent_start_idx);
      }
    } else {
      if (this->KT_is_short((KT *)NULL) || ndrop != 0 ||
          (nfirst < 9999999 && nfirst >= 1)) {
        if (ndrop != 0 || (nfirst < 9999999 && nfirst >= 1)) {
          if (vocabname.length()) {
            fprintf(
                stderr,
                "ndrop or nfirst may not be specified with vocabin. Exit.\n");
            exit(-1);
          }
        }
        NgramCounts_t<int, CT> nc_tmp(1, 0, 500000);
        io::Stream datain(this->m_data_name, "r");
        nc_tmp.count(datain.file, true);
        datain.close();
        fprintf(stderr, "Shrinking\n");
        if (!this->KT_is_short((KT *)NULL))
          nc_tmp.shrink(ndrop, nfirst);
        else {
          nc_tmp.shrink(0, std::min(65534, nfirst));
        }
        nc_tmp.vocab->copy_vocab_to(this->vocab);
        datain.open(this->m_data_name, "r");
        if (!datastorage) {
          this->input_data_size = moc->InitializeCountsFromText(
              datain.file, &(this->vocab), false, order, sent_start_symbol);
        } else {
          datastorage->read(datain.file, this->vocab);
          int sent_start_idx = this->vocab.word_index(sent_start_symbol);
          if (sent_start_idx == 0 && sent_start_symbol.length()) {
            fprintf(stderr, "No sentence start %s(len %d) in vocab, exit.\n",
                    sent_start_symbol.c_str(), (int)sent_start_symbol.length());
            exit(-1);
          }
          this->input_data_size = moc->InitializeCountsFromStorage(
              datastorage, order, sent_start_idx);
        }
      } else {
        io::Stream datain(this->m_data_name, "r");
        this->input_data_size = moc->InitializeCountsFromText(
            datain.file, &(this->vocab), true, order, sent_start_symbol);
        datain.close();
        if (datastorage) {
          datain.open(this->m_data_name, "r");
          datastorage->read(datain.file, this->vocab);
        }
      }
    }
  }

  this->set_order(moc->order());

  if (read_counts != 1) {
    fprintf(stderr, "Estimating bo counts\n");
    estimate_bo_counts(true);
    fprintf(stderr, "Estimating nzer counts\n");
    this->estimate_nzer_counts();
  }

  this->m_optistorage = new Storage<KT>;
  if (this->m_opti_name.length()) {
    fprintf(stderr, "Reading optisource\n");
    io::Stream in(this->m_opti_name, "r");
    this->m_optistorage->read(in.file, this->vocab);
  }
  fprintf(stderr, "Optistorage size %ld\n", (long)this->m_optistorage->size());

  if (sent_start_symbol.size())
    this->m_sent_boundary = this->vocab.word_index(sent_start_symbol);

  this->model_cost_scale = log2(this->vocab.num_words()) + 2 * 10;
}

template <typename KT, typename CT> InterKn_t<KT, CT>::~InterKn_t() {
  delete this->m_optistorage;
  delete moc;
}

template <typename KT, typename ICT>
void InterKn_int_disc<KT, ICT>::set_order(int o) {
  int old_order = this->m_order;
  this->m_order = o;
  m_discount.resize(this->m_order + 1);
  for (int i = old_order + 1; i <= this->m_order; i++) {
    m_discount[i] = m_discount[old_order];
  }
  initialize_minmax();
}

template <typename KT, typename ICT>
void InterKn_int_disc3<KT, ICT>::set_order(int o) {
  int old_order = this->m_order;
  this->m_order = o;
  m_discount.resize(this->m_order + 1);
  for (int i = old_order + 1; i <= this->m_order; i++) {
    m_discount[i] = m_discount[old_order];
  }
  initialize_minmax();
}

template <typename KT, typename CT>
void InterKn_t<KT, CT>::counts2lm(TreeGram *lm) {
  /*********************************************/
  /* Put smoothed model into ngram            */
  /*******************************************/

  remove_zeroprob_grams();
  add_zeroprob_grams();
  add_counts_for_backoffs();

  lm->set_type(TreeGram::INTERPOLATED);
  TreeGram::Gram gr;
  float prob, coeff;
  std::vector<KT> v(1);
  CT num;

  this->vocab.copy_vocab_to(*lm);

  indextype nodes = 0;
  for (int o = 1; o <= this->m_order; o++) {
    // fprintf(stderr,"O%d siex %d\n", o, moc->order_size(o));
    nodes += moc->order_size(o);
  }
  lm->reserve_nodes(nodes);
  // fprintf(stderr,"Reserved %d nodes\n", nodes);

  for (int o = 1; o <= this->m_order; o++) {
    // fprintf(stderr,"Adding o %d\n",o);
    bool breaker = true;
    v.resize(o);
    gr.resize(o);
    GramSorter gramsorter(o, moc->order_size(o));
    moc->StepCountsOrder(true, o, &v[0], &num);
    while (moc->StepCountsOrder(false, o, &v[0], &num)) {
      prob = kn_prob(o, &v[0], num);
      coeff = kn_coeff(o + 1, &v[0]);
      for (int i = 0; i < o; i++)
        gr[i] = v[i];
      // fprintf(stderr,"to sorter: %.4f
      // ",safelogprob(prob));print_indices(stderr, v); fprintf(stderr,"
      // %.4f\n", safelogprob(coeff));
      gramsorter.add_gram(gr, safelogprob(prob), safelogprob2(coeff));
      breaker = false;
    }
    if (breaker)
      break;
    gramsorter.sort();
    for (size_t i = 0; i < gramsorter.num_grams(); i++) {
      GramSorter::Data data = gramsorter.data(i);
      gr = gramsorter.gram(i);
      // fprintf(stderr,"adding ");print_indices(stderr, gr);fprintf(stderr,"
      // %.4f %.4f\n", data.log_prob, data.back_off);
      lm->add_gram(gr, data.log_prob, data.back_off);
    }
  }
  lm->finalize();
}

template <typename KT, typename CT>
void InterKn_t<KT, CT>::counts2ascii(FILE *out) {
  moc->WriteCounts(out);
}

template <typename KT, typename CT>
float InterKn_t<KT, CT>::evaluate(std::vector<float> &discounts) {
  float tmp = m_eval_cache->getvalue(&discounts[0]);
  if (tmp != 1.0)
    return tmp;

  float debug_overthrow = flatv2disc(discounts);
  fprintf(stderr, ": ");
  this->re_estimate_needed();
  double lp = logprob_datastorage(*(this->m_optistorage));
  // double lp=logprob_file(m_opti_name);
  lp += debug_overthrow;
  fprintf(stderr, "%g\n", lp / log10(2.0));
  m_eval_cache->setvalue(&discounts[0], lp);
  return (lp);
}

template <typename KT, typename CT>
double InterKn_t<KT, CT>::logprob_file(const char *name) {
  double logprob = 0.0;

  std::vector<KT> indices;
  indices.reserve(this->m_order);
  io::Stream f(name, "r");
  char w[MAX_WLEN + 1];
  long nwords = 0;
  while (fscanf(f.file, MAX_WLEN_FMT_STRING, w) == 1) {
    const KT idx = this->m_ng->word_index(w);
    if (idx == this->m_sent_boundary) {
      indices.clear();
      indices.push_back(idx);
      continue;
    } else if (indices.size() < this->m_order)
      indices.push_back(idx);
    else {
      for (int i = 0; i < this->m_order - 1; i++)
        indices[i] = indices[i + 1];
      indices[this->m_order - 1] = idx;
      // fprintf(stderr,"read %s(%d)\n",w,indices[m_order-1]);
    }
    if (indices.back() == 0)
      continue; // Do not optimize for unks...
    nwords++;
    // fprintf(stderr,"%s %.1g\n",w,tableprob(indices));
    logprob += safelogprob(tableprob(indices));
  }
  f.close();
  return (-logprob / nwords);
}

template <typename KT, typename CT>
double InterKn_t<KT, CT>::logprob_datastorage(const Storage<KT> &data) {
  double logprob = 0.0;

  std::vector<KT> indices;
  indices.reserve(this->m_order);
  size_t nwords = 0;
  size_t target = 0;

  if (this->m_sent_boundary == -1) {
    for (size_t j = 0; j < data.size(); j++) {
      if (target < this->m_order) {
        target++;
        indices.resize(target);
      }
      for (int i = 1; i <= target; i++) {
        indices[target - i] = data.data(j + 1 - i);
      }
      // print_indices(stderr, indices);
      // fprintf(stderr," %.1g\n",tableprob(indices));
      if (indices.back() == 0)
        continue; // Do not optimize for unks...
      nwords++;
      logprob += safelogprob(tableprob(indices));
    }
    return (-logprob / nwords);
  }
  // Ugly code duplicaton, needs fixing
  for (size_t j = 0; j < data.size(); j++) {
    const KT idx = data.data(j);
    if (idx == this->m_sent_boundary) {
      indices.clear();
      indices.push_back(idx);
      continue;
    } else if (indices.size() < this->m_order)
      indices.push_back(idx);
    else {
      for (int i = 0; i < this->m_order - 1; i++)
        indices[i] = indices[i + 1];
      indices[this->m_order - 1] = idx;
      // fprintf(stderr,"read %s(%d)\n",w,indices[m_order-1]);
    }
    if (indices.back() == 0)
      continue; // Do not optimize for unks...
    nwords++;
    // print_indices(indices);
    // fprintf(stderr,"%.1g\n",tableprob(indices));
    logprob += safelogprob(tableprob(indices));
  }
  return (-logprob / nwords);
}

template <typename KT, typename CT>
double InterKn_t<KT, CT>::tableprob(std::vector<KT> &indices) {
  double prob = 0.0;
  KT *iptr;

  // fprintf(stderr,"looking
  // ");print_indices(stderr,indices);fprintf(stderr,"\n");
  const int looptill = std::min(indices.size(), (size_t)this->m_order);
  for (int n = 1; n <= looptill; n++) {
    iptr = &(indices.back()) - n + 1;
    if (n > 1) {
      // fprintf(stderr,"prob %.4f * coeff %.4f = ", prob, kn_coeff(n, iptr));
      prob *= kn_coeff(n, iptr);
      // fprintf(stderr,"%.4f\n", prob);
    }
    // fprintf(stderr,"oldprob %.4f + prob %.4f = ", prob, kn_prob(n, iptr));
    prob += kn_prob(n, iptr);
    // fprintf(stderr,"%.4f\n",prob);
  }
  // fprintf(stderr,"vc: return %e\n",prob);
  assert(prob >= -1e-03 && prob <= 1.001);
  return (prob);
}

template <typename KT, typename CT>
void InterKn_t<KT, CT>::find_coeffs(float brak, float precision,
                                    float lin_precision) {
  fprintf(stderr, "Gots opti size %ld\n", this->m_optistorage->size());
  if (this->m_optistorage->size() == 0) {
    fprintf(stderr, "Skipping numerical parameter optimization, no "
                    "optimization set specified\n");
    return;
  }
  std::vector<float> searchstart;
  disc2flatv(searchstart);

  /* This is set up for evaluate() the cache the already obtained values */
  m_eval_cache = new sikMatrix<float, float>(searchstart.size(), 1000, 1.0);
  QFit qfit(lin_precision, precision, this);
  qfit.set_initial_point(searchstart);
  qfit.set_minimum(this->m_minvals);
  qfit.set_maximum(this->m_maxvals);

  // std::vector<float> brak(searchstart.size(), brak);
  // qfit.set_searchstartlim(brak);

  searchstart = qfit.minimize(1000 * searchstart.size());
  delete m_eval_cache;

  /* Construct the full model to memory */
  fprintf(stderr, "Optimal discounts: ");
  flatv2disc(searchstart);
  fprintf(stderr, "\n");
  this->re_estimate_needed();
}

template <typename KT, typename CT>
void InterKn_t<KT, CT>::estimate_bo_counts(bool zero_sentstart) {
  if (this->m_absolute_discounting) {
    estimate_bo_counts_absolute_discounting();
    return;
  }

  /* Estimate backoff and lower order counts            */
  std::vector<KT> v(this->m_order);
  CT value;

  if (this->m_sent_boundary < 0) {
    for (int o = this->m_order; o >= 1; o--) {
      moc->StepCountsOrder(true, o, &v[0], &value);
      while (moc->StepCountsOrder(false, o, &v[0], &value)) {
        moc->IncrementBackoffDen(o, &v[0], value);
        if (o > 1) {
          value = std::min(value, m_new_treshold);
          moc->IncrementCount(o - 1, &v[1], value);
        }
      }
    }
    return;
  }

  // Ugly way of doing things...
  for (int o = this->m_order; o >= 1; o--) {
    moc->StepCountsOrder(true, o, &v[0], &value);
    while (moc->StepCountsOrder(false, o, &v[0], &value)) {
      bool flag = false;
      for (int i = 1; i < o; i++) {
        if (v[i] == this->m_sent_boundary) {
          moc->DeleteCurrentST(o);
          flag = true;
          break;
        }
      }
      if (o > 1) {
        if (v[0] == this->m_sent_boundary)
          moc->IncrementCount(o - 1, &v[0], value);
        moc->IncrementCount(o - 1, &v[1], std::min(value, m_new_treshold));
      }
      if (!flag)
        moc->IncrementBackoffDen(o, &v[0], value);
    }
  }
  remove_zeroprob_grams();

  if (zero_sentstart)
    remove_sent_start_prob();
}

template <typename KT, typename CT>
void InterKn_t<KT, CT>::estimate_bo_counts_absolute_discounting(
    bool zero_sentstart) {
  /* Estimate backoff and lower order counts            */
  std::vector<KT> v(this->m_order);
  CT value;

  if (this->m_sent_boundary < 0) {
    for (int o = this->m_order; o >= 1; o--) {
      moc->StepCountsOrder(true, o, &v[0], &value);
      while (moc->StepCountsOrder(false, o, &v[0], &value)) {
        moc->IncrementBackoffDen(o, &v[0], value);
        if (o > 1)
          moc->IncrementCount(o - 1, &v[1], value);
      }
    }
    return;
  }

  // Ugly way of doing things...
  for (int o = this->m_order; o >= 1; o--) {
    moc->StepCountsOrder(true, o, &v[0], &value);
    while (moc->StepCountsOrder(false, o, &v[0], &value)) {
      bool flag = false;
      for (int i = 1; i < o; i++) {
        if (v[i] == this->m_sent_boundary) {
          moc->DeleteCurrentST(o);
          flag = true;
          break;
        }
      }
      if (o > 1)
        moc->IncrementCount(o - 1, &v[1], value);
      if (!flag)
        moc->IncrementBackoffDen(o, &v[0], value);
    }
  }

  if (zero_sentstart)
    remove_sent_start_prob();
}

template <typename KT, typename CT>
void InterKn_t<KT, CT>::remove_zeroprob_grams() {
  for (int o = moc->m_counts.size() - 1; o >= 2; o--) {
    sikMatrix<KT, CT> *m = moc->m_counts[o];
    for (indextype i = 0; i < m->num_entries(); i++)
      if (*(m->Idx2Valp(i)) < 1e-3)
        RemoveEntryIdx(m->m, i--);
  }
}

// Functions for the n-gram growing algorithm

//// DEBUG FUNCTIONS

template <typename KT, typename ICT>
void Varigram_t<KT, ICT>::printmatrix_bo(
    sikMatrix<KT, typename MultiOrderCounts<KT, ICT>::bo_3c> *m) {
  typename MultiOrderCounts<KT, ICT>::bo_3c *value;
  for (indextype i = 0; i < m->num_entries(); i++) {
    print_indices(m->Idx2Keyp(i), m->dims);
    value = m->Idx2Valp(i);
    fprintf(stderr, "=%d(%d %d %d)\n", value->den, value->nzer[0],
            value->nzer[1], value->nzer[2]);
  }
}

template <typename KT, typename ICT>
void Varigram_t<KT, ICT>::printmatrix_bo(
    sikMatrix<KT, typename MultiOrderCounts<KT, ICT>::bo_c> *m) {
  std::vector<KT> idx(m->dims);
  typename MultiOrderCounts<KT, ICT>::bo_c value;
  m->stepthrough(true, &idx[0], &value);
  while (m->stepthrough(false, &idx[0], &value)) {
    print_indices(idx);
    fprintf(stderr, "=%d(%d)\n", value.den, value.nzer);
  }
}
///////////////

Varigram::Varigram(bool use_3nzer, bool absolute_dc)
    : absolute(absolute_dc), m_use_3nzer(use_3nzer), m_datacost_scale(1),
      m_datacost_scale2(0), m_ngram_prune_target(0), m_max_order(INT_MAX),
      m_small_memory(false) {}

template <typename KT, typename ICT>
Varigram_t<KT, ICT>::Varigram_t(bool use_3nzero, bool absolute_dc)
    : Varigram(use_3nzero, absolute_dc), m_kn(NULL), m_initial_ng(NULL),
      m_data(NULL) {}

template <typename KT, typename ICT> Varigram_t<KT, ICT>::~Varigram_t() {
  delete m_kn;
  delete m_initial_ng;
  delete m_data;
}

template <typename KT, typename ICT>
void Varigram_t<KT, ICT>::initialize(std::string infilename, indextype hashsize,
                                     int ndrop, int nfirst,
                                     std::string optiname, std::string clhist,
                                     bool smallmem, std::string vocabname) {
  io::Stream in;
  io::Stream::verbose = true;
  m_data = new Storage_t<KT, ICT>;
  Storage_t<KT, ICT> *datastorage_tmp;
  m_small_memory = smallmem;
  if (m_small_memory)
    datastorage_tmp = NULL;
  else
    datastorage_tmp = m_data;

  fprintf(stderr, "Creating unigram model, ");
  if (hashsize)
    fprintf(stderr, "hashize %d, ", hashsize);
  /* Construct the unigram model, choose ther right template arguments*/
  if (!m_use_3nzer) {
    fprintf(stderr, "using ...\n");
    m_kn = new InterKn_int_disc<KT, ICT>(
        absolute, infilename, vocabname, optiname, false, 1, ndrop, nfirst,
        datastorage_tmp, infilename, clhist, hashsize);
  } else {
    fprintf(stderr, "using modified ...\n");
    m_kn = new InterKn_int_disc3<KT, ICT>(
        absolute, infilename, vocabname, optiname, false, 1, ndrop, nfirst,
        datastorage_tmp, infilename, clhist, hashsize);
  }

  // Setting discounting
  m_kn->init_disc(0.71);
  m_vocab = &(m_kn->vocab);
  // m_kn->use_ehist_pruning(m_data->size());
  fprintf(stderr, "done\n");

  m_infilename = infilename;

  if (clhist.size())
    set_clear_symbol(clhist);
}

template <typename KT, typename ICT>
void Varigram_t<KT, ICT>::grow(int iter2_lim) {
  std::vector<KT> new_history;
  std::vector<int> accepted_mem;
  int tot_acc = 0;
  int iter = 0, iter2 = 0;
  int cur_order;
  int accepted;

  std::vector<sikMatrix<KT, ICT> *> *sik_c;
  sik_c = &(m_kn->moc->m_counts);

  // int update_coeff_counter=1;
  while (iter2 < iter2_lim) {
    int old_hist_size = -1;
    while (true) {
      accepted = 0;
      iter++;
      while (true) {
        if (!m_kn->MocNextVector(new_history) ||
            new_history.size() >= m_max_order)
          goto LAST;

        if (new_history.size() != old_hist_size) {
          if (old_hist_size > 0) {
            m_kn->set_leaveoneout_discounts(old_hist_size);
          }
          // update_coeff_counter++;
          cur_order = new_history.size();
          sikMatrix<KT, ICT> *curref1 = (*sik_c)[new_history.size()],
                             *curref2 = NULL;
          if (iter2 > 0)
            curref2 = (*sik_c)[new_history.size()];
          if (!m_small_memory)
            m_data->initialize_fast_search_lists(new_history.size() + 1,
                                                 curref1, curref2);
          else
            m_data->init_fsl_file(new_history.size() + 1, curref1, m_infilename,
                                  m_vocab);
        }
        old_hist_size = new_history.size();

        if (reestimate_with_history(new_history)) {
          accepted++;
          // if (update_coeff_counter) update_coeff_counter++;
          // if (!(accepted%1000)) {
          //  fprintf(stderr,"New:");
          //  for (size_t j=0;j<new_history.size();j++)
          //    fprintf(stderr," %s",m_kn->vocab.word(new_history[j]).c_str());
          //  fprintf(stderr,"\n");
          //}
        } // else fprintf(stderr,".");

#if 0
	if (update_coeff_counter==10000) {
	  update_coeff_counter=0;
	  /* Re-estimate smoothing, set lower limit of smoothing 
	     with a little data */
	  m_kn->find_coeffs(0.1, 0.7, 1.3); //0.3, 0.5, 1.0
	  if (!m_use_3nzer) {
	    std::vector<float> d;
	    m_kn_ic->disc2flatv(d);
	    for (int i=0;i<d.size();i++) {
	      if (m_kn_ic->MocOrderSize(i+1)<10000 && d[i]<0.20) d[i]=0.20;
	    }
	    m_kn_ic->flatv2disc(d);
	  } else {
	    std::vector<float> d;
	    m_kn_i3c->disc2flatv(d);
	    for (int i=0;i<d.size();i++) 
	      if (m_kn_i3c->MocOrderSize(i+1)<10000 && d[i]<0.20) d[i]=0.20;
	    m_kn_i3c->flatv2disc(d);
	  }
	}
#endif
      }
      accepted_mem.push_back(accepted);
      tot_acc += accepted;
      for (int i = 0; i < accepted_mem.size(); i++)
        fprintf(stderr, "Round %d: accepted %d\n", i, accepted_mem[i]);
    }
  LAST:
    // for (int i=1;i<=m_kn->order();i++) {
    //  fprintf(stderr,"matrix order %d:\n",i);
    //  m_kn->print_matrix(i);
    //}
    if (old_hist_size > 0) {
      m_kn->set_leaveoneout_discounts(old_hist_size);
    }
    m_kn->find_coeffs(0.007, 1e-1, 5e-2);

    // Ugly fix follows, make this look nicer....
    if (m_kn->get_sentence_boundary_symbol() > 0) {
      m_kn->remove_sent_start_prob();
    }

    prune();
    tot_acc += accepted;
    fprintf(stderr, "%d iterations, %d accepted\n", iter, tot_acc);
    iter2++;
  }
}

template <typename KT, typename ICT> void Varigram_t<KT, ICT>::prune() {
  if (m_ngram_prune_target) {
    double cur_scale = m_datacost_scale;

    int round = 0;
    indextype prev_num_grams = m_kn->num_grams();
    double prev_scale = cur_scale * 2;

    while (double(m_kn->num_grams()) > double(m_ngram_prune_target) * 1.03) {

      if (round == 0) {
        fprintf(stderr, "Currently %d ngrams. First prune with E=D=%.5f\n",
                m_kn->num_grams(), cur_scale);
        m_kn->prune_model(cur_scale, 1, m_small_memory ? NULL : m_data);
        ++round;
        continue;
      }

      double scale_diff = cur_scale - prev_scale;
      indextype gram_diff = prev_num_grams - m_kn->num_grams();

      fprintf(stderr,
              "Previous round increased E from %.4f to %.4f and this pruned "
              "the model from %d to %d ngrams\n",
              prev_scale, cur_scale, prev_num_grams, m_kn->num_grams());
      fprintf(stderr, "I still need to remove %d grams\n",
              m_kn->num_grams() - m_ngram_prune_target);

      double increase = (double)(m_kn->num_grams() - m_ngram_prune_target) /
                        (double)(gram_diff);
      fprintf(stderr,
              "Without limits I would increase E with %.4f (which is %.4f %%) "
              "to %.4f\n",
              increase * scale_diff, (increase * scale_diff) / cur_scale,
              cur_scale + (increase * scale_diff));

      prev_scale = cur_scale;
      prev_num_grams = m_kn->num_grams();

      cur_scale = std::max(
          std::min(cur_scale + (increase * scale_diff), cur_scale * 1.5),
          cur_scale * 1.05);

      fprintf(stderr,
              "With limits I increase E with %.4f (which is %.4f %%) to %.4f\n",
              cur_scale - prev_scale, (cur_scale - prev_scale) / prev_scale,
              cur_scale);

      m_kn->prune_model(cur_scale, 1, m_small_memory ? NULL : m_data);
    }

    fprintf(stderr, "Finally, %d grams, which is %.4f %% off target\n",
            m_kn->num_grams(),
            (double)(m_ngram_prune_target - m_kn->num_grams()) /
                (double)(m_ngram_prune_target));
    if (double(m_kn->num_grams()) < double(m_ngram_prune_target) * 0.97) {
      fprintf(stderr,
              "WARNING: we pruned a bit too much! Increase D and run model "
              "training again to get the desired amount of n-grams\n");
    }
  } else {
    m_kn->prune_model(m_datacost_scale2, 1, m_small_memory ? NULL : m_data);
  }
  m_kn->find_coeffs(0.007, 8e-2, 5e-2);
}

template <typename KT, typename ICT>
bool Varigram_t<KT, ICT>::reestimate_with_history(std::vector<KT> &v) {
  /* Modify the counts for the new history */
  bool accepted;
  int idx;
  ICT val;
  std::map<KT, ICT> m_new_c;

#if 0
  fprintf(stderr,"WHIST: ");print_indices(stderr, v); //fprintf(stderr,"\n");
  fprintf(stderr," [");
  for (int i=0;i<v.size();i++) {
    fprintf(stderr," %s", m_vocab->word(v[i]).c_str());
  }
  fprintf(stderr,"]\n");
#endif

  /* Collect relevant statistics */

  // FIXME: This is inefficient, the map has already been built in Storage.hh
  ICT new_c_sum = 0;
  m_data->fast_search_next(&v, &idx, &val);
  m_data->fast_search_next(NULL, &idx, &val);
  while (idx >= 0) {
    m_new_c[idx] += val;
    new_c_sum += val;
    m_data->fast_search_next(NULL, &idx, &val);
  }

  if (m_new_c.size() == 0) {
    // fprintf(stderr,"No hits found\n");
    return (false);
  }

  m_kn->MocResetCaches();
  // m_kn->print_matrix(v.size());
  double delta1 = modify_model(m_new_c, v, 1.0 / new_c_sum);
  double delta2 = m_new_c.size() * m_kn->model_cost_scale;

  const long size = m_kn->num_grams();
  const long origsize = size - m_new_c.size();
  delta2 += size * log2(size) - origsize * log2(origsize);
  delta2 *= m_datacost_scale;

  double delta = delta1 + delta2;

  // fprintf(stderr,"sizes %ld %ld\n", (long) m_new_c->size(), (long) origsize);
  // fprintf(stderr,"scales %g %g\n", m_datacost_scale, m_kn->model_cost_scale);
  // fprintf(stderr,"delta %.1f+%.1f=%.1f",delta1,delta2,delta);

  // fprintf(stderr,"hiorder mod: ");
  // m_kn->print_matrix(v.size()+1);
  if (delta < 0 /*&& delta1 < -10*m_datacost_scale*/) {
    // fprintf(stderr,"\t Accepted.\n");
    accepted = true;
  } else {
    // fprintf(stderr,"\tRejected.\n");
    m_kn->MocUndoCached();
    accepted = false;
  }
  // fprintf(stderr,"hiorder : ");
  // m_kn->print_matrix(v.size()+1);
  return accepted;
}

template <typename KT, typename ICT>
double Varigram_t<KT, ICT>::modify_model(std::map<KT, ICT> &m,
                                         const std::vector<KT> &indices,
                                         const float ml_norm) {
  ICT val;
  int order = indices.size() + 1;
  double logprobdelta = 0.0;

  /* Check, if m_kn->m_order should be increased */
  if (m_kn->order() < order)
    m_kn->set_order(order);

  std::vector<KT> ind(indices);
  ind.resize(ind.size() + 1);

  /* Get the current coding cost */

  // fprintf(stderr,"MAPPED %d\n", m.size());
  float ml_safelogprob = 0.0;

  typename std::map<KT, ICT>::iterator it = m.begin();
  for (; it != m.end(); it++) {
    ind.back() = it->first;
    // fprintf(stderr,"Probbing ");print_indices(ind);fprintf(stderr,"\n");
    logprobdelta += safelogprob(m_kn->tableprob(ind)) * it->second;
    // ml safelogprob enables earlier pruning and makes things slightly
    // faster
    ml_safelogprob += safelogprob(it->second * ml_norm) * it->second;
  }

  const long origsize = m_kn->num_grams();
  const long size = origsize + m.size();
  // fprintf(stderr,"lpdelta %g + ml_safelogpro %g, +dscale1 %g *( mc %g +
  // sdelta %zd + s1 %g-s2 %g\n", logprobdelta, -ml_safelogprob,
  // m_datacost_scale, m_kn->model_cost_scale, m.size(), size*log2(size),
  // origsize*log2(origsize));
  if ((logprobdelta - ml_safelogprob) +
          m_datacost_scale * (m_kn->model_cost_scale * m.size() +
                              size * log2(size) - origsize * log2(origsize)) >=
      0)
    return 1;

  if (!m_use_3nzer) {
    for (it = m.begin(); it != m.end(); it++) {
      MultiOrderCounts<KT, ICT> *moc = m_kn->moc;
      ind.back() = it->first;
      val = it->second;
#if 0
      fprintf(stderr,"Counts %d\n", order);
      m_kn->print_matrix(order); 
      if (order>=2) {
      	fprintf(stderr,"Counts %d\n", order-1);
      	 m_kn->print_matrix(order-1);
      }
      //if (moc->m_backoffs.size()>order) {
      //fprintf(stderr,"backoffs %d\n", order);
//	printmatrix_bo(moc->m_backoffs[order]);
      //     }
      /* Increment the new higher order counts */
      fprintf(stderr,"Increment ");
      print_indices(ind);
      fprintf(stderr,"= %ld\n", (long) val);
#endif
      /* Increment the new higher order counts */
      ICT debug = moc->IncrementCountCache(order, &ind[0], val);
      assert(debug >= 0);
      moc->IncrementBackoffCacheDen(order, &ind[0], val);

      /* The nonzero counts for new order */
      moc->IncrementBackoffCacheNzer(order, &ind[0], 1);

      if (!absolute) {
        if (val > 1) {
          /* In KN smoothing, lower orders are only incremented for
             unseen contexts */
          ICT ic = moc->IncrementCountCache(order - 1, &ind[1], -val + 1);
          if (ic == -val + 1) {
            moc->IncrementCountCache(order - 1, &ind[1], -ic);
            continue;
          } else
            assert(ic > 0);
          moc->IncrementBackoffCacheDen(order - 1, &ind[1], -val + 1);
        }
      }
    }
  } else {
    MultiOrderCounts<KT, ICT> *moc = m_kn->moc;
    for (it = m.begin(); it != m.end(); it++) {
      ind.back() = it->first;
      val = it->second;
#if 0
      fprintf(stderr,"Counts %d\n", order);
      m_kn->print_matrix(order); 
      if (order>=2) {
      	fprintf(stderr,"Counts %d\n", order-1);
      	 m_kn->print_matrix(order-1);
      }
      if (moc->m_backoffs.size()>order) {
	fprintf(stderr,"backoffs %d\n", order);
	printmatrix_bo(moc->m_backoffs[order]);
      }

      typename MultiOrderCounts<KT, ICT>::bo_3c bo_tmp;
      moc->GetBackoff(1, NULL, &bo_tmp);
      fprintf(stderr,"bo 1 : %d\n", bo_tmp.den);
      fprintf(stderr,"nz 1,1 : %d\n", bo_tmp.nzer[0]);
      fprintf(stderr,"nz 1,2 : %d\n", bo_tmp.nzer[1]);
      fprintf(stderr,"nz 1,3+ : %d\n", bo_tmp.nzer[2]);
      
      /* Increment the new higher order counts */
      fprintf(stderr,"Increment ");
      print_indices(ind);
      fprintf(stderr,"= %d\n",val);
#endif
      moc->IncrementCountCache(order, &ind[0], val);
      moc->IncrementBackoffCacheDen(order, &ind[0], val);

      /* The nonzero counts for new order */
      moc->IncrementBackoffCacheNzer(order, &ind[0], std::min(val - 1, (ICT)2),
                                     1);
      // DEBUGGGG
      // printmatrix_bo(moc->m_backoffs[order]);

      if (!absolute) {
        /* In KN smoothing, lower orders are only incremented for
           unseen contexts */
        if (val > 1) {
          const ICT i = moc->IncrementCountCache(order - 1, &ind[1], -val + 1);
          if (i == -val + 1) {
            moc->IncrementCountCache(order - 1, &ind[1], -i);
            continue;
          } else
            assert(i > 0);
          moc->IncrementBackoffCacheDen(order - 1, &ind[1], -val + 1);
          // fprintf(stderr,"decrement BO %d\n",-val+1);

          // int orig_count = moc->GetCount(order-1, &ind[1]);

          if (i > 2)
            continue;
          /* Did the nz class change ? */
          ICT nzer[3];
          nzer[0] = 0;
          nzer[1] = 0;
          nzer[2] = 0;
          if (i == 1) {
            nzer[0] = 1;
            if (val == 2)
              nzer[1] = -1;
            else
              nzer[2] = -1;
            moc->IncrementBackoffCacheNzer(order - 1, &ind[1], nzer);
            continue;
          }
          nzer[1] = 1;
          nzer[2] = -1;
          moc->IncrementBackoffCacheNzer(order - 1, &ind[1], nzer);
        }
      }
    }
  }
  /* Get the current coding cost */
  for (it = m.begin(); it != m.end(); it++) {
    ind.back() = it->first;
    logprobdelta -= safelogprob(m_kn->tableprob(ind)) * it->second;
  }
  return (logprobdelta);
}

template <typename KT, typename ICT>
void Varigram_t<KT, ICT>::write(FILE *out, bool arpa) {
  TreeGram t;
  m_kn->counts2lm(&t);

  if (!arpa) {
    t.write(out, true);
    return;
  }
  t.write(out);
}

template <typename KT, typename ICT>
void Varigram_t<KT, ICT>::get_unigram_counts(std::string &infilename, int ndrop,
                                             int nfirst, int *type) {
  io::Stream in(infilename, "r", io::REOPENABLE);
  m_initial_ng->count(in.file, true);
  if (nfirst <= 0)
    nfirst = 10000000;
  m_initial_ng->shrink(ndrop, nfirst);
}

template <typename KT, typename ICT>
void Varigram_t<KT, ICT>::get_unigram_counts(std::string &infilename, int ndrop,
                                             int nfirst, unsigned short *type) {
  io::Stream in;
  in.open(infilename, "r", io::REOPENABLE);
  NgramCounts_t<int, ICT> tmp_counts(1, 0, 500000);
  tmp_counts.count(in.file, true);
  if (nfirst <= 0 || nfirst > 64000)
    nfirst = 64000;
  tmp_counts.shrink(ndrop, nfirst);
  in.close();
  io::Stream::verbose = false;

  tmp_counts.vocab->copy_vocab_to(*(m_initial_ng->vocab));
  in.open(infilename, "r", io::REOPENABLE);
  m_initial_ng->count(in.file, false);
  in.close();
}

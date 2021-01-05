// Easily modifiable presentation of n-gram language model
#include "HashGram.hh"
#include "def.hh"
#include "str.hh"

template <typename KT> HashGram_t<KT>::~HashGram_t() {
  for (size_t i = 1; i < probs.size(); i++)
    delete probs[i];
  for (size_t i = 1; i < backoffs.size(); i++)
    delete backoffs[i];
}

template <typename KT> void HashGram_t<KT>::read_real(FILE *file) {
  ArpaReader areader(this);
  bool interpolated;
  std::string line;

  areader.read_header(file, interpolated, line);
  if (interpolated) {
    m_type = INTERPOLATED;
  }

  m_order = areader.counts.size();
  probs.resize(m_order + 1);
  backoffs.resize(m_order + 1);

  std::vector<int> tmp_gram(1);
  float log_prob, back_off;
  int prev_order = 0;
  while (areader.next_gram(file, line, tmp_gram, log_prob, back_off)) {
    int order = tmp_gram.size();
    if (order > prev_order) {
      probs[order] = new sikMatrix<KT, float>(order, areader.counts[order - 1],
                                              MINLOGPROB);
      backoffs[order] =
          new sikMatrix<KT, float>(order, areader.counts[order - 1], 0.0);
      prev_order = order;
    }
    // Inefficiency due to abstracting ArpaReader out
    std::vector<KT> gram(tmp_gram.begin(), tmp_gram.end());
    if (log_prob > MINLOGPROB) {
      probs[order]->setvalue(&gram[0], log_prob);
    }
    if (back_off < 0) {
      backoffs[order]->setvalue(&gram[0], back_off);
    }
  }
}

template <typename KT>
void HashGram_t<KT>::write_real(FILE *out, std::string field_separator) {
  std::vector<std::string> strbuf;
  strbuf.reserve(100000);
  std::vector<int> num_grams(m_order + 1, 0);

  for (int o = 1; o < m_order; o++) {
    std::vector<KT> gram(o);
    float logprob, bo;
    backoffs[o]->stepthrough(true, &gram[0], &logprob);
    while (backoffs[o]->stepthrough(false, &gram[0], &bo)) {
      if (bo >= 0.0)
        continue;
      FindEntry(probs[o]->m, (byte *)(&gram[0]), 1);
    }
  }

  for (int o = 1; o <= m_order; o++) {
    std::vector<KT> gram(o);
    float logprob, bo;
    probs[o]->stepthrough(true, &gram[0], &logprob);
    while (probs[o]->stepthrough(false, &gram[0], &logprob)) {
      bo = backoffs[o]->getvalue(&gram[0]);
      if (m_print_zerograms || logprob > MINLOGPROB || bo < 0)
        num_grams[o]++;
    }
  }

  if (m_type == INTERPOLATED)
    fprintf(out, "\\interpolated\n");
  fprintf(out, "\\data\\\n");

  for (int i = 1; i <= m_order; i++) {
    fprintf(out, "ngram %d=%d\n", i, num_grams[i]);
  }

  for (int o = 1; o <= m_order; o++) {
    fprintf(out, "\n\\%d-grams:\n", o);
    std::vector<KT> gram(o);
    float logprob, bo;
    probs[o]->stepthrough(true, &gram[0], &logprob);
    while (probs[o]->stepthrough(false, &gram[0], &logprob)) {
      bo = backoffs[o]->getvalue(&gram[0]);
      if (!m_print_zerograms && (logprob <= MINLOGPROB && bo >= 0))
        continue;
      fprintf(out, "%.4f", logprob);
      fprintf(out, "%s%s", field_separator.c_str(), word(gram[0]).c_str());
      for (int i = 1; i < o; i++)
        fprintf(out, " %s", word(gram[i]).c_str());
      if (bo < 0)
        fprintf(out, "%s%.4f", field_separator.c_str(), bo);
      fprintf(out, "\n");
    }
  }
  fprintf(out, "\\end\\\n");
}

template <typename KT>
float HashGram_t<KT>::log_prob_bo_helper(const std::vector<KT> &gram) {
  const int looptill = std::min(gram.size(), (size_t)m_order);
  int n = looptill;
  m_last_order = 1;
  float log_prob = 0.0;
  const KT *const gram_ptr = &gram[gram.size() - looptill];

  while (true) {
    // fprintf(stderr, "looptill %d, gramidx %d, ", looptill,looptill-n);
    // print_indices(&gram_ptr[looptill-n],n);
    const float prob = probs[n]->getvalue(&gram_ptr[looptill - n]);

    // Full gram found?
    if (prob > MINLOGPROB) {
      // fprintf(stderr,"full gram n=%d, prob %.4f(%.4f), limit
      // %.4f\n",n,prob,log_prob+prob,MINLOGPROB);
      log_prob += prob;
      m_last_order = n;
      break;
    } else if (n == 1) {
      log_prob += MINLOGPROB;
    }

    if (n == 1)
      break;
    log_prob += backoffs[n - 1]->getvalue(&gram_ptr[looptill - n]);
    n--;
  }
  return log_prob;
}

template <typename KT>
float HashGram_t<KT>::log_prob_i_helper(const std::vector<KT> &gram) {
  float prob = 0.0;
  m_last_order = 0;

  const int looptill = std::min(gram.size(), (size_t)m_order);
  for (int n = 1; n <= looptill; n++) {
    if (n > 1)
      prob *= pow(10, backoffs[n - 1]->getvalue(&gram[gram.size() - n]));

    float p = probs[n]->getvalue(&gram[gram.size() - n]);
    if (p <= MINLOGPROB && n > 1) {
      continue;
    }
    m_last_order = n;
    prob += pow(10, p);
  }
  return ((float)safelogprob(prob));
}

// Overloading stuff
template <>
inline float
HashGram_t<unsigned short>::log_prob_bo(const std::vector<int> &gram) {
  std::vector<unsigned short> g(gram.size());
  for (size_t i = 0; i < gram.size(); i++)
    g[i] = gram[i];
  return (log_prob_bo_helper(g));
}

template <>
inline float
HashGram_t<unsigned short>::log_prob_i(const std::vector<int> &gram) {
  std::vector<unsigned short> g(gram.size());
  for (size_t i = 0; i < gram.size(); i++)
    g[i] = gram[i];
  return (log_prob_i_helper(g));
}

template <>
inline float HashGram_t<int>::log_prob_bo(const std::vector<int> &gram) {
  return log_prob_bo_helper(gram);
}

template <>
inline float HashGram_t<int>::log_prob_i(const std::vector<int> &gram) {
  return log_prob_i_helper(gram);
}

template <>
inline float HashGram_t<unsigned short>::log_prob_bo(
    const std::vector<unsigned short> &gram) {
  return log_prob_bo_helper(gram);
}

template <>
inline float HashGram_t<unsigned short>::log_prob_i(
    const std::vector<unsigned short> &gram) {
  return log_prob_i_helper(gram);
}

template <>
inline float
HashGram_t<int>::log_prob_bo(const std::vector<unsigned short> &gram) {
  assert(false);
  return 0;
}

template <>
inline float
HashGram_t<int>::log_prob_i(const std::vector<unsigned short> &gram) {
  assert(false);
  return 0;
}

template <typename KT> void HashGram_t<KT>::prune(float treshold) {
  // As in Stolcke's paper
  assert(m_type == BACKOFF);
  bool verbose = 0;
  std::vector<std::vector<float>> hist_cumuprobs;

  for (int o = m_order; o > 1; o--) { // Does not prune unigrams...
    if (verbose)
      fprintf(stderr, "Processing order %d\n", o);
    std::vector<KT> gram(o);
    ///////////////////////////////////////////////////////////////
    // 1. Cache the total prob mass given to backoff =
    //                               1-total unbackoffed mass
    twofloat f2(0.0, 0.0);
    sikMatrix<KT, twofloat> bo_pcache(o - 1, probs[o - 1]->num_entries(), f2);
    {
      float gram_logprob;
      std::vector<KT> g2(o - 1);
      probs[o]->stepthrough(true, &gram[0], &gram_logprob);
      while (probs[o]->stepthrough(false, &gram[0], &gram_logprob)) {
        f2.float1 = pow(10, gram_logprob);
        for (int i = 1; i < o; i++)
          g2[i - 1] = gram[i];
        f2.float2 = pow(10, log_prob(g2));
        bo_pcache.increment(&gram[0], f2);
      }
    }
    if (verbose)
      fprintf(stderr, "Cached backoff mass.\n");

    float gram_logprob;
    probs[o]->stepthrough(true, &gram[0], &gram_logprob);
    while (probs[o]->stepthrough(false, &gram[0], &gram_logprob)) {
      float trunc_logprob;
      {
        std::vector<KT> g2(o - 1);
        for (int i = 1; i < o; i++)
          g2[i - 1] = gram[i];
        trunc_logprob = log_prob(g2);
        if (verbose) {
          fprintf(stderr, "processing ");
          print_indices(stderr, gram);
          fprintf(stderr, ", gram lp %.4f, trunc %.4f\n", gram_logprob,
                  trunc_logprob);
        }
      }
      ///////////////////////////////////////////////////////////
      // 2. Mariginal history probabilities p(h) = p(h1)p(h2|h1)...
      float hist_log_prob = 0.0;
      // This is the correct way
      std::vector<KT> g2(gram);
      for (int o2 = o - 1; o2 >= 1; o2--) {
        g2.resize(o2);
        hist_log_prob += log_prob(g2);
      }

      if (verbose)
        fprintf(stderr, "histprob %.4f, ", hist_log_prob);

      ////////////////////////////////////////////////////////
      // 3. The values of the a(h) and a(h')
      float alpha, alpha2;
      float bo_psum;
      float new_bo_psum;
      float new_bo_psum2;
      {
        twofloat f2 = bo_pcache.getvalue(&gram[0]);
        bo_psum = f2.float1;
        const float bo_psum2 = f2.float2;
        if (verbose)
          fprintf(stderr, "bo_psum %.4f, bo_psum2 %.4f\n", bo_psum, bo_psum2);
        alpha = (1.0 - bo_psum) / (1.0 - bo_psum2);
        new_bo_psum = bo_psum - pow(10, gram_logprob);
        new_bo_psum2 = bo_psum2 - pow(10, trunc_logprob);
        alpha2 = (1.0 - new_bo_psum) / (1.0 - new_bo_psum2);
      }
      if (verbose)
        fprintf(stderr, "a=%.4f, a'=%.4f\n", alpha, alpha2);
      //////////////////////////////////////////////////////////////
      // 4. Difference in entropy =
      //      -p(h) { p(w|h) [ log p(w|h') + log a'(h) - log p(w|h)]
      //                + [ log a(h') - log a(h) ] * backoffed_mass
      float diff;
      float log_alpha2;
      {
        const float log_alpha = (float)safelogprob(alpha);
        const float hist_prob = (float)pow(10, hist_log_prob);
        const float gram_prob = (float)pow(10, gram_logprob);
        log_alpha2 = (float)safelogprob(alpha2);

        diff = -hist_prob *
               (gram_prob * (trunc_logprob + log_alpha2 - gram_logprob) +
                (1 - bo_psum) * (log_alpha2 - log_alpha));
        if (verbose) {
          fprintf(stderr, "diff = %.4f -p(h) * ", -hist_prob);
          fprintf(stderr,
                  "( %4f p(w|h) * ( %.4f logp(w|h') + %.4f log a'(h) - %.4f "
                  "log p(w|h)) ",
                  gram_prob, trunc_logprob, log_alpha2, gram_logprob);
          fprintf(stderr, "%.4f bomass * ( %.4flog a(h') - %.4f loga(h)))\n",
                  1 - bo_psum, log_alpha2, log_alpha);
          fprintf(stderr, " = diff %.4g (%.4g)\n", diff, treshold);
        }
      }
      if (diff > treshold)
        continue;

      /////////////////////////////////////////////////////////////
      // 5. Modify the model
      {
        if (verbose)
          fprintf(stderr, "rejecting\n");

        // MINLOGPROB-1 is not the default value MINLOGPROB. This prevents
        // the model from reordering the internal ordering and breaking
        // steptrough(). Really bad interface design...
        probs[o]->setvalue(&gram[0], MINLOGPROB - 1);

        if (verbose)
          fprintf(stderr, "Setting backoff %d to %.4f\n", gram[0], log_alpha2);
        backoffs[o - 1]->setvalue(&gram[0], log_alpha2);
        twofloat f2(new_bo_psum, new_bo_psum2);
        bo_pcache.setvalue(&gram[0], f2);
      }
    }
  }
  remove_empty_grams();
}

template <typename KT> void HashGram_t<KT>::remove_empty_grams() {
  for (int o = m_order; o > 1; o--) {
    std::vector<KT> gram(o);
    for (indextype i = probs[o]->num_entries() - 1; i >= 0; i--) {
      // fprintf(stderr,"Checking ");
      // print_indices(probs[o]->Idx2Keyp(i), o);
      // fprintf(stderr," =%.4f ",*(probs[o]->Idx2Valp(i)));
      // fprintf(stderr," (%.4f)\n",
      // backoffs[o]->getvalue(probs[o]->Idx2Keyp(i)));
      if (*(probs[o]->Idx2Valp(i)) <= MINLOGPROB) {
        probs[o]->setvalue(probs[o]->Idx2Keyp(i), MINLOGPROB);
      }
    }
  }
}

template <typename KT> void HashGram_t<KT>::add_zeroprob_grams() {
  for (int o = m_order; o >= 2; o--) {
    std::vector<KT> gram(o);
    float gram_logprob;
    const float inc = 0.0;
    m_print_zerograms = true;

    probs[o]->stepthrough(true, &gram[0], &gram_logprob);
    while (probs[o]->stepthrough(false, &gram[0], &gram_logprob)) {
      // Don't delete even if default value
      probs[o - 1]->increment_wo_del(&gram[0], inc);
    }
  }
}

template <typename KT>
void HashGram_t<KT>::normalize_and_set_bo(std::vector<KT> &prefix,
                                          float explicit_probsum,
                                          float full_bo_probsum) {
  float new_bo =
      safelogprob((1.0f - explicit_probsum) / (1.0f - full_bo_probsum));
  backoffs[prefix.size()]->setvalue(&prefix[0], new_bo);
  // fprintf(stderr, "Bo for prefix");
  // for (auto i: prefix) {
  //         fprintf(stderr, " %d", i);
  // }
  // fprintf(stderr, " %f\n", new_bo);
}

template <typename KT> void HashGram_t<KT>::renormalize_backoffs(int order) {
  std::vector<KT> gram(order);
  std::vector<KT> cur_prefix(order - 1);
  for (int i = 0; i < order - 1; i++) {
    cur_prefix[i] = -1;
  }

  float gram_logprob;
  float explicit_probsum = 0.0f, full_bo_probsum = 0.0f;
  probs[order]->ordered_stepthrough(true, &gram[0], &gram_logprob);
  while (probs[order]->ordered_stepthrough(false, &gram[0], &gram_logprob)) {
    bool different = false;
    if (cur_prefix.size() && cur_prefix[0] == -1) {
      // first pass, uninitialized, so init here
      for (int i = 0; i < order - 1; i++) {
        cur_prefix[i] = gram[i];
      }
    }
    for (int i = 0; i < order - 1; i++) {
      if (cur_prefix[i] != gram[i]) {
        different = true;
        break;
      }
    }

    if (different) {
      // We have accumulated all prob mass for this prefix,
      // we can now re-estimate the backoff coeff
      normalize_and_set_bo(cur_prefix, explicit_probsum, full_bo_probsum);

      // reset the prefix comparison
      for (int i = 0; i < order - 1; i++) {
        cur_prefix[i] = gram[i];
      }
      explicit_probsum = 0.0f;
      full_bo_probsum = 0.0f;
    }
    std::vector<KT> current_postfix(gram.begin() + 1, gram.end());
    explicit_probsum += pow(10, gram_logprob);
    if (order > 1) {
      full_bo_probsum += pow(10, log_prob_bo(current_postfix));
    }
    // fprintf(stderr, "Probsums for prefix");
    // for (auto i: cur_prefix) {
    //         fprintf(stderr, " %d", i);
    // }
    // fprintf(stderr, " = %.8f %.8f\n", explicit_probsum, full_bo_probsum);
  }
  // Remember to sum the last prefix too !
  if (order > 1) {
    normalize_and_set_bo(cur_prefix, explicit_probsum, full_bo_probsum);
  } else if (explicit_probsum > 0.001 && explicit_probsum < 1.0) {
    // Redistribute the remaining unigram mass equally accross unigrams
    const float uniprob = (1.0f - explicit_probsum) / probs[1]->num_entries();
    probs[1]->stepthrough(true, &gram[0], &gram_logprob);
    while (probs[1]->stepthrough(false, &gram[0], &gram_logprob)) {
      probs[1]->setvalue(&gram[0],
                         safelogprob(pow(10, gram_logprob) + uniprob));
    }
  } else {
    fprintf(stderr,
            "Cannot smooth unigrams, explicit_probsum %f not between 0.001 and 1.0.\n",
            explicit_probsum);
  }
}

template <typename KT>
void HashGram_t<KT>::fake_interpolate(HashGram_t<KT> &other, float lambda) {
  if (m_type != BACKOFF || other.m_type != BACKOFF) {
    fprintf(stderr, "Cannot interpolate non-backoff arpa models. Use arpa2arpa "
                    "to convert the models first!\n");
    exit(-1);
  }

  // Assumes identical vocabs for this and other !
  auto max_order = std::max(m_order, other.order());
  // fprintf(stderr, "Maxorder %d\n", max_order);
  probs.resize(max_order + 1);
  backoffs.resize(max_order + 1);
  while (true) {
    if (m_order == max_order)
      break;
    m_order += 1;
    probs[m_order] = new sikMatrix<KT, float>(
        m_order, other.probs[m_order]->num_entries(), MINLOGPROB);
    backoffs[m_order] = new sikMatrix<KT, float>(
        m_order, other.backoffs[m_order]->num_entries(), 0.0);
  }

  for (int o = max_order; o >= 1; o--) {
    std::vector<KT> gram(o);
    float gram_logprob;
    probs[o]->stepthrough(true, &gram[0], &gram_logprob);
    while (probs[o]->stepthrough(false, &gram[0], &gram_logprob)) {
      probs[o]->setvalue(
          &gram[0],
          safelogprob(lambda * pow(10, gram_logprob) +
                      (1 - lambda) * pow(10, other.log_prob_bo(gram))));
    }

    if (o <= other.order()) {
      float gram_logprob;
      other.probs[o]->stepthrough(true, &gram[0], &gram_logprob);
      while (other.probs[o]->stepthrough(false, &gram[0], &gram_logprob)) {
        if (probs[o]->getvalue(&gram[0]) > MINLOGPROB) {
          // These ngrams have been handled in the prev loop

          continue;
        }
        probs[o]->setvalue(&gram[0],
                           safelogprob(lambda * pow(10, log_prob_bo(gram)) +
                                       (1 - lambda) * pow(10, gram_logprob)));
      }
    }
  }
  for (int o = 1; o <= max_order; o++) {
    renormalize_backoffs(o);
  }
}

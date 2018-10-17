// Functions for calculating perplexity
#include "PerplexityFuncs.hh"
#include "HashGram.hh"
#include "TreeGram.hh"
#include "def.hh"
#include <algorithm>

Perplexity::Perplexity(const char *lm_name) {
  init_variables();
  m_lm.reset(new TreeGram());
  m_lm->set_oov("<UNK>");
  m_tmpstring_length = 100;
  m_tmpstring = (char *)malloc(m_tmpstring_length * sizeof(char));
  m_print_unk_warn = false;
  m_skip_unk_prob = true;

  // try {  Commented out: Let the exceptions rise to surface...
  io::Stream in(lm_name, "rb");
  m_lm->read(in.file, 1);
  //}  catch (std::exception &e) {
  // std::cerr << e.what() << std::endl;
  // exit(1);
  //}
  init_special_symbols("", "", "");
}

Perplexity::Perplexity(std::shared_ptr<NGram> lm, const std::string ccs_name,
                       const std::string wb_name, const std::string mb_name,
                       const std::string unk_symbol, bool skip_unk_prob) {
  init_variables();
  m_lm = lm;
  m_skip_unk_prob = skip_unk_prob;

  if (unk_symbol.length())
    m_lm->set_oov(unk_symbol);
  else
    m_lm->set_oov("<UNK>");
  init_special_symbols(ccs_name, wb_name, mb_name);
}

/* This is the old type of constructor */
Perplexity::Perplexity(const std::string lm_name, const int type,
                       const std::string ccs_name, const std::string wb_name,
                       const std::string mb_name, const std::string unk_symbol,
                       const int hashgram, bool skip_unk_prob) {
  init_variables();
  if (hashgram)
    if (hashgram > 0)
      m_lm.reset(new HashGram_t<int>());
    else
      m_lm.reset(new HashGram_t<unsigned short>());
  else
    m_lm.reset(new TreeGram());
  m_skip_unk_prob = skip_unk_prob;

  // fprintf(stderr,"UNK SYMBOL: %s\n",unk_symbol);
  if (unk_symbol.length())
    m_lm->set_oov(unk_symbol);
  else
    m_lm->set_oov("<UNK>");

  try {
    io::Stream in(lm_name, "r");
    m_lm->read(in.file, type);
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    exit(1);
  }
  init_special_symbols(ccs_name, wb_name, mb_name);
}

void Perplexity::init_variables() {
  m_use_unk = 0;
  m_use_ccs = 0;
  m_num_unks = 0;
  m_num_tunks = 0;
  m_num_ccs = 0;
  m_num_pwords = 0;
  m_num_ptokens = 0;
  m_unkflag = 0;
  m_print_unk_warn = true;
  m_logprob = 0.0;
  m_perplexity = 0.0;
  m_token_perplexity = 0.0;
  m_lm2.reset();
  m_alpha = 0.5;
  m_init_hist = m_cur_init_hist = 0;
  m_num_sent_ends = 0;
  m_cw_lpsum = 0.0;
  m_word_logprob = 0.0;
}

void Perplexity::init_special_symbols(const std::string ccs_name,
                                      const std::string wb_name,
                                      const std::string mb_name) {
  if (ccs_name.length()) {
    std::cerr << "Reading ccs";
    find_indices(ccs_name, ccs_vector);
    std::cerr << ", found " << ccs_vector.size() << " context cues."
              << std::endl;
  }

  if (wb_name.length()) {
    std::cerr << "Reading wb";
    find_indices(wb_name, wb_vector);
    std::cerr << ", found " << wb_vector.size() << " word break symbols."
              << std::endl;
    m_wb_type = LISTED;
  } else if (mb_name.length()) {
    std::cerr << "Reading mb";
    load_mbs(mb_name, mb_vector);
    std::cerr << ", found " << mb_vector.size() << " morph break expressions."
              << std::endl;
    m_wb_type = MB_LISTED;
  } else
    m_wb_type = EVERYTIME;

  // std::cerr << "Ngram order: " << m_lm->order() << std::endl;

  /* Space for the hit rates, zeroth is for the requesed amount */
  ngram_hits.resize(m_lm->order() + 1, 0);
}

void Perplexity::find_indices(const std::string fname, std::vector<int> &vec) {
  char buf[1000];

  io::Stream in(fname, "r");
  int temp;
  while (true) {
    int i = fscanf(in.file, "%s", buf);
    if (!i || i == EOF)
      break;

    vec.push_back(temp = m_lm->word_index(buf));
    if (!temp) {
      fprintf(stderr, "Warning, %s not in the lm !\n", buf);
    }
  }
}

void Perplexity::load_mbs(const std::string fname,
                          std::vector<std::string> &vec) {
  char buf[1000];
  io::Stream in(fname, "r");
  while (true) {
    int i = fscanf(in.file, "%s", buf);
    if (!i || i == EOF)
      break;
    std::string sbuf(buf);
    vec.push_back(sbuf);
  }
}

float Perplexity::sentence_logprob(const char *sentence_in) {
  float lpsum = 0.0, foo;
  char *sentence = strdup(sentence_in);
  char *cptr = strtok(sentence, " ");

  if (strlen(sentence) > m_tmpstring_length) {
    m_tmpstring_length = strlen(sentence);
    m_tmpstring =
        (char *)realloc(m_tmpstring, m_tmpstring_length * sizeof(char));
  }

  while (cptr) {
    sscanf(cptr, "%s", m_tmpstring);
    lpsum += logprob(m_tmpstring, foo);
    cptr = strtok(NULL, " ");
  }
  free(sentence);
  return -lpsum;
}

float Perplexity::logprob(const char *word, float &cur_word_lp) {
  // fprintf(stderr, "\"%s\":\n", word);
  cur_word_lp = 0.0;
  if (m_cur_init_hist && m_cur_init_hist == m_init_hist)
    history.clear();
  else if (history.size() == m_lm->order())
    history.pop_front();

  int idx = m_lm->word_index(word);
  if (m_cur_init_hist > 0) {
    if (strncmp("<s>", word, 3) && is_wb(idx)) {
      m_num_pwords++;
    }
    m_cur_init_hist--;
    history.push_back(idx);
    m_lm->set_last_order(0);
    return (0.0);
  }

  if (!strncmp("</s>", word, 4)) {
    m_cur_init_hist = m_init_hist;
    m_num_sent_ends++;
  }

  if (!idx) { /* UNK */
    if (m_print_unk_warn)
      fprintf(stderr, "Unknown token %s\n", word);
    m_num_tunks++;
    if (m_wb_type == MB_LISTED && !is_mb(word)) {
      m_num_unks++;
      cur_word_lp = m_cw_lpsum;
      m_cw_lpsum = 0.0;
    } else if (m_wb_type == EVERYTIME) {
      m_num_unks++;
    } else {
      m_unkflag = 1;
    }
    if (m_skip_unk_prob) {
      history.push_back(0);
      m_lm->set_last_order(0);
      return (0.0);
    }
  }

  if (std::find(ccs_vector.begin(), ccs_vector.end(), idx) !=
      ccs_vector.end()) {
    /* Context cue */
    m_num_ccs++;
    history.push_back(idx);
    m_lm->set_last_order(0);
    return (0.0);
  }

  history.push_back(idx);
  float lp = m_lm->log_prob(history);
  if (m_lm2) {
    float lp2 = m_lm2->log_prob(history);
    lp = (float)safelogprob(m_alpha * pow(10, lp) +
                            (1 - m_alpha) * pow(10, lp2));
  }
  ngram_hits[0]++;
  ngram_hits[m_lm->last_order()]++;

  m_num_ptokens++;
  m_cw_lpsum += lp;
  if (is_wb(idx)) {
    if (m_unkflag) {
      m_unkflag = 0;
      m_num_unks++;
      if (!m_skip_unk_prob) {
        cur_word_lp = m_cw_lpsum;
        m_num_pwords++;
      }
    } else {
      cur_word_lp = m_cw_lpsum;
      m_num_pwords++;
    }
    m_cw_lpsum = 0.0;
  }
  return (lp);
}

void Perplexity::print_hitrates(FILE *out) {
  fprintf(out, "\nNgram hit rates:\n");
  for (int i = 1; i <= m_lm->order(); i++)
    fprintf(out, "%d: %.3f\n", i, (100.0 * ngram_hits[i] / ngram_hits[0]));
}

void Perplexity::reset_hitrates() {
  ngram_hits.clear();
  ngram_hits.resize(m_lm->order() + 1, 0);
}

int Perplexity::get_hitorder(int i) {
  if (i >= (int)ngram_hits.size())
    return 0;
  return ngram_hits[i];
}

bool Perplexity::is_wb(int idx) {
  if (m_wb_type == EVERYTIME)
    return (true);

  if (m_wb_type == LISTED &&
      (std::find(wb_vector.begin(), wb_vector.end(), idx) != wb_vector.end()))
    return (true);

  if (m_wb_type == MB_LISTED) {
    const std::string &w = m_lm->word(idx);
    bool mb = is_mb(w);
    if (mb) {
      return (false); // there is morph boundary preceeding/following the token
    } else {
      return (true); // token is followed/preceeded by a word break
    }
  }
  return (false);
}

bool Perplexity::is_mb(std::string w) {
  if (m_wb_type != MB_LISTED)
    return (false);
  bool mb = false;
  int wl = w.length();
  std::vector<std::string>::iterator itr;
  for (itr = mb_vector.begin(); itr != mb_vector.end(); ++itr) {
    int l = (*itr).length();
    if ((*itr)[0] == '^') {
      if ((*itr).substr(1) == w.substr(0, l - 1)) {
        mb = true;
        break;
      }
    } else if ((*itr)[l - 1] == '$') {
      if (wl >= l - 1) {
        if ((*itr).substr(0, l - 1) == w.substr(wl - l + 1)) {
          mb = true;
          break;
        }
      }
    }
  }
  return (mb);
}

double Perplexity::logprob_file(FILE *in, FILE *out, const int interval) {
  char word[MAX_WLEN + 1];
  int read = 0;
  float interval_sum = 0.0;
  float cur_full_word_lp;
  const float log2coeff = 1 / log10(2.0f);

  while (1) {
    int fsc = fscanf(in, MAX_WLEN_FMT_STRING, word);
    if (!fsc || fsc == EOF)
      break;
    double lp = logprob(word, cur_full_word_lp);
    read++;
    if (out) {
      // fprintf(out,"lp %.4f isum %.4f\n", lp, interval_sum);
      if (interval == 1)
        fprintf(out, "%s %g, order %d (cfw_lp %.4f)\n", word, pow(10, lp),
                m_lm->last_order(), cur_full_word_lp);
      // else fprintf(out,"%s *\n",word);
      else {
        interval_sum += (float)lp;
        if (read % interval == 0) {
          fprintf(out, "%.4f\n", interval_sum * log2coeff);
          interval_sum = 0.0;
        }
      }
    }
    m_word_logprob += cur_full_word_lp;
    // fprintf(stdout,"mwl %.4f , cfwl %.4f\n", m_word_logprob,
    // cur_full_word_lp);
    m_logprob += lp;
  }
  m_word_logprob += m_cw_lpsum;
  // fprintf(stdout, "mwl %.4f, mcwlp %.4f\n", m_word_logprob, m_cw_lpsum);
  m_perplexity_wo_sent_ends = m_word_logprob / (m_num_pwords - m_num_sent_ends);
  m_perplexity = m_word_logprob / m_num_pwords;
  // fprintf(stdout,"Setting m_perplexity to %.4f / %d = %.4f\n",
  //	  m_word_logprob, m_num_pwords, m_perplexity);
  m_token_perplexity = m_logprob / m_num_ptokens;
  // fprintf(stdout,"Setting m_token_perplexity to %.4f / %d = %.4f\n",
  //	  m_logprob, m_num_ptokens, m_token_perplexity);

  if (m_wb_type == EVERYTIME) {
    m_num_unks = m_num_tunks;
    m_num_tunks = 0;
  }

  return (m_logprob);
}

double Perplexity::print_results(FILE *out) {
  if (m_skip_unk_prob) {
    fprintf(out, "\nDropped:   ");
    fprintf(out, "%d UNKS, %.2f %%\n", m_num_unks,
            100.0 * m_num_unks / (m_num_pwords + m_num_unks + m_num_ccs));
    if (m_wb_type != EVERYTIME)
      fprintf(out, "           %d TUNKS, %.2f %%\n", m_num_tunks,
              100.0 * m_num_tunks / (m_num_ptokens + m_num_tunks + m_num_ccs));
  } else {
    fprintf(out, "\nFound:     ");
    fprintf(out, "%d UNKS, %.2f %%\n", m_num_unks,
            100.0 * m_num_unks / (m_num_pwords + m_num_ccs));
    if (m_wb_type != EVERYTIME)
      fprintf(out, "           %d TUNKS, %.2f %%\n", m_num_tunks,
              100.0 * m_num_tunks / (m_num_ptokens + m_num_ccs));
  }
  if (ccs_vector.size()) {
    if (m_skip_unk_prob) {
      fprintf(out, "           %d context cues\n", m_num_ccs);
    } else {
      fprintf(out, "Dropped:   %d context cues\n", m_num_ccs);
    }
  }
  fprintf(out, "Processed: %d words\n", m_num_pwords);
  if (m_wb_type != EVERYTIME)
    fprintf(out, "           %d tokens\n", m_num_ptokens);

  if (m_skip_unk_prob) {
    fprintf(out, "Total:     %d words\n",
            m_num_pwords + m_num_unks + m_num_ccs);
    if (m_wb_type != EVERYTIME)
      fprintf(out, "           %d tokens\n",
              m_num_ptokens + m_num_tunks + m_num_ccs);
  } else {
    fprintf(out, "Total:     %d words\n", m_num_pwords + m_num_ccs);
    if (m_wb_type != EVERYTIME)
      fprintf(out, "           %d tokens\n", m_num_ptokens + m_num_ccs);
  }

  fprintf(out, "\nLogprob %.6f\n", m_word_logprob);
  fprintf(out, "Perplexity %.2f (- %dth root) = %.3f bits\n",
          pow(10, -m_perplexity), m_num_pwords, -m_perplexity / log10(2.0));
  fprintf(out,
          "Perplexity (sentence ends not in normalization) %.2f (- %dth root) "
          "= %.3f bits\n",
          pow(10, -m_perplexity_wo_sent_ends), m_num_pwords - m_num_sent_ends,
          -m_perplexity_wo_sent_ends / log10(2.0));

  if (m_num_pwords != m_num_ptokens) {
    fprintf(out, "\nTokenwise logprob %.6f\n", m_logprob);
    fprintf(out, "Tokenwise perplexity %.2f (- %dth root) = %.3f bits\n",
            pow(10, -m_token_perplexity), m_num_ptokens,
            -m_token_perplexity / log10(2.0));
  }
  return -m_perplexity / log10(2.0);
}

void Perplexity::set_interpolation(std::string lm_name) {
  auto tmp_lm = std::make_shared<HashGram_t<int> >();

  // Add all words from m_lm to tmp_lm
  std::vector<int> gram(1);
  for (int i = 0; i < m_lm->num_words(); i++) {
    gram[0] = tmp_lm->add_word(m_lm->word(i));
    // tmp_lm->probs[1]->setvalue(&gram[0], -59);
  }

  io::Stream in(lm_name, "r");
  tmp_lm->read(in.file, false);

  // Add all words from tmp_lm to m_lm
  for (int i = 0; i < tmp_lm->num_words(); i++) {
    m_lm->add_word(tmp_lm->word(i));
  }

  if (m_lm->order() >= tmp_lm->order())
    m_lm2 = tmp_lm;
  else {
    m_lm2 = m_lm;
    m_lm = tmp_lm;
    m_alpha = 1 - m_alpha;
  }

  assert(m_lm->num_words() == m_lm2->num_words());
  for (int i = 0; i < m_lm->num_words(); i++) {
    assert(m_lm->word(i) == m_lm2->word(i));
  }

  // io::Stream fh1("out1.arpa", "w");
  // io::Stream fh2("out2.arpa", "w");
  // m_lm->write(fh1.file);
  // m_lm2->write(fh2.file);
}

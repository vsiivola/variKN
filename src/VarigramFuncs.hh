// Functions for the n-gram growing algorithm
#ifndef VARIGRAMFUNCS_HH
#define VARIGRAMFUNCS_HH

#include "InterKn.hh"
#include <map>

class Varigram {
public:
  inline Varigram(bool use_3nzer, bool abslute = false);
  virtual ~Varigram() {}
  inline void set_datacost_scale(double f) { m_datacost_scale = f; }
  inline void set_datacost_scale2(double f) { m_datacost_scale2 = f; }
  inline void set_ngram_prune_target(indextype i) { m_ngram_prune_target = i; }
  inline void set_max_order(int i) { m_max_order = i; }

  virtual void initialize(std::string infilename, indextype hashsize, int ndrop,
                          int nfirst, std::string optiname, std::string clhist,
                          bool smallmem, std::string vocabname = "") = 0;
  virtual void grow(int iter2_lim = 1) = 0;
  virtual void write_narpa(FILE *out) = 0;
  virtual void write_debug_counts(FILE *out) = 0;
  virtual void write(FILE *out, bool arpa) = 0;
  virtual void set_clear_symbol(std::string s) = 0;
  virtual void set_zeroprobgrams(bool) = 0;
  virtual void set_cutoffs(std::vector<int> v) = 0;
  virtual void set_discard_unks(bool x) = 0;
  virtual void set_all_discounts(float x) = 0;
  bool absolute;
  void write_vocab(FILE *out) { m_vocab->write(out); }

  inline void write_file(std::string lmname, bool arpa) {
    io::Stream out(lmname, "w");
    write(out.file, arpa);
    out.close();
  }

protected:
  bool m_use_3nzer;
  float m_datacost_scale;
  float m_datacost_scale2;
  indextype m_ngram_prune_target;
  int m_max_order;
  std::string m_infilename;
  Vocabulary *m_vocab;
  bool m_small_memory;
};

template <typename KT, typename ICT> class Varigram_t : public Varigram {
public:
  Varigram_t(bool use_3nzero, bool absolute = false);
  ~Varigram_t();
  void initialize(std::string infilename, indextype hashsize, int ndrop,
                  int nfirst, std::string optiname, std::string clhist,
                  bool smallmem, std::string vocabname = "");
  void grow(int iter2_lim = 1);
  inline void write_narpa(FILE *out) { m_kn->counts2asciilm(out); }
  inline void write_debug_counts(FILE *out) { m_kn->counts2ascii(out); }
  inline void set_zeroprobgrams(bool x) { m_kn->zeroprobgrams = x; }
  void write(FILE *out, bool arpa);

  inline void set_clear_symbol(std::string s) {
    assert(m_kn);
    m_data->clear_lm_history = m_vocab->word_index(s);
    if (!m_data->clear_lm_history) {
      fprintf(stderr,
              "No \"<s>\" in history, --clear_history cannot be used. Exit.\n");
      exit(-1);
    }
    m_kn->set_sentence_boundary_symbol(s);
  }
  void set_cutoffs(std::vector<int> v) { m_kn->cutoffs = v; }
  void set_discard_unks(bool x) { m_kn->discard_ngrams_with_unk = x; }

  void set_all_discounts(float x) { m_kn->init_disc(x); }

private:
  InterKn_t<KT, ICT> *m_kn;
  NgramCounts_t<KT, ICT> *m_initial_ng;

  bool reestimate_with_history(std::vector<KT> &history);
  double modify_model(std::map<KT, ICT> &new_c, const std::vector<KT> &v,
                      const float ml_norm);
  Storage_t<KT, ICT> *m_data;
  void get_unigram_counts(std::string &infilename, int ndrop, int nfirst,
                          int *type);
  void get_unigram_counts(std::string &infilename, int ndrop, int nfirst,
                          unsigned short *type);

  void
  printmatrix_bo(sikMatrix<KT, typename MultiOrderCounts<KT, ICT>::bo_3c> *m);
  void
  printmatrix_bo(sikMatrix<KT, typename MultiOrderCounts<KT, ICT>::bo_c> *m);

  void prune();
};

#include "VarigramFuncs_tmpl.hh"
#endif

// The main library for the n-gram model estimation
#ifndef INTERKN_HH
#define INTERKN_HH

#include <algorithm>

#include <assert.h>
#include <math.h>

#include "MultiOrderCounts.hh"
#include "NgramCounts.hh"
#include "QFit.hh"
#include "Storage.hh"
#include "TreeGram.hh"
#include "io.hh"

class InterKn : public QFitEvaluator /*for evaluate() */ {
public:
  inline InterKn(const bool absolute_discounting, const std::string &dataname,
                 const std::string &optiname, const std::string &prunesource)
      : zeroprobgrams(true), prune_with_real_counts(false),
        discard_cutoffs(false), discard_ngrams_with_unk(false), m_order(0),
        m_sent_boundary(-1), m_ehist_estimate(0),
        m_absolute_discounting(absolute_discounting), m_data_name(dataname),
        m_opti_name(optiname), m_prunedata_name(prunesource) {}
  virtual ~InterKn() {}
  virtual void create_model(float prunetreshold = -1.0) = 0;
  virtual void counts2lm(TreeGram *lm) = 0;
  virtual void counts2asciilm(FILE *out) = 0;
  virtual void counts2ascii(FILE *out) = 0; // For debugging

  inline int order() { return m_order; }
  inline void use_ehist_pruning(int x) { m_ehist_estimate = x; }
  virtual void set_order(int o) = 0;

  virtual void estimate_bo_counts(bool zerosent = false) = 0;
  virtual void estimate_nzer_counts() = 0;
  virtual float evaluate(std::vector<float> &discounts) = 0;
  virtual void find_coeffs(float brak, float precision,
                           float lin_precision) = 0;
  virtual void init_disc(float x) = 0;
  virtual void set_leaveoneout_discounts(int order) = 0;
  virtual indextype num_grams() = 0;

  /* Added for templatized func, too much casting otherwise */
  virtual void MocResetCaches() = 0;
  virtual void MocUndoCached() = 0;
  virtual indextype MocOrderSize(const int o) = 0;
  virtual void write_counts(FILE *f) = 0;
  Vocabulary vocab;

  inline void set_sentence_boundary_symbol(std::string s) {
    m_sent_boundary = vocab.word_index(s);
    if (!m_sent_boundary) {
      fprintf(stderr,
              "No \"<s>\" in history, --clear_history cannot be used. Exit.\n");
      exit(-1);
    }
  }
  inline int get_sentence_boundary_symbol() { return (m_sent_boundary); }

  virtual void set_treshold(float x) { assert(false); }

  typedef float disc;
  class disc3 {
    // FIXME: typedef to std::vector<float> ?
    float disc[3];

  public:
    inline float &operator[](const int i) { return (disc[i]); }
  };

  bool zeroprobgrams;
  bool prune_with_real_counts;

  virtual int debug() { return (-40); }
  virtual void print_matrix(int o) {}
  size_t input_data_size;

  std::vector<int> cutoffs;
  inline int cutoff(int x) {
    if (!cutoffs.size())
      return 0;
    if (x > cutoffs.size())
      return cutoffs.back();
    assert(x > 0);
    return cutoffs[x - 1];
  }
  bool discard_cutoffs;
  bool discard_ngrams_with_unk;
  float model_cost_scale;

protected:
  // NgramCounts *m_ng;
  int m_order;
  virtual inline void re_estimate_needed() {}
  int m_sent_boundary;
  int m_ehist_estimate;
  bool m_absolute_discounting;
  std::string m_data_name, m_opti_name, m_prunedata_name;
  virtual void initialize_minmax() = 0;
  std::vector<float> m_minvals, m_maxvals;
};

template <typename KT> class InterKn_k : public InterKn {
public:
  InterKn_k(const bool abs, const std::string &dataname,
            const std::string &optiname, const std::string &prunesource)
      : InterKn(abs, dataname, optiname, prunesource), m_optistorage(NULL) {}
  virtual ~InterKn_k() {}
  virtual double tableprob(std::vector<KT> &indices) = 0;
  virtual bool MocNextVector(std::vector<KT> &v) = 0;
  virtual void remove_sent_start_prob() = 0;
  // virtual void prune_model(float treshold, bool recorrect_kn, Storage<KT>
  // *real_counts)=0;

protected:
  Storage<KT> *m_optistorage;
  inline bool KT_is_short(int *x) { return (false); }
  inline bool KT_is_short(unsigned short *x) { return (true); }
};

template <typename KT, typename CT> class InterKn_t : public InterKn_k<KT> {
public:
  InterKn_t(const bool abs, const std::string &datasource,
            const std::string &optisource, const std::string &prunesource)
      : InterKn_k<KT>(abs, datasource, optisource, prunesource),
        m_new_treshold(1), m_ori_treshold(1), m_eval_cache(NULL) {}
  ~InterKn_t();
  void constructor_helper(const std::string &vocabname, const int read_counts,
                          const int order, const int ndrop, const int nfirst,
                          Storage_t<KT, CT> *datastorage,
                          const indextype hashsize,
                          const std::string &sent_boundary);
  virtual void estimate_bo_counts(bool zerosent = false);
  void estimate_bo_counts_absolute_discounting(bool zerosent = false);
  void counts2lm(TreeGram *lm);
  virtual void counts2asciilm(FILE *out);
  void counts2ascii(FILE *out); // For debugging
  float evaluate(std::vector<float> &discounts);

  double tableprob(std::vector<KT> &indices);
  void find_coeffs(float brak = -0.1, float precision = 1e-3,
                   float lin_precision = 2e-2);

  double logprob_file(const char *name);
  double logprob_datastorage(const Storage<KT> &data);
  // double model_MDL_cost();
  void clear_lm_sentence_boundaries();

  MultiOrderCounts<KT, CT> *moc;
  virtual void add_zeroprob_grams() = 0;

  /* Added after templatization */
  inline bool MocNextVector(std::vector<KT> &v) { return moc->NextVector(v); }
  inline void MocUndoCached() { moc->UndoCached(); }
  inline indextype MocOrderSize(const int o) { return moc->order_size(o); }
  inline void MocResetCaches() { moc->ResetCaches(); }

  inline indextype num_grams() {
    indextype n_grams = 0;
    for (int i = 1; i <= moc->order(); i++) {
      n_grams += moc->order_size(i);
      // fprintf(stderr,"numg %d=%d\n",i,moc->order_size(i));
    }
    return (n_grams);
  }

  void create_model(float prunetreshold);
  void remove_zeroprob_grams();
  virtual void add_counts_for_backoffs() = 0;
  virtual void remove_sent_start_prob() { assert(false); }
  inline void write_counts(FILE *f) { moc->WriteCounts(f); }

  virtual void prune_model(float treshold, bool recorrect_kn,
                           Storage_t<KT, CT> *real_counts) = 0;

protected:
  inline float kn_prob(const int order, const KT *i);
  virtual float kn_prob(const int order, const KT *i, const CT num) = 0;
  virtual float kn_coeff(const int order, const KT *i) = 0;

  virtual void disc2flatv(std::vector<float> &v) = 0;
  virtual float flatv2disc(std::vector<float> &v) = 0;
  std::vector<float>
  calculate_leaveoneout_discounts(int order, std::vector<float> cur_disc);

  CT m_new_treshold;
  CT m_ori_treshold; // Ugly hack...

  sikMatrix<float, float> *m_eval_cache;
  // NgramCounts_t<KT, CT> *m_ng_typed; // Just to avoid casting

  // virtual inline sikMatrix<KT, CT> *get_ct_matrix(int o, CT *foo, DT *bar)
  // {assert(false);return(0);}
  template <typename BOT> void add_counts_for_backoffs_fbase(BOT *);
  template <typename BOT> void add_zeroprob_grams_fbase(BOT *);
  template <typename BOT>
  void prune_model_fbase(float treshold, bool recorrect_kn,
                         Storage_t<KT, CT> *real_counts, BOT *dummy);
  virtual void prune_gram(std::vector<KT> &v, CT num, bool recorrect_kn,
                          MultiOrderCounts_counter_types::bo_c<CT> *dummy) {
    assert(false);
  }
  virtual void prune_gram(std::vector<KT> &v, CT num, bool recorrect_kn,
                          MultiOrderCounts_counter_types::bo_3c<CT> *dummy) {
    assert(false);
  }
};

template <typename KT, typename ICT>
class InterKn_int_disc : public InterKn_t<KT, ICT> {
public:
  InterKn_int_disc(const bool abs, const std::string data,
                   const std::string vocab, const std::string optisource,
                   const int read_counts, const int order, const int ndrop,
                   const int nfirst, Storage_t<KT, ICT> *datastorage,
                   const std::string prunedata_name,
                   const std::string sent_boundary,
                   const indextype hashsize = 3000000);
  virtual inline void init_disc(float x);
  virtual void estimate_nzer_counts();
  virtual inline void prune_model(float treshold, bool recorrect_kn,
                                  Storage_t<KT, ICT> *real_counts) {
    this->prune_model_fbase(treshold, recorrect_kn, real_counts,
                            (MultiOrderCounts_counter_types::bo_c<ICT> *)NULL);
  }
  virtual float kn_prob(const int order, const KT *i, const ICT num);
  virtual float kn_coeff(const int order, const KT *i);

  virtual inline void remove_sent_start_prob() {
    remove_sent_start_prob_fbase(
        (MultiOrderCounts_counter_types::bo_c<ICT> *)NULL);
  }
  virtual inline void add_counts_for_backoffs() {
    InterKn_t<KT, ICT>::add_counts_for_backoffs_fbase(
        (MultiOrderCounts_counter_types::bo_c<ICT> *)NULL);
  }
  virtual inline void add_zeroprob_grams() {
    InterKn_t<KT, ICT>::add_zeroprob_grams_fbase(
        (MultiOrderCounts_counter_types::bo_c<ICT> *)NULL);
  }
  void print_matrix(int o) {
    if (o < this->moc->m_counts.size() && o >= 1)
      this->moc->m_counts[o]->printmatrix();
  }
  inline void disc2flatv(std::vector<float> &v);
  inline float flatv2disc(std::vector<float> &v);
  virtual inline void set_leaveoneout_discounts(int order) {
    std::vector<float> cur_discount(1, m_discount[order]);
    m_discount[order] =
        this->calculate_leaveoneout_discounts(order, cur_discount)[0];
  }

  void set_order(int o);
  std::vector<InterKn::disc> m_discount;

protected:
  template <typename BOT> void remove_sent_start_prob_fbase(BOT *dummy);
  virtual void prune_gram(std::vector<KT> &v, ICT num, bool recorrect_kn,
                          MultiOrderCounts_counter_types::bo_c<ICT> *bo);

  virtual void initialize_minmax() {
    this->m_minvals.resize(this->m_order, 0.0);
    this->m_maxvals.resize(this->m_order, 1.0);
  }
};

template <typename KT, typename ICT>
class InterKn_int_disc3 : public InterKn_t<KT, ICT> {
public:
  InterKn_int_disc3(const bool abs, const std::string data,
                    const std::string vocab, const std::string optisource,
                    const int read_counts, const int order, const int ndrop,
                    const int nfirst, Storage_t<KT, ICT> *datastorage,
                    const std::string prunedata_name,
                    const std::string sent_boundary,
                    const indextype hashsize = 3000000);

  inline void init_disc(float x);
  void estimate_nzer_counts();

  // inline sikMatrix <KT, int> *get_ct_matrix(int o, int *foo, InterKn::disc3
  // *bar) {return this->moc->m_counts[o];}

  inline void prune_model(float treshold, bool recorrect_kn,
                          Storage_t<KT, ICT> *real_counts) {
    this->prune_model_fbase(treshold, recorrect_kn, real_counts,
                            (MultiOrderCounts_counter_types::bo_3c<ICT> *)NULL);
  }
  float kn_prob(const int order, const KT *i, const ICT num);
  virtual inline float kn_coeff(const int order, const KT *i) {
    return kn_coeff_3nzer(order, i,
                          (MultiOrderCounts_counter_types::bo_3c<ICT> *)NULL);
  }
  int debug() {
    KT i = 1;
    return (this->moc->GetBackoffDen(1, &i));
  }
  void print_matrix(int o) {
    if (o < this->moc->m_counts.size() && o >= 1)
      this->moc->m_counts[o]->printmatrix();
  }
  virtual inline void remove_sent_start_prob() {
    remove_sent_start_prob_fbase(
        (MultiOrderCounts_counter_types::bo_3c<ICT> *)NULL);
  }
  inline void disc2flatv(std::vector<float> &v);
  inline float flatv2disc(std::vector<float> &v);

  virtual inline void set_leaveoneout_discounts(int order) {
    std::vector<float> cur_discount(&m_discount[order][0],
                                    &m_discount[order][0] + 3);
    std::vector<float> loo_d(
        this->calculate_leaveoneout_discounts(order, cur_discount));
    for (int i = 0; i < 3; i++) {
      m_discount[order][i] = loo_d[i];
    }
  }

  void set_order(int o);
  std::vector<InterKn::disc3> m_discount;

  virtual inline void add_counts_for_backoffs() {
    this->add_counts_for_backoffs_fbase(
        (MultiOrderCounts_counter_types::bo_3c<ICT> *)NULL);
  }
  virtual inline void add_zeroprob_grams() {
    this->add_zeroprob_grams_fbase(
        (MultiOrderCounts_counter_types::bo_3c<ICT> *)NULL);
  }

protected:
  virtual void prune_gram(std::vector<KT> &v, ICT num, bool recorrect_kn,
                          MultiOrderCounts_counter_types::bo_3c<ICT> *dummy);
  template <typename BOT>
  float kn_coeff_3nzer(const int order, const KT *i, const BOT *dummy = NULL);
  template <typename BOT> void remove_sent_start_prob_fbase(BOT *dummy);

  virtual void initialize_minmax() {
    this->m_minvals.resize(3 * this->m_order, 0.0);
    this->m_maxvals.resize(3 * this->m_order);
    for (int i = 0; i < this->m_order; i++) {
      this->m_maxvals[3 * i] = 1;
      this->m_maxvals[3 * i + 1] = 2;
      this->m_maxvals[3 * i + 2] = 3;
    }
  }
};

#include "InterKn_tmpl.hh"
#endif

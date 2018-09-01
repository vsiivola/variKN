// Library for storing and modifying the n-gram counts
#ifndef MULTIORDERCOUNTS_HH
#define MULTIORDERCOUNTS_HH

#include "InterKn.hh"
#include "Storage.hh"
#include "sikMatrix.hh"
#include <algorithm>
#include <cassert>
#include <vector>

namespace MultiOrderCounts_counter_types {
// Different smoothing schemes defined here (bo_*)
template <typename CT> struct bo_c {  // back_off-counts
  CT den;                             // denominators
  CT nzer;                            // non-zero counts
  CT prune_den;                       // for pruning
  bo_c operator+=(const bo_c other) { // defined for sikMatrix.increment()
    den += other.den;
    nzer += other.nzer;
    prune_den += other.prune_den;
    return *this;
  }

  bo_c operator-=(const bo_c other) { // defined for UndoCached()
    den -= other.den;
    nzer -= other.nzer;
    prune_den -= other.prune_den;
    return *this;
  }
};

template <typename CT> struct bo_c_fp { // back_off-counts
  CT den;                               // denominators
  float nzer;                           // non-zero counts
  CT lost;                              // for pruning
  CT lost_den;

  bo_c_fp operator+=(const bo_c_fp other) { // defined for sikMatrix.increment()
    den += other.den;
    nzer += other.nzer;
    lost += other.lost;
    lost_den += other.lost_den;
    return *this;
  }

  bo_c_fp operator-=(const bo_c_fp other) { // defined for UndoCached()
    den -= other.den;
    nzer -= other.nzer;
    lost -= other.lost;
    lost_den -= other.lost_den;
    return *this;
  }
};

template <typename CT> struct bo_3c {
  // bo_3c() : den(0), nzer({0, 0, 0}), prune_den(0) {}
  CT den;
  CT nzer[3];
  CT prune_den;

  bo_3c operator+=(const bo_3c other) { // defined for sikMatrix.increment()
    den += other.den;
    nzer[0] += other.nzer[0];
    nzer[1] += other.nzer[1];
    nzer[2] += other.nzer[2];
    prune_den += other.prune_den;
    return *this;
  }
  bo_3c operator-=(const bo_3c other) { // defined for UndoCached
    den -= other.den;
    nzer[0] -= other.nzer[0];
    nzer[1] -= other.nzer[1];
    nzer[2] -= other.nzer[2];
    prune_den -= other.prune_den;
    return *this;
  }
};

template <typename CT> struct bo_3c_fp {
  CT den;
  CT nzer[3];
  CT prune_den;
  CT prune_den_left[3];
  CT prune_den_den[3];

  bo_3c_fp
  operator+=(const bo_3c_fp other) { // defined for sikMatrix.increment()
    den += other.den;
    nzer[0] += other.nzer[0];
    nzer[1] += other.nzer[1];
    nzer[2] += other.nzer[2];
    prune_den += other.prune_den;
    prune_den_left[0] += other.prune_den_left[0];
    prune_den_left[0] += other.prune_den_left[1];
    prune_den_left[0] += other.prune_den_left[2];
    prune_den_den[0] += other.prune_den_den[0];
    prune_den_den[1] += other.prune_den_den[1];
    prune_den_den[2] += other.prune_den_den[2];
    return *this;
  }
  bo_3c_fp operator-=(const bo_3c_fp other) { // defined for UndoCached
    den -= other.den;
    nzer[0] -= other.nzer[0];
    nzer[1] -= other.nzer[1];
    nzer[2] -= other.nzer[2];
    prune_den += other.prune_den;
    prune_den_left[0] -= other.prune_den_left[0];
    prune_den_left[0] -= other.prune_den_left[1];
    prune_den_left[0] -= other.prune_den_left[2];
    prune_den_den[0] -= other.prune_den_den[0];
    prune_den_den[1] -= other.prune_den_den[1];
    prune_den_den[2] -= other.prune_den_den[2];
    return *this;
  }
};
}; // namespace MultiOrderCounts_counter_types

template <typename KT, typename CT, typename BOT>
class MultiOrderCounts_typed_interfaces {
public:
  virtual ~MultiOrderCounts_typed_interfaces() {}
  virtual void GetBackoff(const int order, const KT *v, BOT *value) {
    assert(false);
  }
  virtual void IncrementBackoffCache(const int order, const KT *v,
                                     const BOT *value) {
    assert(false);
  }
  virtual void IncrementBackoff(const int order, const KT *v,
                                const BOT *value) {
    assert(false);
  }
  virtual void IncrementBackoff(const std::vector<KT> &v, const BOT *value) {
    assert(false);
  }
  virtual void *StepBackoffsOrder(const bool init, const int order, KT *indices,
                                  BOT *value) {
    assert(false);
    return (NULL);
  }
  virtual void zero_bo(BOT &z) { assert(false); }
};

template <typename KT, typename CT>
class MultiOrderCounts
    : public MultiOrderCounts_typed_interfaces<
          KT, CT, MultiOrderCounts_counter_types::bo_c<CT>>,
      public MultiOrderCounts_typed_interfaces<
          KT, CT, MultiOrderCounts_counter_types::bo_c_fp<CT>>,
      public MultiOrderCounts_typed_interfaces<
          KT, CT, MultiOrderCounts_counter_types::bo_3c<CT>>,
      public MultiOrderCounts_typed_interfaces<
          KT, CT, MultiOrderCounts_counter_types::bo_3c_fp<CT>> {
public:
  MultiOrderCounts()
      : vocabsize(1000000), hashsize(0), m_cur_order(1), m_cur_ng(0) {}
  virtual ~MultiOrderCounts();

  typedef MultiOrderCounts_counter_types::bo_c<CT> bo_c;
  typedef MultiOrderCounts_counter_types::bo_c_fp<CT> bo_c_fp;
  typedef MultiOrderCounts_counter_types::bo_3c<CT> bo_3c;
  typedef MultiOrderCounts_counter_types::bo_3c_fp<CT> bo_3c_fp;

  virtual int order() = 0;
  inline indextype order_size(int o) {
    if (o <= order())
      return (m_counts[o]->num_entries());
    return (0);
  }
  virtual inline indextype bo_order_size(int o) = 0;

  virtual void WriteCounts(FILE *out) = 0;
  virtual void ReadCounts(FILE *in) = 0;
  virtual void RemoveOrder(int order) = 0;
  long InitializeCountsFromText(FILE *in, Vocabulary *vocab,
                                const bool grow_vocab, const int read_order,
                                const std::string &sent_start_sym);
  long InitializeCountsFromStorage(Storage_t<KT, CT> *data,
                                   const int read_order,
                                   const int sent_start_idx);
  void UseAsCounts(sikMatrix<KT, CT> *mat);
  bool NextVector(std::vector<KT> &v);
  void RandomVector(std::vector<KT> &v);
  void RandomVector(std::vector<KT> &v, Storage_t<KT, CT> &data);
  virtual void clear_derived_counts() = 0;
  virtual void RemoveDefaultBackoffs() = 0;

  inline void *StepCountsOrder(const bool init, const int order, KT *indices,
                               CT *value);

  inline void DeleteCurrentST(const int order) {
    m_counts[order]->delete_current_st();
  }
  inline void *OrderedStepCountsOrder(const bool init, const int order,
                                      KT *indices, CT *value);
  // inline void RemoveEmptyNodes(const int order)
  // {RemoveEmptyNodes(order,0,0);}

  /* Manipulation of counts other than nzer and prune_dennn*/
  inline CT GetCount(const std::vector<KT> &v);
  inline CT GetCount(const int order, const KT *v);
  inline void SetCount(const std::vector<KT> &v, const CT value) {
    SetCount(v.size(), &v[0], value);
  }
  inline void SetCount(const int order, const KT *v, const CT value);
  inline CT IncrementCount(std::vector<KT> &v, const CT value) {
    return (IncrementCount(v.size(), &v[0], value));
  }
  virtual CT GetBackoffDen(const int order, const KT *v) = 0;
  virtual CT GetBackoffNzer(const int order, const KT *v) {
    assert(false);
    return (CT)-1;
  }
  inline CT GetBackoffDen(const std::vector<KT> &v);
  inline CT GetBackoffNzer(const std::vector<KT> &v);

  inline CT IncrementCount(const std::vector<KT> &v, const CT value);
  inline CT IncrementCount(const int order, const KT *v, const CT value);
  virtual void IncrementBackoffDen(const int order, const KT *v,
                                   const CT value) = 0;
  virtual void IncrementBackoffNzer(const int order, const KT *v,
                                    const CT value) {
    assert(false);
  }
  virtual void IncrementBackoffNzer(const int order, const KT *v, const int pos,
                                    const CT value) {
    assert(false);
  }
  /* Cached manipulations */
  virtual void ResetCaches() = 0;
  virtual void UndoCached() = 0;
  CT IncrementCountCache(const int order, const KT *v, const CT value);
  virtual void IncrementBackoffCacheDen(const int order, const KT *v,
                                        const CT value) = 0;
  int vocabsize;
  indextype hashsize;

  std::vector<sikMatrix<KT, CT> *> m_counts; // This should be private

  virtual void IncrementBackoffCacheNzer(const int order, const KT *v,
                                         const CT value) {
    assert(false);
  }
  virtual void IncrementBackoffCacheNzer(const std::vector<KT> &v,
                                         const CT value) {
    assert(false);
  }
  virtual void IncrementBackoffCacheNzer(const int order, const KT *v,
                                         const CT *value) {
    assert(false);
  }
  virtual void IncrementBackoffCacheNzer(const int order, const KT *v,
                                         const int pos, const CT value) {
    assert(false);
  }

#define USE_INTERFACE(x)                                                       \
  using MultiOrderCounts_typed_interfaces<KT, CT, bo_c>::x;                    \
  using MultiOrderCounts_typed_interfaces<KT, CT, bo_c_fp>::x;                 \
  using MultiOrderCounts_typed_interfaces<KT, CT, bo_3c>::x;                   \
  using MultiOrderCounts_typed_interfaces<KT, CT, bo_3c_fp>::x

  USE_INTERFACE(GetBackoff);
  USE_INTERFACE(IncrementBackoffCache);
  USE_INTERFACE(IncrementBackoff);
  USE_INTERFACE(StepBackoffsOrder);
  USE_INTERFACE(zero_bo);
  inline void write_num(FILE *out, const int val) { fprintf(out, "%d", val); }
  inline void write_num(FILE *out, const unsigned int val) {
    fprintf(out, "%u", val);
  }
  inline void write_num(FILE *out, const long val) { fprintf(out, "%ld", val); }
  inline void write_num(FILE *out, const float val) {
    fprintf(out, "%.4f", val);
  }
  inline void write_num(FILE *out, const double val) {
    fprintf(out, "%.4f", val);
  }
  inline void read_num(int *val, const std::string *s, bool *ok) {
    *val = str::str2long(s, ok);
  }

  inline void read_num(long *val, const std::string *s, bool *ok) {
    *val = str::str2long(s, ok);
  }

  inline void read_num(unsigned int *val, const std::string *s, bool *ok) {
    *val = str::str2long(s, ok);
  }

  inline void read_num(float *val, const std::string *s, bool *ok) {
    *val = str::str2float(s, ok);
  }
  inline void read_num(double *val, const std::string *s, bool *ok) {
    *val = str::str2float(s, ok);
  }

  virtual void clear_nzer(int o) = 0;
  virtual void clear_lden(int o) { assert(false); }

protected:
  std::vector<int> m_do_not_delete;
  void allocate_matrices_counts(int o);

  /* For NextVector() */
  int m_cur_order;
  indextype m_cur_ng;

  /* Low-level matrix library functions */
  CT Increment_wo_del(struct matrix *m, const KT *indices, const CT value);

  struct c_cache_t {
    int order;
    CT val;
    indextype index;
  };

  std::vector<c_cache_t> c_cache;
  std::vector<indextype> min_cc_cache;
};

template <typename KT, typename CT, typename BOT>
class MultiOrderCounts_Generic_BOT : public MultiOrderCounts<KT, CT> {
public:
  MultiOrderCounts_Generic_BOT();
  ~MultiOrderCounts_Generic_BOT();

  inline int order() { return MultiOrderCounts<KT, CT>::m_counts.size() - 1; }
  void WriteCounts(FILE *out);
  void ReadCounts(FILE *in);

  void RemoveOrder(int order);
  // void RemoveEmptyNodes(const int order, const indextype start, const
  // indextype bo_start);
  void clear_derived_counts();
  void clear_nzer(int o);

  inline void SetBackoff(const int order, const KT *v, const BOT *value);
  inline void SetBackoff(const std::vector<KT> &v, const BOT *value);
  void GetBackoff(const int order, const KT *v, BOT *value);
  inline CT GetBackoffDen(const std::vector<KT> &v);
  CT GetBackoffDen(const int order, const KT *v);
  void IncrementBackoff(const int order, const KT *v, const BOT *value);
  void IncrementBackoff(const std::vector<KT> &v, const BOT *value);
  void IncrementBackoffDen(const int order, const KT *v, const CT value);

  void ResetCaches();
  void UndoCached();
  void IncrementBackoffCache(const int order, const KT *v, const BOT *value);
  inline void IncrementBackoffCacheDen(const int order, const KT *v,
                                       const CT value);

  virtual inline indextype bo_order_size(int o) {
    if (o < m_backoffs.size())
      return (m_backoffs[o]->num_entries());
    return (0);
  }

  std::vector<sikMatrix<KT, BOT> *> m_backoffs; // This should be private

  void *StepBackoffsOrder(const bool init, const int order, KT *indices,
                          BOT *value);

  inline void RemoveDefaultBackoffs() {
    for (int o = order(); o >= 2; o--) {
      RemoveDefaultValues(m_backoffs[o]->m);
    }
  }

  // These functions are in "wrong" place. Real classes wrap aroud them
  void IncrementBackoffNzer_1nzer(const int order, const KT *v, const CT value);
  inline CT GetBackoffNzer_1nzer(const int order, const KT *v);
  inline void IncrementBackoffCacheNzer_1nzer(const int order, const KT *v,
                                              const CT value);
  inline void GetBackoffNzer_3nzer(const int order, const KT *v, CT *res);
  inline CT GetBackoffNzer_3nzer(const int order, const KT *v, const int which);
  inline void IncrementBackoffNzer_3nzer(const int order, const KT *v,
                                         const int pos, const CT value);
  inline void IncrementBackoffCacheNzer_3nzer(const int order, const KT *v,
                                              const CT *value);
  inline void IncrementBackoffCacheNzer_3nzer(const int order, const KT *v,
                                              const int pos, const CT value);
  inline void zero_bo(typename MultiOrderCounts<KT, CT>::bo_c &z) {
    z = m_bb_init;
  }
  inline void zero_bo(typename MultiOrderCounts<KT, CT>::bo_c_fp &z) {
    z = m_bb_fp_init;
  }
  inline void zero_bo(typename MultiOrderCounts<KT, CT>::bo_3c &z) {
    z = m_3bb_init;
  }
  inline void zero_bo(typename MultiOrderCounts<KT, CT>::bo_3c_fp &z) {
    z = m_3bb_fp_init;
  }

protected:
  virtual void WriteCounts_BOhelper(FILE *out, BOT *bo) = 0;
  virtual void ReadCounts_BOhelper(BOT *bo, std::string *s, bool *ok) = 0;
  void allocate_matrices_backoffs(int o);
  BOT m_uni_bo;

  typename MultiOrderCounts<KT, CT>::bo_c m_bb_init;
  typename MultiOrderCounts<KT, CT>::bo_c_fp m_bb_fp_init;
  typename MultiOrderCounts<KT, CT>::bo_3c m_3bb_init;
  typename MultiOrderCounts<KT, CT>::bo_3c_fp m_3bb_fp_init;

  struct bo_cache_t {
    int order;
    BOT bo;
    indextype index;
  };

  std::vector<bo_cache_t> bo_cache;
  std::vector<indextype> min_bo_cache;

  inline void zero_nz(typename MultiOrderCounts<KT, CT>::bo_c *z) {
    z->nzer = 0;
  }
  inline void zero_nz(typename MultiOrderCounts<KT, CT>::bo_c_fp *z) {
    z->nzer = 0;
  }
  inline void zero_nz(typename MultiOrderCounts<KT, CT>::bo_3c *z) {
    z->nzer[0] = 0;
    z->nzer[1] = 0;
    z->nzer[2] = 0;
  }
  inline void zero_nz(typename MultiOrderCounts<KT, CT>::bo_3c_fp *z) {
    z->nzer[0] = 0;
    z->nzer[1] = 0;
    z->nzer[2] = 0;
  }
  BOT bo_init;
};

template <typename KT, typename CT>
class MultiOrderCounts_1nzer
    : public MultiOrderCounts_Generic_BOT<
          KT, CT, typename MultiOrderCounts<KT, CT>::bo_c> {
public:
  typedef typename MultiOrderCounts<KT, CT>::bo_c bo_c;
  inline void IncrementBackoffNzer(const int order, const KT *v,
                                   const CT value) {
    this->IncrementBackoffNzer_1nzer(order, v, value);
  }
  inline CT GetBackoffNzer(const int order, const KT *v) {
    return this->GetBackoffNzer_1nzer(order, v);
  }
  inline void IncrementBackoffCacheNzer(const int order, const KT *v,
                                        const CT value) {
    this->IncrementBackoffCacheNzer_1nzer(order, v, value);
  }

private:
  inline void WriteCounts_BOhelper(FILE *out, bo_c *bo) {
    this->write_num(out, bo->nzer);
    fprintf(out, " ");
    this->write_num(out, bo->den);
    fprintf(out, " ");
    this->write_num(out, bo->prune_den);
  }

  inline void ReadCounts_BOhelper(bo_c *bo, std::string *s, bool *ok) {
    bo->nzer = str::str2long(s, ok);
    bo->den = str::str2long(s + 1, ok);
    bo->prune_den = str::str2long(s + 2, ok);
  }
};

template <typename KT, typename CT>
class MultiOrderCounts_1nzer_fp
    : public MultiOrderCounts_Generic_BOT<
          KT, CT, typename MultiOrderCounts<KT, CT>::bo_c_fp> {
public:
  typedef typename MultiOrderCounts<KT, CT>::bo_c_fp bo_c_fp;
  inline void IncrementBackoffNzer(const int order, const KT *v,
                                   const CT value) {
    IncrementBackoffNzer_1nzer(order, v, value);
  }
  inline CT GetBackoffNzer(const int order, const KT *v) {
    // return GetBackoffNzer_1nzer(order, v);
    assert(false);
  }
  inline void IncrementBackoffCacheNzer(const int order, const KT *v,
                                        const CT value) {
    IncrementBackoffCacheNzer_1nzer(order, v, value);
  }

private:
  void clear_lden(const int o);
  inline void WriteCounts_BOhelper(FILE *out, bo_c_fp *bo) { assert(false); }
  inline void ReadCounts_BOhelper(bo_c_fp *bo, std::string *s, bool *ok) {
    assert(false);
  }
};

template <typename KT, typename CT>
class MultiOrderCounts_3nzer
    : public MultiOrderCounts_Generic_BOT<
          KT, CT, typename MultiOrderCounts<KT, CT>::bo_3c> {
public:
  typedef typename MultiOrderCounts<KT, CT>::bo_3c bo_3c;

  inline void GetBackoffNzer(const int order, const KT *v, CT *res) {
    GetBackoffNzer_3nzer(order, v, res);
  }

  inline CT GetBackoffNzer(const int order, const KT *v, const int which) {
    return GetBackoffNzer_3nzer(order, v, which);
  }
  inline void IncrementBackoffNzer(const int order, const KT *v, const int pos,
                                   const CT value) {
    this->IncrementBackoffNzer_3nzer(order, v, pos, value);
  }
  // void RemoveBackDefs(struct matrix *m);

  /* Caching functions (for fast deletion of recently added entried */
  inline void IncrementBackoffCacheNzer(const int order, const KT *v,
                                        const CT *value) {
    this->IncrementBackoffCacheNzer_3nzer(order, v, value);
  }

  inline void IncrementBackoffCacheNzer(const int order, const KT *v,
                                        const int pos, const CT value) {
    this->IncrementBackoffCacheNzer_3nzer(order, v, pos, value);
  }

private:
  inline void WriteCounts_BOhelper(FILE *out, bo_3c *bo) {
    fprintf(out, "%ld %ld %ld %ld %ld", (long)bo->nzer[0], (long)bo->nzer[1],
            (long)bo->nzer[2], (long)bo->den, (long)bo->prune_den);
  }
  inline void ReadCounts_BOhelper(bo_3c *bo, std::string *s, bool *ok) {
    bo->nzer[0] = str::str2long(s, ok);
    // fprintf(stderr,"nz %s=%d(%d)\n",s->c_str(), bo->nzer[0], ok);
    bo->nzer[1] = str::str2long(s + 1, ok);
    // fprintf(stderr,"nz %s=%d(%d)\n",(s+1)->c_str(), bo->nzer[1], ok);
    bo->nzer[2] = str::str2long(s + 2, ok);
    // fprintf(stderr,"nz %s=%d(%d)\n",(s+2)->c_str(), bo->nzer[2], ok);
    bo->den = str::str2long(s + 3, ok);
    // fprintf(stderr,"den %s=%d(%d)\n",(s+3)->c_str(), bo->den, ok);
    bo->prune_den = str::str2long(s + 4, ok);
  }
};

template <typename KT, typename CT>
class MultiOrderCounts_3nzer_fp
    : public MultiOrderCounts_Generic_BOT<
          KT, CT, typename MultiOrderCounts<KT, CT>::bo_3c_fp> {
  typedef typename MultiOrderCounts<KT, CT>::bo_3c_fp bo_3c_fp;
  inline void IncrementBackoffNzer(const int order, const KT *v, const int pos,
                                   const CT value) {
    IncrementBackoffNzer_3nzer(order, v, pos, value);
  }

private:
  inline void WriteCounts_BOhelper(FILE *out, bo_3c_fp *bo) { assert(false); }
  inline void ReadCounts_BOhelper(bo_3c_fp *bo, std::string *s, bool *ok) {
    assert(false);
  }
};

#include "MultiOrderCounts_tmpl.hh"
#endif

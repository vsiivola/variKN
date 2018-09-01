// Helper library for storing and modifying the n-gram counts
#ifndef NGRAMCOUNTS_HH
#define NGRAMCOUNTS_HH

#include "Storage.hh"
#include "Vocabulary.hh"
#include "def.hh"
#include "sikMatrix.hh"

class NgramCounts {
public:
  inline NgramCounts() : vocab(new Vocabulary), m_destruct_vocab(true) {}
  virtual inline ~NgramCounts() {
    if (m_destruct_vocab)
      delete vocab;
  }
  virtual void read(FILE *, FILE *) = 0;
  virtual void shrink(float ndrop, int nfirst) = 0;
  virtual int order() = 0;
  virtual void read_vocab(FILE *vocab) = 0;
  virtual long count(FILE *, bool grow_vocab = true) = 0;
  inline void use_vocabulary(Vocabulary *V) {
    if (m_destruct_vocab)
      delete vocab;
    m_destruct_vocab = false;
    vocab = V;
  }
  Vocabulary *vocab;

protected:
  int m_max_vocab;

  /* Needed for reading and writing from template */
  inline int ascii2num(char *c, int *a) { return (atoi(c)); }
  inline long ascii2num(char *c, long *a) {
    return (strtol(c, NULL, 10));
  } // FIXME: This can be broken
  // inline unsigned int ascii2num(char *c, unsigned int *a) {return(atoi(c));}
  inline float ascii2num(char *c, float *a) { return ((float)strtod(c, NULL)); }
  inline float ascii2num(char *c, double *a) {
    return ((float)strtod(c, NULL));
  }
  inline void write_num(FILE *out, const int val) { fprintf(out, "%d", val); }
  inline void write_num(FILE *out, const long val) { fprintf(out, "%ld", val); }
  inline void write_num(FILE *out, const float val) {
    fprintf(out, "%.4f", val);
  }
  bool m_destruct_vocab;
};

template <typename K, typename V> class NgramCounts_t : public NgramCounts {
public:
  NgramCounts_t(const int n, const int max_vocab, indextype hashsize);
  ~NgramCounts_t();
  long count(FILE *, bool grow_vocab = true);
  long count(Storage_t<K, V> *data);
  /* These function names are identical to those of
     vocabulary class, this is highly confusing. Should be fixed */
  void read(FILE *, FILE *vocab);
  void read_vocab(FILE *vocab);
  void write(FILE *, FILE *vocab, bool sort);

  void shrink(float ndrop, int nfirst);
  inline int order() { return m_ds.size(); }
  // bool next(bool init,std::vector<int> &vec, V &value);
  sikMatrix<K, V> *counts;
  // void convert_clustered(ClusterMap<K> *clmap);

private:
  std::vector<K> m_ds;

  /* For shrink */
  struct sortstruct {
    sortstruct();
    int oldidx;
    V count;
    std::string wstring;
    bool operator<(const sortstruct &other) const { // defined for sort()
      return (count > other.count);
    }
  };
};

template <typename K, typename V>
NgramCounts_t<K, V>::sortstruct::sortstruct() : oldidx(-1), count(0) {}

#include "NgramCounts_tmpl.hh"
#endif

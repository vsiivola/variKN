// Copyright (C) 2007  Vesa Siivola. 
// See licence.txt for the terms of distribution.

// Easily modifiable presentation of n-gram language model
#ifndef HASHGRAM_HH
#define HASHGRAM_HH

#include <cstdio>
#include "NGram.hh"
#include "sikMatrix.hh"
#include "ArpaReader.hh"

class HashGram : public NGram, public ArpaReader {
public:
  HashGram () : m_print_zerograms(false) {}
  virtual ~HashGram() {}
  virtual void read(FILE *in, bool binary=false) = 0;
  virtual void write(FILE *out, bool binary=false) = 0;
  virtual void remove_empty_grams()=0;
  virtual void prune(float treshold)=0;
  virtual void add_zeroprob_grams() = 0;

protected:
  struct twofloat {
    twofloat(float a, float b): float1(a), float2(b) {}
    twofloat() {}
    float float1;
    float float2;
    twofloat operator+=(const twofloat other) {
      float1+=other.float1;
      float2+=other.float2;
      return *this;
    }
  };

  bool m_print_zerograms;
};

template <typename KT>
class HashGram_t : public HashGram {
public:
  ~HashGram_t();
  inline void read(FILE *in, bool binary=false) {
    if (binary) {
      fprintf(stderr,"HashGram: Reading binary format unsupported. Exit.\n");
      exit(-1);
    }
    read_real(in);
  }
  void read_real(FILE *in);
  void write(FILE *out, bool binary=false) {
    if (binary) {
      fprintf(stderr,"HashGram: Writing binary format unsupported. Exit.\n");
      exit(-1);
    }
    write_real(out);
  }
  void write_real(FILE *out);
  void remove_empty_grams();
  void prune(float treshold);
  void add_zeroprob_grams();

private:
  std::vector<sikMatrix<KT, float> *> probs;
  std::vector<sikMatrix<KT, float> *> backoffs;

  inline float log_prob_bo(const std::vector<int> &gram); // backoff, default
  inline float log_prob_bo_cl(const std::vector<int> &gram); // Clustered backoff
  inline float log_prob_i(const std::vector<int> &gram); // Interpolated

  inline float log_prob_bo(const std::vector<unsigned short> &gram);
  inline float log_prob_i(const std::vector<unsigned short> &gram); 
  float log_prob_bo_helper(const std::vector<KT> &gram); 
  float log_prob_i_helper(const std::vector<KT> &gram);

  inline float log_prob_bo(const Gram &gram) {
    std::vector<KT> g(gram.size());
    for (size_t i=0;i<gram.size();i++) g[i]=gram[i];
    return(log_prob_bo(g));
  }

  inline float log_prob_i(const Gram &gram){
    std::vector<KT> g(gram.size());
    for (size_t i=0;i<gram.size();i++) g[i]=gram[i];
    return(log_prob_i(g));
  }

};

#include "HashGram_tmpl.hh"
#endif

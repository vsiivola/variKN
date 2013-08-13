// Present and interpolated model through the simple NGram interface
#ifndef INTERTREEGRAM_HH
#define INTERTREEGRAM_HH

#include "TreeGram.hh"
#include "ArpaReader.hh"

class InterTreeGram : public NGram {
public:
  InterTreeGram ( std::vector< std::string >, std::vector<float> );
  ~InterTreeGram ( );

  float log_prob(const Gram &gram);

  // NGram.hh wants us to implement these, but these are actually not needed
  void read(FILE *, bool) { assert(false); }
  void write(FILE *, bool, std::string) { assert(false); }
  float log_prob_bo(const std::vector<int> &gram) { assert(false); } // backoff, default
  float log_prob_i(const std::vector<int> &gram) { assert(false); } // Interpolated

  float log_prob_bo(const std::vector<unsigned short> &gram) { assert(false); } // backoff, default
  float log_prob_i(const std::vector<unsigned short> &gram) { assert(false); } // Interpolated

  inline float log_prob_bo(const Gram &gram) { log_prob(gram); } // Keep this version lean and mean
  float log_prob_i(const Gram &gram) { assert(false); } // Interpolated

  void fetch_bigram_list(int, std::vector<int>&, std::vector<float>&); // For speech recognition LM lookahead
  void test_write(std::string fname, int idx);
private:
  std::vector<TreeGram *> m_models;
  std::vector<float> m_coeffs;
};
#endif

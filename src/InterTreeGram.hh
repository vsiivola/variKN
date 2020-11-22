// Present and interpolated model through the simple NGram interface
#ifndef INTERTREEGRAM_HH
#define INTERTREEGRAM_HH

#include <stdexcept>
#include "ArpaReader.hh"
#include "TreeGram.hh"

class NotImplemented : public std::logic_error {
public:
  NotImplemented() : std::logic_error("Function not yet implemented"){};
};

class InterTreeGram : public NGram {
public:
  InterTreeGram(std::vector<std::string>, std::vector<float>);
  ~InterTreeGram();

  float log_prob(const Gram &gram);

  // NGram.hh wants us to implement these, but these are actually not needed
  void read(FILE *, bool) { throw NotImplemented(); }
  void write(FILE *, bool, std::string) { throw NotImplemented(); }
  float log_prob_bo(const std::vector<int> &gram) {
    throw NotImplemented();
  } // backoff, default
  float log_prob_i(const std::vector<int> &gram) {
    throw NotImplemented();
  } // Interpolated

  float log_prob_bo(const std::vector<unsigned short> &gram) {
    throw NotImplemented();
  } // backoff, default
  float log_prob_i(const std::vector<unsigned short> &gram) {
    throw NotImplemented();
  } // Interpolated

  inline float log_prob_bo(const Gram &gram) {
    return log_prob(gram);
  } // Keep this version lean and mean
  float log_prob_i(const Gram &gram) { throw NotImplemented(); } // Interpolated

  void fetch_bigram_list(
      int, std::vector<int> &,
      std::vector<float> &); // For speech recognition LM lookahead
  void test_write(std::string fname, int idx);

private:
  std::vector<TreeGram *> m_models;
  std::vector<float> m_coeffs;
};
#endif

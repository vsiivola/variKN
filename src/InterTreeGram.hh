// Copyright (C) 2007  Vesa Siivola. 
// See licence.txt for the terms of distribution.

// Present and interpolated model through the simple NGram interface
#ifndef INTERTREEGRAM_HH
#define INTERTREEGRAM_HH

#include "TreeGram.hh"
#include "ArpaReader.hh"

class InterTreeGram : public NGram {
  //friend TreeGram;
  
public:
  InterTreeGram ( std::vector< std::string > );
  ~InterTreeGram ( );
  void read(FILE *, bool) { assert(false); }
  void write(FILE *, bool) { assert(false); }
  float log_prob_bo(const std::vector<int> &gram) { assert(false); } // backoff, default
  float log_prob_i(const std::vector<int> &gram) { assert(false); } // Interpolated

  float log_prob_bo(const std::vector<unsigned short> &gram) { assert(false); } // backoff, default
  float log_prob_i(const std::vector<unsigned short> &gram) { assert(false); } // Interpolated

  float log_prob_bo(const Gram &gram) { assert(false); } // Keep this version lean and mean
  float log_prob_i(const Gram &gram) { assert(false); } // Interpolated

private:
  void initialize_models(std::vector< std::string > );

  std::vector<TreeGram *> m_models;
};
#endif

#include <cstdio>

#include "InterTreeGram.hh"

InterTreeGram::InterTreeGram(std::vector< std::string > lm_names) {
  initialize_vocab(lm_names);

}

void InterTreeGram::initialize_vocab(std::vector< std::string > lm_names) {
  for ( std::vector<std::string>::iterator it = lm_names.begin(); it != lm_names.end(); it ++) {
    fprintf(stdout, "INTERTREEGRAM %s\n", it->c_str());
  }
}

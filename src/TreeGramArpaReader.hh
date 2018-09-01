// Routines for reading and writing arpa format files from and to the
// internal prefix tree format.
#ifndef TREEGRAMARPAREADER_HH
#define TREEGRAMARPAREADER_HH

#include "TreeGram.hh"
#include <stdio.h>

class TreeGramArpaReader {
public:
  void read(FILE *file, TreeGram *tree_gram, bool add_missing_unigrams = false);
  void write(FILE *file, TreeGram *tree_gram,
             std::string field_separator = " ");
  void write_interpolated(FILE *file, TreeGram *treegram,
                          std::string field_separator = " ");
};

#endif /* TREEGRAMARPAREADER_HH */

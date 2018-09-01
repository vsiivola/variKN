// Copyright (C) 2007  Vesa Siivola.
// See licence.txt for the terms of distribution.

// Get the highest order n-gram counts from corpus.
#include "NgramCounts.hh"
#include "conf.hh"
#include "io.hh"

int main(int argc, char **argv) {
  conf::Config config;
  config("Usage: find_ngram_counts [OPTIONS] textin countsout\nExtracts n-gram "
         "counts from given text.\n")('n', "norder=INT", "arg must", "3",
                                      "n-gram order")(
      'B', "vocabin=FILE", "arg", "", "specify vocabulary")(
      'b', "vocabout=FILE", "arg", "",
      "if chosen, a vocabulary file will be written and the output count file "
      "will only contain the word indices")('H', "hashsize=INT", "arg", "0",
                                            "size of the reserved hash table")(
      's', "sort", "", "sort output n-grams to the ordering of the vocabulary");
  config.parse(argc, argv, 2);

  io::Stream::verbose = true;
  io::Stream in(config.arguments[0], "r");
  io::Stream out(config.arguments[1], "w");
  io::Stream vocab, vocabin;
  if (config["vocabout"].specified)
    vocab.open(config["vocabout"].get_str(), "w");

  int n = config["norder"].get_int();
  int hashs = config["hashsize"].get_int();
  bool sort = config["sort"].specified;
  if (sort)
    fprintf(stderr, " Will sort n-grams\n");

  NgramCounts_t<int, int> nc(n, 0, hashs);
  if (config["vocabin"].specified) {
    vocabin.open(config["vocabin"].get_str(), "r");
    nc.read_vocab(vocabin.file);
  }

  nc.count(in.file, !vocabin.file); // vocabin used as boolean
  nc.write(out.file, vocab.file, sort);

  // Some closes missing
  in.close();
  out.close();
  vocab.close();
  vocabin.close();
}

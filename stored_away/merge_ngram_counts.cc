// Copyright (C) 2007  Vesa Siivola.
// See licence.txt for the terms of distribution.

// Program to merge n-gram count files
#include "NgramCounts.hh"
#include "conf.hh"
#include "io.hh"
#include "str.hh"

int main(int argc, char **argv) {
  conf::Config config;
  config("Usage: merge_ngram_counts [OPTIONS] list_of_input_files "
         "output_counts\nMerges several count files.\n")(
      'b', "vocabout=FILE", "arg", "",
      "if chosen, a vocabulary file will be written and the output count file "
      "will only contain the word indices")(
      's', "sort", "", "sort output n-grams to the ordering of the vocabulary")(
      'n', "norder", "arg must", "3", "order of the n-grams");
  config.parse(argc, argv, 2);

  io::Stream::verbose = true;
  io::Stream filelist(config.arguments[0], "r");
  io::Stream out(config.arguments[1], "w");
  io::Stream vocab;
  if (config["vocabout"].specified)
    vocab.open(config["vocabout"].get_str(), "w");
  int n = config["norder"].get_int();
  bool sort = config["sort"].specified;
  NgramCounts_t<int, int> nc(n, 0, 0);

  /* main loop, read in everything */
  std::string s;
  char *sbuf, *cptr;

  io::Stream cur_voc, cur_tri;
  while (str::read_line(&s, filelist.file, true)) {
    sbuf = strdup(s.c_str()); // Ugly, fix to use str::split()
    fprintf(stderr, "read sbuf %s\n", sbuf);
    /* Open files */
    cptr = strtok(sbuf, " ");
    cur_tri.open(cptr, "r");
    cptr = strtok(NULL, " ");
    if (cptr)
      cur_voc.open(cptr, "r");
    free(sbuf);
    nc.read(cur_tri.file, cur_voc.file);
    cur_tri.close();
    cur_voc.close();
  }
  nc.write(out.file, vocab.file, sort);
}

// Copyright (C) 2007  Vesa Siivola.
// See licence.txt for the terms of distribution.

// An unfinished framework for performing entropy pruning
#include "HashGram.hh"
#include "conf.hh"
#include "io.hh"

int main(int argc, char *argv[]) {
  conf::Config config;
  config("Usage: prune arpain arpaout\nPerforms entropy based pruning.\n")(
      't', "treshold=float", "arg must", "", "Pruning treshold.")(
      's', "smallvocab", "", "",
      "Vocabulary is less than 65000 entries. Saves some memory.");
  config.parse(argc, argv, 2);
  const bool smallvocab = config["smallvocab"].specified;

  io::Stream::verbose = true;
  io::Stream in(config.arguments[0], "r");
  io::Stream out(config.arguments[1], "w");

  HashGram *hg;
  if (smallvocab)
    hg = new HashGram_t<unsigned short>;
  else
    hg = new HashGram_t<int>;

  fprintf(stderr, "Reading\n");
  hg->read(in.file);
  in.close();

  fprintf(stderr, "Pruning\n");
  hg->prune(config["treshold"].get_double());

  fprintf(stderr, "Writing\n");
  hg->write(out.file);
  out.close();
}

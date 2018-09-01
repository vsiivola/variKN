// Copyright (C) 2007  Vesa Siivola.
// See licence.txt for the terms of distribution.

// This program converts arpa firmat language models to binary format.
#include "TreeGram.hh"
#include "conf.hh"
#include "io.hh"

int main(int argc, char **argv) {
  conf::Config config;
  config("Usage: bin2arpa bin_in arpa_out\nConverts binary LMs to arpa.\n");
  config.parse(argc, argv, 2, true);

  io::Stream in(config.arguments.at(0), "r");
  io::Stream out(config.arguments.at(1), "w");

  TreeGram ng;

  fprintf(stderr, "Reading\n");
  ng.read(in.file, true);
  in.close();

  fprintf(stderr, "Writing\n");
  ng.write(out.file);
  out.close();
}

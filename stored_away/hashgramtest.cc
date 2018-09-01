// Copyright (C) 2007  Vesa Siivola.
// See licence.txt for the terms of distribution.

// Test the write and read functions of hashgrams
#include "HashGram.hh"
#include "conf.hh"
#include "io.hh"

int main(int argc, char *argv[]) {
  conf::Config config;
  config("Usage:  hashgramtest in out\nTesting.\n");
  config.parse(argc, argv, 2);

  io::Stream::verbose = true;
  io::Stream in(config.arguments[0], "r");
  io::Stream out(config.arguments[1], "w");

  HashGram_t<int> hg;
  hg.read(in.file);
  in.close();
  fprintf(stderr, "Writing\n");
  hg.write(out.file);
  out.close();
}

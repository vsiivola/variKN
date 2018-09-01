// Copyright (C) 2007  Vesa Siivola.
// See licence.txt for the terms of distribution.

// Test ReadCounts and WriteCounts
#include "InterKn.hh"
#include "conf.hh"
#include "io.hh"
#include "str.hh"

int main(int argc, char **argv) {
  conf::Config config;
  config("Usage: moc_writetest counts out\n")(
      '3', "3nzer", "", "",
      "Use 3 discounts per order instead of one. Takes a bit more memory "
      "during model estimation.");

  config.parse(argc, argv, 2);
  const bool use_3nzer = config["3nzer"].specified;

  ClusterMap<int> clmap;
  io::Stream::verbose = true;
  io::Stream in(config.arguments[0], "r");
  io::Stream out(config.arguments[1], "w");

  MultiOrderCounts<int, int> *moc;
  if (use_3nzer)
    moc = new MultiOrderCounts_3nzer<int, int>;
  else
    moc = new MultiOrderCounts_1nzer<int, int>;
  moc->hashsize = 3000000;
  moc->ReadCounts(in.file);
  moc->WriteCounts(out.file);
}

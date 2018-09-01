#include "TreeGramArpaReader.hh"
#include "conf.hh"
#include "io.hh"

int main(int argc, char **argv) {
  conf::Config config;
  config("Usage: arpa2bin arpa_in bin_out\nConverts arpa LMs to binary.\n");
  config.parse(argc, argv, 2, true);

  io::Stream::verbose = true;
  io::Stream in(config.arguments[0], "r");
  io::Stream out(config.arguments[1], "w");

  TreeGramArpaReader areader;
  TreeGram ng;

  fprintf(stderr, "Reading\n");
  areader.read(in.file, &ng);
  in.close();

  if (ng.get_type() == NGram::INTERPOLATED) {
    fprintf(stderr, "Converting to backoff\n");
    ng.convert_to_backoff();
  }

  fprintf(stderr, "Writing\n");
  ng.write(out.file, true);
  out.close();
}

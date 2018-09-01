// This program converts the extended ARPA formats used by some tools
// to regular standard ARPA language models.
#include "TreeGram.hh"
#include "conf.hh"
#include "io.hh"
//#include "HashGram.hh"

int main(int argc, char **argv) {
  conf::Config config;
  config("Usage:  arpa2arpa nonstandard_in standard_out\nConverts interpolated "
         "arpa to backoff arpa.\n")('t', "tabs", "", "",
                                    "Use tabs instead of space between the "
                                    "fields when writing the ARPA file.");
  config.parse(argc, argv, 2, true);

  std::string field_separator = config["tabs"].specified ? "\t" : " ";

  io::Stream::verbose = true;
  io::Stream in(config.arguments[0], "r");
  io::Stream out(config.arguments[1], "w");

  TreeGram ng;
  // HashGram_t<int> ng;

  fprintf(stderr, "Reading\n");
  ng.read(in.file);
  in.close();

  fprintf(stderr, "Writing\n");
  ng.write(out.file, false, field_separator);
  out.close();
}

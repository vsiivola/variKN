// This program adds n-grams to the model so that it can be written as a
// full prefix tree ARPA model. This is required for compability with some
// other tools. Also, the binary format used requires this.
#include <memory>
#include "HashGram.hh"
#include "conf.hh"
#include "io.hh"

int main(int argc, char *argv[]) {
  conf::Config config;
  config(
      "Usage: add_zeroprob_grams arpain arpaout\nAdds grams for treegram.\n")(
      's', "smallvocab", "", "",
      "Vocabulary is less than 65000 entries. Saves some memory.");
  config.parse(argc, argv, 2);

  const bool smallvocab = config["smallvocab"].specified;

  io::Stream::verbose = true;
  io::Stream in(config.arguments[0], "r");
  io::Stream out(config.arguments[1], "w");

  std::unique_ptr<HashGram> hg(
      smallvocab ? std::unique_ptr<HashGram>(new HashGram_t<unsigned short>())
                 : std::unique_ptr<HashGram>(new HashGram_t<int>()));
  fprintf(stderr, "Reading\n");
  hg->read(in.file);
  in.close();

  fprintf(stderr, "Adding zpgs\n");
  hg->add_zeroprob_grams();

  fprintf(stderr, "Writing\n");
  hg->write(out.file);
  out.close();
}

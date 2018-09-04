// Program to grow an n-gram model
#include <memory>
#include "VarigramFuncs.hh"
#include "conf.hh"
#include "io.hh"

int main(int argc, char **argv) {
  conf::Config config;
  config("Usage: varigram_kn [OPTIONS] textin LM_out\nProduces variable span "
         "n-gram from given text\n")(
      'o', "opti=FILE", "arg", "",
      "The devel set for optimizing discount parameters. If not set, use "
      "leave-one-out discount estimates.")(
      'n', "norder=INT", "arg", "0",
      "Maximal order included in the model (default unrestricted)")(
      'D', "dscale=FLOAT", "arg", "-1.0", "Model size scale factor")(
      'E', "dscale2=FLOAT", "arg", "0",
      "Model size scaling during pruning step (default no pruning=0)")(
      'a', "arpa", "", "", "Output arpa instead of binary LM")(
      'x', "narpa", "", "",
      "Output nonstandard interpolated arpa instead of binary LM (uses less "
      "memory)")('3', "3nzer", "", "",
                 "Use 3 discounts per order instead of one. Takes a bit more "
                 "memory during the model estimation.")(
      'f', "nfirst=INT", "arg", "-1",
      "Number of most common words to be included")(
      'd', "ndrop=INT", "arg", "0",
      "Drop all words with less than ndrop occurances. If both nfirst and "
      "ndrop options are specified, the tighter bound is taken")(
      'B', "vocabin=FILE", "arg", "", "Restrict the vocabulary to given words")

      ('H', "hashsize=INT", "arg", "0",
       "Size of the reserved hash table. Speed vs. memory consumption.")(
          's', "smallvocab", "", "",
          "Vocabulary is less than 65000 entries. Saves some memory.")(
          'S', "smallmem", "", "",
          "Do not cache the data. Saves some memory, but slows down.")(
          'C', "clear_history", "", "",
          "Clear LM history on each start of sentence tag (<s>).")
      //('X', "debug_counts=FILE", "arg", "", "Write debug counts.")
      ('Z', "dontaddzpg", "", "",
       "This flag prevents the addition of zeroprob grams to fill out the "
       "treegram structure.")
      //('l', "clusterfile=FILE", "arg", "", "Uses clustering specified in
      // clusterfile.")
      ('A', "absolute", "", "",
       "Use absolute discounting. Several options won't work with this.")
      //('i', "iter=INT", "arg", "1", "Iterate growing and pruning iter times.")
      ('W', "write_counts=FILE", "arg", "",
       "Write resulting count matrices to FILE.")(
          'U', "write_vocab=FILE", "arg", "",
          "Write resulting vocabulary FILE.")(
          'O', "cutoffs=\"val1 val2 ... valN\"", "arg", "",
          "Use the specified cutoffs. The last value is used for all higher "
          "order n-grams.")('N', "discard_unks", "", "",
                            "Remove n-grams containing OOV words.")(
          'L', "longint", "", "",
          "Store counts in a long int type. Needed for big training sets.")(
          'V', "numngramstarget=INT", "arg", "0",
          "Scale model down until there are less than V*1.03 ngrams in the "
          "model")('F', "forcedisc=FLOAT", "arg", "-1.0",
                   "Set all discounts to the given value.");

  config.parse(argc, argv, 2, true);

  const int nfirst = config["nfirst"].get_int();
  const int max_order = config["norder"].get_int();
  const int ndrop = config["ndrop"].get_int();
  const indextype hashs = config["hashsize"].get_int();
  const bool narpa = config["narpa"].specified;
  const bool arpa = config["arpa"].specified;
  const bool use_3nzer = config["3nzer"].specified;
  const std::string optiname(config["opti"].get_str());
  const float dscale = std::max(0.00001, config["dscale"].get_double());
  const float dscale2 = config["dscale2"].get_double();
  const int ngram_prune_target = config["numngramstarget"].get_int();
  const bool smallvocab = config["smallvocab"].specified;
  const bool smallmem = config["smallmem"].specified;
  const bool zpg = config["dontaddzpg"].specified;
  const int iter = 1; // config["iter"].get_int();
  const bool absolute = config["absolute"].specified;
  const std::string countsout(config["write_counts"].get_str());
  const std::string vocabout(config["write_vocab"].get_str());
  const bool discard_unks = config["discard_unks"].specified;
  const bool longint = config["longint"].specified;
  const float force_disc = config["forcedisc"].get_double();

  std::string vocabname;
  if (config["vocabin"].specified)
    vocabname = config["vocabin"].get_str();

  bool ok = true;
  std::vector<int> cutoffs =
      str::long_vec<int>(config["cutoffs"].get_str(), &ok);
  if (!ok) {
    fprintf(stderr, "Error parsing cutoffs, exit\n");
    exit(-1);
  }

  io::Stream vocab_outf, counts_outf;
  if (vocabout.size())
    vocab_outf.open(vocabout, "w");
  if (countsout.size())
    counts_outf.open(countsout, "w");

  std::string infilename(config.arguments.at(0));
  io::Stream::verbose = true;
  if (infilename == "-") {
    fprintf(stderr, "varigram_kn might need to scan the input several times.  "
                    "Cannot read the main data from stdin \"-\". Exit.\n");
    exit(-1);
  }
  io::Stream out(config.arguments.at(1), "w");
  io::Stream::verbose = false;

  std::unique_ptr<Varigram> vg;
  if (!smallvocab)
    if (!longint)
      vg.reset(new Varigram_t<int, int>(use_3nzer, absolute));
    else
      vg.reset(new Varigram_t<int, long>(use_3nzer, absolute));
  else if (!longint)
    vg.reset(new Varigram_t<unsigned short, int>(use_3nzer, absolute));
  else
    vg.reset(new Varigram_t<unsigned short, long>(use_3nzer, absolute));

  if (dscale > 0.0) {
    vg->set_datacost_scale(dscale);
  }
  if (dscale2 > 0.0) {
    vg->set_datacost_scale2(dscale2); // use also pruning
  }
  if (ngram_prune_target > 0) {
    vg->set_ngram_prune_target(ngram_prune_target);
  }
  if (max_order != 0) {
    vg->set_max_order(max_order);
  }

  try {
    if (config["clear_history"].specified)
      vg->initialize(infilename, hashs, ndrop, nfirst, optiname, "<s>",
                     smallmem, vocabname);
    else
      vg->initialize(infilename, hashs, ndrop, nfirst, optiname, "", smallmem,
                     vocabname);

    if (cutoffs.size())
      vg->set_cutoffs(cutoffs);
    vg->set_discard_unks(discard_unks);
    vg->grow(iter);

    if (force_disc >= 0.0f) {
      fprintf(stderr, "Forcing all discounts to %g.\n", force_disc);
      vg->set_all_discounts(force_disc);
    }

    if (vocab_outf.file) {
      vg->write_vocab(vocab_outf.file);
      vocab_outf.close();
    }

    if (counts_outf.file) {
      vg->write_debug_counts(counts_outf.file);
      counts_outf.close();
    }

    if (!narpa) {
      vg->write(out.file, arpa);
    } else {
      if (zpg) {
        vg->set_zeroprobgrams(false);
      }
      vg->write_narpa(out.file);
    }
    out.close();
  } catch (std::exception &e) {
    fprintf(stderr, "%s\n", e.what());
    exit(1);
  }
}

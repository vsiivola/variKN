// This program creates first creates a full  unpruned language model of
// given order. The model can be then pruned, if some options are given.
#include <memory>
#include "InterKn.hh"
#include "conf.hh"
#include "io.hh"

int main(int argc, char **argv) {
  conf::Config config;
  config("Usage: counts2kn [OPTIONS] text_in lm_out.\nCreates an interpolated "
         "KN language model based on the give counts\n")(
      'o', "opti=FILE", "arg", "",
      "The devel set for optimizing discount parameters. If not set, use "
      "leave-one-out discount estimates.")(
      'n', "norder=INT", "arg", "0", "Order of the counts in the count file")(
      'a', "arpa", "", "", "Output arpa instead of binary LM")(
      'x', "narpa", "", "",
      "Output nonstandard interpolated arpa instead of binary LM (uses less "
      "memory)")('3', "3nzer", "", "",
                 "Use 3 discounts per order instead of one. Takes a bit more "
                 "memory during model estimation.")(
      'B', "vocabin=FILE", "arg", "",
      "Input count file contains word indices. This is the corresponding "
      "vocabulary")(
      'H', "hashsize=INT", "arg", "0",
      "Size of the reserved hash table. Speed vs. memory consumption.")(
      'c', "counts", "", "",
      "The first arguments is a counts file with already KN-smoothed counts, "
      "not a text file to be read. Disables options -f and -d")(
      'Z', "h_counts", "", "",
      "The first arguments is a counts file with only highest order counts, "
      "not a text file to be read. Disables options -f and -d")(
      'f', "nfirst=INT", "arg", "99999999",
      "Number of most common words to be included")(
      'd', "ndrop=INT", "arg", "0",
      "Drop all words with less than ndrop occurances. If both nfirst and "
      "ndrop options are specified, the tighter bound is taken")(
      's', "smallvocab", "", "",
      "Vocabulary is less than 65000 entries. Saves some memory.")(
      'p', "prunetreshold", "arg", "-1",
      "Prune out the n-grams, for which the score does not exceed the "
      "treshold. Default: no pruning=0.")(
      'A', "absolute", "", "",
      "Use absolute discounting. Several options won't work with this.")(
      'C', "clear_history", "", "",
      "Clear LM history on each start of sentence tag (<s>).")
      //('e', "ehist", "", "", "Debug ehist pruning.")
      //('l', "clusterfile=FILE", "arg", "", "Uses clustering specified in
      // clusterfile.")
      ('r', "prune_with_real_counts", "", "",
       "Use real counts and not type counts when deciding which grams to use. "
       "Slightly better results, uses more memory.")(
          'R', "real_counts_data=FILE", "arg", "",
          "if -r is specified with -c or -Z, the original training data must "
          "be specified.")('W', "write_counts=FILE", "arg", "",
                           "Write resulting count matrices to FILE.")(
          'U', "write_vocab=FILE", "arg", "",
          "Write resulting vocabulary FILE.")
      //('u', "skip_unks", "arg", "", "Do not estimate any probability for
      // unknow tokens")
      ('O', "cutoffs=\"val1 val2 ... valN\"", "arg", "",
       "Use the specified cutoffs. The last value is used for all higher order "
       "n-grams.")(
          'D', "discard_cutoffs", "", "",
          "Simply throw away cutoff grams, do not adjust lower order counts.")(
          'N', "discard_unks", "", "", "Remove n-grams containing OOV words.")(
          'L', "longint", "", "",
          "Store counts in a long int type. Needed for big training sets.");
  config.parse(argc, argv, 2, true);

  std::string optiname(config["opti"].get_str());
  const int n = config["norder"].get_int();
  const indextype hashs = config["hashsize"].get_int();
  const bool narpa = config["narpa"].specified;
  const bool arpa = config["arpa"].specified;
  const bool use_3nzer = config["3nzer"].specified;
  const bool absolute = config["absolute"].specified;
  const int ndrop = config["ndrop"].get_int();
  const int nfirst = config["nfirst"].get_int();
  const bool smallvocab = config["smallvocab"].specified;
  const float prunetreshold = config["prunetreshold"].get_double();
  // const bool ehist=config["ehist"].specified;
  const bool prune_with_real_counts =
      config["prune_with_real_counts"].specified;
  const bool discard_cutoffs = config["discard_cutoffs"].specified;
  const bool discard_unks = config["discard_unks"].specified;
  const bool longint = config["longint"].specified;
  const std::string countsout(config["write_counts"].get_str());
  const std::string vocabout(config["write_vocab"].get_str());
  const std::string rcfile(config["real_counts_data"].get_str());
  // const bool skip_unks = config["skip_unks"].specified;
  // assert(!skip_unks); // FIXME: implement this
  bool ok = true;
  std::vector<int> cutoffs =
      str::long_vec<int>(config["cutoffs"].get_str(), &ok);
  if (!ok) {
    fprintf(stderr, "Error parsing cutoffs, exit\n");
    exit(-1);
  }

  int counts_in = 0;
  if (config["counts"].specified)
    counts_in = 1;
  if (config["h_counts"].specified) {
    if (counts_in == 1) {
      fprintf(stderr, "Cannot specify both -counts and -h_counts. Exit.\n");
      exit(-1);
    }
    counts_in = -1;
  }
  if (counts_in != 1 && n == 0) {
    fprintf(stderr, "Must specify order -n.\nExit.\n");
    exit(-1);
  }

  io::Stream::verbose = true;
  std::string dataname = config.arguments.at(0);
  if (dataname == "-") {
    fprintf(stderr, "counts2kn might need to scan the input several times.  "
                    "Cannot read the main data from stdin \"-\". Exit.\n");
    exit(-1);
  }
  io::Stream out(config.arguments.at(1), "w");
  std::string vocabname = config["vocabin"].get_str();

  std::string prunedata_name;
  if (counts_in && prune_with_real_counts && rcfile.size() == 0) {
    fprintf(stderr, "Using counts from file: pruning with real counts needs -R "
                    "traindata. Exit.\n");
    exit(-1);
  }
  if (rcfile.size())
    prunedata_name = rcfile;
  else
    prunedata_name = dataname;
  TreeGram lm;

  std::string ss_sym;
  if (config["clear_history"].specified) {
    fprintf(stderr, "clearing history\n");
    ss_sym = "<s>";
  }

  try {
    /* Construct the base kn-smoothed model */
    fprintf(stderr, "Estimating counts\n");
    std::unique_ptr<InterKn> kn;
    bool init_disc = true;

    /* Parse the arguments, create the right kind of model*/
    if (!use_3nzer && !smallvocab) {
      if (!longint)
        kn.reset(new InterKn_int_disc<int, int>(
            absolute, dataname, vocabname, optiname, counts_in, n, ndrop,
            nfirst, nullptr, prunedata_name, ss_sym, hashs));
      else
        kn.reset(new InterKn_int_disc<int, long>(
            absolute, dataname, vocabname, optiname, counts_in, n, ndrop,
            nfirst, nullptr, prunedata_name, ss_sym, hashs));
    } else if (!use_3nzer && smallvocab) {
      if (!longint)
        kn.reset(new InterKn_int_disc<unsigned short, int>(
            absolute, dataname, vocabname, optiname, counts_in, n, ndrop,
            nfirst, nullptr, prunedata_name, ss_sym, hashs));
      else
        kn.reset(new InterKn_int_disc<unsigned short, long>(
            absolute, dataname, vocabname, optiname, counts_in, n, ndrop,
            nfirst, nullptr, prunedata_name, ss_sym, hashs));
    } else if (use_3nzer && !smallvocab) {
      if (!longint)
        kn.reset(new InterKn_int_disc3<int, int>(
            absolute, dataname, vocabname, optiname, counts_in, n, ndrop,
            nfirst, nullptr, prunedata_name, ss_sym, hashs));
      else
        kn.reset(new InterKn_int_disc3<int, long>(
            absolute, dataname, vocabname, optiname, counts_in, n, ndrop,
            nfirst, nullptr, prunedata_name, ss_sym, hashs));
    } else if (use_3nzer && smallvocab) {
      if (!longint)
        kn.reset(new InterKn_int_disc3<unsigned short, int>(
            absolute, dataname, vocabname, optiname, counts_in, n, ndrop,
            nfirst, nullptr, prunedata_name, ss_sym, hashs));
      else
        kn.reset(new InterKn_int_disc3<unsigned short, long>(
            absolute, dataname, vocabname, optiname, counts_in, n, ndrop,
            nfirst, nullptr, prunedata_name, ss_sym, hashs));
    }
    kn->prune_with_real_counts = prune_with_real_counts;
    fprintf(stderr, "The model will use ");
    if (smallvocab)
      fprintf(stderr, "small vocabulary (<65534) and ");
    else
      fprintf(stderr, "large vocabulary and ");
    if (use_3nzer)
      fprintf(stderr, "modified ");
    if (absolute)
      fprintf(stderr, "absolute discounting.\n");
    else
      fprintf(stderr, "Kneser-Ney smoothing.\n");

    if (init_disc) {
      kn->init_disc(0.71);
    }

    // if (ehist) kn->use_ehist_pruning(kn->input_data_size);
    kn->cutoffs = cutoffs;
    kn->discard_cutoffs = discard_cutoffs;
    kn->discard_ngrams_with_unk = discard_unks;
    kn->create_model(std::max((float)0.0, prunetreshold));

    fprintf(stderr, "Writing model\n");
    if (vocabout.size()) {
      io::Stream out(vocabout, "w");
      kn->vocab.write(out.file);
    }

    if (countsout.size()) {
      io::Stream out(countsout, "w");
      kn->write_counts(out.file);
    }

    if (!narpa) {
      kn->counts2lm(&lm);
      lm.write(out.file, !arpa);
    } else {
      kn->counts2asciilm(out.file);
    }
    out.close();
  } catch (std::exception &e) {
    fprintf(stderr, "%s\n", e.what());
    exit(1);
  }
}

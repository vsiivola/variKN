// Calculate the perplexity of the test corpus given a language model
#include "InterTreeGram.hh"
#include "PerplexityFuncs.hh"
#include "conf.hh"
#include <cmath>

int main(int argc, char *argv[]) {
  conf::Config config;
  config("Usage: perplexity [OPTIONS] text_in results_out\nCounts the "
         "perplexity of the test text relative to given model\n")(
      'a', "arpa=FILE", "arg", "", "Arpa format language model")(
      'A', "bin=FILE", "arg", "",
      "Binary format LM. Either arpa or binary model must be specified")(
      'C', "ccs=FILE", "arg", "", "Context cue file. One token per line.")(
      'W', "wb=FILE", "arg", "", "Word break symbol file. One token per line.")(
      'X', "mb=FILE", "arg", "",
      "Morph boundary prefix/postfix file. One token per line.")(
      'u', "unk=STRING", "arg", "",
      "Unk symbol (defaul <UNK>, case sensitive)")(
      'U', "unkwarn", "", "", "Warn if unknown tokens are seen")(
      'o', "includeunks", "", "",
      "Include unknown tokens in perplexity calculations.")(
      'i', "interpolate=FILE", "arg", "", "Interpolate with given arpa LM.")(
      'I', "inter_coeff=FLOAT", "arg", "-1",
      "Interpolation coefficient. The interpolated model will be weighted by "
      "coeff whereas the main model will be weighted by 1.0-coeff.")(
      't', "init_hist=INT", "arg", "",
      "Take n first tokens after \"</s>\" as initial LM history (no "
      "probabilities assigned)")('f', "freegram", "", "",
                                 "No prefix requirements on the n-gram model.")(
      's', "smallvocab", "", "",
      "Vocabulary is less than 65000 entries. Saves some memory.")(
      'S', "probstream=FILE", "arg", "",
      "Write probability of each symbol to FILE")
      //('N',"stream_interval=INT","arg","","Sum every N tokens for each stream
      // output")
      ;
  config.parse(argc, argv, 2, true);
  int lm_type = 0;

  std::string lm_name;
  if (config["arpa"].specified)
    lm_name = config["arpa"].get_str();
  else {
    if (!config["bin"].specified) {
      fprintf(stderr, "Either -arpa or -bin needs to be specified. Exit.\n");
      exit(-1);
    }
    if (config["interpolate"].specified) {
      fprintf(
          stderr,
          "Interpolation of binary models is not implemented (yet?). Exit.\n");
      exit(-1);
    }
    lm_name = config["bin"].get_str();
    lm_type = 1;
  }

  std::string ccs_name, wb_name, mb_name, unk_symbol;
  if (config["ccs"].specified)
    ccs_name = config["ccs"].get_str();
  if (config["wb"].specified)
    wb_name = config["wb"].get_str();
  if (config["mb"].specified)
    mb_name = config["mb"].get_str();
  if (config["unk"].specified)
    unk_symbol = config["unk"].get_str();
  io::Stream stream_out;
  if (config["probstream"].specified)
    stream_out.open(config["probstream"].get_str(), "w");
  bool unkwarn = config["unkwarn"].specified;
  bool skip_unks = !config["includeunks"].specified;
  int freegram = 0;
  if (config["freegram"].specified) {
    if (config["smallvocab"].specified)
      freegram = -1;
    else
      freegram = 1;
  }

  int stream_interval = 1;
  int init_hist = 0;
  if (config["init_hist"].specified) {
    if (config["ccs"].specified) {
      fprintf(stderr, "Can't specify both --init_hist and --ccs. Exit\n");
      exit(-1);
    }
    init_hist = config["init_hist"].get_int();
  }
  io::Stream::verbose = true;
  io::Stream txtin(config.arguments[0], "r");
  io::Stream out(config.arguments[1], "w");
  io::Stream::verbose = false;

  if (wb_name.length() && mb_name.length()) {
    fprintf(stderr, "Cannot specify both -mb and -wb. Remove one. Exit\n");
    exit(-1);
  }

  fprintf(stderr, "lm %s", lm_name.c_str());
  if (ccs_name.length())
    fprintf(stderr, ", ccs %s", ccs_name.c_str());
  if (wb_name.length())
    fprintf(stderr, ", wb_name %s", wb_name.c_str());
  if (mb_name.length())
    fprintf(stderr, ", mb_name %s", mb_name.c_str());
  if (unk_symbol.length())
    fprintf(stderr, ", unk_symbol %s", unk_symbol.c_str());
  fprintf(stderr, "\n");

  std::unique_ptr<Perplexity> lm;
  bool knows_hitrates = true;
  if (config["interpolate"].specified) {
    // Interpolate the old way, using PerplexityFuncs
    if (config["freegram"].specified) {
      lm.reset(new Perplexity(lm_name, lm_type, ccs_name, wb_name, mb_name,
                              unk_symbol, true, skip_unks));
      lm->set_interpolation(config["interpolate"].get_str());
      if (config["inter_coeff"].specified) {
        lm->set_alpha(config["inter_coeff"].get_double());
      }
    } else {
      // Interpolate the new way, using InterTreeGram
      knows_hitrates = false;
      std::vector<std::string> lm_names;
      lm_names.push_back(lm_name);
      lm_names.push_back(config["interpolate"].get_str());

      std::vector<float> coeffs;
      float coeff = 0.5;
      if (config["inter_coeff"].specified) {
        coeff = config["inter_coeff"].get_double();
      }
      coeffs.push_back(1.0 - coeff);
      coeffs.push_back(coeff);
      std::shared_ptr<NGram> itg(new InterTreeGram(lm_names, coeffs));
      lm.reset(new Perplexity(itg, ccs_name, wb_name, mb_name, unk_symbol,
                              skip_unks));
    }
  } else {
    lm.reset(new Perplexity(lm_name, lm_type, ccs_name, wb_name, mb_name,
                            unk_symbol, freegram, skip_unks));
    if (config["inter_coeff"].specified) {
      fprintf(stderr,
              "Only on lm specified, cannot set interpolation coff. Exit\n");
      exit(-1);
    }
  }
  lm->set_unk_warn(unkwarn);

  if (init_hist != 0) {
    lm->set_init_hist(init_hist);
  }
  if (stream_out.file) {
    lm->logprob_file(txtin.file, stream_out.file, stream_interval);
  } else {
    lm->logprob_file(txtin.file, nullptr, 1);
  }
  lm->print_results(out.file);
  if (knows_hitrates) // FIXME: This doesn't work yet for InterTreeGrams
    lm->print_hitrates(out.file);
}

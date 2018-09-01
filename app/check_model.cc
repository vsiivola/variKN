/* Check that the lm is ok, that is for given ngram,
   probabilities sum to 1.0 */

#include "HashGram.hh"
#include "TreeGram.hh"
#include "conf.hh"
#include "def.hh"
#include "io.hh"
#include <cmath>

NGram *tg;

float check_gram(const std::deque<int> &g, bool verbose) {
  std::deque<int> gram(g);
  gram.push_back(-1);
  double prob = 0.0;
  double max_prob = 0.0, max_prob2 = 0.0;
  double tmp_prob;
  int max_idx = -1, max_idx2 = -1;
  for (int i = 0; i < tg->num_words(); i++) {
    gram.back() = i;
    tmp_prob = pow(10, tg->log_prob(gram));
    if (verbose) {
      fprintf(stderr, "%s: %d -> %g (order %d) ", tg->word(i).c_str(), i,
              tmp_prob, tg->last_order());
      // print_indices(gram);
      fprintf(stderr, "\n");
    }
    if (tmp_prob > max_prob) {
      max_prob2 = max_prob;
      max_idx2 = max_idx;
      max_prob = tmp_prob;
      max_idx = i;
    } else if (tmp_prob > max_prob2) {
      max_prob2 = tmp_prob;
      max_idx2 = i;
    }
    prob += tmp_prob;
  }
  if (verbose) {
    fprintf(stderr, "Sum %lg\n", prob);
    fprintf(stderr, "Most probable %s, %lg\n", tg->word(max_idx).c_str(),
            max_prob);
    fprintf(stderr, "Second most probable %s, %lg\n",
            tg->word(max_idx2).c_str(), max_prob2);
  }
  return (prob);
}

int main(int argc, char **argv) {
  conf::Config config;
  config(
      "Use: check_model model_in text_in\nChecks that everything sums to 1.\n")(
      'a', "arpa", "", "", "Arpa model instead of binary model")(
      'd', "datafile", "", "",
      "second argument is a filename, sum ovaer all grams in the file")(
      'f', "freegram", "", "", "Freeprefix arpa model instead of binary model");

  config.parse(argc, argv, 2, true);

  io::Stream::verbose = true;
  io::Stream in(config.arguments[0], "r");
  io::Stream txtin(config.arguments[1], "r");
  bool binary = !config["arpa"].specified;
  bool freegram = config["freegram"].specified;
  bool datafile = config["datafile"].specified;

  if (freegram)
    tg = new HashGram_t<int>;
  else
    tg = new TreeGram;

  tg->read(in.file, binary);
  in.close();

  if (!datafile) {
    std::deque<int> gram;
    char tmpbuf[MAX_WLEN];
    while (fscanf(txtin.file, "%s", tmpbuf) != EOF)
      gram.push_back(tg->word_index(tmpbuf));
    check_gram(gram, true);
    txtin.close();
    exit(0);
  }

  char word[MAX_WLEN];
  std::deque<int> history;
  int errors = 0;
  int num_words = 0;
  while (true) {
    const int fsc = fscanf(txtin.file, "%s", word);
    if (!fsc || fsc == EOF)
      break;
    num_words++;
    if (history.size() == tg->order())
      history.pop_front();
    history.push_back(tg->word_index(word));
    /*
    fprintf(stderr,"Checking [");
    for (int i=0;i<history.size();i++)
      fprintf(stderr," %s", tg->word(history[i]).c_str());
    fprintf(stderr," ]\n");
    */
    const float psum = check_gram(history, false);
    if (psum > 1.001 || psum < 0.99) {
      fprintf(stderr, "Sum failed for ");
      for (int i = 0; i < history.size(); i++) {
        fprintf(stderr, "%s ", tg->word(history[i]).c_str());
      }
      fprintf(stderr, "= %g\n", psum);
      errors++;
    }
  }
  fprintf(stderr, "Done, %d/%d errors found.\n", errors, num_words);
  txtin.close();
  exit(0);
}

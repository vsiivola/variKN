// This program prints out the number of n-grams in the language model.
#include "TreeGram.hh"
#include "conf.hh"
#include "io.hh"
#include "str.hh"
#include <cassert>
#include <cmath>
#include <cstdio>

const float SIMILIMI = 1e-6;
const float PERLIM = 0.005;

void simple_count(const conf::Config &config, io::Stream *in, io::Stream *out) {
  /**********************************************************************/
  // The simple vesion, just read the counts of the arpa file
  /***********************************************************************/
  if (!config["arpa"].specified) {
    fprintf(stderr, "Sorry, the --simple flag only works with arpa format "
                    "models currently. Fix me.\n");
    exit(-1);
  }

  std::string sbuf;
  while (str::read_line(&sbuf, in->file, true) && "\\data\\" != sbuf)
    ;
  int i = 1;
  int i2;
  int count;
  std::vector<int> counts;
  while (str::read_line(&sbuf, in->file, true) &&
         sscanf(sbuf.c_str(), "ngram %d=%d", &i2, &count) == 2) {
    if (i2 != i) {
      fprintf(stderr, "read error, exit\n");
      exit(-1);
    }
    counts.push_back(count);
    i++;
  }
  if (sbuf.length() > 0) {
    fprintf(stderr, "read error,exit\n");
    exit(-1);
  }
  in->close();

  int total = 0;
  for (int i = 0; i < counts.size(); i++) {
    fprintf(out->file, "%d-grams: %d\n", i + 1, counts[i]);
    total += counts[i];
  }
  fprintf(out->file, "--\ntotal %d n-grams\n", total);
  out->close();
}

void dummyless_count(TreeGram *ng, const io::Stream *out) {
  int total_grams = ng->gram_count(1), order_grams;
  fprintf(out->file, "1-grams: %d\n", total_grams);

  TreeGram::Gram indices, prefixindices;
  TreeGram::Iterator iter(nullptr);
  for (int o = 2; o <= ng->order(); o++) {
    order_grams = 0;
    indices.resize(o);
    prefixindices.resize(o - 1);
    iter.reset(ng);
    while (iter.next_order(o)) {
      indices[0] = iter.node(1).word;
      for (int j = 2; j <= o; j++) {
        indices[j - 1] = iter.node(j).word;
        prefixindices[j - 2] = indices[j - 1];
      }

      if (iter.node().back_off < -SIMILIMI) {
        // fprintf(stderr,"accepted, bo %g\n",iter.node().back_off);
        order_grams++;
        continue;
      }
      const float lp_orig = ng->log_prob(indices);
      const float lp_prefix = ng->log_prob(prefixindices);
      // fprintf(stderr,"lp %g, lpp %g, pdif %2.1f, ",lp_orig, lp_prefix,
      // fabs((lp_orig-lp_prefix)/lp_orig));
      if (fabs((lp_orig - lp_prefix) / lp_orig) < PERLIM) {
        // fprintf(stderr,"rejected\n");
        continue;
      }
      // fprintf(stderr,"accepted\n");
      order_grams++;
    }
    fprintf(out->file, "%d-grams: %d\n", o, order_grams);
    total_grams += order_grams;
  }
  fprintf(out->file, "Total grams: %d\n", total_grams);
}

int main(int argc, char **argv) {
  conf::Config config;
  config("Usage:  arpasize in.bin out\nShows the number of grams modeled in an "
         "arpa file. Grams added for bookkeeping are discarded from these "
         "counts.\n")(
      'S', "simple", "", "",
      "Simply use the reported counts and dont's check for dummy grams.")(
      'a', "arpa", "", "", "language model is in arpa format");
  config.parse(argc, argv, 2, true);

  io::Stream::verbose = true;
  io::Stream in(config.arguments[0], "r");
  io::Stream out(config.arguments[1], "w");

  if (!config["simple"].specified) {
    TreeGram ng;
    fprintf(stderr, "Reading model\n");
    ng.read(in.file, !config["arpa"].specified);
    in.close();

    dummyless_count(&ng, &out);
    out.close();
    return (0);
  }
  simple_count(config, &in, &out);
  return (0);
}

#define BOOST_TEST_MODULE varikn_unit_tests
#include <boost/test/included/unit_test.hpp>

#include <cstdio>

#include <HashGram.hh>
#include <InterKn.hh>
#include <InterTreeGram.hh>
#include <PerplexityFuncs.hh>
#include <TreeGram.hh>
#include <VarigramFuncs.hh>

void create_lm(std::string dataname, std::string vocabname,
               std::string optiname, int n, std::string outfname) {
  InterKn_int_disc3<unsigned short, int> kn(false, dataname, vocabname,
                                            optiname, 0, n, 0, 99999999, NULL,
                                            dataname, "", 100);
  kn.init_disc(0.71);
  kn.create_model(0.0);
  io::Stream out(outfname, "w");
  kn.counts2asciilm(out.file);
  out.close();
}

void create_varigram_lm(std::string dataname, std::string vocabname,
                        std::string optiname, int n, std::string outfname,
                        float dscale1 = 1.0, float dscale2 = 2.0) {
  Varigram_t<unsigned short, int> vg(false, false);
  vg.set_datacost_scale(dscale1);
  vg.set_datacost_scale2(dscale2);
  vg.initialize(dataname, 100, 0, 9999999, optiname, "", false, vocabname);
  vg.grow(1);
  io::Stream out(outfname, "w");
  vg.write(out.file, true);
  out.close();
}

float itg_perplexity(std::string model1, std::string model2,
                     std::string infname, bool use_hashgram = false,
                     float alpha = 0.5) {
  std::vector<std::string> lm_names;
  lm_names.push_back(model1);
  lm_names.push_back(model2);

  std::vector<float> coeffs;
  coeffs.push_back(alpha);
  coeffs.push_back(1.0 - alpha);

  auto itg = std::make_shared<InterTreeGram>(lm_names, coeffs);
  io::Stream txtin(infname, "r");
  io::Stream out("-", "w");

  Perplexity lm(itg, "", "", "", "");
  lm.logprob_file(txtin.file, NULL);
  return lm.print_results(out.file);
}

float perplexity(std::string model1, std::string model2, std::string infname,
                 bool use_hashgram = false, float alpha = 0.5) {
  io::Stream txtin(infname, "r");
  io::Stream out("-", "w");

  Perplexity lm(model1, 0, "", "", "", "", false);
  if (model2 != "") {
    lm.set_interpolation(model2);
    lm.set_alpha(alpha);
  }
  lm.logprob_file(txtin.file, NULL);
  return lm.print_results(out.file);
}

void create_simple_models(std::string datadir) {
  float h;

  create_lm(datadir + "/ax20.txt", datadir + "/vocab.txt",
            datadir + "/ax20.txt", 3, datadir + "/a.arpa");
  fprintf(stdout, "Perplexity of a against a:\n");
  h = perplexity(datadir + "/a.arpa", "", datadir + "/ax20.txt", false);
  fprintf(stderr, "h1 %f\n", h);
  BOOST_REQUIRE(h < 0.01);
  h = perplexity(datadir + "/a.arpa", "", datadir + "/ax20.txt", true);
  fprintf(stderr, "h1.1 %f\n", h);
  BOOST_REQUIRE(h < 0.01);

  fprintf(stdout, "Perplexity of b against a:\n");
  create_varigram_lm(datadir + "/bx20.txt", datadir + "/vocab.txt",
                     datadir + "/bx20.txt", 3, datadir + "/b.arpa");
  h = perplexity(datadir + "/b.arpa", "", datadir + "/ax20.txt");
  fprintf(stderr, "h2 %f\n", h);
  BOOST_REQUIRE(h > 10.0);

  fprintf(stdout, "Perplexity of ab against ab:\n");
  create_lm(datadir + "/abx20.txt", "", datadir + "/abx20.txt", 3,
            datadir + "/ab.arpa");
  h = perplexity(datadir + "/ab.arpa", "", datadir + "/abx20.txt");
  fprintf(stderr, "h2.1 %f\n", h);
  BOOST_REQUIRE(h < 0.01);
  create_varigram_lm(datadir + "/abx20.txt", "", datadir + "/abx20.txt", 3,
                     datadir + "/ab-vg.arpa", 0.01, 0.02);
  h = perplexity(datadir + "/ab-vg.arpa", "", datadir + "/abx20.txt");
  fprintf(stderr, "h2.2 %f\n", h);
  BOOST_REQUIRE(h < 0.05);
}

void create_models_without_opti(std::string datadir) {
  float h;
  create_lm(datadir + "/ax20.txt", datadir + "/vocab.txt", "", 3,
            datadir + "/a-noopti.arpa");
  fprintf(stdout, "Perplexity of a against a:\n");
  h = perplexity(datadir + "/a-noopti.arpa", "", datadir + "/ax20.txt", false);
  fprintf(stderr, "h1-no_opti %f\n", h);
  fprintf(stdout, "Perplexity of b against a:\n");

  create_varigram_lm(datadir + "/bx20.txt", datadir + "/vocab.txt", "", 3,
                     datadir + "/b-noopti.arpa");
  h = perplexity(datadir + "/b-noopti.arpa", "", datadir + "/ax20.txt");
  fprintf(stderr, "h2-no_opti %f\n", h);
}

void test_interpolation(std::string datadir) {
  float h;

  fprintf(stdout, "Perplexity of b against b:\n");
  h = perplexity(datadir + "/b.arpa", "", datadir + "/bx20.txt");
  fprintf(stderr, "h3 %f\n", h);
  BOOST_REQUIRE(h < 0.01);

  fprintf(stdout, "Perplexity of ab interp 0.0:\n");
  h = perplexity(datadir + "/a.arpa", datadir + "/b.arpa",
                 datadir + "/abx20.txt", false, 0.01);
  fprintf(stderr, "h4 %f\n", h);

  fprintf(stdout, "Perplexity of ab interp 0.5:\n");
  h = perplexity(datadir + "/a.arpa", datadir + "/b.arpa",
                 datadir + "/abx20.txt", false, 0.5);
  float h2 = itg_perplexity(datadir + "/a.arpa", datadir + "/b.arpa",
                            datadir + "/abx20.txt", false, 0.5);
  fprintf(stderr, "h5 %f %f\n", h, fabs(h - 1.0));
  BOOST_REQUIRE(fabs(h - 1.0) < 0.01);
  BOOST_REQUIRE(fabs(h - h2) < 0.01);

  fprintf(stdout, "Perplexity of ab interp 1.0:\n");
  h = perplexity(datadir + "/a.arpa", datadir + "/b.arpa",
                 datadir + "/abx20.txt", false, 0.99);
  fprintf(stderr, "h6 %f\n", h);
  // BOOST_REQUIRE( abs(lp) < 0.01 );

  fprintf(stdout, "Perplexity of ab interp 0.0 against a:\n");
  h = perplexity(datadir + "/ab.arpa", datadir + "/b.arpa",
                 datadir + "/abx20.txt", false, 1.0);
  h2 = itg_perplexity(datadir + "/ab.arpa", datadir + "/b.arpa",
                      datadir + "/abx20.txt", false, 1.0);
  fprintf(stderr, "h7 %f\n", h);
  BOOST_REQUIRE(fabs(h) < 0.01);
  BOOST_REQUIRE(fabs(h - h2) < 0.01);
}

void test_interpolated_different_vocabs(std::string datadir) {
  create_lm(datadir + "/ax20.txt", "", datadir + "/ax20.txt", 3,
            datadir + "/a-novocab.arpa");
  create_varigram_lm(datadir + "/bx20.txt", "", datadir + "/bx20.txt", 3,
                     datadir + "/b-novocab.arpa");
  fprintf(stdout, "Perplexity of ab interp 0.5 (novocab):\n");
  float h = perplexity(datadir + "/a-novocab.arpa", datadir + "/b-novocab.arpa",
                       datadir + "/abx20.txt", false, 0.5);
  fprintf(stderr, "novocab %f\n", h);
  BOOST_REQUIRE(fabs(h - 1.0) < 0.01);
}

void test_intertreegram(std::string datadir) {
  std::vector<std::string> lm_names;
  lm_names.push_back(std::string(datadir + "/b-novocab.arpa"));
  lm_names.push_back(std::string(datadir + "/a-novocab.arpa"));

  std::vector<float> coeffs;
  coeffs.push_back(0.5);
  coeffs.push_back(0.5);
  // coeffs.push_back(1.0);

  InterTreeGram itg(lm_names, coeffs);
  itg.test_write("foo-itg1.arpa", 0);
  itg.test_write("foo-itg2.arpa", 1);
}

void test_hashinterpolate_write(std::string datadir) {
  fprintf(stderr, "hashinterpolate\n");
  std::vector<std::string> lm_names;
  HashGram_t<int> main_model, second_model;
  main_model.set_oov("<UNK>");
  std::string first_backoff_model(datadir + "/a-novocab-backoff.arpa");
  { // convert to backoff arpa
    TreeGram ng;
    io::Stream firstmodel_in(datadir + "/a-novocab.arpa", "r");
    io::Stream firstmodel_out(first_backoff_model, "w");
    ng.read(firstmodel_in.file);
    firstmodel_in.close();
    ng.write(firstmodel_out.file);
    firstmodel_out.close();
  }

  std::string second_backoff_model(datadir + "/b-novocab-backoff.arpa");
  { // convert to backoff arpa
    TreeGram ng;
    io::Stream secondmodel_in(datadir + "/b-novocab.arpa", "r");
    io::Stream secondmodel_out(second_backoff_model, "w");
    ng.read(secondmodel_in.file);
    secondmodel_in.close();
    ng.write(secondmodel_out.file);
    secondmodel_out.close();
  }
  io::Stream first_backoff_model_in(first_backoff_model, "r");
  main_model.read(first_backoff_model_in.file, false);
  first_backoff_model_in.close();
  io::Stream second_backoff_model_in(second_backoff_model, "r");
  second_model.set_oov("<UNK>");
  main_model.copy_vocab_to(second_model);
  second_model.read(second_backoff_model_in.file, false);
  second_model.copy_vocab_to(main_model);
  main_model.fake_interpolate(second_model, 0.5);
  std::string fakename(datadir + "/fakeinterpolate.arpa");
  io::Stream out(fakename, "w");
  main_model.write(out.file);
  out.close();

  float h = perplexity(fakename, "", datadir + "/abx20.txt");
  float h2 =
      perplexity(datadir + "/a-novocab.arpa", datadir + "/b-novocab.arpa",
                 datadir + "/abx20.txt", false, 0.5);
  fprintf(stderr, "fake interpolate %f vs real interpolate %f\n", h, h2);
  BOOST_REQUIRE(fabs(h - h2) < 0.01);

  // Check that thing sum to one
  {
    std::vector<int> gram0{1, 1, 0};
    std::vector<int> gram1{1, 1, 1};
    std::vector<int> gram2{1, 1, 2};
    float sum = pow(10, main_model.log_prob(gram0)) +
                pow(10, main_model.log_prob(gram1)) +
                pow(10, main_model.log_prob(gram2));
    fprintf(stderr, "Hashfakeinterpolate 3gram sum %f\n", sum);
    BOOST_REQUIRE(fabs(sum - 1.0) < 0.01);
  }

  // Check that thing sum to one
  {
    std::vector<int> gram0{1, 0};
    std::vector<int> gram1{1, 1};
    std::vector<int> gram2{1, 2};
    float sum = pow(10, main_model.log_prob(gram0)) +
                pow(10, main_model.log_prob(gram1)) +
                pow(10, main_model.log_prob(gram2));
    fprintf(stderr, "Hashfakeinterpolate 2gram sum %f\n", sum);
    BOOST_REQUIRE(fabs(sum - 1.0) < 0.01);
  }

  // Check that thing sum to one
  {
    std::vector<int> gram0{2, 0};
    std::vector<int> gram1{2, 1};
    std::vector<int> gram2{2, 2};
    float sum = pow(10, main_model.log_prob(gram0)) +
                pow(10, main_model.log_prob(gram1)) +
                pow(10, main_model.log_prob(gram2));
    fprintf(stderr, "Hashfakeinterpolate 2gramB sum %f\n", sum);
    BOOST_REQUIRE(fabs(sum - 1.0) < 0.01);
  }

  // Check that thing sum to one
  {
    std::vector<int> gram0{0};
    std::vector<int> gram1{1};
    std::vector<int> gram2{2};
    float sum = pow(10, main_model.log_prob(gram0)) +
                pow(10, main_model.log_prob(gram1)) +
                pow(10, main_model.log_prob(gram2));
    fprintf(stderr, "Hashfakeinterpolate 1gram sum %f\n", sum);
    BOOST_REQUIRE(fabs(sum - 1.0) < 0.01);
  }
}

BOOST_AUTO_TEST_CASE(test_create_simple_models) {
  std::string datadir(boost::unit_test::framework::master_test_suite().argv[1]);
  create_simple_models(datadir);
  create_models_without_opti(datadir);
  test_interpolation(datadir);
  test_interpolated_different_vocabs(datadir);
  test_intertreegram(datadir);
  test_hashinterpolate_write(datadir);
  // FIXME: Write unit test to read and write a arpa file with treegram and
  // hashgram
}

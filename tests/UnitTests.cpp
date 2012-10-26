#include <boost/test/minimal.hpp>
#include <cstdio>

#include <InterKn.hh>
#include <VarigramFuncs.hh>
#include <PerplexityFuncs.hh>
#include <InterTreeGram.hh>

void create_lm(std::string dataname, std::string vocabname, std::string optiname, int n,
               std::string outfname) {
  InterKn_int_disc3<unsigned short, int> kn(false, dataname, vocabname, optiname, 0, n, 0, 99999999,
                                            NULL, dataname, "", 100);
  kn.init_disc(0.71);
  kn.create_model(0.0);
  io::Stream out(outfname,"w");
  kn.counts2asciilm(out.file);
  out.close();
}

void create_varigram_lm(std::string dataname, std::string vocabname, std::string optiname, int n,
                        std::string outfname) {
  Varigram_t<unsigned short, int> vg(false, false);
  vg.set_datacost_scale(1);
  vg.set_datacost_scale2(2);
  vg.initialize(dataname, 100, 0, 9999999, optiname, "", false, vocabname);
  vg.grow(1);
  io::Stream out(outfname,"w");
  vg.write(out.file, true);
  out.close();
}

float perplexity(std::string model1, std::string model2, std::string infname, bool use_hashgram=false, float alpha=0.5) {
  io::Stream txtin(infname,"r");
  io::Stream out("-","w");
  
  Perplexity lm(model1, 0, "", "", "", false, false);
  if (model2 != "") {
    lm.set_interpolation(model2);
    lm.set_alpha(alpha);
  }
  lm.logprob_file(txtin.file, NULL);
  return lm.print_results(out.file);
}

void create_simple_models(std::string datadir) {
  float h;

  create_lm(datadir+"/ax20.txt", datadir+"/vocab.txt", datadir+"/ax20.txt", 3, datadir+"/a.arpa");
  fprintf(stdout, "Perplexity of a against a:\n");
  h = perplexity(datadir+"/a.arpa", "", datadir+"/ax20.txt", false);
  fprintf(stderr,"h1 %f\n", h );
  BOOST_REQUIRE( h > -0.01 );
  h = perplexity(datadir+"/a.arpa", "", datadir+"/ax20.txt", true);
  fprintf(stderr,"h1.1 %f\n", h );
  BOOST_REQUIRE( h > -0.01 );

  fprintf(stdout, "Perplexity of b against a:\n");
  create_varigram_lm(datadir+"/bx20.txt", datadir+"/vocab.txt", datadir+"/bx20.txt", 3, datadir+"/b.arpa");
  h = perplexity(datadir+"/b.arpa", "", datadir+"/ax20.txt");
  fprintf(stderr,"h2 %f\n", h );
  BOOST_REQUIRE( h > 10.0 );
}

void test_interpolation(std::string datadir) {
  float h;

  fprintf(stdout, "Perplexity of b against b:\n");
  h = perplexity(datadir+"/b.arpa", "", datadir+"/bx20.txt");
  fprintf(stderr,"h3 %f\n", h );
  BOOST_REQUIRE( h > -0.01 );

  fprintf(stdout, "Perplexity of ab interp 0.0:\n");
  h = perplexity(datadir+"/a.arpa", datadir+"/b.arpa", datadir+"/abx20.txt", 0.01);
  fprintf(stderr,"h4 %f\n", h );

  fprintf(stdout, "Perplexity of ab interp 0.5:\n");
  h = perplexity(datadir+"/a.arpa", datadir+"/b.arpa", datadir+"/abx20.txt", 0.5);
  fprintf(stderr,"h5 %f %f\n", h, fabs(h-1.0) );
  BOOST_REQUIRE( fabs(h-1.0) < 0.01 );

  fprintf(stdout, "Perplexity of ab interp 1.0:\n");
  h = perplexity(datadir+"/a.arpa", datadir+"/b.arpa", datadir+"/abx20.txt", 0.99);
  fprintf(stderr,"h6 %f\n", h );
  // BOOST_REQUIRE( abs(lp) < 0.01 );

  fprintf(stdout, "Perplexity of ab interp 0.0 against a:\n");
  h = perplexity(datadir+"/a.arpa", datadir+"/b.arpa", datadir+"/ax20.txt", 1.00);
  fprintf(stderr,"h7 %f\n", h );
  BOOST_REQUIRE( h > -0.01 );
}

void test_interpolated_different_vocabs(std::string datadir) {
  create_lm(datadir+"/ax20.txt", datadir+"/vocab.txt", datadir+"/ax20.txt", 3, datadir+"/a-novocab.arpa");
  create_varigram_lm(datadir+"/bx20.txt", "", datadir+"/bx20.txt", 3, datadir+"/b-novocab.arpa");
  fprintf(stdout, "Perplexity of ab interp 0.5 (novocab):\n");
  float h = perplexity(datadir+"/a-novocab.arpa", datadir+"/b-novocab.arpa", datadir+"/abx20.txt", 0.5);
  fprintf(stderr,"novocab %f\n", h );
  BOOST_REQUIRE( fabs(h-1.0) < 0.01 );
}

void test_intertreegram(std::string datadir) {
  std::vector< std::string > lm_names;
  lm_names.push_back(std::string(datadir+"/a-novocab.arpa"));
  lm_names.push_back(std::string(datadir+"/b-novocab.arpa"));
  
  InterTreeGram itg(lm_names);
   
}

int test_main( int argc, char *argv[] )             // note the name!
{
  // FIXME: Use the BOOST unit test framework properly
  fprintf(stderr, "Running tests\n"); 
  fprintf(stderr, "datadir: %s\n", argv[1]);
  std::string datadir(argv[1]);
  
  // FIXME: Write unit test to read and write a arpa file with treegram and hashgram
  create_simple_models(datadir);
  test_interpolation(datadir);
  test_interpolated_different_vocabs(datadir);
  test_intertreegram(datadir);

  return 0;
  /*
    // six ways to detect and report the same error:
    BOOST_CHECK( add( 2,2 ) == 4 );        // #1 continues on error
    BOOST_REQUIRE( add( 2,2 ) == 4 );      // #2 throws on error
    if( add( 2,2 ) != 4 )
      BOOST_ERROR( "Ouch..." );            // #3 continues on error
    if( add( 2,2 ) != 4 )
      BOOST_FAIL( "Ouch..." );             // #4 throws on error
    if( add( 2,2 ) != 4 ) throw "Oops..."; // #5 throws on error

    return add( 2, 2 ) == 4 ? 0 : 1;       // #6 returns error code
  */
}


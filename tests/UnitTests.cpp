#include <boost/test/minimal.hpp>
#include <stdio.h>

#include <InterKn.hh>
#include <VarigramFuncs.hh>
#include <PerplexityFuncs.hh>

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

float perplexity(std::string model1, std::string model2, std::string infname, float alpha=0.5) {
  io::Stream txtin(infname,"r");
  io::Stream out("-","w");
  
  Perplexity lm(model1, 0, "", "", "", false, false);
  if (model2 != "") {
    lm.set_interpolation(model2);
    lm.set_alpha(alpha);
  }
  float logprob = lm.logprob_file(txtin.file, NULL);
  lm.print_results(out.file);
  return logprob;
}

int test_main( int argc, char *argv[] )             // note the name!
{
  fprintf(stderr, "Running tests\n"); 
  std::string datadir(argv[1]);
  fprintf(stderr, "datadir: %s\n", argv[1]);

  create_lm(datadir+"/ax20.txt", datadir+"/vocab.txt", datadir+"/ax20.txt", 3, datadir+"/a.arpa");

  fprintf(stdout, "Perplexity of a against a:\n");
  float lp = perplexity(datadir+"/a.arpa", "", datadir+"/ax20.txt");
  fprintf(stderr,"lp1 %f\n", lp );
  BOOST_REQUIRE( abs(lp) < 0.01 );

  fprintf(stdout, "Perplexity of b against a:\n");
  create_varigram_lm(datadir+"/bx20.txt", datadir+"/vocab.txt", datadir+"/bx20.txt", 3, datadir+"/b.arpa");
  lp = perplexity(datadir+"/b.arpa", "", datadir+"/ax20.txt");
  fprintf(stderr,"lp2 %f\n", lp );
  BOOST_REQUIRE( lp < -60.0 );

  fprintf(stdout, "Perplexity of b against b:\n");
  lp = perplexity(datadir+"/b.arpa", "", datadir+"/bx20.txt");
  fprintf(stderr,"lp3 %f\n", lp );
  BOOST_REQUIRE( abs(lp) < 0.01 );

  fprintf(stdout, "Perplexity of ab interp 0.0:\n");
  lp = perplexity(datadir+"/a.arpa", datadir+"/b.arpa", datadir+"/abx20.txt", 0);
  fprintf(stderr,"lp4 %f\n", lp );
  // BOOST_REQUIRE( abs(lp) < 0.01 );

  fprintf(stdout, "Perplexity of ab interp 0.5:\n");
  lp = perplexity(datadir+"/a.arpa", datadir+"/b.arpa", datadir+"/abx20.txt", 0.5);
  fprintf(stderr,"lp5 %f\n", lp );
  // BOOST_REQUIRE( abs(lp) < 0.01 );

  fprintf(stdout, "Perplexity of ab interp 1.0:\n");
  lp = perplexity(datadir+"/a.arpa", datadir+"/b.arpa", datadir+"/abx20.txt", 1.0);
  fprintf(stderr,"lp6 %f\n", lp );
  // BOOST_REQUIRE( abs(lp) < 0.01 );


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


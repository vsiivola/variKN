#include <boost/test/minimal.hpp>
#include <stdio.h>

#include <InterKn.hh>
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

int test_main( int argc, char *argv[] )             // note the name!
{
  fprintf(stderr, "Running tests\n"); 
  std::string datadir(argv[1]);
  fprintf(stderr, "datadir: %s\n", argv[1]);
  create_lm(datadir+"/ax20.txt", datadir+"/vocab.txt", datadir+"/ax20.txt", 3, datadir+"/a.arpa");
  
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


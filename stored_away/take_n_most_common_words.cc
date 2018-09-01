// Copyright (C) 2007  Vesa Siivola.
// See licence.txt for the terms of distribution.

// Program to remove the least common words from counts file
#include "NgramCounts.hh"
#include "conf.hh"
#include "io.hh"

int main(int argc, char **argv) {
  conf::Config config;
  config(
      "Usage: take_n_most_common_words counts_in counts_out\nDrops the least "
      "common words and trims the counts and vocab files accordingly.\n")(
      'n', "norder=INT", "arg must", "3",
      "n-gram order")('B', "vocabin=FILE", "arg", "",
                      "Input count file contains indices, this is the "
                      "corresponding vocabulary")(
      'b', "vocabout=FILE", "arg", "",
      "Write indices to output counts. Write the corresponding vocabulary to "
      "this file")('f', "nfirst=INT", "arg", "99999999",
                   "Number of most common words to be included")(
      'd', "ndrop=FLOAT", "arg", "0",
      "Drop all words with less than ndrop occurances. If both nfirst and "
      "ndrop options are specified, the tighter bound is taken");

  config.parse(argc, argv, 2);
  int nfirst = config["nfirst"].get_int();
  int n = config["norder"].get_int();
  int ndrop = (int)(n * config["ndrop"].get_double());

  io::Stream::verbose = true;
  io::Stream trigramin(config.arguments.at(0), "r");
  io::Stream trigramout(config.arguments.at(1), "w");

  io::Stream vocabin, vocabout;
  if (config["vocabin"].specified)
    vocabin.open(config["vocabin"].get_str(), "r");
  if (config["vocabout"].specified)
    vocabout.open(config["vocabout"].get_str(), "w");

  if (!ndrop && nfirst == 99999999) {
    fprintf(stderr, "Must specify either -ndrop or -nfirst.\nExit.\n");
    exit(0);
  }

  NgramCounts_t<int, int> nc(n, 0, 0);
  nc.read(trigramin.file, vocabin.file);
  nc.shrink(ndrop, nfirst);
  nc.write(trigramout.file, vocabout.file, true);
}

#if 0

/* The old code, probably more memory efficient */ 
int main (int argc, char **argv) {
  struct fdesc vocabin, vocabout, trigramin, trigramout;
  char sbuf[MAX_WLEN];
  int count;
  int *backwards;
  int i;
  std::vector<int> ngram;
  int nfirst,ndrop;
  int num_ngrams=0;
  struct matrix *trigrams;

  //struct sortstruct *wordtable;
  //wordtable=(struct sortstruct *) malloc(sizeof(struct sortstruct *)*MAX_VOCAB);
  std::vector<struct sortstruct> wordtable;

  if (pi_extract_parameter(argc,argv,"-h",OPTION2)) {
    fprintf(stderr,"This program will drop the ngrams having the most uncommon words and trim the vocab file accordingly. Options:\n -vocabin original vocabulary\n -vocabout output trimmed vocabulary\n -triin original ngram file\n -triout trimmed ngram file\n -nfirst number of most common words to take\n -ndrop drop words with less occurances than this\n If both -nfirst and -ndrop are specified, the tighter bound is taken.\n");
    exit(0);
  }

  nfirst=pi_oatoi(pi_extract_parameter(argc,argv,"-nfirst",OPTION),99999999);

  vocabin=pi_open_file(argc,argv,"-vocabin",ALWAYS,"r",1);
  vocabout=pi_open_file(argc,argv,"-vocabout",ALWAYS,"w",1);
  trigramin=pi_open_file(argc,argv,"-triin",ALWAYS,"r",1);
  trigramout=pi_open_file(argc,argv,"-triout",ALWAYS,"w",1);
  int n=pi_oatoi(pi_extract_parameter(argc,argv,"-n",OPTION),3);
  ndrop=(int) (n*pi_ostrtod(pi_extract_parameter(argc,argv,"-ndrop",
						OPTION),0.0));
  if (!ndrop && nfirst == 99999999) {
    fprintf(stderr,"Must specify either -ndrop or -nfirst.\nExit.\n");
    exit(0);
  }
  ngram.resize(n);

  /* Read vocab, reserve index 0 for UNK */
  struct sortstruct ss;
  ss.oldidx=-1;
  ss.count=-1;
  /*wordtable[0].wstring=strdup("!UNKNOWN");*/
  while (fgets(sbuf,MAX_WLEN,vocabin.p)) {
    if (*sbuf=='#') continue;
    pi_chomp(sbuf);
    ss.wstring=strdup(sbuf);
    wordtable.push_back(ss);
  }
  pi_fclose(vocabin);

  if (nfirst>wordtable.size() && ndrop<1) {
    fprintf(stderr,"Requested more words than there is in the original vocabulary.\n Nothing to do. Exit.\n");
    exit(0);
  }

  /* Read word counts */
  while (fgets(sbuf,MAX_WLEN,trigramin.p)) {
    char *cptr;
    cptr=strtok(sbuf," ");
    for (i=0;i<n;i++) {
      if (sscanf(cptr,"%d",&(ngram[i]))!=1) {
	fprintf(stderr,"Error reading ngram file\n");
	exit(-1);
      }
      cptr=strtok(NULL," ");
    }
    if (sscanf(cptr,"%d",&count)!=1) {
      fprintf(stderr,"Error reading count from ngram file\n");
      exit(-1);
    }

    for (i=0;i<n;i++) {
      wordtable[ngram[i]].oldidx=ngram[i];
      wordtable[ngram[i]].count+=count;
    }
    num_ngrams++;
  }

  /* Sort words */
  qsort(&wordtable[1],wordtable.size()-1,sizeof(struct sortstruct), countcmp);

  /* Create a backwards map */
  backwards=(int *) calloc(wordtable.size(),sizeof(int));
  for (i=1;i<=wordtable.size();i++) backwards[i]=0; /* All words unknown */
  for (i=1;(i<=nfirst && wordtable[i].count>ndrop);i++) { 
    backwards[wordtable[i].oldidx]=i; /* Add known words */
    /*fprintf(stderr,"back %d -> %d\n",wordtable[i].oldidx,i);*/
  }
  nfirst=i-1;

  /* Print new vocabulary */
  fprintf(vocabout.p,"## Vocab generated by take_n_most_common_words, indice 0 in ngrams for UNK\n##\n##\n## Includes %d words ##\n",nfirst);
  for (i=1;i<=nfirst;i++) {
    fprintf(vocabout.p,"%s\n",wordtable[i].wstring);
  }
  /*fprintf(vocabout.p,"<s>\n");*/
  pi_fclose(vocabout);

  /* Reread trigram file, print remapped trigrams */
  pi_fclose(trigramin);
  trigramin=pi_open_file(argc,argv,"-triin",ALWAYS,"r",0);

  for (i=0;i<n;i++)
    ngram[i]=nfirst+1;
  trigrams=CreateMatrixI(i,&ngram[0],num_ngrams+1,0);
  /* Read word counts */
  while (fgets(sbuf,MAX_WLEN,trigramin.p)) {
    char *cptr;
    cptr=strtok(sbuf," ");
    for (i=0;i<n;i++) {
      if (sscanf(cptr,"%d",&ngram[i])!=1) {
	fprintf(stderr,"Error reading ngram file\n");
	exit(-1);
      }
      cptr=strtok(NULL," ");
    }
    if (sscanf(cptr,"%d",&count)!=1) {
      fprintf(stderr,"Error reading count from ngram file\n");
      exit(-1);
    }

    for (i=0;i<n;i++) 
      ngram[i]=backwards[ngram[i]];
    IncrementI(trigrams,&ngram[0],count);
  }

  OrderedStepThrough(trigrams,&ngram[0],&count);
  while (OrderedStepThrough(NULL,&ngram[0],&count)) {
    for (i=0;i<n;i++) 
      fprintf(trigramout.p,"%d ",ngram[i]); 
    fprintf(trigramout.p,"%d\n",count);
  }
  pi_fclose(trigramin);
  pi_fclose(trigramout);
}

#endif

// Copyright (C) 2007  Vesa Siivola. 
// See licence.txt for the terms of distribution.

// Easily modifiable presentation of n-gram language model
#include "HashGram.hh"
#include "str.hh"
#include "def.hh"

template <typename KT>
HashGram_t<KT>::~HashGram_t() {
  for (size_t i=1;i<probs.size();i++) delete probs[i];
  for (size_t i=1;i<backoffs.size();i++) delete backoffs[i];
}

void HashGram::read_error()
{
  fprintf(stderr, "HashGram::read(): error on line %d\n", m_lineno);
  exit(1);
}

template <typename KT>
void HashGram_t<KT>::read_real(FILE *file) {
  std::string line;
  std::vector<std::string> vec;

  // Just for efficiency
  line.reserve(128); 
  vec.reserve(16);

  bool ok = true;

  m_lineno = 0;

  // Find header
  while (1) {
    ok = str::read_line(&line, file, true);
    m_lineno++;

    if (!ok) {
      fprintf(stderr, "HashGram::read(): "
	      "error on line %d while waiting \\data\\", m_lineno);
      exit(1);
    }

    if (line == "\\interpolated") m_type=INTERPOLATED;

    if (line == "\\data\\")
      break;
  }

  // Read header
  int order = 1;
  int number_of_nodes = 0;
  int max_order_count = 0;
  while (1) {
    ok = str::read_line(&line, file, true);
    m_lineno++;

    if (!ok) {
      fprintf(stderr, "HashGram::read(): "
	      "error on line %d while reading counts", m_lineno);
      exit(1);
    }
    
    // Header ends in a \-command
    if (line[0] == '\\')
      break;

    // Skip empty lines
    if (line.find_first_not_of(" \t\n") == line.npos)
      continue;

    // All non-empty header lines must be ngram counts
    if (line.substr(0, 6) != "ngram ")
      read_error();
    {
      std::string tmp(line.substr(6));
      str::split(&tmp, "=", false, &vec);
    }
    if (vec.size() != 2)
      read_error();

    int count = atoi(vec[1].c_str());
    if (count > max_order_count)
      max_order_count = count;
    number_of_nodes += count;
    m_counts.push_back(count);
    
    if (atoi(vec[0].c_str()) != order || m_counts.back() < 0)
      read_error();
    order++;
  }

  m_order=m_counts.size();
  probs.resize(m_order+1);
  backoffs.resize(m_order+1);
  // Read ngrams order by order
  for (order = 1; order <= m_counts.size(); order++) {
    fprintf(stderr,"Reserving %d grams for order %d\n", m_counts[order-1], order);
    probs[order]=new sikMatrix<KT, float>(order,m_counts[order-1], MINLOGPROB);
    backoffs[order]=new sikMatrix<KT, float>(order,m_counts[order-1], 0.0);

    // We must always have the correct header line at this point
    if (line[0] != '\\') {
      fprintf(stderr, "HashGram::read(): "
	      "\\%d-grams expected on line %d\n", order, m_lineno);
      exit(1);
    }
    str::clean(&line, " \t");
    str::split(&line, "-", false, &vec);

    if (atoi(vec[0].substr(1).c_str()) != order || vec[1] != "grams:") {
      fprintf(stderr, "HashGram::read(): "
	      "unexpected command on line %d: %s\n", m_lineno, line.c_str());
      exit(1);
    }

    // Read the grams of each order into the sorter
    std::vector<KT> gram(order);
    for (int w = 0; w < m_counts[order-1]; w++) {

      // Read and split the line
      if (!str::read_line(&line, file))
	read_error();
      str::clean(&line, " \t\n");
      m_lineno++;

      // Ignore empty lines
      if (line.find_first_not_of(" \t\n") == line.npos) {
	w--;
	continue;
      }

      str::split(&line, " \t", true, &vec);

      // Check the number of columns on the line
      if (vec.size() < order + 1 || vec.size() > order + 2) {
	fprintf(stderr, "HashGram::read(): "
		"%d columns on line %d\n", (int) vec.size(), m_lineno);
	exit(1);
      }
      if (order == m_counts.size() && vec.size() != order + 1)
	fprintf(stderr, "WARNING: %d columns on line %d\n", (int) vec.size(), 
		m_lineno);

      // FIXME: should we deny new words in higher order ngrams?

      // Parse log-probability, back-off weight and word indices
      // FIXME: check the conversion of floats
      float log_prob = strtod(vec[0].c_str(), NULL);
      float back_off = 0;
      if (vec.size() == order + 2)
	back_off = strtod(vec[order + 1].c_str(), NULL);

      // Add the gram to sorter
      //fprintf(stderr,"add gram [");
      for (int i = 0; i < order; i++) {
	gram[i] = add_word(vec[i + 1]);
      }
      if (log_prob>MINLOGPROB) {
	probs[order]->setvalue(&gram[0],log_prob);
	//fprintf(stderr, "%f ",log_prob); 
      }
      //print_indices(gram); 
      if (back_off<0) {
	backoffs[order]->setvalue(&gram[0],back_off); 
	//fprintf(stderr," %f",back_off);
      }
      //fprintf(stderr,"\n");
    }

    // Skip empty lines before the next order.
    while (1) {
      if (!str::read_line(&line, file, true)) {
	if (ferror(file))
	  read_error();
	if (feof(file))
	  break;
      }
      m_lineno++;

      if (line.find_first_not_of(" \t\n") != line.npos)
	break;
    }
  }
}

template <typename KT>
void HashGram_t<KT>::write_real(FILE *out) {
  std::vector<std::string> strbuf;
  strbuf.reserve(100000);
  std::vector<int> num_grams(m_order+1,0);

  for (int o=1; o<m_order; o++) {
    std::vector<KT> gram(o);
    float logprob, bo;
    backoffs[o]->stepthrough(true, &gram[0], &logprob);     
    while (backoffs[o]->stepthrough(false,&gram[0],&bo)) {
      if (bo>=0.0) continue;
      FindEntry(probs[o]->m, (byte *) (&gram[0]), 1);
    }
  }

  for (int o=1;o<=m_order;o++) {
    std::vector<KT> gram(o);
    float logprob, bo;
    probs[o]->stepthrough(true, &gram[0], &logprob);
    while (probs[o]->stepthrough(false,&gram[0],&logprob)) {
      bo = backoffs[o]->getvalue(&gram[0]);
      if (m_print_zerograms || logprob>MINLOGPROB || bo < 0 )
	num_grams[o]++;
    }
  }

  if (m_type==INTERPOLATED) fprintf(out,"\\interpolated\n");
  fprintf(out,"\\data\\\n");

  for (int i=1;i<=m_order;i++) {
    fprintf(out,"ngram %d=%d\n",i,num_grams[i]);
  }

  for (int o=1;o<=m_order;o++) {
    fprintf(out,"\n\\%d-grams:\n",o);
    std::vector<KT> gram(o);
    float logprob, bo;
    probs[o]->stepthrough(true, &gram[0], &logprob);
    while (probs[o]->stepthrough(false,&gram[0],&logprob)) {
      bo = backoffs[o]->getvalue(&gram[0]);
      if (!m_print_zerograms && (logprob<=MINLOGPROB && bo >= 0 )) continue;
      fprintf(out,"%.4f",logprob);
      for (int i=0;i<o;i++) fprintf(out," %s",word(gram[i]).c_str());
      if (bo <0) 
	fprintf(out, " %.4f",bo);
      fprintf(out,"\n");
    }

  }
  fprintf(out,"\\end\\\n");
}

template <typename KT>
float HashGram_t<KT>::log_prob_bo_helper(const std::vector<KT> &gram)
{
  const int looptill=std::min(gram.size(),(size_t) m_order);
  int n=looptill;
  m_last_order=1;
  float log_prob = 0.0;
  const KT* const gram_ptr = &gram[gram.size()-looptill];

  while (true) {
    //fprintf(stderr, "looptill %d, gramidx %d, ", looptill,looptill-n);
    //print_indices(&gram_ptr[looptill-n],n);
    const float prob=probs[n]->getvalue(&gram_ptr[looptill-n]);

    // Full gram found?
    if (prob>MINLOGPROB) {
      //fprintf(stderr,"full gram n=%d, prob %.4f(%.4f), limit %.4f\n",n,prob,log_prob+prob,MINLOGPROB); 
      log_prob += prob;
      m_last_order = n;
      break;
    } else if (n==1) {
      log_prob += MINLOGPROB;
    }
    
    if (n==1) break;
    log_prob += backoffs[n-1]->getvalue(&gram_ptr[looptill-n]);
    n--;
  }
  return log_prob;
}

template <typename KT>
float HashGram_t<KT>::log_prob_i_helper(const std::vector<KT> &gram) {
  float prob=0.0;
  m_last_order=0;

  const int looptill=std::min(gram.size(),(size_t) m_order);
  for (int n=1;n<=looptill;n++) {
    if (n>1) prob *= pow(10, backoffs[n-1]->getvalue(&gram[gram.size()-n]));
    
    float p = probs[n]->getvalue(&gram[gram.size()-n]);
    if (p<=MINLOGPROB && n>1) {
        continue;
    }
    m_last_order=n;
    prob+=pow(10,p);
  }
  return((float)safelogprob(prob));
}

// Overloading stuff
template <> inline 
float HashGram_t<unsigned short>::log_prob_bo(const std::vector<int> &gram) {
  std::vector<unsigned short> g(gram.size());
  for (size_t i=0;i<gram.size();i++) 
    g[i]=gram[i];
  return(log_prob_bo_helper(g));
}

template <> inline 
float HashGram_t<unsigned short>::log_prob_i(const std::vector<int> &gram) {
  std::vector<unsigned short> g(gram.size());
  for (size_t i=0;i<gram.size();i++) 
    g[i]=gram[i];
  return (log_prob_i_helper(g));
}

template <> inline 
float HashGram_t<int>::log_prob_bo(const std::vector<int> &gram) {
  return log_prob_bo_helper(gram);
}

template <> inline 
float HashGram_t<int>::log_prob_i(const std::vector<int> &gram) {
  return log_prob_i_helper(gram);
}

template <> inline 
float HashGram_t<unsigned short>::log_prob_bo(const std::vector<unsigned short> &gram) {
  return log_prob_bo_helper(gram);
}

template <> inline 
float HashGram_t<unsigned short>::log_prob_i(const std::vector<unsigned short> &gram) {
  return log_prob_i_helper(gram);
}

template <> inline 
float HashGram_t<int>::log_prob_bo(const std::vector<unsigned short> &gram) {
  assert(false);
  return 0;
}

template <> inline 
float HashGram_t<int>::log_prob_i(const std::vector<unsigned short> &gram) {
  assert(false);
  return 0;
}

template <typename KT>
void HashGram_t<KT>::prune(float treshold) {
  // As in Stolcke's paper
  assert(m_type==BACKOFF);
  bool verbose=0;
  std::vector<std::vector<float> > hist_cumuprobs;

  for (int o=m_order;o>1;o--) { // Does not prune unigrams...
    if (verbose) fprintf(stderr,"Processing order %d\n",o);
    std::vector<KT> gram(o);
    ///////////////////////////////////////////////////////////////
    // 1. Cache the total prob mass given to backoff = 
    //                               1-total unbackoffed mass
    twofloat f2(0.0,0.0);
    sikMatrix<KT, twofloat> bo_pcache(o-1,probs[o-1]->num_entries(), f2);
    {
      float gram_logprob;
      std::vector<KT> g2(o-1);
      probs[o]->stepthrough(true, &gram[0], &gram_logprob);
      while (probs[o]->stepthrough(false, &gram[0], &gram_logprob)) {
	f2.float1=pow(10,gram_logprob);
	for (int i=1;i<o;i++) g2[i-1]=gram[i];
	f2.float2=pow(10,log_prob(g2));
	bo_pcache.increment(&gram[0], f2);
      }
    }
    if (verbose) fprintf(stderr,"Cached backoff mass.\n");
    
    float gram_logprob;
    probs[o]->stepthrough(true, &gram[0], &gram_logprob);
    while (probs[o]->stepthrough(false, &gram[0], &gram_logprob)) {
      float trunc_logprob;
      {
	std::vector<KT> g2(o-1);
	for (int i=1;i<o;i++) g2[i-1]=gram[i];
	trunc_logprob=log_prob(g2);
	if (verbose) {
	  fprintf(stderr,"processing ");
	  print_indices(stderr, gram);
	  fprintf(stderr,", gram lp %.4f, trunc %.4f\n", gram_logprob, trunc_logprob);
	}
      }
      ///////////////////////////////////////////////////////////
      // 2. Mariginal history probabilities p(h) = p(h1)p(h2|h1)...
      float hist_log_prob=0.0;
      // This is the correct way
      std::vector<KT> g2(gram);
      for (int o2=o-1;o2>=1;o2--) {
	g2.resize(o2);
	hist_log_prob += log_prob(g2);
      }

      if (verbose) fprintf(stderr,"histprob %.4f, ",hist_log_prob);

      ////////////////////////////////////////////////////////
      // 3. The values of the a(h) and a(h')
      float alpha, alpha2;
      float bo_psum;
      float new_bo_psum;
      float new_bo_psum2;
      {
	twofloat f2=bo_pcache.getvalue(&gram[0]);
	bo_psum = f2.float1;
	const float bo_psum2 = f2.float2;
	if (verbose) 
	  fprintf(stderr,"bo_psum %.4f, bo_psum2 %.4f\n",bo_psum, bo_psum2);
	alpha = (1.0 - bo_psum) / ( 1.0 - bo_psum2);
	new_bo_psum = bo_psum - pow(10, gram_logprob);
	new_bo_psum2 = bo_psum2 - pow(10, trunc_logprob);
	alpha2 = (1.0 - new_bo_psum) / ( 1.0 - new_bo_psum2);
      }
      if (verbose) fprintf(stderr,"a=%.4f, a'=%.4f\n",alpha, alpha2);
      //////////////////////////////////////////////////////////////
      // 4. Difference in entropy = 
      //      -p(h) { p(w|h) [ log p(w|h') + log a'(h) - log p(w|h)] 
      //                + [ log a(h') - log a(h) ] * backoffed_mass
      float diff;
      float log_alpha2;
      { 
	const float log_alpha=(float)safelogprob(alpha);
	const float hist_prob=(float)pow(10, hist_log_prob);
	const float gram_prob=(float)pow(10, gram_logprob);
	log_alpha2=(float)safelogprob(alpha2);

	diff = - hist_prob
	  * (gram_prob * (trunc_logprob + log_alpha2 - gram_logprob)
	    + (1-bo_psum) *(log_alpha2-log_alpha));
	if (verbose) {
	  fprintf(stderr,"diff = %.4f -p(h) * ", -hist_prob);
	  fprintf(stderr,"( %4f p(w|h) * ( %.4f logp(w|h') + %.4f log a'(h) - %.4f log p(w|h)) ", gram_prob, trunc_logprob, log_alpha2, gram_logprob);
	  fprintf(stderr,"%.4f bomass * ( %.4flog a(h') - %.4f loga(h)))\n", 1-bo_psum, log_alpha2, log_alpha);
	  fprintf(stderr," = diff %.4g (%.4g)\n", diff, treshold);
	}
      }
      if (diff > treshold) continue;
      
      /////////////////////////////////////////////////////////////
      // 5. Modify the model
      {
	if (verbose) fprintf(stderr,"rejecting\n");

	// MINLOGPROB-1 is not the default value MINLOGPROB. This prevents
	// the model from reordering the internal ordering and breaking
	// steptrough(). Really bad interface design...
	probs[o]->setvalue(&gram[0],MINLOGPROB-1); 

	if (verbose) 
	  fprintf(stderr,"Setting backoff %d to %.4f\n", gram[0], log_alpha2);
	backoffs[o-1]->setvalue(&gram[0], log_alpha2);
	twofloat f2(new_bo_psum, new_bo_psum2);
	bo_pcache.setvalue(&gram[0], f2);
      }
    }
  }
  remove_empty_grams();
}

template <typename KT>
void HashGram_t<KT>::remove_empty_grams() {
  for (int o=m_order;o>1;o--) {
    std::vector<KT> gram(o);
    for (indextype i=probs[o]->num_entries()-1;i>=0;i--) {
      //fprintf(stderr,"Checking "); 
      //print_indices(probs[o]->Idx2Keyp(i), o);
      //fprintf(stderr," =%.4f ",*(probs[o]->Idx2Valp(i)));
      //fprintf(stderr," (%.4f)\n", backoffs[o]->getvalue(probs[o]->Idx2Keyp(i)));
      if (*(probs[o]->Idx2Valp(i)) <= MINLOGPROB) {
	probs[o]->setvalue(probs[o]->Idx2Keyp(i), MINLOGPROB);
      }
    }
  }
}

template <typename KT>
void HashGram_t<KT>::add_zeroprob_grams() {
  for (int o=m_order;o>=2;o--) { 
    std::vector<KT> gram(o);
    float gram_logprob;
    const float inc=0.0;
    m_print_zerograms=true;

    probs[o]->stepthrough(true, &gram[0], &gram_logprob);
    while (probs[o]->stepthrough(false, &gram[0], &gram_logprob)) {
      //Don't delete even if default value
      probs[o-1]->increment_wo_del(&gram[0],inc); 
    }
  }
}


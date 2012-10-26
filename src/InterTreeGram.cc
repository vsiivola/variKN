#include "io.hh"
#include "def.hh"
#include "TreeGramArpaReader.hh"
#include "InterTreeGram.hh"

InterTreeGram::InterTreeGram(std::vector< std::string > lm_names, std::vector<float> coeffs) {
  if (lm_names.size() != coeffs.size()) {
    fprintf(stderr, "InterTreeGram::InterTreeGram: There must be as many interpolation coeffs as there are LMs. Exit.\n");
    exit(1);
  }

  float coeff_sum=0.0;
  for(std::vector<float>::iterator j=coeffs.begin();j!=coeffs.end();++j) {
    coeff_sum += *j;
  }

  if (coeff_sum < 0.99 || coeff_sum>1.01) {
    fprintf(stderr, "InterTreeGram::InterTreeGram: Interpolation coeffs must sum to 1 (!=%f). Exit.\n", coeff_sum);
    exit(1);
  }
  m_coeffs = coeffs;

  // Combine vocab from all models
  for ( std::vector<std::string>::iterator it = lm_names.begin(); it != lm_names.end(); it ++) {
    fprintf(stdout, "INTERTREEGRAM %s\n", it->c_str());
    ArpaReader areader(this);
    io::Stream lm_in(*it, "r");
    
    std::string line;
    bool dummy;
    areader.read_header(lm_in.file, dummy, line);

    std::vector<int> tmp_gram(1);
    float log_prob, back_off;
    // This adds all words to vocab and discards the other information
    while( areader.next_gram(lm_in.file, line, tmp_gram, log_prob, back_off) && tmp_gram.size()==1) {}
  }

  // Initialize all models
  for ( std::vector<std::string>::iterator it = lm_names.begin(); it != lm_names.end(); it ++) {
    io::Stream lm_in(*it, "r");
    TreeGram *lm = new TreeGram;

    // Copy vocab;
    for (int i=0; i<num_words(); i++) {
      lm->add_word(word(i));
    }
    int real_num_words = num_words();
    //copy_vocab_to(*lm);
    assert(lm->num_words() == real_num_words);

    // Read n_grams
    TreeGramArpaReader tgar;
    tgar.read(lm_in.file, lm);
    m_models.push_back(lm);
    assert(lm->num_words() == real_num_words);
  }
}

InterTreeGram::~InterTreeGram(void) {
  for (std::vector<TreeGram *>::iterator j=m_models.begin();j!=m_models.end();++j) {
    delete *j;
  }
}

float InterTreeGram::log_prob(std::vector<int> &gram) {
  double prob=0.0;
  for (int i=0; i<m_models.size(); i++) {
    prob += m_coeffs[i] * pow(10, m_models[i]->log_prob(gram));
  }
  return safelogprob(prob);
}

void InterTreeGram::test_write(std::string fname, int idx) {
  io::Stream lm_out(fname, "w");
  TreeGramArpaReader tga;

  tga.write(lm_out.file, m_models[idx]);
}
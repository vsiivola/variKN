#include "io.hh"
#include "TreeGramArpaReader.hh"
#include "InterTreeGram.hh"

InterTreeGram::InterTreeGram(std::vector< std::string > lm_names) {
  initialize_models(lm_names);
}

InterTreeGram::~InterTreeGram(void) {
  for (std::vector<TreeGram *>::iterator j=m_models.begin();j!=m_models.end();++j) {
    delete *j;
  }
}

void InterTreeGram::initialize_models(std::vector< std::string > lm_names) {
  // Combine vocab from all models
  for ( std::vector<std::string>::iterator it = lm_names.begin(); it != lm_names.end(); it ++) {
    fprintf(stdout, "INTERTREEGRAM %s\n", it->c_str());
    ArpaReader areader(this);
    io::Stream lm_in(*it, "r");
    
    std::string line;
    bool dummy;
    areader.read_header(lm_in.file, dummy, line);
  }

  // Initialize all models
  for ( std::vector<std::string>::iterator it = lm_names.begin(); it != lm_names.end(); it ++) {
    io::Stream lm_in(*it, "r");
    TreeGram *lm = new TreeGram;

    // Copy vocab;
    for (int i=0; i<num_words(); i++) {
      lm->add_word(word(i));
    }

    // Read n_grams
    TreeGramArpaReader tgar;
    tgar.read(lm_in.file, lm);
    m_models.push_back(lm);
  }
}

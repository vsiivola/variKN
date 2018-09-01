// Routines for reading arpa format files to the internal prefix tree
// format.
#include "ArpaReader.hh"
#include "def.hh"
#include "str.hh"
#include <stdlib.h>

void ArpaReader::read_error() {
  fprintf(stderr, "ArpaReader::read(): error on line %d\n", m_lineno);
  exit(1);
}

void ArpaReader::read_header(FILE *file, bool &interpolated,
                             std::string &line) {
  int order;
  // Just for efficiency
  m_vec.reserve(16);

  bool ok = true;
  interpolated = false;

  m_lineno = 0;

  // Find header
  while (1) {
    ok = str::read_line(&line, file, true);
    m_lineno++;

    if (!ok) {
      fprintf(stderr,
              "ArpaReader::read(): "
              "error on line %d while waiting \\data\\",
              m_lineno);
      exit(1);
    }

    // iARPA is used by the IRSTLM toolkit
    if (line == "\\interpolated" || line == "iARPA")
      interpolated = true;

    if (line == "\\data\\")
      break;
  }

  // Read header
  order = 1;
  int max_order_count = 0;
  while (1) {
    ok = str::read_line(&line, file, true);
    m_lineno++;

    if (!ok) {
      fprintf(stderr,
              "ArpaReader::read(): "
              "error on line %d while reading counts",
              m_lineno);
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
      str::split(&tmp, "=", false, &m_vec);
    }
    if (m_vec.size() != 2)
      read_error();

    int count = atoi(m_vec[1].c_str());
    if (count > max_order_count)
      max_order_count = count;
    counts.push_back(count);

    if (atoi(m_vec[0].c_str()) != order || counts.back() < 0)
      read_error();
    order++;
  }
}

bool ArpaReader::next_gram(FILE *file, std::string &line,
                           std::vector<int> &gram, float &log_prob,
                           float &back_off) {
  // Read ngrams order by order
  while (m_read_order == 0 || m_gram_num >= counts[m_read_order - 1]) {
    m_gram_num = 0;
    m_read_order++;

    // Skip empty lines before the next order.
    bool skip_empty_lines = line != "\\1-grams:";
    while (skip_empty_lines) {
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

    // We must always have the correct header line at this point
    if (m_read_order > counts.size()) {
      if (line != "\\end\\") {
        fprintf(stderr,
                "ArpaReader::next_gram():"
                "expected end, got '%s' on line %d\n",
                line.c_str(), m_lineno);
        exit(1);
      }
      return false;
    }

    fprintf(stderr, "Found %d grams for order %d\n", counts[m_read_order - 1],
            m_read_order);

    if (line[0] != '\\') {
      fprintf(stderr,
              "ArpaReader::next_gram(): "
              "\\%d-grams expected on line %d\n",
              m_read_order, m_lineno);
      exit(1);
    }

    str::clean(&line, " \t");
    str::split(&line, "-", false, &m_vec);

    if (atoi(m_vec[0].substr(1).c_str()) != m_read_order ||
        m_vec[1] != "grams:") {
      fprintf(stderr,
              "ArpaReader::next_gram(): "
              "unexpected command on line %d: %s\n",
              m_lineno, line.c_str());
      exit(1);
    }
  }

  gram.resize(m_read_order);
  // Read and split the line
  while (true) {
    if (!str::read_line(&line, file))
      read_error();
    str::clean(&line, " \t\n");
    m_lineno++;

    // Ignore empty lines
    if (line.find_first_not_of(" \t\n") == line.npos) {
      continue;
    }
    break;
  }

  str::split(&line, " \t", true, &m_vec);

  // Check the number of columns on the line
  if (m_vec.size() < m_read_order + 1 || m_vec.size() > m_read_order + 2) {
    fprintf(stderr,
            "ArpaReader::next_gram(): "
            "%d columns on line %d\n",
            (int)m_vec.size(), m_lineno);
    exit(1);
  }
  if (m_read_order == counts.size() && m_vec.size() != m_read_order + 1) {
    fprintf(stderr, "WARNING: %d columns on line %d\n", (int)m_vec.size(),
            m_lineno);
  }

  // FIXME: should we deny new words in higher order ngrams?

  // Parse log-probability, back-off weight and word indices
  // FIXME: check the conversion of floats
  log_prob = strtod(m_vec[0].c_str(), NULL);
  back_off = 0;
  if (m_vec.size() == m_read_order + 2)
    back_off = strtod(m_vec[m_read_order + 1].c_str(), NULL);

  // Add the gram to sorter
  // fprintf(stderr,"add gram [");
  for (int i = 0; i < m_read_order; i++) {
    gram[i] = m_vocab->add_word(m_vec[i + 1]);
  }

  m_gram_num++;
  return true;
}

#ifndef DEF_HH
#define DEF_HH

#include <deque>
#include <math.h>
#include <stdio.h>
#include <vector>

const int MAX_WLEN = 1000;
#define MAX_WLEN_FMT_STRING                                                    \
  "%1000s" // For limiting fscanf width, must be modified to match MAX_WLEN
const double MINLOGPROB = -60;
const double MINPROB = 1e-60;

inline double safelogprob(double x) {
  if (x > MINPROB)
    return (log10(x));
  return (MINLOGPROB);
}

inline double safelogprob2(double x) {
  if (x > MINPROB)
    return (log10(x));
  if (x < 0)
    return (0.0);
  return (MINLOGPROB);
}

///////////////////////////////////////////////////
// For debug

inline void print_indices(FILE *o, const int *const v, const int dim) {
  fprintf(o, "[");
  for (int i = 0; i < dim; i++)
    fprintf(o, " %d", v[i]);
  fprintf(o, " ]");
}

inline void print_indices(FILE *o, const unsigned short *v, const int dim) {
  fprintf(o, "[");
  for (int i = 0; i < dim; i++)
    fprintf(o, " %d", v[i]);
  fprintf(o, " ]");
}

inline void print_indices(FILE *o, const std::vector<int> &v) {
  print_indices(o, &v[0], v.size());
}

inline void print_indices(FILE *o, const std::vector<unsigned short> &v) {
  print_indices(o, &v[0], v.size());
}

inline void print_findices(FILE *o, const float *v, const int dim) {
  fprintf(stderr, "[");
  for (int i = 0; i < dim; i++)
    fprintf(o, " %.4f", v[i]);
  fprintf(stderr, " ]");
}

inline void print_indices(FILE *o, const std::vector<float> &v) {
  print_findices(o, &v[0], v.size());
}

typedef std::deque<int> Gram;

inline void print_indices(FILE *o, const Gram &v) {
  fprintf(o, "[");
  for (size_t i = 0; i < v.size(); i++)
    fprintf(o, " %d", v[i]);
  fprintf(o, " ]");
}
#endif

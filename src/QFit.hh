// Numerical search, quadratic fit
#ifndef QFIT_HH
#define QFIT_HH

#include <cassert>
#include <vector>

class QFitEvaluator {
public:
  virtual float evaluate(std::vector<float> &) = 0;
  virtual ~QFitEvaluator(){};
};

class QFit {
public:
  inline QFit(float ltol, float atol, QFitEvaluator *qe)
      : m_ltol(ltol), m_atol(atol), m_qe(qe) {}
  std::vector<float> minimize(int maxiter);
  inline void set_initial_point(std::vector<float> v) { m_initial_point = v; }
  inline void set_minimum(std::vector<float> v) { m_low_limits = v; }
  inline void set_maximum(std::vector<float> v) { m_high_limits = v; }
  inline void set_searchstartlim(std::vector<float> v) { m_sslim = v; }

private:
  float m_ltol;
  float m_atol;
  QFitEvaluator *m_qe;
  std::vector<float> m_initial_point;
  std::vector<float> m_res;
  std::vector<float> m_low_limits;
  std::vector<float> m_high_limits;
  std::vector<float> m_sslim;
  void check_limits(int biter, int c, float &l1, float &l2, float &l3,
                    float &t1, float &t2, float &t3);
};

#endif

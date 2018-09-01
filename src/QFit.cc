// Numerical search, quadratic fit
#include "QFit.hh"
#include "def.hh"
#include <math.h>

const float SMALLNUM = 1e-5f;

std::vector<float> QFit::minimize(int maxiter) {
  // fprintf(stderr, "vsizes %zd %zd %zd\n", m_initial_point.size(),
  // m_low_limits.size(), m_high_limits.size());

  assert(m_initial_point.size() > 0);
  assert(m_initial_point.size() == m_low_limits.size() &&
         m_initial_point.size() == m_high_limits.size());

  m_res = m_initial_point;
  float l1 = 0, l2, l3 = m_atol + 1;
  int iter = 0, biter = 0;
  float tolsum = m_atol + 1;
  float utol = 10 * m_ltol;
  while (iter < maxiter) {
    if (biter == 1)
      utol = m_ltol;
    biter++;
    // fprintf(stderr,"tolsum=%g <=> m_atol %g\n",tolsum, m_atol);
    if (tolsum < m_atol) {
      // fprintf(stderr," Out1 : %g\n", l3-l1);
      break;
    }
    tolsum = 0;

    for (int c = 0; c < m_res.size(); c++) {
      // Crude setting for the search limits, this should be improved
      float t1, t2, t3;
      check_limits(biter, c, l1, l2, l3, t1, t2, t3);
      float l = 0, t = 0;

      while (iter < maxiter) {
        iter++;
        // fprintf(stderr,"%d: ",iter);
        if (l3 - l1 < utol) {
          tolsum = tolsum + l3 - l1;
          // fprintf(stderr," Out2 : %g\n", l3-l1);
          break;
        }
        if (t2 >= t1 || t2 >= t3) {
          // Function is not quasiconvex, change search limits:
          // fprintf(stderr,"not qconvex");
          // fprintf(stderr,", orig lim %g < %g <%g", l1,l2,l3);
          if (t2 >= t3) {
            // fprintf(stderr," a\n");
            l1 = l2;
            t1 = l2;
          } else {
            // fprintf(stderr," b\n");
            l3 = l2;
            t3 = t2;
          }
          l2 = l1 + (l3 - l1) / 2;
          // fprintf(stderr,"new lim %g<%g<%g\n", l1,l2,l3);
          m_res[c] = l2;
          t2 = m_qe->evaluate(m_res);
          // fprintf(stderr,"targets %g<%g<%g\n", t1,t2,t3);
        } else {
          // Fit a quadratic function and find the minimum
          // fprintf(stderr,"qfit\n");
          l = (-(t3 - t2) * (l1 * l1 - l2 * l2) +
               (t1 - t2) * (l3 * l3 - l2 * l2)) /
              ((t1 - t2) * (l3 - l2) - (t3 - t2) * (l1 - l2)) / 2;
          // fprintf(stderr, "(-(t3 %g-t2 %g)*(l1 %g*l1-l2 %g*l2)+(t1 %g-t2
          // %g)*(l3%g*l3-l2 %g*l2))/((t1 %g-t2 %g)*(l3 %g-l2 %g)-(t3 %g-t2
          // %g)*(l1 %g-l2 %g))/2=l %g\n", t3, t2, l1, l2,
          // t1,t2,l3,l2,t1,t2,l3,l2,t3,t2,l1,l2,l);
          m_res[c] = l;
          // print_indices(stderr,m_res);fprintf(stderr,"\n");
          t = m_qe->evaluate(m_res);
          if (fabs(l - l2) < utol / 4) {
            if (l2 - l1 < l3 - l2)
              l = l + utol / 2;
            else
              l = l - utol / 2;
          }
          m_res[c] = l;
          // print_indices(stderr,m_res);fprintf(stderr,"\n");
          t = m_qe->evaluate(m_res);

          // Update the limits
          if (l > l2) {
            if (t >= t2) {
              l3 = l;
              t3 = t;
            } else {
              l1 = l2;
              t1 = t2;
              l2 = l;
              t2 = t;
            }
          } else if (l < l2) {
            if (t >= t2) {
              l1 = l;
              t1 = t;
            } else {
              l3 = l2;
              t3 = t2;
              l2 = l;
              t2 = t;
            }
          }
        }
        // fprintf(stderr,"new limits  %g < %g < %g\n", l1, l2, l3);
        // fprintf(stderr,"with values %g < %g < %g\n", t1, t2, t3);
      }
      m_res[c] = l2;
      tolsum = tolsum + l3 - l1;
    }
  }
  return (m_res);
}

void QFit::check_limits(int biter, int c, float &l1, float &l2, float &l3,
                        float &t1, float &t2, float &t3) {

  l1 = m_res[c] - (m_res[c] - m_low_limits[c]) / pow(2.0f, biter - 1);
  l2 = m_res[c];
  l3 = m_res[c] + (m_high_limits[c] - m_res[c]) / pow(2.0f, biter - 1);

  if (m_sslim.size()) { // if initial limits for the search are given
    assert(m_sslim.size() == m_res.size());
    if (l1 < m_res[c] - m_sslim[c])
      l1 = m_res[c] - m_sslim[c];
    if (l3 > m_res[c] + m_sslim[c])
      l3 = m_res[c] + m_sslim[c];
  }

  if (l2 - l1 < SMALLNUM)
    l2 = l1 + SMALLNUM;
  if (l3 - l2 < SMALLNUM)
    l2 = l3 - SMALLNUM;

  // fprintf(stderr,"init_limits %g<%g<%g\n",
  // m_low_limits[c],m_res[c],m_high_limits[c]); fprintf(stderr,"start_limits
  // %g<%g<%g\n", l1,l2,l3);
  while (true) {
    m_res[c] = l1;
    t1 = m_qe->evaluate(m_res);
    m_res[c] = l2;
    t2 = m_qe->evaluate(m_res);
    m_res[c] = l3;
    t3 = m_qe->evaluate(m_res);
    if (t1 < t2) {
      float tmp = l1;
      l1 = std::max(m_low_limits[c], l2 - 2 * (l2 - l1));
      if (l1 == tmp)
        break;
      l3 = l2;
      l2 = tmp;
      // fprintf(stderr,"mid1 %g<%g<%g\n", l1,l2,l3);
      continue;
    }
    if (t3 < t2) {
      float tmp = l3;
      l3 = std::min(m_high_limits[c], l2 + 2 * (l3 - l2));
      if (l3 == tmp)
        break;
      l1 = l2;
      l2 = tmp;
      // fprintf(stderr,"mid2 %g<%g<%g\n", l1,l2,l3);
      continue;
    }
    break;
  }
  // fprintf(stderr,"end_limits %g<%g<%g\n", l1,l2,l3);
}

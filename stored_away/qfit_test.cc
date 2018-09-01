// Copyright (C) 2007  Vesa Siivola.
// See licence.txt for the terms of distribution.

// Test the numerical search

/*
% This is a Matlab test script to produce the "correct" results:
i_vec=[1:1:100];
j_vec=[1:1:100];

A=zeros(length(i_vec), length(j_vec));

for x=i_vec,
  for y=j_vec,
    xx=x/100;
    yy=y/100;
    A(x,y)=17*xx^2+20*yy^2-15*yy-15*xx+19*(xx-0.5)*(yy-0.5);
  end
end
[f ff]=min(A);

%imagesc(A);
surf(A);

[i,j] = min(A); [i2,j2]= min(i);
[j(j2) j2]/100
% (0.51, 0.37)
min(min(A))
% -6.065
*/

#include "QFit.hh"

class Funkkis : public QFitEvaluator {
public:
  inline float evaluate(std::vector<float> &v) {
    fprintf(stderr, "Evaluating at [%f %f]:", v[0], v[1]);
    float res = (17 * v[0] * v[0] + 20 * v[1] * v[1] - 15 * v[1] - 15 * v[0] +
                 19 * (v[0] - 0.5) * (v[1] - 0.5));
    fprintf(stderr, "%g\n", res);
    return (res);
  }
};

int main(int argc, char **argv) {
  Funkkis f;
  QFit qfit(1e-4, 0.9e-4, &f);

  std::vector<float> res(2, 1.0);
  qfit.set_initial_point(res);
  res.clear();
  res.resize(2, -0.5);
  qfit.set_minimum(res);
  res.clear();
  res.resize(2, 1.5);
  qfit.set_maximum(res);

  fprintf(stderr, "Minimization\n");
  res = qfit.minimize(10000);

  fprintf(stderr, "Output of results\n");
  fprintf(stdout, "Minimum: [");
  for (int i = 0; i < res.size(); i++) {
    fprintf(stdout, "%g ", res[i]);
  }
  fprintf(stdout, "]\n");
  fprintf(stdout, "Value at minimum: %g\n", f.evaluate(res));
}

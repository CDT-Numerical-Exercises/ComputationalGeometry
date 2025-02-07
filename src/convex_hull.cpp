#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "convex_hull.h"

double calculate_D(const gsl_vector *p, const gsl_vector *q,
                   const gsl_vector *r) {
  const double px = gsl_vector_get(p, 0);
  const double py = gsl_vector_get(p, 1);
  const double qx = gsl_vector_get(q, 0);
  const double qy = gsl_vector_get(q, 1);
  const double rx = gsl_vector_get(r, 0);
  const double ry = gsl_vector_get(r, 1);

  const double D = (qx*ry - rx*qy) - px*(ry - qy) + py*(rx - qx);

  return D;
}

double calculate_D(const gsl_matrix *X, const size_t pi, const size_t qi,
                   const size_t ri) {
  const gsl_vector_const_view p = gsl_matrix_const_row(X, pi);
  const gsl_vector_const_view q = gsl_matrix_const_row(X, qi);
  const gsl_vector_const_view r = gsl_matrix_const_row(X, ri);
  return calculate_D(&p.vector, &q.vector, &r.vector);
}


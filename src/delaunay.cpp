#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <unordered_set>
#include <boost/container_hash/hash.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "delaunay.h"

const double R3 = sqrt(3);
const double IR3 = 1.0/sqrt(3);

union Supertriangle {
  double verts[6];
  struct {
    double Ax;
    double Ay;
    double Bx;
    double By;
    double Cx;
    double Cy;
  };
};

// outputs the vertices of the supertriangle to the supertriangle matrix.
// this should be a pre-allocated 3x2 matrix
int make_supertriangle(const gsl_matrix *data, gsl_matrix *supertriangle, double delta) {
  if (supertriangle->size1 != 3 && supertriangle->size2 != 2) {
    GSL_ERROR("supertriangle should be 3x2", GSL_EINVAL);
    return -1;
  }

  double xmin = gsl_matrix_get(data, 0, 0);
  double xmax = xmin;
  for (int row = 1; row < data->size1; ++row) {
    double v = gsl_matrix_get(data, row, 0);
    if (v < xmin) xmin = v;
    if (v > xmax) xmax = v;
  }

  double ymin = gsl_matrix_get(data, 0, 1);
  double ymax = ymin;
  for (int row = 1; row < data->size1; ++row) {
    double v = gsl_matrix_get(data, row, 1);
    if (v < ymin) ymin = v;
    if (v > ymax) ymax = v;
  }

  // find supertriangle vertices
  Supertriangle S;
  S.Cy = ymin;
  S.By = ymin;
  S.Ax = (xmin + xmax)/2.0;
  S.Bx = xmin - IR3*(ymax - ymin);
  S.Ay = (S.Ax - S.Bx)*R3 + ymin;
  S.Cx = 2*S.Ax - S.Bx;

  // add delta
  S.Ay += delta;
  S.Bx -= R3*delta/2;
  S.By -= delta/2;
  S.Cx += R3*delta/2;
  S.Cy -= delta/2;

  // make a view into the supertriangle array and copy it onto the
  // output
  gsl_matrix_const_view S_view = gsl_matrix_const_view_array(S.verts, 3, 2);
  gsl_matrix_memcpy(supertriangle, &S_view.matrix);
  return 0;
}

// calculates the circumcircle for the 3 specified points, and outputs
// the radius to r and the centre to X.
// points should be a 3x2 matrix
// X should be a pre-allocated vector of length 2.
int find_circumcircle(const gsl_matrix *points, double &r, gsl_vector *X) {
  if (points->size1 != 3 && points->size2 != 2) {
    GSL_ERROR("points should be 3x2", GSL_EINVAL);
    return -1;
  }
  if (X->size != 2) {
    GSL_ERROR("X should be of length 2", GSL_EINVAL);
    return -2;
  }

  const double x1 = gsl_matrix_get(points, 0, 0);
  const double x2 = gsl_matrix_get(points, 1, 0);
  const double x3 = gsl_matrix_get(points, 2, 0);
  const double y1 = gsl_matrix_get(points, 0, 1);
  const double y2 = gsl_matrix_get(points, 1, 1);
  const double y3 = gsl_matrix_get(points, 2, 1);

  // calculate the norm^2 terms (these can then be reused)
  const double n1 = x1*x1 + y1*y1;
  const double n2 = x2*x2 + y2*y2;
  const double n3 = x3*x3 + y3*y3;

  // calculate the determinants
  const double a = x1*(y2 - y3) - y1*(x2 - x3) + (x2*y3 - x3*y2);
  const double bx = -( n1*(y2 - y3) - y1*(n2 - n3) + (n2*y3 - n3*y2) );
  const double by = n1*(x2 - x3) - x1*(n2 - n3) + (n2*x3 - n3*x2);
  const double c = -( n1*(x2*y3 - x3*y2) - x1*(n2*y3 - n3*y2) + y1*(n2*x3 - n3*x2) );

  r = sqrt(bx*bx + by*by - 4*a*c)/(2*abs(a));
  gsl_vector_set(X, 0, (-bx)/(2*a));
  gsl_vector_set(X, 1, (-by)/(2*a));

  return 0;
}

bool is_in_circumcircle(const gsl_vector *point, const double r, const gsl_vector *X) {
  double p[2];
  gsl_vector_view p_view = gsl_vector_view_array(p, 2);
  gsl_vector_memcpy(&p_view.vector, point);
  gsl_vector_sub(&p_view.vector, X);
  return gsl_blas_dnrm2(&p_view.vector) <= r;
}

bool is_in_circumcircle(const gsl_vector *point, const gsl_matrix *verts) {
  double r;
  double X[2];
  gsl_vector_view X_view = gsl_vector_view_array(X, 2);
  find_circumcircle(verts, r, &X_view.vector);
  return is_in_circumcircle(point, r, &X_view.vector);
}

/*
  Define everything we need to use Edge in an unordered set
*/

bool Edge::operator==(const Edge b) const {
  return ((A == b.A && B == b.B) || (A == b.B && B == b.A));
}

// computes a hash that's independent of the order the vertices are
// specified in. This allows Edge(a, b) to be recognised as the same
// thing as Edge(b, a).
std::size_t hash_value(Edge const &p) {
  size_t seed = 0;
  if (p.A > p.B) {
    boost::hash_combine(seed, p.A);
    boost::hash_combine(seed, p.B);
    return seed;
  }
  boost::hash_combine(seed, p.B);
  boost::hash_combine(seed, p.A);
  return seed;
}

/*
  ===============
*/


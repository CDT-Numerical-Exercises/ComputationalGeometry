#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H 1

#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double calculate_D(const gsl_vector *p, const gsl_vector *q, const gsl_vector *r);
double calculate_D(const gsl_matrix *X, const size_t pi, const size_t qi, const size_t ri);

std::vector<size_t> convex_hull(const gsl_matrix *verts);

#endif

#ifndef DELAUNAY_H
#define DELAUNAY_H 1

#include <gsl/gsl_matrix.h>

int make_supertriangle(const gsl_matrix *data, gsl_matrix *supertriangle, double delta = 1);
int find_circumcircle(const gsl_matrix *points, double &r, gsl_vector *X);

#endif

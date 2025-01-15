#ifndef DELAUNAY_H
#define DELAUNAY_H 1

#include <gsl/gsl_matrix.h>

int make_supertriangle(const gsl_matrix *data, gsl_matrix *supertriangle);

#endif

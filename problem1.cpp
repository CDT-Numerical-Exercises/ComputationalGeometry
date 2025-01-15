#include <iostream>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gnuplot-iostream/gnuplot-iostream.h>

#include "helpers.h"
#include "delaunay.h"

const double data[4*2] = {
  0, 0,
  1, 3,
  2, 1,
  5, 4
};

int main() {
  const gsl_matrix_const_view d_view = gsl_matrix_const_view_array(data, 4, 2);
  double st[6];
  gsl_matrix_view supertriangle_view = gsl_matrix_view_array(st, 3, 2);
  make_supertriangle(&d_view.matrix, &supertriangle_view.matrix);
  print_matrix(supertriangle_view);

  Gnuplot gp;
  gp << "set yrange[-1:10]\nset xrange[-3:8]\n";
  gp << "plot '-' with circles, '-' with lines\n";
  for (int row = 0; row < d_view.matrix.size1; ++row) {
    for (int col = 0; col < d_view.matrix.size2 - 1; ++col) {
      gp << gsl_matrix_get(&d_view.matrix, row, col) << " ";
    }
    gp << gsl_matrix_get(&d_view.matrix, row, d_view.matrix.size2 - 1);
    gp << "\n";
  }
  gp << "e\n";
  
  for (int row = 0; row < supertriangle_view.matrix.size1; ++row) {
    for (int col = 0; col < supertriangle_view.matrix.size2 - 1; ++col) {
      gp << gsl_matrix_get(&supertriangle_view.matrix, row, col) << " ";
    }
    gp << gsl_matrix_get(&supertriangle_view.matrix, row, supertriangle_view.matrix.size2 - 1);
    gp << "\n";
  }
  gp << gsl_matrix_get(&supertriangle_view.matrix, 0, 0) << " " << gsl_matrix_get(&supertriangle_view.matrix, 0, 1) << "\n";
  gp << "e\n";
}

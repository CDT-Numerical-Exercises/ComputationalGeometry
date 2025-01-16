#include <iostream>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gnuplot-iostream/gnuplot-iostream.h>

#include "helpers.h"
#include "delaunay.h"

const double data[3*2 + 2*2] = {
  0, 0,
  1, 1,
  -0.25, 1.6,
  0, 0,
  1, 1
};

int main() {
  gsl_vector *X = gsl_vector_alloc(2);
  double r;

  // check three different permutations
  {
    const gsl_matrix_const_view d_view = gsl_matrix_const_view_array(data, 3, 2);

    find_circumcircle(&d_view.matrix, r, X);

    std::cout << "Centre: ";
    print_vector(X);
    std::cout << "r: " << r << std::endl;
  }
  {
    const gsl_matrix_const_view d_view = gsl_matrix_const_view_array(data+2, 3, 2);

    find_circumcircle(&d_view.matrix, r, X);

    std::cout << "Centre: ";
    print_vector(X);
    std::cout << "r: " << r << std::endl;
  }
  {
    const gsl_matrix_const_view d_view = gsl_matrix_const_view_array(data+4, 3, 2);

    find_circumcircle(&d_view.matrix, r, X);

    std::cout << "Centre: ";
    print_vector(X);
    std::cout << "r: " << r << std::endl;
  }

  // plot it
  const gsl_matrix_const_view d_view = gsl_matrix_const_view_array(data, 3, 2);
  Gnuplot gp;
  gp << "set size square\n";
  gp << "set yrange[-2:2]\nset xrange[-2:2]\n";
  // std::cout << "set object 1 circle at " << gsl_vector_get(X, 0) << "," << gsl_vector_get(X, 1) << " size scr " << r << " fc rgb 'navy'\n" << std::endl;
  gp << "set object 1 circle at " << gsl_vector_get(X, 0) << "," << gsl_vector_get(X, 1) << " size " << r << " fc rgb 'red'\n";
  gp << "plot '-' with circles\n";
  for (int row = 0; row < d_view.matrix.size1; ++row) {
    for (int col = 0; col < d_view.matrix.size2 - 1; ++col) {
      gp << gsl_matrix_get(&d_view.matrix, row, col) << " ";
    }
    gp << gsl_matrix_get(&d_view.matrix, row, d_view.matrix.size2 - 1);
    gp << "\n";
  }
  // add the centre
  gp << gsl_vector_get(X, 0) << " " << gsl_vector_get(X, 1) << "\n";
  gp << "e\n";

  gsl_vector_free(X);
}

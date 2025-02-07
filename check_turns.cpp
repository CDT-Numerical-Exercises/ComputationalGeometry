#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "convex_hull.h"

const double data[9*2] = {
  0, 0,
  0, 1,
  -1, 2,
  1, 2,
  -1, 0.5,
  1, 0.5,

  1, 1,
  1.5, 2,
  2, 1.5
};

int main() {
  const gsl_matrix_const_view d_view = gsl_matrix_const_view_array(data, 9, 2);

  using namespace std;

  cout << "Up and Left: " << calculate_D(&d_view.matrix, 0, 1, 2) << endl;
  cout << "             " << (calculate_D(&d_view.matrix, 0, 1, 2) > 0) << endl;
  cout << "Up and Right: " << calculate_D(&d_view.matrix, 0, 1, 3) << endl;
  cout << "              " << (calculate_D(&d_view.matrix, 0, 1, 3) < 0) << endl;
  cout << "Down and Left: " << calculate_D(&d_view.matrix, 0, 1, 4) << endl;
  cout << "               " << (calculate_D(&d_view.matrix, 0, 1, 4) > 0) << endl;
  cout << "Down and Right: " << calculate_D(&d_view.matrix, 0, 1, 5) << endl;
  cout << "                " << (calculate_D(&d_view.matrix, 0, 1, 5) < 0) << endl;

  cout << "Example 1: " << calculate_D(&d_view.matrix, 0, 6, 7) << endl;
  cout << "           " << (calculate_D(&d_view.matrix, 0, 6, 7) > 0) << endl;
  cout << "Example 2: " << calculate_D(&d_view.matrix, 0, 6, 8) << endl;
  cout << "           " << (calculate_D(&d_view.matrix, 0, 6, 8) < 0) << endl;
}

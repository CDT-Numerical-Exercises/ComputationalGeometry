#define _USE_MATH_DEFINES
#include <cmath>
#include <gsl/gsl_matrix.h>

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
  double ymax = xmin;
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

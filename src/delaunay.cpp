#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <boost/container_hash/hash.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gnuplot-iostream/gnuplot-iostream.h>
#include <filesystem>

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

int get_point(const size_t i, const gsl_matrix *verts, const gsl_matrix *supertriangle,
               gsl_vector *X) {
  const size_t N = verts->size1;
  if (i < N) {
    gsl_vector_set(X, 0, gsl_matrix_get(verts, i, 0));
    gsl_vector_set(X, 1, gsl_matrix_get(verts, i, 1));
    return 0;
  }
  if (supertriangle == NULL) {
    return -1;
  }
  const size_t j = i - N;
  gsl_vector_set(X, 0, gsl_matrix_get(supertriangle, j, 0));
  gsl_vector_set(X, 1, gsl_matrix_get(supertriangle, j, 1));
  return 0;
}

int Triangle::to_matrix(const gsl_matrix *verts,
                        const gsl_matrix *supertriangle, gsl_matrix *X) const {
  double v[2];
  gsl_vector_view v_view = gsl_vector_view_array(v, 2);

  gsl_vector_view row = gsl_matrix_row(X, 0);
  int ret = get_point(A, verts, supertriangle, &v_view.vector);
  if (ret != 0) return ret;
  gsl_vector_memcpy(&row.vector, &v_view.vector);

  row = gsl_matrix_row(X, 1);
  ret = get_point(B, verts, supertriangle, &v_view.vector);
  if (ret != 0) return ret;
  gsl_vector_memcpy(&row.vector, &v_view.vector);

  row = gsl_matrix_row(X, 2);
  ret = get_point(C, verts, supertriangle, &v_view.vector);
  if (ret != 0) return ret;
  gsl_vector_memcpy(&row.vector, &v_view.vector);

  return 0;
}

int Triangle::to_matrix(const gsl_matrix *verts, gsl_matrix *X) const {
  return Triangle::to_matrix(verts, NULL, X);
}

// basically a shorthand function to let us test whether a point is
// inside the triangle without having to manually call to_matrix first
// also handles the cleanup by stack-allocating everything
bool Triangle::is_in_circumcircle(const gsl_vector *point, const gsl_matrix *verts,
                        const gsl_matrix *supertriangle) const {
  double V[6];
  gsl_matrix_view V_view = gsl_matrix_view_array(V, 3, 2);
  to_matrix(verts, supertriangle, &V_view.matrix);
  return ::is_in_circumcircle(point, &V_view.matrix);
}

void Triangle::to_edges(Edge &a, Edge &b, Edge &c) const {
  a = Edge(A, B);
  b = Edge(B, C);
  c = Edge(C, A);
}

void save_triangulation(const std::filesystem::path fn, const size_t i, const gsl_matrix *verts, const gsl_matrix *supertriangle, const std::vector<Triangle> &triangles) {
  #ifdef ROBUST_BOUNDS
  // figure out the bounds from the supertriangle
  double xmin = gsl_matrix_get(supertriangle, 0, 0);
  double xmax = xmin;
  for (int row = 1; row < supertriangle->size1; ++row) {
    double v = gsl_matrix_get(supertriangle, row, 0);
    if (v < xmin) xmin = v;
    if (v > xmax) xmax = v;
  }

  double ymin = gsl_matrix_get(supertriangle, 0, 1);
  double ymax = ymin;
  for (int row = 1; row < supertriangle->size1; ++row) {
    double v = gsl_matrix_get(supertriangle, row, 1);
    if (v < ymin) ymin = v;
    if (v > ymax) ymax = v;
  }
  #else
  // faster way of getting the bounds is to use the known structure of
  // the supertriangle matrix
  // this is a little less robusts, but should be a bit faster
  const double xmin = gsl_matrix_get(supertriangle, 1, 0);
  const double xmax = gsl_matrix_get(supertriangle, 2, 0);
  const double ymin = gsl_matrix_get(supertriangle, 1, 1);
  const double ymax = gsl_matrix_get(supertriangle, 0, 1);
  #endif
  
  Gnuplot gp;
  gp << "set terminal pngcairo size 350,262 enhanced font 'Verdana,10'\n";
  gp << "set output " << fn << "\n"; // be warned -- this may not be sanitised
  gp << "set yrange[" << ymin << ":" << ymax << "]\n";
  gp << "set xrange[" << xmin << ":" << xmax << "]\n";

  gp << "set key off\n";
  gp << "plot '-' with circles, "; // for the points

  // add the edges to an unordered set
  // this will make plotting easier
  std::unordered_set<Edge> edges;
  for (const Triangle &tri : triangles) {
    Edge A; Edge B; Edge C;
    tri.to_edges(A, B, C);
    edges.insert(A);
    edges.insert(B);
    edges.insert(C);
  }

  // add a curve for every edge
  const size_t N_edges = edges.size();
  for (size_t i = 1; i < N_edges; ++i) {
    gp << "'-' with lines, ";
  }
  gp << "'-' with lines\n";

  double X[2];
  gsl_vector_view X_view = gsl_vector_view_array(X, 2);
  
  for (int row = 0; row <= i; ++row) {
    get_point(row, verts, supertriangle, &X_view.vector);
    gp << gsl_vector_get(&X_view.vector, 0) << " ";
    gp << gsl_vector_get(&X_view.vector, 1) << "\n";
  }
  gp << "e\n";

  // draw the lines
  for (const Edge &e : edges) {
    get_point(e.A, verts, supertriangle, &X_view.vector);
    gp << gsl_vector_get(&X_view.vector, 0) << " ";
    gp << gsl_vector_get(&X_view.vector, 1) << "\n";
    get_point(e.B, verts, supertriangle, &X_view.vector);
    gp << gsl_vector_get(&X_view.vector, 0) << " ";
    gp << gsl_vector_get(&X_view.vector, 1) << "\ne\n";
  }
}

// Takes a gsl matrix of vertices (must be Nx2) and determines the
// Delaunay triangulation of these vertices.
// Returns a C++ vector containing the indices of the vertices forming
// each triangle. The length of the vector will be a multiple of
// three, and each number corresponds to a row of the input
// matrix. Each group of three numbers in the vector points to three
// vertices forming a triangle in the Delaunay triangulation.
std::vector<Triangle>
delaunay_triangulate(const gsl_matrix *verts) {
  std::vector<Triangle> triangles;

  double supertriangle[6];
  gsl_matrix_view st_view = gsl_matrix_view_array(supertriangle, 3, 2);
  make_supertriangle(verts, &st_view.matrix);

  // add the supertriangle to the triangles list
  const double N = verts->size1;
  triangles.push_back(Triangle(N, N+1, N+2));

  // go through all points, including the supertriangle points
  double point[2];
  gsl_vector_view p_view = gsl_vector_view_array(point, 2);
  gsl_vector *p = &p_view.vector; // not entirely safe, but p_view will never go out of scope
  size_t frame = 0;
  for (size_t i = 0; i < (N+3); ++i) {
    // save the triangulation in its current state
    char buf[50];
    snprintf(buf, sizeof(buf), "frame%03lu.png", frame);
    save_triangulation(buf, i, verts, &st_view.matrix, triangles);
    ++frame;
    
    // retrieve the point
    get_point(i, verts, &st_view.matrix, p);
    
    std::unordered_map<Edge,size_t> edges;

    // find all triangles the point lies inside, and add their edges
    // to the edges set
    std::vector<size_t> to_remove;
    for (size_t j = 0; j < triangles.size(); ++j) {
      const Triangle &tri = triangles.at(j);
      if (tri.is_in_circumcircle(p, verts, &st_view.matrix)) {
        Edge A; Edge B; Edge C;
        tri.to_edges(A, B, C);
        // this construction is safe
        // if the key does not exist, it will be initialised to 0 first
        // https://stackoverflow.com/a/48844516
        edges[A] += 1;
        edges[B] += 1;
        edges[C] += 1;
        // mark this triangle for removal
        to_remove.insert(to_remove.begin(), j);
      }
    }

    // remove all triangles
    // we've added them to the front, so we can safely just iterate
    for (size_t j : to_remove) {
      auto it = triangles.begin();
      std::advance(it, j);
      triangles.erase(it);
    }

    snprintf(buf, sizeof(buf), "frame%03lu.png", frame);
    save_triangulation(buf, i, verts, &st_view.matrix, triangles);
    ++frame;

    // construct new triangles from each edge + the point
    for (auto edge : edges) {
      if (edge.second == 1) {
        const Triangle tri(edge.first.A, edge.first.B, i);
        triangles.push_back(tri);
      }
    }

    snprintf(buf, sizeof(buf), "frame%03lu.png", frame);
    save_triangulation(buf, i, verts, &st_view.matrix, triangles);
    ++frame;
  }

  // drop any triangles that contain the supertriangle vertices
  std::vector<size_t> to_remove;
  for (size_t j = 0; j < triangles.size(); ++j) {
    const Triangle &tri = triangles.at(j);
    if (tri.A >= N || tri.B >= N || tri.C >= N) {
      to_remove.insert(to_remove.begin(), j);
    }
  }
  for (size_t j : to_remove) {
    auto it = triangles.begin();
    std::advance(it, j);
    triangles.erase(it);
  }

  char buf[50];
  snprintf(buf, sizeof(buf), "frame%03lu.png", frame);
  save_triangulation(buf, N-1, verts, &st_view.matrix, triangles);

  return triangles;
}

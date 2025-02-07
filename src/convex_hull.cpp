#include <vector>
#include <algorithm>
#include <filesystem>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gnuplot-iostream/gnuplot-iostream.h>

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

std::vector<size_t> argsort_by_x(const gsl_matrix *verts) {
  // create a vector of indices
  std::vector<size_t> args;
  for (size_t i = 0; i < verts->size1; ++i) {
    args.push_back(i);
  }

  // sort the indices based on the x coordinate of the corresponding
  // vertex
  std::sort(args.begin(), args.end(), [verts](const size_t &a, const size_t &b){
    return gsl_matrix_get(verts, a, 0) < gsl_matrix_get(verts, b, 0); });

  return args;
}

void save_convex_hull(const std::filesystem::path fn, const gsl_matrix *verts,
                      const std::vector<size_t> Lupper) {
  Gnuplot gp;
  gp << "set terminal pngcairo size 350,262 enhanced font 'Verdana,10'\n";
  gp << "set output " << fn << "\n"; // be warned -- this may not be sanitised
  gp << "set key off\n";
  gp << "plot '-' with circles, "; // for the points
  gp << "'-' with lines\n"; // for the hull

  // plot the points
  for (int row = 0; row < verts->size1; ++row) {
    for (int col = 0; col < verts->size2 - 1; ++col) {
      gp << gsl_matrix_get(verts, row, col) << " ";
    }
    gp << gsl_matrix_get(verts, row, verts->size2 - 1);
    gp << "\n";
  }
  gp << "e\n";

  // draw the hull
  for (const size_t &i : Lupper) {
    gsl_vector_const_view row_view = gsl_matrix_const_row(verts, i);
    gp << gsl_vector_get(&row_view.vector, 0) << " " << gsl_vector_get(&row_view.vector, 1) << "\n";
  }
  gp << "e\n";
}

void save_convex_hull(const std::filesystem::path fn, const gsl_matrix *verts,
                      const std::vector<size_t> Lupper,
                      const std::vector<size_t> Llower) {
  std::vector<size_t> L;
  L.reserve(Lupper.size() + Llower.size());
  L.insert(L.end(), Lupper.begin(), Lupper.end());
  L.insert(L.end(), Llower.begin(), Llower.end());
  save_convex_hull(fn, verts, L);
}

// returns a vector of integers corresponding to the indices of the
// vertices in the input matrix. These are in order to create a
// clockwise convex hull.
template <bool save_anim>
std::vector<size_t> convex_hull(const gsl_matrix *verts, const std::filesystem::path frame_dir) {
  char buf[50];
  size_t frame = 0;
  
  std::vector<size_t> points = argsort_by_x(verts);
  const size_t n_points = points.size();

  // iterate starting from i = 2
  std::vector<size_t> Lupper = { points[0], points[1] };
  if constexpr (save_anim) {
    snprintf(buf, sizeof(buf), "frame%03lu.png", frame);
    save_convex_hull(frame_dir / buf, verts, Lupper);
    ++frame;
  }
  for (size_t i = 2; i < points.size(); ++i) {
    Lupper.push_back(points[i]);

    if constexpr (save_anim) {
      snprintf(buf, sizeof(buf), "frame%03lu.png", frame);
      save_convex_hull(frame_dir / buf, verts, Lupper);
      ++frame;
    }

    // remove points until the last turn is a right turn, or we don't
    // have enough points to define a right turn.
    // right turn is when D < 0
    size_t current_length = Lupper.size();
    while (current_length > 2 && calculate_D(verts, Lupper[current_length - 3], Lupper[current_length - 2], Lupper[current_length - 1]) > 0) {
      // remove the middle of the last 3 points
      // i.e. the second to last point
      --current_length;
      auto it = std::prev(std::prev(Lupper.end()));
      Lupper.erase(it);
    }
  }

  if constexpr (save_anim) {
    snprintf(buf, sizeof(buf), "frame%03lu.png", frame);
    save_convex_hull(frame_dir / buf, verts, Lupper);
    ++frame;
  }

  // iterate backwards
  std::vector<size_t> Llower = { points[n_points - 1], points[n_points - 2] };
  if constexpr (save_anim) {
    snprintf(buf, sizeof(buf), "frame%03lu.png", frame);
    save_convex_hull(frame_dir / buf, verts, Lupper, Llower);
    ++frame;
  }
  for (size_t i = n_points - 2; i > 0; --i) {
    Llower.push_back(points[i-1]);

    if constexpr (save_anim) {
      snprintf(buf, sizeof(buf), "frame%03lu.png", frame);
      save_convex_hull(frame_dir / buf, verts, Lupper, Llower);
      ++frame;
    }

    size_t current_length = Llower.size();
    while (current_length > 2 && calculate_D(verts, Llower[current_length - 3], Llower[current_length - 2], Llower[current_length - 1]) > 0) {
      // remove the middle of the last 3 points
      // i.e. the second to last point
      --current_length;
      auto it = std::prev(std::prev(Llower.end()));
      Llower.erase(it);
    }
  }

  if constexpr (save_anim) {
    snprintf(buf, sizeof(buf), "frame%03lu.png", frame);
    save_convex_hull(frame_dir / buf, verts, Lupper, Llower);
    ++frame;
  }

  // remove the first and last points
  Llower.erase(Llower.begin());
  Llower.erase(std::prev(Llower.end()));

  // glue the vectors together
  Lupper.reserve(Lupper.size() + Llower.size());
  Lupper.insert(Lupper.end(), Llower.begin(), Llower.end());

  if constexpr (save_anim) {
    std::vector<size_t> L(Lupper);
    L.push_back(L[0]);
    snprintf(buf, sizeof(buf), "frame%03lu.png", frame);
    save_convex_hull(frame_dir / buf, verts, L);
    ++frame;
  }

  return Lupper;
}

std::vector<size_t> convex_hull(const gsl_matrix *verts) {
  return convex_hull<false>(verts, "");
}

std::vector<size_t> convex_hull(const gsl_matrix *verts, const std::filesystem::path frame_dir) {
  // create the directory if it doesn't exist
  if (!std::filesystem::exists(frame_dir)) {
    std::filesystem::create_directory(frame_dir);
  }
  return convex_hull<true>(verts, frame_dir);
}

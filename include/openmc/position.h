#ifndef OPENMC_POSITION_H
#define OPENMC_POSITION_H

#include <cmath> // for sqrt
#include <iostream>
#include <stdexcept> // for out_of_range

#include "openmc/array.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Type representing a position in Cartesian coordinates
//==============================================================================

typedef double double4_t __attribute__ ((vector_size (4 * sizeof(double))));

struct Position {
  // Constructors
  Position() = default;
  Position(double4_t r_) : r {r_} {}; 
  Position(double x_, double y_, double z_) : r {x_, y_, z_, 0.0} {};
  // Position(double x_, double y_, double z_) : x {x_}, y {y_}, z {z_} {};
  Position(const double xyz[]) : r {xyz[0], xyz[1], xyz[2], 0.0} {};
  // Position(const double xyz[]) : x {xyz[0]}, y {xyz[1]}, z {xyz[2]} {};
  Position(const vector<double>& xyz) : r {xyz[0], xyz[1], xyz[2], 0.0} {};
  // Position(const vector<double>& xyz) : x {xyz[0]}, y {xyz[1]}, z {xyz[2]} {};
  Position(const array<double, 3>& xyz) : r {xyz[0], xyz[1], xyz[2]} {};
  // Position(const array<double, 3>& xyz) : x {xyz[0]}, y {xyz[1]}, z {xyz[2]} {};

  // Unary operators
  Position& operator+=(Position);
  Position& operator+=(double);
  Position& operator-=(Position);
  Position& operator-=(double);
  Position& operator*=(Position);
  Position& operator*=(double);
  Position& operator/=(Position);
  Position& operator/=(double);
  Position operator-() const;

  const double& operator[](int i) const
  {
      return r[i];
    // switch (i) {
    // case 0:
    //   return x;
    // case 1:
    //   return y;
    // case 2:
    //   return z;
    // default:
    //   throw std::out_of_range {"Index in Position must be between 0 and 2."};
    // }
  }
  double& operator[](int i)
  {
    return r[i];
    // switch (i) {
    // case 0:
    //   return x;
    // case 1:
    //   return y;
    // case 2:
    //   return z;
    // default:
    //   throw std::out_of_range {"Index in Position must be between 0 and 2."};
    // }
  }

  // Access to x, y, or z by compile time known index (specializations below)
  template<int i>
  const double& get() const
  {
    throw std::out_of_range {"Index in Position must be between 0 and 2."};
  }
  template<int i>
  double& get()
  {
    throw std::out_of_range {"Index in Position must be between 0 and 2."};
  }

  // Other member functions

  //! Dot product of two vectors
  //! \param[in] other Vector to take dot product with
  //! \result Resulting dot product
  inline double dot(Position other) const
  {
    double4_t dot_r = r * other.r;
    return dot_r[0] + dot_r[1] + dot_r[2];
    // return x * other.x + y * other.y + z * other.z;
  }
  inline double norm() const { 
      double4_t square = r * r;
      return std::sqrt(square[0] + square[1] + square[2]); 
  }
  // inline double norm() const { return std::sqrt(x * x + y * y + z * z); }

  //! Reflect a direction across a normal vector
  //! \param[in] other Vector to reflect across
  //! \result Reflected vector
  Position reflect(Position n) const;

  //! Rotate the position based on a rotation matrix
  Position rotate(const vector<double>& rotation) const;

  // Data members
  // double x = 0.;
  // double y = 0.;
  // double z = 0.;
  double4_t r = { 0.0, 0.0, 0.0, 0.0 };
};

// Compile-time known member index access functions
template<>
inline const double& Position::get<0>() const
{
  return r[0];
  // return x;
}
template<>
inline const double& Position::get<1>() const
{
  return r[1];
  // return y;
}
template<>
inline const double& Position::get<2>() const
{
  return r[2];
  // return z;
}
template<>
inline double& Position::get<0>()
{
  return x;
}
template<>
inline double& Position::get<1>()
{
  return y;
}
template<>
inline double& Position::get<2>()
{
  return z;
}

// Binary operators
inline Position operator+(Position a, Position b)
{
  return a += b;
}
inline Position operator+(Position a, double b)
{
  return a += b;
}
inline Position operator+(double a, Position b)
{
  return b += a;
}

inline Position operator-(Position a, Position b)
{
  return a -= b;
}
inline Position operator-(Position a, double b)
{
  return a -= b;
}
inline Position operator-(double a, Position b)
{
  return b -= a;
}

inline Position operator*(Position a, Position b)
{
  return a *= b;
}
inline Position operator*(Position a, double b)
{
  return a *= b;
}
inline Position operator*(double a, Position b)
{
  return b *= a;
}

inline Position operator/(Position a, Position b)
{
  return a /= b;
}
inline Position operator/(Position a, double b)
{
  return a /= b;
}
inline Position operator/(double a, Position b)
{
  return b /= a;
}

inline Position Position::reflect(Position n) const
{
  const double projection = n.dot(*this);
  const double magnitude = n.dot(n);
  n *= (2.0 * projection / magnitude);
  return *this - n;
}

inline bool operator==(Position a, Position b)
{
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

inline bool operator!=(Position a, Position b)
{
  return a.x != b.x || a.y != b.y || a.z != b.z;
}

std::ostream& operator<<(std::ostream& os, Position a);

//==============================================================================
//! Type representing a vector direction in Cartesian coordinates
//==============================================================================

using Direction = Position;

} // namespace openmc

#endif // OPENMC_POSITION_H

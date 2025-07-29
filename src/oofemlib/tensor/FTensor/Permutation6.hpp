#pragma once

#include "error_when_instantiated.hpp"

namespace FTensor
{
template <int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
struct Permutation6 {
  Permutation6() = default;

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N2, N3, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, i, j, k, l, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N2, N3, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, i, j, k, m, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N2, N4, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, i, j, k, m, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N2, N4, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, i, j, k, n, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N2, N5, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim2, Dim5, Dim4, Dim3, i, j, k, n, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N2, N5, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N3, N2, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim3, Dim2, Dim5, Dim4, i, j, l, k, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N3, N2, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, i, j, l, m, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N3, N4, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, i, j, l, m, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N3, N4, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, i, j, l, n, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N3, N5, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim3, Dim5, Dim4, Dim2, i, j, l, n, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N3, N5, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N4, N2, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim4, Dim2, Dim5, Dim3, i, j, m, k, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N4, N2, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim4, Dim3, Dim2, Dim5, i, j, m, l, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N4, N3, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim4, Dim3, Dim5, Dim2, i, j, m, l, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N4, N3, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, i, j, m, n, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N4, N5, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim4, Dim5, Dim3, Dim2, i, j, m, n, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N4, N5, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, i, j, n, k, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N5, N2, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim5, Dim2, Dim4, Dim3, i, j, n, k, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N5, N2, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim5, Dim3, Dim2, Dim4, i, j, n, l, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N5, N3, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim5, Dim3, Dim4, Dim2, i, j, n, l, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N5, N3, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim5, Dim4, Dim2, Dim3, i, j, n, m, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N5, N4, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim1, Dim5, Dim4, Dim3, Dim2, i, j, n, m, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N1, N5, N4, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N1, N3, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim1, Dim3, Dim5, Dim4, i, k, j, l, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N1, N3, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim1, Dim4, Dim3, Dim5, i, k, j, m, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N1, N4, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim1, Dim4, Dim5, Dim3, i, k, j, m, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N1, N4, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim1, Dim5, Dim3, Dim4, i, k, j, n, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N1, N5, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim1, Dim5, Dim4, Dim3, i, k, j, n, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N1, N5, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N3, N1, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim3, Dim1, Dim5, Dim4, i, k, l, j, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N3, N1, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, i, k, l, m, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N3, N4, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, i, k, l, m, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N3, N4, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, i, k, l, n, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N3, N5, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim3, Dim5, Dim4, Dim1, i, k, l, n, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N3, N5, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N4, N1, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim4, Dim1, Dim5, Dim3, i, k, m, j, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N4, N1, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim4, Dim3, Dim1, Dim5, i, k, m, l, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N4, N3, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim4, Dim3, Dim5, Dim1, i, k, m, l, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N4, N3, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, i, k, m, n, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N4, N5, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim4, Dim5, Dim3, Dim1, i, k, m, n, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N4, N5, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, i, k, n, j, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N5, N1, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim5, Dim1, Dim4, Dim3, i, k, n, j, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N5, N1, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim5, Dim3, Dim1, Dim4, i, k, n, l, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N5, N3, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim5, Dim3, Dim4, Dim1, i, k, n, l, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N5, N3, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim5, Dim4, Dim1, Dim3, i, k, n, m, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N5, N4, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim2, Dim5, Dim4, Dim3, Dim1, i, k, n, m, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N2, N5, N4, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N1, N2, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim1, Dim2, Dim5, Dim4, i, l, j, k, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N1, N2, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim1, Dim4, Dim2, Dim5, i, l, j, m, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N1, N4, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim1, Dim4, Dim5, Dim2, i, l, j, m, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N1, N4, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim1, Dim5, Dim2, Dim4, i, l, j, n, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N1, N5, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim1, Dim5, Dim4, Dim2, i, l, j, n, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N1, N5, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim2, Dim1, Dim4, Dim5, i, l, k, j, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N2, N1, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim2, Dim1, Dim5, Dim4, i, l, k, j, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N2, N1, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim2, Dim4, Dim1, Dim5, i, l, k, m, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N2, N4, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim2, Dim4, Dim5, Dim1, i, l, k, m, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N2, N4, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim2, Dim5, Dim1, Dim4, i, l, k, n, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N2, N5, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim2, Dim5, Dim4, Dim1, i, l, k, n, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N2, N5, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N4, N1, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim4, Dim1, Dim5, Dim2, i, l, m, j, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N4, N1, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim4, Dim2, Dim1, Dim5, i, l, m, k, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N4, N2, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim4, Dim2, Dim5, Dim1, i, l, m, k, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N4, N2, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, i, l, m, n, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N4, N5, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim4, Dim5, Dim2, Dim1, i, l, m, n, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N4, N5, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, i, l, n, j, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N5, N1, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim5, Dim1, Dim4, Dim2, i, l, n, j, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N5, N1, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim5, Dim2, Dim1, Dim4, i, l, n, k, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N5, N2, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim5, Dim2, Dim4, Dim1, i, l, n, k, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N5, N2, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim5, Dim4, Dim1, Dim2, i, l, n, m, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N5, N4, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim3, Dim5, Dim4, Dim2, Dim1, i, l, n, m, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N3, N5, N4, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N1, N2, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim1, Dim2, Dim5, Dim3, i, m, j, k, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N1, N2, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim1, Dim3, Dim2, Dim5, i, m, j, l, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N1, N3, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim1, Dim3, Dim5, Dim2, i, m, j, l, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N1, N3, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim1, Dim5, Dim2, Dim3, i, m, j, n, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N1, N5, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim1, Dim5, Dim3, Dim2, i, m, j, n, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N1, N5, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim2, Dim1, Dim3, Dim5, i, m, k, j, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N2, N1, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim2, Dim1, Dim5, Dim3, i, m, k, j, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N2, N1, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim2, Dim3, Dim1, Dim5, i, m, k, l, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N2, N3, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim2, Dim3, Dim5, Dim1, i, m, k, l, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N2, N3, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim2, Dim5, Dim1, Dim3, i, m, k, n, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N2, N5, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim2, Dim5, Dim3, Dim1, i, m, k, n, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N2, N5, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim3, Dim1, Dim2, Dim5, i, m, l, j, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N3, N1, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim3, Dim1, Dim5, Dim2, i, m, l, j, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N3, N1, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim3, Dim2, Dim1, Dim5, i, m, l, k, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N3, N2, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim3, Dim2, Dim5, Dim1, i, m, l, k, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N3, N2, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim3, Dim5, Dim1, Dim2, i, m, l, n, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N3, N5, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim3, Dim5, Dim2, Dim1, i, m, l, n, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N3, N5, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, i, m, n, j, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N5, N1, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim5, Dim1, Dim3, Dim2, i, m, n, j, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N5, N1, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim5, Dim2, Dim1, Dim3, i, m, n, k, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N5, N2, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim5, Dim2, Dim3, Dim1, i, m, n, k, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N5, N2, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim5, Dim3, Dim1, Dim2, i, m, n, l, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N5, N3, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim4, Dim5, Dim3, Dim2, Dim1, i, m, n, l, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N4, N5, N3, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N1, N2, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim1, Dim2, Dim4, Dim3, i, n, j, k, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N1, N2, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim1, Dim3, Dim2, Dim4, i, n, j, l, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N1, N3, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim1, Dim3, Dim4, Dim2, i, n, j, l, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N1, N3, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim1, Dim4, Dim2, Dim3, i, n, j, m, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N1, N4, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim1, Dim4, Dim3, Dim2, i, n, j, m, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N1, N4, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim2, Dim1, Dim3, Dim4, i, n, k, j, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N2, N1, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim2, Dim1, Dim4, Dim3, i, n, k, j, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N2, N1, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim2, Dim3, Dim1, Dim4, i, n, k, l, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N2, N3, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim2, Dim3, Dim4, Dim1, i, n, k, l, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N2, N3, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim2, Dim4, Dim1, Dim3, i, n, k, m, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N2, N4, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim2, Dim4, Dim3, Dim1, i, n, k, m, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N2, N4, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim3, Dim1, Dim2, Dim4, i, n, l, j, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N3, N1, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim3, Dim1, Dim4, Dim2, i, n, l, j, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N3, N1, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim3, Dim2, Dim1, Dim4, i, n, l, k, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N3, N2, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim3, Dim2, Dim4, Dim1, i, n, l, k, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N3, N2, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim3, Dim4, Dim1, Dim2, i, n, l, m, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N3, N4, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim3, Dim4, Dim2, Dim1, i, n, l, m, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N3, N4, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim4, Dim1, Dim2, Dim3, i, n, m, j, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N4, N1, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim4, Dim1, Dim3, Dim2, i, n, m, j, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N4, N1, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim4, Dim2, Dim1, Dim3, i, n, m, k, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N4, N2, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim4, Dim2, Dim3, Dim1, i, n, m, k, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N4, N2, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim4, Dim3, Dim1, Dim2, i, n, m, l, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N4, N3, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim0, Dim5, Dim4, Dim3, Dim2, Dim1, i, n, m, l, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N0, N5, N4, N3, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim2, Dim3, Dim4, Dim5, j, i, k, l, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N2, N3, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim2, Dim3, Dim5, Dim4, j, i, k, l, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N2, N3, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim2, Dim4, Dim3, Dim5, j, i, k, m, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N2, N4, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim2, Dim4, Dim5, Dim3, j, i, k, m, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N2, N4, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim2, Dim5, Dim3, Dim4, j, i, k, n, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N2, N5, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim2, Dim5, Dim4, Dim3, j, i, k, n, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N2, N5, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim3, Dim2, Dim4, Dim5, j, i, l, k, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N3, N2, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim3, Dim2, Dim5, Dim4, j, i, l, k, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N3, N2, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim3, Dim4, Dim2, Dim5, j, i, l, m, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N3, N4, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim3, Dim4, Dim5, Dim2, j, i, l, m, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N3, N4, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim3, Dim5, Dim2, Dim4, j, i, l, n, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N3, N5, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim3, Dim5, Dim4, Dim2, j, i, l, n, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N3, N5, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim4, Dim2, Dim3, Dim5, j, i, m, k, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N4, N2, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim4, Dim2, Dim5, Dim3, j, i, m, k, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N4, N2, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim4, Dim3, Dim2, Dim5, j, i, m, l, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N4, N3, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim4, Dim3, Dim5, Dim2, j, i, m, l, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N4, N3, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim4, Dim5, Dim2, Dim3, j, i, m, n, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N4, N5, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim4, Dim5, Dim3, Dim2, j, i, m, n, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N4, N5, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim5, Dim2, Dim3, Dim4, j, i, n, k, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N5, N2, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim5, Dim2, Dim4, Dim3, j, i, n, k, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N5, N2, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim5, Dim3, Dim2, Dim4, j, i, n, l, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N5, N3, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim5, Dim3, Dim4, Dim2, j, i, n, l, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N5, N3, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim5, Dim4, Dim2, Dim3, j, i, n, m, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N5, N4, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim0, Dim5, Dim4, Dim3, Dim2, j, i, n, m, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N0, N5, N4, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N0, N3, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim0, Dim3, Dim5, Dim4, j, k, i, l, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N0, N3, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim0, Dim4, Dim3, Dim5, j, k, i, m, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N0, N4, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim0, Dim4, Dim5, Dim3, j, k, i, m, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N0, N4, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim0, Dim5, Dim3, Dim4, j, k, i, n, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N0, N5, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim0, Dim5, Dim4, Dim3, j, k, i, n, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N0, N5, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N3, N0, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim3, Dim0, Dim5, Dim4, j, k, l, i, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N3, N0, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, j, k, l, m, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N3, N4, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, j, k, l, m, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N3, N4, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, j, k, l, n, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N3, N5, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim3, Dim5, Dim4, Dim0, j, k, l, n, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N3, N5, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N4, N0, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim4, Dim0, Dim5, Dim3, j, k, m, i, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N4, N0, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim4, Dim3, Dim0, Dim5, j, k, m, l, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N4, N3, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim4, Dim3, Dim5, Dim0, j, k, m, l, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N4, N3, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, j, k, m, n, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N4, N5, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim4, Dim5, Dim3, Dim0, j, k, m, n, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N4, N5, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, j, k, n, i, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N5, N0, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim5, Dim0, Dim4, Dim3, j, k, n, i, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N5, N0, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim5, Dim3, Dim0, Dim4, j, k, n, l, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N5, N3, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim5, Dim3, Dim4, Dim0, j, k, n, l, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N5, N3, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim5, Dim4, Dim0, Dim3, j, k, n, m, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N5, N4, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim2, Dim5, Dim4, Dim3, Dim0, j, k, n, m, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N2, N5, N4, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N0, N2, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim0, Dim2, Dim5, Dim4, j, l, i, k, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N0, N2, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim0, Dim4, Dim2, Dim5, j, l, i, m, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N0, N4, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim0, Dim4, Dim5, Dim2, j, l, i, m, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N0, N4, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim0, Dim5, Dim2, Dim4, j, l, i, n, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N0, N5, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim0, Dim5, Dim4, Dim2, j, l, i, n, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N0, N5, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim2, Dim0, Dim4, Dim5, j, l, k, i, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N2, N0, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim2, Dim0, Dim5, Dim4, j, l, k, i, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N2, N0, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim2, Dim4, Dim0, Dim5, j, l, k, m, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N2, N4, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim2, Dim4, Dim5, Dim0, j, l, k, m, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N2, N4, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim2, Dim5, Dim0, Dim4, j, l, k, n, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N2, N5, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim2, Dim5, Dim4, Dim0, j, l, k, n, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N2, N5, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N4, N0, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim4, Dim0, Dim5, Dim2, j, l, m, i, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N4, N0, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim4, Dim2, Dim0, Dim5, j, l, m, k, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N4, N2, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim4, Dim2, Dim5, Dim0, j, l, m, k, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N4, N2, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, j, l, m, n, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N4, N5, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim4, Dim5, Dim2, Dim0, j, l, m, n, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N4, N5, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, j, l, n, i, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N5, N0, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim5, Dim0, Dim4, Dim2, j, l, n, i, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N5, N0, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim5, Dim2, Dim0, Dim4, j, l, n, k, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N5, N2, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim5, Dim2, Dim4, Dim0, j, l, n, k, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N5, N2, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim5, Dim4, Dim0, Dim2, j, l, n, m, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N5, N4, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim3, Dim5, Dim4, Dim2, Dim0, j, l, n, m, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N3, N5, N4, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N0, N2, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim0, Dim2, Dim5, Dim3, j, m, i, k, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N0, N2, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim0, Dim3, Dim2, Dim5, j, m, i, l, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N0, N3, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim0, Dim3, Dim5, Dim2, j, m, i, l, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N0, N3, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim0, Dim5, Dim2, Dim3, j, m, i, n, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N0, N5, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim0, Dim5, Dim3, Dim2, j, m, i, n, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N0, N5, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim2, Dim0, Dim3, Dim5, j, m, k, i, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N2, N0, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim2, Dim0, Dim5, Dim3, j, m, k, i, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N2, N0, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim2, Dim3, Dim0, Dim5, j, m, k, l, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N2, N3, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim2, Dim3, Dim5, Dim0, j, m, k, l, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N2, N3, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim2, Dim5, Dim0, Dim3, j, m, k, n, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N2, N5, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim2, Dim5, Dim3, Dim0, j, m, k, n, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N2, N5, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim3, Dim0, Dim2, Dim5, j, m, l, i, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N3, N0, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim3, Dim0, Dim5, Dim2, j, m, l, i, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N3, N0, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim3, Dim2, Dim0, Dim5, j, m, l, k, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N3, N2, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim3, Dim2, Dim5, Dim0, j, m, l, k, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N3, N2, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim3, Dim5, Dim0, Dim2, j, m, l, n, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N3, N5, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim3, Dim5, Dim2, Dim0, j, m, l, n, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N3, N5, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, j, m, n, i, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N5, N0, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim5, Dim0, Dim3, Dim2, j, m, n, i, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N5, N0, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim5, Dim2, Dim0, Dim3, j, m, n, k, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N5, N2, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim5, Dim2, Dim3, Dim0, j, m, n, k, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N5, N2, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim5, Dim3, Dim0, Dim2, j, m, n, l, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N5, N3, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim4, Dim5, Dim3, Dim2, Dim0, j, m, n, l, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N4, N5, N3, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N0, N2, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim0, Dim2, Dim4, Dim3, j, n, i, k, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N0, N2, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim0, Dim3, Dim2, Dim4, j, n, i, l, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N0, N3, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim0, Dim3, Dim4, Dim2, j, n, i, l, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N0, N3, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim0, Dim4, Dim2, Dim3, j, n, i, m, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N0, N4, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim0, Dim4, Dim3, Dim2, j, n, i, m, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N0, N4, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim2, Dim0, Dim3, Dim4, j, n, k, i, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N2, N0, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim2, Dim0, Dim4, Dim3, j, n, k, i, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N2, N0, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim2, Dim3, Dim0, Dim4, j, n, k, l, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N2, N3, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim2, Dim3, Dim4, Dim0, j, n, k, l, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N2, N3, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim2, Dim4, Dim0, Dim3, j, n, k, m, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N2, N4, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim2, Dim4, Dim3, Dim0, j, n, k, m, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N2, N4, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim3, Dim0, Dim2, Dim4, j, n, l, i, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N3, N0, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim3, Dim0, Dim4, Dim2, j, n, l, i, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N3, N0, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim3, Dim2, Dim0, Dim4, j, n, l, k, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N3, N2, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim3, Dim2, Dim4, Dim0, j, n, l, k, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N3, N2, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim3, Dim4, Dim0, Dim2, j, n, l, m, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N3, N4, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim3, Dim4, Dim2, Dim0, j, n, l, m, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N3, N4, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim4, Dim0, Dim2, Dim3, j, n, m, i, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N4, N0, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim4, Dim0, Dim3, Dim2, j, n, m, i, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N4, N0, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim4, Dim2, Dim0, Dim3, j, n, m, k, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N4, N2, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim4, Dim2, Dim3, Dim0, j, n, m, k, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N4, N2, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim4, Dim3, Dim0, Dim2, j, n, m, l, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N4, N3, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim1, Dim5, Dim4, Dim3, Dim2, Dim0, j, n, m, l, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N1, N5, N4, N3, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim1, Dim3, Dim4, Dim5, k, i, j, l, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N1, N3, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim1, Dim3, Dim5, Dim4, k, i, j, l, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N1, N3, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim1, Dim4, Dim3, Dim5, k, i, j, m, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N1, N4, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim1, Dim4, Dim5, Dim3, k, i, j, m, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N1, N4, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim1, Dim5, Dim3, Dim4, k, i, j, n, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N1, N5, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim1, Dim5, Dim4, Dim3, k, i, j, n, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N1, N5, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim3, Dim1, Dim4, Dim5, k, i, l, j, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N3, N1, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim3, Dim1, Dim5, Dim4, k, i, l, j, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N3, N1, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim3, Dim4, Dim1, Dim5, k, i, l, m, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N3, N4, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim3, Dim4, Dim5, Dim1, k, i, l, m, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N3, N4, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim3, Dim5, Dim1, Dim4, k, i, l, n, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N3, N5, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim3, Dim5, Dim4, Dim1, k, i, l, n, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N3, N5, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim4, Dim1, Dim3, Dim5, k, i, m, j, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N4, N1, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim4, Dim1, Dim5, Dim3, k, i, m, j, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N4, N1, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim4, Dim3, Dim1, Dim5, k, i, m, l, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N4, N3, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim4, Dim3, Dim5, Dim1, k, i, m, l, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N4, N3, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim4, Dim5, Dim1, Dim3, k, i, m, n, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N4, N5, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim4, Dim5, Dim3, Dim1, k, i, m, n, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N4, N5, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim5, Dim1, Dim3, Dim4, k, i, n, j, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N5, N1, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim5, Dim1, Dim4, Dim3, k, i, n, j, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N5, N1, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim5, Dim3, Dim1, Dim4, k, i, n, l, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N5, N3, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim5, Dim3, Dim4, Dim1, k, i, n, l, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N5, N3, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim5, Dim4, Dim1, Dim3, k, i, n, m, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N5, N4, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim0, Dim5, Dim4, Dim3, Dim1, k, i, n, m, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N0, N5, N4, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim0, Dim3, Dim4, Dim5, k, j, i, l, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N0, N3, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim0, Dim3, Dim5, Dim4, k, j, i, l, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N0, N3, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim0, Dim4, Dim3, Dim5, k, j, i, m, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N0, N4, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim0, Dim4, Dim5, Dim3, k, j, i, m, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N0, N4, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim0, Dim5, Dim3, Dim4, k, j, i, n, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N0, N5, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim0, Dim5, Dim4, Dim3, k, j, i, n, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N0, N5, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim3, Dim0, Dim4, Dim5, k, j, l, i, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N3, N0, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim3, Dim0, Dim5, Dim4, k, j, l, i, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N3, N0, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim3, Dim4, Dim0, Dim5, k, j, l, m, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N3, N4, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim3, Dim4, Dim5, Dim0, k, j, l, m, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N3, N4, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim3, Dim5, Dim0, Dim4, k, j, l, n, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N3, N5, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim3, Dim5, Dim4, Dim0, k, j, l, n, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N3, N5, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim4, Dim0, Dim3, Dim5, k, j, m, i, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N4, N0, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim4, Dim0, Dim5, Dim3, k, j, m, i, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N4, N0, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim4, Dim3, Dim0, Dim5, k, j, m, l, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N4, N3, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim4, Dim3, Dim5, Dim0, k, j, m, l, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N4, N3, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim4, Dim5, Dim0, Dim3, k, j, m, n, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N4, N5, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim4, Dim5, Dim3, Dim0, k, j, m, n, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N4, N5, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim5, Dim0, Dim3, Dim4, k, j, n, i, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N5, N0, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim5, Dim0, Dim4, Dim3, k, j, n, i, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N5, N0, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim5, Dim3, Dim0, Dim4, k, j, n, l, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N5, N3, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim5, Dim3, Dim4, Dim0, k, j, n, l, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N5, N3, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim5, Dim4, Dim0, Dim3, k, j, n, m, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N5, N4, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim1, Dim5, Dim4, Dim3, Dim0, k, j, n, m, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N1, N5, N4, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N0, N1, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim0, Dim1, Dim5, Dim4, k, l, i, j, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N0, N1, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim0, Dim4, Dim1, Dim5, k, l, i, m, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N0, N4, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim0, Dim4, Dim5, Dim1, k, l, i, m, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N0, N4, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim0, Dim5, Dim1, Dim4, k, l, i, n, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N0, N5, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim0, Dim5, Dim4, Dim1, k, l, i, n, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N0, N5, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim1, Dim0, Dim4, Dim5, k, l, j, i, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N1, N0, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim1, Dim0, Dim5, Dim4, k, l, j, i, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N1, N0, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim1, Dim4, Dim0, Dim5, k, l, j, m, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N1, N4, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim1, Dim4, Dim5, Dim0, k, l, j, m, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N1, N4, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim1, Dim5, Dim0, Dim4, k, l, j, n, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N1, N5, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim1, Dim5, Dim4, Dim0, k, l, j, n, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N1, N5, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N4, N0, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, k, l, m, i, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N4, N0, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, k, l, m, j, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N4, N1, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, k, l, m, j, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N4, N1, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, k, l, m, n, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N4, N5, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, k, l, m, n, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N4, N5, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, k, l, n, i, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N5, N0, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim5, Dim0, Dim4, Dim1, k, l, n, i, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N5, N0, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim5, Dim1, Dim0, Dim4, k, l, n, j, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N5, N1, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim5, Dim1, Dim4, Dim0, k, l, n, j, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N5, N1, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim5, Dim4, Dim0, Dim1, k, l, n, m, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N5, N4, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim3, Dim5, Dim4, Dim1, Dim0, k, l, n, m, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N3, N5, N4, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N0, N1, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim0, Dim1, Dim5, Dim3, k, m, i, j, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N0, N1, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim0, Dim3, Dim1, Dim5, k, m, i, l, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N0, N3, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim0, Dim3, Dim5, Dim1, k, m, i, l, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N0, N3, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim0, Dim5, Dim1, Dim3, k, m, i, n, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N0, N5, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim0, Dim5, Dim3, Dim1, k, m, i, n, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N0, N5, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim1, Dim0, Dim3, Dim5, k, m, j, i, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N1, N0, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim1, Dim0, Dim5, Dim3, k, m, j, i, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N1, N0, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim1, Dim3, Dim0, Dim5, k, m, j, l, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N1, N3, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim1, Dim3, Dim5, Dim0, k, m, j, l, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N1, N3, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim1, Dim5, Dim0, Dim3, k, m, j, n, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N1, N5, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim1, Dim5, Dim3, Dim0, k, m, j, n, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N1, N5, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim3, Dim0, Dim1, Dim5, k, m, l, i, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N3, N0, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim3, Dim0, Dim5, Dim1, k, m, l, i, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N3, N0, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim3, Dim1, Dim0, Dim5, k, m, l, j, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N3, N1, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim3, Dim1, Dim5, Dim0, k, m, l, j, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N3, N1, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim3, Dim5, Dim0, Dim1, k, m, l, n, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N3, N5, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim3, Dim5, Dim1, Dim0, k, m, l, n, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N3, N5, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, k, m, n, i, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N5, N0, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim5, Dim0, Dim3, Dim1, k, m, n, i, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N5, N0, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim5, Dim1, Dim0, Dim3, k, m, n, j, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N5, N1, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim5, Dim1, Dim3, Dim0, k, m, n, j, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N5, N1, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim5, Dim3, Dim0, Dim1, k, m, n, l, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N5, N3, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim4, Dim5, Dim3, Dim1, Dim0, k, m, n, l, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N4, N5, N3, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N0, N1, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim0, Dim1, Dim4, Dim3, k, n, i, j, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N0, N1, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim0, Dim3, Dim1, Dim4, k, n, i, l, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N0, N3, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim0, Dim3, Dim4, Dim1, k, n, i, l, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N0, N3, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim0, Dim4, Dim1, Dim3, k, n, i, m, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N0, N4, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim0, Dim4, Dim3, Dim1, k, n, i, m, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N0, N4, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim1, Dim0, Dim3, Dim4, k, n, j, i, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N1, N0, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim1, Dim0, Dim4, Dim3, k, n, j, i, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N1, N0, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim1, Dim3, Dim0, Dim4, k, n, j, l, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N1, N3, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim1, Dim3, Dim4, Dim0, k, n, j, l, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N1, N3, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim1, Dim4, Dim0, Dim3, k, n, j, m, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N1, N4, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim1, Dim4, Dim3, Dim0, k, n, j, m, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N1, N4, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim3, Dim0, Dim1, Dim4, k, n, l, i, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N3, N0, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim3, Dim0, Dim4, Dim1, k, n, l, i, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N3, N0, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim3, Dim1, Dim0, Dim4, k, n, l, j, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N3, N1, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim3, Dim1, Dim4, Dim0, k, n, l, j, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N3, N1, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim3, Dim4, Dim0, Dim1, k, n, l, m, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N3, N4, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim3, Dim4, Dim1, Dim0, k, n, l, m, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N3, N4, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim4, Dim0, Dim1, Dim3, k, n, m, i, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N4, N0, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim4, Dim0, Dim3, Dim1, k, n, m, i, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N4, N0, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim4, Dim1, Dim0, Dim3, k, n, m, j, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N4, N1, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim4, Dim1, Dim3, Dim0, k, n, m, j, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N4, N1, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim4, Dim3, Dim0, Dim1, k, n, m, l, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N4, N3, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim2, Dim5, Dim4, Dim3, Dim1, Dim0, k, n, m, l, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N2, N5, N4, N3, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim1, Dim2, Dim4, Dim5, l, i, j, k, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N1, N2, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim1, Dim2, Dim5, Dim4, l, i, j, k, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N1, N2, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim1, Dim4, Dim2, Dim5, l, i, j, m, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N1, N4, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim1, Dim4, Dim5, Dim2, l, i, j, m, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N1, N4, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim1, Dim5, Dim2, Dim4, l, i, j, n, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N1, N5, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim1, Dim5, Dim4, Dim2, l, i, j, n, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N1, N5, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim2, Dim1, Dim4, Dim5, l, i, k, j, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N2, N1, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim2, Dim1, Dim5, Dim4, l, i, k, j, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N2, N1, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim2, Dim4, Dim1, Dim5, l, i, k, m, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N2, N4, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim2, Dim4, Dim5, Dim1, l, i, k, m, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N2, N4, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim2, Dim5, Dim1, Dim4, l, i, k, n, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N2, N5, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim2, Dim5, Dim4, Dim1, l, i, k, n, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N2, N5, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim4, Dim1, Dim2, Dim5, l, i, m, j, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N4, N1, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim4, Dim1, Dim5, Dim2, l, i, m, j, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N4, N1, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim4, Dim2, Dim1, Dim5, l, i, m, k, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N4, N2, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim4, Dim2, Dim5, Dim1, l, i, m, k, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N4, N2, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim4, Dim5, Dim1, Dim2, l, i, m, n, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N4, N5, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim4, Dim5, Dim2, Dim1, l, i, m, n, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N4, N5, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim5, Dim1, Dim2, Dim4, l, i, n, j, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N5, N1, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim5, Dim1, Dim4, Dim2, l, i, n, j, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N5, N1, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim5, Dim2, Dim1, Dim4, l, i, n, k, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N5, N2, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim5, Dim2, Dim4, Dim1, l, i, n, k, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N5, N2, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim5, Dim4, Dim1, Dim2, l, i, n, m, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N5, N4, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim0, Dim5, Dim4, Dim2, Dim1, l, i, n, m, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N0, N5, N4, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim0, Dim2, Dim4, Dim5, l, j, i, k, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N0, N2, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim0, Dim2, Dim5, Dim4, l, j, i, k, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N0, N2, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim0, Dim4, Dim2, Dim5, l, j, i, m, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N0, N4, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim0, Dim4, Dim5, Dim2, l, j, i, m, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N0, N4, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim0, Dim5, Dim2, Dim4, l, j, i, n, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N0, N5, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim0, Dim5, Dim4, Dim2, l, j, i, n, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N0, N5, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim2, Dim0, Dim4, Dim5, l, j, k, i, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N2, N0, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim2, Dim0, Dim5, Dim4, l, j, k, i, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N2, N0, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim2, Dim4, Dim0, Dim5, l, j, k, m, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N2, N4, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim2, Dim4, Dim5, Dim0, l, j, k, m, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N2, N4, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim2, Dim5, Dim0, Dim4, l, j, k, n, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N2, N5, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim2, Dim5, Dim4, Dim0, l, j, k, n, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N2, N5, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim4, Dim0, Dim2, Dim5, l, j, m, i, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N4, N0, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim4, Dim0, Dim5, Dim2, l, j, m, i, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N4, N0, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim4, Dim2, Dim0, Dim5, l, j, m, k, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N4, N2, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim4, Dim2, Dim5, Dim0, l, j, m, k, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N4, N2, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim4, Dim5, Dim0, Dim2, l, j, m, n, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N4, N5, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim4, Dim5, Dim2, Dim0, l, j, m, n, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N4, N5, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim5, Dim0, Dim2, Dim4, l, j, n, i, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N5, N0, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim5, Dim0, Dim4, Dim2, l, j, n, i, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N5, N0, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim5, Dim2, Dim0, Dim4, l, j, n, k, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N5, N2, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim5, Dim2, Dim4, Dim0, l, j, n, k, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N5, N2, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim5, Dim4, Dim0, Dim2, l, j, n, m, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N5, N4, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim1, Dim5, Dim4, Dim2, Dim0, l, j, n, m, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N1, N5, N4, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim0, Dim1, Dim4, Dim5, l, k, i, j, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N0, N1, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim0, Dim1, Dim5, Dim4, l, k, i, j, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N0, N1, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim0, Dim4, Dim1, Dim5, l, k, i, m, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N0, N4, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim0, Dim4, Dim5, Dim1, l, k, i, m, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N0, N4, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim0, Dim5, Dim1, Dim4, l, k, i, n, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N0, N5, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim0, Dim5, Dim4, Dim1, l, k, i, n, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N0, N5, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim1, Dim0, Dim4, Dim5, l, k, j, i, m, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N1, N0, N4, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim1, Dim0, Dim5, Dim4, l, k, j, i, n, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N1, N0, N5, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim1, Dim4, Dim0, Dim5, l, k, j, m, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N1, N4, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim1, Dim4, Dim5, Dim0, l, k, j, m, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N1, N4, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim1, Dim5, Dim0, Dim4, l, k, j, n, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N1, N5, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim1, Dim5, Dim4, Dim0, l, k, j, n, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N1, N5, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim4, Dim0, Dim1, Dim5, l, k, m, i, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N4, N0, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim4, Dim0, Dim5, Dim1, l, k, m, i, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N4, N0, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim4, Dim1, Dim0, Dim5, l, k, m, j, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N4, N1, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim4, Dim1, Dim5, Dim0, l, k, m, j, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N4, N1, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim4, Dim5, Dim0, Dim1, l, k, m, n, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N4, N5, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim4, Dim5, Dim1, Dim0, l, k, m, n, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N4, N5, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim5, Dim0, Dim1, Dim4, l, k, n, i, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N5, N0, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim5, Dim0, Dim4, Dim1, l, k, n, i, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N5, N0, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim5, Dim1, Dim0, Dim4, l, k, n, j, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N5, N1, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim5, Dim1, Dim4, Dim0, l, k, n, j, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N5, N1, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim5, Dim4, Dim0, Dim1, l, k, n, m, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N5, N4, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim2, Dim5, Dim4, Dim1, Dim0, l, k, n, m, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N2, N5, N4, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N0, N1, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim0, Dim1, Dim5, Dim2, l, m, i, j, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N0, N1, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim0, Dim2, Dim1, Dim5, l, m, i, k, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N0, N2, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim0, Dim2, Dim5, Dim1, l, m, i, k, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N0, N2, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim0, Dim5, Dim1, Dim2, l, m, i, n, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N0, N5, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim0, Dim5, Dim2, Dim1, l, m, i, n, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N0, N5, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim1, Dim0, Dim2, Dim5, l, m, j, i, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N1, N0, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim1, Dim0, Dim5, Dim2, l, m, j, i, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N1, N0, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim1, Dim2, Dim0, Dim5, l, m, j, k, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N1, N2, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim1, Dim2, Dim5, Dim0, l, m, j, k, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N1, N2, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim1, Dim5, Dim0, Dim2, l, m, j, n, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N1, N5, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim1, Dim5, Dim2, Dim0, l, m, j, n, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N1, N5, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim2, Dim0, Dim1, Dim5, l, m, k, i, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N2, N0, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim2, Dim0, Dim5, Dim1, l, m, k, i, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N2, N0, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim2, Dim1, Dim0, Dim5, l, m, k, j, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N2, N1, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim2, Dim1, Dim5, Dim0, l, m, k, j, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N2, N1, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim2, Dim5, Dim0, Dim1, l, m, k, n, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N2, N5, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim2, Dim5, Dim1, Dim0, l, m, k, n, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N2, N5, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, l, m, n, i, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N5, N0, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim5, Dim0, Dim2, Dim1, l, m, n, i, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N5, N0, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim5, Dim1, Dim0, Dim2, l, m, n, j, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N5, N1, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim5, Dim1, Dim2, Dim0, l, m, n, j, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N5, N1, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim5, Dim2, Dim0, Dim1, l, m, n, k, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N5, N2, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim4, Dim5, Dim2, Dim1, Dim0, l, m, n, k, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N4, N5, N2, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N0, N1, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim0, Dim1, Dim4, Dim2, l, n, i, j, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N0, N1, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim0, Dim2, Dim1, Dim4, l, n, i, k, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N0, N2, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim0, Dim2, Dim4, Dim1, l, n, i, k, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N0, N2, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim0, Dim4, Dim1, Dim2, l, n, i, m, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N0, N4, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim0, Dim4, Dim2, Dim1, l, n, i, m, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N0, N4, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim1, Dim0, Dim2, Dim4, l, n, j, i, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N1, N0, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim1, Dim0, Dim4, Dim2, l, n, j, i, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N1, N0, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim1, Dim2, Dim0, Dim4, l, n, j, k, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N1, N2, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim1, Dim2, Dim4, Dim0, l, n, j, k, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N1, N2, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim1, Dim4, Dim0, Dim2, l, n, j, m, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N1, N4, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim1, Dim4, Dim2, Dim0, l, n, j, m, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N1, N4, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim2, Dim0, Dim1, Dim4, l, n, k, i, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N2, N0, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim2, Dim0, Dim4, Dim1, l, n, k, i, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N2, N0, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim2, Dim1, Dim0, Dim4, l, n, k, j, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N2, N1, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim2, Dim1, Dim4, Dim0, l, n, k, j, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N2, N1, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim2, Dim4, Dim0, Dim1, l, n, k, m, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N2, N4, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim2, Dim4, Dim1, Dim0, l, n, k, m, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N2, N4, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim4, Dim0, Dim1, Dim2, l, n, m, i, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N4, N0, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim4, Dim0, Dim2, Dim1, l, n, m, i, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N4, N0, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim4, Dim1, Dim0, Dim2, l, n, m, j, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N4, N1, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim4, Dim1, Dim2, Dim0, l, n, m, j, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N4, N1, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim4, Dim2, Dim0, Dim1, l, n, m, k, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N4, N2, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim3, Dim5, Dim4, Dim2, Dim1, Dim0, l, n, m, k, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N3, N5, N4, N2, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim1, Dim2, Dim3, Dim5, m, i, j, k, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N1, N2, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim1, Dim2, Dim5, Dim3, m, i, j, k, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N1, N2, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim1, Dim3, Dim2, Dim5, m, i, j, l, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N1, N3, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim1, Dim3, Dim5, Dim2, m, i, j, l, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N1, N3, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim1, Dim5, Dim2, Dim3, m, i, j, n, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N1, N5, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim1, Dim5, Dim3, Dim2, m, i, j, n, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N1, N5, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim2, Dim1, Dim3, Dim5, m, i, k, j, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N2, N1, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim2, Dim1, Dim5, Dim3, m, i, k, j, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N2, N1, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim2, Dim3, Dim1, Dim5, m, i, k, l, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N2, N3, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim2, Dim3, Dim5, Dim1, m, i, k, l, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N2, N3, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim2, Dim5, Dim1, Dim3, m, i, k, n, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N2, N5, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim2, Dim5, Dim3, Dim1, m, i, k, n, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N2, N5, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim3, Dim1, Dim2, Dim5, m, i, l, j, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N3, N1, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim3, Dim1, Dim5, Dim2, m, i, l, j, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N3, N1, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim3, Dim2, Dim1, Dim5, m, i, l, k, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N3, N2, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim3, Dim2, Dim5, Dim1, m, i, l, k, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N3, N2, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim3, Dim5, Dim1, Dim2, m, i, l, n, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N3, N5, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim3, Dim5, Dim2, Dim1, m, i, l, n, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N3, N5, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim5, Dim1, Dim2, Dim3, m, i, n, j, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N5, N1, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim5, Dim1, Dim3, Dim2, m, i, n, j, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N5, N1, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim5, Dim2, Dim1, Dim3, m, i, n, k, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N5, N2, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim5, Dim2, Dim3, Dim1, m, i, n, k, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N5, N2, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim5, Dim3, Dim1, Dim2, m, i, n, l, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N5, N3, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim0, Dim5, Dim3, Dim2, Dim1, m, i, n, l, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N0, N5, N3, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim0, Dim2, Dim3, Dim5, m, j, i, k, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N0, N2, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim0, Dim2, Dim5, Dim3, m, j, i, k, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N0, N2, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim0, Dim3, Dim2, Dim5, m, j, i, l, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N0, N3, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim0, Dim3, Dim5, Dim2, m, j, i, l, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N0, N3, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim0, Dim5, Dim2, Dim3, m, j, i, n, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N0, N5, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim0, Dim5, Dim3, Dim2, m, j, i, n, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N0, N5, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim2, Dim0, Dim3, Dim5, m, j, k, i, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N2, N0, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim2, Dim0, Dim5, Dim3, m, j, k, i, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N2, N0, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim2, Dim3, Dim0, Dim5, m, j, k, l, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N2, N3, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim2, Dim3, Dim5, Dim0, m, j, k, l, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N2, N3, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim2, Dim5, Dim0, Dim3, m, j, k, n, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N2, N5, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim2, Dim5, Dim3, Dim0, m, j, k, n, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N2, N5, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim3, Dim0, Dim2, Dim5, m, j, l, i, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N3, N0, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim3, Dim0, Dim5, Dim2, m, j, l, i, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N3, N0, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim3, Dim2, Dim0, Dim5, m, j, l, k, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N3, N2, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim3, Dim2, Dim5, Dim0, m, j, l, k, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N3, N2, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim3, Dim5, Dim0, Dim2, m, j, l, n, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N3, N5, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim3, Dim5, Dim2, Dim0, m, j, l, n, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N3, N5, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim5, Dim0, Dim2, Dim3, m, j, n, i, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N5, N0, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim5, Dim0, Dim3, Dim2, m, j, n, i, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N5, N0, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim5, Dim2, Dim0, Dim3, m, j, n, k, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N5, N2, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim5, Dim2, Dim3, Dim0, m, j, n, k, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N5, N2, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim5, Dim3, Dim0, Dim2, m, j, n, l, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N5, N3, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim1, Dim5, Dim3, Dim2, Dim0, m, j, n, l, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N1, N5, N3, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim0, Dim1, Dim3, Dim5, m, k, i, j, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N0, N1, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim0, Dim1, Dim5, Dim3, m, k, i, j, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N0, N1, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim0, Dim3, Dim1, Dim5, m, k, i, l, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N0, N3, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim0, Dim3, Dim5, Dim1, m, k, i, l, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N0, N3, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim0, Dim5, Dim1, Dim3, m, k, i, n, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N0, N5, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim0, Dim5, Dim3, Dim1, m, k, i, n, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N0, N5, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim1, Dim0, Dim3, Dim5, m, k, j, i, l, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N1, N0, N3, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim1, Dim0, Dim5, Dim3, m, k, j, i, n, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N1, N0, N5, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim1, Dim3, Dim0, Dim5, m, k, j, l, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N1, N3, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim1, Dim3, Dim5, Dim0, m, k, j, l, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N1, N3, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim1, Dim5, Dim0, Dim3, m, k, j, n, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N1, N5, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim1, Dim5, Dim3, Dim0, m, k, j, n, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N1, N5, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim3, Dim0, Dim1, Dim5, m, k, l, i, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N3, N0, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim3, Dim0, Dim5, Dim1, m, k, l, i, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N3, N0, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim3, Dim1, Dim0, Dim5, m, k, l, j, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N3, N1, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim3, Dim1, Dim5, Dim0, m, k, l, j, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N3, N1, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim3, Dim5, Dim0, Dim1, m, k, l, n, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N3, N5, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim3, Dim5, Dim1, Dim0, m, k, l, n, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N3, N5, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim5, Dim0, Dim1, Dim3, m, k, n, i, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N5, N0, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim5, Dim0, Dim3, Dim1, m, k, n, i, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N5, N0, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim5, Dim1, Dim0, Dim3, m, k, n, j, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N5, N1, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim5, Dim1, Dim3, Dim0, m, k, n, j, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N5, N1, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim5, Dim3, Dim0, Dim1, m, k, n, l, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N5, N3, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim2, Dim5, Dim3, Dim1, Dim0, m, k, n, l, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N2, N5, N3, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim0, Dim1, Dim2, Dim5, m, l, i, j, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N0, N1, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim0, Dim1, Dim5, Dim2, m, l, i, j, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N0, N1, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim0, Dim2, Dim1, Dim5, m, l, i, k, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N0, N2, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim0, Dim2, Dim5, Dim1, m, l, i, k, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N0, N2, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim0, Dim5, Dim1, Dim2, m, l, i, n, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N0, N5, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim0, Dim5, Dim2, Dim1, m, l, i, n, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N0, N5, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim1, Dim0, Dim2, Dim5, m, l, j, i, k, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N1, N0, N2, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim1, Dim0, Dim5, Dim2, m, l, j, i, n, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N1, N0, N5, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim1, Dim2, Dim0, Dim5, m, l, j, k, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N1, N2, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim1, Dim2, Dim5, Dim0, m, l, j, k, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N1, N2, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim1, Dim5, Dim0, Dim2, m, l, j, n, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N1, N5, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim1, Dim5, Dim2, Dim0, m, l, j, n, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N1, N5, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim2, Dim0, Dim1, Dim5, m, l, k, i, j, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N2, N0, N1, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim2, Dim0, Dim5, Dim1, m, l, k, i, n, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N2, N0, N5, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim2, Dim1, Dim0, Dim5, m, l, k, j, i, n> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N2, N1, N0, N5 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim2, Dim1, Dim5, Dim0, m, l, k, j, n, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N2, N1, N5, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim2, Dim5, Dim0, Dim1, m, l, k, n, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N2, N5, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim2, Dim5, Dim1, Dim0, m, l, k, n, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N2, N5, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim5, Dim0, Dim1, Dim2, m, l, n, i, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N5, N0, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim5, Dim0, Dim2, Dim1, m, l, n, i, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N5, N0, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim5, Dim1, Dim0, Dim2, m, l, n, j, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N5, N1, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim5, Dim1, Dim2, Dim0, m, l, n, j, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N5, N1, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim5, Dim2, Dim0, Dim1, m, l, n, k, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N5, N2, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim3, Dim5, Dim2, Dim1, Dim0, m, l, n, k, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N3, N5, N2, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N0, N1, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim0, Dim1, Dim3, Dim2, m, n, i, j, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N0, N1, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim0, Dim2, Dim1, Dim3, m, n, i, k, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N0, N2, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim0, Dim2, Dim3, Dim1, m, n, i, k, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N0, N2, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim0, Dim3, Dim1, Dim2, m, n, i, l, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N0, N3, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim0, Dim3, Dim2, Dim1, m, n, i, l, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N0, N3, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim1, Dim0, Dim2, Dim3, m, n, j, i, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N1, N0, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim1, Dim0, Dim3, Dim2, m, n, j, i, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N1, N0, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim1, Dim2, Dim0, Dim3, m, n, j, k, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N1, N2, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim1, Dim2, Dim3, Dim0, m, n, j, k, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N1, N2, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim1, Dim3, Dim0, Dim2, m, n, j, l, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N1, N3, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim1, Dim3, Dim2, Dim0, m, n, j, l, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N1, N3, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim2, Dim0, Dim1, Dim3, m, n, k, i, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N2, N0, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim2, Dim0, Dim3, Dim1, m, n, k, i, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N2, N0, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim2, Dim1, Dim0, Dim3, m, n, k, j, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N2, N1, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim2, Dim1, Dim3, Dim0, m, n, k, j, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N2, N1, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim2, Dim3, Dim0, Dim1, m, n, k, l, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N2, N3, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim2, Dim3, Dim1, Dim0, m, n, k, l, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N2, N3, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim3, Dim0, Dim1, Dim2, m, n, l, i, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N3, N0, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim3, Dim0, Dim2, Dim1, m, n, l, i, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N3, N0, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim3, Dim1, Dim0, Dim2, m, n, l, j, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N3, N1, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim3, Dim1, Dim2, Dim0, m, n, l, j, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N3, N1, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim3, Dim2, Dim0, Dim1, m, n, l, k, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N3, N2, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim4, Dim5, Dim3, Dim2, Dim1, Dim0, m, n, l, k, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N4, N5, N3, N2, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim1, Dim2, Dim3, Dim4, n, i, j, k, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N1, N2, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim1, Dim2, Dim4, Dim3, n, i, j, k, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N1, N2, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim1, Dim3, Dim2, Dim4, n, i, j, l, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N1, N3, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim1, Dim3, Dim4, Dim2, n, i, j, l, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N1, N3, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim1, Dim4, Dim2, Dim3, n, i, j, m, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N1, N4, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim1, Dim4, Dim3, Dim2, n, i, j, m, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N1, N4, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim2, Dim1, Dim3, Dim4, n, i, k, j, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N2, N1, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim2, Dim1, Dim4, Dim3, n, i, k, j, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N2, N1, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim2, Dim3, Dim1, Dim4, n, i, k, l, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N2, N3, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim2, Dim3, Dim4, Dim1, n, i, k, l, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N2, N3, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim2, Dim4, Dim1, Dim3, n, i, k, m, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N2, N4, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim2, Dim4, Dim3, Dim1, n, i, k, m, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N2, N4, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim3, Dim1, Dim2, Dim4, n, i, l, j, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N3, N1, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim3, Dim1, Dim4, Dim2, n, i, l, j, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N3, N1, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim3, Dim2, Dim1, Dim4, n, i, l, k, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N3, N2, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim3, Dim2, Dim4, Dim1, n, i, l, k, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N3, N2, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim3, Dim4, Dim1, Dim2, n, i, l, m, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N3, N4, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim3, Dim4, Dim2, Dim1, n, i, l, m, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N3, N4, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim4, Dim1, Dim2, Dim3, n, i, m, j, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N4, N1, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim4, Dim1, Dim3, Dim2, n, i, m, j, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N4, N1, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim4, Dim2, Dim1, Dim3, n, i, m, k, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N4, N2, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim4, Dim2, Dim3, Dim1, n, i, m, k, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N4, N2, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim4, Dim3, Dim1, Dim2, n, i, m, l, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N4, N3, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim0, Dim4, Dim3, Dim2, Dim1, n, i, m, l, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N0, N4, N3, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim0, Dim2, Dim3, Dim4, n, j, i, k, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N0, N2, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim0, Dim2, Dim4, Dim3, n, j, i, k, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N0, N2, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim0, Dim3, Dim2, Dim4, n, j, i, l, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N0, N3, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim0, Dim3, Dim4, Dim2, n, j, i, l, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N0, N3, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim0, Dim4, Dim2, Dim3, n, j, i, m, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N0, N4, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim0, Dim4, Dim3, Dim2, n, j, i, m, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N0, N4, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim2, Dim0, Dim3, Dim4, n, j, k, i, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N2, N0, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim2, Dim0, Dim4, Dim3, n, j, k, i, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N2, N0, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim2, Dim3, Dim0, Dim4, n, j, k, l, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N2, N3, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim2, Dim3, Dim4, Dim0, n, j, k, l, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N2, N3, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim2, Dim4, Dim0, Dim3, n, j, k, m, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N2, N4, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim2, Dim4, Dim3, Dim0, n, j, k, m, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N2, N4, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim3, Dim0, Dim2, Dim4, n, j, l, i, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N3, N0, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim3, Dim0, Dim4, Dim2, n, j, l, i, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N3, N0, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim3, Dim2, Dim0, Dim4, n, j, l, k, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N3, N2, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim3, Dim2, Dim4, Dim0, n, j, l, k, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N3, N2, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim3, Dim4, Dim0, Dim2, n, j, l, m, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N3, N4, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim3, Dim4, Dim2, Dim0, n, j, l, m, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N3, N4, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim4, Dim0, Dim2, Dim3, n, j, m, i, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N4, N0, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim4, Dim0, Dim3, Dim2, n, j, m, i, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N4, N0, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim4, Dim2, Dim0, Dim3, n, j, m, k, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N4, N2, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim4, Dim2, Dim3, Dim0, n, j, m, k, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N4, N2, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim4, Dim3, Dim0, Dim2, n, j, m, l, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N4, N3, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim1, Dim4, Dim3, Dim2, Dim0, n, j, m, l, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N1, N4, N3, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim0, Dim1, Dim3, Dim4, n, k, i, j, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N0, N1, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim0, Dim1, Dim4, Dim3, n, k, i, j, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N0, N1, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim0, Dim3, Dim1, Dim4, n, k, i, l, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N0, N3, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim0, Dim3, Dim4, Dim1, n, k, i, l, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N0, N3, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim0, Dim4, Dim1, Dim3, n, k, i, m, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N0, N4, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim0, Dim4, Dim3, Dim1, n, k, i, m, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N0, N4, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim1, Dim0, Dim3, Dim4, n, k, j, i, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N1, N0, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim1, Dim0, Dim4, Dim3, n, k, j, i, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N1, N0, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim1, Dim3, Dim0, Dim4, n, k, j, l, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N1, N3, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim1, Dim3, Dim4, Dim0, n, k, j, l, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N1, N3, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim1, Dim4, Dim0, Dim3, n, k, j, m, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N1, N4, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim1, Dim4, Dim3, Dim0, n, k, j, m, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N1, N4, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim3, Dim0, Dim1, Dim4, n, k, l, i, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N3, N0, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim3, Dim0, Dim4, Dim1, n, k, l, i, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N3, N0, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim3, Dim1, Dim0, Dim4, n, k, l, j, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N3, N1, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim3, Dim1, Dim4, Dim0, n, k, l, j, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N3, N1, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim3, Dim4, Dim0, Dim1, n, k, l, m, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N3, N4, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim3, Dim4, Dim1, Dim0, n, k, l, m, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N3, N4, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim4, Dim0, Dim1, Dim3, n, k, m, i, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N4, N0, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim4, Dim0, Dim3, Dim1, n, k, m, i, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N4, N0, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim4, Dim1, Dim0, Dim3, n, k, m, j, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N4, N1, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim4, Dim1, Dim3, Dim0, n, k, m, j, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N4, N1, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim4, Dim3, Dim0, Dim1, n, k, m, l, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N4, N3, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim2, Dim4, Dim3, Dim1, Dim0, n, k, m, l, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N2, N4, N3, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim0, Dim1, Dim2, Dim4, n, l, i, j, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N0, N1, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim0, Dim1, Dim4, Dim2, n, l, i, j, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N0, N1, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim0, Dim2, Dim1, Dim4, n, l, i, k, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N0, N2, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim0, Dim2, Dim4, Dim1, n, l, i, k, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N0, N2, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim0, Dim4, Dim1, Dim2, n, l, i, m, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N0, N4, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim0, Dim4, Dim2, Dim1, n, l, i, m, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N0, N4, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim1, Dim0, Dim2, Dim4, n, l, j, i, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N1, N0, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim1, Dim0, Dim4, Dim2, n, l, j, i, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N1, N0, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim1, Dim2, Dim0, Dim4, n, l, j, k, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N1, N2, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim1, Dim2, Dim4, Dim0, n, l, j, k, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N1, N2, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim1, Dim4, Dim0, Dim2, n, l, j, m, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N1, N4, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim1, Dim4, Dim2, Dim0, n, l, j, m, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N1, N4, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim2, Dim0, Dim1, Dim4, n, l, k, i, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N2, N0, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim2, Dim0, Dim4, Dim1, n, l, k, i, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N2, N0, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim2, Dim1, Dim0, Dim4, n, l, k, j, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N2, N1, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim2, Dim1, Dim4, Dim0, n, l, k, j, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N2, N1, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim2, Dim4, Dim0, Dim1, n, l, k, m, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N2, N4, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim2, Dim4, Dim1, Dim0, n, l, k, m, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N2, N4, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim4, Dim0, Dim1, Dim2, n, l, m, i, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N4, N0, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim4, Dim0, Dim2, Dim1, n, l, m, i, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N4, N0, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim4, Dim1, Dim0, Dim2, n, l, m, j, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N4, N1, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim4, Dim1, Dim2, Dim0, n, l, m, j, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N4, N1, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim4, Dim2, Dim0, Dim1, n, l, m, k, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N4, N2, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim3, Dim4, Dim2, Dim1, Dim0, n, l, m, k, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N3, N4, N2, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim0, Dim1, Dim2, Dim3, n, m, i, j, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N0, N1, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim0, Dim1, Dim3, Dim2, n, m, i, j, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N0, N1, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim0, Dim2, Dim1, Dim3, n, m, i, k, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N0, N2, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim0, Dim2, Dim3, Dim1, n, m, i, k, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N0, N2, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim0, Dim3, Dim1, Dim2, n, m, i, l, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N0, N3, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim0, Dim3, Dim2, Dim1, n, m, i, l, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N0, N3, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim1, Dim0, Dim2, Dim3, n, m, j, i, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N1, N0, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim1, Dim0, Dim3, Dim2, n, m, j, i, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N1, N0, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim1, Dim2, Dim0, Dim3, n, m, j, k, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N1, N2, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim1, Dim2, Dim3, Dim0, n, m, j, k, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N1, N2, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim1, Dim3, Dim0, Dim2, n, m, j, l, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N1, N3, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim1, Dim3, Dim2, Dim0, n, m, j, l, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N1, N3, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim2, Dim0, Dim1, Dim3, n, m, k, i, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N2, N0, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim2, Dim0, Dim3, Dim1, n, m, k, i, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N2, N0, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim2, Dim1, Dim0, Dim3, n, m, k, j, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N2, N1, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim2, Dim1, Dim3, Dim0, n, m, k, j, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N2, N1, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim2, Dim3, Dim0, Dim1, n, m, k, l, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N2, N3, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim2, Dim3, Dim1, Dim0, n, m, k, l, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N2, N3, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim3, Dim0, Dim1, Dim2, n, m, l, i, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N3, N0, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim3, Dim0, Dim2, Dim1, n, m, l, i, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N3, N0, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim3, Dim1, Dim0, Dim2, n, m, l, j, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N3, N1, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim3, Dim1, Dim2, Dim0, n, m, l, j, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N3, N1, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim3, Dim2, Dim0, Dim1, n, m, l, k, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N3, N2, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor6_Expr<B, U, Dim5, Dim4, Dim3, Dim2, Dim1, Dim0, n, m, l, k, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    return rhs( N5, N4, N3, N2, N1, N0 );
  }

// Catch all version for incompatible indices
  template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3, int Dim1_4, int Dim1_5, char i1, char j1, char k1, char l1, char m1, char n1>
  U eval( const Tensor6_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, Dim1_4, Dim1_5, i1, j1, k1, l1, m1, n1> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    static_assert( error_when_instantiated<B>(), "Incompatible indices" );
    return rhs( 0, 0, 0, 0, 0, 0 );
  }
};

}

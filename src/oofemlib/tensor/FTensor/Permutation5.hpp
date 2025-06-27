#pragma once

#include "error_when_instantiated.hpp"

namespace FTensor
{
template <int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
struct Permutation5 {
  Permutation5() = default;

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N1, N2, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim1, Dim2, Dim4, Dim3, i, j, k, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N1, N2, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim1, Dim3, Dim2, Dim4, i, j, l, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N1, N3, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim1, Dim3, Dim4, Dim2, i, j, l, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N1, N3, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim1, Dim4, Dim2, Dim3, i, j, m, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N1, N4, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim1, Dim4, Dim3, Dim2, i, j, m, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N1, N4, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim2, Dim1, Dim3, Dim4, i, k, j, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N2, N1, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim2, Dim1, Dim4, Dim3, i, k, j, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N2, N1, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim2, Dim3, Dim1, Dim4, i, k, l, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N2, N3, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim2, Dim3, Dim4, Dim1, i, k, l, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N2, N3, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim2, Dim4, Dim1, Dim3, i, k, m, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N2, N4, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim2, Dim4, Dim3, Dim1, i, k, m, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N2, N4, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim3, Dim1, Dim2, Dim4, i, l, j, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N3, N1, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim3, Dim1, Dim4, Dim2, i, l, j, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N3, N1, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim3, Dim2, Dim1, Dim4, i, l, k, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N3, N2, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim3, Dim2, Dim4, Dim1, i, l, k, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N3, N2, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim3, Dim4, Dim1, Dim2, i, l, m, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N3, N4, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim3, Dim4, Dim2, Dim1, i, l, m, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N3, N4, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim4, Dim1, Dim2, Dim3, i, m, j, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N4, N1, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim4, Dim1, Dim3, Dim2, i, m, j, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N4, N1, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim4, Dim2, Dim1, Dim3, i, m, k, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N4, N2, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim4, Dim2, Dim3, Dim1, i, m, k, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N4, N2, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim4, Dim3, Dim1, Dim2, i, m, l, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N4, N3, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim0, Dim4, Dim3, Dim2, Dim1, i, m, l, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N0, N4, N3, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N0, N2, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim0, Dim2, Dim4, Dim3, j, i, k, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N0, N2, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim0, Dim3, Dim2, Dim4, j, i, l, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N0, N3, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim0, Dim3, Dim4, Dim2, j, i, l, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N0, N3, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim0, Dim4, Dim2, Dim3, j, i, m, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N0, N4, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim0, Dim4, Dim3, Dim2, j, i, m, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N0, N4, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim2, Dim0, Dim3, Dim4, j, k, i, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N2, N0, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim2, Dim0, Dim4, Dim3, j, k, i, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N2, N0, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim2, Dim3, Dim0, Dim4, j, k, l, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N2, N3, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim2, Dim3, Dim4, Dim0, j, k, l, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N2, N3, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim2, Dim4, Dim0, Dim3, j, k, m, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N2, N4, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim2, Dim4, Dim3, Dim0, j, k, m, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N2, N4, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim3, Dim0, Dim2, Dim4, j, l, i, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N3, N0, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim3, Dim0, Dim4, Dim2, j, l, i, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N3, N0, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim3, Dim2, Dim0, Dim4, j, l, k, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N3, N2, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim3, Dim2, Dim4, Dim0, j, l, k, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N3, N2, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim3, Dim4, Dim0, Dim2, j, l, m, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N3, N4, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim3, Dim4, Dim2, Dim0, j, l, m, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N3, N4, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim4, Dim0, Dim2, Dim3, j, m, i, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N4, N0, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim4, Dim0, Dim3, Dim2, j, m, i, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N4, N0, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim4, Dim2, Dim0, Dim3, j, m, k, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N4, N2, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim4, Dim2, Dim3, Dim0, j, m, k, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N4, N2, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim4, Dim3, Dim0, Dim2, j, m, l, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N4, N3, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim1, Dim4, Dim3, Dim2, Dim0, j, m, l, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N1, N4, N3, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N0, N1, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim0, Dim1, Dim4, Dim3, k, i, j, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N0, N1, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim0, Dim3, Dim1, Dim4, k, i, l, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N0, N3, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim0, Dim3, Dim4, Dim1, k, i, l, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N0, N3, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim0, Dim4, Dim1, Dim3, k, i, m, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N0, N4, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim0, Dim4, Dim3, Dim1, k, i, m, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N0, N4, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim1, Dim0, Dim3, Dim4, k, j, i, l, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N1, N0, N3, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim1, Dim0, Dim4, Dim3, k, j, i, m, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N1, N0, N4, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim1, Dim3, Dim0, Dim4, k, j, l, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N1, N3, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim1, Dim3, Dim4, Dim0, k, j, l, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N1, N3, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim1, Dim4, Dim0, Dim3, k, j, m, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N1, N4, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim1, Dim4, Dim3, Dim0, k, j, m, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N1, N4, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim3, Dim0, Dim1, Dim4, k, l, i, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N3, N0, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim3, Dim0, Dim4, Dim1, k, l, i, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N3, N0, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim3, Dim1, Dim0, Dim4, k, l, j, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N3, N1, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim3, Dim1, Dim4, Dim0, k, l, j, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N3, N1, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim3, Dim4, Dim0, Dim1, k, l, m, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N3, N4, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim3, Dim4, Dim1, Dim0, k, l, m, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N3, N4, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim4, Dim0, Dim1, Dim3, k, m, i, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N4, N0, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim4, Dim0, Dim3, Dim1, k, m, i, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N4, N0, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim4, Dim1, Dim0, Dim3, k, m, j, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N4, N1, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim4, Dim1, Dim3, Dim0, k, m, j, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N4, N1, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim4, Dim3, Dim0, Dim1, k, m, l, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N4, N3, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim2, Dim4, Dim3, Dim1, Dim0, k, m, l, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N2, N4, N3, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N0, N1, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim0, Dim1, Dim4, Dim2, l, i, j, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N0, N1, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim0, Dim2, Dim1, Dim4, l, i, k, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N0, N2, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim0, Dim2, Dim4, Dim1, l, i, k, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N0, N2, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim0, Dim4, Dim1, Dim2, l, i, m, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N0, N4, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim0, Dim4, Dim2, Dim1, l, i, m, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N0, N4, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim1, Dim0, Dim2, Dim4, l, j, i, k, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N1, N0, N2, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim1, Dim0, Dim4, Dim2, l, j, i, m, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N1, N0, N4, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim1, Dim2, Dim0, Dim4, l, j, k, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N1, N2, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim1, Dim2, Dim4, Dim0, l, j, k, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N1, N2, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim1, Dim4, Dim0, Dim2, l, j, m, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N1, N4, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim1, Dim4, Dim2, Dim0, l, j, m, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N1, N4, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim2, Dim0, Dim1, Dim4, l, k, i, j, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N2, N0, N1, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim2, Dim0, Dim4, Dim1, l, k, i, m, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N2, N0, N4, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim2, Dim1, Dim0, Dim4, l, k, j, i, m> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N2, N1, N0, N4 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim2, Dim1, Dim4, Dim0, l, k, j, m, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N2, N1, N4, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim2, Dim4, Dim0, Dim1, l, k, m, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N2, N4, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim2, Dim4, Dim1, Dim0, l, k, m, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N2, N4, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim4, Dim0, Dim1, Dim2, l, m, i, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N4, N0, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim4, Dim0, Dim2, Dim1, l, m, i, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N4, N0, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim4, Dim1, Dim0, Dim2, l, m, j, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N4, N1, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim4, Dim1, Dim2, Dim0, l, m, j, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N4, N1, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim4, Dim2, Dim0, Dim1, l, m, k, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N4, N2, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim3, Dim4, Dim2, Dim1, Dim0, l, m, k, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N3, N4, N2, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N0, N1, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim0, Dim1, Dim3, Dim2, m, i, j, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N0, N1, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim0, Dim2, Dim1, Dim3, m, i, k, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N0, N2, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim0, Dim2, Dim3, Dim1, m, i, k, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N0, N2, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim0, Dim3, Dim1, Dim2, m, i, l, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N0, N3, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim0, Dim3, Dim2, Dim1, m, i, l, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N0, N3, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim1, Dim0, Dim2, Dim3, m, j, i, k, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N1, N0, N2, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim1, Dim0, Dim3, Dim2, m, j, i, l, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N1, N0, N3, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim1, Dim2, Dim0, Dim3, m, j, k, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N1, N2, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim1, Dim2, Dim3, Dim0, m, j, k, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N1, N2, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim1, Dim3, Dim0, Dim2, m, j, l, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N1, N3, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim1, Dim3, Dim2, Dim0, m, j, l, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N1, N3, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim2, Dim0, Dim1, Dim3, m, k, i, j, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N2, N0, N1, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim2, Dim0, Dim3, Dim1, m, k, i, l, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N2, N0, N3, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim2, Dim1, Dim0, Dim3, m, k, j, i, l> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N2, N1, N0, N3 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim2, Dim1, Dim3, Dim0, m, k, j, l, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N2, N1, N3, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim2, Dim3, Dim0, Dim1, m, k, l, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N2, N3, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim2, Dim3, Dim1, Dim0, m, k, l, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N2, N3, N1, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim3, Dim0, Dim1, Dim2, m, l, i, j, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N3, N0, N1, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim3, Dim0, Dim2, Dim1, m, l, i, k, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N3, N0, N2, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim3, Dim1, Dim0, Dim2, m, l, j, i, k> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N3, N1, N0, N2 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim3, Dim1, Dim2, Dim0, m, l, j, k, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N3, N1, N2, N0 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim3, Dim2, Dim0, Dim1, m, l, k, i, j> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N3, N2, N0, N1 );
  }

  template <class B, class U>
  U eval( const Tensor5_Expr<B, U, Dim4, Dim3, Dim2, Dim1, Dim0, m, l, k, j, i> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    return rhs( N4, N3, N2, N1, N0 );
  }

    // Catch all version for incompatible indices
  template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3, int Dim1_4, char i1, char j1, char k1, char l1, char m1>
  U eval( const Tensor5_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, Dim1_4, i1, j1, k1, l1, m1> &rhs,
      const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    static_assert( error_when_instantiated<B>(), "Incompatible indices" );
    return rhs( 0, 0, 0, 0, 0 );
  }
};

}

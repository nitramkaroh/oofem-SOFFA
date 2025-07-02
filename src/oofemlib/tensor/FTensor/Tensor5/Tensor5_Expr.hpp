/* Declare a wrapper class for generic rank 5 Tensor expressions. */

#pragma once

#include "Tensor5_minus_Tensor5.hpp"
#include "Tensor5_plus_Tensor5.hpp"
#include "Tensor5_times_Tensor1_single.hpp"
#include "Tensor5_times_Tensor2_single.hpp"
#include "Tensor5_times_Tensor2_double.hpp"
//#include "Tensor5_times_Tensor2_symmetric.hpp"
#include "Tensor5_times_Tensor3_triple.hpp"
#include "Tensor5_times_generic.hpp"

#include "../permute.hpp"

namespace FTensor
{
  template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i,
            char j, char k, char l, char m>
  class Tensor5_Expr
  {
    A iter;

  public:
    Tensor5_Expr(const A &a) : iter(a) {}
    T operator()(const int N1, const int N2, const int N3, const int N4, const int N5) const
    {
      return iter(N1, N2, N3, N4, N5);
    }
  };

  template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i,
            char j, char k, char l, char m>
  class Tensor5_Expr<Tensor5<A, Dim0, Dim1, Dim2, Dim3, Dim4>, T, Dim0, Dim1, Dim2,
                     Dim3, Dim4, i, j, k, l, m>
  {
    Tensor5<A, Dim0, Dim1, Dim2, Dim3, Dim4> &iter;

  public:
    Tensor5_Expr(Tensor5<A, Dim0, Dim1, Dim2, Dim3, Dim4> &a) : iter(a) {}
    T &operator()(const int N1, const int N2, const int N3, const int N4, const int N5)
    {
      return iter(N1, N2, N3, N4, N5);
    }
    T operator()(const int N1, const int N2, const int N3, const int N4, const int N5) const
    {
      return iter(N1, N2, N3, N4, N5);
    }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3, int Dim1_4,
              char i_1, char j_1, char k_1, char l_1, char m_1>
    auto &equals(const Tensor5_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, Dim1_4, i_1,
                                    j_1, k_1, l_1, m_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          for(int kk = 0; kk < Dim2; ++kk)
            for(int ll = 0; ll < Dim3; ++ll)
              for ( int mm = 0; mm < Dim4; ++mm )
              {
                iter(ii, jj, kk, ll, mm) = permute(*this, rhs, ii, jj, kk, ll, mm);
              }
      return *this;
    }

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3, int Dim1_4,
              char i_1, char j_1, char k_1, char l_1, char m_1>
    auto &operator=(const Tensor5_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, Dim1_4,
                                       i_1, j_1, k_1, l_1, m_1> &rhs)
    {
      return equals(rhs);
    }

    auto &
    operator=(const Tensor5_Expr<Tensor5<A, Dim0, Dim1, Dim2, Dim3, Dim4>, T, Dim0,
                                 Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &rhs)
    {
      return equals(rhs);
    }

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3, int Dim1_4,
              char i_1, char j_1, char k_1, char l_1, char m_1>
    auto &operator+=(const Tensor5_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, Dim1_4,
                                        i_1, j_1, k_1, l_1, m_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          for(int kk = 0; kk < Dim2; ++kk)
            for(int ll = 0; ll < Dim3; ++ll)
              for ( int mm = 0; mm < Dim4; ++mm )
              {
                iter(ii, jj, kk, ll, mm) += permute(*this, rhs, ii, jj, kk, ll, mm);
              }
      return *this;
    }

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3, int Dim1_4,
              char i_1, char j_1, char k_1, char l_1, char m_1>
    auto &operator-=(const Tensor5_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, Dim1_4,
                                        i_1, j_1, k_1, l_1, m_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          for(int kk = 0; kk < Dim2; ++kk)
            for(int ll = 0; ll < Dim3; ++ll)
              for ( int mm = 0; mm < Dim4; ++mm )
              {
                iter(ii, jj, kk, ll, mm) += permute(*this, rhs, ii, jj, kk, ll, mm);
              }
      return *this;
    }

    template <class U>
    Tensor5_Expr<Tensor5<A, Dim0, Dim1, Dim2, Dim3, Dim4>, T, Dim0,
        Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &
    operator=( const U &u )
    {
      for ( int ii = 0; ii < Dim0; ++ii )
        for ( int jj = 0; jj < Dim1; ++jj )
          for ( int kk = 0; kk < Dim2; ++kk )
            for ( int ll = 0; ll < Dim3; ++ll )
              for ( int mm = 0; mm < Dim4; ++mm ) {
                iter( ii, jj, kk, ll, mm ) = u;
              }
      return *this;
    }

  };
}

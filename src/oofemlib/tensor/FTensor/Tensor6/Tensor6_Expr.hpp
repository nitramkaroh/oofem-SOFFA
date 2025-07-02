/* Declare a wrapper class for generic rank 6 Tensor expressions. */

#pragma once

#include "Tensor6_minus_Tensor6.hpp"
#include "Tensor6_plus_Tensor6.hpp"
#include "Tensor6_times_Tensor1_single.hpp"
#include "Tensor6_times_Tensor2_single.hpp"
#include "Tensor6_times_Tensor2_double.hpp"
#include "Tensor6_times_Tensor2_symmetric.hpp"
#include "Tensor6_times_Tensor3_triple.hpp"
#include "Tensor6_times_generic.hpp"

#include "../permute.hpp"

namespace FTensor
{
  template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i,
            char j, char k, char l, char m, char n>
  class Tensor6_Expr
  {
    A iter;

  public:
    Tensor6_Expr(const A &a) : iter(a) {}
    T operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
    {
      return iter(N1, N2, N3, N4, N5, N6);
    }
  };

  template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i,
            char j, char k, char l, char m, char n>
  class Tensor6_Expr<Tensor6<A, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5>, T, Dim0, Dim1, Dim2,
                     Dim3, Dim4, Dim5, i, j, k, l, m, n>
  {
    Tensor6<A, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5> &iter;

  public:
    Tensor6_Expr(Tensor6<A, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5> &a) : iter(a) {}
    T &operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6)
    {
      return iter(N1, N2, N3, N4, N5. N6);
    }
    T operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
    {
      return iter(N1, N2, N3, N4, N5, N6);
    }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3, int Dim1_4, int Dim1_5,
              char i_1, char j_1, char k_1, char l_1, char m_1, char n_1>
    auto &equals(const Tensor6_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, Dim1_4, Dim1_5, i_1,
                                    j_1, k_1, l_1, m_1, n_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          for(int kk = 0; kk < Dim2; ++kk)
            for(int ll = 0; ll < Dim3; ++ll)
              for ( int mm = 0; mm < Dim4; ++mm )
                for ( int nn = 0; nn < Dim5; ++nn )
                {
                  iter( ii, jj, kk, ll, mm, nn ) = permute( *this, rhs, ii, jj, kk, ll, mm, nn );
                }
      return *this;
    }

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3, int Dim1_4, int Dim1_5,
              char i_1, char j_1, char k_1, char l_1, char m_1, char n_1>
    auto &operator=(const Tensor6_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, Dim1_4, Dim1_5,
                                       i_1, j_1, k_1, l_1, m_1, n_1> &rhs)
    {
      return equals(rhs);
    }

    auto &
    operator=(const Tensor6_Expr<Tensor6<A, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5>, T, Dim0,
                                 Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &rhs)
    {
      return equals(rhs);
    }

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3, int Dim1_4, int Dim1_5,
              char i_1, char j_1, char k_1, char l_1, char m_1, char n_1>
    auto &operator+=(const Tensor6_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, Dim1_4, Dim1_5,
                                        i_1, j_1, k_1, l_1, m_1, n_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          for(int kk = 0; kk < Dim2; ++kk)
            for(int ll = 0; ll < Dim3; ++ll)
              for ( int mm = 0; mm < Dim4; ++mm )
                for ( int nn = 0; nn < Dim5; ++nn ) 
                {
                  iter( ii, jj, kk, ll, mm, nn ) += permute( *this, rhs, ii, jj, kk, ll, mm, nn );
                }
      return *this;
    }

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3, int Dim1_4, int Dim1_5,
              char i_1, char j_1, char k_1, char l_1, char m_1, char n_1>
    auto &operator-=(const Tensor6_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, Dim1_4, Dim1_5,
                                        i_1, j_1, k_1, l_1, m_1, n_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          for(int kk = 0; kk < Dim2; ++kk)
            for(int ll = 0; ll < Dim3; ++ll)
              for ( int mm = 0; mm < Dim4; ++mm )
                for ( int nn = 0; nn < Dim5; ++nn )
                {
                  iter( ii, jj, kk, ll, mm, nn ) += permute( *this, rhs, ii, jj, kk, ll, mm, nn );
                }
      return *this;
    }

    template <class U>
    Tensor6_Expr<Tensor6<A, Dim0, Dim1, Dim2, D
      im3, Dim4, Dim5>, T, Dim0,
        Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &
    operator=( const U &u )
    {
      for ( int ii = 0; ii < Dim0; ++ii )
        for ( int jj = 0; jj < Dim1; ++jj )
          for ( int kk = 0; kk < Dim2; ++kk )
            for ( int ll = 0; ll < Dim3; ++ll )
              for ( int mm = 0; mm < Dim4; ++mm )
                for ( int nn = 0; nn < Dim5; ++nn ) {
                  iter( ii, jj, kk, ll, mm, nn ) = u;
                }
      return *this;
    }

  };
}

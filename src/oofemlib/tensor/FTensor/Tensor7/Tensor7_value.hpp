/* A general version, not for pointers. */

#pragma once

#include "Tensor7_contracted.hpp"

#include <ostream>
#include "floatmatrixf.h"
#ifdef FTENSOR_DEBUG
#include <sstream>
#include <stdexcept>
#endif

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
            int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5, int Tensor_Dim6>
  class Tensor7
  {
  protected:
    T data[Tensor_Dim0][Tensor_Dim1][Tensor_Dim2][Tensor_Dim3][Tensor_Dim4][Tensor_Dim5][Tensor_Dim6];

  public:
    /* Initialization operators */
    template <class... U> Tensor7(U... d) : data{d...}
    {
      static_assert(sizeof...(d) == sizeof(data) / sizeof(T),
                    "Incorrect number of Arguments. Constructor should "
                    "initialize the entire Tensor");
    }

    Tensor7() {}
    /* There are two operator(int,int,int,int,int)'s, one for non-consts
       that lets you change the value, and one for consts that
       doesn't. */
    T &operator()( const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const int N7 )
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0
          || N3 >= Tensor_Dim2 || N3 < 0 || N4 >= Tensor_Dim3 || N4 < 0 || N5 >= Tensor_Dim4 || N5 < 0
          || N6 >= Tensor_Dim5 || N6 < 0 || N7 >= Tensor_Dim6|| N7 < 0 )
        {
          std::stringstream s;
          s << "Bad index in Tensor7<T," << Tensor_Dim0 << "," << Tensor_Dim1
            << "," << Tensor_Dim2 << "," << Tensor_Dim3 << "," << Tensor_Dim4 << "," << Tensor_Dim5 << "," << Tensor_Dim6 << ">.operator(" << N1
            << "," << N2 << "," << N3 << "," << N4 "," << N5 << "," << N6 "," << N7 << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[N1][N2][N3][N4][N5][N6][N7];
    }

    T operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const int N7) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0
          || N3 >= Tensor_Dim2 || N3 < 0 || N4 >= Tensor_Dim3 || N4 < 0 || N5 >= Tensor_Dim4 || N5 < 0
          || N6 >= Tensor_Dim5 || N6 < 0 || N7 >= Tensor_Dim6 || N7 < 0 )
        {
          std::stringstream s;
        s << "Bad index in Tensor7<T," << Tensor_Dim0 << "," << Tensor_Dim1
            << "," << Tensor_Dim2 << "," << Tensor_Dim3 << "," << Tensor_Dim4 "," << Tensor_Dim5 "," << Tensor_Dim6 << ">.operator(" << N1
            << "," << N2 << "," << N3 << "," << N4 "," << N5 << "," << N6 << "," << N7 << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[N1][N2][N3][N4][N5][N6][N7];
    }

    /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. */

    template <char i, char j, char k, char l, char m, char n, char o, int Dim0, int Dim1, int Dim2,
              int Dim3, int Dim4, int Dim5, int Dim6>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2
            && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim6 ),
      Tensor7_Expr<Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>,
                   T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>,
               const Index<k, Dim2>, const Index<l, Dim3>, 
               const Index<m, Dim4>, const Index<n, Dim5>,
               const Index<o, Dim6> )
    {
      return Tensor7_Expr<
        Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T,
        Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>(*this);
    };

    template <char i, char j, char k, char l, char m, char n, char o, int Dim0, int Dim1, int Dim2,
              int Dim3, int Dim4, int Dim5, int Dim6>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2
            && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim6 ),
      Tensor7_Expr<
        const Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>,
        T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>,
               const Index<k, Dim2>, const Index<l, Dim3>,
               const Index<m, Dim4>, const Index<n, Dim5>,
               const Index<o, Dim6> ) const
    {
      return Tensor7_Expr<
        const Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>,
        T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>(*this);
    };

    /* These operators are for internal contractions, resulting in a Tensor5.
     * For example something like A(i,i,k,l,m,n,o) */

  // (i, i, k, l, m, n, o, )
    template <char i, char k, char l, char m, char n, char o, int Dim01, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6>
    auto operator()( const Index<i, Dim01>, const Index<i, Dim01>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<n, Dim5>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim01 && Tensor_Dim1 >= Dim01 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_01<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim01>;
      return Tensor5_Expr<TensorExpr, T, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>( TensorExpr( *this ) );
    };

    // (i, j, i, l, m, n, o, )
    template <char i, char j, char l, char m, char n, char o, int Dim02, int Dim1, int Dim3, int Dim4, int Dim5, int Dim6>
    auto operator()( const Index<i, Dim02>, const Index<j, Dim1>, const Index<i, Dim02>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<n, Dim5>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim02 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim02 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_02<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim02>;
      return Tensor5_Expr<TensorExpr, T, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>( TensorExpr( *this ) );
    };

    // (i, j, k, i, m, n, o, )
    template <char i, char j, char k, char m, char n, char o, int Dim03, int Dim1, int Dim2, int Dim4, int Dim5, int Dim6>
    auto operator()( const Index<i, Dim03>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<i, Dim03>, const Index<m, Dim4>, const Index<n, Dim5>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim03 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim03 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_03<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim03>;
      return Tensor5_Expr<TensorExpr, T, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>( TensorExpr( *this ) );
    };

    // (i, j, k, l, i, n, o, )
    template <char i, char j, char k, char l, char n, char o, int Dim04, int Dim1, int Dim2, int Dim3, int Dim5, int Dim6>
    auto operator()( const Index<i, Dim04>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<i, Dim04>, const Index<n, Dim5>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim04 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim04 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_04<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim04>;
      return Tensor5_Expr<TensorExpr, T, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, i, o, )
    template <char i, char j, char k, char l, char m, char o, int Dim05, int Dim1, int Dim2, int Dim3, int Dim4, int Dim6>
    auto operator()( const Index<i, Dim05>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<i, Dim05>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim05 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim05 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_05<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim05>;
      return Tensor5_Expr<TensorExpr, T, Dim1, Dim2, Dim3, Dim4, Dim6, j, k, l, m, o>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, n, i, )
    template <char i, char j, char k, char l, char m, char n, int Dim06, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5>
    auto operator()( const Index<i, Dim06>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<n, Dim5>, const Index<i, Dim06> ) const
    {
      static_assert( Tensor_Dim0 >= Dim06 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim06,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_06<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim06>;
      return Tensor5_Expr<TensorExpr, T, Dim1, Dim2, Dim3, Dim4, Dim5, j, k, l, m, n>( TensorExpr( *this ) );
    };

    // (i, j, j, l, m, n, o, )
    template <char i, char j, char l, char m, char n, char o, int Dim0, int Dim12, int Dim3, int Dim4, int Dim5, int Dim6>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim12>, const Index<j, Dim12>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<n, Dim5>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim12 && Tensor_Dim2 >= Dim12 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_12<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim12>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>( TensorExpr( *this ) );
    };

    // (i, j, k, j, m, n, o, )
    template <char i, char j, char k, char m, char n, char o, int Dim0, int Dim13, int Dim2, int Dim4, int Dim5, int Dim6>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim13>, const Index<k, Dim2>, const Index<j, Dim13>, const Index<m, Dim4>, const Index<n, Dim5>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim13 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim13 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_13<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim13>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>( TensorExpr( *this ) );
    };

    // (i, j, k, l, j, n, o, )
    template <char i, char j, char k, char l, char n, char o, int Dim0, int Dim14, int Dim2, int Dim3, int Dim5, int Dim6>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim14>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<j, Dim14>, const Index<n, Dim5>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim14 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim14 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_14<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim14>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, j, o, )
    template <char i, char j, char k, char l, char m, char o, int Dim0, int Dim15, int Dim2, int Dim3, int Dim4, int Dim6>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim15>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<j, Dim15>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim15 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim15 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_15<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim15>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim2, Dim3, Dim4, Dim6, i, k, l, m, o>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, n, j, )
    template <char i, char j, char k, char l, char m, char n, int Dim0, int Dim16, int Dim2, int Dim3, int Dim4, int Dim5>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim16>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<n, Dim5>, const Index<j, Dim16> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim16 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim16,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_16<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim16>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim2, Dim3, Dim4, Dim5, i, k, l, m, n>( TensorExpr( *this ) );
    };

    // (i, j, k, k, m, n, o, )
    template <char i, char j, char k, char m, char n, char o, int Dim0, int Dim1, int Dim23, int Dim4, int Dim5, int Dim6>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim23>, const Index<k, Dim23>, const Index<m, Dim4>, const Index<n, Dim5>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim23 && Tensor_Dim3 >= Dim23 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_23<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim23>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>( TensorExpr( *this ) );
    };

    // (i, j, k, l, k, n, o, )
    template <char i, char j, char k, char l, char n, char o, int Dim0, int Dim1, int Dim24, int Dim3, int Dim5, int Dim6>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim24>, const Index<l, Dim3>, const Index<k, Dim24>, const Index<n, Dim5>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim24 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim24 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_24<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim24>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, k, o, )
    template <char i, char j, char k, char l, char m, char o, int Dim0, int Dim1, int Dim25, int Dim3, int Dim4, int Dim6>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim25>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<k, Dim25>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim25 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim25 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_25<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim25>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim1, Dim3, Dim4, Dim6, i, j, l, m, o>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, n, k, )
    template <char i, char j, char k, char l, char m, char n, int Dim0, int Dim1, int Dim26, int Dim3, int Dim4, int Dim5>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim26>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<n, Dim5>, const Index<k, Dim26> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim26 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim26,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_26<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim26>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim1, Dim3, Dim4, Dim5, i, j, l, m, n>( TensorExpr( *this ) );
    };

    // (i, j, k, l, l, n, o, )
    template <char i, char j, char k, char l, char n, char o, int Dim0, int Dim1, int Dim2, int Dim34, int Dim5, int Dim6>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim34>, const Index<l, Dim34>, const Index<n, Dim5>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim34 && Tensor_Dim4 >= Dim34 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_34<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim34>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, l, o, )
    template <char i, char j, char k, char l, char m, char o, int Dim0, int Dim1, int Dim2, int Dim35, int Dim4, int Dim6>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim35>, const Index<m, Dim4>, const Index<l, Dim35>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim35 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim35 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_35<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim35>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim1, Dim2, Dim4, Dim6, i, j, k, m, o>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, n, l, )
    template <char i, char j, char k, char l, char m, char n, int Dim0, int Dim1, int Dim2, int Dim36, int Dim4, int Dim5>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim36>, const Index<m, Dim4>, const Index<n, Dim5>, const Index<l, Dim36> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim36 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim36,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_36<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim36>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim1, Dim2, Dim4, Dim5, i, j, k, m, n>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, m, o, )
    template <char i, char j, char k, char l, char m, char o, int Dim0, int Dim1, int Dim2, int Dim3, int Dim45, int Dim6>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim45>, const Index<m, Dim45>, const Index<o, Dim6> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim45 && Tensor_Dim5 >= Dim45 && Tensor_Dim6 >= Dim6,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_45<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim45>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim1, Dim2, Dim3, Dim6, i, j, k, l, o>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, n, m, )
    template <char i, char j, char k, char l, char m, char n, int Dim0, int Dim1, int Dim2, int Dim3, int Dim46, int Dim5>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim46>, const Index<n, Dim5>, const Index<m, Dim46> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim46 && Tensor_Dim5 >= Dim5 && Tensor_Dim6 >= Dim46,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_46<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim46>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim1, Dim2, Dim3, Dim5, i, j, k, l, n>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, n, n, )
    template <char i, char j, char k, char l, char m, char n, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim56>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<n, Dim56>, const Index<n, Dim56> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim56 && Tensor_Dim6 >= Dim56,
          "Incompatible indices" );

      using TensorExpr = Tensor7_contracted_56<
          Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6>, T, Dim56>;
      return Tensor5_Expr<TensorExpr, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>( TensorExpr( *this ) );
    };







    /* This is for expressions where a number is used for one slot, and
       an index for another, yielding a Tensor3_Expr.  The non-const
       versions don't actually create a Tensor3_number_rhs_[0123] object.
       They create a Tensor3_Expr directly, which provides the
       appropriate indexing operators.  The const versions do create a
       Tensor3_number_[0123]. */
    // TODO
  };



}

/// JSON compatible output

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
            int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5, int Tensor_Dim6>
  std::ostream &Tensor7_0000001(
    std::ostream &os,
    const Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6> &t,
    const int &iterator0, const int &iterator1, const int &iterator2, const int &iterator3, const int &iterator4, const int &iterator5 )
  {
    os << '[';
    for(int i = 0; i < Tensor_Dim6 - 1; ++i)
      {
        os << t(iterator0, iterator1, iterator2, iterator3, iterator4, iterator5, i);
        os << ',';
      }
    os << t(iterator0, iterator1, iterator2, iterator3, iterator4, iterator5, Tensor_Dim6 - 1);
    os << ']';

    return os;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
      int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5, int Tensor_Dim6>
  std::ostream &Tensor7_0000010(
      std::ostream &os,
      const Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6> &t,
      const int &iterator0, const int &iterator1, const int &iterator2, const int &iterator3, const int &iterator4 )
  {
    os << '[';
    for ( int i = 0; i < Tensor_Dim5 - 1; ++i ) {
      FTensor::Tensor7_0000001( os, t, iterator0, iterator1, iterator2, iterator3, iterator4, i );
      os << ',';
    }
    FTensor::Tensor7_0000001( os, t, iterator0, iterator1, iterator2, iterator3, iterator4, Tensor_Dim5 - 1 );
    os << ']';

    return os;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
            int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5, int Tensor_Dim6>
  std::ostream &Tensor7_0000100(
    std::ostream &os,
    const Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6> &t,
    const int &iterator0, const int &iterator1, const int &iterator2, const int &iterator3)
  {
    os << '[';
    for(int i = 0; i < Tensor_Dim4 - 1; ++i)
      {
        FTensor::Tensor7_0000010(os, t, iterator0, iterator1, iterator2, iterator3, i);
        os << ',';
      }
    FTensor::Tensor7_0000010(os, t, iterator0, iterator1, iterator2, Tensor_Dim4 - 1);
    os << ']';

    return os;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
            int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5, int Tensor_Dim6>
  std::ostream &Tensor7_0001000(
    std::ostream &os,
    const Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6> &t,
    const int &iterator0, const int &iterator1, const int &iterator2)
  {
    os << '[';
    for(int i = 0; i < Tensor_Dim3 - 1; ++i)
      {
        FTensor::Tensor7_0000100(os, t, iterator0, iterator1, iterator2, i);
        os << ',';
      }
    FTensor::Tensor7_0000100(os, t, iterator0, iterator1, iterator2, Tensor_Dim3 - 1);
    os << ']';

    return os;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
    int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5, int Tensor_Dim6>
std::ostream &Tensor7_0010000(
    std::ostream &os,
    const Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6> &t,
    const int &iterator0, const int &iterator1 )
{
  os << '[';
  for ( int i = 0; i < Tensor_Dim2 - 1; ++i ) {
    FTensor::Tensor7_0001000( os, t, iterator0, iterator1, i );
    os << ',';
  }
  FTensor::Tensor7_0001000( os, t, iterator0, iterator1, Tensor_Dim2 - 1 );
  os << ']';

  return os;
}

 template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
    int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5, int Tensor_Dim6>
std::ostream &Tensor7_0100000(
    std::ostream &os,
    const Tensor7<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6> &t,
    const int &iterator0 )
{
  os << '[';
  for ( int i = 0; i < Tensor_Dim1 - 1; ++i ) {
    FTensor::Tensor7_0010000( os, t, iterator0, i );
    os << ',';
  }
  FTensor::Tensor7_0010000( os, t, iterator0, Tensor_Dim1 - 1 );
  os << ']';

  return os;
}

template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
          int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5, int Tensor_Dim6>
std::ostream &operator<<(std::ostream &os,
          const FTensor::Tensor7<T, Tensor_Dim0, Tensor_Dim1,Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5, Tensor_Dim6> &t)
{
  os << '[';
  for(int i = 0; i < Tensor_Dim0 - 1; ++i)
    {
      FTensor::Tensor7_0100000(os, t, i);
      os << ',';
    }
  FTensor::Tensor7_0100000(os, t, Tensor_Dim0 - 1);
  os << ']';

  return os;
}

} //namespace FTensor

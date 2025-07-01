/* A general version, not for pointers. */

#pragma once

#include "Tensor6_contracted.hpp"

#include <ostream>
#include "floatmatrixf.h"
#ifdef FTENSOR_DEBUG
#include <sstream>
#include <stdexcept>
#endif

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
            int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5>
  class Tensor6
  {
  protected:
    T data[Tensor_Dim0][Tensor_Dim1][Tensor_Dim2][Tensor_Dim3][Tensor_Dim4][Tensor_Dim5];

  public:
    /* Initialization operators */
    template <class... U> Tensor6(U... d) : data{d...}
    {
      static_assert(sizeof...(d) == sizeof(data) / sizeof(T),
                    "Incorrect number of Arguments. Constructor should "
                    "initialize the entire Tensor");
    }

    Tensor6() {}
    /* There are two operator(int,int,int,int,int)'s, one for non-consts
       that lets you change the value, and one for consts that
       doesn't. */
    T &operator()( const int N1, const int N2, const int N3, const int N4, const int N5, const int N6 )
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0
          || N3 >= Tensor_Dim2 || N3 < 0 || N4 >= Tensor_Dim3 || N4 < 0 || N5 >= Tensor_Dim4 || N5 < 0
          || N6 >= Tensor_Dim5 || N6 < 0 )
        {
          std::stringstream s;
          s << "Bad index in Tensor6<T," << Tensor_Dim0 << "," << Tensor_Dim1
            << "," << Tensor_Dim2 << "," << Tensor_Dim3 << "," << Tensor_Dim4 << Tensor_Dim5 << ">.operator(" << N1
            << "," << N2 << "," << N3 << "," << N4 "," << N5 << "," << N6 << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[N1][N2][N3][N4][N5][N6];
    }

    T operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0
          || N3 >= Tensor_Dim2 || N3 < 0 || N4 >= Tensor_Dim3 || N4 < 0 || N5 >= Tensor_Dim4 || N5 < 0
          || N6 >= Tensor_Dim5 || N6 < 0 )
        {
          std::stringstream s;
        s << "Bad index in Tensor6<T," << Tensor_Dim0 << "," << Tensor_Dim1
          << "," << Tensor_Dim2 << "," << Tensor_Dim3 << "," << Tensor_Dim4 << ">.operator(" << N1
          << "," << N2 << "," << N3 << "," << N4 "," << N5 << "," << N6 << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[N1][N2][N3][N4][N5][N6];
    }

    /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. */

    template <char i, char j, char k, char l, char m, char n, int Dim0, int Dim1, int Dim2,
              int Dim3, int Dim4, int Dim5>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2
            && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 ),
      Tensor6_Expr<Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>,
                   T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>,
               const Index<k, Dim2>, const Index<l, Dim3>, 
               const Index<m, Dim4>, const Index<n, Dim5> )
    {
      return Tensor6_Expr<
        Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T,
        Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>(*this);
    };

    template <char i, char j, char k, char l, char m, char n, int Dim0, int Dim1, int Dim2,
              int Dim3, int Dim4, int Dim5>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2
            && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5 ),
      Tensor6_Expr<
        const Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>,
        T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>,
        const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<n, Dim5> ) const
    {
      return Tensor6_Expr<
        const Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>,
        T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>(*this);
    };

    /* These operators are for internal contractions, resulting in a Tensor4.
     * For example something like A(i,i,k,l,m,n) */

  // (i, i, k, l, m, n, )
    template <char i, char k, char l, char m, char n, int Dim01, int Dim2, int Dim3, int Dim4, int Dim5>
    auto operator()( const Index<i, Dim01>, const Index<i, Dim01>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<n, Dim5> ) const
    {
      static_assert( Tensor_Dim0 >= Dim01 && Tensor_Dim1 >= Dim01 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_01<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim01>;
      return Tensor4_Expr<TensorExpr, T, Dim2, Dim3, Dim4, Dim5, k, l, m, n>( TensorExpr( *this ) );
    };

    // (i, j, i, l, m, n, )
    template <char i, char j, char l, char m, char n, int Dim02, int Dim1, int Dim3, int Dim4, int Dim5>
    auto operator()( const Index<i, Dim02>, const Index<j, Dim1>, const Index<i, Dim02>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<n, Dim5> ) const
    {
      static_assert( Tensor_Dim0 >= Dim02 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim02 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_02<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim02>;
      return Tensor4_Expr<TensorExpr, T, Dim1, Dim3, Dim4, Dim5, j, l, m, n>( TensorExpr( *this ) );
    };

    // (i, j, k, i, m, n, )
    template <char i, char j, char k, char m, char n, int Dim03, int Dim1, int Dim2, int Dim4, int Dim5>
    auto operator()( const Index<i, Dim03>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<i, Dim03>, const Index<m, Dim4>, const Index<n, Dim5> ) const
    {
      static_assert( Tensor_Dim0 >= Dim03 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim03 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_03<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim03>;
      return Tensor4_Expr<TensorExpr, T, Dim1, Dim2, Dim4, Dim5, j, k, m, n>( TensorExpr( *this ) );
    };

    // (i, j, k, l, i, n, )
    template <char i, char j, char k, char l, char n, int Dim04, int Dim1, int Dim2, int Dim3, int Dim5>
    auto operator()( const Index<i, Dim04>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<i, Dim04>, const Index<n, Dim5> ) const
    {
      static_assert( Tensor_Dim0 >= Dim04 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim04 && Tensor_Dim5 >= Dim5,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_04<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim04>;
      return Tensor4_Expr<TensorExpr, T, Dim1, Dim2, Dim3, Dim5, j, k, l, n>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, i, )
    template <char i, char j, char k, char l, char m, int Dim05, int Dim1, int Dim2, int Dim3, int Dim4>
    auto operator()( const Index<i, Dim05>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<i, Dim05> ) const
    {
      static_assert( Tensor_Dim0 >= Dim05 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim05,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_05<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim05>;
      return Tensor4_Expr<TensorExpr, T, Dim1, Dim2, Dim3, Dim4, j, k, l, m>( TensorExpr( *this ) );
    };

    // (i, j, j, l, m, n, )
    template <char i, char j, char l, char m, char n, int Dim0, int Dim12, int Dim3, int Dim4, int Dim5>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim12>, const Index<j, Dim12>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<n, Dim5> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim12 && Tensor_Dim2 >= Dim12 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_12<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim12>;
      return Tensor4_Expr<TensorExpr, T, Dim0, Dim3, Dim4, Dim5, i, l, m, n>( TensorExpr( *this ) );
    };

    // (i, j, k, j, m, n, )
    template <char i, char j, char k, char m, char n, int Dim0, int Dim13, int Dim2, int Dim4, int Dim5>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim13>, const Index<k, Dim2>, const Index<j, Dim13>, const Index<m, Dim4>, const Index<n, Dim5> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim13 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim13 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_13<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim13>;
      return Tensor4_Expr<TensorExpr, T, Dim0, Dim2, Dim4, Dim5, i, k, m, n>( TensorExpr( *this ) );
    };

    // (i, j, k, l, j, n, )
    template <char i, char j, char k, char l, char n, int Dim0, int Dim14, int Dim2, int Dim3, int Dim5>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim14>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<j, Dim14>, const Index<n, Dim5> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim14 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim14 && Tensor_Dim5 >= Dim5,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_14<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim14>;
      return Tensor4_Expr<TensorExpr, T, Dim0, Dim2, Dim3, Dim5, i, k, l, n>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, j, )
    template <char i, char j, char k, char l, char m, int Dim0, int Dim15, int Dim2, int Dim3, int Dim4>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim15>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<j, Dim15> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim15 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim15,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_15<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim15>;
      return Tensor4_Expr<TensorExpr, T, Dim0, Dim2, Dim3, Dim4, i, k, l, m>( TensorExpr( *this ) );
    };

    // (i, j, k, k, m, n, )
    template <char i, char j, char k, char m, char n, int Dim0, int Dim1, int Dim23, int Dim4, int Dim5>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim23>, const Index<k, Dim23>, const Index<m, Dim4>, const Index<n, Dim5> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim23 && Tensor_Dim3 >= Dim23 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim5,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_23<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim23>;
      return Tensor4_Expr<TensorExpr, T, Dim0, Dim1, Dim4, Dim5, i, j, m, n>( TensorExpr( *this ) );
    };

    // (i, j, k, l, k, n, )
    template <char i, char j, char k, char l, char n, int Dim0, int Dim1, int Dim24, int Dim3, int Dim5>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim24>, const Index<l, Dim3>, const Index<k, Dim24>, const Index<n, Dim5> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim24 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim24 && Tensor_Dim5 >= Dim5,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_24<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim24>;
      return Tensor4_Expr<TensorExpr, T, Dim0, Dim1, Dim3, Dim5, i, j, l, n>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, k, )
    template <char i, char j, char k, char l, char m, int Dim0, int Dim1, int Dim25, int Dim3, int Dim4>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim25>, const Index<l, Dim3>, const Index<m, Dim4>, const Index<k, Dim25> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim25 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim25,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_25<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim25>;
      return Tensor4_Expr<TensorExpr, T, Dim0, Dim1, Dim3, Dim4, i, j, l, m>( TensorExpr( *this ) );
    };

    // (i, j, k, l, l, n, )
    template <char i, char j, char k, char l, char n, int Dim0, int Dim1, int Dim2, int Dim34, int Dim5>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim34>, const Index<l, Dim34>, const Index<n, Dim5> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim34 && Tensor_Dim4 >= Dim34 && Tensor_Dim5 >= Dim5,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_34<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim34>;
      return Tensor4_Expr<TensorExpr, T, Dim0, Dim1, Dim2, Dim5, i, j, k, n>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, l, )
    template <char i, char j, char k, char l, char m, int Dim0, int Dim1, int Dim2, int Dim35, int Dim4>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim35>, const Index<m, Dim4>, const Index<l, Dim35> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim35 && Tensor_Dim4 >= Dim4 && Tensor_Dim5 >= Dim35,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_35<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim35>;
      return Tensor4_Expr<TensorExpr, T, Dim0, Dim1, Dim2, Dim4, i, j, k, m>( TensorExpr( *this ) );
    };

    // (i, j, k, l, m, m, )
    template <char i, char j, char k, char l, char m, int Dim0, int Dim1, int Dim2, int Dim3, int Dim45>
    auto operator()( const Index<i, Dim0>, const Index<j, Dim1>, const Index<k, Dim2>, const Index<l, Dim3>, const Index<m, Dim45>, const Index<m, Dim45> ) const
    {
      static_assert( Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3 && Tensor_Dim4 >= Dim45 && Tensor_Dim5 >= Dim45,
          "Incompatible indices" );

      using TensorExpr = Tensor6_contracted_45<
          Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5>, T, Dim45>;
      return Tensor4_Expr<TensorExpr, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l>( TensorExpr( *this ) );
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
            int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5>
  std::ostream &Tensor6_000001(
    std::ostream &os,
    const Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5> &t,
    const int &iterator0, const int &iterator1, const int &iterator2, const int &iterator3, const int &iterator4)
  {
    os << '[';
    for(int i = 0; i < Tensor_Dim5 - 1; ++i)
      {
        os << t(iterator0, iterator1, iterator2, iterator3, iterator4, i);
        os << ',';
      }
    os << t(iterator0, iterator1, iterator2, iterator3, iterator4, Tensor_Dim5 - 1);
    os << ']';

    return os;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
            int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5>
  std::ostream &Tensor6_000010(
    std::ostream &os,
    const Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5> &t,
    const int &iterator0, const int &iterator1, const int &iterator2, const int &iterator3)
  {
    os << '[';
    for(int i = 0; i < Tensor_Dim4 - 1; ++i)
      {
        FTensor::Tensor6_000001(os, t, iterator0, iterator1, iterator2, iterator3, i);
        os << ',';
      }
    FTensor::Tensor6_000001(os, t, iterator0, iterator1, iterator2, Tensor_Dim4 - 1);
    os << ']';

    return os;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
            int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5>
  std::ostream &Tensor6_000100(
    std::ostream &os,
    const Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5> &t,
    const int &iterator0, const int &iterator1, const int &iterator2)
  {
    os << '[';
    for(int i = 0; i < Tensor_Dim3 - 1; ++i)
      {
        FTensor::Tensor6_000010(os, t, iterator0, iterator1, iterator2, i);
        os << ',';
      }
    FTensor::Tensor6_000010(os, t, iterator0, iterator1, iterator2, Tensor_Dim3 - 1);
    os << ']';

    return os;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
    int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5>
std::ostream &Tensor6_001000(
    std::ostream &os,
    const Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5> &t,
    const int &iterator0, const int &iterator1 )
{
  os << '[';
  for ( int i = 0; i < Tensor_Dim2 - 1; ++i ) {
    FTensor::Tensor6_000100( os, t, iterator0, iterator1, i );
    os << ',';
  }
  FTensor::Tensor6_000100( os, t, iterator0, iterator1, Tensor_Dim2 - 1 );
  os << ']';

  return os;
}

 template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
    int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5>
std::ostream &Tensor6_010000(
    std::ostream &os,
    const Tensor6<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5> &t,
    const int &iterator0 )
{
  os << '[';
  for ( int i = 0; i < Tensor_Dim1 - 1; ++i ) {
    FTensor::Tensor6_001000( os, t, iterator0, i );
    os << ',';
  }
  FTensor::Tensor6_001000( os, t, iterator0, Tensor_Dim1 - 1 );
  os << ']';

  return os;
}

template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
          int Tensor_Dim3, int Tensor_Dim4, int Tensor_Dim5>
std::ostream &operator<<(std::ostream &os,
                         const FTensor::Tensor6<T, Tensor_Dim0, Tensor_Dim1,
                                                Tensor_Dim2, Tensor_Dim3, Tensor_Dim4, Tensor_Dim5> &t)
{
  os << '[';
  for(int i = 0; i < Tensor_Dim0 - 1; ++i)
    {
      FTensor::Tensor6_010000(os, t, i);
      os << ',';
    }
  FTensor::Tensor6_010000(os, t, Tensor_Dim0 - 1);
  os << ']';

  return os;
}

} //namespace FTensor

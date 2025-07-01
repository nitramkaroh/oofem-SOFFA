/* Multiplies a Tensor6 with a generic, yielding a Tensor6. */

#pragma once

namespace FTensor
{
  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
  class Tensor6_times_generic
  {
    const Tensor6_Expr<A, T, Dim0, Dim1,  Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
    const U d;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
    {
      return iterA(N1, N2, N3, N4, N5, N6) * d;
    }

    Tensor6_times_generic(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
                          const U &d0)
        : iterA(a), d(d0)
    {}
  };

  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
  Tensor6_Expr<Tensor6_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>,
               typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>
  operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a, const U &d0)
  {
    using TensorExpr = Tensor6_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
    return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i,
                        j, k, l, m, n>(TensorExpr(a, d0));
  }

  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
  Tensor6_Expr<Tensor6_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>,
               typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>
  operator*(const U &d0, const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
  {
    using TensorExpr = Tensor6_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
    return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i,
                        j, k, l, m, n>(TensorExpr(a, d0));
  }
}

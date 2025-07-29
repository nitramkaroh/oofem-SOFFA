/* Multiplies a Tensor5 with a generic, yielding a Tensor5. */

#pragma once

namespace FTensor
{
  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
  class Tensor5_times_generic
  {
    const Tensor5_Expr<A, T, Dim0, Dim1,  Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
    const U d;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5) const
    {
      return iterA(N1, N2, N3, N4, N5) * d;
    }

    Tensor5_times_generic(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
                          const U &d0)
        : iterA(a), d(d0)
    {}
  };

  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
  Tensor5_Expr<Tensor5_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>,
               typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>
  operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a, const U &d0)
  {
    using TensorExpr = Tensor5_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
    return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i,
                        j, k, l, m>(TensorExpr(a, d0));
  }

  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
  Tensor5_Expr<Tensor5_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>,
               typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>
  operator*(const U &d0, const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
  {
    using TensorExpr = Tensor5_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
    return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i,
                        j, k, l, m>(TensorExpr(a, d0));
  }
}

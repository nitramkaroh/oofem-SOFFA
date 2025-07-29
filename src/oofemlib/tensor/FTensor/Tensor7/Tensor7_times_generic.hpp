/* Multiplies a Tensor7 with a generic, yielding a Tensor7. */

#pragma once

namespace FTensor
{
  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
  class Tensor7_times_generic
  {
    const Tensor7_Expr<A, T, Dim0, Dim1,  Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
    const U d;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const int N7) const
    {
      return iterA(N1, N2, N3, N4, N5, N6, N7) * d;
    }

    Tensor7_times_generic(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
                          const U &d0)
        : iterA(a), d(d0)
    {}
  };

  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
  Tensor7_Expr<Tensor7_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>,
               typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>
  operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a, const U &d0)
  {
    using TensorExpr = Tensor7_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
    return Tensor7_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i,
                        j, k, l, m, n, o>(TensorExpr(a, d0));
  }

  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
  Tensor7_Expr<Tensor7_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>,
               typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>
  operator*(const U &d0, const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
  {
    using TensorExpr = Tensor7_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
    return Tensor7_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i,
                        j, k, l, m, n, o>(TensorExpr(a, d0));
  }
}

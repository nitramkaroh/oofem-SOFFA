#pragma once

namespace FTensor
{
// A(i, j, k, l, )*B(m, n, o, ) -> Tensor7
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
class Tensor4_times_Tensor3_outer
{
  Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
  Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> iterB;

public:
  Tensor4_times_Tensor3_outer(
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
    const Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> &b)
      : iterA(a), iterB(b)
  {}
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const int N7) const
  {
    return iterA(N1, N2, N3, N4) * iterB(N5, N6, N7);
  }
};

// A(i, j, k, l, )*B(m, n, o, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor7_Expr<
  Tensor4_times_Tensor3_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>
operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
          const Tensor3_Expr<A, T, Dim4, Dim5, Dim6, m, n, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor7_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>(TensorExpr(a, b));
}

// B(m, n, o, )*A(i, j, k, l, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor7_Expr<
  Tensor4_times_Tensor3_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>
operator*(const Tensor3_Expr<A, T, Dim4, Dim5, Dim6, m, n, o> &b,
          const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

} //namespace FTensor

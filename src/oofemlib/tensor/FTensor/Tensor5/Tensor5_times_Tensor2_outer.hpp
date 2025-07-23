#pragma once

namespace FTensor
{
// A(i, j, k, l, m, )*B(n, o, ) -> Tensor7
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
class Tensor5_times_Tensor2_outer
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor2_Expr<B, U, Dim5, Dim6, n, o> iterB;

public:
  Tensor5_times_Tensor2_outer(
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
    const Tensor2_Expr<B, U, Dim5, Dim6, n, o> &b)
      : iterA(a), iterB(b)
  {}
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const int N7) const
  {
    return iterA(N1, N2, N3, N4, N5) * iterB(N6, N7);
  }
};

// A(i, j, k, l, m, )*B(n, o, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor7_Expr<
  Tensor5_times_Tensor2_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>
operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
          const Tensor2_Expr<B, U, Dim5, Dim6, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor2_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor7_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>(TensorExpr(a, b));
}

// B(n, o, )*A(i, j, k, l, m, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor7_Expr<
  Tensor5_times_Tensor2_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>
operator*(const Tensor2_Expr<B, U, Dim5, Dim6, n, o> &b,
          const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

} //namespace FTensor

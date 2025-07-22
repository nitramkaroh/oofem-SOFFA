#pragma once

namespace FTensor
{
// A(i, j, k, l, m, n, )*B(o, ) -> Tensor7
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
class Tensor6_times_Tensor1_outer
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor1_Expr<B, U, Dim6, o> iterB;

public:
  Tensor6_times_Tensor1_outer(
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor1_Expr<B, U, Dim6, o> &b)
      : iterA(a), iterB(b)
  {}
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const int N7) const
  {
    return iterA(N1, N2, N3, N4, N5, N6) * iterB(N7);
  }
};

// A(i, j, k, l, m, n, )*B(o, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor7_Expr<
  Tensor6_times_Tensor1_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>
operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
          const Tensor1_Expr<A, T, Dim6, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor1_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor7_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>(TensorExpr(a, b));
}

// B(o, )*A(i, j, k, l, m, n, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor7_Expr<
  Tensor6_times_Tensor1_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>
operator*(const Tensor1_Expr<A, T, Dim6, o> &b,
          const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

} //namespace FTensor

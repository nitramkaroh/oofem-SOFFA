#pragma once

namespace FTensor
{
// A(i, j, k, l, m, )*B(n, ) -> Tensor6
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
class Tensor5_times_Tensor1_outer
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor1_Expr<B, U, Dim5, n> iterB;

public:
  Tensor5_times_Tensor1_outer(
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
    const Tensor1_Expr<B, U, Dim5, n> &b)
      : iterA(a), iterB(b)
  {}
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    return iterA(N1, N2, N3, N4, N5) * iterB(N6);
  }
};

// A(i, j, k, l, m, )*B(n, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor6_Expr<
  Tensor5_times_Tensor1_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>
operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
          const Tensor1_Expr<A, T, Dim5, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor1_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>(TensorExpr(a, b));
}

// B(n, )*A(i, j, k, l, m, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor6_Expr<
  Tensor5_times_Tensor1_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>
operator*(const Tensor1_Expr<A, T, Dim5, n> &b,
          const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

} //namespace FTensor

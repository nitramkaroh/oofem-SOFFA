#pragma once

namespace FTensor
{
// A(i, j, k, l, )*B(m, n, ) -> Tensor6
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
class Tensor4_times_Tensor2_outer
{
  Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
  Tensor2_Expr<B, U, Dim4, Dim5, m, n> iterB;

public:
  Tensor4_times_Tensor2_outer(
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
    const Tensor2_Expr<B, U, Dim4, Dim5, m, n> &b)
      : iterA(a), iterB(b)
  {}
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    return iterA(N1, N2, N3, N4) * iterB(N5, N6);
  }
};

// A(i, j, k, l, )*B(m, n, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor6_Expr<
  Tensor4_times_Tensor2_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>
operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
          const Tensor2_Expr<B, U, Dim4, Dim5, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor2_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>(TensorExpr(a, b));
}

// B(m, n, )*A(i, j, k, l, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor6_Expr<
  Tensor4_times_Tensor2_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>
operator*(const Tensor2_Expr<B, U, Dim4, Dim5, m, n> &b,
          const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

} //namespace FTensor

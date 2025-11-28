#pragma once

namespace FTensor
{
// A(i, j, k, )*B(l, m, ) -> Tensor5
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
class Tensor3_times_Tensor2_outer
{
  Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
  Tensor2_Expr<B, U, Dim3, Dim4, l, m> iterB;

public:
  Tensor3_times_Tensor2_outer(
    const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
    const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b)
      : iterA(a), iterB(b)
  {}
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5) const
  {
    return iterA(N1, N2, N3) * iterB(N4, N5);
  }
};

// A(i, j, k, )*B(l, m, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor5_Expr<
  Tensor3_times_Tensor2_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>
operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
          const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor2_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>(TensorExpr(a, b));
}

// B(l, m, )*A(i, j, k, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor5_Expr<
  Tensor3_times_Tensor2_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>
operator*(const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b,
          const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
{
  return a * b;
}

} //namespace FTensor

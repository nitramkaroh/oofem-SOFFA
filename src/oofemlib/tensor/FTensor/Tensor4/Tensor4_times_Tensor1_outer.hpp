#pragma once

namespace FTensor
{
// A(i, j, k, l, )*B(m, ) -> Tensor5
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
class Tensor4_times_Tensor1_outer
{
  Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
  Tensor1_Expr<B, U, Dim4, m> iterB;

public:
  Tensor4_times_Tensor1_outer(
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
    const Tensor1_Expr<B, U, Dim4, m> &b)
      : iterA(a), iterB(b)
  {}
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5) const
  {
    return iterA(N1, N2, N3, N4) * iterB(N5);
  }
};

// A(i, j, k, l, )*B(m, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor5_Expr<
  Tensor4_times_Tensor1_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>
operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
          const Tensor1_Expr<A, T, Dim4, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor1_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>(TensorExpr(a, b));
}

// B(m, )*A(i, j, k, l, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor5_Expr<
  Tensor4_times_Tensor1_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>
operator*(const Tensor1_Expr<A, T, Dim4, m> &b,
          const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

} //namespace FTensor

#pragma once

namespace FTensor
{
// A(i, j, k, )*B(l, ) -> Tensor4
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
class Tensor3_times_Tensor1_outer
{
  Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
  Tensor1_Expr<B, U, Dim3, l> iterB;

public:
  Tensor3_times_Tensor1_outer(
    const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
    const Tensor1_Expr<B, U, Dim3, l> &b)
      : iterA(a), iterB(b)
  {}
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4) const
  {
    return iterA(N1, N2, N3) * iterB(N4);
  }
};

// A(i, j, k, )*B(l, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
Tensor4_Expr<
  Tensor3_times_Tensor1_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>
operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
          const Tensor1_Expr<B, U, Dim3, l> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor1_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
}

// B(l, )*A(i, j, k, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
Tensor4_Expr<
  Tensor3_times_Tensor1_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>
operator*(const Tensor1_Expr<B, U, Dim3, l> &b,
          const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
{
  return a * b;
}

} //namespace FTensor

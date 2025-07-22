#pragma once

namespace FTensor
{
// A(i, j, k, )*B(l, m, n, ) -> Tensor6
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
class Tensor3_times_Tensor3_outer
{
  Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
  Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> iterB;

public:
  Tensor3_times_Tensor3_outer(
    const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> &b)
      : iterA(a), iterB(b)
  {}
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    return iterA(N1, N2, N3) * iterB(N4, N5, N6);
  }
};

// A(i, j, k, )*B(l, m, n, ) operator:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor6_Expr<
  Tensor3_times_Tensor3_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>,
  typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>
operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
          const Tensor3_Expr<A, T, Dim3, Dim4, Dim5, l, m, n> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_outer<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>(TensorExpr(a, b));
}

} //namespace FTensor

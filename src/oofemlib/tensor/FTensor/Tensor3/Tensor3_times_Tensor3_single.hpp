#pragma once

namespace FTensor
{
//Universal A(i, j, k, )*B(l, m, n, ) single contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n, int DimA, int DimB, int DimX, int DimY, int DimZ, char a, char b, char x, char y, char z>
class Tensor3_times_Tensor3_single
{
  Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
  Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      {
        // Permutation is where the indices get checked.
        result += Permutation3<DimA, DimB, DimX, a, b, x>().eval(
          iterA, N1, N2, xx) * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx, N3, N4);
      }
    return result;
  }

  Tensor3_times_Tensor3_single(
    const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &iter_a,
    const Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, )*B(i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim4, i, l, m> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim0, Dim3, Dim4, i, j, k, i, l, m, Dim1, Dim2, Dim0, Dim3, Dim4, j, k, i, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>(TensorExpr(a, b));
}

//A(i, j, k, )*B(l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim4, l, i, m> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, i, j, k, l, i, m, Dim1, Dim2, Dim0, Dim3, Dim4, j, k, i, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>(TensorExpr(a, b));
}

//A(i, j, k, )*B(l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim0, l, m, i> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, i, j, k, l, m, i, Dim1, Dim2, Dim0, Dim3, Dim4, j, k, i, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>(TensorExpr(a, b));
}

//A(i, j, k, )*B(j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim4, j, l, m> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim1, Dim3, Dim4, i, j, k, j, l, m, Dim0, Dim2, Dim1, Dim3, Dim4, i, k, j, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>(TensorExpr(a, b));
}

//A(i, j, k, )*B(l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim4, l, j, m> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, i, j, k, l, j, m, Dim0, Dim2, Dim1, Dim3, Dim4, i, k, j, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>(TensorExpr(a, b));
}

//A(i, j, k, )*B(l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim1, l, m, j> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, i, j, k, l, m, j, Dim0, Dim2, Dim1, Dim3, Dim4, i, k, j, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>(TensorExpr(a, b));
}

//A(i, j, k, )*B(k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim2, Dim3, Dim4, i, j, k, k, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>(TensorExpr(a, b));
}

//A(i, j, k, )*B(l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim4, l, k, m> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, i, j, k, l, k, m, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>(TensorExpr(a, b));
}

//A(i, j, k, )*B(l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim2, l, m, k> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, i, j, k, l, m, k, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>(TensorExpr(a, b));
}

} //namespace FTensor

#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, )*B(n, o, ) single contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o, int DimA, int DimB, int DimC, int DimD, int DimX, int DimY, char a, char b, char c, char d, char x, char y>
class Tensor5_times_Tensor2_single
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor2_Expr<B, U, Dim5, Dim6, n, o> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      {
        // Permutation is where the indices get checked.
        result += Permutation5<DimA, DimB, DimC, DimD, DimX, a, b, c, d, x>().eval(
          iterA, N1, N2, N3, N4, xx) * Permutation2<DimX, DimY, x, y>().eval(iterB, xx, N5);
      }
    return result;
  }

  Tensor5_times_Tensor2_single(
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &iter_a,
    const Tensor2_Expr<B, U, Dim5, Dim6, n, o> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, )*B(i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor2_Expr<B, U, Dim0, Dim5, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, i, j, k, l, m, i, n, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, j, k, l, m, i, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, j, k, l, m, n>(TensorExpr(a, b));
}

//B(i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor2_Expr<B, U, Dim0, Dim5, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor2_Expr<B, U, Dim5, Dim0, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, i, j, k, l, m, n, i, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, j, k, l, m, i, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, j, k, l, m, n>(TensorExpr(a, b));
}

//B(n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor2_Expr<B, U, Dim5, Dim0, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor2_Expr<B, U, Dim1, Dim5, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, i, j, k, l, m, j, n, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, i, k, l, m, j, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, i, k, l, m, n>(TensorExpr(a, b));
}

//B(j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor2_Expr<B, U, Dim1, Dim5, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor2_Expr<B, U, Dim5, Dim1, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, i, j, k, l, m, n, j, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, i, k, l, m, j, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, i, k, l, m, n>(TensorExpr(a, b));
}

//B(n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor2_Expr<B, U, Dim5, Dim1, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor2_Expr<B, U, Dim2, Dim5, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, i, j, k, l, m, k, n, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, i, j, l, m, k, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, i, j, l, m, n>(TensorExpr(a, b));
}

//B(k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor2_Expr<B, U, Dim2, Dim5, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor2_Expr<B, U, Dim5, Dim2, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, i, j, k, l, m, n, k, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, i, j, l, m, k, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, i, j, l, m, n>(TensorExpr(a, b));
}

//B(n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor2_Expr<B, U, Dim5, Dim2, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor2_Expr<B, U, Dim3, Dim5, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, i, j, k, l, m, l, n, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, i, j, k, m, l, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, i, j, k, m, n>(TensorExpr(a, b));
}

//B(l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor2_Expr<B, U, Dim3, Dim5, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor2_Expr<B, U, Dim5, Dim3, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, i, j, k, l, m, n, l, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, i, j, k, m, l, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, i, j, k, m, n>(TensorExpr(a, b));
}

//B(n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor2_Expr<B, U, Dim5, Dim3, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor2_Expr<B, U, Dim4, Dim5, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, i, j, k, l, m, m, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, i, j, k, l, n>(TensorExpr(a, b));
}

//B(m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor2_Expr<B, U, Dim4, Dim5, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor2_Expr<B, U, Dim5, Dim4, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, i, j, k, l, m, n, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, i, j, k, l, n>(TensorExpr(a, b));
}

//B(n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor2_Expr<B, U, Dim5, Dim4, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

} //namespace FTensor

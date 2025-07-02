#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, n, )*B(o, p, q, ) double contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q, int DimA, int DimB, int DimC, int DimD, int DimX, int DimY, int DimZ, char a, char b, char c, char d, char x, char y, char z>
class Tensor6_times_Tensor3_double
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor3_Expr<B, U, Dim6, Dim7, Dim8, o, p, q> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        {
          // Permutation is where the indices get checked.
          result += Permutation6<DimA, DimB, DimC, DimD, DimX, DimY, a, b, c, d, x, y>().eval(
            iterA, N1, N2, N3, N4, xx, yy) * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx, yy, N5);
        }
    return result;
  }

  Tensor6_times_Tensor3_double(
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &iter_a,
    const Tensor3_Expr<B, U, Dim6, Dim7, Dim8, o, p, q> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, n, )*B(j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim6, j, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, i, j, k, l, m, n, j, i, o, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, k, l, m, n, j, i, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(j, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim6, j, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim6, i, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, i, j, k, l, m, n, i, j, o, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, k, l, m, n, i, j, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(i, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim6, i, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim6, Dim0, j, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, i, j, k, l, m, n, j, o, i, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, k, l, m, n, j, o, i>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(j, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim6, Dim0, j, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim6, Dim1, i, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, i, j, k, l, m, n, i, o, j, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, k, l, m, n, i, o, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(i, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim6, Dim1, i, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim1, Dim0, o, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, i, j, k, l, m, n, o, j, i, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, k, l, m, n, o, j, i>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(o, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim1, Dim0, o, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim0, Dim1, o, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, i, j, k, l, m, n, o, i, j, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, k, l, m, n, o, i, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(o, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim0, Dim1, o, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim6, k, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, i, j, k, l, m, n, k, i, o, Dim1, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, j, l, m, n, k, i, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(k, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim6, k, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim6, i, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, i, j, k, l, m, n, i, k, o, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, j, l, m, n, i, k, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(i, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim6, i, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim6, Dim0, k, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, i, j, k, l, m, n, k, o, i, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, j, l, m, n, k, o, i>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(k, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim6, Dim0, k, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim6, Dim2, i, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, i, j, k, l, m, n, i, o, k, Dim1, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, j, l, m, n, i, o, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(i, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim6, Dim2, i, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim2, Dim0, o, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, i, j, k, l, m, n, o, k, i, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, j, l, m, n, o, k, i>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(o, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim2, Dim0, o, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim0, Dim2, o, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, i, j, k, l, m, n, o, i, k, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, j, l, m, n, o, i, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(o, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim0, Dim2, o, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim6, l, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim6, i, j, k, l, m, n, l, i, o, Dim1, Dim2, Dim4, Dim5, Dim3, Dim0, Dim6, j, k, m, n, l, i, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(l, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim6, l, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim6, i, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim6, i, j, k, l, m, n, i, l, o, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, j, k, m, n, i, l, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(i, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim6, i, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim6, Dim0, l, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim0, i, j, k, l, m, n, l, o, i, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim0, j, k, m, n, l, o, i>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(l, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim6, Dim0, l, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim6, Dim3, i, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim3, i, j, k, l, m, n, i, o, l, Dim1, Dim2, Dim4, Dim5, Dim0, Dim6, Dim3, j, k, m, n, i, o, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(i, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim6, Dim3, i, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim3, Dim0, o, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, i, j, k, l, m, n, o, l, i, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, Dim0, j, k, m, n, o, l, i>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(o, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim3, Dim0, o, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim0, Dim3, o, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, i, j, k, l, m, n, o, i, l, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, j, k, m, n, o, i, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(o, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim0, Dim3, o, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim6, m, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim6, i, j, k, l, m, n, m, i, o, Dim1, Dim2, Dim3, Dim5, Dim4, Dim0, Dim6, j, k, l, n, m, i, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(m, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim6, m, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim6, i, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim6, i, j, k, l, m, n, i, m, o, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, j, k, l, n, i, m, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(i, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim6, i, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim6, Dim0, m, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim0, i, j, k, l, m, n, m, o, i, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim0, j, k, l, n, m, o, i>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(m, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim6, Dim0, m, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim6, Dim4, i, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim4, i, j, k, l, m, n, i, o, m, Dim1, Dim2, Dim3, Dim5, Dim0, Dim6, Dim4, j, k, l, n, i, o, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(i, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim6, Dim4, i, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim4, Dim0, o, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, i, j, k, l, m, n, o, m, i, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, Dim0, j, k, l, n, o, m, i>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(o, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim4, Dim0, o, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim0, Dim4, o, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, i, j, k, l, m, n, o, i, m, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, j, k, l, n, o, i, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(o, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim0, Dim4, o, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim6, n, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim6, i, j, k, l, m, n, n, i, o, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, j, k, l, m, n, i, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, j, k, l, m, o>(TensorExpr(a, b));
}

//B(n, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim6, n, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim6, i, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim6, i, j, k, l, m, n, i, n, o, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, j, k, l, m, i, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, j, k, l, m, o>(TensorExpr(a, b));
}

//B(i, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim6, i, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim0, n, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim0, i, j, k, l, m, n, n, o, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, j, k, l, m, n, o, i>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, j, k, l, m, o>(TensorExpr(a, b));
}

//B(n, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim0, n, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim6, Dim5, i, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim5, i, j, k, l, m, n, i, o, n, Dim1, Dim2, Dim3, Dim4, Dim0, Dim6, Dim5, j, k, l, m, i, o, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, j, k, l, m, o>(TensorExpr(a, b));
}

//B(i, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim6, Dim5, i, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim5, Dim0, o, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, i, j, k, l, m, n, o, n, i, Dim1, Dim2, Dim3, Dim4, Dim6, Dim5, Dim0, j, k, l, m, o, n, i>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, j, k, l, m, o>(TensorExpr(a, b));
}

//B(o, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim5, Dim0, o, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim0, Dim5, o, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, i, j, k, l, m, n, o, i, n, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, j, k, l, m, o, i, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, j, k, l, m, o>(TensorExpr(a, b));
}

//B(o, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim0, Dim5, o, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim6, k, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, i, j, k, l, m, n, k, j, o, Dim0, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, i, l, m, n, k, j, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(k, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim6, k, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim6, j, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, i, j, k, l, m, n, j, k, o, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, i, l, m, n, j, k, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(j, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim6, j, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim6, Dim1, k, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, i, j, k, l, m, n, k, o, j, Dim0, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, i, l, m, n, k, o, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(k, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim6, Dim1, k, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim6, Dim2, j, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, i, j, k, l, m, n, j, o, k, Dim0, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, i, l, m, n, j, o, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(j, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim6, Dim2, j, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim2, Dim1, o, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, i, j, k, l, m, n, o, k, j, Dim0, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, i, l, m, n, o, k, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(o, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim2, Dim1, o, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim1, Dim2, o, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, i, j, k, l, m, n, o, j, k, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, i, l, m, n, o, j, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(o, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim1, Dim2, o, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim6, l, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim6, i, j, k, l, m, n, l, j, o, Dim0, Dim2, Dim4, Dim5, Dim3, Dim1, Dim6, i, k, m, n, l, j, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(l, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim6, l, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim6, j, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim6, i, j, k, l, m, n, j, l, o, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, i, k, m, n, j, l, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(j, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim6, j, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim6, Dim1, l, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim1, i, j, k, l, m, n, l, o, j, Dim0, Dim2, Dim4, Dim5, Dim3, Dim6, Dim1, i, k, m, n, l, o, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(l, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim6, Dim1, l, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim6, Dim3, j, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim3, i, j, k, l, m, n, j, o, l, Dim0, Dim2, Dim4, Dim5, Dim1, Dim6, Dim3, i, k, m, n, j, o, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(j, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim6, Dim3, j, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim3, Dim1, o, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, i, j, k, l, m, n, o, l, j, Dim0, Dim2, Dim4, Dim5, Dim6, Dim3, Dim1, i, k, m, n, o, l, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(o, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim3, Dim1, o, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim1, Dim3, o, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, i, j, k, l, m, n, o, j, l, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, i, k, m, n, o, j, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(o, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim1, Dim3, o, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim6, m, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim6, i, j, k, l, m, n, m, j, o, Dim0, Dim2, Dim3, Dim5, Dim4, Dim1, Dim6, i, k, l, n, m, j, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(m, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim6, m, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim6, j, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim6, i, j, k, l, m, n, j, m, o, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, i, k, l, n, j, m, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(j, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim6, j, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim6, Dim1, m, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim1, i, j, k, l, m, n, m, o, j, Dim0, Dim2, Dim3, Dim5, Dim4, Dim6, Dim1, i, k, l, n, m, o, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(m, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim6, Dim1, m, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim6, Dim4, j, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim4, i, j, k, l, m, n, j, o, m, Dim0, Dim2, Dim3, Dim5, Dim1, Dim6, Dim4, i, k, l, n, j, o, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(j, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim6, Dim4, j, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim4, Dim1, o, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, i, j, k, l, m, n, o, m, j, Dim0, Dim2, Dim3, Dim5, Dim6, Dim4, Dim1, i, k, l, n, o, m, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(o, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim4, Dim1, o, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim1, Dim4, o, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, i, j, k, l, m, n, o, j, m, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, i, k, l, n, o, j, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(o, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim1, Dim4, o, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim6, n, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim6, i, j, k, l, m, n, n, j, o, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, i, k, l, m, n, j, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, i, k, l, m, o>(TensorExpr(a, b));
}

//B(n, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim6, n, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim6, j, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim6, i, j, k, l, m, n, j, n, o, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, i, k, l, m, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, i, k, l, m, o>(TensorExpr(a, b));
}

//B(j, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim6, j, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim1, n, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim1, i, j, k, l, m, n, n, o, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, i, k, l, m, n, o, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, i, k, l, m, o>(TensorExpr(a, b));
}

//B(n, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim1, n, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim6, Dim5, j, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim5, i, j, k, l, m, n, j, o, n, Dim0, Dim2, Dim3, Dim4, Dim1, Dim6, Dim5, i, k, l, m, j, o, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, i, k, l, m, o>(TensorExpr(a, b));
}

//B(j, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim6, Dim5, j, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim5, Dim1, o, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, i, j, k, l, m, n, o, n, j, Dim0, Dim2, Dim3, Dim4, Dim6, Dim5, Dim1, i, k, l, m, o, n, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, i, k, l, m, o>(TensorExpr(a, b));
}

//B(o, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim5, Dim1, o, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim1, Dim5, o, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, i, j, k, l, m, n, o, j, n, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, i, k, l, m, o, j, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, i, k, l, m, o>(TensorExpr(a, b));
}

//B(o, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim1, Dim5, o, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim6, l, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim6, i, j, k, l, m, n, l, k, o, Dim0, Dim1, Dim4, Dim5, Dim3, Dim2, Dim6, i, j, m, n, l, k, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(l, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim6, l, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim6, k, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim6, i, j, k, l, m, n, k, l, o, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, i, j, m, n, k, l, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(k, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim6, k, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim6, Dim2, l, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim2, i, j, k, l, m, n, l, o, k, Dim0, Dim1, Dim4, Dim5, Dim3, Dim6, Dim2, i, j, m, n, l, o, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(l, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim6, Dim2, l, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim6, Dim3, k, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim3, i, j, k, l, m, n, k, o, l, Dim0, Dim1, Dim4, Dim5, Dim2, Dim6, Dim3, i, j, m, n, k, o, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(k, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim6, Dim3, k, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim3, Dim2, o, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, i, j, k, l, m, n, o, l, k, Dim0, Dim1, Dim4, Dim5, Dim6, Dim3, Dim2, i, j, m, n, o, l, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(o, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim3, Dim2, o, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim2, Dim3, o, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, i, j, k, l, m, n, o, k, l, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, i, j, m, n, o, k, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(o, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim2, Dim3, o, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim6, m, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim6, i, j, k, l, m, n, m, k, o, Dim0, Dim1, Dim3, Dim5, Dim4, Dim2, Dim6, i, j, l, n, m, k, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(m, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim6, m, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim6, k, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim6, i, j, k, l, m, n, k, m, o, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, i, j, l, n, k, m, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(k, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim6, k, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim6, Dim2, m, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim2, i, j, k, l, m, n, m, o, k, Dim0, Dim1, Dim3, Dim5, Dim4, Dim6, Dim2, i, j, l, n, m, o, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(m, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim6, Dim2, m, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim6, Dim4, k, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim4, i, j, k, l, m, n, k, o, m, Dim0, Dim1, Dim3, Dim5, Dim2, Dim6, Dim4, i, j, l, n, k, o, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(k, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim6, Dim4, k, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim4, Dim2, o, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, i, j, k, l, m, n, o, m, k, Dim0, Dim1, Dim3, Dim5, Dim6, Dim4, Dim2, i, j, l, n, o, m, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(o, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim4, Dim2, o, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim2, Dim4, o, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, i, j, k, l, m, n, o, k, m, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, i, j, l, n, o, k, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(o, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim2, Dim4, o, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim6, n, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim6, i, j, k, l, m, n, n, k, o, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, i, j, l, m, n, k, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, i, j, l, m, o>(TensorExpr(a, b));
}

//B(n, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim6, n, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim6, k, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim6, i, j, k, l, m, n, k, n, o, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, i, j, l, m, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, i, j, l, m, o>(TensorExpr(a, b));
}

//B(k, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim6, k, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim2, n, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim2, i, j, k, l, m, n, n, o, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, i, j, l, m, n, o, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, i, j, l, m, o>(TensorExpr(a, b));
}

//B(n, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim2, n, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim6, Dim5, k, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim5, i, j, k, l, m, n, k, o, n, Dim0, Dim1, Dim3, Dim4, Dim2, Dim6, Dim5, i, j, l, m, k, o, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, i, j, l, m, o>(TensorExpr(a, b));
}

//B(k, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim6, Dim5, k, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim5, Dim2, o, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, i, j, k, l, m, n, o, n, k, Dim0, Dim1, Dim3, Dim4, Dim6, Dim5, Dim2, i, j, l, m, o, n, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, i, j, l, m, o>(TensorExpr(a, b));
}

//B(o, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim5, Dim2, o, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim2, Dim5, o, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, i, j, k, l, m, n, o, k, n, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, i, j, l, m, o, k, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, i, j, l, m, o>(TensorExpr(a, b));
}

//B(o, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim2, Dim5, o, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim6, m, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim6, i, j, k, l, m, n, m, l, o, Dim0, Dim1, Dim2, Dim5, Dim4, Dim3, Dim6, i, j, k, n, m, l, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(m, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim6, m, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim6, l, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim6, i, j, k, l, m, n, l, m, o, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, i, j, k, n, l, m, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(l, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim6, l, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim6, Dim3, m, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim3, i, j, k, l, m, n, m, o, l, Dim0, Dim1, Dim2, Dim5, Dim4, Dim6, Dim3, i, j, k, n, m, o, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(m, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim6, Dim3, m, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim6, Dim4, l, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim4, i, j, k, l, m, n, l, o, m, Dim0, Dim1, Dim2, Dim5, Dim3, Dim6, Dim4, i, j, k, n, l, o, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(l, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim6, Dim4, l, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim4, Dim3, o, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, i, j, k, l, m, n, o, m, l, Dim0, Dim1, Dim2, Dim5, Dim6, Dim4, Dim3, i, j, k, n, o, m, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(o, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim4, Dim3, o, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim3, Dim4, o, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, i, j, k, l, m, n, o, l, m, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, i, j, k, n, o, l, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(o, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim3, Dim4, o, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim6, n, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim6, i, j, k, l, m, n, n, l, o, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, i, j, k, m, n, l, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, i, j, k, m, o>(TensorExpr(a, b));
}

//B(n, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim6, n, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim6, l, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim6, i, j, k, l, m, n, l, n, o, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, i, j, k, m, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, i, j, k, m, o>(TensorExpr(a, b));
}

//B(l, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim6, l, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim3, n, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim3, i, j, k, l, m, n, n, o, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, i, j, k, m, n, o, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, i, j, k, m, o>(TensorExpr(a, b));
}

//B(n, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim3, n, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim6, Dim5, l, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim5, i, j, k, l, m, n, l, o, n, Dim0, Dim1, Dim2, Dim4, Dim3, Dim6, Dim5, i, j, k, m, l, o, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, i, j, k, m, o>(TensorExpr(a, b));
}

//B(l, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim6, Dim5, l, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim5, Dim3, o, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, i, j, k, l, m, n, o, n, l, Dim0, Dim1, Dim2, Dim4, Dim6, Dim5, Dim3, i, j, k, m, o, n, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, i, j, k, m, o>(TensorExpr(a, b));
}

//B(o, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim5, Dim3, o, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim3, Dim5, o, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, i, j, k, l, m, n, o, l, n, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, i, j, k, m, o, l, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, i, j, k, m, o>(TensorExpr(a, b));
}

//B(o, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim3, Dim5, o, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim6, n, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim6, i, j, k, l, m, n, n, m, o, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, i, j, k, l, n, m, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, i, j, k, l, o>(TensorExpr(a, b));
}

//B(n, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim6, n, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim6, i, j, k, l, m, n, m, n, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, i, j, k, l, o>(TensorExpr(a, b));
}

//B(m, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim4, n, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim4, i, j, k, l, m, n, n, o, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, i, j, k, l, n, o, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, i, j, k, l, o>(TensorExpr(a, b));
}

//B(n, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim4, n, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim6, Dim5, m, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim5, i, j, k, l, m, n, m, o, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim5, i, j, k, l, m, o, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, i, j, k, l, o>(TensorExpr(a, b));
}

//B(m, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim6, Dim5, m, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim5, Dim4, o, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, i, j, k, l, m, n, o, n, m, Dim0, Dim1, Dim2, Dim3, Dim6, Dim5, Dim4, i, j, k, l, o, n, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, i, j, k, l, o>(TensorExpr(a, b));
}

//B(o, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim5, Dim4, o, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim4, Dim5, o, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, i, j, k, l, m, n, o, m, n, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, i, j, k, l, o, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, i, j, k, l, o>(TensorExpr(a, b));
}

//B(o, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim4, Dim5, o, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

} //namespace FTensor

#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, n, )*B(o, p, q, ) triple contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q, int DimA, int DimB, int DimC, int DimX, int DimY, int DimZ, char a, char b, char c, char x, char y, char z>
class Tensor6_times_Tensor3_triple
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor3_Expr<B, U, Dim6, Dim7, Dim8, o, p, q> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        for(int zz = 0; zz < DimZ; ++zz)
          {
            // Permutation is where the indices get checked.
            result += Permutation6<DimA, DimB, DimC, DimX, DimY, DimZ, a, b, c, x, y, z>().eval(
              iterA, N1, N2, N3, xx, yy, zz) * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx, yy, zz);
          }
    return result;
  }

  Tensor6_times_Tensor3_triple(
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &iter_a,
    const Tensor3_Expr<B, U, Dim6, Dim7, Dim8, o, p, q> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, n, )*B(k, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim0, i, j, k, l, m, n, k, j, i, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, l, m, n, i, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(k, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim0, i, j, k, l, m, n, j, k, i, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, l, m, n, i, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(j, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim1, i, j, k, l, m, n, k, i, j, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, l, m, n, i, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(k, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim2, i, j, k, l, m, n, j, i, k, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, l, m, n, i, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(j, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim1, i, j, k, l, m, n, i, k, j, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, l, m, n, i, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(i, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, i, j, k, l, m, n, i, j, k, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, l, m, n, i, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(i, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim0, l, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim0, i, j, k, l, m, n, l, j, i, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, k, m, n, i, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(l, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim0, l, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim0, j, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim0, i, j, k, l, m, n, j, l, i, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, k, m, n, i, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(j, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim0, j, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim1, l, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim1, i, j, k, l, m, n, l, i, j, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, k, m, n, i, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(l, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim1, l, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim3, j, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim3, i, j, k, l, m, n, j, i, l, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, k, m, n, i, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(j, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim3, j, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim1, i, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim1, i, j, k, l, m, n, i, l, j, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, k, m, n, i, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(i, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim1, i, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim3, i, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim3, i, j, k, l, m, n, i, j, l, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, k, m, n, i, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(i, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim3, i, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim0, m, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim0, i, j, k, l, m, n, m, j, i, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, k, l, n, i, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(m, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim0, m, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim0, j, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim0, i, j, k, l, m, n, j, m, i, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, k, l, n, i, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(j, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim0, j, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim1, m, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim1, i, j, k, l, m, n, m, i, j, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, k, l, n, i, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(m, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim1, m, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim4, j, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim4, i, j, k, l, m, n, j, i, m, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, k, l, n, i, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(j, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim4, j, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim1, i, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim1, i, j, k, l, m, n, i, m, j, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, k, l, n, i, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(i, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim1, i, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim4, i, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim4, i, j, k, l, m, n, i, j, m, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, k, l, n, i, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(i, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim4, i, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim0, n, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim0, i, j, k, l, m, n, n, j, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, k, l, m>(TensorExpr(a, b));
}

//B(n, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim0, n, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim0, j, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim0, i, j, k, l, m, n, j, n, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, k, l, m>(TensorExpr(a, b));
}

//B(j, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim0, j, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim1, n, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim1, i, j, k, l, m, n, n, i, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, k, l, m>(TensorExpr(a, b));
}

//B(n, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim1, n, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim5, j, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim5, i, j, k, l, m, n, j, i, n, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, k, l, m>(TensorExpr(a, b));
}

//B(j, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim5, j, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim1, i, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim1, i, j, k, l, m, n, i, n, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, k, l, m>(TensorExpr(a, b));
}

//B(i, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim1, i, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim5, i, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim5, i, j, k, l, m, n, i, j, n, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, k, l, m>(TensorExpr(a, b));
}

//B(i, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim5, i, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim0, l, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim0, i, j, k, l, m, n, l, k, i, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, j, m, n, i, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(l, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim0, l, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim0, k, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim0, i, j, k, l, m, n, k, l, i, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, j, m, n, i, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(k, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim0, k, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim2, l, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim2, i, j, k, l, m, n, l, i, k, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, j, m, n, i, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(l, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim2, l, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim3, k, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim3, i, j, k, l, m, n, k, i, l, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, j, m, n, i, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(k, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim3, k, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim2, i, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim2, i, j, k, l, m, n, i, l, k, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, j, m, n, i, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(i, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim2, i, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim3, i, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim3, i, j, k, l, m, n, i, k, l, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, j, m, n, i, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(i, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim3, i, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim0, m, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim0, i, j, k, l, m, n, m, k, i, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, j, l, n, i, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(m, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim0, m, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim0, k, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim0, i, j, k, l, m, n, k, m, i, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, j, l, n, i, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(k, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim0, k, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim2, m, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim2, i, j, k, l, m, n, m, i, k, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, j, l, n, i, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(m, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim2, m, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim4, k, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim4, i, j, k, l, m, n, k, i, m, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, j, l, n, i, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(k, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim4, k, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim2, i, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim2, i, j, k, l, m, n, i, m, k, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, j, l, n, i, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(i, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim2, i, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim4, i, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim4, i, j, k, l, m, n, i, k, m, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, j, l, n, i, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(i, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim4, i, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim0, n, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim0, i, j, k, l, m, n, n, k, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, j, l, m>(TensorExpr(a, b));
}

//B(n, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim0, n, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim0, k, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim0, i, j, k, l, m, n, k, n, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, j, l, m>(TensorExpr(a, b));
}

//B(k, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim0, k, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim2, n, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim2, i, j, k, l, m, n, n, i, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, j, l, m>(TensorExpr(a, b));
}

//B(n, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim2, n, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim5, k, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim5, i, j, k, l, m, n, k, i, n, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, j, l, m>(TensorExpr(a, b));
}

//B(k, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim5, k, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim2, i, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim2, i, j, k, l, m, n, i, n, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, j, l, m>(TensorExpr(a, b));
}

//B(i, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim2, i, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim5, i, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim5, i, j, k, l, m, n, i, k, n, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, j, l, m>(TensorExpr(a, b));
}

//B(i, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim5, i, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim0, m, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim0, i, j, k, l, m, n, m, l, i, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, j, k, n, i, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(m, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim0, m, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim0, l, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim0, i, j, k, l, m, n, l, m, i, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, j, k, n, i, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(l, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim0, l, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim3, m, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim3, i, j, k, l, m, n, m, i, l, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, j, k, n, i, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(m, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim3, m, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim4, l, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim4, i, j, k, l, m, n, l, i, m, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, j, k, n, i, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(l, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim4, l, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim3, i, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim3, i, j, k, l, m, n, i, m, l, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, j, k, n, i, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(i, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim3, i, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim4, i, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim4, i, j, k, l, m, n, i, l, m, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, j, k, n, i, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(i, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim4, i, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim0, n, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim0, i, j, k, l, m, n, n, l, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, j, k, m>(TensorExpr(a, b));
}

//B(n, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim0, n, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim0, l, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim0, i, j, k, l, m, n, l, n, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, j, k, m>(TensorExpr(a, b));
}

//B(l, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim0, l, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim3, n, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim3, i, j, k, l, m, n, n, i, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, j, k, m>(TensorExpr(a, b));
}

//B(n, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim3, n, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim5, l, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim5, i, j, k, l, m, n, l, i, n, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, j, k, m>(TensorExpr(a, b));
}

//B(l, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim5, l, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim3, i, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim3, i, j, k, l, m, n, i, n, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, j, k, m>(TensorExpr(a, b));
}

//B(i, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim3, i, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim5, i, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim5, i, j, k, l, m, n, i, l, n, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, j, k, m>(TensorExpr(a, b));
}

//B(i, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim5, i, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim0, n, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim0, i, j, k, l, m, n, n, m, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, j, k, l>(TensorExpr(a, b));
}

//B(n, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim0, n, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim0, m, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim0, i, j, k, l, m, n, m, n, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, j, k, l>(TensorExpr(a, b));
}

//B(m, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim0, m, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim4, n, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim4, i, j, k, l, m, n, n, i, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, j, k, l>(TensorExpr(a, b));
}

//B(n, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim4, n, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim5, m, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim5, i, j, k, l, m, n, m, i, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, j, k, l>(TensorExpr(a, b));
}

//B(m, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim5, m, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim4, i, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim4, i, j, k, l, m, n, i, n, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, j, k, l>(TensorExpr(a, b));
}

//B(i, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim4, i, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim5, i, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim5, i, j, k, l, m, n, i, m, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, j, k, l>(TensorExpr(a, b));
}

//B(i, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim5, i, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim1, l, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim1, i, j, k, l, m, n, l, k, j, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, i, m, n, j, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(l, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim1, l, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim1, k, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim1, i, j, k, l, m, n, k, l, j, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, i, m, n, j, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(k, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim1, k, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim2, l, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim2, i, j, k, l, m, n, l, j, k, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, i, m, n, j, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(l, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim2, l, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim3, k, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim3, i, j, k, l, m, n, k, j, l, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, i, m, n, j, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(k, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim3, k, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim2, j, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim2, i, j, k, l, m, n, j, l, k, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, i, m, n, j, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(j, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim2, j, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim3, j, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim3, i, j, k, l, m, n, j, k, l, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, i, m, n, j, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(j, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim3, j, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim1, m, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim1, i, j, k, l, m, n, m, k, j, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, i, l, n, j, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(m, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim1, m, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim1, k, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim1, i, j, k, l, m, n, k, m, j, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, i, l, n, j, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(k, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim1, k, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim2, m, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim2, i, j, k, l, m, n, m, j, k, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, i, l, n, j, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(m, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim2, m, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim4, k, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim4, i, j, k, l, m, n, k, j, m, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, i, l, n, j, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(k, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim4, k, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim2, j, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim2, i, j, k, l, m, n, j, m, k, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, i, l, n, j, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(j, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim2, j, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim4, j, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim4, i, j, k, l, m, n, j, k, m, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, i, l, n, j, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(j, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim4, j, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim1, n, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim1, i, j, k, l, m, n, n, k, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, i, l, m>(TensorExpr(a, b));
}

//B(n, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim1, n, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim1, k, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim1, i, j, k, l, m, n, k, n, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, i, l, m>(TensorExpr(a, b));
}

//B(k, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim1, k, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim2, n, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim2, i, j, k, l, m, n, n, j, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, i, l, m>(TensorExpr(a, b));
}

//B(n, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim2, n, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim5, k, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim5, i, j, k, l, m, n, k, j, n, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, i, l, m>(TensorExpr(a, b));
}

//B(k, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim5, k, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim2, j, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim2, i, j, k, l, m, n, j, n, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, i, l, m>(TensorExpr(a, b));
}

//B(j, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim2, j, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim5, j, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim5, i, j, k, l, m, n, j, k, n, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, i, l, m>(TensorExpr(a, b));
}

//B(j, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim5, j, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim1, m, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim1, i, j, k, l, m, n, m, l, j, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, i, k, n, j, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(m, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim1, m, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim1, l, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim1, i, j, k, l, m, n, l, m, j, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, i, k, n, j, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(l, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim1, l, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim3, m, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim3, i, j, k, l, m, n, m, j, l, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, i, k, n, j, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(m, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim3, m, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim4, l, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim4, i, j, k, l, m, n, l, j, m, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, i, k, n, j, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(l, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim4, l, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim3, j, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim3, i, j, k, l, m, n, j, m, l, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, i, k, n, j, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(j, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim3, j, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim4, j, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim4, i, j, k, l, m, n, j, l, m, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, i, k, n, j, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(j, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim4, j, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim1, n, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim1, i, j, k, l, m, n, n, l, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, i, k, m>(TensorExpr(a, b));
}

//B(n, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim1, n, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim1, l, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim1, i, j, k, l, m, n, l, n, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, i, k, m>(TensorExpr(a, b));
}

//B(l, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim1, l, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim3, n, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim3, i, j, k, l, m, n, n, j, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, i, k, m>(TensorExpr(a, b));
}

//B(n, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim3, n, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim5, l, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim5, i, j, k, l, m, n, l, j, n, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, i, k, m>(TensorExpr(a, b));
}

//B(l, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim5, l, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim3, j, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim3, i, j, k, l, m, n, j, n, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, i, k, m>(TensorExpr(a, b));
}

//B(j, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim3, j, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim5, j, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim5, i, j, k, l, m, n, j, l, n, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, i, k, m>(TensorExpr(a, b));
}

//B(j, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim5, j, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim1, n, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim1, i, j, k, l, m, n, n, m, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, i, k, l>(TensorExpr(a, b));
}

//B(n, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim1, n, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim1, m, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim1, i, j, k, l, m, n, m, n, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, i, k, l>(TensorExpr(a, b));
}

//B(m, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim1, m, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim4, n, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim4, i, j, k, l, m, n, n, j, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, i, k, l>(TensorExpr(a, b));
}

//B(n, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim4, n, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim5, m, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim5, i, j, k, l, m, n, m, j, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, i, k, l>(TensorExpr(a, b));
}

//B(m, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim5, m, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim4, j, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim4, i, j, k, l, m, n, j, n, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, i, k, l>(TensorExpr(a, b));
}

//B(j, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim4, j, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim5, j, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim5, i, j, k, l, m, n, j, m, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, i, k, l>(TensorExpr(a, b));
}

//B(j, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim5, j, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim2, m, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim2, i, j, k, l, m, n, m, l, k, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, i, j, n, k, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(m, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim2, m, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim2, l, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim2, i, j, k, l, m, n, l, m, k, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, i, j, n, k, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(l, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim2, l, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim3, m, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim3, i, j, k, l, m, n, m, k, l, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, i, j, n, k, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(m, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim3, m, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim4, l, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim4, i, j, k, l, m, n, l, k, m, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, i, j, n, k, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(l, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim4, l, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim3, k, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim3, i, j, k, l, m, n, k, m, l, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, i, j, n, k, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(k, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim3, k, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim4, i, j, k, l, m, n, k, l, m, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, i, j, n, k, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(k, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim2, n, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim2, i, j, k, l, m, n, n, l, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, i, j, m>(TensorExpr(a, b));
}

//B(n, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim2, n, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim2, l, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim2, i, j, k, l, m, n, l, n, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, i, j, m>(TensorExpr(a, b));
}

//B(l, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim2, l, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim3, n, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim3, i, j, k, l, m, n, n, k, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, i, j, m>(TensorExpr(a, b));
}

//B(n, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim3, n, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim5, l, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim5, i, j, k, l, m, n, l, k, n, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, i, j, m>(TensorExpr(a, b));
}

//B(l, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim5, l, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim3, k, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim3, i, j, k, l, m, n, k, n, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, i, j, m>(TensorExpr(a, b));
}

//B(k, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim3, k, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim5, k, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim5, i, j, k, l, m, n, k, l, n, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, i, j, m>(TensorExpr(a, b));
}

//B(k, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim5, k, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim2, n, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim2, i, j, k, l, m, n, n, m, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, i, j, l>(TensorExpr(a, b));
}

//B(n, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim2, n, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim2, m, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim2, i, j, k, l, m, n, m, n, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, i, j, l>(TensorExpr(a, b));
}

//B(m, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim2, m, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim4, n, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim4, i, j, k, l, m, n, n, k, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, i, j, l>(TensorExpr(a, b));
}

//B(n, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim4, n, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim5, m, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim5, i, j, k, l, m, n, m, k, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, i, j, l>(TensorExpr(a, b));
}

//B(m, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim5, m, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim4, k, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim4, i, j, k, l, m, n, k, n, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, i, j, l>(TensorExpr(a, b));
}

//B(k, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim4, k, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim5, k, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim5, i, j, k, l, m, n, k, m, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, i, j, l>(TensorExpr(a, b));
}

//B(k, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim5, k, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim3, n, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim3, i, j, k, l, m, n, n, m, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, i, j, k>(TensorExpr(a, b));
}

//B(n, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim3, n, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim3, m, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim3, i, j, k, l, m, n, m, n, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, i, j, k>(TensorExpr(a, b));
}

//B(m, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim3, m, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim4, n, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim4, i, j, k, l, m, n, n, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, i, j, k>(TensorExpr(a, b));
}

//B(n, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim4, n, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim5, m, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim5, i, j, k, l, m, n, m, l, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, i, j, k>(TensorExpr(a, b));
}

//B(m, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim5, m, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim4, l, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim4, i, j, k, l, m, n, l, n, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, i, j, k>(TensorExpr(a, b));
}

//B(l, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim4, l, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim5, i, j, k, l, m, n, l, m, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, i, j, k>(TensorExpr(a, b));
}

//B(l, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

} //namespace FTensor

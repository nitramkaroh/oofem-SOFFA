#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, n, )*B(o, p, q, r, ) triple contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, int Dim9, char i, char j, char k, char l, char m, char n, char o, char p, char q, char r, int DimA, int DimB, int DimC, int DimX, int DimY, int DimZ, int DimW, char a, char b, char c, char x, char y, char z, char w>
class Tensor6_times_Tensor4_triple
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim9, o, p, q, r> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        for(int zz = 0; zz < DimZ; ++zz)
          {
            // Permutation is where the indices get checked.
            result += Permutation6<DimA, DimB, DimC, DimX, DimY, DimZ, a, b, c, x, y, z>().eval(
              iterA, N1, N2, N3, xx, yy, zz) * Permutation4<DimX, DimY, DimZ, DimW, x, y, z, w>().eval(iterB, xx, yy, zz, N4);
          }
    return result;
  }

  Tensor6_times_Tensor4_triple(
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &iter_a,
    const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim9, o, p, q, r> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, n, )*B(k, j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim6, k, j, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim0, Dim6, i, j, k, l, m, n, k, j, i, o, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(k, j, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim6, k, j, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim6, j, k, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim0, Dim6, i, j, k, l, m, n, j, k, i, o, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(j, k, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim6, j, k, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim6, k, i, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim1, Dim6, i, j, k, l, m, n, k, i, j, o, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(k, i, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim6, k, i, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim6, j, i, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim2, Dim6, i, j, k, l, m, n, j, i, k, o, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(j, i, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim6, j, i, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim6, i, k, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim1, Dim6, i, j, k, l, m, n, i, k, j, o, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(i, k, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim6, i, k, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim6, i, j, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, i, j, k, l, m, n, i, j, k, o, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(i, j, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim6, i, j, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim6, Dim0, k, j, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, Dim0, i, j, k, l, m, n, k, j, o, i, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(k, j, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim6, Dim0, k, j, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim6, Dim0, j, k, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim0, i, j, k, l, m, n, j, k, o, i, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(j, k, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim6, Dim0, j, k, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim6, Dim1, k, i, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, Dim1, i, j, k, l, m, n, k, i, o, j, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(k, i, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim6, Dim1, k, i, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim6, Dim2, j, i, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, Dim2, i, j, k, l, m, n, j, i, o, k, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(j, i, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim6, Dim2, j, i, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim6, Dim1, i, k, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim1, i, j, k, l, m, n, i, k, o, j, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(i, k, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim6, Dim1, i, k, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim6, Dim2, i, j, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim2, i, j, k, l, m, n, i, j, o, k, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(i, j, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim6, Dim2, i, j, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim1, Dim0, k, o, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, Dim0, i, j, k, l, m, n, k, o, j, i, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(k, o, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim1, Dim0, k, o, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim2, Dim0, j, o, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, Dim0, i, j, k, l, m, n, j, o, k, i, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(j, o, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim2, Dim0, j, o, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim0, Dim1, k, o, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, Dim1, i, j, k, l, m, n, k, o, i, j, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(k, o, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim0, Dim1, k, o, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim0, Dim2, j, o, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, Dim2, i, j, k, l, m, n, j, o, i, k, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(j, o, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim0, Dim2, j, o, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim2, Dim1, i, o, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, Dim1, i, j, k, l, m, n, i, o, k, j, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(i, o, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim2, Dim1, i, o, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim1, Dim2, i, o, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, Dim2, i, j, k, l, m, n, i, o, j, k, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(i, o, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim1, Dim2, i, o, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim1, Dim0, o, k, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim0, i, j, k, l, m, n, o, k, j, i, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(o, k, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim1, Dim0, o, k, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim2, Dim0, o, j, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim0, i, j, k, l, m, n, o, j, k, i, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(o, j, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim2, Dim0, o, j, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim0, Dim1, o, k, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim1, i, j, k, l, m, n, o, k, i, j, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(o, k, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim0, Dim1, o, k, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim0, Dim2, o, j, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim2, i, j, k, l, m, n, o, j, i, k, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(o, j, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim0, Dim2, o, j, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim2, Dim1, o, i, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim1, i, j, k, l, m, n, o, i, k, j, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(o, i, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim2, Dim1, o, i, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim1, Dim2, o, i, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim2, i, j, k, l, m, n, o, i, j, k, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(o, i, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim1, Dim2, o, i, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim6, l, j, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim0, Dim6, i, j, k, l, m, n, l, j, i, o, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(l, j, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim6, l, j, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim6, j, l, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim0, Dim6, i, j, k, l, m, n, j, l, i, o, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(j, l, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim6, j, l, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim6, l, i, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim1, Dim6, i, j, k, l, m, n, l, i, j, o, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(l, i, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim6, l, i, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim6, j, i, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim3, Dim6, i, j, k, l, m, n, j, i, l, o, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(j, i, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim6, j, i, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim6, i, l, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim1, Dim6, i, j, k, l, m, n, i, l, j, o, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(i, l, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim6, i, l, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim6, i, j, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, i, j, k, l, m, n, i, j, l, o, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(i, j, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim6, i, j, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim6, Dim0, l, j, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim6, Dim0, i, j, k, l, m, n, l, j, o, i, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(l, j, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim6, Dim0, l, j, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim6, Dim0, j, l, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim6, Dim0, i, j, k, l, m, n, j, l, o, i, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(j, l, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim6, Dim0, j, l, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim6, Dim1, l, i, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim6, Dim1, i, j, k, l, m, n, l, i, o, j, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(l, i, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim6, Dim1, l, i, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim6, Dim3, j, i, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, Dim3, i, j, k, l, m, n, j, i, o, l, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(j, i, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim6, Dim3, j, i, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim6, Dim1, i, l, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim6, Dim1, i, j, k, l, m, n, i, l, o, j, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(i, l, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim6, Dim1, i, l, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim6, Dim3, i, j, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim3, i, j, k, l, m, n, i, j, o, l, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(i, j, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim6, Dim3, i, j, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim1, Dim0, l, o, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim1, Dim0, i, j, k, l, m, n, l, o, j, i, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(l, o, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim1, Dim0, l, o, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim3, Dim0, j, o, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim3, Dim0, i, j, k, l, m, n, j, o, l, i, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(j, o, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim3, Dim0, j, o, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim0, Dim1, l, o, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim0, Dim1, i, j, k, l, m, n, l, o, i, j, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(l, o, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim0, Dim1, l, o, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim0, Dim3, j, o, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, Dim3, i, j, k, l, m, n, j, o, i, l, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(j, o, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim0, Dim3, j, o, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim3, Dim1, i, o, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim3, Dim1, i, j, k, l, m, n, i, o, l, j, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(i, o, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim3, Dim1, i, o, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim1, Dim3, i, o, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, Dim3, i, j, k, l, m, n, i, o, j, l, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(i, o, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim1, Dim3, i, o, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim1, Dim0, o, l, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim0, i, j, k, l, m, n, o, l, j, i, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(o, l, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim1, Dim0, o, l, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim3, Dim0, o, j, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim0, i, j, k, l, m, n, o, j, l, i, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(o, j, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim3, Dim0, o, j, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim0, Dim1, o, l, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim1, i, j, k, l, m, n, o, l, i, j, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(o, l, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim0, Dim1, o, l, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim0, Dim3, o, j, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim3, i, j, k, l, m, n, o, j, i, l, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(o, j, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim0, Dim3, o, j, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim3, Dim1, o, i, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim1, i, j, k, l, m, n, o, i, l, j, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(o, i, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim3, Dim1, o, i, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim1, Dim3, o, i, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim3, i, j, k, l, m, n, o, i, j, l, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(o, i, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim1, Dim3, o, i, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim6, m, j, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim0, Dim6, i, j, k, l, m, n, m, j, i, o, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(m, j, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim6, m, j, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim6, j, m, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim0, Dim6, i, j, k, l, m, n, j, m, i, o, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(j, m, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim6, j, m, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim6, m, i, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim1, Dim6, i, j, k, l, m, n, m, i, j, o, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(m, i, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim6, m, i, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim6, j, i, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim4, Dim6, i, j, k, l, m, n, j, i, m, o, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(j, i, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim6, j, i, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim6, i, m, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim1, Dim6, i, j, k, l, m, n, i, m, j, o, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(i, m, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim6, i, m, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim6, i, j, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim4, Dim6, i, j, k, l, m, n, i, j, m, o, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(i, j, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim6, i, j, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim6, Dim0, m, j, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim6, Dim0, i, j, k, l, m, n, m, j, o, i, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(m, j, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim6, Dim0, m, j, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim6, Dim0, j, m, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim6, Dim0, i, j, k, l, m, n, j, m, o, i, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(j, m, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim6, Dim0, j, m, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim6, Dim1, m, i, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim6, Dim1, i, j, k, l, m, n, m, i, o, j, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(m, i, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim6, Dim1, m, i, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim6, Dim4, j, i, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, Dim4, i, j, k, l, m, n, j, i, o, m, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(j, i, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim6, Dim4, j, i, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim6, Dim1, i, m, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim6, Dim1, i, j, k, l, m, n, i, m, o, j, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(i, m, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim6, Dim1, i, m, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim6, Dim4, i, j, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim4, i, j, k, l, m, n, i, j, o, m, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(i, j, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim6, Dim4, i, j, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim1, Dim0, m, o, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim1, Dim0, i, j, k, l, m, n, m, o, j, i, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(m, o, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim1, Dim0, m, o, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim4, Dim0, j, o, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim4, Dim0, i, j, k, l, m, n, j, o, m, i, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(j, o, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim4, Dim0, j, o, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim0, Dim1, m, o, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim0, Dim1, i, j, k, l, m, n, m, o, i, j, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(m, o, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim0, Dim1, m, o, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim0, Dim4, j, o, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, Dim4, i, j, k, l, m, n, j, o, i, m, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(j, o, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim0, Dim4, j, o, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim4, Dim1, i, o, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim4, Dim1, i, j, k, l, m, n, i, o, m, j, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(i, o, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim4, Dim1, i, o, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim1, Dim4, i, o, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, Dim4, i, j, k, l, m, n, i, o, j, m, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(i, o, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim1, Dim4, i, o, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim1, Dim0, o, m, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim0, i, j, k, l, m, n, o, m, j, i, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(o, m, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim1, Dim0, o, m, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim4, Dim0, o, j, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim0, i, j, k, l, m, n, o, j, m, i, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(o, j, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim4, Dim0, o, j, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim0, Dim1, o, m, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim1, i, j, k, l, m, n, o, m, i, j, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(o, m, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim0, Dim1, o, m, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim0, Dim4, o, j, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim4, i, j, k, l, m, n, o, j, i, m, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(o, j, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim0, Dim4, o, j, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim4, Dim1, o, i, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim1, i, j, k, l, m, n, o, i, m, j, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(o, i, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim4, Dim1, o, i, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim1, Dim4, o, i, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim4, i, j, k, l, m, n, o, i, j, m, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(o, i, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim1, Dim4, o, i, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim6, n, j, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim0, Dim6, i, j, k, l, m, n, n, j, i, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(n, j, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim6, n, j, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim6, j, n, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim0, Dim6, i, j, k, l, m, n, j, n, i, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(j, n, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim6, j, n, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim6, n, i, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim1, Dim6, i, j, k, l, m, n, n, i, j, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(n, i, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim6, n, i, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim6, j, i, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim5, Dim6, i, j, k, l, m, n, j, i, n, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(j, i, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim6, j, i, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim6, i, n, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim1, Dim6, i, j, k, l, m, n, i, n, j, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(i, n, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim6, i, n, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim6, i, j, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim5, Dim6, i, j, k, l, m, n, i, j, n, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(i, j, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim6, i, j, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim0, n, j, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim6, Dim0, i, j, k, l, m, n, n, j, o, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(n, j, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim0, n, j, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim0, j, n, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim6, Dim0, i, j, k, l, m, n, j, n, o, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(j, n, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim0, j, n, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim1, n, i, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim6, Dim1, i, j, k, l, m, n, n, i, o, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(n, i, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim1, n, i, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim6, Dim5, j, i, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, Dim5, i, j, k, l, m, n, j, i, o, n, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(j, i, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim6, Dim5, j, i, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim1, i, n, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim6, Dim1, i, j, k, l, m, n, i, n, o, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(i, n, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim1, i, n, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim6, Dim5, i, j, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim5, i, j, k, l, m, n, i, j, o, n, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(i, j, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim6, Dim5, i, j, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim0, n, o, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim1, Dim0, i, j, k, l, m, n, n, o, j, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(n, o, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim0, n, o, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim5, Dim0, j, o, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim5, Dim0, i, j, k, l, m, n, j, o, n, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(j, o, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim5, Dim0, j, o, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim1, n, o, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim0, Dim1, i, j, k, l, m, n, n, o, i, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(n, o, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim1, n, o, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim0, Dim5, j, o, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, Dim5, i, j, k, l, m, n, j, o, i, n, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(j, o, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim0, Dim5, j, o, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim5, Dim1, i, o, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim5, Dim1, i, j, k, l, m, n, i, o, n, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(i, o, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim5, Dim1, i, o, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim1, Dim5, i, o, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, Dim5, i, j, k, l, m, n, i, o, j, n, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(i, o, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim1, Dim5, i, o, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim1, Dim0, o, n, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim0, i, j, k, l, m, n, o, n, j, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(o, n, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim1, Dim0, o, n, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim5, Dim0, o, j, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim0, i, j, k, l, m, n, o, j, n, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(o, j, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim5, Dim0, o, j, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim0, Dim1, o, n, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim1, i, j, k, l, m, n, o, n, i, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(o, n, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim0, Dim1, o, n, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim0, Dim5, o, j, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim5, i, j, k, l, m, n, o, j, i, n, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(o, j, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim0, Dim5, o, j, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim5, Dim1, o, i, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim1, i, j, k, l, m, n, o, i, n, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(o, i, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim5, Dim1, o, i, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim1, Dim5, o, i, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim5, i, j, k, l, m, n, o, i, j, n, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(o, i, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim1, Dim5, o, i, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim6, l, k, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim0, Dim6, i, j, k, l, m, n, l, k, i, o, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(l, k, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim6, l, k, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim6, k, l, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim0, Dim6, i, j, k, l, m, n, k, l, i, o, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(k, l, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim6, k, l, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim6, l, i, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim2, Dim6, i, j, k, l, m, n, l, i, k, o, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(l, i, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim6, l, i, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim6, k, i, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim3, Dim6, i, j, k, l, m, n, k, i, l, o, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(k, i, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim6, k, i, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim6, i, l, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim2, Dim6, i, j, k, l, m, n, i, l, k, o, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(i, l, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim6, i, l, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim6, i, k, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, i, j, k, l, m, n, i, k, l, o, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(i, k, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim6, i, k, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim6, Dim0, l, k, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim6, Dim0, i, j, k, l, m, n, l, k, o, i, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(l, k, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim6, Dim0, l, k, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim6, Dim0, k, l, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim6, Dim0, i, j, k, l, m, n, k, l, o, i, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(k, l, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim6, Dim0, k, l, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim6, Dim2, l, i, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim6, Dim2, i, j, k, l, m, n, l, i, o, k, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(l, i, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim6, Dim2, l, i, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim6, Dim3, k, i, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, Dim3, i, j, k, l, m, n, k, i, o, l, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(k, i, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim6, Dim3, k, i, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim6, Dim2, i, l, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim6, Dim2, i, j, k, l, m, n, i, l, o, k, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(i, l, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim6, Dim2, i, l, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim6, Dim3, i, k, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim3, i, j, k, l, m, n, i, k, o, l, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(i, k, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim6, Dim3, i, k, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim2, Dim0, l, o, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim2, Dim0, i, j, k, l, m, n, l, o, k, i, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(l, o, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim2, Dim0, l, o, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim3, Dim0, k, o, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim3, Dim0, i, j, k, l, m, n, k, o, l, i, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(k, o, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim3, Dim0, k, o, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim0, Dim2, l, o, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim0, Dim2, i, j, k, l, m, n, l, o, i, k, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(l, o, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim0, Dim2, l, o, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim0, Dim3, k, o, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, Dim3, i, j, k, l, m, n, k, o, i, l, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(k, o, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim0, Dim3, k, o, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim3, Dim2, i, o, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim3, Dim2, i, j, k, l, m, n, i, o, l, k, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(i, o, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim3, Dim2, i, o, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim2, Dim3, i, o, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, Dim3, i, j, k, l, m, n, i, o, k, l, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(i, o, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim2, Dim3, i, o, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim2, Dim0, o, l, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim0, i, j, k, l, m, n, o, l, k, i, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(o, l, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim2, Dim0, o, l, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim3, Dim0, o, k, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim0, i, j, k, l, m, n, o, k, l, i, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(o, k, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim3, Dim0, o, k, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim0, Dim2, o, l, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim2, i, j, k, l, m, n, o, l, i, k, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(o, l, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim0, Dim2, o, l, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim0, Dim3, o, k, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim3, i, j, k, l, m, n, o, k, i, l, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(o, k, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim0, Dim3, o, k, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim3, Dim2, o, i, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim2, i, j, k, l, m, n, o, i, l, k, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(o, i, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim3, Dim2, o, i, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim2, Dim3, o, i, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim3, i, j, k, l, m, n, o, i, k, l, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(o, i, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim2, Dim3, o, i, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim6, m, k, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim0, Dim6, i, j, k, l, m, n, m, k, i, o, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(m, k, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim6, m, k, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim6, k, m, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim0, Dim6, i, j, k, l, m, n, k, m, i, o, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(k, m, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim6, k, m, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim6, m, i, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim2, Dim6, i, j, k, l, m, n, m, i, k, o, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(m, i, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim6, m, i, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim6, k, i, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim4, Dim6, i, j, k, l, m, n, k, i, m, o, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(k, i, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim6, k, i, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim6, i, m, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim2, Dim6, i, j, k, l, m, n, i, m, k, o, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(i, m, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim6, i, m, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim6, i, k, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim4, Dim6, i, j, k, l, m, n, i, k, m, o, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(i, k, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim6, i, k, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim6, Dim0, m, k, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim6, Dim0, i, j, k, l, m, n, m, k, o, i, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(m, k, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim6, Dim0, m, k, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim6, Dim0, k, m, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim6, Dim0, i, j, k, l, m, n, k, m, o, i, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(k, m, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim6, Dim0, k, m, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim6, Dim2, m, i, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim6, Dim2, i, j, k, l, m, n, m, i, o, k, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(m, i, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim6, Dim2, m, i, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim6, Dim4, k, i, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, Dim4, i, j, k, l, m, n, k, i, o, m, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(k, i, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim6, Dim4, k, i, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim6, Dim2, i, m, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim6, Dim2, i, j, k, l, m, n, i, m, o, k, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(i, m, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim6, Dim2, i, m, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim6, Dim4, i, k, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim4, i, j, k, l, m, n, i, k, o, m, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(i, k, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim6, Dim4, i, k, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim2, Dim0, m, o, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim2, Dim0, i, j, k, l, m, n, m, o, k, i, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(m, o, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim2, Dim0, m, o, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim4, Dim0, k, o, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim4, Dim0, i, j, k, l, m, n, k, o, m, i, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(k, o, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim4, Dim0, k, o, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim0, Dim2, m, o, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim0, Dim2, i, j, k, l, m, n, m, o, i, k, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(m, o, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim0, Dim2, m, o, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim0, Dim4, k, o, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, Dim4, i, j, k, l, m, n, k, o, i, m, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(k, o, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim0, Dim4, k, o, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim4, Dim2, i, o, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim4, Dim2, i, j, k, l, m, n, i, o, m, k, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(i, o, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim4, Dim2, i, o, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim2, Dim4, i, o, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, Dim4, i, j, k, l, m, n, i, o, k, m, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(i, o, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim2, Dim4, i, o, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim2, Dim0, o, m, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim0, i, j, k, l, m, n, o, m, k, i, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(o, m, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim2, Dim0, o, m, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim4, Dim0, o, k, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim0, i, j, k, l, m, n, o, k, m, i, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(o, k, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim4, Dim0, o, k, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim0, Dim2, o, m, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim2, i, j, k, l, m, n, o, m, i, k, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(o, m, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim0, Dim2, o, m, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim0, Dim4, o, k, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim4, i, j, k, l, m, n, o, k, i, m, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(o, k, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim0, Dim4, o, k, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim4, Dim2, o, i, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim2, i, j, k, l, m, n, o, i, m, k, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(o, i, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim4, Dim2, o, i, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim2, Dim4, o, i, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim4, i, j, k, l, m, n, o, i, k, m, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(o, i, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim2, Dim4, o, i, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim6, n, k, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim0, Dim6, i, j, k, l, m, n, n, k, i, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(n, k, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim6, n, k, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim6, k, n, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim0, Dim6, i, j, k, l, m, n, k, n, i, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(k, n, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim6, k, n, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim6, n, i, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim2, Dim6, i, j, k, l, m, n, n, i, k, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(n, i, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim6, n, i, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim6, k, i, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim5, Dim6, i, j, k, l, m, n, k, i, n, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(k, i, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim6, k, i, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim6, i, n, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim2, Dim6, i, j, k, l, m, n, i, n, k, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(i, n, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim6, i, n, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim6, i, k, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim5, Dim6, i, j, k, l, m, n, i, k, n, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(i, k, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim6, i, k, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim0, n, k, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim6, Dim0, i, j, k, l, m, n, n, k, o, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(n, k, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim0, n, k, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim0, k, n, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim6, Dim0, i, j, k, l, m, n, k, n, o, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(k, n, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim0, k, n, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim2, n, i, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim6, Dim2, i, j, k, l, m, n, n, i, o, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(n, i, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim2, n, i, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim6, Dim5, k, i, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, Dim5, i, j, k, l, m, n, k, i, o, n, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(k, i, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim6, Dim5, k, i, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim2, i, n, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim6, Dim2, i, j, k, l, m, n, i, n, o, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(i, n, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim2, i, n, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim6, Dim5, i, k, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim5, i, j, k, l, m, n, i, k, o, n, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(i, k, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim6, Dim5, i, k, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim0, n, o, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim2, Dim0, i, j, k, l, m, n, n, o, k, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(n, o, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim0, n, o, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim5, Dim0, k, o, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim5, Dim0, i, j, k, l, m, n, k, o, n, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(k, o, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim5, Dim0, k, o, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim2, n, o, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim0, Dim2, i, j, k, l, m, n, n, o, i, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(n, o, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim2, n, o, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim0, Dim5, k, o, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, Dim5, i, j, k, l, m, n, k, o, i, n, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(k, o, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim0, Dim5, k, o, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim5, Dim2, i, o, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim5, Dim2, i, j, k, l, m, n, i, o, n, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(i, o, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim5, Dim2, i, o, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim2, Dim5, i, o, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, Dim5, i, j, k, l, m, n, i, o, k, n, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(i, o, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim2, Dim5, i, o, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim2, Dim0, o, n, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim0, i, j, k, l, m, n, o, n, k, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(o, n, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim2, Dim0, o, n, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim5, Dim0, o, k, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim0, i, j, k, l, m, n, o, k, n, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(o, k, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim5, Dim0, o, k, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim0, Dim2, o, n, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim2, i, j, k, l, m, n, o, n, i, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(o, n, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim0, Dim2, o, n, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim0, Dim5, o, k, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim5, i, j, k, l, m, n, o, k, i, n, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(o, k, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim0, Dim5, o, k, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim5, Dim2, o, i, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim2, i, j, k, l, m, n, o, i, n, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(o, i, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim5, Dim2, o, i, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim2, Dim5, o, i, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim5, i, j, k, l, m, n, o, i, k, n, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(o, i, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim2, Dim5, o, i, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim6, m, l, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim0, Dim6, i, j, k, l, m, n, m, l, i, o, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(m, l, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim6, m, l, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim6, l, m, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim0, Dim6, i, j, k, l, m, n, l, m, i, o, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(l, m, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim6, l, m, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim6, m, i, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim3, Dim6, i, j, k, l, m, n, m, i, l, o, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(m, i, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim6, m, i, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim6, l, i, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim4, Dim6, i, j, k, l, m, n, l, i, m, o, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(l, i, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim6, l, i, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim6, i, m, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim3, Dim6, i, j, k, l, m, n, i, m, l, o, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(i, m, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim6, i, m, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim6, i, l, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim4, Dim6, i, j, k, l, m, n, i, l, m, o, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(i, l, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim6, i, l, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim6, Dim0, m, l, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim6, Dim0, i, j, k, l, m, n, m, l, o, i, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(m, l, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim6, Dim0, m, l, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim6, Dim0, l, m, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim6, Dim0, i, j, k, l, m, n, l, m, o, i, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(l, m, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim6, Dim0, l, m, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim6, Dim3, m, i, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim6, Dim3, i, j, k, l, m, n, m, i, o, l, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(m, i, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim6, Dim3, m, i, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim6, Dim4, l, i, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim6, Dim4, i, j, k, l, m, n, l, i, o, m, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(l, i, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim6, Dim4, l, i, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim6, Dim3, i, m, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim6, Dim3, i, j, k, l, m, n, i, m, o, l, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(i, m, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim6, Dim3, i, m, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim6, Dim4, i, l, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim6, Dim4, i, j, k, l, m, n, i, l, o, m, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(i, l, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim6, Dim4, i, l, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim3, Dim0, m, o, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim3, Dim0, i, j, k, l, m, n, m, o, l, i, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(m, o, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim3, Dim0, m, o, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim4, Dim0, l, o, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim4, Dim0, i, j, k, l, m, n, l, o, m, i, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(l, o, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim4, Dim0, l, o, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim0, Dim3, m, o, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim0, Dim3, i, j, k, l, m, n, m, o, i, l, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(m, o, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim0, Dim3, m, o, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim0, Dim4, l, o, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim0, Dim4, i, j, k, l, m, n, l, o, i, m, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(l, o, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim0, Dim4, l, o, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim4, Dim3, i, o, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim4, Dim3, i, j, k, l, m, n, i, o, m, l, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(i, o, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim4, Dim3, i, o, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim3, Dim4, i, o, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim3, Dim4, i, j, k, l, m, n, i, o, l, m, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(i, o, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim3, Dim4, i, o, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim3, Dim0, o, m, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim0, i, j, k, l, m, n, o, m, l, i, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(o, m, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim3, Dim0, o, m, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim4, Dim0, o, l, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim0, i, j, k, l, m, n, o, l, m, i, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(o, l, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim4, Dim0, o, l, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim0, Dim3, o, m, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim3, i, j, k, l, m, n, o, m, i, l, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(o, m, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim0, Dim3, o, m, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim0, Dim4, o, l, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim4, i, j, k, l, m, n, o, l, i, m, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(o, l, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim0, Dim4, o, l, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim4, Dim3, o, i, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim3, i, j, k, l, m, n, o, i, m, l, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(o, i, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim4, Dim3, o, i, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim3, Dim4, o, i, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim4, i, j, k, l, m, n, o, i, l, m, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(o, i, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim3, Dim4, o, i, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim6, n, l, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim0, Dim6, i, j, k, l, m, n, n, l, i, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(n, l, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim6, n, l, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim6, l, n, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim0, Dim6, i, j, k, l, m, n, l, n, i, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(l, n, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim6, l, n, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim6, n, i, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim3, Dim6, i, j, k, l, m, n, n, i, l, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(n, i, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim6, n, i, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim6, l, i, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim5, Dim6, i, j, k, l, m, n, l, i, n, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(l, i, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim6, l, i, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim6, i, n, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim3, Dim6, i, j, k, l, m, n, i, n, l, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(i, n, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim6, i, n, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim6, i, l, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim5, Dim6, i, j, k, l, m, n, i, l, n, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(i, l, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim6, i, l, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim0, n, l, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim6, Dim0, i, j, k, l, m, n, n, l, o, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(n, l, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim0, n, l, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim0, l, n, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim6, Dim0, i, j, k, l, m, n, l, n, o, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(l, n, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim0, l, n, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim3, n, i, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim6, Dim3, i, j, k, l, m, n, n, i, o, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(n, i, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim3, n, i, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim6, Dim5, l, i, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim6, Dim5, i, j, k, l, m, n, l, i, o, n, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(l, i, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim6, Dim5, l, i, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim3, i, n, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim6, Dim3, i, j, k, l, m, n, i, n, o, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(i, n, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim3, i, n, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim6, Dim5, i, l, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim6, Dim5, i, j, k, l, m, n, i, l, o, n, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(i, l, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim6, Dim5, i, l, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim0, n, o, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim3, Dim0, i, j, k, l, m, n, n, o, l, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(n, o, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim0, n, o, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim5, Dim0, l, o, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim5, Dim0, i, j, k, l, m, n, l, o, n, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(l, o, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim5, Dim0, l, o, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim3, n, o, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim0, Dim3, i, j, k, l, m, n, n, o, i, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(n, o, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim3, n, o, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim0, Dim5, l, o, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim0, Dim5, i, j, k, l, m, n, l, o, i, n, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(l, o, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim0, Dim5, l, o, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim5, Dim3, i, o, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim5, Dim3, i, j, k, l, m, n, i, o, n, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(i, o, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim5, Dim3, i, o, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim3, Dim5, i, o, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim3, Dim5, i, j, k, l, m, n, i, o, l, n, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(i, o, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim3, Dim5, i, o, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim3, Dim0, o, n, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim0, i, j, k, l, m, n, o, n, l, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(o, n, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim3, Dim0, o, n, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim5, Dim0, o, l, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim0, i, j, k, l, m, n, o, l, n, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(o, l, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim5, Dim0, o, l, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim0, Dim3, o, n, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim3, i, j, k, l, m, n, o, n, i, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(o, n, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim0, Dim3, o, n, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim0, Dim5, o, l, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim5, i, j, k, l, m, n, o, l, i, n, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(o, l, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim0, Dim5, o, l, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim5, Dim3, o, i, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim3, i, j, k, l, m, n, o, i, n, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(o, i, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim5, Dim3, o, i, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim3, Dim5, o, i, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim5, i, j, k, l, m, n, o, i, l, n, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(o, i, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim3, Dim5, o, i, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim6, n, m, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim0, Dim6, i, j, k, l, m, n, n, m, i, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(n, m, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim6, n, m, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim6, m, n, i, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim0, Dim6, i, j, k, l, m, n, m, n, i, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(m, n, i, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim6, m, n, i, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim6, n, i, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim4, Dim6, i, j, k, l, m, n, n, i, m, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(n, i, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim6, n, i, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim6, m, i, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim5, Dim6, i, j, k, l, m, n, m, i, n, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(m, i, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim6, m, i, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim6, i, n, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim4, Dim6, i, j, k, l, m, n, i, n, m, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(i, n, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim6, i, n, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim6, i, m, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim5, Dim6, i, j, k, l, m, n, i, m, n, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(i, m, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim6, i, m, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim0, n, m, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim6, Dim0, i, j, k, l, m, n, n, m, o, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(n, m, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim0, n, m, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim0, m, n, o, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim6, Dim0, i, j, k, l, m, n, m, n, o, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(m, n, o, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim0, m, n, o, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim4, n, i, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim6, Dim4, i, j, k, l, m, n, n, i, o, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(n, i, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim4, n, i, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim6, Dim5, m, i, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim6, Dim5, i, j, k, l, m, n, m, i, o, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(m, i, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim6, Dim5, m, i, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim4, i, n, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim6, Dim4, i, j, k, l, m, n, i, n, o, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(i, n, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim4, i, n, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim6, Dim5, i, m, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim6, Dim5, i, j, k, l, m, n, i, m, o, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(i, m, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim6, Dim5, i, m, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim0, n, o, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim4, Dim0, i, j, k, l, m, n, n, o, m, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(n, o, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim0, n, o, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim5, Dim0, m, o, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim5, Dim0, i, j, k, l, m, n, m, o, n, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(m, o, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim5, Dim0, m, o, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim4, n, o, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim0, Dim4, i, j, k, l, m, n, n, o, i, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(n, o, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim4, n, o, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim0, Dim5, m, o, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim0, Dim5, i, j, k, l, m, n, m, o, i, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(m, o, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim0, Dim5, m, o, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim5, Dim4, i, o, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim5, Dim4, i, j, k, l, m, n, i, o, n, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(i, o, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim5, Dim4, i, o, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim4, Dim5, i, o, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim4, Dim5, i, j, k, l, m, n, i, o, m, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(i, o, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim4, Dim5, i, o, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim4, Dim0, o, n, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim0, i, j, k, l, m, n, o, n, m, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(o, n, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim4, Dim0, o, n, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim5, Dim0, o, m, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim0, i, j, k, l, m, n, o, m, n, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(o, m, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim5, Dim0, o, m, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim0, Dim4, o, n, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim4, i, j, k, l, m, n, o, n, i, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(o, n, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim0, Dim4, o, n, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim0, Dim5, o, m, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim5, i, j, k, l, m, n, o, m, i, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(o, m, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim0, Dim5, o, m, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim5, Dim4, o, i, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim4, i, j, k, l, m, n, o, i, n, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(o, i, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim5, Dim4, o, i, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim4, Dim5, o, i, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim5, i, j, k, l, m, n, o, i, m, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(o, i, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim4, Dim5, o, i, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim6, l, k, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim1, Dim6, i, j, k, l, m, n, l, k, j, o, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(l, k, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim6, l, k, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim6, k, l, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim1, Dim6, i, j, k, l, m, n, k, l, j, o, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(k, l, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim6, k, l, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim6, l, j, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim2, Dim6, i, j, k, l, m, n, l, j, k, o, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(l, j, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim6, l, j, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim6, k, j, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim3, Dim6, i, j, k, l, m, n, k, j, l, o, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(k, j, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim6, k, j, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim6, j, l, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim2, Dim6, i, j, k, l, m, n, j, l, k, o, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(j, l, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim6, j, l, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim6, j, k, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, j, k, l, m, n, j, k, l, o, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(j, k, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim6, j, k, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim6, Dim1, l, k, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim6, Dim1, i, j, k, l, m, n, l, k, o, j, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(l, k, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim6, Dim1, l, k, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim6, Dim1, k, l, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim6, Dim1, i, j, k, l, m, n, k, l, o, j, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(k, l, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim6, Dim1, k, l, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim6, Dim2, l, j, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim6, Dim2, i, j, k, l, m, n, l, j, o, k, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(l, j, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim6, Dim2, l, j, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim6, Dim3, k, j, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, Dim3, i, j, k, l, m, n, k, j, o, l, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(k, j, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim6, Dim3, k, j, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim6, Dim2, j, l, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim6, Dim2, i, j, k, l, m, n, j, l, o, k, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(j, l, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim6, Dim2, j, l, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim6, Dim3, j, k, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim3, i, j, k, l, m, n, j, k, o, l, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(j, k, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim6, Dim3, j, k, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim2, Dim1, l, o, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim2, Dim1, i, j, k, l, m, n, l, o, k, j, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(l, o, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim2, Dim1, l, o, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim3, Dim1, k, o, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim3, Dim1, i, j, k, l, m, n, k, o, l, j, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(k, o, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim3, Dim1, k, o, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim1, Dim2, l, o, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim1, Dim2, i, j, k, l, m, n, l, o, j, k, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(l, o, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim1, Dim2, l, o, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim1, Dim3, k, o, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, Dim3, i, j, k, l, m, n, k, o, j, l, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(k, o, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim1, Dim3, k, o, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim3, Dim2, j, o, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim3, Dim2, i, j, k, l, m, n, j, o, l, k, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(j, o, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim3, Dim2, j, o, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim2, Dim3, j, o, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, Dim3, i, j, k, l, m, n, j, o, k, l, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(j, o, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim2, Dim3, j, o, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim2, Dim1, o, l, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim1, i, j, k, l, m, n, o, l, k, j, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(o, l, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim2, Dim1, o, l, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim3, Dim1, o, k, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim1, i, j, k, l, m, n, o, k, l, j, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(o, k, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim3, Dim1, o, k, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim1, Dim2, o, l, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim2, i, j, k, l, m, n, o, l, j, k, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(o, l, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim1, Dim2, o, l, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim1, Dim3, o, k, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim3, i, j, k, l, m, n, o, k, j, l, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(o, k, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim1, Dim3, o, k, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim3, Dim2, o, j, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim2, i, j, k, l, m, n, o, j, l, k, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(o, j, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim3, Dim2, o, j, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim2, Dim3, o, j, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim3, i, j, k, l, m, n, o, j, k, l, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(o, j, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim2, Dim3, o, j, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim6, m, k, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim1, Dim6, i, j, k, l, m, n, m, k, j, o, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(m, k, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim6, m, k, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim6, k, m, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim1, Dim6, i, j, k, l, m, n, k, m, j, o, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(k, m, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim6, k, m, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim6, m, j, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim2, Dim6, i, j, k, l, m, n, m, j, k, o, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(m, j, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim6, m, j, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim6, k, j, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim4, Dim6, i, j, k, l, m, n, k, j, m, o, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(k, j, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim6, k, j, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim6, j, m, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim2, Dim6, i, j, k, l, m, n, j, m, k, o, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(j, m, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim6, j, m, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim6, j, k, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim4, Dim6, i, j, k, l, m, n, j, k, m, o, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(j, k, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim6, j, k, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim6, Dim1, m, k, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim6, Dim1, i, j, k, l, m, n, m, k, o, j, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(m, k, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim6, Dim1, m, k, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim6, Dim1, k, m, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim6, Dim1, i, j, k, l, m, n, k, m, o, j, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(k, m, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim6, Dim1, k, m, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim6, Dim2, m, j, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim6, Dim2, i, j, k, l, m, n, m, j, o, k, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(m, j, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim6, Dim2, m, j, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim6, Dim4, k, j, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, Dim4, i, j, k, l, m, n, k, j, o, m, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(k, j, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim6, Dim4, k, j, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim6, Dim2, j, m, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim6, Dim2, i, j, k, l, m, n, j, m, o, k, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(j, m, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim6, Dim2, j, m, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim6, Dim4, j, k, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim4, i, j, k, l, m, n, j, k, o, m, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(j, k, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim6, Dim4, j, k, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim2, Dim1, m, o, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim2, Dim1, i, j, k, l, m, n, m, o, k, j, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(m, o, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim2, Dim1, m, o, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim4, Dim1, k, o, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim4, Dim1, i, j, k, l, m, n, k, o, m, j, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(k, o, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim4, Dim1, k, o, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim1, Dim2, m, o, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim1, Dim2, i, j, k, l, m, n, m, o, j, k, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(m, o, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim1, Dim2, m, o, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim1, Dim4, k, o, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, Dim4, i, j, k, l, m, n, k, o, j, m, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(k, o, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim1, Dim4, k, o, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim4, Dim2, j, o, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim4, Dim2, i, j, k, l, m, n, j, o, m, k, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(j, o, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim4, Dim2, j, o, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim2, Dim4, j, o, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, Dim4, i, j, k, l, m, n, j, o, k, m, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(j, o, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim2, Dim4, j, o, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim2, Dim1, o, m, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim1, i, j, k, l, m, n, o, m, k, j, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(o, m, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim2, Dim1, o, m, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim4, Dim1, o, k, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim1, i, j, k, l, m, n, o, k, m, j, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(o, k, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim4, Dim1, o, k, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim1, Dim2, o, m, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim2, i, j, k, l, m, n, o, m, j, k, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(o, m, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim1, Dim2, o, m, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim1, Dim4, o, k, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim4, i, j, k, l, m, n, o, k, j, m, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(o, k, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim1, Dim4, o, k, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim4, Dim2, o, j, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim2, i, j, k, l, m, n, o, j, m, k, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(o, j, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim4, Dim2, o, j, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim2, Dim4, o, j, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim4, i, j, k, l, m, n, o, j, k, m, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(o, j, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim2, Dim4, o, j, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim6, n, k, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim1, Dim6, i, j, k, l, m, n, n, k, j, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(n, k, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim6, n, k, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim6, k, n, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim1, Dim6, i, j, k, l, m, n, k, n, j, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(k, n, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim6, k, n, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim6, n, j, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim2, Dim6, i, j, k, l, m, n, n, j, k, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(n, j, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim6, n, j, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim6, k, j, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim5, Dim6, i, j, k, l, m, n, k, j, n, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(k, j, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim6, k, j, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim6, j, n, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim2, Dim6, i, j, k, l, m, n, j, n, k, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(j, n, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim6, j, n, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim6, j, k, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim5, Dim6, i, j, k, l, m, n, j, k, n, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(j, k, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim6, j, k, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim1, n, k, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim6, Dim1, i, j, k, l, m, n, n, k, o, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(n, k, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim1, n, k, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim1, k, n, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim6, Dim1, i, j, k, l, m, n, k, n, o, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(k, n, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim1, k, n, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim2, n, j, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim6, Dim2, i, j, k, l, m, n, n, j, o, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(n, j, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim2, n, j, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim6, Dim5, k, j, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, Dim5, i, j, k, l, m, n, k, j, o, n, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(k, j, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim6, Dim5, k, j, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim2, j, n, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim6, Dim2, i, j, k, l, m, n, j, n, o, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(j, n, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim2, j, n, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim6, Dim5, j, k, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim5, i, j, k, l, m, n, j, k, o, n, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(j, k, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim6, Dim5, j, k, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim1, n, o, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim2, Dim1, i, j, k, l, m, n, n, o, k, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(n, o, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim1, n, o, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim5, Dim1, k, o, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim5, Dim1, i, j, k, l, m, n, k, o, n, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(k, o, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim5, Dim1, k, o, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim2, n, o, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim1, Dim2, i, j, k, l, m, n, n, o, j, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(n, o, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim2, n, o, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim1, Dim5, k, o, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, Dim5, i, j, k, l, m, n, k, o, j, n, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(k, o, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim1, Dim5, k, o, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim5, Dim2, j, o, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim5, Dim2, i, j, k, l, m, n, j, o, n, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(j, o, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim5, Dim2, j, o, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim2, Dim5, j, o, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, Dim5, i, j, k, l, m, n, j, o, k, n, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(j, o, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim2, Dim5, j, o, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim2, Dim1, o, n, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim1, i, j, k, l, m, n, o, n, k, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(o, n, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim2, Dim1, o, n, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim5, Dim1, o, k, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim1, i, j, k, l, m, n, o, k, n, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(o, k, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim5, Dim1, o, k, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim1, Dim2, o, n, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim2, i, j, k, l, m, n, o, n, j, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(o, n, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim1, Dim2, o, n, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim1, Dim5, o, k, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim5, i, j, k, l, m, n, o, k, j, n, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(o, k, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim1, Dim5, o, k, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim5, Dim2, o, j, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim2, i, j, k, l, m, n, o, j, n, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(o, j, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim5, Dim2, o, j, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim2, Dim5, o, j, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim5, i, j, k, l, m, n, o, j, k, n, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(o, j, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim2, Dim5, o, j, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim6, m, l, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim1, Dim6, i, j, k, l, m, n, m, l, j, o, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(m, l, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim6, m, l, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim6, l, m, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim1, Dim6, i, j, k, l, m, n, l, m, j, o, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(l, m, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim6, l, m, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim6, m, j, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim3, Dim6, i, j, k, l, m, n, m, j, l, o, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(m, j, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim6, m, j, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim6, l, j, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim4, Dim6, i, j, k, l, m, n, l, j, m, o, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(l, j, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim6, l, j, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim6, j, m, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim3, Dim6, i, j, k, l, m, n, j, m, l, o, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(j, m, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim6, j, m, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim6, j, l, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim4, Dim6, i, j, k, l, m, n, j, l, m, o, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(j, l, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim6, j, l, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim6, Dim1, m, l, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim6, Dim1, i, j, k, l, m, n, m, l, o, j, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(m, l, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim6, Dim1, m, l, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim6, Dim1, l, m, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim6, Dim1, i, j, k, l, m, n, l, m, o, j, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(l, m, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim6, Dim1, l, m, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim6, Dim3, m, j, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim6, Dim3, i, j, k, l, m, n, m, j, o, l, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(m, j, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim6, Dim3, m, j, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim6, Dim4, l, j, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim6, Dim4, i, j, k, l, m, n, l, j, o, m, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(l, j, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim6, Dim4, l, j, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim6, Dim3, j, m, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim6, Dim3, i, j, k, l, m, n, j, m, o, l, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(j, m, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim6, Dim3, j, m, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim6, Dim4, j, l, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim6, Dim4, i, j, k, l, m, n, j, l, o, m, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(j, l, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim6, Dim4, j, l, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim3, Dim1, m, o, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim3, Dim1, i, j, k, l, m, n, m, o, l, j, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(m, o, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim3, Dim1, m, o, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim4, Dim1, l, o, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim4, Dim1, i, j, k, l, m, n, l, o, m, j, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(l, o, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim4, Dim1, l, o, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim1, Dim3, m, o, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim1, Dim3, i, j, k, l, m, n, m, o, j, l, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(m, o, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim1, Dim3, m, o, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim1, Dim4, l, o, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim1, Dim4, i, j, k, l, m, n, l, o, j, m, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(l, o, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim1, Dim4, l, o, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim4, Dim3, j, o, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim4, Dim3, i, j, k, l, m, n, j, o, m, l, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(j, o, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim4, Dim3, j, o, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim3, Dim4, j, o, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim3, Dim4, i, j, k, l, m, n, j, o, l, m, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(j, o, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim3, Dim4, j, o, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim3, Dim1, o, m, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim1, i, j, k, l, m, n, o, m, l, j, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(o, m, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim3, Dim1, o, m, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim4, Dim1, o, l, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim1, i, j, k, l, m, n, o, l, m, j, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(o, l, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim4, Dim1, o, l, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim1, Dim3, o, m, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim3, i, j, k, l, m, n, o, m, j, l, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(o, m, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim1, Dim3, o, m, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim1, Dim4, o, l, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim4, i, j, k, l, m, n, o, l, j, m, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(o, l, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim1, Dim4, o, l, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim4, Dim3, o, j, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim3, i, j, k, l, m, n, o, j, m, l, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(o, j, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim4, Dim3, o, j, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim3, Dim4, o, j, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim4, i, j, k, l, m, n, o, j, l, m, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(o, j, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim3, Dim4, o, j, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim6, n, l, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim1, Dim6, i, j, k, l, m, n, n, l, j, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(n, l, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim6, n, l, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim6, l, n, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim1, Dim6, i, j, k, l, m, n, l, n, j, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(l, n, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim6, l, n, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim6, n, j, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim3, Dim6, i, j, k, l, m, n, n, j, l, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(n, j, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim6, n, j, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim6, l, j, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim5, Dim6, i, j, k, l, m, n, l, j, n, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(l, j, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim6, l, j, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim6, j, n, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim3, Dim6, i, j, k, l, m, n, j, n, l, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(j, n, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim6, j, n, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim6, j, l, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim5, Dim6, i, j, k, l, m, n, j, l, n, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(j, l, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim6, j, l, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim1, n, l, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim6, Dim1, i, j, k, l, m, n, n, l, o, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(n, l, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim1, n, l, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim1, l, n, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim6, Dim1, i, j, k, l, m, n, l, n, o, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(l, n, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim1, l, n, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim3, n, j, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim6, Dim3, i, j, k, l, m, n, n, j, o, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(n, j, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim3, n, j, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim6, Dim5, l, j, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim6, Dim5, i, j, k, l, m, n, l, j, o, n, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(l, j, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim6, Dim5, l, j, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim3, j, n, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim6, Dim3, i, j, k, l, m, n, j, n, o, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(j, n, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim3, j, n, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim6, Dim5, j, l, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim6, Dim5, i, j, k, l, m, n, j, l, o, n, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(j, l, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim6, Dim5, j, l, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim1, n, o, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim3, Dim1, i, j, k, l, m, n, n, o, l, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(n, o, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim1, n, o, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim5, Dim1, l, o, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim5, Dim1, i, j, k, l, m, n, l, o, n, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(l, o, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim5, Dim1, l, o, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim3, n, o, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim1, Dim3, i, j, k, l, m, n, n, o, j, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(n, o, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim3, n, o, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim1, Dim5, l, o, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim1, Dim5, i, j, k, l, m, n, l, o, j, n, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(l, o, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim1, Dim5, l, o, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim5, Dim3, j, o, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim5, Dim3, i, j, k, l, m, n, j, o, n, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(j, o, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim5, Dim3, j, o, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim3, Dim5, j, o, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim3, Dim5, i, j, k, l, m, n, j, o, l, n, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(j, o, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim3, Dim5, j, o, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim3, Dim1, o, n, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim1, i, j, k, l, m, n, o, n, l, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(o, n, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim3, Dim1, o, n, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim5, Dim1, o, l, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim1, i, j, k, l, m, n, o, l, n, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(o, l, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim5, Dim1, o, l, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim1, Dim3, o, n, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim3, i, j, k, l, m, n, o, n, j, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(o, n, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim1, Dim3, o, n, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim1, Dim5, o, l, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim5, i, j, k, l, m, n, o, l, j, n, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(o, l, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim1, Dim5, o, l, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim5, Dim3, o, j, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim3, i, j, k, l, m, n, o, j, n, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(o, j, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim5, Dim3, o, j, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim3, Dim5, o, j, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim5, i, j, k, l, m, n, o, j, l, n, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(o, j, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim3, Dim5, o, j, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim6, n, m, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim1, Dim6, i, j, k, l, m, n, n, m, j, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(n, m, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim6, n, m, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim6, m, n, j, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim1, Dim6, i, j, k, l, m, n, m, n, j, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(m, n, j, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim6, m, n, j, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim6, n, j, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim4, Dim6, i, j, k, l, m, n, n, j, m, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(n, j, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim6, n, j, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim6, m, j, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim5, Dim6, i, j, k, l, m, n, m, j, n, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(m, j, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim6, m, j, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim6, j, n, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim4, Dim6, i, j, k, l, m, n, j, n, m, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(j, n, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim6, j, n, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim6, j, m, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim5, Dim6, i, j, k, l, m, n, j, m, n, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(j, m, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim6, j, m, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim1, n, m, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim6, Dim1, i, j, k, l, m, n, n, m, o, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(n, m, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim1, n, m, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim1, m, n, o, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim6, Dim1, i, j, k, l, m, n, m, n, o, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(m, n, o, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim1, m, n, o, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim4, n, j, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim6, Dim4, i, j, k, l, m, n, n, j, o, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(n, j, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim4, n, j, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim6, Dim5, m, j, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim6, Dim5, i, j, k, l, m, n, m, j, o, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(m, j, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim6, Dim5, m, j, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim4, j, n, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim6, Dim4, i, j, k, l, m, n, j, n, o, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(j, n, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim4, j, n, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim6, Dim5, j, m, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim6, Dim5, i, j, k, l, m, n, j, m, o, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(j, m, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim6, Dim5, j, m, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim1, n, o, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim4, Dim1, i, j, k, l, m, n, n, o, m, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(n, o, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim1, n, o, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim5, Dim1, m, o, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim5, Dim1, i, j, k, l, m, n, m, o, n, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(m, o, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim5, Dim1, m, o, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim4, n, o, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim1, Dim4, i, j, k, l, m, n, n, o, j, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(n, o, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim4, n, o, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim1, Dim5, m, o, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim1, Dim5, i, j, k, l, m, n, m, o, j, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(m, o, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim1, Dim5, m, o, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim5, Dim4, j, o, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim5, Dim4, i, j, k, l, m, n, j, o, n, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(j, o, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim5, Dim4, j, o, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim4, Dim5, j, o, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim4, Dim5, i, j, k, l, m, n, j, o, m, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(j, o, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim4, Dim5, j, o, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim4, Dim1, o, n, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim1, i, j, k, l, m, n, o, n, m, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(o, n, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim4, Dim1, o, n, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim5, Dim1, o, m, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim1, i, j, k, l, m, n, o, m, n, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(o, m, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim5, Dim1, o, m, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim1, Dim4, o, n, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim4, i, j, k, l, m, n, o, n, j, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(o, n, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim1, Dim4, o, n, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim1, Dim5, o, m, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim5, i, j, k, l, m, n, o, m, j, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(o, m, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim1, Dim5, o, m, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim5, Dim4, o, j, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim4, i, j, k, l, m, n, o, j, n, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(o, j, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim5, Dim4, o, j, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim4, Dim5, o, j, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim5, i, j, k, l, m, n, o, j, m, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(o, j, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim4, Dim5, o, j, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim6, m, l, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim2, Dim6, i, j, k, l, m, n, m, l, k, o, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(m, l, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim6, m, l, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim6, l, m, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim2, Dim6, i, j, k, l, m, n, l, m, k, o, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(l, m, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim6, l, m, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim6, m, k, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim3, Dim6, i, j, k, l, m, n, m, k, l, o, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(m, k, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim6, m, k, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim6, l, k, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim4, Dim6, i, j, k, l, m, n, l, k, m, o, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(l, k, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim6, l, k, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim6, k, m, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim3, Dim6, i, j, k, l, m, n, k, m, l, o, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(k, m, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim6, k, m, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim6, k, l, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, k, l, m, n, k, l, m, o, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(k, l, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim6, k, l, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim6, Dim2, m, l, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim6, Dim2, i, j, k, l, m, n, m, l, o, k, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(m, l, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim6, Dim2, m, l, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim6, Dim2, l, m, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim6, Dim2, i, j, k, l, m, n, l, m, o, k, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(l, m, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim6, Dim2, l, m, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim6, Dim3, m, k, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim6, Dim3, i, j, k, l, m, n, m, k, o, l, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(m, k, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim6, Dim3, m, k, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim6, Dim4, l, k, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim6, Dim4, i, j, k, l, m, n, l, k, o, m, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(l, k, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim6, Dim4, l, k, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim6, Dim3, k, m, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim6, Dim3, i, j, k, l, m, n, k, m, o, l, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(k, m, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim6, Dim3, k, m, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim6, Dim4, k, l, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim6, Dim4, i, j, k, l, m, n, k, l, o, m, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(k, l, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim6, Dim4, k, l, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim3, Dim2, m, o, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim3, Dim2, i, j, k, l, m, n, m, o, l, k, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(m, o, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim3, Dim2, m, o, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim4, Dim2, l, o, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim4, Dim2, i, j, k, l, m, n, l, o, m, k, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(l, o, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim4, Dim2, l, o, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim2, Dim3, m, o, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim2, Dim3, i, j, k, l, m, n, m, o, k, l, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(m, o, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim2, Dim3, m, o, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim2, Dim4, l, o, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim2, Dim4, i, j, k, l, m, n, l, o, k, m, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(l, o, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim2, Dim4, l, o, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim4, Dim3, k, o, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim4, Dim3, i, j, k, l, m, n, k, o, m, l, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(k, o, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim4, Dim3, k, o, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim3, Dim4, k, o, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim3, Dim4, i, j, k, l, m, n, k, o, l, m, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(k, o, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim3, Dim4, k, o, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim3, Dim2, o, m, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim2, i, j, k, l, m, n, o, m, l, k, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(o, m, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim3, Dim2, o, m, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim4, Dim2, o, l, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim2, i, j, k, l, m, n, o, l, m, k, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(o, l, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim4, Dim2, o, l, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim2, Dim3, o, m, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim3, i, j, k, l, m, n, o, m, k, l, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(o, m, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim2, Dim3, o, m, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim2, Dim4, o, l, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim4, i, j, k, l, m, n, o, l, k, m, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(o, l, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim2, Dim4, o, l, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim4, Dim3, o, k, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim3, i, j, k, l, m, n, o, k, m, l, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(o, k, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim4, Dim3, o, k, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim3, Dim4, o, k, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim4, i, j, k, l, m, n, o, k, l, m, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(o, k, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim3, Dim4, o, k, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim6, n, l, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim2, Dim6, i, j, k, l, m, n, n, l, k, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(n, l, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim6, n, l, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim6, l, n, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim2, Dim6, i, j, k, l, m, n, l, n, k, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(l, n, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim6, l, n, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim6, n, k, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim3, Dim6, i, j, k, l, m, n, n, k, l, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(n, k, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim6, n, k, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim6, l, k, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim5, Dim6, i, j, k, l, m, n, l, k, n, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(l, k, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim6, l, k, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim6, k, n, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim3, Dim6, i, j, k, l, m, n, k, n, l, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(k, n, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim6, k, n, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim6, k, l, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim5, Dim6, i, j, k, l, m, n, k, l, n, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(k, l, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim6, k, l, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim2, n, l, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim6, Dim2, i, j, k, l, m, n, n, l, o, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(n, l, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim2, n, l, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim2, l, n, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim6, Dim2, i, j, k, l, m, n, l, n, o, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(l, n, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim2, l, n, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim3, n, k, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim6, Dim3, i, j, k, l, m, n, n, k, o, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(n, k, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim3, n, k, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim6, Dim5, l, k, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim6, Dim5, i, j, k, l, m, n, l, k, o, n, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(l, k, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim6, Dim5, l, k, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim3, k, n, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim6, Dim3, i, j, k, l, m, n, k, n, o, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(k, n, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim3, k, n, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim6, Dim5, k, l, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim6, Dim5, i, j, k, l, m, n, k, l, o, n, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(k, l, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim6, Dim5, k, l, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim2, n, o, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim3, Dim2, i, j, k, l, m, n, n, o, l, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(n, o, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim2, n, o, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim5, Dim2, l, o, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim5, Dim2, i, j, k, l, m, n, l, o, n, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(l, o, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim5, Dim2, l, o, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim3, n, o, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim2, Dim3, i, j, k, l, m, n, n, o, k, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(n, o, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim3, n, o, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim2, Dim5, l, o, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim2, Dim5, i, j, k, l, m, n, l, o, k, n, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(l, o, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim2, Dim5, l, o, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim5, Dim3, k, o, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim5, Dim3, i, j, k, l, m, n, k, o, n, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(k, o, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim5, Dim3, k, o, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim3, Dim5, k, o, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim3, Dim5, i, j, k, l, m, n, k, o, l, n, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(k, o, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim3, Dim5, k, o, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim3, Dim2, o, n, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim2, i, j, k, l, m, n, o, n, l, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(o, n, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim3, Dim2, o, n, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim5, Dim2, o, l, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim2, i, j, k, l, m, n, o, l, n, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(o, l, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim5, Dim2, o, l, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim2, Dim3, o, n, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim3, i, j, k, l, m, n, o, n, k, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(o, n, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim2, Dim3, o, n, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim2, Dim5, o, l, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim5, i, j, k, l, m, n, o, l, k, n, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(o, l, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim2, Dim5, o, l, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim5, Dim3, o, k, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim3, i, j, k, l, m, n, o, k, n, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(o, k, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim5, Dim3, o, k, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim3, Dim5, o, k, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim5, i, j, k, l, m, n, o, k, l, n, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(o, k, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim3, Dim5, o, k, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim6, n, m, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim2, Dim6, i, j, k, l, m, n, n, m, k, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(n, m, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim6, n, m, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim6, m, n, k, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim2, Dim6, i, j, k, l, m, n, m, n, k, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(m, n, k, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim6, m, n, k, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim6, n, k, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim4, Dim6, i, j, k, l, m, n, n, k, m, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(n, k, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim6, n, k, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim6, m, k, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim5, Dim6, i, j, k, l, m, n, m, k, n, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(m, k, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim6, m, k, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim6, k, n, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim4, Dim6, i, j, k, l, m, n, k, n, m, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(k, n, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim6, k, n, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim6, k, m, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim5, Dim6, i, j, k, l, m, n, k, m, n, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(k, m, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim6, k, m, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim2, n, m, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim6, Dim2, i, j, k, l, m, n, n, m, o, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(n, m, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim2, n, m, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim2, m, n, o, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim6, Dim2, i, j, k, l, m, n, m, n, o, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(m, n, o, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim2, m, n, o, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim4, n, k, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim6, Dim4, i, j, k, l, m, n, n, k, o, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(n, k, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim4, n, k, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim6, Dim5, m, k, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim6, Dim5, i, j, k, l, m, n, m, k, o, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(m, k, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim6, Dim5, m, k, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim4, k, n, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim6, Dim4, i, j, k, l, m, n, k, n, o, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(k, n, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim4, k, n, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim6, Dim5, k, m, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim6, Dim5, i, j, k, l, m, n, k, m, o, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(k, m, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim6, Dim5, k, m, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim2, n, o, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim4, Dim2, i, j, k, l, m, n, n, o, m, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(n, o, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim2, n, o, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim5, Dim2, m, o, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim5, Dim2, i, j, k, l, m, n, m, o, n, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(m, o, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim5, Dim2, m, o, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim4, n, o, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim2, Dim4, i, j, k, l, m, n, n, o, k, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(n, o, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim4, n, o, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim2, Dim5, m, o, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim2, Dim5, i, j, k, l, m, n, m, o, k, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(m, o, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim2, Dim5, m, o, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim5, Dim4, k, o, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim5, Dim4, i, j, k, l, m, n, k, o, n, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(k, o, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim5, Dim4, k, o, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim4, Dim5, k, o, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim4, Dim5, i, j, k, l, m, n, k, o, m, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(k, o, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim4, Dim5, k, o, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim4, Dim2, o, n, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim2, i, j, k, l, m, n, o, n, m, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(o, n, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim4, Dim2, o, n, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim5, Dim2, o, m, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim2, i, j, k, l, m, n, o, m, n, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(o, m, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim5, Dim2, o, m, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim2, Dim4, o, n, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim4, i, j, k, l, m, n, o, n, k, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(o, n, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim2, Dim4, o, n, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim2, Dim5, o, m, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim5, i, j, k, l, m, n, o, m, k, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(o, m, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim2, Dim5, o, m, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim5, Dim4, o, k, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim4, i, j, k, l, m, n, o, k, n, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(o, k, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim5, Dim4, o, k, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim4, Dim5, o, k, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim5, i, j, k, l, m, n, o, k, m, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(o, k, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim4, Dim5, o, k, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim6, n, m, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim3, Dim6, i, j, k, l, m, n, n, m, l, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(n, m, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim6, n, m, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim6, m, n, l, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim3, Dim6, i, j, k, l, m, n, m, n, l, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(m, n, l, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim6, m, n, l, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim6, n, l, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim4, Dim6, i, j, k, l, m, n, n, l, m, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(n, l, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim6, n, l, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim6, m, l, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim5, Dim6, i, j, k, l, m, n, m, l, n, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(m, l, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim6, m, l, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim6, l, n, m, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim4, Dim6, i, j, k, l, m, n, l, n, m, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(l, n, m, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim6, l, n, m, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim6, l, m, n, o> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, l, m, n, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(l, m, n, o, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim6, l, m, n, o> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim3, n, m, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim6, Dim3, i, j, k, l, m, n, n, m, o, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(n, m, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim3, n, m, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim3, m, n, o, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim6, Dim3, i, j, k, l, m, n, m, n, o, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(m, n, o, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim3, m, n, o, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim4, n, l, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim6, Dim4, i, j, k, l, m, n, n, l, o, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(n, l, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim4, n, l, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim6, Dim5, m, l, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim6, Dim5, i, j, k, l, m, n, m, l, o, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(m, l, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim6, Dim5, m, l, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim4, l, n, o, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim6, Dim4, i, j, k, l, m, n, l, n, o, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(l, n, o, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim4, l, n, o, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim6, Dim5, l, m, o, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim6, Dim5, i, j, k, l, m, n, l, m, o, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(l, m, o, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim6, Dim5, l, m, o, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim3, n, o, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim4, Dim3, i, j, k, l, m, n, n, o, m, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(n, o, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim3, n, o, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim5, Dim3, m, o, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim5, Dim3, i, j, k, l, m, n, m, o, n, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(m, o, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim5, Dim3, m, o, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim4, n, o, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim3, Dim4, i, j, k, l, m, n, n, o, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(n, o, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim4, n, o, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim3, Dim5, m, o, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim3, Dim5, i, j, k, l, m, n, m, o, l, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(m, o, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim3, Dim5, m, o, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim5, Dim4, l, o, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim5, Dim4, i, j, k, l, m, n, l, o, n, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(l, o, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim5, Dim4, l, o, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim4, Dim5, l, o, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim4, Dim5, i, j, k, l, m, n, l, o, m, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(l, o, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim4, Dim5, l, o, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim4, Dim3, o, n, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim3, i, j, k, l, m, n, o, n, m, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(o, n, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim4, Dim3, o, n, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim5, Dim3, o, m, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim3, i, j, k, l, m, n, o, m, n, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(o, m, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim5, Dim3, o, m, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim3, Dim4, o, n, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim4, i, j, k, l, m, n, o, n, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(o, n, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim3, Dim4, o, n, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim3, Dim5, o, m, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim5, i, j, k, l, m, n, o, m, l, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(o, m, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim3, Dim5, o, m, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim5, Dim4, o, l, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim4, i, j, k, l, m, n, o, l, n, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(o, l, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim5, Dim4, o, l, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim4, Dim5, o, l, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim5, i, j, k, l, m, n, o, l, m, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(o, l, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim4, Dim5, o, l, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

} //namespace FTensor

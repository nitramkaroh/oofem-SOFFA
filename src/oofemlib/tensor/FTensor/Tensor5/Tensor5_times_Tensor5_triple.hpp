#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, )*B(n, o, p, q, r, ) triple contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, int Dim9, char i, char j, char k, char l, char m, char n, char o, char p, char q, char r, int DimA, int DimB, int DimX, int DimY, int DimZ, int DimW, int DimV, char a, char b, char x, char y, char z, char w, char v>
class Tensor5_times_Tensor5_triple
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim9, n, o, p, q, r> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        for(int zz = 0; zz < DimZ; ++zz)
          {
            // Permutation is where the indices get checked.
            result += Permutation5<DimA, DimB, DimX, DimY, DimZ, a, b, x, y, z>().eval(
              iterA, N1, N2, xx, yy, zz) * Permutation5<DimX, DimY, DimZ, DimW, DimV, x, y, z, w, v>().eval(iterB, xx, yy, zz, N3, N4);
          }
    return result;
  }

  Tensor5_times_Tensor5_triple(
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &iter_a,
    const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim9, n, o, p, q, r> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, )*B(k, j, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim1, Dim0, Dim5, Dim6, k, j, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim0, Dim5, Dim6, i, j, k, l, m, k, j, i, n, o, Dim3, Dim4, Dim2, Dim1, Dim0, Dim5, Dim6, l, m, k, j, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, k, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim2, Dim0, Dim5, Dim6, j, k, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim0, Dim5, Dim6, i, j, k, l, m, j, k, i, n, o, Dim3, Dim4, Dim1, Dim2, Dim0, Dim5, Dim6, l, m, j, k, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, i, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim0, Dim1, Dim5, Dim6, k, i, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim1, Dim5, Dim6, i, j, k, l, m, k, i, j, n, o, Dim3, Dim4, Dim2, Dim0, Dim1, Dim5, Dim6, l, m, k, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, i, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim0, Dim2, Dim5, Dim6, j, i, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim2, Dim5, Dim6, i, j, k, l, m, j, i, k, n, o, Dim3, Dim4, Dim1, Dim0, Dim2, Dim5, Dim6, l, m, j, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, k, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim2, Dim1, Dim5, Dim6, i, k, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim1, Dim5, Dim6, i, j, k, l, m, i, k, j, n, o, Dim3, Dim4, Dim0, Dim2, Dim1, Dim5, Dim6, l, m, i, k, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, j, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, l, m, i, j, k, n, o, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, Dim6, l, m, i, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, j, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim1, Dim5, Dim0, Dim6, k, j, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, Dim0, Dim6, i, j, k, l, m, k, j, n, i, o, Dim3, Dim4, Dim2, Dim1, Dim5, Dim0, Dim6, l, m, k, j, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, k, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim2, Dim5, Dim0, Dim6, j, k, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, Dim0, Dim6, i, j, k, l, m, j, k, n, i, o, Dim3, Dim4, Dim1, Dim2, Dim5, Dim0, Dim6, l, m, j, k, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, i, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim0, Dim5, Dim1, Dim6, k, i, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, Dim1, Dim6, i, j, k, l, m, k, i, n, j, o, Dim3, Dim4, Dim2, Dim0, Dim5, Dim1, Dim6, l, m, k, i, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, i, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim0, Dim5, Dim2, Dim6, j, i, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, Dim2, Dim6, i, j, k, l, m, j, i, n, k, o, Dim3, Dim4, Dim1, Dim0, Dim5, Dim2, Dim6, l, m, j, i, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, k, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim2, Dim5, Dim1, Dim6, i, k, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, Dim1, Dim6, i, j, k, l, m, i, k, n, j, o, Dim3, Dim4, Dim0, Dim2, Dim5, Dim1, Dim6, l, m, i, k, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, j, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim1, Dim5, Dim2, Dim6, i, j, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim2, Dim6, i, j, k, l, m, i, j, n, k, o, Dim3, Dim4, Dim0, Dim1, Dim5, Dim2, Dim6, l, m, i, j, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, j, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim1, Dim5, Dim6, Dim0, k, j, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, Dim6, Dim0, i, j, k, l, m, k, j, n, o, i, Dim3, Dim4, Dim2, Dim1, Dim5, Dim6, Dim0, l, m, k, j, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, k, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim2, Dim5, Dim6, Dim0, j, k, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim0, i, j, k, l, m, j, k, n, o, i, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim0, l, m, j, k, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, i, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim0, Dim5, Dim6, Dim1, k, i, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, Dim6, Dim1, i, j, k, l, m, k, i, n, o, j, Dim3, Dim4, Dim2, Dim0, Dim5, Dim6, Dim1, l, m, k, i, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, i, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim0, Dim5, Dim6, Dim2, j, i, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, Dim6, Dim2, i, j, k, l, m, j, i, n, o, k, Dim3, Dim4, Dim1, Dim0, Dim5, Dim6, Dim2, l, m, j, i, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, k, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim2, Dim5, Dim6, Dim1, i, k, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim1, i, j, k, l, m, i, k, n, o, j, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim1, l, m, i, k, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, j, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim1, Dim5, Dim6, Dim2, i, j, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim2, i, j, k, l, m, i, j, n, o, k, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim2, l, m, i, j, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim1, Dim0, Dim6, k, n, j, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, Dim0, Dim6, i, j, k, l, m, k, n, j, i, o, Dim3, Dim4, Dim2, Dim5, Dim1, Dim0, Dim6, l, m, k, n, j, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim2, Dim0, Dim6, j, n, k, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, Dim0, Dim6, i, j, k, l, m, j, n, k, i, o, Dim3, Dim4, Dim1, Dim5, Dim2, Dim0, Dim6, l, m, j, n, k, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim0, Dim1, Dim6, k, n, i, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, Dim1, Dim6, i, j, k, l, m, k, n, i, j, o, Dim3, Dim4, Dim2, Dim5, Dim0, Dim1, Dim6, l, m, k, n, i, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim0, Dim2, Dim6, j, n, i, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, Dim2, Dim6, i, j, k, l, m, j, n, i, k, o, Dim3, Dim4, Dim1, Dim5, Dim0, Dim2, Dim6, l, m, j, n, i, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim2, Dim1, Dim6, i, n, k, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, Dim1, Dim6, i, j, k, l, m, i, n, k, j, o, Dim3, Dim4, Dim0, Dim5, Dim2, Dim1, Dim6, l, m, i, n, k, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim1, Dim2, Dim6, i, n, j, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, Dim2, Dim6, i, j, k, l, m, i, n, j, k, o, Dim3, Dim4, Dim0, Dim5, Dim1, Dim2, Dim6, l, m, i, n, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim1, Dim6, Dim0, k, n, j, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, Dim6, Dim0, i, j, k, l, m, k, n, j, o, i, Dim3, Dim4, Dim2, Dim5, Dim1, Dim6, Dim0, l, m, k, n, j, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim2, Dim6, Dim0, j, n, k, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, Dim6, Dim0, i, j, k, l, m, j, n, k, o, i, Dim3, Dim4, Dim1, Dim5, Dim2, Dim6, Dim0, l, m, j, n, k, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim0, Dim6, Dim1, k, n, i, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, Dim6, Dim1, i, j, k, l, m, k, n, i, o, j, Dim3, Dim4, Dim2, Dim5, Dim0, Dim6, Dim1, l, m, k, n, i, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim0, Dim6, Dim2, j, n, i, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, Dim6, Dim2, i, j, k, l, m, j, n, i, o, k, Dim3, Dim4, Dim1, Dim5, Dim0, Dim6, Dim2, l, m, j, n, i, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim2, Dim6, Dim1, i, n, k, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, Dim6, Dim1, i, j, k, l, m, i, n, k, o, j, Dim3, Dim4, Dim0, Dim5, Dim2, Dim6, Dim1, l, m, i, n, k, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim1, Dim6, Dim2, i, n, j, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, Dim6, Dim2, i, j, k, l, m, i, n, j, o, k, Dim3, Dim4, Dim0, Dim5, Dim1, Dim6, Dim2, l, m, i, n, j, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim1, Dim0, k, n, o, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim1, Dim0, i, j, k, l, m, k, n, o, j, i, Dim3, Dim4, Dim2, Dim5, Dim6, Dim1, Dim0, l, m, k, n, o, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim2, Dim0, j, n, o, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim2, Dim0, i, j, k, l, m, j, n, o, k, i, Dim3, Dim4, Dim1, Dim5, Dim6, Dim2, Dim0, l, m, j, n, o, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim0, Dim1, k, n, o, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim0, Dim1, i, j, k, l, m, k, n, o, i, j, Dim3, Dim4, Dim2, Dim5, Dim6, Dim0, Dim1, l, m, k, n, o, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim0, Dim2, j, n, o, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim0, Dim2, i, j, k, l, m, j, n, o, i, k, Dim3, Dim4, Dim1, Dim5, Dim6, Dim0, Dim2, l, m, j, n, o, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim2, Dim1, i, n, o, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim2, Dim1, i, j, k, l, m, i, n, o, k, j, Dim3, Dim4, Dim0, Dim5, Dim6, Dim2, Dim1, l, m, i, n, o, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim1, Dim2, i, n, o, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim1, Dim2, i, j, k, l, m, i, n, o, j, k, Dim3, Dim4, Dim0, Dim5, Dim6, Dim1, Dim2, l, m, i, n, o, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim1, Dim0, Dim6, n, k, j, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim0, Dim6, i, j, k, l, m, n, k, j, i, o, Dim3, Dim4, Dim5, Dim2, Dim1, Dim0, Dim6, l, m, n, k, j, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim2, Dim0, Dim6, n, j, k, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim0, Dim6, i, j, k, l, m, n, j, k, i, o, Dim3, Dim4, Dim5, Dim1, Dim2, Dim0, Dim6, l, m, n, j, k, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim0, Dim1, Dim6, n, k, i, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim1, Dim6, i, j, k, l, m, n, k, i, j, o, Dim3, Dim4, Dim5, Dim2, Dim0, Dim1, Dim6, l, m, n, k, i, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim0, Dim2, Dim6, n, j, i, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim2, Dim6, i, j, k, l, m, n, j, i, k, o, Dim3, Dim4, Dim5, Dim1, Dim0, Dim2, Dim6, l, m, n, j, i, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim2, Dim1, Dim6, n, i, k, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim1, Dim6, i, j, k, l, m, n, i, k, j, o, Dim3, Dim4, Dim5, Dim0, Dim2, Dim1, Dim6, l, m, n, i, k, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim1, Dim2, Dim6, n, i, j, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, i, j, k, l, m, n, i, j, k, o, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim6, l, m, n, i, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim1, Dim6, Dim0, n, k, j, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, Dim0, i, j, k, l, m, n, k, j, o, i, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, Dim0, l, m, n, k, j, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim2, Dim6, Dim0, n, j, k, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim0, i, j, k, l, m, n, j, k, o, i, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim0, l, m, n, j, k, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim0, Dim6, Dim1, n, k, i, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, Dim1, i, j, k, l, m, n, k, i, o, j, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, Dim1, l, m, n, k, i, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim0, Dim6, Dim2, n, j, i, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, Dim2, i, j, k, l, m, n, j, i, o, k, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, Dim2, l, m, n, j, i, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim2, Dim6, Dim1, n, i, k, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim1, i, j, k, l, m, n, i, k, o, j, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim1, l, m, n, i, k, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim1, Dim6, Dim2, n, i, j, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim2, i, j, k, l, m, n, i, j, o, k, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim2, l, m, n, i, j, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim1, Dim0, n, k, o, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, Dim0, i, j, k, l, m, n, k, o, j, i, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, Dim0, l, m, n, k, o, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim2, Dim0, n, j, o, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, Dim0, i, j, k, l, m, n, j, o, k, i, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, Dim0, l, m, n, j, o, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim0, Dim1, n, k, o, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, Dim1, i, j, k, l, m, n, k, o, i, j, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, Dim1, l, m, n, k, o, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim0, Dim2, n, j, o, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, Dim2, i, j, k, l, m, n, j, o, i, k, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, Dim2, l, m, n, j, o, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim2, Dim1, n, i, o, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, Dim1, i, j, k, l, m, n, i, o, k, j, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, Dim1, l, m, n, i, o, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim1, Dim2, n, i, o, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, Dim2, i, j, k, l, m, n, i, o, j, k, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, Dim2, l, m, n, i, o, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim1, Dim0, n, o, k, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim0, i, j, k, l, m, n, o, k, j, i, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim0, l, m, n, o, k, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim2, Dim0, n, o, j, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim0, i, j, k, l, m, n, o, j, k, i, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim0, l, m, n, o, j, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim0, Dim1, n, o, k, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim1, i, j, k, l, m, n, o, k, i, j, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim1, l, m, n, o, k, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim0, Dim2, n, o, j, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim2, i, j, k, l, m, n, o, j, i, k, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim2, l, m, n, o, j, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim2, Dim1, n, o, i, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim1, i, j, k, l, m, n, o, i, k, j, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim1, l, m, n, o, i, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim1, Dim2, n, o, i, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim2, i, j, k, l, m, n, o, i, j, k, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim2, l, m, n, o, i, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, j, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim1, Dim0, Dim5, Dim6, l, j, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim0, Dim5, Dim6, i, j, k, l, m, l, j, i, n, o, Dim2, Dim4, Dim3, Dim1, Dim0, Dim5, Dim6, k, m, l, j, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, l, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim3, Dim0, Dim5, Dim6, j, l, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim0, Dim5, Dim6, i, j, k, l, m, j, l, i, n, o, Dim2, Dim4, Dim1, Dim3, Dim0, Dim5, Dim6, k, m, j, l, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, i, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim0, Dim1, Dim5, Dim6, l, i, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim1, Dim5, Dim6, i, j, k, l, m, l, i, j, n, o, Dim2, Dim4, Dim3, Dim0, Dim1, Dim5, Dim6, k, m, l, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, i, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim0, Dim3, Dim5, Dim6, j, i, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim3, Dim5, Dim6, i, j, k, l, m, j, i, l, n, o, Dim2, Dim4, Dim1, Dim0, Dim3, Dim5, Dim6, k, m, j, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, l, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim3, Dim1, Dim5, Dim6, i, l, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim1, Dim5, Dim6, i, j, k, l, m, i, l, j, n, o, Dim2, Dim4, Dim0, Dim3, Dim1, Dim5, Dim6, k, m, i, l, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, j, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, k, l, m, i, j, l, n, o, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, Dim6, k, m, i, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, j, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim1, Dim5, Dim0, Dim6, l, j, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, Dim0, Dim6, i, j, k, l, m, l, j, n, i, o, Dim2, Dim4, Dim3, Dim1, Dim5, Dim0, Dim6, k, m, l, j, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, l, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim3, Dim5, Dim0, Dim6, j, l, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, Dim0, Dim6, i, j, k, l, m, j, l, n, i, o, Dim2, Dim4, Dim1, Dim3, Dim5, Dim0, Dim6, k, m, j, l, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, i, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim0, Dim5, Dim1, Dim6, l, i, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, Dim1, Dim6, i, j, k, l, m, l, i, n, j, o, Dim2, Dim4, Dim3, Dim0, Dim5, Dim1, Dim6, k, m, l, i, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, i, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim0, Dim5, Dim3, Dim6, j, i, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, Dim3, Dim6, i, j, k, l, m, j, i, n, l, o, Dim2, Dim4, Dim1, Dim0, Dim5, Dim3, Dim6, k, m, j, i, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, l, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim3, Dim5, Dim1, Dim6, i, l, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, Dim1, Dim6, i, j, k, l, m, i, l, n, j, o, Dim2, Dim4, Dim0, Dim3, Dim5, Dim1, Dim6, k, m, i, l, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, j, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim1, Dim5, Dim3, Dim6, i, j, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim3, Dim6, i, j, k, l, m, i, j, n, l, o, Dim2, Dim4, Dim0, Dim1, Dim5, Dim3, Dim6, k, m, i, j, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, j, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim1, Dim5, Dim6, Dim0, l, j, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, Dim6, Dim0, i, j, k, l, m, l, j, n, o, i, Dim2, Dim4, Dim3, Dim1, Dim5, Dim6, Dim0, k, m, l, j, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, l, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim3, Dim5, Dim6, Dim0, j, l, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, Dim6, Dim0, i, j, k, l, m, j, l, n, o, i, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim0, k, m, j, l, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, i, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim0, Dim5, Dim6, Dim1, l, i, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, Dim6, Dim1, i, j, k, l, m, l, i, n, o, j, Dim2, Dim4, Dim3, Dim0, Dim5, Dim6, Dim1, k, m, l, i, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, i, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim0, Dim5, Dim6, Dim3, j, i, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, Dim6, Dim3, i, j, k, l, m, j, i, n, o, l, Dim2, Dim4, Dim1, Dim0, Dim5, Dim6, Dim3, k, m, j, i, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, l, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim3, Dim5, Dim6, Dim1, i, l, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, Dim6, Dim1, i, j, k, l, m, i, l, n, o, j, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim1, k, m, i, l, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, j, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim1, Dim5, Dim6, Dim3, i, j, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim3, i, j, k, l, m, i, j, n, o, l, Dim2, Dim4, Dim0, Dim1, Dim5, Dim6, Dim3, k, m, i, j, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim1, Dim0, Dim6, l, n, j, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, Dim0, Dim6, i, j, k, l, m, l, n, j, i, o, Dim2, Dim4, Dim3, Dim5, Dim1, Dim0, Dim6, k, m, l, n, j, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim3, Dim0, Dim6, j, n, l, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, Dim0, Dim6, i, j, k, l, m, j, n, l, i, o, Dim2, Dim4, Dim1, Dim5, Dim3, Dim0, Dim6, k, m, j, n, l, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim0, Dim1, Dim6, l, n, i, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, Dim1, Dim6, i, j, k, l, m, l, n, i, j, o, Dim2, Dim4, Dim3, Dim5, Dim0, Dim1, Dim6, k, m, l, n, i, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim0, Dim3, Dim6, j, n, i, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, Dim3, Dim6, i, j, k, l, m, j, n, i, l, o, Dim2, Dim4, Dim1, Dim5, Dim0, Dim3, Dim6, k, m, j, n, i, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim3, Dim1, Dim6, i, n, l, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, Dim1, Dim6, i, j, k, l, m, i, n, l, j, o, Dim2, Dim4, Dim0, Dim5, Dim3, Dim1, Dim6, k, m, i, n, l, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim1, Dim3, Dim6, i, n, j, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, Dim3, Dim6, i, j, k, l, m, i, n, j, l, o, Dim2, Dim4, Dim0, Dim5, Dim1, Dim3, Dim6, k, m, i, n, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim1, Dim6, Dim0, l, n, j, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, Dim6, Dim0, i, j, k, l, m, l, n, j, o, i, Dim2, Dim4, Dim3, Dim5, Dim1, Dim6, Dim0, k, m, l, n, j, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim3, Dim6, Dim0, j, n, l, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, Dim6, Dim0, i, j, k, l, m, j, n, l, o, i, Dim2, Dim4, Dim1, Dim5, Dim3, Dim6, Dim0, k, m, j, n, l, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim0, Dim6, Dim1, l, n, i, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, Dim6, Dim1, i, j, k, l, m, l, n, i, o, j, Dim2, Dim4, Dim3, Dim5, Dim0, Dim6, Dim1, k, m, l, n, i, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim0, Dim6, Dim3, j, n, i, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, Dim6, Dim3, i, j, k, l, m, j, n, i, o, l, Dim2, Dim4, Dim1, Dim5, Dim0, Dim6, Dim3, k, m, j, n, i, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim3, Dim6, Dim1, i, n, l, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, Dim6, Dim1, i, j, k, l, m, i, n, l, o, j, Dim2, Dim4, Dim0, Dim5, Dim3, Dim6, Dim1, k, m, i, n, l, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim1, Dim6, Dim3, i, n, j, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, Dim6, Dim3, i, j, k, l, m, i, n, j, o, l, Dim2, Dim4, Dim0, Dim5, Dim1, Dim6, Dim3, k, m, i, n, j, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim1, Dim0, l, n, o, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim1, Dim0, i, j, k, l, m, l, n, o, j, i, Dim2, Dim4, Dim3, Dim5, Dim6, Dim1, Dim0, k, m, l, n, o, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim3, Dim0, j, n, o, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim3, Dim0, i, j, k, l, m, j, n, o, l, i, Dim2, Dim4, Dim1, Dim5, Dim6, Dim3, Dim0, k, m, j, n, o, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim0, Dim1, l, n, o, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim0, Dim1, i, j, k, l, m, l, n, o, i, j, Dim2, Dim4, Dim3, Dim5, Dim6, Dim0, Dim1, k, m, l, n, o, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim0, Dim3, j, n, o, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim0, Dim3, i, j, k, l, m, j, n, o, i, l, Dim2, Dim4, Dim1, Dim5, Dim6, Dim0, Dim3, k, m, j, n, o, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim3, Dim1, i, n, o, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim3, Dim1, i, j, k, l, m, i, n, o, l, j, Dim2, Dim4, Dim0, Dim5, Dim6, Dim3, Dim1, k, m, i, n, o, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim1, Dim3, i, n, o, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim1, Dim3, i, j, k, l, m, i, n, o, j, l, Dim2, Dim4, Dim0, Dim5, Dim6, Dim1, Dim3, k, m, i, n, o, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim1, Dim0, Dim6, n, l, j, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim0, Dim6, i, j, k, l, m, n, l, j, i, o, Dim2, Dim4, Dim5, Dim3, Dim1, Dim0, Dim6, k, m, n, l, j, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim3, Dim0, Dim6, n, j, l, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim0, Dim6, i, j, k, l, m, n, j, l, i, o, Dim2, Dim4, Dim5, Dim1, Dim3, Dim0, Dim6, k, m, n, j, l, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim0, Dim1, Dim6, n, l, i, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim1, Dim6, i, j, k, l, m, n, l, i, j, o, Dim2, Dim4, Dim5, Dim3, Dim0, Dim1, Dim6, k, m, n, l, i, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim0, Dim3, Dim6, n, j, i, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim3, Dim6, i, j, k, l, m, n, j, i, l, o, Dim2, Dim4, Dim5, Dim1, Dim0, Dim3, Dim6, k, m, n, j, i, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim3, Dim1, Dim6, n, i, l, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim1, Dim6, i, j, k, l, m, n, i, l, j, o, Dim2, Dim4, Dim5, Dim0, Dim3, Dim1, Dim6, k, m, n, i, l, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim1, Dim3, Dim6, n, i, j, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, i, j, k, l, m, n, i, j, l, o, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, Dim6, k, m, n, i, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim1, Dim6, Dim0, n, l, j, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim6, Dim0, i, j, k, l, m, n, l, j, o, i, Dim2, Dim4, Dim5, Dim3, Dim1, Dim6, Dim0, k, m, n, l, j, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim3, Dim6, Dim0, n, j, l, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim6, Dim0, i, j, k, l, m, n, j, l, o, i, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim0, k, m, n, j, l, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim0, Dim6, Dim1, n, l, i, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim6, Dim1, i, j, k, l, m, n, l, i, o, j, Dim2, Dim4, Dim5, Dim3, Dim0, Dim6, Dim1, k, m, n, l, i, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim0, Dim6, Dim3, n, j, i, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, Dim3, i, j, k, l, m, n, j, i, o, l, Dim2, Dim4, Dim5, Dim1, Dim0, Dim6, Dim3, k, m, n, j, i, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim3, Dim6, Dim1, n, i, l, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim6, Dim1, i, j, k, l, m, n, i, l, o, j, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim1, k, m, n, i, l, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim1, Dim6, Dim3, n, i, j, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim3, i, j, k, l, m, n, i, j, o, l, Dim2, Dim4, Dim5, Dim0, Dim1, Dim6, Dim3, k, m, n, i, j, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim1, Dim0, n, l, o, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim1, Dim0, i, j, k, l, m, n, l, o, j, i, Dim2, Dim4, Dim5, Dim3, Dim6, Dim1, Dim0, k, m, n, l, o, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim3, Dim0, n, j, o, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim3, Dim0, i, j, k, l, m, n, j, o, l, i, Dim2, Dim4, Dim5, Dim1, Dim6, Dim3, Dim0, k, m, n, j, o, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim0, Dim1, n, l, o, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim0, Dim1, i, j, k, l, m, n, l, o, i, j, Dim2, Dim4, Dim5, Dim3, Dim6, Dim0, Dim1, k, m, n, l, o, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim0, Dim3, n, j, o, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, Dim3, i, j, k, l, m, n, j, o, i, l, Dim2, Dim4, Dim5, Dim1, Dim6, Dim0, Dim3, k, m, n, j, o, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim3, Dim1, n, i, o, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim3, Dim1, i, j, k, l, m, n, i, o, l, j, Dim2, Dim4, Dim5, Dim0, Dim6, Dim3, Dim1, k, m, n, i, o, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim1, Dim3, n, i, o, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, Dim3, i, j, k, l, m, n, i, o, j, l, Dim2, Dim4, Dim5, Dim0, Dim6, Dim1, Dim3, k, m, n, i, o, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim1, Dim0, n, o, l, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim0, i, j, k, l, m, n, o, l, j, i, Dim2, Dim4, Dim5, Dim6, Dim3, Dim1, Dim0, k, m, n, o, l, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim3, Dim0, n, o, j, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim0, i, j, k, l, m, n, o, j, l, i, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim0, k, m, n, o, j, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim0, Dim1, n, o, l, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim1, i, j, k, l, m, n, o, l, i, j, Dim2, Dim4, Dim5, Dim6, Dim3, Dim0, Dim1, k, m, n, o, l, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim0, Dim3, n, o, j, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim3, i, j, k, l, m, n, o, j, i, l, Dim2, Dim4, Dim5, Dim6, Dim1, Dim0, Dim3, k, m, n, o, j, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim3, Dim1, n, o, i, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim1, i, j, k, l, m, n, o, i, l, j, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim1, k, m, n, o, i, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim1, Dim3, n, o, i, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim3, i, j, k, l, m, n, o, i, j, l, Dim2, Dim4, Dim5, Dim6, Dim0, Dim1, Dim3, k, m, n, o, i, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, j, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim1, Dim0, Dim5, Dim6, m, j, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim0, Dim5, Dim6, i, j, k, l, m, m, j, i, n, o, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, Dim6, k, l, m, j, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, m, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim4, Dim0, Dim5, Dim6, j, m, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim0, Dim5, Dim6, i, j, k, l, m, j, m, i, n, o, Dim2, Dim3, Dim1, Dim4, Dim0, Dim5, Dim6, k, l, j, m, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, i, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim0, Dim1, Dim5, Dim6, m, i, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim1, Dim5, Dim6, i, j, k, l, m, m, i, j, n, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, i, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim0, Dim4, Dim5, Dim6, j, i, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim4, Dim5, Dim6, i, j, k, l, m, j, i, m, n, o, Dim2, Dim3, Dim1, Dim0, Dim4, Dim5, Dim6, k, l, j, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, m, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim4, Dim1, Dim5, Dim6, i, m, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim1, Dim5, Dim6, i, j, k, l, m, i, m, j, n, o, Dim2, Dim3, Dim0, Dim4, Dim1, Dim5, Dim6, k, l, i, m, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, j, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, k, l, m, i, j, m, n, o, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, Dim6, k, l, i, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, j, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim1, Dim5, Dim0, Dim6, m, j, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim5, Dim0, Dim6, i, j, k, l, m, m, j, n, i, o, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, Dim6, k, l, m, j, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, m, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim4, Dim5, Dim0, Dim6, j, m, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim5, Dim0, Dim6, i, j, k, l, m, j, m, n, i, o, Dim2, Dim3, Dim1, Dim4, Dim5, Dim0, Dim6, k, l, j, m, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, i, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim0, Dim5, Dim1, Dim6, m, i, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim5, Dim1, Dim6, i, j, k, l, m, m, i, n, j, o, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, Dim6, k, l, m, i, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, i, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim0, Dim5, Dim4, Dim6, j, i, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, Dim4, Dim6, i, j, k, l, m, j, i, n, m, o, Dim2, Dim3, Dim1, Dim0, Dim5, Dim4, Dim6, k, l, j, i, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, m, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim4, Dim5, Dim1, Dim6, i, m, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim5, Dim1, Dim6, i, j, k, l, m, i, m, n, j, o, Dim2, Dim3, Dim0, Dim4, Dim5, Dim1, Dim6, k, l, i, m, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, j, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim1, Dim5, Dim4, Dim6, i, j, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim4, Dim6, i, j, k, l, m, i, j, n, m, o, Dim2, Dim3, Dim0, Dim1, Dim5, Dim4, Dim6, k, l, i, j, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, j, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim1, Dim5, Dim6, Dim0, m, j, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim5, Dim6, Dim0, i, j, k, l, m, m, j, n, o, i, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim0, k, l, m, j, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, m, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim4, Dim5, Dim6, Dim0, j, m, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim5, Dim6, Dim0, i, j, k, l, m, j, m, n, o, i, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim0, k, l, j, m, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, i, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim0, Dim5, Dim6, Dim1, m, i, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim5, Dim6, Dim1, i, j, k, l, m, m, i, n, o, j, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim1, k, l, m, i, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, i, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim0, Dim5, Dim6, Dim4, j, i, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, Dim6, Dim4, i, j, k, l, m, j, i, n, o, m, Dim2, Dim3, Dim1, Dim0, Dim5, Dim6, Dim4, k, l, j, i, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, m, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim4, Dim5, Dim6, Dim1, i, m, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim5, Dim6, Dim1, i, j, k, l, m, i, m, n, o, j, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim1, k, l, i, m, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, j, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim1, Dim5, Dim6, Dim4, i, j, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim4, i, j, k, l, m, i, j, n, o, m, Dim2, Dim3, Dim0, Dim1, Dim5, Dim6, Dim4, k, l, i, j, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim1, Dim0, Dim6, m, n, j, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim1, Dim0, Dim6, i, j, k, l, m, m, n, j, i, o, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, k, l, m, n, j, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim4, Dim0, Dim6, j, n, m, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim4, Dim0, Dim6, i, j, k, l, m, j, n, m, i, o, Dim2, Dim3, Dim1, Dim5, Dim4, Dim0, Dim6, k, l, j, n, m, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim0, Dim1, Dim6, m, n, i, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim0, Dim1, Dim6, i, j, k, l, m, m, n, i, j, o, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, k, l, m, n, i, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim0, Dim4, Dim6, j, n, i, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, Dim4, Dim6, i, j, k, l, m, j, n, i, m, o, Dim2, Dim3, Dim1, Dim5, Dim0, Dim4, Dim6, k, l, j, n, i, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim4, Dim1, Dim6, i, n, m, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim4, Dim1, Dim6, i, j, k, l, m, i, n, m, j, o, Dim2, Dim3, Dim0, Dim5, Dim4, Dim1, Dim6, k, l, i, n, m, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim1, Dim4, Dim6, i, n, j, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, Dim4, Dim6, i, j, k, l, m, i, n, j, m, o, Dim2, Dim3, Dim0, Dim5, Dim1, Dim4, Dim6, k, l, i, n, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim1, Dim6, Dim0, m, n, j, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim1, Dim6, Dim0, i, j, k, l, m, m, n, j, o, i, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, k, l, m, n, j, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim4, Dim6, Dim0, j, n, m, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim4, Dim6, Dim0, i, j, k, l, m, j, n, m, o, i, Dim2, Dim3, Dim1, Dim5, Dim4, Dim6, Dim0, k, l, j, n, m, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim0, Dim6, Dim1, m, n, i, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim0, Dim6, Dim1, i, j, k, l, m, m, n, i, o, j, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, k, l, m, n, i, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim0, Dim6, Dim4, j, n, i, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, Dim6, Dim4, i, j, k, l, m, j, n, i, o, m, Dim2, Dim3, Dim1, Dim5, Dim0, Dim6, Dim4, k, l, j, n, i, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim4, Dim6, Dim1, i, n, m, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim4, Dim6, Dim1, i, j, k, l, m, i, n, m, o, j, Dim2, Dim3, Dim0, Dim5, Dim4, Dim6, Dim1, k, l, i, n, m, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim1, Dim6, Dim4, i, n, j, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, Dim6, Dim4, i, j, k, l, m, i, n, j, o, m, Dim2, Dim3, Dim0, Dim5, Dim1, Dim6, Dim4, k, l, i, n, j, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim1, Dim0, m, n, o, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim1, Dim0, i, j, k, l, m, m, n, o, j, i, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, k, l, m, n, o, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim4, Dim0, j, n, o, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim4, Dim0, i, j, k, l, m, j, n, o, m, i, Dim2, Dim3, Dim1, Dim5, Dim6, Dim4, Dim0, k, l, j, n, o, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim0, Dim1, m, n, o, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim0, Dim1, i, j, k, l, m, m, n, o, i, j, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, k, l, m, n, o, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim0, Dim4, j, n, o, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim0, Dim4, i, j, k, l, m, j, n, o, i, m, Dim2, Dim3, Dim1, Dim5, Dim6, Dim0, Dim4, k, l, j, n, o, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim4, Dim1, i, n, o, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim4, Dim1, i, j, k, l, m, i, n, o, m, j, Dim2, Dim3, Dim0, Dim5, Dim6, Dim4, Dim1, k, l, i, n, o, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim1, Dim4, i, n, o, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim1, Dim4, i, j, k, l, m, i, n, o, j, m, Dim2, Dim3, Dim0, Dim5, Dim6, Dim1, Dim4, k, l, i, n, o, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim1, Dim0, Dim6, n, m, j, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim0, Dim6, i, j, k, l, m, n, m, j, i, o, Dim2, Dim3, Dim5, Dim4, Dim1, Dim0, Dim6, k, l, n, m, j, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim4, Dim0, Dim6, n, j, m, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim0, Dim6, i, j, k, l, m, n, j, m, i, o, Dim2, Dim3, Dim5, Dim1, Dim4, Dim0, Dim6, k, l, n, j, m, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim0, Dim1, Dim6, n, m, i, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim1, Dim6, i, j, k, l, m, n, m, i, j, o, Dim2, Dim3, Dim5, Dim4, Dim0, Dim1, Dim6, k, l, n, m, i, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim0, Dim4, Dim6, n, j, i, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim4, Dim6, i, j, k, l, m, n, j, i, m, o, Dim2, Dim3, Dim5, Dim1, Dim0, Dim4, Dim6, k, l, n, j, i, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim4, Dim1, Dim6, n, i, m, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim1, Dim6, i, j, k, l, m, n, i, m, j, o, Dim2, Dim3, Dim5, Dim0, Dim4, Dim1, Dim6, k, l, n, i, m, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim1, Dim4, Dim6, n, i, j, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim4, Dim6, i, j, k, l, m, n, i, j, m, o, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, Dim6, k, l, n, i, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim1, Dim6, Dim0, n, m, j, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim6, Dim0, i, j, k, l, m, n, m, j, o, i, Dim2, Dim3, Dim5, Dim4, Dim1, Dim6, Dim0, k, l, n, m, j, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim4, Dim6, Dim0, n, j, m, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim6, Dim0, i, j, k, l, m, n, j, m, o, i, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim0, k, l, n, j, m, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim0, Dim6, Dim1, n, m, i, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim6, Dim1, i, j, k, l, m, n, m, i, o, j, Dim2, Dim3, Dim5, Dim4, Dim0, Dim6, Dim1, k, l, n, m, i, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim0, Dim6, Dim4, n, j, i, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, Dim4, i, j, k, l, m, n, j, i, o, m, Dim2, Dim3, Dim5, Dim1, Dim0, Dim6, Dim4, k, l, n, j, i, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim4, Dim6, Dim1, n, i, m, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim6, Dim1, i, j, k, l, m, n, i, m, o, j, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim1, k, l, n, i, m, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim1, Dim6, Dim4, n, i, j, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim4, i, j, k, l, m, n, i, j, o, m, Dim2, Dim3, Dim5, Dim0, Dim1, Dim6, Dim4, k, l, n, i, j, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim1, Dim0, n, m, o, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim1, Dim0, i, j, k, l, m, n, m, o, j, i, Dim2, Dim3, Dim5, Dim4, Dim6, Dim1, Dim0, k, l, n, m, o, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim4, Dim0, n, j, o, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim4, Dim0, i, j, k, l, m, n, j, o, m, i, Dim2, Dim3, Dim5, Dim1, Dim6, Dim4, Dim0, k, l, n, j, o, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim0, Dim1, n, m, o, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim0, Dim1, i, j, k, l, m, n, m, o, i, j, Dim2, Dim3, Dim5, Dim4, Dim6, Dim0, Dim1, k, l, n, m, o, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim0, Dim4, n, j, o, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, Dim4, i, j, k, l, m, n, j, o, i, m, Dim2, Dim3, Dim5, Dim1, Dim6, Dim0, Dim4, k, l, n, j, o, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim4, Dim1, n, i, o, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim4, Dim1, i, j, k, l, m, n, i, o, m, j, Dim2, Dim3, Dim5, Dim0, Dim6, Dim4, Dim1, k, l, n, i, o, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim1, Dim4, n, i, o, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, Dim4, i, j, k, l, m, n, i, o, j, m, Dim2, Dim3, Dim5, Dim0, Dim6, Dim1, Dim4, k, l, n, i, o, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim1, Dim0, n, o, m, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim0, i, j, k, l, m, n, o, m, j, i, Dim2, Dim3, Dim5, Dim6, Dim4, Dim1, Dim0, k, l, n, o, m, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim4, Dim0, n, o, j, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim0, i, j, k, l, m, n, o, j, m, i, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim0, k, l, n, o, j, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim0, Dim1, n, o, m, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim1, i, j, k, l, m, n, o, m, i, j, Dim2, Dim3, Dim5, Dim6, Dim4, Dim0, Dim1, k, l, n, o, m, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim0, Dim4, n, o, j, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim4, i, j, k, l, m, n, o, j, i, m, Dim2, Dim3, Dim5, Dim6, Dim1, Dim0, Dim4, k, l, n, o, j, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim4, Dim1, n, o, i, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim1, i, j, k, l, m, n, o, i, m, j, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim1, k, l, n, o, i, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim1, Dim4, n, o, i, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim4, i, j, k, l, m, n, o, i, j, m, Dim2, Dim3, Dim5, Dim6, Dim0, Dim1, Dim4, k, l, n, o, i, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, k, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim2, Dim0, Dim5, Dim6, l, k, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim0, Dim5, Dim6, i, j, k, l, m, l, k, i, n, o, Dim1, Dim4, Dim3, Dim2, Dim0, Dim5, Dim6, j, m, l, k, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, l, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim3, Dim0, Dim5, Dim6, k, l, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim0, Dim5, Dim6, i, j, k, l, m, k, l, i, n, o, Dim1, Dim4, Dim2, Dim3, Dim0, Dim5, Dim6, j, m, k, l, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, i, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim0, Dim2, Dim5, Dim6, l, i, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim2, Dim5, Dim6, i, j, k, l, m, l, i, k, n, o, Dim1, Dim4, Dim3, Dim0, Dim2, Dim5, Dim6, j, m, l, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, i, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim0, Dim3, Dim5, Dim6, k, i, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim3, Dim5, Dim6, i, j, k, l, m, k, i, l, n, o, Dim1, Dim4, Dim2, Dim0, Dim3, Dim5, Dim6, j, m, k, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, l, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim3, Dim2, Dim5, Dim6, i, l, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim2, Dim5, Dim6, i, j, k, l, m, i, l, k, n, o, Dim1, Dim4, Dim0, Dim3, Dim2, Dim5, Dim6, j, m, i, l, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, k, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim3, Dim5, Dim6, i, j, k, l, m, i, k, l, n, o, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, Dim6, j, m, i, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, k, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim2, Dim5, Dim0, Dim6, l, k, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, Dim0, Dim6, i, j, k, l, m, l, k, n, i, o, Dim1, Dim4, Dim3, Dim2, Dim5, Dim0, Dim6, j, m, l, k, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, l, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim3, Dim5, Dim0, Dim6, k, l, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, Dim0, Dim6, i, j, k, l, m, k, l, n, i, o, Dim1, Dim4, Dim2, Dim3, Dim5, Dim0, Dim6, j, m, k, l, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, i, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim0, Dim5, Dim2, Dim6, l, i, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, Dim2, Dim6, i, j, k, l, m, l, i, n, k, o, Dim1, Dim4, Dim3, Dim0, Dim5, Dim2, Dim6, j, m, l, i, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, i, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim0, Dim5, Dim3, Dim6, k, i, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, Dim3, Dim6, i, j, k, l, m, k, i, n, l, o, Dim1, Dim4, Dim2, Dim0, Dim5, Dim3, Dim6, j, m, k, i, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, l, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim3, Dim5, Dim2, Dim6, i, l, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, Dim2, Dim6, i, j, k, l, m, i, l, n, k, o, Dim1, Dim4, Dim0, Dim3, Dim5, Dim2, Dim6, j, m, i, l, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, k, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim2, Dim5, Dim3, Dim6, i, k, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, Dim3, Dim6, i, j, k, l, m, i, k, n, l, o, Dim1, Dim4, Dim0, Dim2, Dim5, Dim3, Dim6, j, m, i, k, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, k, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim2, Dim5, Dim6, Dim0, l, k, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, Dim6, Dim0, i, j, k, l, m, l, k, n, o, i, Dim1, Dim4, Dim3, Dim2, Dim5, Dim6, Dim0, j, m, l, k, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, l, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim3, Dim5, Dim6, Dim0, k, l, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, Dim6, Dim0, i, j, k, l, m, k, l, n, o, i, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim0, j, m, k, l, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, i, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim0, Dim5, Dim6, Dim2, l, i, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, Dim6, Dim2, i, j, k, l, m, l, i, n, o, k, Dim1, Dim4, Dim3, Dim0, Dim5, Dim6, Dim2, j, m, l, i, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, i, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim0, Dim5, Dim6, Dim3, k, i, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, Dim6, Dim3, i, j, k, l, m, k, i, n, o, l, Dim1, Dim4, Dim2, Dim0, Dim5, Dim6, Dim3, j, m, k, i, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, l, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim3, Dim5, Dim6, Dim2, i, l, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, Dim6, Dim2, i, j, k, l, m, i, l, n, o, k, Dim1, Dim4, Dim0, Dim3, Dim5, Dim6, Dim2, j, m, i, l, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, k, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim2, Dim5, Dim6, Dim3, i, k, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim3, i, j, k, l, m, i, k, n, o, l, Dim1, Dim4, Dim0, Dim2, Dim5, Dim6, Dim3, j, m, i, k, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim2, Dim0, Dim6, l, n, k, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, Dim0, Dim6, i, j, k, l, m, l, n, k, i, o, Dim1, Dim4, Dim3, Dim5, Dim2, Dim0, Dim6, j, m, l, n, k, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim3, Dim0, Dim6, k, n, l, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, Dim0, Dim6, i, j, k, l, m, k, n, l, i, o, Dim1, Dim4, Dim2, Dim5, Dim3, Dim0, Dim6, j, m, k, n, l, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim0, Dim2, Dim6, l, n, i, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, Dim2, Dim6, i, j, k, l, m, l, n, i, k, o, Dim1, Dim4, Dim3, Dim5, Dim0, Dim2, Dim6, j, m, l, n, i, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim0, Dim3, Dim6, k, n, i, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, Dim3, Dim6, i, j, k, l, m, k, n, i, l, o, Dim1, Dim4, Dim2, Dim5, Dim0, Dim3, Dim6, j, m, k, n, i, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim3, Dim2, Dim6, i, n, l, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, Dim2, Dim6, i, j, k, l, m, i, n, l, k, o, Dim1, Dim4, Dim0, Dim5, Dim3, Dim2, Dim6, j, m, i, n, l, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim2, Dim3, Dim6, i, n, k, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, Dim3, Dim6, i, j, k, l, m, i, n, k, l, o, Dim1, Dim4, Dim0, Dim5, Dim2, Dim3, Dim6, j, m, i, n, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim2, Dim6, Dim0, l, n, k, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, Dim6, Dim0, i, j, k, l, m, l, n, k, o, i, Dim1, Dim4, Dim3, Dim5, Dim2, Dim6, Dim0, j, m, l, n, k, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim3, Dim6, Dim0, k, n, l, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, Dim6, Dim0, i, j, k, l, m, k, n, l, o, i, Dim1, Dim4, Dim2, Dim5, Dim3, Dim6, Dim0, j, m, k, n, l, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim0, Dim6, Dim2, l, n, i, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, Dim6, Dim2, i, j, k, l, m, l, n, i, o, k, Dim1, Dim4, Dim3, Dim5, Dim0, Dim6, Dim2, j, m, l, n, i, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim0, Dim6, Dim3, k, n, i, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, Dim6, Dim3, i, j, k, l, m, k, n, i, o, l, Dim1, Dim4, Dim2, Dim5, Dim0, Dim6, Dim3, j, m, k, n, i, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim3, Dim6, Dim2, i, n, l, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, Dim6, Dim2, i, j, k, l, m, i, n, l, o, k, Dim1, Dim4, Dim0, Dim5, Dim3, Dim6, Dim2, j, m, i, n, l, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim2, Dim6, Dim3, i, n, k, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, Dim6, Dim3, i, j, k, l, m, i, n, k, o, l, Dim1, Dim4, Dim0, Dim5, Dim2, Dim6, Dim3, j, m, i, n, k, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim2, Dim0, l, n, o, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim2, Dim0, i, j, k, l, m, l, n, o, k, i, Dim1, Dim4, Dim3, Dim5, Dim6, Dim2, Dim0, j, m, l, n, o, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim3, Dim0, k, n, o, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim3, Dim0, i, j, k, l, m, k, n, o, l, i, Dim1, Dim4, Dim2, Dim5, Dim6, Dim3, Dim0, j, m, k, n, o, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim0, Dim2, l, n, o, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim0, Dim2, i, j, k, l, m, l, n, o, i, k, Dim1, Dim4, Dim3, Dim5, Dim6, Dim0, Dim2, j, m, l, n, o, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim0, Dim3, k, n, o, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim0, Dim3, i, j, k, l, m, k, n, o, i, l, Dim1, Dim4, Dim2, Dim5, Dim6, Dim0, Dim3, j, m, k, n, o, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim3, Dim2, i, n, o, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim3, Dim2, i, j, k, l, m, i, n, o, l, k, Dim1, Dim4, Dim0, Dim5, Dim6, Dim3, Dim2, j, m, i, n, o, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim2, Dim3, i, n, o, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim2, Dim3, i, j, k, l, m, i, n, o, k, l, Dim1, Dim4, Dim0, Dim5, Dim6, Dim2, Dim3, j, m, i, n, o, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim2, Dim0, Dim6, n, l, k, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim0, Dim6, i, j, k, l, m, n, l, k, i, o, Dim1, Dim4, Dim5, Dim3, Dim2, Dim0, Dim6, j, m, n, l, k, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim3, Dim0, Dim6, n, k, l, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim0, Dim6, i, j, k, l, m, n, k, l, i, o, Dim1, Dim4, Dim5, Dim2, Dim3, Dim0, Dim6, j, m, n, k, l, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim0, Dim2, Dim6, n, l, i, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim2, Dim6, i, j, k, l, m, n, l, i, k, o, Dim1, Dim4, Dim5, Dim3, Dim0, Dim2, Dim6, j, m, n, l, i, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim0, Dim3, Dim6, n, k, i, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim3, Dim6, i, j, k, l, m, n, k, i, l, o, Dim1, Dim4, Dim5, Dim2, Dim0, Dim3, Dim6, j, m, n, k, i, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim3, Dim2, Dim6, n, i, l, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim2, Dim6, i, j, k, l, m, n, i, l, k, o, Dim1, Dim4, Dim5, Dim0, Dim3, Dim2, Dim6, j, m, n, i, l, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim2, Dim3, Dim6, n, i, k, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, i, j, k, l, m, n, i, k, l, o, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, Dim6, j, m, n, i, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim2, Dim6, Dim0, n, l, k, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim6, Dim0, i, j, k, l, m, n, l, k, o, i, Dim1, Dim4, Dim5, Dim3, Dim2, Dim6, Dim0, j, m, n, l, k, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim3, Dim6, Dim0, n, k, l, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim6, Dim0, i, j, k, l, m, n, k, l, o, i, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim0, j, m, n, k, l, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim0, Dim6, Dim2, n, l, i, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim6, Dim2, i, j, k, l, m, n, l, i, o, k, Dim1, Dim4, Dim5, Dim3, Dim0, Dim6, Dim2, j, m, n, l, i, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim0, Dim6, Dim3, n, k, i, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, Dim3, i, j, k, l, m, n, k, i, o, l, Dim1, Dim4, Dim5, Dim2, Dim0, Dim6, Dim3, j, m, n, k, i, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim3, Dim6, Dim2, n, i, l, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim6, Dim2, i, j, k, l, m, n, i, l, o, k, Dim1, Dim4, Dim5, Dim0, Dim3, Dim6, Dim2, j, m, n, i, l, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim2, Dim6, Dim3, n, i, k, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim3, i, j, k, l, m, n, i, k, o, l, Dim1, Dim4, Dim5, Dim0, Dim2, Dim6, Dim3, j, m, n, i, k, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim2, Dim0, n, l, o, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim2, Dim0, i, j, k, l, m, n, l, o, k, i, Dim1, Dim4, Dim5, Dim3, Dim6, Dim2, Dim0, j, m, n, l, o, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim3, Dim0, n, k, o, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim3, Dim0, i, j, k, l, m, n, k, o, l, i, Dim1, Dim4, Dim5, Dim2, Dim6, Dim3, Dim0, j, m, n, k, o, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim0, Dim2, n, l, o, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim0, Dim2, i, j, k, l, m, n, l, o, i, k, Dim1, Dim4, Dim5, Dim3, Dim6, Dim0, Dim2, j, m, n, l, o, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim0, Dim3, n, k, o, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, Dim3, i, j, k, l, m, n, k, o, i, l, Dim1, Dim4, Dim5, Dim2, Dim6, Dim0, Dim3, j, m, n, k, o, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim3, Dim2, n, i, o, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim3, Dim2, i, j, k, l, m, n, i, o, l, k, Dim1, Dim4, Dim5, Dim0, Dim6, Dim3, Dim2, j, m, n, i, o, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim2, Dim3, n, i, o, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, Dim3, i, j, k, l, m, n, i, o, k, l, Dim1, Dim4, Dim5, Dim0, Dim6, Dim2, Dim3, j, m, n, i, o, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim2, Dim0, n, o, l, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim0, i, j, k, l, m, n, o, l, k, i, Dim1, Dim4, Dim5, Dim6, Dim3, Dim2, Dim0, j, m, n, o, l, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim3, Dim0, n, o, k, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim0, i, j, k, l, m, n, o, k, l, i, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim0, j, m, n, o, k, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim0, Dim2, n, o, l, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim2, i, j, k, l, m, n, o, l, i, k, Dim1, Dim4, Dim5, Dim6, Dim3, Dim0, Dim2, j, m, n, o, l, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim0, Dim3, n, o, k, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim3, i, j, k, l, m, n, o, k, i, l, Dim1, Dim4, Dim5, Dim6, Dim2, Dim0, Dim3, j, m, n, o, k, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim3, Dim2, n, o, i, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim2, i, j, k, l, m, n, o, i, l, k, Dim1, Dim4, Dim5, Dim6, Dim0, Dim3, Dim2, j, m, n, o, i, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim2, Dim3, n, o, i, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim3, i, j, k, l, m, n, o, i, k, l, Dim1, Dim4, Dim5, Dim6, Dim0, Dim2, Dim3, j, m, n, o, i, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, k, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim2, Dim0, Dim5, Dim6, m, k, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim0, Dim5, Dim6, i, j, k, l, m, m, k, i, n, o, Dim1, Dim3, Dim4, Dim2, Dim0, Dim5, Dim6, j, l, m, k, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, m, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim4, Dim0, Dim5, Dim6, k, m, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim0, Dim5, Dim6, i, j, k, l, m, k, m, i, n, o, Dim1, Dim3, Dim2, Dim4, Dim0, Dim5, Dim6, j, l, k, m, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, i, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim0, Dim2, Dim5, Dim6, m, i, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim2, Dim5, Dim6, i, j, k, l, m, m, i, k, n, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, i, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim0, Dim4, Dim5, Dim6, k, i, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim4, Dim5, Dim6, i, j, k, l, m, k, i, m, n, o, Dim1, Dim3, Dim2, Dim0, Dim4, Dim5, Dim6, j, l, k, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, m, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim4, Dim2, Dim5, Dim6, i, m, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim2, Dim5, Dim6, i, j, k, l, m, i, m, k, n, o, Dim1, Dim3, Dim0, Dim4, Dim2, Dim5, Dim6, j, l, i, m, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, k, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim4, Dim5, Dim6, i, j, k, l, m, i, k, m, n, o, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, Dim6, j, l, i, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, k, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim2, Dim5, Dim0, Dim6, m, k, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim5, Dim0, Dim6, i, j, k, l, m, m, k, n, i, o, Dim1, Dim3, Dim4, Dim2, Dim5, Dim0, Dim6, j, l, m, k, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, m, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim4, Dim5, Dim0, Dim6, k, m, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim5, Dim0, Dim6, i, j, k, l, m, k, m, n, i, o, Dim1, Dim3, Dim2, Dim4, Dim5, Dim0, Dim6, j, l, k, m, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, i, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim0, Dim5, Dim2, Dim6, m, i, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim5, Dim2, Dim6, i, j, k, l, m, m, i, n, k, o, Dim1, Dim3, Dim4, Dim0, Dim5, Dim2, Dim6, j, l, m, i, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, i, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim0, Dim5, Dim4, Dim6, k, i, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, Dim4, Dim6, i, j, k, l, m, k, i, n, m, o, Dim1, Dim3, Dim2, Dim0, Dim5, Dim4, Dim6, j, l, k, i, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, m, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim4, Dim5, Dim2, Dim6, i, m, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim5, Dim2, Dim6, i, j, k, l, m, i, m, n, k, o, Dim1, Dim3, Dim0, Dim4, Dim5, Dim2, Dim6, j, l, i, m, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, k, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim2, Dim5, Dim4, Dim6, i, k, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, Dim4, Dim6, i, j, k, l, m, i, k, n, m, o, Dim1, Dim3, Dim0, Dim2, Dim5, Dim4, Dim6, j, l, i, k, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, k, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim2, Dim5, Dim6, Dim0, m, k, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim5, Dim6, Dim0, i, j, k, l, m, m, k, n, o, i, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim0, j, l, m, k, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, m, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim4, Dim5, Dim6, Dim0, k, m, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim5, Dim6, Dim0, i, j, k, l, m, k, m, n, o, i, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim0, j, l, k, m, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, i, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim0, Dim5, Dim6, Dim2, m, i, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim5, Dim6, Dim2, i, j, k, l, m, m, i, n, o, k, Dim1, Dim3, Dim4, Dim0, Dim5, Dim6, Dim2, j, l, m, i, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, i, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim0, Dim5, Dim6, Dim4, k, i, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, Dim6, Dim4, i, j, k, l, m, k, i, n, o, m, Dim1, Dim3, Dim2, Dim0, Dim5, Dim6, Dim4, j, l, k, i, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, m, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim4, Dim5, Dim6, Dim2, i, m, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim5, Dim6, Dim2, i, j, k, l, m, i, m, n, o, k, Dim1, Dim3, Dim0, Dim4, Dim5, Dim6, Dim2, j, l, i, m, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, k, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim2, Dim5, Dim6, Dim4, i, k, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim4, i, j, k, l, m, i, k, n, o, m, Dim1, Dim3, Dim0, Dim2, Dim5, Dim6, Dim4, j, l, i, k, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim2, Dim0, Dim6, m, n, k, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim2, Dim0, Dim6, i, j, k, l, m, m, n, k, i, o, Dim1, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, j, l, m, n, k, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim4, Dim0, Dim6, k, n, m, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim4, Dim0, Dim6, i, j, k, l, m, k, n, m, i, o, Dim1, Dim3, Dim2, Dim5, Dim4, Dim0, Dim6, j, l, k, n, m, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim0, Dim2, Dim6, m, n, i, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim0, Dim2, Dim6, i, j, k, l, m, m, n, i, k, o, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, j, l, m, n, i, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim0, Dim4, Dim6, k, n, i, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, Dim4, Dim6, i, j, k, l, m, k, n, i, m, o, Dim1, Dim3, Dim2, Dim5, Dim0, Dim4, Dim6, j, l, k, n, i, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim4, Dim2, Dim6, i, n, m, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim4, Dim2, Dim6, i, j, k, l, m, i, n, m, k, o, Dim1, Dim3, Dim0, Dim5, Dim4, Dim2, Dim6, j, l, i, n, m, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim2, Dim4, Dim6, i, n, k, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, Dim4, Dim6, i, j, k, l, m, i, n, k, m, o, Dim1, Dim3, Dim0, Dim5, Dim2, Dim4, Dim6, j, l, i, n, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim2, Dim6, Dim0, m, n, k, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim2, Dim6, Dim0, i, j, k, l, m, m, n, k, o, i, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, j, l, m, n, k, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim4, Dim6, Dim0, k, n, m, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim4, Dim6, Dim0, i, j, k, l, m, k, n, m, o, i, Dim1, Dim3, Dim2, Dim5, Dim4, Dim6, Dim0, j, l, k, n, m, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim0, Dim6, Dim2, m, n, i, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim0, Dim6, Dim2, i, j, k, l, m, m, n, i, o, k, Dim1, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, j, l, m, n, i, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim0, Dim6, Dim4, k, n, i, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, Dim6, Dim4, i, j, k, l, m, k, n, i, o, m, Dim1, Dim3, Dim2, Dim5, Dim0, Dim6, Dim4, j, l, k, n, i, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim4, Dim6, Dim2, i, n, m, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim4, Dim6, Dim2, i, j, k, l, m, i, n, m, o, k, Dim1, Dim3, Dim0, Dim5, Dim4, Dim6, Dim2, j, l, i, n, m, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim2, Dim6, Dim4, i, n, k, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, Dim6, Dim4, i, j, k, l, m, i, n, k, o, m, Dim1, Dim3, Dim0, Dim5, Dim2, Dim6, Dim4, j, l, i, n, k, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim2, Dim0, m, n, o, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim2, Dim0, i, j, k, l, m, m, n, o, k, i, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, j, l, m, n, o, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim4, Dim0, k, n, o, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim4, Dim0, i, j, k, l, m, k, n, o, m, i, Dim1, Dim3, Dim2, Dim5, Dim6, Dim4, Dim0, j, l, k, n, o, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim0, Dim2, m, n, o, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim0, Dim2, i, j, k, l, m, m, n, o, i, k, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, j, l, m, n, o, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim0, Dim4, k, n, o, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim0, Dim4, i, j, k, l, m, k, n, o, i, m, Dim1, Dim3, Dim2, Dim5, Dim6, Dim0, Dim4, j, l, k, n, o, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim4, Dim2, i, n, o, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim4, Dim2, i, j, k, l, m, i, n, o, m, k, Dim1, Dim3, Dim0, Dim5, Dim6, Dim4, Dim2, j, l, i, n, o, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim2, Dim4, i, n, o, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim2, Dim4, i, j, k, l, m, i, n, o, k, m, Dim1, Dim3, Dim0, Dim5, Dim6, Dim2, Dim4, j, l, i, n, o, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim2, Dim0, Dim6, n, m, k, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim0, Dim6, i, j, k, l, m, n, m, k, i, o, Dim1, Dim3, Dim5, Dim4, Dim2, Dim0, Dim6, j, l, n, m, k, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim4, Dim0, Dim6, n, k, m, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim0, Dim6, i, j, k, l, m, n, k, m, i, o, Dim1, Dim3, Dim5, Dim2, Dim4, Dim0, Dim6, j, l, n, k, m, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim0, Dim2, Dim6, n, m, i, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim2, Dim6, i, j, k, l, m, n, m, i, k, o, Dim1, Dim3, Dim5, Dim4, Dim0, Dim2, Dim6, j, l, n, m, i, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim0, Dim4, Dim6, n, k, i, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim4, Dim6, i, j, k, l, m, n, k, i, m, o, Dim1, Dim3, Dim5, Dim2, Dim0, Dim4, Dim6, j, l, n, k, i, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim4, Dim2, Dim6, n, i, m, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim2, Dim6, i, j, k, l, m, n, i, m, k, o, Dim1, Dim3, Dim5, Dim0, Dim4, Dim2, Dim6, j, l, n, i, m, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim2, Dim4, Dim6, n, i, k, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim4, Dim6, i, j, k, l, m, n, i, k, m, o, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, Dim6, j, l, n, i, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim2, Dim6, Dim0, n, m, k, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim6, Dim0, i, j, k, l, m, n, m, k, o, i, Dim1, Dim3, Dim5, Dim4, Dim2, Dim6, Dim0, j, l, n, m, k, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim4, Dim6, Dim0, n, k, m, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim6, Dim0, i, j, k, l, m, n, k, m, o, i, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim0, j, l, n, k, m, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim0, Dim6, Dim2, n, m, i, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim6, Dim2, i, j, k, l, m, n, m, i, o, k, Dim1, Dim3, Dim5, Dim4, Dim0, Dim6, Dim2, j, l, n, m, i, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim0, Dim6, Dim4, n, k, i, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, Dim4, i, j, k, l, m, n, k, i, o, m, Dim1, Dim3, Dim5, Dim2, Dim0, Dim6, Dim4, j, l, n, k, i, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim4, Dim6, Dim2, n, i, m, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim6, Dim2, i, j, k, l, m, n, i, m, o, k, Dim1, Dim3, Dim5, Dim0, Dim4, Dim6, Dim2, j, l, n, i, m, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim2, Dim6, Dim4, n, i, k, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim4, i, j, k, l, m, n, i, k, o, m, Dim1, Dim3, Dim5, Dim0, Dim2, Dim6, Dim4, j, l, n, i, k, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim2, Dim0, n, m, o, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim2, Dim0, i, j, k, l, m, n, m, o, k, i, Dim1, Dim3, Dim5, Dim4, Dim6, Dim2, Dim0, j, l, n, m, o, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim4, Dim0, n, k, o, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim4, Dim0, i, j, k, l, m, n, k, o, m, i, Dim1, Dim3, Dim5, Dim2, Dim6, Dim4, Dim0, j, l, n, k, o, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim0, Dim2, n, m, o, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim0, Dim2, i, j, k, l, m, n, m, o, i, k, Dim1, Dim3, Dim5, Dim4, Dim6, Dim0, Dim2, j, l, n, m, o, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim0, Dim4, n, k, o, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, Dim4, i, j, k, l, m, n, k, o, i, m, Dim1, Dim3, Dim5, Dim2, Dim6, Dim0, Dim4, j, l, n, k, o, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim4, Dim2, n, i, o, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim4, Dim2, i, j, k, l, m, n, i, o, m, k, Dim1, Dim3, Dim5, Dim0, Dim6, Dim4, Dim2, j, l, n, i, o, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim2, Dim4, n, i, o, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, Dim4, i, j, k, l, m, n, i, o, k, m, Dim1, Dim3, Dim5, Dim0, Dim6, Dim2, Dim4, j, l, n, i, o, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim2, Dim0, n, o, m, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim0, i, j, k, l, m, n, o, m, k, i, Dim1, Dim3, Dim5, Dim6, Dim4, Dim2, Dim0, j, l, n, o, m, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim4, Dim0, n, o, k, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim0, i, j, k, l, m, n, o, k, m, i, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim0, j, l, n, o, k, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim0, Dim2, n, o, m, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim2, i, j, k, l, m, n, o, m, i, k, Dim1, Dim3, Dim5, Dim6, Dim4, Dim0, Dim2, j, l, n, o, m, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim0, Dim4, n, o, k, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim4, i, j, k, l, m, n, o, k, i, m, Dim1, Dim3, Dim5, Dim6, Dim2, Dim0, Dim4, j, l, n, o, k, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim4, Dim2, n, o, i, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim2, i, j, k, l, m, n, o, i, m, k, Dim1, Dim3, Dim5, Dim6, Dim0, Dim4, Dim2, j, l, n, o, i, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim2, Dim4, n, o, i, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim4, i, j, k, l, m, n, o, i, k, m, Dim1, Dim3, Dim5, Dim6, Dim0, Dim2, Dim4, j, l, n, o, i, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, l, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim3, Dim0, Dim5, Dim6, m, l, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim0, Dim5, Dim6, i, j, k, l, m, m, l, i, n, o, Dim1, Dim2, Dim4, Dim3, Dim0, Dim5, Dim6, j, k, m, l, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, m, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim4, Dim0, Dim5, Dim6, l, m, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim0, Dim5, Dim6, i, j, k, l, m, l, m, i, n, o, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, j, k, l, m, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, i, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim0, Dim3, Dim5, Dim6, m, i, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim3, Dim5, Dim6, i, j, k, l, m, m, i, l, n, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, i, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim0, Dim4, Dim5, Dim6, l, i, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim4, Dim5, Dim6, i, j, k, l, m, l, i, m, n, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, m, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim4, Dim3, Dim5, Dim6, i, m, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim3, Dim5, Dim6, i, j, k, l, m, i, m, l, n, o, Dim1, Dim2, Dim0, Dim4, Dim3, Dim5, Dim6, j, k, i, m, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, i, l, m, n, o, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, Dim6, j, k, i, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, l, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim3, Dim5, Dim0, Dim6, m, l, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim5, Dim0, Dim6, i, j, k, l, m, m, l, n, i, o, Dim1, Dim2, Dim4, Dim3, Dim5, Dim0, Dim6, j, k, m, l, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, m, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim4, Dim5, Dim0, Dim6, l, m, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim5, Dim0, Dim6, i, j, k, l, m, l, m, n, i, o, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, j, k, l, m, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, i, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim0, Dim5, Dim3, Dim6, m, i, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim5, Dim3, Dim6, i, j, k, l, m, m, i, n, l, o, Dim1, Dim2, Dim4, Dim0, Dim5, Dim3, Dim6, j, k, m, i, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, i, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim0, Dim5, Dim4, Dim6, l, i, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, Dim4, Dim6, i, j, k, l, m, l, i, n, m, o, Dim1, Dim2, Dim3, Dim0, Dim5, Dim4, Dim6, j, k, l, i, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, m, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim4, Dim5, Dim3, Dim6, i, m, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim5, Dim3, Dim6, i, j, k, l, m, i, m, n, l, o, Dim1, Dim2, Dim0, Dim4, Dim5, Dim3, Dim6, j, k, i, m, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, l, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim3, Dim5, Dim4, Dim6, i, l, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, Dim4, Dim6, i, j, k, l, m, i, l, n, m, o, Dim1, Dim2, Dim0, Dim3, Dim5, Dim4, Dim6, j, k, i, l, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, l, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim3, Dim5, Dim6, Dim0, m, l, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim5, Dim6, Dim0, i, j, k, l, m, m, l, n, o, i, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim0, j, k, m, l, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, m, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim4, Dim5, Dim6, Dim0, l, m, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim5, Dim6, Dim0, i, j, k, l, m, l, m, n, o, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, j, k, l, m, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, i, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim0, Dim5, Dim6, Dim3, m, i, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim5, Dim6, Dim3, i, j, k, l, m, m, i, n, o, l, Dim1, Dim2, Dim4, Dim0, Dim5, Dim6, Dim3, j, k, m, i, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, i, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim0, Dim5, Dim6, Dim4, l, i, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, Dim6, Dim4, i, j, k, l, m, l, i, n, o, m, Dim1, Dim2, Dim3, Dim0, Dim5, Dim6, Dim4, j, k, l, i, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, m, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim4, Dim5, Dim6, Dim3, i, m, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim5, Dim6, Dim3, i, j, k, l, m, i, m, n, o, l, Dim1, Dim2, Dim0, Dim4, Dim5, Dim6, Dim3, j, k, i, m, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, l, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim3, Dim5, Dim6, Dim4, i, l, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, Dim6, Dim4, i, j, k, l, m, i, l, n, o, m, Dim1, Dim2, Dim0, Dim3, Dim5, Dim6, Dim4, j, k, i, l, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim3, Dim0, Dim6, m, n, l, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim3, Dim0, Dim6, i, j, k, l, m, m, n, l, i, o, Dim1, Dim2, Dim4, Dim5, Dim3, Dim0, Dim6, j, k, m, n, l, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim4, Dim0, Dim6, l, n, m, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim4, Dim0, Dim6, i, j, k, l, m, l, n, m, i, o, Dim1, Dim2, Dim3, Dim5, Dim4, Dim0, Dim6, j, k, l, n, m, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim0, Dim3, Dim6, m, n, i, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim0, Dim3, Dim6, i, j, k, l, m, m, n, i, l, o, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, j, k, m, n, i, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim0, Dim4, Dim6, l, n, i, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, Dim4, Dim6, i, j, k, l, m, l, n, i, m, o, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, j, k, l, n, i, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim4, Dim3, Dim6, i, n, m, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim4, Dim3, Dim6, i, j, k, l, m, i, n, m, l, o, Dim1, Dim2, Dim0, Dim5, Dim4, Dim3, Dim6, j, k, i, n, m, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim3, Dim4, Dim6, i, n, l, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, Dim4, Dim6, i, j, k, l, m, i, n, l, m, o, Dim1, Dim2, Dim0, Dim5, Dim3, Dim4, Dim6, j, k, i, n, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim3, Dim6, Dim0, m, n, l, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim3, Dim6, Dim0, i, j, k, l, m, m, n, l, o, i, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim0, j, k, m, n, l, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim4, Dim6, Dim0, l, n, m, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim4, Dim6, Dim0, i, j, k, l, m, l, n, m, o, i, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim0, j, k, l, n, m, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim0, Dim6, Dim3, m, n, i, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim0, Dim6, Dim3, i, j, k, l, m, m, n, i, o, l, Dim1, Dim2, Dim4, Dim5, Dim0, Dim6, Dim3, j, k, m, n, i, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim0, Dim6, Dim4, l, n, i, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, Dim6, Dim4, i, j, k, l, m, l, n, i, o, m, Dim1, Dim2, Dim3, Dim5, Dim0, Dim6, Dim4, j, k, l, n, i, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim4, Dim6, Dim3, i, n, m, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim4, Dim6, Dim3, i, j, k, l, m, i, n, m, o, l, Dim1, Dim2, Dim0, Dim5, Dim4, Dim6, Dim3, j, k, i, n, m, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim3, Dim6, Dim4, i, n, l, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, Dim6, Dim4, i, j, k, l, m, i, n, l, o, m, Dim1, Dim2, Dim0, Dim5, Dim3, Dim6, Dim4, j, k, i, n, l, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim3, Dim0, m, n, o, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim3, Dim0, i, j, k, l, m, m, n, o, l, i, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, Dim0, j, k, m, n, o, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim4, Dim0, l, n, o, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim4, Dim0, i, j, k, l, m, l, n, o, m, i, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, Dim0, j, k, l, n, o, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim0, Dim3, m, n, o, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim0, Dim3, i, j, k, l, m, m, n, o, i, l, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, j, k, m, n, o, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim0, Dim4, l, n, o, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim0, Dim4, i, j, k, l, m, l, n, o, i, m, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, j, k, l, n, o, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim4, Dim3, i, n, o, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim4, Dim3, i, j, k, l, m, i, n, o, m, l, Dim1, Dim2, Dim0, Dim5, Dim6, Dim4, Dim3, j, k, i, n, o, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim3, Dim4, i, n, o, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim3, Dim4, i, j, k, l, m, i, n, o, l, m, Dim1, Dim2, Dim0, Dim5, Dim6, Dim3, Dim4, j, k, i, n, o, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim3, Dim0, Dim6, n, m, l, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim0, Dim6, i, j, k, l, m, n, m, l, i, o, Dim1, Dim2, Dim5, Dim4, Dim3, Dim0, Dim6, j, k, n, m, l, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim4, Dim0, Dim6, n, l, m, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim0, Dim6, i, j, k, l, m, n, l, m, i, o, Dim1, Dim2, Dim5, Dim3, Dim4, Dim0, Dim6, j, k, n, l, m, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim0, Dim3, Dim6, n, m, i, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim3, Dim6, i, j, k, l, m, n, m, i, l, o, Dim1, Dim2, Dim5, Dim4, Dim0, Dim3, Dim6, j, k, n, m, i, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim0, Dim4, Dim6, n, l, i, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim4, Dim6, i, j, k, l, m, n, l, i, m, o, Dim1, Dim2, Dim5, Dim3, Dim0, Dim4, Dim6, j, k, n, l, i, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim4, Dim3, Dim6, n, i, m, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim3, Dim6, i, j, k, l, m, n, i, m, l, o, Dim1, Dim2, Dim5, Dim0, Dim4, Dim3, Dim6, j, k, n, i, m, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim3, Dim4, Dim6, n, i, l, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim4, Dim6, i, j, k, l, m, n, i, l, m, o, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, Dim6, j, k, n, i, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim3, Dim6, Dim0, n, m, l, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim6, Dim0, i, j, k, l, m, n, m, l, o, i, Dim1, Dim2, Dim5, Dim4, Dim3, Dim6, Dim0, j, k, n, m, l, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim4, Dim6, Dim0, n, l, m, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim6, Dim0, i, j, k, l, m, n, l, m, o, i, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim0, j, k, n, l, m, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim0, Dim6, Dim3, n, m, i, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim6, Dim3, i, j, k, l, m, n, m, i, o, l, Dim1, Dim2, Dim5, Dim4, Dim0, Dim6, Dim3, j, k, n, m, i, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim0, Dim6, Dim4, n, l, i, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim6, Dim4, i, j, k, l, m, n, l, i, o, m, Dim1, Dim2, Dim5, Dim3, Dim0, Dim6, Dim4, j, k, n, l, i, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim4, Dim6, Dim3, n, i, m, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim6, Dim3, i, j, k, l, m, n, i, m, o, l, Dim1, Dim2, Dim5, Dim0, Dim4, Dim6, Dim3, j, k, n, i, m, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim3, Dim6, Dim4, n, i, l, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim6, Dim4, i, j, k, l, m, n, i, l, o, m, Dim1, Dim2, Dim5, Dim0, Dim3, Dim6, Dim4, j, k, n, i, l, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim3, Dim0, n, m, o, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim3, Dim0, i, j, k, l, m, n, m, o, l, i, Dim1, Dim2, Dim5, Dim4, Dim6, Dim3, Dim0, j, k, n, m, o, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim4, Dim0, n, l, o, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim4, Dim0, i, j, k, l, m, n, l, o, m, i, Dim1, Dim2, Dim5, Dim3, Dim6, Dim4, Dim0, j, k, n, l, o, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim0, Dim3, n, m, o, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim0, Dim3, i, j, k, l, m, n, m, o, i, l, Dim1, Dim2, Dim5, Dim4, Dim6, Dim0, Dim3, j, k, n, m, o, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim0, Dim4, n, l, o, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim0, Dim4, i, j, k, l, m, n, l, o, i, m, Dim1, Dim2, Dim5, Dim3, Dim6, Dim0, Dim4, j, k, n, l, o, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim4, Dim3, n, i, o, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim4, Dim3, i, j, k, l, m, n, i, o, m, l, Dim1, Dim2, Dim5, Dim0, Dim6, Dim4, Dim3, j, k, n, i, o, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim3, Dim4, n, i, o, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim3, Dim4, i, j, k, l, m, n, i, o, l, m, Dim1, Dim2, Dim5, Dim0, Dim6, Dim3, Dim4, j, k, n, i, o, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim3, Dim0, n, o, m, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim0, i, j, k, l, m, n, o, m, l, i, Dim1, Dim2, Dim5, Dim6, Dim4, Dim3, Dim0, j, k, n, o, m, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim4, Dim0, n, o, l, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim0, i, j, k, l, m, n, o, l, m, i, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim0, j, k, n, o, l, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim0, Dim3, n, o, m, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim3, i, j, k, l, m, n, o, m, i, l, Dim1, Dim2, Dim5, Dim6, Dim4, Dim0, Dim3, j, k, n, o, m, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim0, Dim4, n, o, l, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim4, i, j, k, l, m, n, o, l, i, m, Dim1, Dim2, Dim5, Dim6, Dim3, Dim0, Dim4, j, k, n, o, l, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim4, Dim3, n, o, i, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim3, i, j, k, l, m, n, o, i, m, l, Dim1, Dim2, Dim5, Dim6, Dim0, Dim4, Dim3, j, k, n, o, i, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim3, Dim4, n, o, i, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim4, i, j, k, l, m, n, o, i, l, m, Dim1, Dim2, Dim5, Dim6, Dim0, Dim3, Dim4, j, k, n, o, i, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, k, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim2, Dim1, Dim5, Dim6, l, k, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim1, Dim5, Dim6, i, j, k, l, m, l, k, j, n, o, Dim0, Dim4, Dim3, Dim2, Dim1, Dim5, Dim6, i, m, l, k, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, l, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim3, Dim1, Dim5, Dim6, k, l, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim1, Dim5, Dim6, i, j, k, l, m, k, l, j, n, o, Dim0, Dim4, Dim2, Dim3, Dim1, Dim5, Dim6, i, m, k, l, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, j, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim1, Dim2, Dim5, Dim6, l, j, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim2, Dim5, Dim6, i, j, k, l, m, l, j, k, n, o, Dim0, Dim4, Dim3, Dim1, Dim2, Dim5, Dim6, i, m, l, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, j, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim1, Dim3, Dim5, Dim6, k, j, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim3, Dim5, Dim6, i, j, k, l, m, k, j, l, n, o, Dim0, Dim4, Dim2, Dim1, Dim3, Dim5, Dim6, i, m, k, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, l, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim3, Dim2, Dim5, Dim6, j, l, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim2, Dim5, Dim6, i, j, k, l, m, j, l, k, n, o, Dim0, Dim4, Dim1, Dim3, Dim2, Dim5, Dim6, i, m, j, l, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, k, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim3, Dim5, Dim6, i, j, k, l, m, j, k, l, n, o, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, Dim6, i, m, j, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, k, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim2, Dim5, Dim1, Dim6, l, k, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, Dim1, Dim6, i, j, k, l, m, l, k, n, j, o, Dim0, Dim4, Dim3, Dim2, Dim5, Dim1, Dim6, i, m, l, k, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, l, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim3, Dim5, Dim1, Dim6, k, l, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, Dim1, Dim6, i, j, k, l, m, k, l, n, j, o, Dim0, Dim4, Dim2, Dim3, Dim5, Dim1, Dim6, i, m, k, l, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, j, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim1, Dim5, Dim2, Dim6, l, j, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, Dim2, Dim6, i, j, k, l, m, l, j, n, k, o, Dim0, Dim4, Dim3, Dim1, Dim5, Dim2, Dim6, i, m, l, j, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, j, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim1, Dim5, Dim3, Dim6, k, j, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, Dim3, Dim6, i, j, k, l, m, k, j, n, l, o, Dim0, Dim4, Dim2, Dim1, Dim5, Dim3, Dim6, i, m, k, j, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, l, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim3, Dim5, Dim2, Dim6, j, l, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, Dim2, Dim6, i, j, k, l, m, j, l, n, k, o, Dim0, Dim4, Dim1, Dim3, Dim5, Dim2, Dim6, i, m, j, l, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, k, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim2, Dim5, Dim3, Dim6, j, k, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, Dim3, Dim6, i, j, k, l, m, j, k, n, l, o, Dim0, Dim4, Dim1, Dim2, Dim5, Dim3, Dim6, i, m, j, k, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, k, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim2, Dim5, Dim6, Dim1, l, k, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, Dim6, Dim1, i, j, k, l, m, l, k, n, o, j, Dim0, Dim4, Dim3, Dim2, Dim5, Dim6, Dim1, i, m, l, k, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, l, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim3, Dim5, Dim6, Dim1, k, l, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, Dim6, Dim1, i, j, k, l, m, k, l, n, o, j, Dim0, Dim4, Dim2, Dim3, Dim5, Dim6, Dim1, i, m, k, l, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, j, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim1, Dim5, Dim6, Dim2, l, j, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, Dim6, Dim2, i, j, k, l, m, l, j, n, o, k, Dim0, Dim4, Dim3, Dim1, Dim5, Dim6, Dim2, i, m, l, j, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, j, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim1, Dim5, Dim6, Dim3, k, j, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, Dim6, Dim3, i, j, k, l, m, k, j, n, o, l, Dim0, Dim4, Dim2, Dim1, Dim5, Dim6, Dim3, i, m, k, j, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, l, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim3, Dim5, Dim6, Dim2, j, l, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, Dim6, Dim2, i, j, k, l, m, j, l, n, o, k, Dim0, Dim4, Dim1, Dim3, Dim5, Dim6, Dim2, i, m, j, l, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, k, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim2, Dim5, Dim6, Dim3, j, k, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim3, i, j, k, l, m, j, k, n, o, l, Dim0, Dim4, Dim1, Dim2, Dim5, Dim6, Dim3, i, m, j, k, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim2, Dim1, Dim6, l, n, k, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, Dim1, Dim6, i, j, k, l, m, l, n, k, j, o, Dim0, Dim4, Dim3, Dim5, Dim2, Dim1, Dim6, i, m, l, n, k, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim3, Dim1, Dim6, k, n, l, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, Dim1, Dim6, i, j, k, l, m, k, n, l, j, o, Dim0, Dim4, Dim2, Dim5, Dim3, Dim1, Dim6, i, m, k, n, l, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim1, Dim2, Dim6, l, n, j, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, Dim2, Dim6, i, j, k, l, m, l, n, j, k, o, Dim0, Dim4, Dim3, Dim5, Dim1, Dim2, Dim6, i, m, l, n, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim1, Dim3, Dim6, k, n, j, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, Dim3, Dim6, i, j, k, l, m, k, n, j, l, o, Dim0, Dim4, Dim2, Dim5, Dim1, Dim3, Dim6, i, m, k, n, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim3, Dim2, Dim6, j, n, l, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, Dim2, Dim6, i, j, k, l, m, j, n, l, k, o, Dim0, Dim4, Dim1, Dim5, Dim3, Dim2, Dim6, i, m, j, n, l, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim2, Dim3, Dim6, j, n, k, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, Dim3, Dim6, i, j, k, l, m, j, n, k, l, o, Dim0, Dim4, Dim1, Dim5, Dim2, Dim3, Dim6, i, m, j, n, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim2, Dim6, Dim1, l, n, k, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, Dim6, Dim1, i, j, k, l, m, l, n, k, o, j, Dim0, Dim4, Dim3, Dim5, Dim2, Dim6, Dim1, i, m, l, n, k, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim3, Dim6, Dim1, k, n, l, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, Dim6, Dim1, i, j, k, l, m, k, n, l, o, j, Dim0, Dim4, Dim2, Dim5, Dim3, Dim6, Dim1, i, m, k, n, l, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim1, Dim6, Dim2, l, n, j, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, Dim6, Dim2, i, j, k, l, m, l, n, j, o, k, Dim0, Dim4, Dim3, Dim5, Dim1, Dim6, Dim2, i, m, l, n, j, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim1, Dim6, Dim3, k, n, j, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, Dim6, Dim3, i, j, k, l, m, k, n, j, o, l, Dim0, Dim4, Dim2, Dim5, Dim1, Dim6, Dim3, i, m, k, n, j, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim3, Dim6, Dim2, j, n, l, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, Dim6, Dim2, i, j, k, l, m, j, n, l, o, k, Dim0, Dim4, Dim1, Dim5, Dim3, Dim6, Dim2, i, m, j, n, l, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim2, Dim6, Dim3, j, n, k, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, Dim6, Dim3, i, j, k, l, m, j, n, k, o, l, Dim0, Dim4, Dim1, Dim5, Dim2, Dim6, Dim3, i, m, j, n, k, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim2, Dim1, l, n, o, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim2, Dim1, i, j, k, l, m, l, n, o, k, j, Dim0, Dim4, Dim3, Dim5, Dim6, Dim2, Dim1, i, m, l, n, o, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim3, Dim1, k, n, o, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim3, Dim1, i, j, k, l, m, k, n, o, l, j, Dim0, Dim4, Dim2, Dim5, Dim6, Dim3, Dim1, i, m, k, n, o, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim1, Dim2, l, n, o, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim1, Dim2, i, j, k, l, m, l, n, o, j, k, Dim0, Dim4, Dim3, Dim5, Dim6, Dim1, Dim2, i, m, l, n, o, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim1, Dim3, k, n, o, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim1, Dim3, i, j, k, l, m, k, n, o, j, l, Dim0, Dim4, Dim2, Dim5, Dim6, Dim1, Dim3, i, m, k, n, o, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim3, Dim2, j, n, o, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim3, Dim2, i, j, k, l, m, j, n, o, l, k, Dim0, Dim4, Dim1, Dim5, Dim6, Dim3, Dim2, i, m, j, n, o, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim2, Dim3, j, n, o, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim2, Dim3, i, j, k, l, m, j, n, o, k, l, Dim0, Dim4, Dim1, Dim5, Dim6, Dim2, Dim3, i, m, j, n, o, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim2, Dim1, Dim6, n, l, k, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim1, Dim6, i, j, k, l, m, n, l, k, j, o, Dim0, Dim4, Dim5, Dim3, Dim2, Dim1, Dim6, i, m, n, l, k, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim3, Dim1, Dim6, n, k, l, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim1, Dim6, i, j, k, l, m, n, k, l, j, o, Dim0, Dim4, Dim5, Dim2, Dim3, Dim1, Dim6, i, m, n, k, l, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim1, Dim2, Dim6, n, l, j, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim2, Dim6, i, j, k, l, m, n, l, j, k, o, Dim0, Dim4, Dim5, Dim3, Dim1, Dim2, Dim6, i, m, n, l, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim1, Dim3, Dim6, n, k, j, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim3, Dim6, i, j, k, l, m, n, k, j, l, o, Dim0, Dim4, Dim5, Dim2, Dim1, Dim3, Dim6, i, m, n, k, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim3, Dim2, Dim6, n, j, l, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim2, Dim6, i, j, k, l, m, n, j, l, k, o, Dim0, Dim4, Dim5, Dim1, Dim3, Dim2, Dim6, i, m, n, j, l, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim2, Dim3, Dim6, n, j, k, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, j, k, l, m, n, j, k, l, o, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, Dim6, i, m, n, j, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim2, Dim6, Dim1, n, l, k, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim6, Dim1, i, j, k, l, m, n, l, k, o, j, Dim0, Dim4, Dim5, Dim3, Dim2, Dim6, Dim1, i, m, n, l, k, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim3, Dim6, Dim1, n, k, l, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim6, Dim1, i, j, k, l, m, n, k, l, o, j, Dim0, Dim4, Dim5, Dim2, Dim3, Dim6, Dim1, i, m, n, k, l, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim1, Dim6, Dim2, n, l, j, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim6, Dim2, i, j, k, l, m, n, l, j, o, k, Dim0, Dim4, Dim5, Dim3, Dim1, Dim6, Dim2, i, m, n, l, j, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim1, Dim6, Dim3, n, k, j, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, Dim3, i, j, k, l, m, n, k, j, o, l, Dim0, Dim4, Dim5, Dim2, Dim1, Dim6, Dim3, i, m, n, k, j, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim3, Dim6, Dim2, n, j, l, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim6, Dim2, i, j, k, l, m, n, j, l, o, k, Dim0, Dim4, Dim5, Dim1, Dim3, Dim6, Dim2, i, m, n, j, l, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim2, Dim6, Dim3, n, j, k, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim3, i, j, k, l, m, n, j, k, o, l, Dim0, Dim4, Dim5, Dim1, Dim2, Dim6, Dim3, i, m, n, j, k, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim2, Dim1, n, l, o, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim2, Dim1, i, j, k, l, m, n, l, o, k, j, Dim0, Dim4, Dim5, Dim3, Dim6, Dim2, Dim1, i, m, n, l, o, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim3, Dim1, n, k, o, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim3, Dim1, i, j, k, l, m, n, k, o, l, j, Dim0, Dim4, Dim5, Dim2, Dim6, Dim3, Dim1, i, m, n, k, o, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim1, Dim2, n, l, o, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim1, Dim2, i, j, k, l, m, n, l, o, j, k, Dim0, Dim4, Dim5, Dim3, Dim6, Dim1, Dim2, i, m, n, l, o, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim1, Dim3, n, k, o, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, Dim3, i, j, k, l, m, n, k, o, j, l, Dim0, Dim4, Dim5, Dim2, Dim6, Dim1, Dim3, i, m, n, k, o, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim3, Dim2, n, j, o, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim3, Dim2, i, j, k, l, m, n, j, o, l, k, Dim0, Dim4, Dim5, Dim1, Dim6, Dim3, Dim2, i, m, n, j, o, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim2, Dim3, n, j, o, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, Dim3, i, j, k, l, m, n, j, o, k, l, Dim0, Dim4, Dim5, Dim1, Dim6, Dim2, Dim3, i, m, n, j, o, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim2, Dim1, n, o, l, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim1, i, j, k, l, m, n, o, l, k, j, Dim0, Dim4, Dim5, Dim6, Dim3, Dim2, Dim1, i, m, n, o, l, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim3, Dim1, n, o, k, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim1, i, j, k, l, m, n, o, k, l, j, Dim0, Dim4, Dim5, Dim6, Dim2, Dim3, Dim1, i, m, n, o, k, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim1, Dim2, n, o, l, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim2, i, j, k, l, m, n, o, l, j, k, Dim0, Dim4, Dim5, Dim6, Dim3, Dim1, Dim2, i, m, n, o, l, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim1, Dim3, n, o, k, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim3, i, j, k, l, m, n, o, k, j, l, Dim0, Dim4, Dim5, Dim6, Dim2, Dim1, Dim3, i, m, n, o, k, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim3, Dim2, n, o, j, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim2, i, j, k, l, m, n, o, j, l, k, Dim0, Dim4, Dim5, Dim6, Dim1, Dim3, Dim2, i, m, n, o, j, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim2, Dim3, n, o, j, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim3, i, j, k, l, m, n, o, j, k, l, Dim0, Dim4, Dim5, Dim6, Dim1, Dim2, Dim3, i, m, n, o, j, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, k, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim2, Dim1, Dim5, Dim6, m, k, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim1, Dim5, Dim6, i, j, k, l, m, m, k, j, n, o, Dim0, Dim3, Dim4, Dim2, Dim1, Dim5, Dim6, i, l, m, k, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, m, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim4, Dim1, Dim5, Dim6, k, m, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim1, Dim5, Dim6, i, j, k, l, m, k, m, j, n, o, Dim0, Dim3, Dim2, Dim4, Dim1, Dim5, Dim6, i, l, k, m, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, j, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim1, Dim2, Dim5, Dim6, m, j, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim2, Dim5, Dim6, i, j, k, l, m, m, j, k, n, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, j, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim1, Dim4, Dim5, Dim6, k, j, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim4, Dim5, Dim6, i, j, k, l, m, k, j, m, n, o, Dim0, Dim3, Dim2, Dim1, Dim4, Dim5, Dim6, i, l, k, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, m, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim4, Dim2, Dim5, Dim6, j, m, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim2, Dim5, Dim6, i, j, k, l, m, j, m, k, n, o, Dim0, Dim3, Dim1, Dim4, Dim2, Dim5, Dim6, i, l, j, m, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, k, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim4, Dim5, Dim6, i, j, k, l, m, j, k, m, n, o, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, Dim6, i, l, j, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, k, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim2, Dim5, Dim1, Dim6, m, k, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim5, Dim1, Dim6, i, j, k, l, m, m, k, n, j, o, Dim0, Dim3, Dim4, Dim2, Dim5, Dim1, Dim6, i, l, m, k, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, m, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim4, Dim5, Dim1, Dim6, k, m, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim5, Dim1, Dim6, i, j, k, l, m, k, m, n, j, o, Dim0, Dim3, Dim2, Dim4, Dim5, Dim1, Dim6, i, l, k, m, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, j, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim1, Dim5, Dim2, Dim6, m, j, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim5, Dim2, Dim6, i, j, k, l, m, m, j, n, k, o, Dim0, Dim3, Dim4, Dim1, Dim5, Dim2, Dim6, i, l, m, j, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, j, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim1, Dim5, Dim4, Dim6, k, j, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, Dim4, Dim6, i, j, k, l, m, k, j, n, m, o, Dim0, Dim3, Dim2, Dim1, Dim5, Dim4, Dim6, i, l, k, j, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, m, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim4, Dim5, Dim2, Dim6, j, m, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim5, Dim2, Dim6, i, j, k, l, m, j, m, n, k, o, Dim0, Dim3, Dim1, Dim4, Dim5, Dim2, Dim6, i, l, j, m, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, k, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim2, Dim5, Dim4, Dim6, j, k, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, Dim4, Dim6, i, j, k, l, m, j, k, n, m, o, Dim0, Dim3, Dim1, Dim2, Dim5, Dim4, Dim6, i, l, j, k, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, k, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim2, Dim5, Dim6, Dim1, m, k, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim5, Dim6, Dim1, i, j, k, l, m, m, k, n, o, j, Dim0, Dim3, Dim4, Dim2, Dim5, Dim6, Dim1, i, l, m, k, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, m, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim4, Dim5, Dim6, Dim1, k, m, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim5, Dim6, Dim1, i, j, k, l, m, k, m, n, o, j, Dim0, Dim3, Dim2, Dim4, Dim5, Dim6, Dim1, i, l, k, m, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, j, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim1, Dim5, Dim6, Dim2, m, j, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim5, Dim6, Dim2, i, j, k, l, m, m, j, n, o, k, Dim0, Dim3, Dim4, Dim1, Dim5, Dim6, Dim2, i, l, m, j, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, j, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim1, Dim5, Dim6, Dim4, k, j, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, Dim6, Dim4, i, j, k, l, m, k, j, n, o, m, Dim0, Dim3, Dim2, Dim1, Dim5, Dim6, Dim4, i, l, k, j, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, m, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim4, Dim5, Dim6, Dim2, j, m, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim5, Dim6, Dim2, i, j, k, l, m, j, m, n, o, k, Dim0, Dim3, Dim1, Dim4, Dim5, Dim6, Dim2, i, l, j, m, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, k, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim2, Dim5, Dim6, Dim4, j, k, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim4, i, j, k, l, m, j, k, n, o, m, Dim0, Dim3, Dim1, Dim2, Dim5, Dim6, Dim4, i, l, j, k, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim2, Dim1, Dim6, m, n, k, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim2, Dim1, Dim6, i, j, k, l, m, m, n, k, j, o, Dim0, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, i, l, m, n, k, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim4, Dim1, Dim6, k, n, m, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim4, Dim1, Dim6, i, j, k, l, m, k, n, m, j, o, Dim0, Dim3, Dim2, Dim5, Dim4, Dim1, Dim6, i, l, k, n, m, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim1, Dim2, Dim6, m, n, j, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim1, Dim2, Dim6, i, j, k, l, m, m, n, j, k, o, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, i, l, m, n, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim1, Dim4, Dim6, k, n, j, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, Dim4, Dim6, i, j, k, l, m, k, n, j, m, o, Dim0, Dim3, Dim2, Dim5, Dim1, Dim4, Dim6, i, l, k, n, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim4, Dim2, Dim6, j, n, m, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim4, Dim2, Dim6, i, j, k, l, m, j, n, m, k, o, Dim0, Dim3, Dim1, Dim5, Dim4, Dim2, Dim6, i, l, j, n, m, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim2, Dim4, Dim6, j, n, k, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, Dim4, Dim6, i, j, k, l, m, j, n, k, m, o, Dim0, Dim3, Dim1, Dim5, Dim2, Dim4, Dim6, i, l, j, n, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim2, Dim6, Dim1, m, n, k, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim2, Dim6, Dim1, i, j, k, l, m, m, n, k, o, j, Dim0, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, i, l, m, n, k, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim4, Dim6, Dim1, k, n, m, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim4, Dim6, Dim1, i, j, k, l, m, k, n, m, o, j, Dim0, Dim3, Dim2, Dim5, Dim4, Dim6, Dim1, i, l, k, n, m, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim1, Dim6, Dim2, m, n, j, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim1, Dim6, Dim2, i, j, k, l, m, m, n, j, o, k, Dim0, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, i, l, m, n, j, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim1, Dim6, Dim4, k, n, j, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, Dim6, Dim4, i, j, k, l, m, k, n, j, o, m, Dim0, Dim3, Dim2, Dim5, Dim1, Dim6, Dim4, i, l, k, n, j, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim4, Dim6, Dim2, j, n, m, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim4, Dim6, Dim2, i, j, k, l, m, j, n, m, o, k, Dim0, Dim3, Dim1, Dim5, Dim4, Dim6, Dim2, i, l, j, n, m, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim2, Dim6, Dim4, j, n, k, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, Dim6, Dim4, i, j, k, l, m, j, n, k, o, m, Dim0, Dim3, Dim1, Dim5, Dim2, Dim6, Dim4, i, l, j, n, k, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim2, Dim1, m, n, o, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim2, Dim1, i, j, k, l, m, m, n, o, k, j, Dim0, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, i, l, m, n, o, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim4, Dim1, k, n, o, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim4, Dim1, i, j, k, l, m, k, n, o, m, j, Dim0, Dim3, Dim2, Dim5, Dim6, Dim4, Dim1, i, l, k, n, o, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim1, Dim2, m, n, o, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim1, Dim2, i, j, k, l, m, m, n, o, j, k, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, i, l, m, n, o, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim1, Dim4, k, n, o, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim1, Dim4, i, j, k, l, m, k, n, o, j, m, Dim0, Dim3, Dim2, Dim5, Dim6, Dim1, Dim4, i, l, k, n, o, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim4, Dim2, j, n, o, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim4, Dim2, i, j, k, l, m, j, n, o, m, k, Dim0, Dim3, Dim1, Dim5, Dim6, Dim4, Dim2, i, l, j, n, o, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim2, Dim4, j, n, o, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim2, Dim4, i, j, k, l, m, j, n, o, k, m, Dim0, Dim3, Dim1, Dim5, Dim6, Dim2, Dim4, i, l, j, n, o, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim2, Dim1, Dim6, n, m, k, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim1, Dim6, i, j, k, l, m, n, m, k, j, o, Dim0, Dim3, Dim5, Dim4, Dim2, Dim1, Dim6, i, l, n, m, k, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim4, Dim1, Dim6, n, k, m, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim1, Dim6, i, j, k, l, m, n, k, m, j, o, Dim0, Dim3, Dim5, Dim2, Dim4, Dim1, Dim6, i, l, n, k, m, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim1, Dim2, Dim6, n, m, j, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim2, Dim6, i, j, k, l, m, n, m, j, k, o, Dim0, Dim3, Dim5, Dim4, Dim1, Dim2, Dim6, i, l, n, m, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim1, Dim4, Dim6, n, k, j, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim4, Dim6, i, j, k, l, m, n, k, j, m, o, Dim0, Dim3, Dim5, Dim2, Dim1, Dim4, Dim6, i, l, n, k, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim4, Dim2, Dim6, n, j, m, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim2, Dim6, i, j, k, l, m, n, j, m, k, o, Dim0, Dim3, Dim5, Dim1, Dim4, Dim2, Dim6, i, l, n, j, m, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim2, Dim4, Dim6, n, j, k, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim4, Dim6, i, j, k, l, m, n, j, k, m, o, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, Dim6, i, l, n, j, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim2, Dim6, Dim1, n, m, k, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim6, Dim1, i, j, k, l, m, n, m, k, o, j, Dim0, Dim3, Dim5, Dim4, Dim2, Dim6, Dim1, i, l, n, m, k, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim4, Dim6, Dim1, n, k, m, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim6, Dim1, i, j, k, l, m, n, k, m, o, j, Dim0, Dim3, Dim5, Dim2, Dim4, Dim6, Dim1, i, l, n, k, m, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim1, Dim6, Dim2, n, m, j, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim6, Dim2, i, j, k, l, m, n, m, j, o, k, Dim0, Dim3, Dim5, Dim4, Dim1, Dim6, Dim2, i, l, n, m, j, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim1, Dim6, Dim4, n, k, j, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, Dim4, i, j, k, l, m, n, k, j, o, m, Dim0, Dim3, Dim5, Dim2, Dim1, Dim6, Dim4, i, l, n, k, j, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim4, Dim6, Dim2, n, j, m, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim6, Dim2, i, j, k, l, m, n, j, m, o, k, Dim0, Dim3, Dim5, Dim1, Dim4, Dim6, Dim2, i, l, n, j, m, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim2, Dim6, Dim4, n, j, k, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim4, i, j, k, l, m, n, j, k, o, m, Dim0, Dim3, Dim5, Dim1, Dim2, Dim6, Dim4, i, l, n, j, k, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim2, Dim1, n, m, o, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim2, Dim1, i, j, k, l, m, n, m, o, k, j, Dim0, Dim3, Dim5, Dim4, Dim6, Dim2, Dim1, i, l, n, m, o, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim4, Dim1, n, k, o, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim4, Dim1, i, j, k, l, m, n, k, o, m, j, Dim0, Dim3, Dim5, Dim2, Dim6, Dim4, Dim1, i, l, n, k, o, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim1, Dim2, n, m, o, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim1, Dim2, i, j, k, l, m, n, m, o, j, k, Dim0, Dim3, Dim5, Dim4, Dim6, Dim1, Dim2, i, l, n, m, o, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim1, Dim4, n, k, o, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, Dim4, i, j, k, l, m, n, k, o, j, m, Dim0, Dim3, Dim5, Dim2, Dim6, Dim1, Dim4, i, l, n, k, o, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim4, Dim2, n, j, o, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim4, Dim2, i, j, k, l, m, n, j, o, m, k, Dim0, Dim3, Dim5, Dim1, Dim6, Dim4, Dim2, i, l, n, j, o, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim2, Dim4, n, j, o, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, Dim4, i, j, k, l, m, n, j, o, k, m, Dim0, Dim3, Dim5, Dim1, Dim6, Dim2, Dim4, i, l, n, j, o, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim2, Dim1, n, o, m, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim1, i, j, k, l, m, n, o, m, k, j, Dim0, Dim3, Dim5, Dim6, Dim4, Dim2, Dim1, i, l, n, o, m, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim4, Dim1, n, o, k, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim1, i, j, k, l, m, n, o, k, m, j, Dim0, Dim3, Dim5, Dim6, Dim2, Dim4, Dim1, i, l, n, o, k, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim1, Dim2, n, o, m, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim2, i, j, k, l, m, n, o, m, j, k, Dim0, Dim3, Dim5, Dim6, Dim4, Dim1, Dim2, i, l, n, o, m, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim1, Dim4, n, o, k, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim4, i, j, k, l, m, n, o, k, j, m, Dim0, Dim3, Dim5, Dim6, Dim2, Dim1, Dim4, i, l, n, o, k, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim4, Dim2, n, o, j, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim2, i, j, k, l, m, n, o, j, m, k, Dim0, Dim3, Dim5, Dim6, Dim1, Dim4, Dim2, i, l, n, o, j, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim2, Dim4, n, o, j, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim4, i, j, k, l, m, n, o, j, k, m, Dim0, Dim3, Dim5, Dim6, Dim1, Dim2, Dim4, i, l, n, o, j, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, l, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim3, Dim1, Dim5, Dim6, m, l, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim1, Dim5, Dim6, i, j, k, l, m, m, l, j, n, o, Dim0, Dim2, Dim4, Dim3, Dim1, Dim5, Dim6, i, k, m, l, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, m, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim4, Dim1, Dim5, Dim6, l, m, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim1, Dim5, Dim6, i, j, k, l, m, l, m, j, n, o, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, i, k, l, m, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, j, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim1, Dim3, Dim5, Dim6, m, j, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim3, Dim5, Dim6, i, j, k, l, m, m, j, l, n, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, j, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim1, Dim4, Dim5, Dim6, l, j, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim4, Dim5, Dim6, i, j, k, l, m, l, j, m, n, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, m, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim4, Dim3, Dim5, Dim6, j, m, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim3, Dim5, Dim6, i, j, k, l, m, j, m, l, n, o, Dim0, Dim2, Dim1, Dim4, Dim3, Dim5, Dim6, i, k, j, m, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, j, l, m, n, o, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, Dim6, i, k, j, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, l, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim3, Dim5, Dim1, Dim6, m, l, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim5, Dim1, Dim6, i, j, k, l, m, m, l, n, j, o, Dim0, Dim2, Dim4, Dim3, Dim5, Dim1, Dim6, i, k, m, l, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, m, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim4, Dim5, Dim1, Dim6, l, m, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim5, Dim1, Dim6, i, j, k, l, m, l, m, n, j, o, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, i, k, l, m, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, j, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim1, Dim5, Dim3, Dim6, m, j, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim5, Dim3, Dim6, i, j, k, l, m, m, j, n, l, o, Dim0, Dim2, Dim4, Dim1, Dim5, Dim3, Dim6, i, k, m, j, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, j, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim1, Dim5, Dim4, Dim6, l, j, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, Dim4, Dim6, i, j, k, l, m, l, j, n, m, o, Dim0, Dim2, Dim3, Dim1, Dim5, Dim4, Dim6, i, k, l, j, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, m, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim4, Dim5, Dim3, Dim6, j, m, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim5, Dim3, Dim6, i, j, k, l, m, j, m, n, l, o, Dim0, Dim2, Dim1, Dim4, Dim5, Dim3, Dim6, i, k, j, m, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, l, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim3, Dim5, Dim4, Dim6, j, l, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, Dim4, Dim6, i, j, k, l, m, j, l, n, m, o, Dim0, Dim2, Dim1, Dim3, Dim5, Dim4, Dim6, i, k, j, l, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, l, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim3, Dim5, Dim6, Dim1, m, l, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim5, Dim6, Dim1, i, j, k, l, m, m, l, n, o, j, Dim0, Dim2, Dim4, Dim3, Dim5, Dim6, Dim1, i, k, m, l, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, m, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim4, Dim5, Dim6, Dim1, l, m, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim5, Dim6, Dim1, i, j, k, l, m, l, m, n, o, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, i, k, l, m, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, j, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim1, Dim5, Dim6, Dim3, m, j, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim5, Dim6, Dim3, i, j, k, l, m, m, j, n, o, l, Dim0, Dim2, Dim4, Dim1, Dim5, Dim6, Dim3, i, k, m, j, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, j, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim1, Dim5, Dim6, Dim4, l, j, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, Dim6, Dim4, i, j, k, l, m, l, j, n, o, m, Dim0, Dim2, Dim3, Dim1, Dim5, Dim6, Dim4, i, k, l, j, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, m, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim4, Dim5, Dim6, Dim3, j, m, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim5, Dim6, Dim3, i, j, k, l, m, j, m, n, o, l, Dim0, Dim2, Dim1, Dim4, Dim5, Dim6, Dim3, i, k, j, m, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, l, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim3, Dim5, Dim6, Dim4, j, l, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, Dim6, Dim4, i, j, k, l, m, j, l, n, o, m, Dim0, Dim2, Dim1, Dim3, Dim5, Dim6, Dim4, i, k, j, l, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim3, Dim1, Dim6, m, n, l, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim3, Dim1, Dim6, i, j, k, l, m, m, n, l, j, o, Dim0, Dim2, Dim4, Dim5, Dim3, Dim1, Dim6, i, k, m, n, l, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim4, Dim1, Dim6, l, n, m, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim4, Dim1, Dim6, i, j, k, l, m, l, n, m, j, o, Dim0, Dim2, Dim3, Dim5, Dim4, Dim1, Dim6, i, k, l, n, m, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim1, Dim3, Dim6, m, n, j, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim1, Dim3, Dim6, i, j, k, l, m, m, n, j, l, o, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, i, k, m, n, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim1, Dim4, Dim6, l, n, j, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, Dim4, Dim6, i, j, k, l, m, l, n, j, m, o, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, i, k, l, n, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim4, Dim3, Dim6, j, n, m, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim4, Dim3, Dim6, i, j, k, l, m, j, n, m, l, o, Dim0, Dim2, Dim1, Dim5, Dim4, Dim3, Dim6, i, k, j, n, m, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim3, Dim4, Dim6, j, n, l, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, Dim4, Dim6, i, j, k, l, m, j, n, l, m, o, Dim0, Dim2, Dim1, Dim5, Dim3, Dim4, Dim6, i, k, j, n, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim3, Dim6, Dim1, m, n, l, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim3, Dim6, Dim1, i, j, k, l, m, m, n, l, o, j, Dim0, Dim2, Dim4, Dim5, Dim3, Dim6, Dim1, i, k, m, n, l, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim4, Dim6, Dim1, l, n, m, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim4, Dim6, Dim1, i, j, k, l, m, l, n, m, o, j, Dim0, Dim2, Dim3, Dim5, Dim4, Dim6, Dim1, i, k, l, n, m, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim1, Dim6, Dim3, m, n, j, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim1, Dim6, Dim3, i, j, k, l, m, m, n, j, o, l, Dim0, Dim2, Dim4, Dim5, Dim1, Dim6, Dim3, i, k, m, n, j, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim1, Dim6, Dim4, l, n, j, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, Dim6, Dim4, i, j, k, l, m, l, n, j, o, m, Dim0, Dim2, Dim3, Dim5, Dim1, Dim6, Dim4, i, k, l, n, j, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim4, Dim6, Dim3, j, n, m, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim4, Dim6, Dim3, i, j, k, l, m, j, n, m, o, l, Dim0, Dim2, Dim1, Dim5, Dim4, Dim6, Dim3, i, k, j, n, m, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim3, Dim6, Dim4, j, n, l, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, Dim6, Dim4, i, j, k, l, m, j, n, l, o, m, Dim0, Dim2, Dim1, Dim5, Dim3, Dim6, Dim4, i, k, j, n, l, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim3, Dim1, m, n, o, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim3, Dim1, i, j, k, l, m, m, n, o, l, j, Dim0, Dim2, Dim4, Dim5, Dim6, Dim3, Dim1, i, k, m, n, o, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim4, Dim1, l, n, o, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim4, Dim1, i, j, k, l, m, l, n, o, m, j, Dim0, Dim2, Dim3, Dim5, Dim6, Dim4, Dim1, i, k, l, n, o, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim1, Dim3, m, n, o, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim1, Dim3, i, j, k, l, m, m, n, o, j, l, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, i, k, m, n, o, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim1, Dim4, l, n, o, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim1, Dim4, i, j, k, l, m, l, n, o, j, m, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, i, k, l, n, o, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim4, Dim3, j, n, o, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim4, Dim3, i, j, k, l, m, j, n, o, m, l, Dim0, Dim2, Dim1, Dim5, Dim6, Dim4, Dim3, i, k, j, n, o, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim3, Dim4, j, n, o, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim3, Dim4, i, j, k, l, m, j, n, o, l, m, Dim0, Dim2, Dim1, Dim5, Dim6, Dim3, Dim4, i, k, j, n, o, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim3, Dim1, Dim6, n, m, l, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim1, Dim6, i, j, k, l, m, n, m, l, j, o, Dim0, Dim2, Dim5, Dim4, Dim3, Dim1, Dim6, i, k, n, m, l, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim4, Dim1, Dim6, n, l, m, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim1, Dim6, i, j, k, l, m, n, l, m, j, o, Dim0, Dim2, Dim5, Dim3, Dim4, Dim1, Dim6, i, k, n, l, m, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim1, Dim3, Dim6, n, m, j, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim3, Dim6, i, j, k, l, m, n, m, j, l, o, Dim0, Dim2, Dim5, Dim4, Dim1, Dim3, Dim6, i, k, n, m, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim1, Dim4, Dim6, n, l, j, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim4, Dim6, i, j, k, l, m, n, l, j, m, o, Dim0, Dim2, Dim5, Dim3, Dim1, Dim4, Dim6, i, k, n, l, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim4, Dim3, Dim6, n, j, m, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim3, Dim6, i, j, k, l, m, n, j, m, l, o, Dim0, Dim2, Dim5, Dim1, Dim4, Dim3, Dim6, i, k, n, j, m, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim3, Dim4, Dim6, n, j, l, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim4, Dim6, i, j, k, l, m, n, j, l, m, o, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, Dim6, i, k, n, j, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim3, Dim6, Dim1, n, m, l, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim6, Dim1, i, j, k, l, m, n, m, l, o, j, Dim0, Dim2, Dim5, Dim4, Dim3, Dim6, Dim1, i, k, n, m, l, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim4, Dim6, Dim1, n, l, m, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim6, Dim1, i, j, k, l, m, n, l, m, o, j, Dim0, Dim2, Dim5, Dim3, Dim4, Dim6, Dim1, i, k, n, l, m, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim1, Dim6, Dim3, n, m, j, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim6, Dim3, i, j, k, l, m, n, m, j, o, l, Dim0, Dim2, Dim5, Dim4, Dim1, Dim6, Dim3, i, k, n, m, j, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim1, Dim6, Dim4, n, l, j, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim6, Dim4, i, j, k, l, m, n, l, j, o, m, Dim0, Dim2, Dim5, Dim3, Dim1, Dim6, Dim4, i, k, n, l, j, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim4, Dim6, Dim3, n, j, m, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim6, Dim3, i, j, k, l, m, n, j, m, o, l, Dim0, Dim2, Dim5, Dim1, Dim4, Dim6, Dim3, i, k, n, j, m, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim3, Dim6, Dim4, n, j, l, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim6, Dim4, i, j, k, l, m, n, j, l, o, m, Dim0, Dim2, Dim5, Dim1, Dim3, Dim6, Dim4, i, k, n, j, l, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim3, Dim1, n, m, o, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim3, Dim1, i, j, k, l, m, n, m, o, l, j, Dim0, Dim2, Dim5, Dim4, Dim6, Dim3, Dim1, i, k, n, m, o, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim4, Dim1, n, l, o, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim4, Dim1, i, j, k, l, m, n, l, o, m, j, Dim0, Dim2, Dim5, Dim3, Dim6, Dim4, Dim1, i, k, n, l, o, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim1, Dim3, n, m, o, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim1, Dim3, i, j, k, l, m, n, m, o, j, l, Dim0, Dim2, Dim5, Dim4, Dim6, Dim1, Dim3, i, k, n, m, o, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim1, Dim4, n, l, o, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim1, Dim4, i, j, k, l, m, n, l, o, j, m, Dim0, Dim2, Dim5, Dim3, Dim6, Dim1, Dim4, i, k, n, l, o, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim4, Dim3, n, j, o, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim4, Dim3, i, j, k, l, m, n, j, o, m, l, Dim0, Dim2, Dim5, Dim1, Dim6, Dim4, Dim3, i, k, n, j, o, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim3, Dim4, n, j, o, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim3, Dim4, i, j, k, l, m, n, j, o, l, m, Dim0, Dim2, Dim5, Dim1, Dim6, Dim3, Dim4, i, k, n, j, o, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim3, Dim1, n, o, m, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim1, i, j, k, l, m, n, o, m, l, j, Dim0, Dim2, Dim5, Dim6, Dim4, Dim3, Dim1, i, k, n, o, m, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim4, Dim1, n, o, l, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim1, i, j, k, l, m, n, o, l, m, j, Dim0, Dim2, Dim5, Dim6, Dim3, Dim4, Dim1, i, k, n, o, l, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim1, Dim3, n, o, m, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim3, i, j, k, l, m, n, o, m, j, l, Dim0, Dim2, Dim5, Dim6, Dim4, Dim1, Dim3, i, k, n, o, m, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim1, Dim4, n, o, l, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim4, i, j, k, l, m, n, o, l, j, m, Dim0, Dim2, Dim5, Dim6, Dim3, Dim1, Dim4, i, k, n, o, l, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim4, Dim3, n, o, j, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim3, i, j, k, l, m, n, o, j, m, l, Dim0, Dim2, Dim5, Dim6, Dim1, Dim4, Dim3, i, k, n, o, j, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim3, Dim4, n, o, j, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim4, i, j, k, l, m, n, o, j, l, m, Dim0, Dim2, Dim5, Dim6, Dim1, Dim3, Dim4, i, k, n, o, j, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, l, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim3, Dim2, Dim5, Dim6, m, l, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim2, Dim5, Dim6, i, j, k, l, m, m, l, k, n, o, Dim0, Dim1, Dim4, Dim3, Dim2, Dim5, Dim6, i, j, m, l, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, m, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim4, Dim2, Dim5, Dim6, l, m, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim2, Dim5, Dim6, i, j, k, l, m, l, m, k, n, o, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, i, j, l, m, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, k, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim2, Dim3, Dim5, Dim6, m, k, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, k, l, m, m, k, l, n, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, k, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim2, Dim4, Dim5, Dim6, l, k, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, k, l, m, l, k, m, n, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, m, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim4, Dim3, Dim5, Dim6, k, m, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim3, Dim5, Dim6, i, j, k, l, m, k, m, l, n, o, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, i, j, k, m, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, k, l, m, n, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, l, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim3, Dim5, Dim2, Dim6, m, l, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim5, Dim2, Dim6, i, j, k, l, m, m, l, n, k, o, Dim0, Dim1, Dim4, Dim3, Dim5, Dim2, Dim6, i, j, m, l, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, m, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim4, Dim5, Dim2, Dim6, l, m, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim5, Dim2, Dim6, i, j, k, l, m, l, m, n, k, o, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, i, j, l, m, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, k, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim2, Dim5, Dim3, Dim6, m, k, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim5, Dim3, Dim6, i, j, k, l, m, m, k, n, l, o, Dim0, Dim1, Dim4, Dim2, Dim5, Dim3, Dim6, i, j, m, k, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, k, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim2, Dim5, Dim4, Dim6, l, k, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, Dim4, Dim6, i, j, k, l, m, l, k, n, m, o, Dim0, Dim1, Dim3, Dim2, Dim5, Dim4, Dim6, i, j, l, k, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, m, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim4, Dim5, Dim3, Dim6, k, m, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim5, Dim3, Dim6, i, j, k, l, m, k, m, n, l, o, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, i, j, k, m, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, l, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim3, Dim5, Dim4, Dim6, k, l, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, Dim4, Dim6, i, j, k, l, m, k, l, n, m, o, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, i, j, k, l, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, l, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim3, Dim5, Dim6, Dim2, m, l, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim5, Dim6, Dim2, i, j, k, l, m, m, l, n, o, k, Dim0, Dim1, Dim4, Dim3, Dim5, Dim6, Dim2, i, j, m, l, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, m, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim4, Dim5, Dim6, Dim2, l, m, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim5, Dim6, Dim2, i, j, k, l, m, l, m, n, o, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, i, j, l, m, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, k, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim2, Dim5, Dim6, Dim3, m, k, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim5, Dim6, Dim3, i, j, k, l, m, m, k, n, o, l, Dim0, Dim1, Dim4, Dim2, Dim5, Dim6, Dim3, i, j, m, k, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, k, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim2, Dim5, Dim6, Dim4, l, k, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, Dim6, Dim4, i, j, k, l, m, l, k, n, o, m, Dim0, Dim1, Dim3, Dim2, Dim5, Dim6, Dim4, i, j, l, k, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, m, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim4, Dim5, Dim6, Dim3, k, m, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim5, Dim6, Dim3, i, j, k, l, m, k, m, n, o, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, i, j, k, m, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, l, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim3, Dim5, Dim6, Dim4, k, l, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, Dim6, Dim4, i, j, k, l, m, k, l, n, o, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, i, j, k, l, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim3, Dim2, Dim6, m, n, l, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim3, Dim2, Dim6, i, j, k, l, m, m, n, l, k, o, Dim0, Dim1, Dim4, Dim5, Dim3, Dim2, Dim6, i, j, m, n, l, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim4, Dim2, Dim6, l, n, m, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim4, Dim2, Dim6, i, j, k, l, m, l, n, m, k, o, Dim0, Dim1, Dim3, Dim5, Dim4, Dim2, Dim6, i, j, l, n, m, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim2, Dim3, Dim6, m, n, k, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim2, Dim3, Dim6, i, j, k, l, m, m, n, k, l, o, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, i, j, m, n, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim2, Dim4, Dim6, l, n, k, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, Dim4, Dim6, i, j, k, l, m, l, n, k, m, o, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, i, j, l, n, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim4, Dim3, Dim6, k, n, m, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim4, Dim3, Dim6, i, j, k, l, m, k, n, m, l, o, Dim0, Dim1, Dim2, Dim5, Dim4, Dim3, Dim6, i, j, k, n, m, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim3, Dim4, Dim6, k, n, l, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, Dim4, Dim6, i, j, k, l, m, k, n, l, m, o, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, i, j, k, n, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim3, Dim6, Dim2, m, n, l, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim3, Dim6, Dim2, i, j, k, l, m, m, n, l, o, k, Dim0, Dim1, Dim4, Dim5, Dim3, Dim6, Dim2, i, j, m, n, l, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim4, Dim6, Dim2, l, n, m, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim4, Dim6, Dim2, i, j, k, l, m, l, n, m, o, k, Dim0, Dim1, Dim3, Dim5, Dim4, Dim6, Dim2, i, j, l, n, m, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim2, Dim6, Dim3, m, n, k, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim2, Dim6, Dim3, i, j, k, l, m, m, n, k, o, l, Dim0, Dim1, Dim4, Dim5, Dim2, Dim6, Dim3, i, j, m, n, k, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim2, Dim6, Dim4, l, n, k, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, Dim6, Dim4, i, j, k, l, m, l, n, k, o, m, Dim0, Dim1, Dim3, Dim5, Dim2, Dim6, Dim4, i, j, l, n, k, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim4, Dim6, Dim3, k, n, m, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim4, Dim6, Dim3, i, j, k, l, m, k, n, m, o, l, Dim0, Dim1, Dim2, Dim5, Dim4, Dim6, Dim3, i, j, k, n, m, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim3, Dim6, Dim4, k, n, l, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, Dim6, Dim4, i, j, k, l, m, k, n, l, o, m, Dim0, Dim1, Dim2, Dim5, Dim3, Dim6, Dim4, i, j, k, n, l, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim3, Dim2, m, n, o, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim3, Dim2, i, j, k, l, m, m, n, o, l, k, Dim0, Dim1, Dim4, Dim5, Dim6, Dim3, Dim2, i, j, m, n, o, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim4, Dim2, l, n, o, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim4, Dim2, i, j, k, l, m, l, n, o, m, k, Dim0, Dim1, Dim3, Dim5, Dim6, Dim4, Dim2, i, j, l, n, o, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim2, Dim3, m, n, o, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim2, Dim3, i, j, k, l, m, m, n, o, k, l, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, i, j, m, n, o, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim2, Dim4, l, n, o, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim2, Dim4, i, j, k, l, m, l, n, o, k, m, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, i, j, l, n, o, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim4, Dim3, k, n, o, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim4, Dim3, i, j, k, l, m, k, n, o, m, l, Dim0, Dim1, Dim2, Dim5, Dim6, Dim4, Dim3, i, j, k, n, o, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim3, Dim4, k, n, o, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim3, Dim4, i, j, k, l, m, k, n, o, l, m, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, i, j, k, n, o, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim3, Dim2, Dim6, n, m, l, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim2, Dim6, i, j, k, l, m, n, m, l, k, o, Dim0, Dim1, Dim5, Dim4, Dim3, Dim2, Dim6, i, j, n, m, l, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim4, Dim2, Dim6, n, l, m, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim2, Dim6, i, j, k, l, m, n, l, m, k, o, Dim0, Dim1, Dim5, Dim3, Dim4, Dim2, Dim6, i, j, n, l, m, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim2, Dim3, Dim6, n, m, k, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim3, Dim6, i, j, k, l, m, n, m, k, l, o, Dim0, Dim1, Dim5, Dim4, Dim2, Dim3, Dim6, i, j, n, m, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim2, Dim4, Dim6, n, l, k, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim4, Dim6, i, j, k, l, m, n, l, k, m, o, Dim0, Dim1, Dim5, Dim3, Dim2, Dim4, Dim6, i, j, n, l, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim4, Dim3, Dim6, n, k, m, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim3, Dim6, i, j, k, l, m, n, k, m, l, o, Dim0, Dim1, Dim5, Dim2, Dim4, Dim3, Dim6, i, j, n, k, m, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim3, Dim4, Dim6, n, k, l, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, k, l, m, n, k, l, m, o, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, Dim6, i, j, n, k, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim3, Dim6, Dim2, n, m, l, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim6, Dim2, i, j, k, l, m, n, m, l, o, k, Dim0, Dim1, Dim5, Dim4, Dim3, Dim6, Dim2, i, j, n, m, l, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim4, Dim6, Dim2, n, l, m, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim6, Dim2, i, j, k, l, m, n, l, m, o, k, Dim0, Dim1, Dim5, Dim3, Dim4, Dim6, Dim2, i, j, n, l, m, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim2, Dim6, Dim3, n, m, k, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim6, Dim3, i, j, k, l, m, n, m, k, o, l, Dim0, Dim1, Dim5, Dim4, Dim2, Dim6, Dim3, i, j, n, m, k, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim2, Dim6, Dim4, n, l, k, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim6, Dim4, i, j, k, l, m, n, l, k, o, m, Dim0, Dim1, Dim5, Dim3, Dim2, Dim6, Dim4, i, j, n, l, k, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim4, Dim6, Dim3, n, k, m, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim6, Dim3, i, j, k, l, m, n, k, m, o, l, Dim0, Dim1, Dim5, Dim2, Dim4, Dim6, Dim3, i, j, n, k, m, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim3, Dim6, Dim4, n, k, l, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim6, Dim4, i, j, k, l, m, n, k, l, o, m, Dim0, Dim1, Dim5, Dim2, Dim3, Dim6, Dim4, i, j, n, k, l, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim3, Dim2, n, m, o, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim3, Dim2, i, j, k, l, m, n, m, o, l, k, Dim0, Dim1, Dim5, Dim4, Dim6, Dim3, Dim2, i, j, n, m, o, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim4, Dim2, n, l, o, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim4, Dim2, i, j, k, l, m, n, l, o, m, k, Dim0, Dim1, Dim5, Dim3, Dim6, Dim4, Dim2, i, j, n, l, o, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim2, Dim3, n, m, o, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim2, Dim3, i, j, k, l, m, n, m, o, k, l, Dim0, Dim1, Dim5, Dim4, Dim6, Dim2, Dim3, i, j, n, m, o, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim2, Dim4, n, l, o, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim2, Dim4, i, j, k, l, m, n, l, o, k, m, Dim0, Dim1, Dim5, Dim3, Dim6, Dim2, Dim4, i, j, n, l, o, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim4, Dim3, n, k, o, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim4, Dim3, i, j, k, l, m, n, k, o, m, l, Dim0, Dim1, Dim5, Dim2, Dim6, Dim4, Dim3, i, j, n, k, o, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim3, Dim4, n, k, o, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim3, Dim4, i, j, k, l, m, n, k, o, l, m, Dim0, Dim1, Dim5, Dim2, Dim6, Dim3, Dim4, i, j, n, k, o, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim3, Dim2, n, o, m, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim2, i, j, k, l, m, n, o, m, l, k, Dim0, Dim1, Dim5, Dim6, Dim4, Dim3, Dim2, i, j, n, o, m, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim4, Dim2, n, o, l, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim2, i, j, k, l, m, n, o, l, m, k, Dim0, Dim1, Dim5, Dim6, Dim3, Dim4, Dim2, i, j, n, o, l, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim2, Dim3, n, o, m, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim3, i, j, k, l, m, n, o, m, k, l, Dim0, Dim1, Dim5, Dim6, Dim4, Dim2, Dim3, i, j, n, o, m, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim2, Dim4, n, o, l, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim4, i, j, k, l, m, n, o, l, k, m, Dim0, Dim1, Dim5, Dim6, Dim3, Dim2, Dim4, i, j, n, o, l, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim4, Dim3, n, o, k, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim3, i, j, k, l, m, n, o, k, m, l, Dim0, Dim1, Dim5, Dim6, Dim2, Dim4, Dim3, i, j, n, o, k, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim3, Dim4, n, o, k, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim4, i, j, k, l, m, n, o, k, l, m, Dim0, Dim1, Dim5, Dim6, Dim2, Dim3, Dim4, i, j, n, o, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

} //namespace FTensor

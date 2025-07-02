#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, n, o, )*B(p, q, r, ) triple contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, int Dim9, char i, char j, char k, char l, char m, char n, char o, char p, char q, char r, int DimA, int DimB, int DimC, int DimD, int DimX, int DimY, int DimZ, char a, char b, char c, char d, char x, char y, char z>
class Tensor7_times_Tensor3_triple
{
  Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
  Tensor3_Expr<B, U, Dim7, Dim8, Dim9, p, q, r> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        for(int zz = 0; zz < DimZ; ++zz)
          {
            // Permutation is where the indices get checked.
            result += Permutation7<DimA, DimB, DimC, DimD, DimX, DimY, DimZ, a, b, c, d, x, y, z>().eval(
              iterA, N1, N2, N3, N4, xx, yy, zz) * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx, yy, zz);
          }
    return result;
  }

  Tensor7_times_Tensor3_triple(
    const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &iter_a,
    const Tensor3_Expr<B, U, Dim7, Dim8, Dim9, p, q, r> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, n, o, )*B(k, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim0, i, j, k, l, m, n, o, k, j, i, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim0, l, m, n, o, k, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(k, j, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim0, i, j, k, l, m, n, o, j, k, i, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim0, l, m, n, o, j, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(j, k, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim1, i, j, k, l, m, n, o, k, i, j, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim1, l, m, n, o, k, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(k, i, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim2, i, j, k, l, m, n, o, j, i, k, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim2, l, m, n, o, j, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(j, i, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim1, i, j, k, l, m, n, o, i, k, j, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim1, l, m, n, o, i, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(i, k, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim2, i, j, k, l, m, n, o, i, j, k, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim2, l, m, n, o, i, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, Dim6, l, m, n, o>(TensorExpr(a, b));
}

//B(i, j, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim0, l, j, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim0, i, j, k, l, m, n, o, l, j, i, Dim2, Dim4, Dim5, Dim6, Dim3, Dim1, Dim0, k, m, n, o, l, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(l, j, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim0, l, j, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim0, j, l, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim0, i, j, k, l, m, n, o, j, l, i, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim0, k, m, n, o, j, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(j, l, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim0, j, l, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim1, l, i, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim1, i, j, k, l, m, n, o, l, i, j, Dim2, Dim4, Dim5, Dim6, Dim3, Dim0, Dim1, k, m, n, o, l, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(l, i, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim1, l, i, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim3, j, i, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim3, i, j, k, l, m, n, o, j, i, l, Dim2, Dim4, Dim5, Dim6, Dim1, Dim0, Dim3, k, m, n, o, j, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(j, i, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim3, j, i, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim1, i, l, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim1, i, j, k, l, m, n, o, i, l, j, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim1, k, m, n, o, i, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(i, l, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim1, i, l, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim3, i, j, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim3, i, j, k, l, m, n, o, i, j, l, Dim2, Dim4, Dim5, Dim6, Dim0, Dim1, Dim3, k, m, n, o, i, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, Dim6, k, m, n, o>(TensorExpr(a, b));
}

//B(i, j, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim3, i, j, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim0, m, j, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim0, i, j, k, l, m, n, o, m, j, i, Dim2, Dim3, Dim5, Dim6, Dim4, Dim1, Dim0, k, l, n, o, m, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(m, j, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim0, m, j, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim0, j, m, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim0, i, j, k, l, m, n, o, j, m, i, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim0, k, l, n, o, j, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(j, m, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim0, j, m, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim1, m, i, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim1, i, j, k, l, m, n, o, m, i, j, Dim2, Dim3, Dim5, Dim6, Dim4, Dim0, Dim1, k, l, n, o, m, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(m, i, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim1, m, i, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim4, j, i, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim4, i, j, k, l, m, n, o, j, i, m, Dim2, Dim3, Dim5, Dim6, Dim1, Dim0, Dim4, k, l, n, o, j, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(j, i, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim4, j, i, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim1, i, m, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim1, i, j, k, l, m, n, o, i, m, j, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim1, k, l, n, o, i, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(i, m, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim1, i, m, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim4, i, j, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim4, i, j, k, l, m, n, o, i, j, m, Dim2, Dim3, Dim5, Dim6, Dim0, Dim1, Dim4, k, l, n, o, i, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, Dim6, k, l, n, o>(TensorExpr(a, b));
}

//B(i, j, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim4, i, j, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim0, n, j, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim0, i, j, k, l, m, n, o, n, j, i, Dim2, Dim3, Dim4, Dim6, Dim5, Dim1, Dim0, k, l, m, o, n, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(n, j, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim0, n, j, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim0, j, n, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim0, i, j, k, l, m, n, o, j, n, i, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim0, k, l, m, o, j, n, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(j, n, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim0, j, n, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim1, n, i, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim1, i, j, k, l, m, n, o, n, i, j, Dim2, Dim3, Dim4, Dim6, Dim5, Dim0, Dim1, k, l, m, o, n, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(n, i, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim1, n, i, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim5, j, i, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim5, i, j, k, l, m, n, o, j, i, n, Dim2, Dim3, Dim4, Dim6, Dim1, Dim0, Dim5, k, l, m, o, j, i, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(j, i, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim5, j, i, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim1, i, n, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim1, i, j, k, l, m, n, o, i, n, j, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim1, k, l, m, o, i, n, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(i, n, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim1, i, n, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim5, i, j, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim5, i, j, k, l, m, n, o, i, j, n, Dim2, Dim3, Dim4, Dim6, Dim0, Dim1, Dim5, k, l, m, o, i, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim6, k, l, m, o>(TensorExpr(a, b));
}

//B(i, j, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim5, i, j, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim1, Dim0, o, j, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim1, Dim0, i, j, k, l, m, n, o, o, j, i, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, k, l, m, n, o, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//B(o, j, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim1, Dim0, o, j, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim6, Dim0, j, o, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim6, Dim0, i, j, k, l, m, n, o, j, o, i, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, k, l, m, n, j, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//B(j, o, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim6, Dim0, j, o, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim0, Dim1, o, i, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim0, Dim1, i, j, k, l, m, n, o, o, i, j, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, k, l, m, n, o, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//B(o, i, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim0, Dim1, o, i, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim6, j, i, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim6, i, j, k, l, m, n, o, j, i, o, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, k, l, m, n, j, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//B(j, i, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim6, j, i, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim6, Dim1, i, o, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim6, Dim1, i, j, k, l, m, n, o, i, o, j, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, k, l, m, n, i, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//B(i, o, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim6, Dim1, i, o, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim6, i, j, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim6, i, j, k, l, m, n, o, i, j, o, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, k, l, m, n, i, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//B(i, j, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim6, i, j, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim0, l, k, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim0, i, j, k, l, m, n, o, l, k, i, Dim1, Dim4, Dim5, Dim6, Dim3, Dim2, Dim0, j, m, n, o, l, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(l, k, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim0, l, k, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim0, k, l, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim0, i, j, k, l, m, n, o, k, l, i, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim0, j, m, n, o, k, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(k, l, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim0, k, l, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim2, l, i, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim2, i, j, k, l, m, n, o, l, i, k, Dim1, Dim4, Dim5, Dim6, Dim3, Dim0, Dim2, j, m, n, o, l, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(l, i, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim2, l, i, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim3, k, i, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim3, i, j, k, l, m, n, o, k, i, l, Dim1, Dim4, Dim5, Dim6, Dim2, Dim0, Dim3, j, m, n, o, k, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(k, i, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim3, k, i, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim2, i, l, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim2, i, j, k, l, m, n, o, i, l, k, Dim1, Dim4, Dim5, Dim6, Dim0, Dim3, Dim2, j, m, n, o, i, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(i, l, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim2, i, l, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim3, i, k, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim3, i, j, k, l, m, n, o, i, k, l, Dim1, Dim4, Dim5, Dim6, Dim0, Dim2, Dim3, j, m, n, o, i, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, Dim6, j, m, n, o>(TensorExpr(a, b));
}

//B(i, k, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim3, i, k, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim0, m, k, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim0, i, j, k, l, m, n, o, m, k, i, Dim1, Dim3, Dim5, Dim6, Dim4, Dim2, Dim0, j, l, n, o, m, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(m, k, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim0, m, k, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim0, k, m, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim0, i, j, k, l, m, n, o, k, m, i, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim0, j, l, n, o, k, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(k, m, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim0, k, m, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim2, m, i, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim2, i, j, k, l, m, n, o, m, i, k, Dim1, Dim3, Dim5, Dim6, Dim4, Dim0, Dim2, j, l, n, o, m, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(m, i, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim2, m, i, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim4, k, i, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim4, i, j, k, l, m, n, o, k, i, m, Dim1, Dim3, Dim5, Dim6, Dim2, Dim0, Dim4, j, l, n, o, k, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(k, i, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim4, k, i, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim2, i, m, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim2, i, j, k, l, m, n, o, i, m, k, Dim1, Dim3, Dim5, Dim6, Dim0, Dim4, Dim2, j, l, n, o, i, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(i, m, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim2, i, m, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim4, i, k, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim4, i, j, k, l, m, n, o, i, k, m, Dim1, Dim3, Dim5, Dim6, Dim0, Dim2, Dim4, j, l, n, o, i, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, Dim6, j, l, n, o>(TensorExpr(a, b));
}

//B(i, k, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim4, i, k, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim0, n, k, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim0, i, j, k, l, m, n, o, n, k, i, Dim1, Dim3, Dim4, Dim6, Dim5, Dim2, Dim0, j, l, m, o, n, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(n, k, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim0, n, k, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim0, k, n, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim0, i, j, k, l, m, n, o, k, n, i, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim0, j, l, m, o, k, n, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(k, n, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim0, k, n, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim2, n, i, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim2, i, j, k, l, m, n, o, n, i, k, Dim1, Dim3, Dim4, Dim6, Dim5, Dim0, Dim2, j, l, m, o, n, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(n, i, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim2, n, i, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim5, k, i, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim5, i, j, k, l, m, n, o, k, i, n, Dim1, Dim3, Dim4, Dim6, Dim2, Dim0, Dim5, j, l, m, o, k, i, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(k, i, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim5, k, i, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim2, i, n, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim2, i, j, k, l, m, n, o, i, n, k, Dim1, Dim3, Dim4, Dim6, Dim0, Dim5, Dim2, j, l, m, o, i, n, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(i, n, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim2, i, n, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim5, i, k, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim5, i, j, k, l, m, n, o, i, k, n, Dim1, Dim3, Dim4, Dim6, Dim0, Dim2, Dim5, j, l, m, o, i, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim6, j, l, m, o>(TensorExpr(a, b));
}

//B(i, k, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim5, i, k, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim2, Dim0, o, k, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim2, Dim0, i, j, k, l, m, n, o, o, k, i, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, j, l, m, n, o, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//B(o, k, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim2, Dim0, o, k, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim6, Dim0, k, o, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim6, Dim0, i, j, k, l, m, n, o, k, o, i, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, j, l, m, n, k, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//B(k, o, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim6, Dim0, k, o, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim0, Dim2, o, i, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim0, Dim2, i, j, k, l, m, n, o, o, i, k, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, j, l, m, n, o, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//B(o, i, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim0, Dim2, o, i, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim6, k, i, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim6, i, j, k, l, m, n, o, k, i, o, Dim1, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, j, l, m, n, k, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//B(k, i, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim6, k, i, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim6, Dim2, i, o, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim6, Dim2, i, j, k, l, m, n, o, i, o, k, Dim1, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, j, l, m, n, i, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//B(i, o, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim6, Dim2, i, o, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim6, i, k, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim6, i, j, k, l, m, n, o, i, k, o, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, j, l, m, n, i, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//B(i, k, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim6, i, k, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim0, m, l, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim0, i, j, k, l, m, n, o, m, l, i, Dim1, Dim2, Dim5, Dim6, Dim4, Dim3, Dim0, j, k, n, o, m, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(m, l, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim0, m, l, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim0, l, m, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim0, i, j, k, l, m, n, o, l, m, i, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim0, j, k, n, o, l, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(l, m, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim0, l, m, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim3, m, i, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim3, i, j, k, l, m, n, o, m, i, l, Dim1, Dim2, Dim5, Dim6, Dim4, Dim0, Dim3, j, k, n, o, m, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(m, i, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim3, m, i, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim4, l, i, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim4, i, j, k, l, m, n, o, l, i, m, Dim1, Dim2, Dim5, Dim6, Dim3, Dim0, Dim4, j, k, n, o, l, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(l, i, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim4, l, i, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim3, i, m, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim3, i, j, k, l, m, n, o, i, m, l, Dim1, Dim2, Dim5, Dim6, Dim0, Dim4, Dim3, j, k, n, o, i, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(i, m, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim3, i, m, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim4, i, l, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim4, i, j, k, l, m, n, o, i, l, m, Dim1, Dim2, Dim5, Dim6, Dim0, Dim3, Dim4, j, k, n, o, i, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, Dim6, j, k, n, o>(TensorExpr(a, b));
}

//B(i, l, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim4, i, l, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim0, n, l, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim0, i, j, k, l, m, n, o, n, l, i, Dim1, Dim2, Dim4, Dim6, Dim5, Dim3, Dim0, j, k, m, o, n, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(n, l, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim0, n, l, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim0, l, n, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim0, i, j, k, l, m, n, o, l, n, i, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim0, j, k, m, o, l, n, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(l, n, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim0, l, n, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim3, n, i, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim3, i, j, k, l, m, n, o, n, i, l, Dim1, Dim2, Dim4, Dim6, Dim5, Dim0, Dim3, j, k, m, o, n, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(n, i, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim3, n, i, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim5, l, i, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim5, i, j, k, l, m, n, o, l, i, n, Dim1, Dim2, Dim4, Dim6, Dim3, Dim0, Dim5, j, k, m, o, l, i, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(l, i, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim5, l, i, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim3, i, n, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim3, i, j, k, l, m, n, o, i, n, l, Dim1, Dim2, Dim4, Dim6, Dim0, Dim5, Dim3, j, k, m, o, i, n, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(i, n, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim3, i, n, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim5, i, l, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim5, i, j, k, l, m, n, o, i, l, n, Dim1, Dim2, Dim4, Dim6, Dim0, Dim3, Dim5, j, k, m, o, i, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim6, j, k, m, o>(TensorExpr(a, b));
}

//B(i, l, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim5, i, l, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim3, Dim0, o, l, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim3, Dim0, i, j, k, l, m, n, o, o, l, i, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, Dim0, j, k, m, n, o, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//B(o, l, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim3, Dim0, o, l, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim6, Dim0, l, o, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim6, Dim0, i, j, k, l, m, n, o, l, o, i, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim0, j, k, m, n, l, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//B(l, o, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim6, Dim0, l, o, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim0, Dim3, o, i, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim0, Dim3, i, j, k, l, m, n, o, o, i, l, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, j, k, m, n, o, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//B(o, i, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim0, Dim3, o, i, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim6, l, i, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim6, i, j, k, l, m, n, o, l, i, o, Dim1, Dim2, Dim4, Dim5, Dim3, Dim0, Dim6, j, k, m, n, l, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//B(l, i, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim6, l, i, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim6, Dim3, i, o, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim6, Dim3, i, j, k, l, m, n, o, i, o, l, Dim1, Dim2, Dim4, Dim5, Dim0, Dim6, Dim3, j, k, m, n, i, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//B(i, o, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim6, Dim3, i, o, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim6, i, l, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim6, i, j, k, l, m, n, o, i, l, o, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, j, k, m, n, i, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//B(i, l, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim6, i, l, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim0, n, m, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim0, i, j, k, l, m, n, o, n, m, i, Dim1, Dim2, Dim3, Dim6, Dim5, Dim4, Dim0, j, k, l, o, n, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(n, m, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim0, n, m, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim0, m, n, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim0, i, j, k, l, m, n, o, m, n, i, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim0, j, k, l, o, m, n, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(m, n, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim0, m, n, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim4, n, i, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim4, i, j, k, l, m, n, o, n, i, m, Dim1, Dim2, Dim3, Dim6, Dim5, Dim0, Dim4, j, k, l, o, n, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(n, i, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim4, n, i, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim5, m, i, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim5, i, j, k, l, m, n, o, m, i, n, Dim1, Dim2, Dim3, Dim6, Dim4, Dim0, Dim5, j, k, l, o, m, i, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(m, i, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim5, m, i, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim4, i, n, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim4, i, j, k, l, m, n, o, i, n, m, Dim1, Dim2, Dim3, Dim6, Dim0, Dim5, Dim4, j, k, l, o, i, n, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(i, n, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim4, i, n, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim5, i, m, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim5, i, j, k, l, m, n, o, i, m, n, Dim1, Dim2, Dim3, Dim6, Dim0, Dim4, Dim5, j, k, l, o, i, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim6, j, k, l, o>(TensorExpr(a, b));
}

//B(i, m, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim5, i, m, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim4, Dim0, o, m, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim4, Dim0, i, j, k, l, m, n, o, o, m, i, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, Dim0, j, k, l, n, o, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>(TensorExpr(a, b));
}

//B(o, m, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim4, Dim0, o, m, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim6, Dim0, m, o, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim6, Dim0, i, j, k, l, m, n, o, m, o, i, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim0, j, k, l, n, m, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>(TensorExpr(a, b));
}

//B(m, o, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim6, Dim0, m, o, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim0, Dim4, o, i, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim0, Dim4, i, j, k, l, m, n, o, o, i, m, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, j, k, l, n, o, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>(TensorExpr(a, b));
}

//B(o, i, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim0, Dim4, o, i, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim6, m, i, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim6, i, j, k, l, m, n, o, m, i, o, Dim1, Dim2, Dim3, Dim5, Dim4, Dim0, Dim6, j, k, l, n, m, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>(TensorExpr(a, b));
}

//B(m, i, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim6, m, i, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim6, Dim4, i, o, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim6, Dim4, i, j, k, l, m, n, o, i, o, m, Dim1, Dim2, Dim3, Dim5, Dim0, Dim6, Dim4, j, k, l, n, i, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>(TensorExpr(a, b));
}

//B(i, o, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim6, Dim4, i, o, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim6, i, m, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim6, i, j, k, l, m, n, o, i, m, o, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, j, k, l, n, i, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>(TensorExpr(a, b));
}

//B(i, m, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim6, i, m, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim5, Dim0, o, n, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim5, Dim0, i, j, k, l, m, n, o, o, n, i, Dim1, Dim2, Dim3, Dim4, Dim6, Dim5, Dim0, j, k, l, m, o, n, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>(TensorExpr(a, b));
}

//B(o, n, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim5, Dim0, o, n, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim0, n, o, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim6, Dim0, i, j, k, l, m, n, o, n, o, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, j, k, l, m, n, o, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>(TensorExpr(a, b));
}

//B(n, o, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim0, n, o, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim0, Dim5, o, i, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim0, Dim5, i, j, k, l, m, n, o, o, i, n, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, j, k, l, m, o, i, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>(TensorExpr(a, b));
}

//B(o, i, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim0, Dim5, o, i, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim6, n, i, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim6, i, j, k, l, m, n, o, n, i, o, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, j, k, l, m, n, i, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>(TensorExpr(a, b));
}

//B(n, i, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim6, n, i, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim6, Dim5, i, o, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim6, Dim5, i, j, k, l, m, n, o, i, o, n, Dim1, Dim2, Dim3, Dim4, Dim0, Dim6, Dim5, j, k, l, m, i, o, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>(TensorExpr(a, b));
}

//B(i, o, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim6, Dim5, i, o, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim6, i, n, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim6, i, j, k, l, m, n, o, i, n, o, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, j, k, l, m, i, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>(TensorExpr(a, b));
}

//B(i, n, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim6, i, n, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim1, l, k, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim1, i, j, k, l, m, n, o, l, k, j, Dim0, Dim4, Dim5, Dim6, Dim3, Dim2, Dim1, i, m, n, o, l, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(l, k, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim1, l, k, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim1, k, l, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim1, i, j, k, l, m, n, o, k, l, j, Dim0, Dim4, Dim5, Dim6, Dim2, Dim3, Dim1, i, m, n, o, k, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(k, l, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim1, k, l, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim2, l, j, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim2, i, j, k, l, m, n, o, l, j, k, Dim0, Dim4, Dim5, Dim6, Dim3, Dim1, Dim2, i, m, n, o, l, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(l, j, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim2, l, j, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim3, k, j, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim3, i, j, k, l, m, n, o, k, j, l, Dim0, Dim4, Dim5, Dim6, Dim2, Dim1, Dim3, i, m, n, o, k, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(k, j, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim3, k, j, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim2, j, l, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim2, i, j, k, l, m, n, o, j, l, k, Dim0, Dim4, Dim5, Dim6, Dim1, Dim3, Dim2, i, m, n, o, j, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(j, l, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim2, j, l, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim3, j, k, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim3, i, j, k, l, m, n, o, j, k, l, Dim0, Dim4, Dim5, Dim6, Dim1, Dim2, Dim3, i, m, n, o, j, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, Dim6, i, m, n, o>(TensorExpr(a, b));
}

//B(j, k, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim3, j, k, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim1, m, k, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim1, i, j, k, l, m, n, o, m, k, j, Dim0, Dim3, Dim5, Dim6, Dim4, Dim2, Dim1, i, l, n, o, m, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(m, k, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim1, m, k, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim1, k, m, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim1, i, j, k, l, m, n, o, k, m, j, Dim0, Dim3, Dim5, Dim6, Dim2, Dim4, Dim1, i, l, n, o, k, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(k, m, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim1, k, m, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim2, m, j, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim2, i, j, k, l, m, n, o, m, j, k, Dim0, Dim3, Dim5, Dim6, Dim4, Dim1, Dim2, i, l, n, o, m, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(m, j, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim2, m, j, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim4, k, j, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim4, i, j, k, l, m, n, o, k, j, m, Dim0, Dim3, Dim5, Dim6, Dim2, Dim1, Dim4, i, l, n, o, k, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(k, j, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim4, k, j, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim2, j, m, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim2, i, j, k, l, m, n, o, j, m, k, Dim0, Dim3, Dim5, Dim6, Dim1, Dim4, Dim2, i, l, n, o, j, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(j, m, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim2, j, m, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim4, j, k, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim4, i, j, k, l, m, n, o, j, k, m, Dim0, Dim3, Dim5, Dim6, Dim1, Dim2, Dim4, i, l, n, o, j, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, Dim6, i, l, n, o>(TensorExpr(a, b));
}

//B(j, k, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim4, j, k, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim1, n, k, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim1, i, j, k, l, m, n, o, n, k, j, Dim0, Dim3, Dim4, Dim6, Dim5, Dim2, Dim1, i, l, m, o, n, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(n, k, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim1, n, k, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim1, k, n, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim1, i, j, k, l, m, n, o, k, n, j, Dim0, Dim3, Dim4, Dim6, Dim2, Dim5, Dim1, i, l, m, o, k, n, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(k, n, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim1, k, n, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim2, n, j, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim2, i, j, k, l, m, n, o, n, j, k, Dim0, Dim3, Dim4, Dim6, Dim5, Dim1, Dim2, i, l, m, o, n, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(n, j, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim2, n, j, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim5, k, j, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim5, i, j, k, l, m, n, o, k, j, n, Dim0, Dim3, Dim4, Dim6, Dim2, Dim1, Dim5, i, l, m, o, k, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(k, j, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim5, k, j, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim2, j, n, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim2, i, j, k, l, m, n, o, j, n, k, Dim0, Dim3, Dim4, Dim6, Dim1, Dim5, Dim2, i, l, m, o, j, n, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(j, n, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim2, j, n, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim5, j, k, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim5, i, j, k, l, m, n, o, j, k, n, Dim0, Dim3, Dim4, Dim6, Dim1, Dim2, Dim5, i, l, m, o, j, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim6, i, l, m, o>(TensorExpr(a, b));
}

//B(j, k, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim5, j, k, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim2, Dim1, o, k, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim2, Dim1, i, j, k, l, m, n, o, o, k, j, Dim0, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, i, l, m, n, o, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//B(o, k, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim2, Dim1, o, k, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim6, Dim1, k, o, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim6, Dim1, i, j, k, l, m, n, o, k, o, j, Dim0, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, i, l, m, n, k, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//B(k, o, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim6, Dim1, k, o, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim1, Dim2, o, j, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim1, Dim2, i, j, k, l, m, n, o, o, j, k, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, i, l, m, n, o, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//B(o, j, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim1, Dim2, o, j, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim6, k, j, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim6, i, j, k, l, m, n, o, k, j, o, Dim0, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, i, l, m, n, k, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//B(k, j, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim6, k, j, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim6, Dim2, j, o, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim6, Dim2, i, j, k, l, m, n, o, j, o, k, Dim0, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, i, l, m, n, j, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//B(j, o, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim6, Dim2, j, o, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim6, j, k, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim6, i, j, k, l, m, n, o, j, k, o, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, i, l, m, n, j, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//B(j, k, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim6, j, k, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim1, m, l, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim1, i, j, k, l, m, n, o, m, l, j, Dim0, Dim2, Dim5, Dim6, Dim4, Dim3, Dim1, i, k, n, o, m, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(m, l, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim1, m, l, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim1, l, m, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim1, i, j, k, l, m, n, o, l, m, j, Dim0, Dim2, Dim5, Dim6, Dim3, Dim4, Dim1, i, k, n, o, l, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(l, m, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim1, l, m, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim3, m, j, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim3, i, j, k, l, m, n, o, m, j, l, Dim0, Dim2, Dim5, Dim6, Dim4, Dim1, Dim3, i, k, n, o, m, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(m, j, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim3, m, j, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim4, l, j, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim4, i, j, k, l, m, n, o, l, j, m, Dim0, Dim2, Dim5, Dim6, Dim3, Dim1, Dim4, i, k, n, o, l, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(l, j, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim4, l, j, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim3, j, m, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim3, i, j, k, l, m, n, o, j, m, l, Dim0, Dim2, Dim5, Dim6, Dim1, Dim4, Dim3, i, k, n, o, j, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(j, m, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim3, j, m, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim4, j, l, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim4, i, j, k, l, m, n, o, j, l, m, Dim0, Dim2, Dim5, Dim6, Dim1, Dim3, Dim4, i, k, n, o, j, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, Dim6, i, k, n, o>(TensorExpr(a, b));
}

//B(j, l, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim4, j, l, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim1, n, l, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim1, i, j, k, l, m, n, o, n, l, j, Dim0, Dim2, Dim4, Dim6, Dim5, Dim3, Dim1, i, k, m, o, n, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(n, l, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim1, n, l, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim1, l, n, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim1, i, j, k, l, m, n, o, l, n, j, Dim0, Dim2, Dim4, Dim6, Dim3, Dim5, Dim1, i, k, m, o, l, n, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(l, n, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim1, l, n, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim3, n, j, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim3, i, j, k, l, m, n, o, n, j, l, Dim0, Dim2, Dim4, Dim6, Dim5, Dim1, Dim3, i, k, m, o, n, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(n, j, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim3, n, j, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim5, l, j, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim5, i, j, k, l, m, n, o, l, j, n, Dim0, Dim2, Dim4, Dim6, Dim3, Dim1, Dim5, i, k, m, o, l, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(l, j, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim5, l, j, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim3, j, n, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim3, i, j, k, l, m, n, o, j, n, l, Dim0, Dim2, Dim4, Dim6, Dim1, Dim5, Dim3, i, k, m, o, j, n, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(j, n, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim3, j, n, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim5, j, l, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim5, i, j, k, l, m, n, o, j, l, n, Dim0, Dim2, Dim4, Dim6, Dim1, Dim3, Dim5, i, k, m, o, j, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim6, i, k, m, o>(TensorExpr(a, b));
}

//B(j, l, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim5, j, l, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim3, Dim1, o, l, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim3, Dim1, i, j, k, l, m, n, o, o, l, j, Dim0, Dim2, Dim4, Dim5, Dim6, Dim3, Dim1, i, k, m, n, o, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//B(o, l, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim3, Dim1, o, l, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim6, Dim1, l, o, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim6, Dim1, i, j, k, l, m, n, o, l, o, j, Dim0, Dim2, Dim4, Dim5, Dim3, Dim6, Dim1, i, k, m, n, l, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//B(l, o, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim6, Dim1, l, o, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim1, Dim3, o, j, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim1, Dim3, i, j, k, l, m, n, o, o, j, l, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, i, k, m, n, o, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//B(o, j, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim1, Dim3, o, j, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim6, l, j, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim6, i, j, k, l, m, n, o, l, j, o, Dim0, Dim2, Dim4, Dim5, Dim3, Dim1, Dim6, i, k, m, n, l, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//B(l, j, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim6, l, j, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim6, Dim3, j, o, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim6, Dim3, i, j, k, l, m, n, o, j, o, l, Dim0, Dim2, Dim4, Dim5, Dim1, Dim6, Dim3, i, k, m, n, j, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//B(j, o, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim6, Dim3, j, o, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim6, j, l, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim6, i, j, k, l, m, n, o, j, l, o, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, i, k, m, n, j, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//B(j, l, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim6, j, l, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim1, n, m, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim1, i, j, k, l, m, n, o, n, m, j, Dim0, Dim2, Dim3, Dim6, Dim5, Dim4, Dim1, i, k, l, o, n, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(n, m, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim1, n, m, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim1, m, n, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim1, i, j, k, l, m, n, o, m, n, j, Dim0, Dim2, Dim3, Dim6, Dim4, Dim5, Dim1, i, k, l, o, m, n, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(m, n, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim1, m, n, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim4, n, j, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim4, i, j, k, l, m, n, o, n, j, m, Dim0, Dim2, Dim3, Dim6, Dim5, Dim1, Dim4, i, k, l, o, n, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(n, j, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim4, n, j, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim5, m, j, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim5, i, j, k, l, m, n, o, m, j, n, Dim0, Dim2, Dim3, Dim6, Dim4, Dim1, Dim5, i, k, l, o, m, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(m, j, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim5, m, j, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim4, j, n, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim4, i, j, k, l, m, n, o, j, n, m, Dim0, Dim2, Dim3, Dim6, Dim1, Dim5, Dim4, i, k, l, o, j, n, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(j, n, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim4, j, n, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim5, j, m, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim5, i, j, k, l, m, n, o, j, m, n, Dim0, Dim2, Dim3, Dim6, Dim1, Dim4, Dim5, i, k, l, o, j, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim6, i, k, l, o>(TensorExpr(a, b));
}

//B(j, m, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim5, j, m, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim4, Dim1, o, m, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim4, Dim1, i, j, k, l, m, n, o, o, m, j, Dim0, Dim2, Dim3, Dim5, Dim6, Dim4, Dim1, i, k, l, n, o, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>(TensorExpr(a, b));
}

//B(o, m, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim4, Dim1, o, m, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim6, Dim1, m, o, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim6, Dim1, i, j, k, l, m, n, o, m, o, j, Dim0, Dim2, Dim3, Dim5, Dim4, Dim6, Dim1, i, k, l, n, m, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>(TensorExpr(a, b));
}

//B(m, o, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim6, Dim1, m, o, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim1, Dim4, o, j, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim1, Dim4, i, j, k, l, m, n, o, o, j, m, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, i, k, l, n, o, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>(TensorExpr(a, b));
}

//B(o, j, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim1, Dim4, o, j, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim6, m, j, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim6, i, j, k, l, m, n, o, m, j, o, Dim0, Dim2, Dim3, Dim5, Dim4, Dim1, Dim6, i, k, l, n, m, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>(TensorExpr(a, b));
}

//B(m, j, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim6, m, j, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim6, Dim4, j, o, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim6, Dim4, i, j, k, l, m, n, o, j, o, m, Dim0, Dim2, Dim3, Dim5, Dim1, Dim6, Dim4, i, k, l, n, j, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>(TensorExpr(a, b));
}

//B(j, o, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim6, Dim4, j, o, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim6, j, m, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim6, i, j, k, l, m, n, o, j, m, o, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, i, k, l, n, j, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>(TensorExpr(a, b));
}

//B(j, m, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim6, j, m, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim5, Dim1, o, n, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim5, Dim1, i, j, k, l, m, n, o, o, n, j, Dim0, Dim2, Dim3, Dim4, Dim6, Dim5, Dim1, i, k, l, m, o, n, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>(TensorExpr(a, b));
}

//B(o, n, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim5, Dim1, o, n, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim1, n, o, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim6, Dim1, i, j, k, l, m, n, o, n, o, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, i, k, l, m, n, o, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>(TensorExpr(a, b));
}

//B(n, o, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim1, n, o, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim1, Dim5, o, j, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim1, Dim5, i, j, k, l, m, n, o, o, j, n, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, i, k, l, m, o, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>(TensorExpr(a, b));
}

//B(o, j, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim1, Dim5, o, j, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim6, n, j, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim6, i, j, k, l, m, n, o, n, j, o, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, i, k, l, m, n, j, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>(TensorExpr(a, b));
}

//B(n, j, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim6, n, j, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim6, Dim5, j, o, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim6, Dim5, i, j, k, l, m, n, o, j, o, n, Dim0, Dim2, Dim3, Dim4, Dim1, Dim6, Dim5, i, k, l, m, j, o, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>(TensorExpr(a, b));
}

//B(j, o, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim6, Dim5, j, o, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim6, j, n, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim6, i, j, k, l, m, n, o, j, n, o, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, i, k, l, m, j, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>(TensorExpr(a, b));
}

//B(j, n, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim6, j, n, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim2, m, l, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim2, i, j, k, l, m, n, o, m, l, k, Dim0, Dim1, Dim5, Dim6, Dim4, Dim3, Dim2, i, j, n, o, m, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(m, l, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim2, m, l, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim2, l, m, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim2, i, j, k, l, m, n, o, l, m, k, Dim0, Dim1, Dim5, Dim6, Dim3, Dim4, Dim2, i, j, n, o, l, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(l, m, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim2, l, m, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim3, m, k, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim3, i, j, k, l, m, n, o, m, k, l, Dim0, Dim1, Dim5, Dim6, Dim4, Dim2, Dim3, i, j, n, o, m, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(m, k, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim3, m, k, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim4, l, k, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim4, i, j, k, l, m, n, o, l, k, m, Dim0, Dim1, Dim5, Dim6, Dim3, Dim2, Dim4, i, j, n, o, l, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(l, k, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim4, l, k, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim3, k, m, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim3, i, j, k, l, m, n, o, k, m, l, Dim0, Dim1, Dim5, Dim6, Dim2, Dim4, Dim3, i, j, n, o, k, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(k, m, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim3, k, m, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim4, i, j, k, l, m, n, o, k, l, m, Dim0, Dim1, Dim5, Dim6, Dim2, Dim3, Dim4, i, j, n, o, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, Dim6, i, j, n, o>(TensorExpr(a, b));
}

//B(k, l, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim2, n, l, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim2, i, j, k, l, m, n, o, n, l, k, Dim0, Dim1, Dim4, Dim6, Dim5, Dim3, Dim2, i, j, m, o, n, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(n, l, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim2, n, l, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim2, l, n, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim2, i, j, k, l, m, n, o, l, n, k, Dim0, Dim1, Dim4, Dim6, Dim3, Dim5, Dim2, i, j, m, o, l, n, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(l, n, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim2, l, n, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim3, n, k, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim3, i, j, k, l, m, n, o, n, k, l, Dim0, Dim1, Dim4, Dim6, Dim5, Dim2, Dim3, i, j, m, o, n, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(n, k, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim3, n, k, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim5, l, k, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim5, i, j, k, l, m, n, o, l, k, n, Dim0, Dim1, Dim4, Dim6, Dim3, Dim2, Dim5, i, j, m, o, l, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(l, k, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim5, l, k, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim3, k, n, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim3, i, j, k, l, m, n, o, k, n, l, Dim0, Dim1, Dim4, Dim6, Dim2, Dim5, Dim3, i, j, m, o, k, n, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(k, n, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim3, k, n, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim5, k, l, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim5, i, j, k, l, m, n, o, k, l, n, Dim0, Dim1, Dim4, Dim6, Dim2, Dim3, Dim5, i, j, m, o, k, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim6, i, j, m, o>(TensorExpr(a, b));
}

//B(k, l, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim5, k, l, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim3, Dim2, o, l, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim3, Dim2, i, j, k, l, m, n, o, o, l, k, Dim0, Dim1, Dim4, Dim5, Dim6, Dim3, Dim2, i, j, m, n, o, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//B(o, l, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim3, Dim2, o, l, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim6, Dim2, l, o, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim6, Dim2, i, j, k, l, m, n, o, l, o, k, Dim0, Dim1, Dim4, Dim5, Dim3, Dim6, Dim2, i, j, m, n, l, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//B(l, o, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim6, Dim2, l, o, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim2, Dim3, o, k, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim2, Dim3, i, j, k, l, m, n, o, o, k, l, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, i, j, m, n, o, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//B(o, k, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim2, Dim3, o, k, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim6, l, k, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim6, i, j, k, l, m, n, o, l, k, o, Dim0, Dim1, Dim4, Dim5, Dim3, Dim2, Dim6, i, j, m, n, l, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//B(l, k, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim6, l, k, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim6, Dim3, k, o, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim6, Dim3, i, j, k, l, m, n, o, k, o, l, Dim0, Dim1, Dim4, Dim5, Dim2, Dim6, Dim3, i, j, m, n, k, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//B(k, o, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim6, Dim3, k, o, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim6, k, l, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim6, i, j, k, l, m, n, o, k, l, o, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, i, j, m, n, k, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//B(k, l, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim6, k, l, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim2, n, m, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim2, i, j, k, l, m, n, o, n, m, k, Dim0, Dim1, Dim3, Dim6, Dim5, Dim4, Dim2, i, j, l, o, n, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(n, m, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim2, n, m, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim2, m, n, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim2, i, j, k, l, m, n, o, m, n, k, Dim0, Dim1, Dim3, Dim6, Dim4, Dim5, Dim2, i, j, l, o, m, n, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(m, n, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim2, m, n, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim4, n, k, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim4, i, j, k, l, m, n, o, n, k, m, Dim0, Dim1, Dim3, Dim6, Dim5, Dim2, Dim4, i, j, l, o, n, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(n, k, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim4, n, k, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim5, m, k, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim5, i, j, k, l, m, n, o, m, k, n, Dim0, Dim1, Dim3, Dim6, Dim4, Dim2, Dim5, i, j, l, o, m, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(m, k, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim5, m, k, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim4, k, n, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim4, i, j, k, l, m, n, o, k, n, m, Dim0, Dim1, Dim3, Dim6, Dim2, Dim5, Dim4, i, j, l, o, k, n, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(k, n, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim4, k, n, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim5, k, m, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim5, i, j, k, l, m, n, o, k, m, n, Dim0, Dim1, Dim3, Dim6, Dim2, Dim4, Dim5, i, j, l, o, k, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim6, i, j, l, o>(TensorExpr(a, b));
}

//B(k, m, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim5, k, m, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim4, Dim2, o, m, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim4, Dim2, i, j, k, l, m, n, o, o, m, k, Dim0, Dim1, Dim3, Dim5, Dim6, Dim4, Dim2, i, j, l, n, o, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>(TensorExpr(a, b));
}

//B(o, m, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim4, Dim2, o, m, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim6, Dim2, m, o, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim6, Dim2, i, j, k, l, m, n, o, m, o, k, Dim0, Dim1, Dim3, Dim5, Dim4, Dim6, Dim2, i, j, l, n, m, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>(TensorExpr(a, b));
}

//B(m, o, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim6, Dim2, m, o, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim2, Dim4, o, k, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim2, Dim4, i, j, k, l, m, n, o, o, k, m, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, i, j, l, n, o, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>(TensorExpr(a, b));
}

//B(o, k, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim2, Dim4, o, k, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim6, m, k, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim6, i, j, k, l, m, n, o, m, k, o, Dim0, Dim1, Dim3, Dim5, Dim4, Dim2, Dim6, i, j, l, n, m, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>(TensorExpr(a, b));
}

//B(m, k, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim6, m, k, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim6, Dim4, k, o, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim6, Dim4, i, j, k, l, m, n, o, k, o, m, Dim0, Dim1, Dim3, Dim5, Dim2, Dim6, Dim4, i, j, l, n, k, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>(TensorExpr(a, b));
}

//B(k, o, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim6, Dim4, k, o, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim6, k, m, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim6, i, j, k, l, m, n, o, k, m, o, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, i, j, l, n, k, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>(TensorExpr(a, b));
}

//B(k, m, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim6, k, m, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim5, Dim2, o, n, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim5, Dim2, i, j, k, l, m, n, o, o, n, k, Dim0, Dim1, Dim3, Dim4, Dim6, Dim5, Dim2, i, j, l, m, o, n, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>(TensorExpr(a, b));
}

//B(o, n, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim5, Dim2, o, n, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim2, n, o, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim6, Dim2, i, j, k, l, m, n, o, n, o, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, i, j, l, m, n, o, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>(TensorExpr(a, b));
}

//B(n, o, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim2, n, o, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim2, Dim5, o, k, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim2, Dim5, i, j, k, l, m, n, o, o, k, n, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, i, j, l, m, o, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>(TensorExpr(a, b));
}

//B(o, k, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim2, Dim5, o, k, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim6, n, k, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim6, i, j, k, l, m, n, o, n, k, o, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, i, j, l, m, n, k, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>(TensorExpr(a, b));
}

//B(n, k, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim6, n, k, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim6, Dim5, k, o, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim6, Dim5, i, j, k, l, m, n, o, k, o, n, Dim0, Dim1, Dim3, Dim4, Dim2, Dim6, Dim5, i, j, l, m, k, o, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>(TensorExpr(a, b));
}

//B(k, o, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim6, Dim5, k, o, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim6, k, n, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim6, i, j, k, l, m, n, o, k, n, o, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, i, j, l, m, k, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>(TensorExpr(a, b));
}

//B(k, n, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim6, k, n, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim3, n, m, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim3, i, j, k, l, m, n, o, n, m, l, Dim0, Dim1, Dim2, Dim6, Dim5, Dim4, Dim3, i, j, k, o, n, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(n, m, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim3, n, m, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim3, m, n, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim3, i, j, k, l, m, n, o, m, n, l, Dim0, Dim1, Dim2, Dim6, Dim4, Dim5, Dim3, i, j, k, o, m, n, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(m, n, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim3, m, n, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim4, n, l, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim4, i, j, k, l, m, n, o, n, l, m, Dim0, Dim1, Dim2, Dim6, Dim5, Dim3, Dim4, i, j, k, o, n, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(n, l, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim4, n, l, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim5, m, l, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim5, i, j, k, l, m, n, o, m, l, n, Dim0, Dim1, Dim2, Dim6, Dim4, Dim3, Dim5, i, j, k, o, m, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(m, l, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim5, m, l, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim4, l, n, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim4, i, j, k, l, m, n, o, l, n, m, Dim0, Dim1, Dim2, Dim6, Dim3, Dim5, Dim4, i, j, k, o, l, n, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(l, n, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim4, l, n, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim5, i, j, k, l, m, n, o, l, m, n, Dim0, Dim1, Dim2, Dim6, Dim3, Dim4, Dim5, i, j, k, o, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim6, i, j, k, o>(TensorExpr(a, b));
}

//B(l, m, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim4, Dim3, o, m, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim4, Dim3, i, j, k, l, m, n, o, o, m, l, Dim0, Dim1, Dim2, Dim5, Dim6, Dim4, Dim3, i, j, k, n, o, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>(TensorExpr(a, b));
}

//B(o, m, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim4, Dim3, o, m, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim6, Dim3, m, o, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim6, Dim3, i, j, k, l, m, n, o, m, o, l, Dim0, Dim1, Dim2, Dim5, Dim4, Dim6, Dim3, i, j, k, n, m, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>(TensorExpr(a, b));
}

//B(m, o, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim6, Dim3, m, o, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim3, Dim4, o, l, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim3, Dim4, i, j, k, l, m, n, o, o, l, m, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, i, j, k, n, o, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>(TensorExpr(a, b));
}

//B(o, l, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim3, Dim4, o, l, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim6, m, l, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim6, i, j, k, l, m, n, o, m, l, o, Dim0, Dim1, Dim2, Dim5, Dim4, Dim3, Dim6, i, j, k, n, m, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>(TensorExpr(a, b));
}

//B(m, l, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim6, m, l, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim6, Dim4, l, o, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim6, Dim4, i, j, k, l, m, n, o, l, o, m, Dim0, Dim1, Dim2, Dim5, Dim3, Dim6, Dim4, i, j, k, n, l, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>(TensorExpr(a, b));
}

//B(l, o, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim6, Dim4, l, o, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim6, l, m, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim6, i, j, k, l, m, n, o, l, m, o, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, i, j, k, n, l, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>(TensorExpr(a, b));
}

//B(l, m, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim6, l, m, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim5, Dim3, o, n, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim5, Dim3, i, j, k, l, m, n, o, o, n, l, Dim0, Dim1, Dim2, Dim4, Dim6, Dim5, Dim3, i, j, k, m, o, n, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>(TensorExpr(a, b));
}

//B(o, n, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim5, Dim3, o, n, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim3, n, o, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim6, Dim3, i, j, k, l, m, n, o, n, o, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, i, j, k, m, n, o, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>(TensorExpr(a, b));
}

//B(n, o, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim3, n, o, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim3, Dim5, o, l, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim3, Dim5, i, j, k, l, m, n, o, o, l, n, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, i, j, k, m, o, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>(TensorExpr(a, b));
}

//B(o, l, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim3, Dim5, o, l, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim6, n, l, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim6, i, j, k, l, m, n, o, n, l, o, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, i, j, k, m, n, l, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>(TensorExpr(a, b));
}

//B(n, l, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim6, n, l, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim6, Dim5, l, o, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim6, Dim5, i, j, k, l, m, n, o, l, o, n, Dim0, Dim1, Dim2, Dim4, Dim3, Dim6, Dim5, i, j, k, m, l, o, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>(TensorExpr(a, b));
}

//B(l, o, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim6, Dim5, l, o, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim6, l, n, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim6, i, j, k, l, m, n, o, l, n, o, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, i, j, k, m, l, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>(TensorExpr(a, b));
}

//B(l, n, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim6, l, n, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim5, Dim4, o, n, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim5, Dim4, i, j, k, l, m, n, o, o, n, m, Dim0, Dim1, Dim2, Dim3, Dim6, Dim5, Dim4, i, j, k, l, o, n, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
}

//B(o, n, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim5, Dim4, o, n, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim4, n, o, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim6, Dim4, i, j, k, l, m, n, o, n, o, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, i, j, k, l, n, o, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
}

//B(n, o, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim4, n, o, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim4, Dim5, o, m, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim4, Dim5, i, j, k, l, m, n, o, o, m, n, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, i, j, k, l, o, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
}

//B(o, m, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim4, Dim5, o, m, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim6, n, m, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim6, i, j, k, l, m, n, o, n, m, o, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, i, j, k, l, n, m, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
}

//B(n, m, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim6, n, m, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim6, Dim5, m, o, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim6, Dim5, i, j, k, l, m, n, o, m, o, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim5, i, j, k, l, m, o, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
}

//B(m, o, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim6, Dim5, m, o, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim6, i, j, k, l, m, n, o, m, n, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
}

//B(m, n, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

} //namespace FTensor

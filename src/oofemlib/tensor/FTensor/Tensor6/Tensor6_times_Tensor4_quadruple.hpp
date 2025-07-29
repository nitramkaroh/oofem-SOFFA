#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, n, )*B(o, p, q, r, ) quadruple contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, int Dim9, char i, char j, char k, char l, char m, char n, char o, char p, char q, char r, int DimA, int DimB, int DimX, int DimY, int DimZ, int DimW, char a, char b, char x, char y, char z, char w>
class Tensor6_times_Tensor4_quadruple
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim9, o, p, q, r> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        for(int zz = 0; zz < DimZ; ++zz)
          for(int ww = 0; ww < DimW; ++ww)
            {
              // Permutation is where the indices get checked.
              result += Permutation6<DimA, DimB, DimX, DimY, DimZ, DimW, a, b, x, y, z, w>().eval(
                iterA, N1, N2, xx, yy, zz, ww) * Permutation4<DimX, DimY, DimZ, DimW, x, y, z, w>().eval(iterB, xx, yy, zz, ww);
            }
    return result;
  }

  Tensor6_times_Tensor4_quadruple(
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &iter_a,
    const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim9, o, p, q, r> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, n, )*B(l, k, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim0, l, k, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim1, Dim0, i, j, k, l, m, n, l, k, j, i, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(l, k, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim0, l, k, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim0, k, l, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim1, Dim0, i, j, k, l, m, n, k, l, j, i, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(k, l, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim0, k, l, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim0, l, j, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim2, Dim0, i, j, k, l, m, n, l, j, k, i, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(l, j, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim0, l, j, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim0, k, j, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim3, Dim0, i, j, k, l, m, n, k, j, l, i, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(k, j, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim0, k, j, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim0, j, l, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim2, Dim0, i, j, k, l, m, n, j, l, k, i, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(j, l, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim0, j, l, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim0, j, k, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim3, Dim0, i, j, k, l, m, n, j, k, l, i, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(j, k, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim0, j, k, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim1, l, k, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim0, Dim1, i, j, k, l, m, n, l, k, i, j, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(l, k, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim1, l, k, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim1, k, l, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim0, Dim1, i, j, k, l, m, n, k, l, i, j, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(k, l, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim1, k, l, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim2, l, j, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim0, Dim2, i, j, k, l, m, n, l, j, i, k, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(l, j, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim2, l, j, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim3, k, j, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim0, Dim3, i, j, k, l, m, n, k, j, i, l, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(k, j, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim3, k, j, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim2, j, l, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim0, Dim2, i, j, k, l, m, n, j, l, i, k, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(j, l, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim2, j, l, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim3, j, k, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim0, Dim3, i, j, k, l, m, n, j, k, i, l, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(j, k, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim3, j, k, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim1, l, i, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim2, Dim1, i, j, k, l, m, n, l, i, k, j, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(l, i, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim1, l, i, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim1, k, i, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim3, Dim1, i, j, k, l, m, n, k, i, l, j, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(k, i, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim1, k, i, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim2, l, i, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim1, Dim2, i, j, k, l, m, n, l, i, j, k, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(l, i, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim2, l, i, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim3, k, i, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim1, Dim3, i, j, k, l, m, n, k, i, j, l, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(k, i, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim3, k, i, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim2, j, i, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim3, Dim2, i, j, k, l, m, n, j, i, l, k, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(j, i, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim2, j, i, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim3, j, i, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim2, Dim3, i, j, k, l, m, n, j, i, k, l, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(j, i, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim3, j, i, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim1, i, l, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim2, Dim1, i, j, k, l, m, n, i, l, k, j, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(i, l, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim1, i, l, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim1, i, k, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim3, Dim1, i, j, k, l, m, n, i, k, l, j, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(i, k, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim1, i, k, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim2, i, l, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim1, Dim2, i, j, k, l, m, n, i, l, j, k, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(i, l, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim2, i, l, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim3, i, k, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim1, Dim3, i, j, k, l, m, n, i, k, j, l, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(i, k, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim3, i, k, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim2, i, j, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim3, Dim2, i, j, k, l, m, n, i, j, l, k, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(i, j, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim2, i, j, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, i, j, k, l, m, n, i, j, k, l, Dim4, Dim5, Dim0, Dim1, Dim2, Dim3, m, n, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim4, Dim5, m, n>(TensorExpr(a, b));
}

//B(i, j, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim0, m, k, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim1, Dim0, i, j, k, l, m, n, m, k, j, i, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(m, k, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim0, m, k, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim0, k, m, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim1, Dim0, i, j, k, l, m, n, k, m, j, i, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(k, m, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim0, k, m, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim0, m, j, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim2, Dim0, i, j, k, l, m, n, m, j, k, i, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(m, j, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim0, m, j, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim0, k, j, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim4, Dim0, i, j, k, l, m, n, k, j, m, i, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(k, j, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim0, k, j, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim0, j, m, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim2, Dim0, i, j, k, l, m, n, j, m, k, i, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(j, m, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim0, j, m, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim0, j, k, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim4, Dim0, i, j, k, l, m, n, j, k, m, i, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(j, k, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim0, j, k, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim1, m, k, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim0, Dim1, i, j, k, l, m, n, m, k, i, j, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(m, k, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim1, m, k, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim1, k, m, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim0, Dim1, i, j, k, l, m, n, k, m, i, j, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(k, m, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim1, k, m, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim2, m, j, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim0, Dim2, i, j, k, l, m, n, m, j, i, k, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(m, j, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim2, m, j, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim4, k, j, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim0, Dim4, i, j, k, l, m, n, k, j, i, m, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(k, j, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim4, k, j, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim2, j, m, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim0, Dim2, i, j, k, l, m, n, j, m, i, k, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(j, m, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim2, j, m, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim4, j, k, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim0, Dim4, i, j, k, l, m, n, j, k, i, m, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(j, k, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim4, j, k, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim1, m, i, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim2, Dim1, i, j, k, l, m, n, m, i, k, j, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(m, i, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim1, m, i, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim1, k, i, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim4, Dim1, i, j, k, l, m, n, k, i, m, j, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(k, i, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim1, k, i, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim2, m, i, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim1, Dim2, i, j, k, l, m, n, m, i, j, k, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(m, i, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim2, m, i, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim4, k, i, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim1, Dim4, i, j, k, l, m, n, k, i, j, m, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(k, i, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim4, k, i, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim2, j, i, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim4, Dim2, i, j, k, l, m, n, j, i, m, k, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(j, i, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim2, j, i, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim4, j, i, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim2, Dim4, i, j, k, l, m, n, j, i, k, m, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(j, i, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim4, j, i, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim1, i, m, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim2, Dim1, i, j, k, l, m, n, i, m, k, j, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(i, m, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim1, i, m, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim1, i, k, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim4, Dim1, i, j, k, l, m, n, i, k, m, j, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(i, k, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim1, i, k, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim2, i, m, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim1, Dim2, i, j, k, l, m, n, i, m, j, k, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(i, m, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim2, i, m, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim4, i, k, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim1, Dim4, i, j, k, l, m, n, i, k, j, m, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(i, k, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim4, i, k, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim2, i, j, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim4, Dim2, i, j, k, l, m, n, i, j, m, k, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(i, j, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim2, i, j, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim4, i, j, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim4, i, j, k, l, m, n, i, j, k, m, Dim3, Dim5, Dim0, Dim1, Dim2, Dim4, l, n, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim5, l, n>(TensorExpr(a, b));
}

//B(i, j, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim4, i, j, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim0, n, k, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim1, Dim0, i, j, k, l, m, n, n, k, j, i, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(n, k, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim0, n, k, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim0, k, n, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim1, Dim0, i, j, k, l, m, n, k, n, j, i, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(k, n, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim0, k, n, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim0, n, j, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim2, Dim0, i, j, k, l, m, n, n, j, k, i, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(n, j, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim0, n, j, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim0, k, j, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim5, Dim0, i, j, k, l, m, n, k, j, n, i, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(k, j, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim0, k, j, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim0, j, n, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim2, Dim0, i, j, k, l, m, n, j, n, k, i, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(j, n, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim0, j, n, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim0, j, k, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim5, Dim0, i, j, k, l, m, n, j, k, n, i, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(j, k, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim0, j, k, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim1, n, k, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim0, Dim1, i, j, k, l, m, n, n, k, i, j, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(n, k, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim1, n, k, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim1, k, n, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim0, Dim1, i, j, k, l, m, n, k, n, i, j, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(k, n, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim1, k, n, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim2, n, j, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim0, Dim2, i, j, k, l, m, n, n, j, i, k, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(n, j, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim2, n, j, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim5, k, j, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim0, Dim5, i, j, k, l, m, n, k, j, i, n, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(k, j, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim5, k, j, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim2, j, n, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim0, Dim2, i, j, k, l, m, n, j, n, i, k, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(j, n, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim2, j, n, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim5, j, k, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim0, Dim5, i, j, k, l, m, n, j, k, i, n, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(j, k, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim5, j, k, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim1, n, i, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim2, Dim1, i, j, k, l, m, n, n, i, k, j, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(n, i, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim1, n, i, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim1, k, i, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim5, Dim1, i, j, k, l, m, n, k, i, n, j, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(k, i, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim1, k, i, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim2, n, i, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim1, Dim2, i, j, k, l, m, n, n, i, j, k, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(n, i, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim2, n, i, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim5, k, i, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim1, Dim5, i, j, k, l, m, n, k, i, j, n, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(k, i, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim5, k, i, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim2, j, i, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim5, Dim2, i, j, k, l, m, n, j, i, n, k, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(j, i, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim2, j, i, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim5, j, i, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim2, Dim5, i, j, k, l, m, n, j, i, k, n, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(j, i, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim5, j, i, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim1, i, n, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim2, Dim1, i, j, k, l, m, n, i, n, k, j, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(i, n, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim1, i, n, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim1, i, k, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim5, Dim1, i, j, k, l, m, n, i, k, n, j, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(i, k, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim1, i, k, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim2, i, n, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim1, Dim2, i, j, k, l, m, n, i, n, j, k, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(i, n, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim2, i, n, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim5, i, k, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim1, Dim5, i, j, k, l, m, n, i, k, j, n, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(i, k, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim5, i, k, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim2, i, j, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim5, Dim2, i, j, k, l, m, n, i, j, n, k, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(i, j, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim2, i, j, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim5, i, j, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, Dim5, i, j, k, l, m, n, i, j, k, n, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//B(i, j, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim5, i, j, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim0, m, l, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim1, Dim0, i, j, k, l, m, n, m, l, j, i, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(m, l, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim0, m, l, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim0, l, m, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim1, Dim0, i, j, k, l, m, n, l, m, j, i, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(l, m, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim0, l, m, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim0, m, j, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim3, Dim0, i, j, k, l, m, n, m, j, l, i, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(m, j, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim0, m, j, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim0, l, j, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim4, Dim0, i, j, k, l, m, n, l, j, m, i, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(l, j, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim0, l, j, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim0, j, m, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim3, Dim0, i, j, k, l, m, n, j, m, l, i, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(j, m, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim0, j, m, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim0, j, l, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim4, Dim0, i, j, k, l, m, n, j, l, m, i, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(j, l, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim0, j, l, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim1, m, l, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim0, Dim1, i, j, k, l, m, n, m, l, i, j, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(m, l, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim1, m, l, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim1, l, m, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim0, Dim1, i, j, k, l, m, n, l, m, i, j, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(l, m, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim1, l, m, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim3, m, j, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim0, Dim3, i, j, k, l, m, n, m, j, i, l, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(m, j, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim3, m, j, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim4, l, j, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim0, Dim4, i, j, k, l, m, n, l, j, i, m, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(l, j, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim4, l, j, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim3, j, m, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim0, Dim3, i, j, k, l, m, n, j, m, i, l, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(j, m, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim3, j, m, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim4, j, l, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim0, Dim4, i, j, k, l, m, n, j, l, i, m, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(j, l, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim4, j, l, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim1, m, i, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim3, Dim1, i, j, k, l, m, n, m, i, l, j, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(m, i, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim1, m, i, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim1, l, i, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim4, Dim1, i, j, k, l, m, n, l, i, m, j, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(l, i, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim1, l, i, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim3, m, i, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim1, Dim3, i, j, k, l, m, n, m, i, j, l, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(m, i, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim3, m, i, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim4, l, i, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim1, Dim4, i, j, k, l, m, n, l, i, j, m, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(l, i, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim4, l, i, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim3, j, i, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim4, Dim3, i, j, k, l, m, n, j, i, m, l, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(j, i, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim3, j, i, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim4, j, i, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim3, Dim4, i, j, k, l, m, n, j, i, l, m, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(j, i, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim4, j, i, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim1, i, m, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim3, Dim1, i, j, k, l, m, n, i, m, l, j, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(i, m, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim1, i, m, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim1, i, l, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim4, Dim1, i, j, k, l, m, n, i, l, m, j, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(i, l, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim1, i, l, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim3, i, m, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim1, Dim3, i, j, k, l, m, n, i, m, j, l, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(i, m, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim3, i, m, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim4, i, l, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim1, Dim4, i, j, k, l, m, n, i, l, j, m, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(i, l, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim4, i, l, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim3, i, j, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim4, Dim3, i, j, k, l, m, n, i, j, m, l, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(i, j, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim3, i, j, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim4, i, j, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim3, Dim4, i, j, k, l, m, n, i, j, l, m, Dim2, Dim5, Dim0, Dim1, Dim3, Dim4, k, n, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim5, k, n>(TensorExpr(a, b));
}

//B(i, j, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim4, i, j, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim0, n, l, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim1, Dim0, i, j, k, l, m, n, n, l, j, i, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(n, l, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim0, n, l, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim0, l, n, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim1, Dim0, i, j, k, l, m, n, l, n, j, i, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(l, n, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim0, l, n, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim0, n, j, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim3, Dim0, i, j, k, l, m, n, n, j, l, i, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(n, j, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim0, n, j, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim0, l, j, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim5, Dim0, i, j, k, l, m, n, l, j, n, i, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(l, j, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim0, l, j, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim0, j, n, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim3, Dim0, i, j, k, l, m, n, j, n, l, i, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(j, n, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim0, j, n, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim0, j, l, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim5, Dim0, i, j, k, l, m, n, j, l, n, i, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(j, l, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim0, j, l, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim1, n, l, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim0, Dim1, i, j, k, l, m, n, n, l, i, j, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(n, l, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim1, n, l, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim1, l, n, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim0, Dim1, i, j, k, l, m, n, l, n, i, j, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(l, n, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim1, l, n, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim3, n, j, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim0, Dim3, i, j, k, l, m, n, n, j, i, l, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(n, j, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim3, n, j, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim5, l, j, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim0, Dim5, i, j, k, l, m, n, l, j, i, n, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(l, j, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim5, l, j, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim3, j, n, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim0, Dim3, i, j, k, l, m, n, j, n, i, l, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(j, n, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim3, j, n, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim5, j, l, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim0, Dim5, i, j, k, l, m, n, j, l, i, n, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(j, l, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim5, j, l, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim1, n, i, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim3, Dim1, i, j, k, l, m, n, n, i, l, j, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(n, i, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim1, n, i, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim1, l, i, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim5, Dim1, i, j, k, l, m, n, l, i, n, j, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(l, i, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim1, l, i, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim3, n, i, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim1, Dim3, i, j, k, l, m, n, n, i, j, l, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(n, i, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim3, n, i, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim5, l, i, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim1, Dim5, i, j, k, l, m, n, l, i, j, n, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(l, i, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim5, l, i, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim3, j, i, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim5, Dim3, i, j, k, l, m, n, j, i, n, l, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(j, i, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim3, j, i, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim5, j, i, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim3, Dim5, i, j, k, l, m, n, j, i, l, n, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(j, i, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim5, j, i, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim1, i, n, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim3, Dim1, i, j, k, l, m, n, i, n, l, j, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(i, n, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim1, i, n, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim1, i, l, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim5, Dim1, i, j, k, l, m, n, i, l, n, j, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(i, l, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim1, i, l, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim3, i, n, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim1, Dim3, i, j, k, l, m, n, i, n, j, l, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(i, n, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim3, i, n, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim5, i, l, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim1, Dim5, i, j, k, l, m, n, i, l, j, n, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(i, l, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim5, i, l, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim3, i, j, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim5, Dim3, i, j, k, l, m, n, i, j, n, l, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(i, j, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim3, i, j, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim5, i, j, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim3, Dim5, i, j, k, l, m, n, i, j, l, n, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//B(i, j, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim5, i, j, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim0, n, m, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim1, Dim0, i, j, k, l, m, n, n, m, j, i, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(n, m, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim0, n, m, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim0, m, n, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim1, Dim0, i, j, k, l, m, n, m, n, j, i, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(m, n, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim0, m, n, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim0, n, j, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim4, Dim0, i, j, k, l, m, n, n, j, m, i, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(n, j, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim0, n, j, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim0, m, j, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim5, Dim0, i, j, k, l, m, n, m, j, n, i, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(m, j, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim0, m, j, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim0, j, n, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim4, Dim0, i, j, k, l, m, n, j, n, m, i, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(j, n, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim0, j, n, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim0, j, m, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim5, Dim0, i, j, k, l, m, n, j, m, n, i, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(j, m, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim0, j, m, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim1, n, m, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim0, Dim1, i, j, k, l, m, n, n, m, i, j, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(n, m, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim1, n, m, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim1, m, n, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim0, Dim1, i, j, k, l, m, n, m, n, i, j, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(m, n, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim1, m, n, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim4, n, j, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim0, Dim4, i, j, k, l, m, n, n, j, i, m, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(n, j, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim4, n, j, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim5, m, j, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim0, Dim5, i, j, k, l, m, n, m, j, i, n, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(m, j, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim5, m, j, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim4, j, n, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim0, Dim4, i, j, k, l, m, n, j, n, i, m, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(j, n, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim4, j, n, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim5, j, m, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim0, Dim5, i, j, k, l, m, n, j, m, i, n, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(j, m, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim5, j, m, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim1, n, i, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim4, Dim1, i, j, k, l, m, n, n, i, m, j, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(n, i, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim1, n, i, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim1, m, i, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim5, Dim1, i, j, k, l, m, n, m, i, n, j, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(m, i, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim1, m, i, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim4, n, i, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim1, Dim4, i, j, k, l, m, n, n, i, j, m, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(n, i, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim4, n, i, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim5, m, i, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim1, Dim5, i, j, k, l, m, n, m, i, j, n, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(m, i, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim5, m, i, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim4, j, i, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim5, Dim4, i, j, k, l, m, n, j, i, n, m, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(j, i, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim4, j, i, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim5, j, i, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim4, Dim5, i, j, k, l, m, n, j, i, m, n, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(j, i, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim5, j, i, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim1, i, n, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim4, Dim1, i, j, k, l, m, n, i, n, m, j, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(i, n, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim1, i, n, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim1, i, m, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim5, Dim1, i, j, k, l, m, n, i, m, n, j, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(i, m, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim1, i, m, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim4, i, n, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim1, Dim4, i, j, k, l, m, n, i, n, j, m, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(i, n, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim4, i, n, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim5, i, m, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim1, Dim5, i, j, k, l, m, n, i, m, j, n, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(i, m, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim5, i, m, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim4, i, j, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim5, Dim4, i, j, k, l, m, n, i, j, n, m, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(i, j, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim4, i, j, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim5, i, j, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim4, Dim5, i, j, k, l, m, n, i, j, m, n, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//B(i, j, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim5, i, j, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim0, m, l, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim2, Dim0, i, j, k, l, m, n, m, l, k, i, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(m, l, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim0, m, l, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim0, l, m, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim2, Dim0, i, j, k, l, m, n, l, m, k, i, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(l, m, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim0, l, m, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim0, m, k, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim3, Dim0, i, j, k, l, m, n, m, k, l, i, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(m, k, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim0, m, k, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim0, l, k, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim4, Dim0, i, j, k, l, m, n, l, k, m, i, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(l, k, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim0, l, k, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim0, k, m, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim3, Dim0, i, j, k, l, m, n, k, m, l, i, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(k, m, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim0, k, m, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim0, k, l, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim4, Dim0, i, j, k, l, m, n, k, l, m, i, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(k, l, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim0, k, l, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim2, m, l, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim0, Dim2, i, j, k, l, m, n, m, l, i, k, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(m, l, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim2, m, l, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim2, l, m, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim0, Dim2, i, j, k, l, m, n, l, m, i, k, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(l, m, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim2, l, m, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim3, m, k, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim0, Dim3, i, j, k, l, m, n, m, k, i, l, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(m, k, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim3, m, k, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim4, l, k, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim0, Dim4, i, j, k, l, m, n, l, k, i, m, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(l, k, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim4, l, k, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim3, k, m, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim0, Dim3, i, j, k, l, m, n, k, m, i, l, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(k, m, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim3, k, m, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim4, k, l, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim0, Dim4, i, j, k, l, m, n, k, l, i, m, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(k, l, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim4, k, l, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim2, m, i, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim3, Dim2, i, j, k, l, m, n, m, i, l, k, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(m, i, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim2, m, i, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim2, l, i, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim4, Dim2, i, j, k, l, m, n, l, i, m, k, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(l, i, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim2, l, i, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim3, m, i, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim2, Dim3, i, j, k, l, m, n, m, i, k, l, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(m, i, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim3, m, i, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim4, l, i, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim2, Dim4, i, j, k, l, m, n, l, i, k, m, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(l, i, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim4, l, i, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim3, k, i, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim4, Dim3, i, j, k, l, m, n, k, i, m, l, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(k, i, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim3, k, i, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim4, k, i, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim3, Dim4, i, j, k, l, m, n, k, i, l, m, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(k, i, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim4, k, i, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim2, i, m, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim3, Dim2, i, j, k, l, m, n, i, m, l, k, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(i, m, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim2, i, m, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim2, i, l, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim4, Dim2, i, j, k, l, m, n, i, l, m, k, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(i, l, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim2, i, l, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim3, i, m, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim2, Dim3, i, j, k, l, m, n, i, m, k, l, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(i, m, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim3, i, m, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim4, i, l, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim2, Dim4, i, j, k, l, m, n, i, l, k, m, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(i, l, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim4, i, l, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim3, i, k, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim4, Dim3, i, j, k, l, m, n, i, k, m, l, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(i, k, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim3, i, k, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim4, i, k, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim3, Dim4, i, j, k, l, m, n, i, k, l, m, Dim1, Dim5, Dim0, Dim2, Dim3, Dim4, j, n, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim5, j, n>(TensorExpr(a, b));
}

//B(i, k, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim4, i, k, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim0, n, l, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim2, Dim0, i, j, k, l, m, n, n, l, k, i, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(n, l, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim0, n, l, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim0, l, n, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim2, Dim0, i, j, k, l, m, n, l, n, k, i, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(l, n, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim0, l, n, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim0, n, k, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim3, Dim0, i, j, k, l, m, n, n, k, l, i, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(n, k, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim0, n, k, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim0, l, k, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim5, Dim0, i, j, k, l, m, n, l, k, n, i, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(l, k, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim0, l, k, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim0, k, n, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim3, Dim0, i, j, k, l, m, n, k, n, l, i, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(k, n, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim0, k, n, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim0, k, l, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim5, Dim0, i, j, k, l, m, n, k, l, n, i, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(k, l, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim0, k, l, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim2, n, l, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim0, Dim2, i, j, k, l, m, n, n, l, i, k, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(n, l, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim2, n, l, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim2, l, n, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim0, Dim2, i, j, k, l, m, n, l, n, i, k, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(l, n, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim2, l, n, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim3, n, k, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim0, Dim3, i, j, k, l, m, n, n, k, i, l, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(n, k, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim3, n, k, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim5, l, k, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim0, Dim5, i, j, k, l, m, n, l, k, i, n, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(l, k, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim5, l, k, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim3, k, n, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim0, Dim3, i, j, k, l, m, n, k, n, i, l, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(k, n, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim3, k, n, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim5, k, l, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim0, Dim5, i, j, k, l, m, n, k, l, i, n, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(k, l, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim5, k, l, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim2, n, i, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim3, Dim2, i, j, k, l, m, n, n, i, l, k, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(n, i, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim2, n, i, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim2, l, i, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim5, Dim2, i, j, k, l, m, n, l, i, n, k, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(l, i, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim2, l, i, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim3, n, i, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim2, Dim3, i, j, k, l, m, n, n, i, k, l, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(n, i, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim3, n, i, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim5, l, i, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim2, Dim5, i, j, k, l, m, n, l, i, k, n, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(l, i, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim5, l, i, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim3, k, i, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim5, Dim3, i, j, k, l, m, n, k, i, n, l, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(k, i, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim3, k, i, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim5, k, i, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim3, Dim5, i, j, k, l, m, n, k, i, l, n, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(k, i, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim5, k, i, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim2, i, n, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim3, Dim2, i, j, k, l, m, n, i, n, l, k, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(i, n, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim2, i, n, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim2, i, l, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim5, Dim2, i, j, k, l, m, n, i, l, n, k, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(i, l, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim2, i, l, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim3, i, n, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim2, Dim3, i, j, k, l, m, n, i, n, k, l, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(i, n, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim3, i, n, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim5, i, l, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim2, Dim5, i, j, k, l, m, n, i, l, k, n, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(i, l, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim5, i, l, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim3, i, k, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim5, Dim3, i, j, k, l, m, n, i, k, n, l, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(i, k, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim3, i, k, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim5, i, k, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim3, Dim5, i, j, k, l, m, n, i, k, l, n, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//B(i, k, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim5, i, k, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim0, n, m, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim2, Dim0, i, j, k, l, m, n, n, m, k, i, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(n, m, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim0, n, m, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim0, m, n, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim2, Dim0, i, j, k, l, m, n, m, n, k, i, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(m, n, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim0, m, n, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim0, n, k, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim4, Dim0, i, j, k, l, m, n, n, k, m, i, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(n, k, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim0, n, k, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim0, m, k, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim5, Dim0, i, j, k, l, m, n, m, k, n, i, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(m, k, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim0, m, k, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim0, k, n, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim4, Dim0, i, j, k, l, m, n, k, n, m, i, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(k, n, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim0, k, n, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim0, k, m, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim5, Dim0, i, j, k, l, m, n, k, m, n, i, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(k, m, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim0, k, m, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim2, n, m, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim0, Dim2, i, j, k, l, m, n, n, m, i, k, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(n, m, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim2, n, m, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim2, m, n, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim0, Dim2, i, j, k, l, m, n, m, n, i, k, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(m, n, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim2, m, n, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim4, n, k, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim0, Dim4, i, j, k, l, m, n, n, k, i, m, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(n, k, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim4, n, k, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim5, m, k, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim0, Dim5, i, j, k, l, m, n, m, k, i, n, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(m, k, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim5, m, k, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim4, k, n, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim0, Dim4, i, j, k, l, m, n, k, n, i, m, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(k, n, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim4, k, n, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim5, k, m, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim0, Dim5, i, j, k, l, m, n, k, m, i, n, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(k, m, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim5, k, m, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim2, n, i, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim4, Dim2, i, j, k, l, m, n, n, i, m, k, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(n, i, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim2, n, i, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim2, m, i, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim5, Dim2, i, j, k, l, m, n, m, i, n, k, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(m, i, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim2, m, i, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim4, n, i, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim2, Dim4, i, j, k, l, m, n, n, i, k, m, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(n, i, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim4, n, i, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim5, m, i, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim2, Dim5, i, j, k, l, m, n, m, i, k, n, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(m, i, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim5, m, i, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim4, k, i, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim5, Dim4, i, j, k, l, m, n, k, i, n, m, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(k, i, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim4, k, i, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim5, k, i, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim4, Dim5, i, j, k, l, m, n, k, i, m, n, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(k, i, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim5, k, i, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim2, i, n, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim4, Dim2, i, j, k, l, m, n, i, n, m, k, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(i, n, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim2, i, n, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim2, i, m, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim5, Dim2, i, j, k, l, m, n, i, m, n, k, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(i, m, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim2, i, m, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim4, i, n, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim2, Dim4, i, j, k, l, m, n, i, n, k, m, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(i, n, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim4, i, n, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim5, i, m, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim2, Dim5, i, j, k, l, m, n, i, m, k, n, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(i, m, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim5, i, m, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim4, i, k, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim5, Dim4, i, j, k, l, m, n, i, k, n, m, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(i, k, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim4, i, k, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim5, i, k, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim4, Dim5, i, j, k, l, m, n, i, k, m, n, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//B(i, k, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim5, i, k, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim0, n, m, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim3, Dim0, i, j, k, l, m, n, n, m, l, i, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(n, m, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim0, n, m, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim0, m, n, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim3, Dim0, i, j, k, l, m, n, m, n, l, i, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(m, n, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim0, m, n, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim0, n, l, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim4, Dim0, i, j, k, l, m, n, n, l, m, i, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(n, l, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim0, n, l, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim0, m, l, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim5, Dim0, i, j, k, l, m, n, m, l, n, i, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(m, l, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim0, m, l, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim0, l, n, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim4, Dim0, i, j, k, l, m, n, l, n, m, i, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(l, n, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim0, l, n, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim0, l, m, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim5, Dim0, i, j, k, l, m, n, l, m, n, i, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(l, m, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim0, l, m, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim3, n, m, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim0, Dim3, i, j, k, l, m, n, n, m, i, l, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(n, m, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim3, n, m, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim3, m, n, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim0, Dim3, i, j, k, l, m, n, m, n, i, l, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(m, n, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim3, m, n, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim4, n, l, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim0, Dim4, i, j, k, l, m, n, n, l, i, m, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(n, l, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim4, n, l, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim5, m, l, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim0, Dim5, i, j, k, l, m, n, m, l, i, n, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(m, l, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim5, m, l, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim4, l, n, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim0, Dim4, i, j, k, l, m, n, l, n, i, m, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(l, n, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim4, l, n, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim5, l, m, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim0, Dim5, i, j, k, l, m, n, l, m, i, n, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(l, m, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim5, l, m, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim3, n, i, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim4, Dim3, i, j, k, l, m, n, n, i, m, l, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(n, i, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim3, n, i, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim3, m, i, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim5, Dim3, i, j, k, l, m, n, m, i, n, l, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(m, i, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim3, m, i, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim4, n, i, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim3, Dim4, i, j, k, l, m, n, n, i, l, m, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(n, i, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim4, n, i, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim5, m, i, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim3, Dim5, i, j, k, l, m, n, m, i, l, n, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(m, i, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim5, m, i, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim4, l, i, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim5, Dim4, i, j, k, l, m, n, l, i, n, m, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(l, i, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim4, l, i, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim5, l, i, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim4, Dim5, i, j, k, l, m, n, l, i, m, n, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(l, i, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim5, l, i, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim3, i, n, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim4, Dim3, i, j, k, l, m, n, i, n, m, l, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(i, n, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim3, i, n, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim3, i, m, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim5, Dim3, i, j, k, l, m, n, i, m, n, l, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(i, m, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim3, i, m, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim4, i, n, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim3, Dim4, i, j, k, l, m, n, i, n, l, m, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(i, n, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim4, i, n, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim5, i, m, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim3, Dim5, i, j, k, l, m, n, i, m, l, n, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(i, m, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim5, i, m, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim4, i, l, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim5, Dim4, i, j, k, l, m, n, i, l, n, m, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(i, l, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim4, i, l, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim5, i, l, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim4, Dim5, i, j, k, l, m, n, i, l, m, n, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, j, k>(TensorExpr(a, b));
}

//B(i, l, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim5, i, l, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim1, m, l, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim2, Dim1, i, j, k, l, m, n, m, l, k, j, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(m, l, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim1, m, l, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim1, l, m, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim2, Dim1, i, j, k, l, m, n, l, m, k, j, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(l, m, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim1, l, m, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim1, m, k, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim3, Dim1, i, j, k, l, m, n, m, k, l, j, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(m, k, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim1, m, k, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim1, l, k, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim4, Dim1, i, j, k, l, m, n, l, k, m, j, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(l, k, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim1, l, k, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim1, k, m, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim3, Dim1, i, j, k, l, m, n, k, m, l, j, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(k, m, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim1, k, m, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim1, k, l, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim4, Dim1, i, j, k, l, m, n, k, l, m, j, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(k, l, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim1, k, l, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim2, m, l, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim1, Dim2, i, j, k, l, m, n, m, l, j, k, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(m, l, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim2, m, l, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim2, l, m, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim1, Dim2, i, j, k, l, m, n, l, m, j, k, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(l, m, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim2, l, m, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim3, m, k, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim1, Dim3, i, j, k, l, m, n, m, k, j, l, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(m, k, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim3, m, k, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim4, l, k, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim1, Dim4, i, j, k, l, m, n, l, k, j, m, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(l, k, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim4, l, k, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim3, k, m, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim1, Dim3, i, j, k, l, m, n, k, m, j, l, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(k, m, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim3, k, m, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim4, k, l, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim1, Dim4, i, j, k, l, m, n, k, l, j, m, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(k, l, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim4, k, l, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim2, m, j, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim3, Dim2, i, j, k, l, m, n, m, j, l, k, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(m, j, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim2, m, j, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim2, l, j, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim4, Dim2, i, j, k, l, m, n, l, j, m, k, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(l, j, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim2, l, j, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim3, m, j, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim2, Dim3, i, j, k, l, m, n, m, j, k, l, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(m, j, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim3, m, j, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim4, l, j, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim2, Dim4, i, j, k, l, m, n, l, j, k, m, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(l, j, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim4, l, j, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim3, k, j, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim4, Dim3, i, j, k, l, m, n, k, j, m, l, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(k, j, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim3, k, j, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim4, k, j, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim3, Dim4, i, j, k, l, m, n, k, j, l, m, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(k, j, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim4, k, j, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim2, j, m, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim3, Dim2, i, j, k, l, m, n, j, m, l, k, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(j, m, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim2, j, m, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim2, j, l, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim4, Dim2, i, j, k, l, m, n, j, l, m, k, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(j, l, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim2, j, l, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim3, j, m, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim2, Dim3, i, j, k, l, m, n, j, m, k, l, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(j, m, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim3, j, m, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim4, j, l, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim2, Dim4, i, j, k, l, m, n, j, l, k, m, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(j, l, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim4, j, l, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim3, j, k, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim4, Dim3, i, j, k, l, m, n, j, k, m, l, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(j, k, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim3, j, k, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim4, j, k, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m, n, j, k, l, m, Dim0, Dim5, Dim1, Dim2, Dim3, Dim4, i, n, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim5, i, n>(TensorExpr(a, b));
}

//B(j, k, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim4, j, k, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim1, n, l, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim2, Dim1, i, j, k, l, m, n, n, l, k, j, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(n, l, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim1, n, l, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim1, l, n, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim2, Dim1, i, j, k, l, m, n, l, n, k, j, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(l, n, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim1, l, n, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim1, n, k, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim3, Dim1, i, j, k, l, m, n, n, k, l, j, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(n, k, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim1, n, k, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim1, l, k, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim5, Dim1, i, j, k, l, m, n, l, k, n, j, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(l, k, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim1, l, k, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim1, k, n, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim3, Dim1, i, j, k, l, m, n, k, n, l, j, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(k, n, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim1, k, n, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim1, k, l, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim5, Dim1, i, j, k, l, m, n, k, l, n, j, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(k, l, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim1, k, l, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim2, n, l, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim1, Dim2, i, j, k, l, m, n, n, l, j, k, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(n, l, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim2, n, l, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim2, l, n, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim1, Dim2, i, j, k, l, m, n, l, n, j, k, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(l, n, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim2, l, n, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim3, n, k, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim1, Dim3, i, j, k, l, m, n, n, k, j, l, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(n, k, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim3, n, k, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim5, l, k, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim1, Dim5, i, j, k, l, m, n, l, k, j, n, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(l, k, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim5, l, k, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim3, k, n, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim1, Dim3, i, j, k, l, m, n, k, n, j, l, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(k, n, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim3, k, n, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim5, k, l, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim1, Dim5, i, j, k, l, m, n, k, l, j, n, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(k, l, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim5, k, l, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim2, n, j, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim3, Dim2, i, j, k, l, m, n, n, j, l, k, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(n, j, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim2, n, j, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim2, l, j, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim5, Dim2, i, j, k, l, m, n, l, j, n, k, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(l, j, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim2, l, j, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim3, n, j, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim2, Dim3, i, j, k, l, m, n, n, j, k, l, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(n, j, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim3, n, j, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim5, l, j, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim2, Dim5, i, j, k, l, m, n, l, j, k, n, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(l, j, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim5, l, j, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim3, k, j, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim5, Dim3, i, j, k, l, m, n, k, j, n, l, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(k, j, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim3, k, j, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim5, k, j, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim3, Dim5, i, j, k, l, m, n, k, j, l, n, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(k, j, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim5, k, j, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim2, j, n, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim3, Dim2, i, j, k, l, m, n, j, n, l, k, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(j, n, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim2, j, n, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim2, j, l, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim5, Dim2, i, j, k, l, m, n, j, l, n, k, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(j, l, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim2, j, l, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim3, j, n, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim2, Dim3, i, j, k, l, m, n, j, n, k, l, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(j, n, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim3, j, n, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim5, j, l, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim2, Dim5, i, j, k, l, m, n, j, l, k, n, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(j, l, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim5, j, l, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim3, j, k, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim5, Dim3, i, j, k, l, m, n, j, k, n, l, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(j, k, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim3, j, k, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim5, j, k, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim3, Dim5, i, j, k, l, m, n, j, k, l, n, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//B(j, k, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim5, j, k, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim1, n, m, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim2, Dim1, i, j, k, l, m, n, n, m, k, j, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(n, m, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim1, n, m, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim1, m, n, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim2, Dim1, i, j, k, l, m, n, m, n, k, j, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(m, n, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim1, m, n, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim1, n, k, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim4, Dim1, i, j, k, l, m, n, n, k, m, j, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(n, k, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim1, n, k, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim1, m, k, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim5, Dim1, i, j, k, l, m, n, m, k, n, j, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(m, k, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim1, m, k, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim1, k, n, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim4, Dim1, i, j, k, l, m, n, k, n, m, j, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(k, n, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim1, k, n, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim1, k, m, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim5, Dim1, i, j, k, l, m, n, k, m, n, j, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(k, m, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim1, k, m, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim2, n, m, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim1, Dim2, i, j, k, l, m, n, n, m, j, k, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(n, m, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim2, n, m, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim2, m, n, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim1, Dim2, i, j, k, l, m, n, m, n, j, k, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(m, n, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim2, m, n, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim4, n, k, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim1, Dim4, i, j, k, l, m, n, n, k, j, m, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(n, k, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim4, n, k, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim5, m, k, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim1, Dim5, i, j, k, l, m, n, m, k, j, n, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(m, k, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim5, m, k, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim4, k, n, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim1, Dim4, i, j, k, l, m, n, k, n, j, m, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(k, n, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim4, k, n, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim5, k, m, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim1, Dim5, i, j, k, l, m, n, k, m, j, n, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(k, m, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim5, k, m, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim2, n, j, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim4, Dim2, i, j, k, l, m, n, n, j, m, k, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(n, j, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim2, n, j, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim2, m, j, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim5, Dim2, i, j, k, l, m, n, m, j, n, k, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(m, j, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim2, m, j, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim4, n, j, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim2, Dim4, i, j, k, l, m, n, n, j, k, m, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(n, j, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim4, n, j, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim5, m, j, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim2, Dim5, i, j, k, l, m, n, m, j, k, n, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(m, j, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim5, m, j, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim4, k, j, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim5, Dim4, i, j, k, l, m, n, k, j, n, m, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(k, j, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim4, k, j, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim5, k, j, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim4, Dim5, i, j, k, l, m, n, k, j, m, n, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(k, j, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim5, k, j, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim2, j, n, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim4, Dim2, i, j, k, l, m, n, j, n, m, k, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(j, n, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim2, j, n, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim2, j, m, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim5, Dim2, i, j, k, l, m, n, j, m, n, k, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(j, m, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim2, j, m, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim4, j, n, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim2, Dim4, i, j, k, l, m, n, j, n, k, m, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(j, n, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim4, j, n, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim5, j, m, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim2, Dim5, i, j, k, l, m, n, j, m, k, n, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(j, m, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim5, j, m, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim4, j, k, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim5, Dim4, i, j, k, l, m, n, j, k, n, m, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(j, k, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim4, j, k, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim5, j, k, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim4, Dim5, i, j, k, l, m, n, j, k, m, n, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//B(j, k, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim5, j, k, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim1, n, m, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim3, Dim1, i, j, k, l, m, n, n, m, l, j, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(n, m, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim1, n, m, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim1, m, n, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim3, Dim1, i, j, k, l, m, n, m, n, l, j, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(m, n, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim1, m, n, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim1, n, l, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim4, Dim1, i, j, k, l, m, n, n, l, m, j, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(n, l, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim1, n, l, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim1, m, l, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim5, Dim1, i, j, k, l, m, n, m, l, n, j, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(m, l, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim1, m, l, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim1, l, n, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim4, Dim1, i, j, k, l, m, n, l, n, m, j, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(l, n, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim1, l, n, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim1, l, m, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim5, Dim1, i, j, k, l, m, n, l, m, n, j, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(l, m, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim1, l, m, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim3, n, m, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim1, Dim3, i, j, k, l, m, n, n, m, j, l, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(n, m, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim3, n, m, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim3, m, n, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim1, Dim3, i, j, k, l, m, n, m, n, j, l, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(m, n, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim3, m, n, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim4, n, l, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim1, Dim4, i, j, k, l, m, n, n, l, j, m, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(n, l, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim4, n, l, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim5, m, l, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim1, Dim5, i, j, k, l, m, n, m, l, j, n, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(m, l, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim5, m, l, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim4, l, n, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim1, Dim4, i, j, k, l, m, n, l, n, j, m, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(l, n, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim4, l, n, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim5, l, m, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim1, Dim5, i, j, k, l, m, n, l, m, j, n, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(l, m, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim5, l, m, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim3, n, j, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim4, Dim3, i, j, k, l, m, n, n, j, m, l, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(n, j, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim3, n, j, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim3, m, j, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim5, Dim3, i, j, k, l, m, n, m, j, n, l, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(m, j, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim3, m, j, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim4, n, j, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim3, Dim4, i, j, k, l, m, n, n, j, l, m, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(n, j, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim4, n, j, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim5, m, j, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim3, Dim5, i, j, k, l, m, n, m, j, l, n, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(m, j, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim5, m, j, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim4, l, j, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim5, Dim4, i, j, k, l, m, n, l, j, n, m, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(l, j, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim4, l, j, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim5, l, j, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim4, Dim5, i, j, k, l, m, n, l, j, m, n, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(l, j, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim5, l, j, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim3, j, n, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim4, Dim3, i, j, k, l, m, n, j, n, m, l, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(j, n, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim3, j, n, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim3, j, m, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim5, Dim3, i, j, k, l, m, n, j, m, n, l, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(j, m, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim3, j, m, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim4, j, n, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim3, Dim4, i, j, k, l, m, n, j, n, l, m, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(j, n, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim4, j, n, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim5, j, m, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim3, Dim5, i, j, k, l, m, n, j, m, l, n, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(j, m, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim5, j, m, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim4, j, l, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim5, Dim4, i, j, k, l, m, n, j, l, n, m, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(j, l, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim4, j, l, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim5, j, l, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim4, Dim5, i, j, k, l, m, n, j, l, m, n, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, i, k>(TensorExpr(a, b));
}

//B(j, l, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim5, j, l, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim2, n, m, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim3, Dim2, i, j, k, l, m, n, n, m, l, k, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(n, m, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim2, n, m, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim2, m, n, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim3, Dim2, i, j, k, l, m, n, m, n, l, k, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(m, n, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim2, m, n, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim2, n, l, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim4, Dim2, i, j, k, l, m, n, n, l, m, k, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(n, l, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim2, n, l, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim2, m, l, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim5, Dim2, i, j, k, l, m, n, m, l, n, k, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(m, l, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim2, m, l, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim2, l, n, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim4, Dim2, i, j, k, l, m, n, l, n, m, k, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(l, n, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim2, l, n, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim2, l, m, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim5, Dim2, i, j, k, l, m, n, l, m, n, k, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(l, m, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim2, l, m, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim3, n, m, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim2, Dim3, i, j, k, l, m, n, n, m, k, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(n, m, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim3, n, m, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim3, m, n, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim2, Dim3, i, j, k, l, m, n, m, n, k, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(m, n, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim3, m, n, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim4, n, l, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim2, Dim4, i, j, k, l, m, n, n, l, k, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(n, l, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim4, n, l, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim5, m, l, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim2, Dim5, i, j, k, l, m, n, m, l, k, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(m, l, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim5, m, l, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim4, l, n, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim2, Dim4, i, j, k, l, m, n, l, n, k, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(l, n, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim4, l, n, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim5, l, m, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim2, Dim5, i, j, k, l, m, n, l, m, k, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(l, m, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim5, l, m, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim3, n, k, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim4, Dim3, i, j, k, l, m, n, n, k, m, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(n, k, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim3, n, k, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim3, m, k, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim5, Dim3, i, j, k, l, m, n, m, k, n, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(m, k, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim3, m, k, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim4, n, k, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim3, Dim4, i, j, k, l, m, n, n, k, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(n, k, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim4, n, k, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim5, m, k, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim3, Dim5, i, j, k, l, m, n, m, k, l, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(m, k, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim5, m, k, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim4, l, k, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim5, Dim4, i, j, k, l, m, n, l, k, n, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(l, k, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim4, l, k, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim5, l, k, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim4, Dim5, i, j, k, l, m, n, l, k, m, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(l, k, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim5, l, k, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim3, k, n, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim4, Dim3, i, j, k, l, m, n, k, n, m, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(k, n, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim3, k, n, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim3, k, m, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim5, Dim3, i, j, k, l, m, n, k, m, n, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(k, m, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim3, k, m, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim4, k, n, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim3, Dim4, i, j, k, l, m, n, k, n, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(k, n, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim4, k, n, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim5, k, m, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim3, Dim5, i, j, k, l, m, n, k, m, l, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(k, m, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim5, k, m, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim4, k, l, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim5, Dim4, i, j, k, l, m, n, k, l, n, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(k, l, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim4, k, l, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim5, k, l, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n, k, l, m, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, i, j>(TensorExpr(a, b));
}

//B(k, l, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim5, k, l, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

} //namespace FTensor

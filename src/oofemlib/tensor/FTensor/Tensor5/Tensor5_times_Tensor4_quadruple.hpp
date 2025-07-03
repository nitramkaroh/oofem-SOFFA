#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, )*B(n, o, p, q, ) quadruple contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q, int DimA, int DimX, int DimY, int DimZ, int DimW, char a, char x, char y, char z, char w>
class Tensor5_times_Tensor4_quadruple
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim8, n, o, p, q> iterB;

public:
  typename promote<T, U>::V operator()(const int N1) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        for(int zz = 0; zz < DimZ; ++zz)
          for(int ww = 0; ww < DimW; ++ww)
            {
              // Permutation is where the indices get checked.
              result += Permutation5<DimA, DimX, DimY, DimZ, DimW, a, x, y, z, w>().eval(
                iterA, N1, xx, yy, zz, ww) * Permutation4<DimX, DimY, DimZ, DimW, x, y, z, w>().eval(iterB, xx, yy, zz, ww);
            }
    return result;
  }

  Tensor5_times_Tensor4_quadruple(
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &iter_a,
    const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim8, n, o, p, q> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, )*B(l, k, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim0, l, k, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim1, Dim0, i, j, k, l, m, l, k, j, i, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(l, k, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim0, l, k, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim0, k, l, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim1, Dim0, i, j, k, l, m, k, l, j, i, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(k, l, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim0, k, l, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim0, l, j, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim2, Dim0, i, j, k, l, m, l, j, k, i, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(l, j, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim0, l, j, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim0, k, j, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim3, Dim0, i, j, k, l, m, k, j, l, i, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(k, j, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim0, k, j, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim0, j, l, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim2, Dim0, i, j, k, l, m, j, l, k, i, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(j, l, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim0, j, l, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim0, j, k, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim3, Dim0, i, j, k, l, m, j, k, l, i, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(j, k, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim0, j, k, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim1, l, k, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim0, Dim1, i, j, k, l, m, l, k, i, j, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(l, k, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim1, l, k, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim1, k, l, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim0, Dim1, i, j, k, l, m, k, l, i, j, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(k, l, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim1, k, l, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim2, l, j, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim0, Dim2, i, j, k, l, m, l, j, i, k, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(l, j, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim2, l, j, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim3, k, j, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim0, Dim3, i, j, k, l, m, k, j, i, l, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(k, j, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim3, k, j, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim2, j, l, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim0, Dim2, i, j, k, l, m, j, l, i, k, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(j, l, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim2, j, l, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim3, j, k, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim0, Dim3, i, j, k, l, m, j, k, i, l, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(j, k, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim3, j, k, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim1, l, i, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim2, Dim1, i, j, k, l, m, l, i, k, j, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(l, i, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim1, l, i, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim1, k, i, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim3, Dim1, i, j, k, l, m, k, i, l, j, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(k, i, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim1, k, i, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim2, l, i, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim1, Dim2, i, j, k, l, m, l, i, j, k, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(l, i, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim2, l, i, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim3, k, i, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim1, Dim3, i, j, k, l, m, k, i, j, l, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(k, i, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim3, k, i, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, i, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim2, j, i, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim3, Dim2, i, j, k, l, m, j, i, l, k, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(j, i, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim2, j, i, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, i, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim3, j, i, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim2, Dim3, i, j, k, l, m, j, i, k, l, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(j, i, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim3, j, i, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim1, i, l, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim2, Dim1, i, j, k, l, m, i, l, k, j, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(i, l, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim1, i, l, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim1, i, k, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim3, Dim1, i, j, k, l, m, i, k, l, j, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(i, k, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim1, i, k, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim2, i, l, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim1, Dim2, i, j, k, l, m, i, l, j, k, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(i, l, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim2, i, l, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim3, i, k, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim1, Dim3, i, j, k, l, m, i, k, j, l, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(i, k, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim3, i, k, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim2, i, j, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim3, Dim2, i, j, k, l, m, i, j, l, k, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(i, j, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim2, i, j, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim2, Dim3, i, j, k, l, m, i, j, k, l, Dim4, Dim0, Dim1, Dim2, Dim3, m, i, j, k, l>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim4, m>(TensorExpr(a, b));
}

//B(i, j, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim0, m, k, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim1, Dim0, i, j, k, l, m, m, k, j, i, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(m, k, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim0, m, k, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim0, k, m, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim1, Dim0, i, j, k, l, m, k, m, j, i, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(k, m, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim0, k, m, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim0, m, j, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim2, Dim0, i, j, k, l, m, m, j, k, i, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(m, j, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim0, m, j, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim0, k, j, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim4, Dim0, i, j, k, l, m, k, j, m, i, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(k, j, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim0, k, j, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim0, j, m, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim2, Dim0, i, j, k, l, m, j, m, k, i, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(j, m, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim0, j, m, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim0, j, k, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim4, Dim0, i, j, k, l, m, j, k, m, i, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(j, k, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim0, j, k, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim1, m, k, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim0, Dim1, i, j, k, l, m, m, k, i, j, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(m, k, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim1, m, k, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim1, k, m, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim0, Dim1, i, j, k, l, m, k, m, i, j, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(k, m, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim1, k, m, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim2, m, j, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim0, Dim2, i, j, k, l, m, m, j, i, k, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(m, j, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim2, m, j, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim4, k, j, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim0, Dim4, i, j, k, l, m, k, j, i, m, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(k, j, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim4, k, j, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim2, j, m, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim0, Dim2, i, j, k, l, m, j, m, i, k, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(j, m, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim2, j, m, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim4, j, k, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim0, Dim4, i, j, k, l, m, j, k, i, m, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(j, k, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim4, j, k, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim1, m, i, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim2, Dim1, i, j, k, l, m, m, i, k, j, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(m, i, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim1, m, i, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim1, k, i, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim4, Dim1, i, j, k, l, m, k, i, m, j, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(k, i, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim1, k, i, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim2, m, i, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim1, Dim2, i, j, k, l, m, m, i, j, k, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(m, i, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim2, m, i, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim4, k, i, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim1, Dim4, i, j, k, l, m, k, i, j, m, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(k, i, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim4, k, i, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim2, j, i, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim4, Dim2, i, j, k, l, m, j, i, m, k, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(j, i, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim2, j, i, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim4, j, i, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim2, Dim4, i, j, k, l, m, j, i, k, m, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(j, i, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim4, j, i, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim1, i, m, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim2, Dim1, i, j, k, l, m, i, m, k, j, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(i, m, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim1, i, m, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim1, i, k, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim4, Dim1, i, j, k, l, m, i, k, m, j, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(i, k, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim1, i, k, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim2, i, m, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim1, Dim2, i, j, k, l, m, i, m, j, k, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(i, m, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim2, i, m, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim4, i, k, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim1, Dim4, i, j, k, l, m, i, k, j, m, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(i, k, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim4, i, k, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim2, i, j, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim4, Dim2, i, j, k, l, m, i, j, m, k, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(i, j, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim2, i, j, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim4, i, j, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim2, Dim4, i, j, k, l, m, i, j, k, m, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim3, l>(TensorExpr(a, b));
}

//B(i, j, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim4, i, j, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim0, m, l, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim1, Dim0, i, j, k, l, m, m, l, j, i, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(m, l, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim0, m, l, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim0, l, m, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim1, Dim0, i, j, k, l, m, l, m, j, i, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(l, m, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim0, l, m, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim0, m, j, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim3, Dim0, i, j, k, l, m, m, j, l, i, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(m, j, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim0, m, j, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim0, l, j, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim4, Dim0, i, j, k, l, m, l, j, m, i, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(l, j, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim0, l, j, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim0, j, m, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim3, Dim0, i, j, k, l, m, j, m, l, i, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(j, m, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim0, j, m, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim0, j, l, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim4, Dim0, i, j, k, l, m, j, l, m, i, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(j, l, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim0, j, l, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim1, m, l, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim0, Dim1, i, j, k, l, m, m, l, i, j, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(m, l, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim1, m, l, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim1, l, m, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim0, Dim1, i, j, k, l, m, l, m, i, j, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(l, m, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim1, l, m, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim3, m, j, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim0, Dim3, i, j, k, l, m, m, j, i, l, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(m, j, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim3, m, j, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim4, l, j, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim0, Dim4, i, j, k, l, m, l, j, i, m, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(l, j, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim4, l, j, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim3, j, m, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim0, Dim3, i, j, k, l, m, j, m, i, l, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(j, m, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim3, j, m, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim4, j, l, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim0, Dim4, i, j, k, l, m, j, l, i, m, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(j, l, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim4, j, l, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim1, m, i, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim3, Dim1, i, j, k, l, m, m, i, l, j, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(m, i, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim1, m, i, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim1, l, i, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim4, Dim1, i, j, k, l, m, l, i, m, j, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(l, i, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim1, l, i, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim3, m, i, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim1, Dim3, i, j, k, l, m, m, i, j, l, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(m, i, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim3, m, i, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim4, l, i, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim1, Dim4, i, j, k, l, m, l, i, j, m, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(l, i, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim4, l, i, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim3, j, i, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim4, Dim3, i, j, k, l, m, j, i, m, l, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(j, i, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim3, j, i, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim4, j, i, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim3, Dim4, i, j, k, l, m, j, i, l, m, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(j, i, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim4, j, i, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim1, i, m, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim3, Dim1, i, j, k, l, m, i, m, l, j, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(i, m, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim1, i, m, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim1, i, l, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim4, Dim1, i, j, k, l, m, i, l, m, j, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(i, l, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim1, i, l, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim3, i, m, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim1, Dim3, i, j, k, l, m, i, m, j, l, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(i, m, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim3, i, m, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim4, i, l, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim1, Dim4, i, j, k, l, m, i, l, j, m, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(i, l, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim4, i, l, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim3, i, j, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim4, Dim3, i, j, k, l, m, i, j, m, l, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(i, j, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim3, i, j, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim4, i, j, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim3, Dim4, i, j, k, l, m, i, j, l, m, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim2, k>(TensorExpr(a, b));
}

//B(i, j, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim4, i, j, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim0, m, l, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim2, Dim0, i, j, k, l, m, m, l, k, i, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(m, l, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim0, m, l, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim0, l, m, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim2, Dim0, i, j, k, l, m, l, m, k, i, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(l, m, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim0, l, m, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim0, m, k, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim3, Dim0, i, j, k, l, m, m, k, l, i, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(m, k, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim0, m, k, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim0, l, k, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim4, Dim0, i, j, k, l, m, l, k, m, i, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(l, k, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim0, l, k, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim0, k, m, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim3, Dim0, i, j, k, l, m, k, m, l, i, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(k, m, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim0, k, m, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim0, k, l, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim4, Dim0, i, j, k, l, m, k, l, m, i, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(k, l, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim0, k, l, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim2, m, l, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim0, Dim2, i, j, k, l, m, m, l, i, k, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(m, l, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim2, m, l, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim2, l, m, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim0, Dim2, i, j, k, l, m, l, m, i, k, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(l, m, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim2, l, m, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim3, m, k, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim0, Dim3, i, j, k, l, m, m, k, i, l, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(m, k, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim3, m, k, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim4, l, k, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim0, Dim4, i, j, k, l, m, l, k, i, m, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(l, k, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim4, l, k, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim3, k, m, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim0, Dim3, i, j, k, l, m, k, m, i, l, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(k, m, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim3, k, m, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim4, k, l, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim0, Dim4, i, j, k, l, m, k, l, i, m, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(k, l, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim4, k, l, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim2, m, i, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim3, Dim2, i, j, k, l, m, m, i, l, k, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(m, i, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim2, m, i, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim2, l, i, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim4, Dim2, i, j, k, l, m, l, i, m, k, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(l, i, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim2, l, i, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim3, m, i, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim2, Dim3, i, j, k, l, m, m, i, k, l, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(m, i, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim3, m, i, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim4, l, i, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim2, Dim4, i, j, k, l, m, l, i, k, m, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(l, i, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim4, l, i, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim3, k, i, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim4, Dim3, i, j, k, l, m, k, i, m, l, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(k, i, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim3, k, i, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim4, k, i, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim3, Dim4, i, j, k, l, m, k, i, l, m, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(k, i, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim4, k, i, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim2, i, m, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim3, Dim2, i, j, k, l, m, i, m, l, k, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(i, m, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim2, i, m, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim2, i, l, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim4, Dim2, i, j, k, l, m, i, l, m, k, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(i, l, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim2, i, l, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim3, i, m, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim2, Dim3, i, j, k, l, m, i, m, k, l, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(i, m, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim3, i, m, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim4, i, l, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim2, Dim4, i, j, k, l, m, i, l, k, m, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(i, l, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim4, i, l, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim3, i, k, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim4, Dim3, i, j, k, l, m, i, k, m, l, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(i, k, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim3, i, k, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim4, i, k, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim3, Dim4, i, j, k, l, m, i, k, l, m, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim1, j>(TensorExpr(a, b));
}

//B(i, k, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim4, i, k, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim1, m, l, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim2, Dim1, i, j, k, l, m, m, l, k, j, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(m, l, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim1, m, l, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim1, l, m, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim2, Dim1, i, j, k, l, m, l, m, k, j, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(l, m, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim1, l, m, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim1, m, k, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim3, Dim1, i, j, k, l, m, m, k, l, j, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(m, k, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim1, m, k, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim1, l, k, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim4, Dim1, i, j, k, l, m, l, k, m, j, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(l, k, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim1, l, k, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim1, k, m, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim3, Dim1, i, j, k, l, m, k, m, l, j, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(k, m, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim1, k, m, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim1, k, l, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim4, Dim1, i, j, k, l, m, k, l, m, j, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(k, l, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim1, k, l, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim2, m, l, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim1, Dim2, i, j, k, l, m, m, l, j, k, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(m, l, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim2, m, l, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim2, l, m, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim1, Dim2, i, j, k, l, m, l, m, j, k, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(l, m, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim2, l, m, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim3, m, k, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim1, Dim3, i, j, k, l, m, m, k, j, l, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(m, k, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim3, m, k, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim4, l, k, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim1, Dim4, i, j, k, l, m, l, k, j, m, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(l, k, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim4, l, k, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim3, k, m, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim1, Dim3, i, j, k, l, m, k, m, j, l, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(k, m, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim3, k, m, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim4, k, l, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim1, Dim4, i, j, k, l, m, k, l, j, m, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(k, l, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim4, k, l, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim2, m, j, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim3, Dim2, i, j, k, l, m, m, j, l, k, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(m, j, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim2, m, j, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim2, l, j, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim4, Dim2, i, j, k, l, m, l, j, m, k, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(l, j, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim2, l, j, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim3, m, j, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim2, Dim3, i, j, k, l, m, m, j, k, l, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(m, j, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim3, m, j, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim4, l, j, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim2, Dim4, i, j, k, l, m, l, j, k, m, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(l, j, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim4, l, j, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim3, k, j, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim4, Dim3, i, j, k, l, m, k, j, m, l, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(k, j, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim3, k, j, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim4, k, j, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim3, Dim4, i, j, k, l, m, k, j, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(k, j, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim4, k, j, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim2, j, m, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim3, Dim2, i, j, k, l, m, j, m, l, k, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(j, m, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim2, j, m, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim2, j, l, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim4, Dim2, i, j, k, l, m, j, l, m, k, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(j, l, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim2, j, l, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim3, j, m, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim2, Dim3, i, j, k, l, m, j, m, k, l, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(j, m, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim3, j, m, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim4, j, l, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim2, Dim4, i, j, k, l, m, j, l, k, m, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(j, l, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim4, j, l, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim3, j, k, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim4, Dim3, i, j, k, l, m, j, k, m, l, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(j, k, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim3, j, k, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim4, j, k, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_quadruple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m, j, k, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor1_Expr<TensorExpr, typename promote<T,U>::V, Dim0, i>(TensorExpr(a, b));
}

//B(j, k, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim4, j, k, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

} //namespace FTensor

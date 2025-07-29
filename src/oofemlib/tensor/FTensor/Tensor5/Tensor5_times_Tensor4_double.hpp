#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, )*B(n, o, p, q, ) double contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q, int DimA, int DimB, int DimC, int DimX, int DimY, int DimZ, int DimW, char a, char b, char c, char x, char y, char z, char w>
class Tensor5_times_Tensor4_double
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim8, n, o, p, q> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        {
          // Permutation is where the indices get checked.
          result += Permutation5<DimA, DimB, DimC, DimX, DimY, a, b, c, x, y>().eval(
            iterA, N1, N2, N3, xx, yy) * Permutation4<DimX, DimY, DimZ, DimW, x, y, z, w>().eval(iterB, xx, yy, N4, N5);
        }
    return result;
  }

  Tensor5_times_Tensor4_double(
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &iter_a,
    const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim8, n, o, p, q> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, )*B(j, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim6, j, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, Dim6, i, j, k, l, m, j, i, n, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(j, i, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim6, j, i, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim6, i, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, i, j, k, l, m, i, j, n, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(i, j, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim6, i, j, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim6, j, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, Dim6, i, j, k, l, m, j, n, i, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(j, n, i, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim6, j, n, i, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim6, i, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, Dim6, i, j, k, l, m, i, n, j, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(i, n, j, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim6, i, n, j, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim0, j, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim0, i, j, k, l, m, j, n, o, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(j, n, o, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim0, j, n, o, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim1, i, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim1, i, j, k, l, m, i, n, o, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(i, n, o, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim1, i, n, o, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim6, n, j, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, i, j, k, l, m, n, j, i, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(n, j, i, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim6, n, j, i, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim6, n, i, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, i, j, k, l, m, n, i, j, o, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(n, i, j, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim6, n, i, j, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim0, n, j, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, i, j, k, l, m, n, j, o, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(n, j, o, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim0, n, j, o, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim1, n, i, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, i, j, k, l, m, n, i, o, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(n, i, o, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim1, n, i, o, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim0, n, o, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, i, j, k, l, m, n, o, j, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(n, o, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim0, n, o, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim1, n, o, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, i, j, k, l, m, n, o, i, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, k, l, m, i, j, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(n, o, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim1, n, o, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim6, k, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, Dim6, i, j, k, l, m, k, i, n, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(k, i, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim6, k, i, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim6, i, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, i, j, k, l, m, i, k, n, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(i, k, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim6, i, k, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim6, k, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, Dim6, i, j, k, l, m, k, n, i, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(k, n, i, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim6, k, n, i, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim6, i, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, Dim6, i, j, k, l, m, i, n, k, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(i, n, k, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim6, i, n, k, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim0, k, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim0, i, j, k, l, m, k, n, o, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(k, n, o, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim0, k, n, o, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim2, i, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim2, i, j, k, l, m, i, n, o, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(i, n, o, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim2, i, n, o, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim6, n, k, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, i, j, k, l, m, n, k, i, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(n, k, i, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim6, n, k, i, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim6, n, i, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, i, j, k, l, m, n, i, k, o, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(n, i, k, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim6, n, i, k, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim0, n, k, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, i, j, k, l, m, n, k, o, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(n, k, o, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim0, n, k, o, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim2, n, i, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, i, j, k, l, m, n, i, o, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(n, i, o, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim2, n, i, o, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim0, n, o, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, i, j, k, l, m, n, o, k, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(n, o, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim0, n, o, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim2, n, o, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, i, j, k, l, m, n, o, i, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, j, l, m, i, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(n, o, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim2, n, o, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim6, l, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, Dim6, i, j, k, l, m, l, i, n, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(l, i, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim6, l, i, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim6, i, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, Dim6, i, j, k, l, m, i, l, n, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(i, l, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim6, i, l, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim6, l, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, Dim6, i, j, k, l, m, l, n, i, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(l, n, i, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim6, l, n, i, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim6, i, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, Dim6, i, j, k, l, m, i, n, l, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(i, n, l, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim6, i, n, l, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim0, l, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim0, i, j, k, l, m, l, n, o, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(l, n, o, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim0, l, n, o, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim3, i, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim3, i, j, k, l, m, i, n, o, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(i, n, o, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim3, i, n, o, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim6, n, l, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim6, i, j, k, l, m, n, l, i, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(n, l, i, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim6, n, l, i, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim6, n, i, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim6, i, j, k, l, m, n, i, l, o, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(n, i, l, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim6, n, i, l, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim0, n, l, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim0, i, j, k, l, m, n, l, o, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(n, l, o, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim0, n, l, o, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim3, n, i, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim3, i, j, k, l, m, n, i, o, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(n, i, o, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim3, n, i, o, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim0, n, o, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, i, j, k, l, m, n, o, l, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(n, o, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim0, n, o, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim3, n, o, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, i, j, k, l, m, n, o, i, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, j, k, m, i, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(n, o, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim3, n, o, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim6, m, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim5, Dim6, i, j, k, l, m, m, i, n, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(m, i, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim6, m, i, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim6, i, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim5, Dim6, i, j, k, l, m, i, m, n, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(i, m, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim6, i, m, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim6, m, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim0, Dim6, i, j, k, l, m, m, n, i, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(m, n, i, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim6, m, n, i, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim6, i, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim4, Dim6, i, j, k, l, m, i, n, m, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(i, n, m, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim6, i, n, m, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim0, m, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim0, i, j, k, l, m, m, n, o, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(m, n, o, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim0, m, n, o, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim4, i, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim4, i, j, k, l, m, i, n, o, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(i, n, o, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim4, i, n, o, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim6, n, m, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim6, i, j, k, l, m, n, m, i, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(n, m, i, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim6, n, m, i, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim6, n, i, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim6, i, j, k, l, m, n, i, m, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(n, i, m, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim6, n, i, m, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim0, n, m, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim0, i, j, k, l, m, n, m, o, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(n, m, o, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim0, n, m, o, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim4, n, i, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim4, i, j, k, l, m, n, i, o, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(n, i, o, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim4, n, i, o, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim0, n, o, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, i, j, k, l, m, n, o, m, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(n, o, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim0, n, o, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim4, n, o, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, i, j, k, l, m, n, o, i, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(n, o, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim4, n, o, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim6, k, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, Dim6, i, j, k, l, m, k, j, n, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(k, j, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim6, k, j, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim6, j, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, j, k, l, m, j, k, n, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(j, k, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim6, j, k, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim6, k, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, Dim6, i, j, k, l, m, k, n, j, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(k, n, j, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim6, k, n, j, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim6, j, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, Dim6, i, j, k, l, m, j, n, k, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(j, n, k, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim6, j, n, k, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim1, k, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim1, i, j, k, l, m, k, n, o, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(k, n, o, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim1, k, n, o, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim2, j, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim2, i, j, k, l, m, j, n, o, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(j, n, o, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim2, j, n, o, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim6, n, k, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, i, j, k, l, m, n, k, j, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(n, k, j, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim6, n, k, j, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim6, n, j, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, i, j, k, l, m, n, j, k, o, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(n, j, k, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim6, n, j, k, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim1, n, k, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, i, j, k, l, m, n, k, o, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(n, k, o, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim1, n, k, o, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim2, n, j, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, i, j, k, l, m, n, j, o, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(n, j, o, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim2, n, j, o, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim1, n, o, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, i, j, k, l, m, n, o, k, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(n, o, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim1, n, o, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim2, n, o, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, i, j, k, l, m, n, o, j, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, i, l, m, j, k, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(n, o, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim2, n, o, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim6, l, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, Dim6, i, j, k, l, m, l, j, n, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(l, j, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim6, l, j, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim6, j, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, Dim6, i, j, k, l, m, j, l, n, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(j, l, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim6, j, l, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim6, l, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, Dim6, i, j, k, l, m, l, n, j, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(l, n, j, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim6, l, n, j, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim6, j, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, Dim6, i, j, k, l, m, j, n, l, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(j, n, l, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim6, j, n, l, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim1, l, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim1, i, j, k, l, m, l, n, o, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(l, n, o, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim1, l, n, o, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim3, j, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim3, i, j, k, l, m, j, n, o, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(j, n, o, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim3, j, n, o, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim6, n, l, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim6, i, j, k, l, m, n, l, j, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(n, l, j, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim6, n, l, j, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim6, n, j, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim6, i, j, k, l, m, n, j, l, o, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(n, j, l, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim6, n, j, l, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim1, n, l, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim1, i, j, k, l, m, n, l, o, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(n, l, o, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim1, n, l, o, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim3, n, j, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim3, i, j, k, l, m, n, j, o, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(n, j, o, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim3, n, j, o, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim1, n, o, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, i, j, k, l, m, n, o, l, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(n, o, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim1, n, o, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim3, n, o, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, i, j, k, l, m, n, o, j, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, i, k, m, j, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(n, o, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim3, n, o, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim6, m, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim5, Dim6, i, j, k, l, m, m, j, n, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(m, j, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim6, m, j, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim6, j, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim5, Dim6, i, j, k, l, m, j, m, n, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(j, m, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim6, j, m, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim6, m, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim1, Dim6, i, j, k, l, m, m, n, j, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(m, n, j, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim6, m, n, j, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim6, j, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim4, Dim6, i, j, k, l, m, j, n, m, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(j, n, m, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim6, j, n, m, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim1, m, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim1, i, j, k, l, m, m, n, o, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(m, n, o, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim1, m, n, o, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim4, j, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim4, i, j, k, l, m, j, n, o, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(j, n, o, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim4, j, n, o, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim6, n, m, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim6, i, j, k, l, m, n, m, j, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(n, m, j, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim6, n, m, j, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim6, n, j, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim6, i, j, k, l, m, n, j, m, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(n, j, m, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim6, n, j, m, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim1, n, m, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim1, i, j, k, l, m, n, m, o, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(n, m, o, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim1, n, m, o, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim4, n, j, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim4, i, j, k, l, m, n, j, o, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(n, j, o, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim4, n, j, o, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim1, n, o, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, i, j, k, l, m, n, o, m, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(n, o, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim1, n, o, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim4, n, o, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, i, j, k, l, m, n, o, j, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(n, o, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim4, n, o, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim6, l, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, Dim6, i, j, k, l, m, l, k, n, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(l, k, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim6, l, k, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim6, k, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, k, l, m, k, l, n, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(k, l, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim6, k, l, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim6, l, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, Dim6, i, j, k, l, m, l, n, k, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(l, n, k, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim6, l, n, k, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim6, k, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, Dim6, i, j, k, l, m, k, n, l, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(k, n, l, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim6, k, n, l, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim2, l, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim2, i, j, k, l, m, l, n, o, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(l, n, o, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim2, l, n, o, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim3, k, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim3, i, j, k, l, m, k, n, o, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(k, n, o, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim3, k, n, o, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim6, n, l, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim6, i, j, k, l, m, n, l, k, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(n, l, k, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim6, n, l, k, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim6, n, k, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim6, i, j, k, l, m, n, k, l, o, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(n, k, l, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim6, n, k, l, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim2, n, l, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim2, i, j, k, l, m, n, l, o, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(n, l, o, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim2, n, l, o, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim3, n, k, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim3, i, j, k, l, m, n, k, o, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(n, k, o, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim3, n, k, o, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim2, n, o, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, i, j, k, l, m, n, o, l, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(n, o, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim2, n, o, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim3, n, o, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, i, j, k, l, m, n, o, k, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, i, j, m, k, l, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(n, o, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim3, n, o, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim6, m, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim5, Dim6, i, j, k, l, m, m, k, n, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(m, k, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim6, m, k, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim6, k, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim5, Dim6, i, j, k, l, m, k, m, n, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(k, m, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim6, k, m, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim6, m, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim2, Dim6, i, j, k, l, m, m, n, k, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(m, n, k, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim6, m, n, k, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim6, k, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim4, Dim6, i, j, k, l, m, k, n, m, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(k, n, m, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim6, k, n, m, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim2, m, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim2, i, j, k, l, m, m, n, o, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(m, n, o, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim2, m, n, o, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim4, k, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim4, i, j, k, l, m, k, n, o, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(k, n, o, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim4, k, n, o, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim6, n, m, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim6, i, j, k, l, m, n, m, k, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(n, m, k, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim6, n, m, k, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim6, n, k, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim6, i, j, k, l, m, n, k, m, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(n, k, m, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim6, n, k, m, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim2, n, m, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim2, i, j, k, l, m, n, m, o, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(n, m, o, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim2, n, m, o, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim4, n, k, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim4, i, j, k, l, m, n, k, o, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(n, k, o, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim4, n, k, o, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim2, n, o, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, i, j, k, l, m, n, o, m, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(n, o, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim2, n, o, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim4, n, o, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, i, j, k, l, m, n, o, k, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(n, o, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim4, n, o, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim6, m, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim5, Dim6, i, j, k, l, m, m, l, n, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(m, l, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim6, m, l, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim6, l, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, l, m, n, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(l, m, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim6, l, m, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim6, m, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim3, Dim6, i, j, k, l, m, m, n, l, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(m, n, l, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim6, m, n, l, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim6, l, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim4, Dim6, i, j, k, l, m, l, n, m, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(l, n, m, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim6, l, n, m, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim3, m, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim3, i, j, k, l, m, m, n, o, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(m, n, o, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim3, m, n, o, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim4, l, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim4, i, j, k, l, m, l, n, o, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(l, n, o, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim4, l, n, o, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim6, n, m, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim6, i, j, k, l, m, n, m, l, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(n, m, l, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim6, n, m, l, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim6, n, l, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim6, i, j, k, l, m, n, l, m, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(n, l, m, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim6, n, l, m, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim3, n, m, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim3, i, j, k, l, m, n, m, o, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(n, m, o, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim3, n, m, o, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim4, n, l, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim4, i, j, k, l, m, n, l, o, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(n, l, o, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim4, n, l, o, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim3, n, o, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, i, j, k, l, m, n, o, m, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(n, o, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim3, n, o, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim4, n, o, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, i, j, k, l, m, n, o, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(n, o, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim4, n, o, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

} //namespace FTensor

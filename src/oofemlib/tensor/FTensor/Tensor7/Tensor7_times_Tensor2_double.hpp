#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, n, o, )*B(p, q, ) double contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q, int DimA, int DimB, int DimC, int DimD, int DimE, int DimX, int DimY, char a, char b, char c, char d, char e, char x, char y>
class Tensor7_times_Tensor2_double
{
  Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
  Tensor2_Expr<B, U, Dim7, Dim8, p, q> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        {
          // Permutation is where the indices get checked.
          result += Permutation7<DimA, DimB, DimC, DimD, DimE, DimX, DimY, a, b, c, d, e, x, y>().eval(
            iterA, N1, N2, N3, N4, N5, xx, yy) * Permutation2<DimX, DimY, x, y>().eval(iterB, xx, yy);
        }
    return result;
  }

  Tensor7_times_Tensor2_double(
    const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &iter_a,
    const Tensor2_Expr<B, U, Dim7, Dim8, p, q> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, n, o, )*B(j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, i, j, k, l, m, n, o, j, i, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, k, l, m, n, o, i, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(j, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, i, j, k, l, m, n, o, i, j, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, k, l, m, n, o, i, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, k, l, m, n, o>(TensorExpr(a, b));
}

//B(i, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim2, Dim0, k, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, i, j, k, l, m, n, o, k, i, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, j, l, m, n, o, i, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(k, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim2, Dim0, k, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim0, Dim2, i, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, i, j, k, l, m, n, o, i, k, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, j, l, m, n, o, i, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, j, l, m, n, o>(TensorExpr(a, b));
}

//B(i, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim0, Dim2, i, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim3, Dim0, l, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, i, j, k, l, m, n, o, l, i, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, j, k, m, n, o, i, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(l, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim3, Dim0, l, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim0, Dim3, i, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, i, j, k, l, m, n, o, i, l, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, j, k, m, n, o, i, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, j, k, m, n, o>(TensorExpr(a, b));
}

//B(i, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim0, Dim3, i, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim4, Dim0, m, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, i, j, k, l, m, n, o, m, i, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, j, k, l, n, o, i, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(m, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim4, Dim0, m, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim0, Dim4, i, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, i, j, k, l, m, n, o, i, m, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, j, k, l, n, o, i, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, j, k, l, n, o>(TensorExpr(a, b));
}

//B(i, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim0, Dim4, i, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim5, Dim0, n, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, i, j, k, l, m, n, o, n, i, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, j, k, l, m, o, i, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, j, k, l, m, o>(TensorExpr(a, b));
}

//B(n, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim5, Dim0, n, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim0, Dim5, i, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, i, j, k, l, m, n, o, i, n, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, j, k, l, m, o, i, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, j, k, l, m, o>(TensorExpr(a, b));
}

//B(i, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim0, Dim5, i, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim6, Dim0, o, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim0, i, j, k, l, m, n, o, o, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, j, k, l, m, n, i, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, j, k, l, m, n>(TensorExpr(a, b));
}

//B(o, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim6, Dim0, o, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim0, Dim6, i, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim6, i, j, k, l, m, n, o, i, o, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, j, k, l, m, n, i, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, j, k, l, m, n>(TensorExpr(a, b));
}

//B(i, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim0, Dim6, i, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, i, j, k, l, m, n, o, k, j, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, i, l, m, n, o, j, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(k, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, i, j, k, l, m, n, o, j, k, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, i, l, m, n, o, j, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, i, l, m, n, o>(TensorExpr(a, b));
}

//B(j, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim3, Dim1, l, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, i, j, k, l, m, n, o, l, j, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, i, k, m, n, o, j, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(l, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim3, Dim1, l, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim1, Dim3, j, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, i, j, k, l, m, n, o, j, l, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, i, k, m, n, o, j, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, i, k, m, n, o>(TensorExpr(a, b));
}

//B(j, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim1, Dim3, j, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim4, Dim1, m, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, i, j, k, l, m, n, o, m, j, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, i, k, l, n, o, j, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(m, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim4, Dim1, m, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim1, Dim4, j, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, i, j, k, l, m, n, o, j, m, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, i, k, l, n, o, j, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, i, k, l, n, o>(TensorExpr(a, b));
}

//B(j, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim1, Dim4, j, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim5, Dim1, n, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, i, j, k, l, m, n, o, n, j, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, i, k, l, m, o, j, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, i, k, l, m, o>(TensorExpr(a, b));
}

//B(n, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim5, Dim1, n, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim1, Dim5, j, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, i, j, k, l, m, n, o, j, n, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, i, k, l, m, o, j, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, i, k, l, m, o>(TensorExpr(a, b));
}

//B(j, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim1, Dim5, j, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim6, Dim1, o, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim1, i, j, k, l, m, n, o, o, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, i, k, l, m, n, j, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, i, k, l, m, n>(TensorExpr(a, b));
}

//B(o, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim6, Dim1, o, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim1, Dim6, j, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim6, i, j, k, l, m, n, o, j, o, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, i, k, l, m, n, j, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, i, k, l, m, n>(TensorExpr(a, b));
}

//B(j, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim1, Dim6, j, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, i, j, k, l, m, n, o, l, k, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, i, j, m, n, o, k, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(l, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, i, j, k, l, m, n, o, k, l, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, i, j, m, n, o, k, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, i, j, m, n, o>(TensorExpr(a, b));
}

//B(k, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim4, Dim2, m, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, i, j, k, l, m, n, o, m, k, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, i, j, l, n, o, k, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(m, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim4, Dim2, m, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim2, Dim4, k, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, i, j, k, l, m, n, o, k, m, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, i, j, l, n, o, k, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, i, j, l, n, o>(TensorExpr(a, b));
}

//B(k, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim2, Dim4, k, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim5, Dim2, n, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, i, j, k, l, m, n, o, n, k, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, i, j, l, m, o, k, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, i, j, l, m, o>(TensorExpr(a, b));
}

//B(n, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim5, Dim2, n, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim2, Dim5, k, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, i, j, k, l, m, n, o, k, n, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, i, j, l, m, o, k, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, i, j, l, m, o>(TensorExpr(a, b));
}

//B(k, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim2, Dim5, k, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim6, Dim2, o, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim2, i, j, k, l, m, n, o, o, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, i, j, l, m, n, k, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, i, j, l, m, n>(TensorExpr(a, b));
}

//B(o, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim6, Dim2, o, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim2, Dim6, k, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim6, i, j, k, l, m, n, o, k, o, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, i, j, l, m, n, k, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, i, j, l, m, n>(TensorExpr(a, b));
}

//B(k, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim2, Dim6, k, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim4, Dim3, m, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, i, j, k, l, m, n, o, m, l, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, i, j, k, n, o, l, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(m, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim4, Dim3, m, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, i, j, k, l, m, n, o, l, m, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, i, j, k, n, o, l, m>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, i, j, k, n, o>(TensorExpr(a, b));
}

//B(l, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim5, Dim3, n, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, i, j, k, l, m, n, o, n, l, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, i, j, k, m, o, l, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, i, j, k, m, o>(TensorExpr(a, b));
}

//B(n, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim5, Dim3, n, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim3, Dim5, l, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, i, j, k, l, m, n, o, l, n, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, i, j, k, m, o, l, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, i, j, k, m, o>(TensorExpr(a, b));
}

//B(l, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim3, Dim5, l, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim6, Dim3, o, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim3, i, j, k, l, m, n, o, o, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, i, j, k, m, n, l, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, i, j, k, m, n>(TensorExpr(a, b));
}

//B(o, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim6, Dim3, o, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim3, Dim6, l, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim6, i, j, k, l, m, n, o, l, o, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, i, j, k, m, n, l, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, i, j, k, m, n>(TensorExpr(a, b));
}

//B(l, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim3, Dim6, l, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim5, Dim4, n, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, i, j, k, l, m, n, o, n, m, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, i, j, k, l, o, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, i, j, k, l, o>(TensorExpr(a, b));
}

//B(n, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim5, Dim4, n, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim4, Dim5, m, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, i, j, k, l, m, n, o, m, n, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, i, j, k, l, o, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, i, j, k, l, o>(TensorExpr(a, b));
}

//B(m, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim4, Dim5, m, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim6, Dim4, o, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim4, i, j, k, l, m, n, o, o, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, i, j, k, l, n, m, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, i, j, k, l, n>(TensorExpr(a, b));
}

//B(o, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim6, Dim4, o, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim4, Dim6, m, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim6, i, j, k, l, m, n, o, m, o, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, i, j, k, l, n, m, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, i, j, k, l, n>(TensorExpr(a, b));
}

//B(m, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim4, Dim6, m, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim6, Dim5, o, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim5, i, j, k, l, m, n, o, o, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>(TensorExpr(a, b));
}

//B(o, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim6, Dim5, o, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim5, Dim6, n, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim6, i, j, k, l, m, n, o, n, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>(TensorExpr(a, b));
}

//B(n, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor2_Expr<B, U, Dim5, Dim6, n, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

} //namespace FTensor

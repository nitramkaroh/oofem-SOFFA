#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, n, o, )*B(p, q, r, s, ) double contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, int Dim9, int Dim10, char i, char j, char k, char l, char m, char n, char o, char p, char q, char r, char s, int DimA, int DimB, int DimC, int DimD, int DimE, int DimX, int DimY, int DimZ, int DimW, char a, char b, char c, char d, char e, char x, char y, char z, char w>
class Tensor7_times_Tensor4_double
{
  Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
  Tensor4_Expr<B, U, Dim7, Dim8, Dim9, Dim10, p, q, r, s> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const int N7) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        {
          // Permutation is where the indices get checked.
          result += Permutation7<DimA, DimB, DimC, DimD, DimE, DimX, DimY, a, b, c, d, e, x, y>().eval(
            iterA, N1, N2, N3, N4, N5, xx, yy) * Permutation4<DimX, DimY, DimZ, DimW, x, y, z, w>().eval(iterB, xx, yy, N6, N7);
        }
    return result;
  }

  Tensor7_times_Tensor4_double(
    const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &iter_a,
    const Tensor4_Expr<B, U, Dim7, Dim8, Dim9, Dim10, p, q, r, s> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, n, o, )*B(j, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim7, Dim8, j, i, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim7, Dim8, i, j, k, l, m, n, o, j, i, p, q, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, k, l, m, n, o, i, j, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, i, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim7, Dim8, j, i, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim7, Dim8, i, j, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, i, j, k, l, m, n, o, i, j, p, q, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, k, l, m, n, o, i, j, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, j, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim7, Dim8, i, j, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim7, Dim0, Dim8, j, p, i, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim0, Dim8, i, j, k, l, m, n, o, j, p, i, q, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, k, l, m, n, o, i, j, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, p, i, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim7, Dim0, Dim8, j, p, i, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim7, Dim1, Dim8, i, p, j, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim1, Dim8, i, j, k, l, m, n, o, i, p, j, q, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, k, l, m, n, o, i, j, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, p, j, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim7, Dim1, Dim8, i, p, j, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim7, Dim8, Dim0, j, p, q, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim8, Dim0, i, j, k, l, m, n, o, j, p, q, i, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, k, l, m, n, o, i, j, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, p, q, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim7, Dim8, Dim0, j, p, q, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim7, Dim8, Dim1, i, p, q, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim8, Dim1, i, j, k, l, m, n, o, i, p, q, j, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, k, l, m, n, o, i, j, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, p, q, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim7, Dim8, Dim1, i, p, q, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim1, Dim0, Dim8, p, j, i, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim0, Dim8, i, j, k, l, m, n, o, p, j, i, q, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, k, l, m, n, o, i, j, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, j, i, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim1, Dim0, Dim8, p, j, i, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim0, Dim1, Dim8, p, i, j, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim1, Dim8, i, j, k, l, m, n, o, p, i, j, q, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, k, l, m, n, o, i, j, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, i, j, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim0, Dim1, Dim8, p, i, j, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim1, Dim8, Dim0, p, j, q, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim8, Dim0, i, j, k, l, m, n, o, p, j, q, i, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, k, l, m, n, o, i, j, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, j, q, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim1, Dim8, Dim0, p, j, q, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim0, Dim8, Dim1, p, i, q, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim8, Dim1, i, j, k, l, m, n, o, p, i, q, j, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, k, l, m, n, o, i, j, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, i, q, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim0, Dim8, Dim1, p, i, q, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim1, Dim0, p, q, j, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim1, Dim0, i, j, k, l, m, n, o, p, q, j, i, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, k, l, m, n, o, i, j, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, j, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim1, Dim0, p, q, j, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim0, Dim1, p, q, i, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim0, Dim1, i, j, k, l, m, n, o, p, q, i, j, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, k, l, m, n, o, i, j, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, i, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim0, Dim1, p, q, i, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim7, Dim8, k, i, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim7, Dim8, i, j, k, l, m, n, o, k, i, p, q, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, j, l, m, n, o, i, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, i, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim7, Dim8, k, i, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim7, Dim8, i, k, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, i, j, k, l, m, n, o, i, k, p, q, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, j, l, m, n, o, i, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, k, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim7, Dim8, i, k, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim7, Dim0, Dim8, k, p, i, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim0, Dim8, i, j, k, l, m, n, o, k, p, i, q, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, j, l, m, n, o, i, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, p, i, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim7, Dim0, Dim8, k, p, i, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim7, Dim2, Dim8, i, p, k, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim2, Dim8, i, j, k, l, m, n, o, i, p, k, q, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, j, l, m, n, o, i, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, p, k, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim7, Dim2, Dim8, i, p, k, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim7, Dim8, Dim0, k, p, q, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim8, Dim0, i, j, k, l, m, n, o, k, p, q, i, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, j, l, m, n, o, i, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, p, q, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim7, Dim8, Dim0, k, p, q, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim7, Dim8, Dim2, i, p, q, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim8, Dim2, i, j, k, l, m, n, o, i, p, q, k, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, j, l, m, n, o, i, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, p, q, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim7, Dim8, Dim2, i, p, q, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim2, Dim0, Dim8, p, k, i, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim0, Dim8, i, j, k, l, m, n, o, p, k, i, q, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, j, l, m, n, o, i, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, k, i, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim2, Dim0, Dim8, p, k, i, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim0, Dim2, Dim8, p, i, k, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim2, Dim8, i, j, k, l, m, n, o, p, i, k, q, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, j, l, m, n, o, i, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, i, k, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim0, Dim2, Dim8, p, i, k, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim2, Dim8, Dim0, p, k, q, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim8, Dim0, i, j, k, l, m, n, o, p, k, q, i, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, j, l, m, n, o, i, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, k, q, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim2, Dim8, Dim0, p, k, q, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim0, Dim8, Dim2, p, i, q, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim8, Dim2, i, j, k, l, m, n, o, p, i, q, k, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, j, l, m, n, o, i, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, i, q, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim0, Dim8, Dim2, p, i, q, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim2, Dim0, p, q, k, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim2, Dim0, i, j, k, l, m, n, o, p, q, k, i, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, j, l, m, n, o, i, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, k, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim2, Dim0, p, q, k, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim0, Dim2, p, q, i, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim0, Dim2, i, j, k, l, m, n, o, p, q, i, k, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, j, l, m, n, o, i, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, i, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim0, Dim2, p, q, i, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim7, Dim8, l, i, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim7, Dim8, i, j, k, l, m, n, o, l, i, p, q, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, j, k, m, n, o, i, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, i, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim7, Dim8, l, i, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim7, Dim8, i, l, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, i, j, k, l, m, n, o, i, l, p, q, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, j, k, m, n, o, i, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, l, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim7, Dim8, i, l, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim7, Dim0, Dim8, l, p, i, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim0, Dim8, i, j, k, l, m, n, o, l, p, i, q, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, j, k, m, n, o, i, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, p, i, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim7, Dim0, Dim8, l, p, i, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim7, Dim3, Dim8, i, p, l, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim3, Dim8, i, j, k, l, m, n, o, i, p, l, q, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, j, k, m, n, o, i, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, p, l, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim7, Dim3, Dim8, i, p, l, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim7, Dim8, Dim0, l, p, q, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim8, Dim0, i, j, k, l, m, n, o, l, p, q, i, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, j, k, m, n, o, i, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, p, q, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim7, Dim8, Dim0, l, p, q, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim7, Dim8, Dim3, i, p, q, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim8, Dim3, i, j, k, l, m, n, o, i, p, q, l, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, j, k, m, n, o, i, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, p, q, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim7, Dim8, Dim3, i, p, q, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim3, Dim0, Dim8, p, l, i, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim0, Dim8, i, j, k, l, m, n, o, p, l, i, q, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, j, k, m, n, o, i, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, l, i, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim3, Dim0, Dim8, p, l, i, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim0, Dim3, Dim8, p, i, l, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim3, Dim8, i, j, k, l, m, n, o, p, i, l, q, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, j, k, m, n, o, i, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, i, l, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim0, Dim3, Dim8, p, i, l, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim3, Dim8, Dim0, p, l, q, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim8, Dim0, i, j, k, l, m, n, o, p, l, q, i, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, j, k, m, n, o, i, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, l, q, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim3, Dim8, Dim0, p, l, q, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim0, Dim8, Dim3, p, i, q, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim8, Dim3, i, j, k, l, m, n, o, p, i, q, l, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, j, k, m, n, o, i, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, i, q, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim0, Dim8, Dim3, p, i, q, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim3, Dim0, p, q, l, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim3, Dim0, i, j, k, l, m, n, o, p, q, l, i, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, j, k, m, n, o, i, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, l, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim3, Dim0, p, q, l, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim0, Dim3, p, q, i, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim0, Dim3, i, j, k, l, m, n, o, p, q, i, l, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, j, k, m, n, o, i, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, i, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim0, Dim3, p, q, i, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim7, Dim8, m, i, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim7, Dim8, i, j, k, l, m, n, o, m, i, p, q, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, j, k, l, n, o, i, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, i, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim7, Dim8, m, i, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim7, Dim8, i, m, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, i, j, k, l, m, n, o, i, m, p, q, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, j, k, l, n, o, i, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(i, m, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim7, Dim8, i, m, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim7, Dim0, Dim8, m, p, i, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim0, Dim8, i, j, k, l, m, n, o, m, p, i, q, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, j, k, l, n, o, i, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, p, i, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim7, Dim0, Dim8, m, p, i, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim7, Dim4, Dim8, i, p, m, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim4, Dim8, i, j, k, l, m, n, o, i, p, m, q, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, j, k, l, n, o, i, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(i, p, m, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim7, Dim4, Dim8, i, p, m, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim7, Dim8, Dim0, m, p, q, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim8, Dim0, i, j, k, l, m, n, o, m, p, q, i, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, j, k, l, n, o, i, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, p, q, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim7, Dim8, Dim0, m, p, q, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim7, Dim8, Dim4, i, p, q, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim8, Dim4, i, j, k, l, m, n, o, i, p, q, m, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, j, k, l, n, o, i, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(i, p, q, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim7, Dim8, Dim4, i, p, q, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim4, Dim0, Dim8, p, m, i, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim0, Dim8, i, j, k, l, m, n, o, p, m, i, q, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, j, k, l, n, o, i, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, m, i, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim4, Dim0, Dim8, p, m, i, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim0, Dim4, Dim8, p, i, m, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim4, Dim8, i, j, k, l, m, n, o, p, i, m, q, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, j, k, l, n, o, i, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, i, m, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim0, Dim4, Dim8, p, i, m, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim4, Dim8, Dim0, p, m, q, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim8, Dim0, i, j, k, l, m, n, o, p, m, q, i, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, j, k, l, n, o, i, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, m, q, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim4, Dim8, Dim0, p, m, q, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim0, Dim8, Dim4, p, i, q, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim8, Dim4, i, j, k, l, m, n, o, p, i, q, m, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, j, k, l, n, o, i, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, i, q, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim0, Dim8, Dim4, p, i, q, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim4, Dim0, p, q, m, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim4, Dim0, i, j, k, l, m, n, o, p, q, m, i, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, j, k, l, n, o, i, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, m, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim4, Dim0, p, q, m, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim0, Dim4, p, q, i, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim0, Dim4, i, j, k, l, m, n, o, p, q, i, m, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, j, k, l, n, o, i, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, i, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim0, Dim4, p, q, i, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim7, Dim8, n, i, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim7, Dim8, i, j, k, l, m, n, o, n, i, p, q, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, Dim8, j, k, l, m, o, i, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, i, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim7, Dim8, n, i, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim7, Dim8, i, n, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim7, Dim8, i, j, k, l, m, n, o, i, n, p, q, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, Dim8, j, k, l, m, o, i, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(i, n, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim7, Dim8, i, n, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim7, Dim0, Dim8, n, p, i, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim0, Dim8, i, j, k, l, m, n, o, n, p, i, q, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, Dim8, j, k, l, m, o, i, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, p, i, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim7, Dim0, Dim8, n, p, i, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim7, Dim5, Dim8, i, p, n, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim5, Dim8, i, j, k, l, m, n, o, i, p, n, q, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, Dim8, j, k, l, m, o, i, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(i, p, n, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim7, Dim5, Dim8, i, p, n, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim7, Dim8, Dim0, n, p, q, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim8, Dim0, i, j, k, l, m, n, o, n, p, q, i, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, Dim8, j, k, l, m, o, i, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, p, q, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim7, Dim8, Dim0, n, p, q, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim7, Dim8, Dim5, i, p, q, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim8, Dim5, i, j, k, l, m, n, o, i, p, q, n, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, Dim8, j, k, l, m, o, i, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(i, p, q, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim7, Dim8, Dim5, i, p, q, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim5, Dim0, Dim8, p, n, i, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim0, Dim8, i, j, k, l, m, n, o, p, n, i, q, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, Dim8, j, k, l, m, o, i, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, n, i, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim5, Dim0, Dim8, p, n, i, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim0, Dim5, Dim8, p, i, n, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim5, Dim8, i, j, k, l, m, n, o, p, i, n, q, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, Dim8, j, k, l, m, o, i, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, i, n, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim0, Dim5, Dim8, p, i, n, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim5, Dim8, Dim0, p, n, q, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim8, Dim0, i, j, k, l, m, n, o, p, n, q, i, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, Dim8, j, k, l, m, o, i, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, n, q, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim5, Dim8, Dim0, p, n, q, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim0, Dim8, Dim5, p, i, q, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim8, Dim5, i, j, k, l, m, n, o, p, i, q, n, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, Dim8, j, k, l, m, o, i, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, i, q, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim0, Dim8, Dim5, p, i, q, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim5, Dim0, p, q, n, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim5, Dim0, i, j, k, l, m, n, o, p, q, n, i, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, Dim8, j, k, l, m, o, i, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, q, n, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim5, Dim0, p, q, n, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim0, Dim5, p, q, i, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim0, Dim5, i, j, k, l, m, n, o, p, q, i, n, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, Dim8, j, k, l, m, o, i, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, q, i, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim0, Dim5, p, q, i, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim7, Dim8, o, i, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim0, Dim7, Dim8, i, j, k, l, m, n, o, o, i, p, q, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, j, k, l, m, n, i, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, j, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(o, i, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim7, Dim8, o, i, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim7, Dim8, i, o, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, i, o, p, q, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, j, k, l, m, n, i, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, j, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(i, o, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim7, Dim8, i, o, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim0, Dim8, o, p, i, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim0, Dim8, i, j, k, l, m, n, o, o, p, i, q, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, j, k, l, m, n, i, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, j, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(o, p, i, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim0, Dim8, o, p, i, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, o, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim7, Dim6, Dim8, i, p, o, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim6, Dim8, i, j, k, l, m, n, o, i, p, o, q, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, j, k, l, m, n, i, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, j, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(i, p, o, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim7, Dim6, Dim8, i, p, o, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim0, o, p, q, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim8, Dim0, i, j, k, l, m, n, o, o, p, q, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, j, k, l, m, n, i, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, j, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(o, p, q, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim0, o, p, q, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, q, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim0, Dim7, Dim8, Dim6, i, p, q, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim8, Dim6, i, j, k, l, m, n, o, i, p, q, o, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, j, k, l, m, n, i, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, j, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(i, p, q, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim7, Dim8, Dim6, i, p, q, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim6, Dim0, Dim8, p, o, i, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim0, Dim8, i, j, k, l, m, n, o, p, o, i, q, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, j, k, l, m, n, i, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, j, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, o, i, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim6, Dim0, Dim8, p, o, i, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, o, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim0, Dim6, Dim8, p, i, o, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim6, Dim8, i, j, k, l, m, n, o, p, i, o, q, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, j, k, l, m, n, i, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, j, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, i, o, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim0, Dim6, Dim8, p, i, o, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim6, Dim8, Dim0, p, o, q, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim8, Dim0, i, j, k, l, m, n, o, p, o, q, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, j, k, l, m, n, i, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, j, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, o, q, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim6, Dim8, Dim0, p, o, q, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, q, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim0, Dim8, Dim6, p, i, q, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim8, Dim6, i, j, k, l, m, n, o, p, i, q, o, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, j, k, l, m, n, i, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, j, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, i, q, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim0, Dim8, Dim6, p, i, q, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim6, Dim0, p, q, o, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim6, Dim0, i, j, k, l, m, n, o, p, q, o, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, j, k, l, m, n, i, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, j, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, q, o, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim6, Dim0, p, q, o, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim0, Dim6, p, q, i, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim0, Dim6, i, j, k, l, m, n, o, p, q, i, o, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, j, k, l, m, n, i, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, j, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, q, i, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim0, Dim6, p, q, i, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim7, Dim8, k, j, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim7, Dim8, i, j, k, l, m, n, o, k, j, p, q, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, l, m, n, o, j, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, j, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim7, Dim8, k, j, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim7, Dim8, j, k, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, j, k, l, m, n, o, j, k, p, q, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, l, m, n, o, j, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, k, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim7, Dim8, j, k, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim7, Dim1, Dim8, k, p, j, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim1, Dim8, i, j, k, l, m, n, o, k, p, j, q, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, l, m, n, o, j, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, p, j, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim7, Dim1, Dim8, k, p, j, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim7, Dim2, Dim8, j, p, k, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim2, Dim8, i, j, k, l, m, n, o, j, p, k, q, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, l, m, n, o, j, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, p, k, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim7, Dim2, Dim8, j, p, k, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim7, Dim8, Dim1, k, p, q, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim8, Dim1, i, j, k, l, m, n, o, k, p, q, j, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, l, m, n, o, j, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, p, q, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim7, Dim8, Dim1, k, p, q, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim7, Dim8, Dim2, j, p, q, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim8, Dim2, i, j, k, l, m, n, o, j, p, q, k, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, l, m, n, o, j, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, p, q, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim7, Dim8, Dim2, j, p, q, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim2, Dim1, Dim8, p, k, j, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim1, Dim8, i, j, k, l, m, n, o, p, k, j, q, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, l, m, n, o, j, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, k, j, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim2, Dim1, Dim8, p, k, j, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim1, Dim2, Dim8, p, j, k, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim2, Dim8, i, j, k, l, m, n, o, p, j, k, q, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, l, m, n, o, j, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, j, k, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim1, Dim2, Dim8, p, j, k, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim2, Dim8, Dim1, p, k, q, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim8, Dim1, i, j, k, l, m, n, o, p, k, q, j, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, l, m, n, o, j, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, k, q, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim2, Dim8, Dim1, p, k, q, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim1, Dim8, Dim2, p, j, q, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim8, Dim2, i, j, k, l, m, n, o, p, j, q, k, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, l, m, n, o, j, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, j, q, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim1, Dim8, Dim2, p, j, q, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim2, Dim1, p, q, k, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim2, Dim1, i, j, k, l, m, n, o, p, q, k, j, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, l, m, n, o, j, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, k, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim2, Dim1, p, q, k, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim1, Dim2, p, q, j, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim1, Dim2, i, j, k, l, m, n, o, p, q, j, k, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, l, m, n, o, j, k, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, j, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim1, Dim2, p, q, j, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim7, Dim8, l, j, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim7, Dim8, i, j, k, l, m, n, o, l, j, p, q, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, k, m, n, o, j, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, j, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim7, Dim8, l, j, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim7, Dim8, j, l, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, j, k, l, m, n, o, j, l, p, q, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, k, m, n, o, j, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, l, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim7, Dim8, j, l, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim7, Dim1, Dim8, l, p, j, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim1, Dim8, i, j, k, l, m, n, o, l, p, j, q, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, k, m, n, o, j, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, p, j, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim7, Dim1, Dim8, l, p, j, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim7, Dim3, Dim8, j, p, l, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim3, Dim8, i, j, k, l, m, n, o, j, p, l, q, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, k, m, n, o, j, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, p, l, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim7, Dim3, Dim8, j, p, l, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim7, Dim8, Dim1, l, p, q, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim8, Dim1, i, j, k, l, m, n, o, l, p, q, j, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, k, m, n, o, j, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, p, q, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim7, Dim8, Dim1, l, p, q, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim7, Dim8, Dim3, j, p, q, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim8, Dim3, i, j, k, l, m, n, o, j, p, q, l, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, k, m, n, o, j, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, p, q, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim7, Dim8, Dim3, j, p, q, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim3, Dim1, Dim8, p, l, j, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim1, Dim8, i, j, k, l, m, n, o, p, l, j, q, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, k, m, n, o, j, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, l, j, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim3, Dim1, Dim8, p, l, j, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim1, Dim3, Dim8, p, j, l, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim3, Dim8, i, j, k, l, m, n, o, p, j, l, q, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, k, m, n, o, j, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, j, l, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim1, Dim3, Dim8, p, j, l, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim3, Dim8, Dim1, p, l, q, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim8, Dim1, i, j, k, l, m, n, o, p, l, q, j, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, k, m, n, o, j, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, l, q, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim3, Dim8, Dim1, p, l, q, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim1, Dim8, Dim3, p, j, q, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim8, Dim3, i, j, k, l, m, n, o, p, j, q, l, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, k, m, n, o, j, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, j, q, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim1, Dim8, Dim3, p, j, q, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim3, Dim1, p, q, l, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim3, Dim1, i, j, k, l, m, n, o, p, q, l, j, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, k, m, n, o, j, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, l, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim3, Dim1, p, q, l, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim1, Dim3, p, q, j, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim1, Dim3, i, j, k, l, m, n, o, p, q, j, l, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, k, m, n, o, j, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, j, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim1, Dim3, p, q, j, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim7, Dim8, m, j, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim7, Dim8, i, j, k, l, m, n, o, m, j, p, q, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, k, l, n, o, j, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, j, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim7, Dim8, m, j, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim7, Dim8, j, m, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, j, k, l, m, n, o, j, m, p, q, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, k, l, n, o, j, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(j, m, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim7, Dim8, j, m, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim7, Dim1, Dim8, m, p, j, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim1, Dim8, i, j, k, l, m, n, o, m, p, j, q, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, k, l, n, o, j, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, p, j, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim7, Dim1, Dim8, m, p, j, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim7, Dim4, Dim8, j, p, m, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim4, Dim8, i, j, k, l, m, n, o, j, p, m, q, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, k, l, n, o, j, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(j, p, m, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim7, Dim4, Dim8, j, p, m, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim7, Dim8, Dim1, m, p, q, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim8, Dim1, i, j, k, l, m, n, o, m, p, q, j, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, k, l, n, o, j, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, p, q, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim7, Dim8, Dim1, m, p, q, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim7, Dim8, Dim4, j, p, q, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim8, Dim4, i, j, k, l, m, n, o, j, p, q, m, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, k, l, n, o, j, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(j, p, q, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim7, Dim8, Dim4, j, p, q, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim4, Dim1, Dim8, p, m, j, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim1, Dim8, i, j, k, l, m, n, o, p, m, j, q, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, k, l, n, o, j, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, m, j, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim4, Dim1, Dim8, p, m, j, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim1, Dim4, Dim8, p, j, m, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim4, Dim8, i, j, k, l, m, n, o, p, j, m, q, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, k, l, n, o, j, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, j, m, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim1, Dim4, Dim8, p, j, m, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim4, Dim8, Dim1, p, m, q, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim8, Dim1, i, j, k, l, m, n, o, p, m, q, j, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, k, l, n, o, j, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, m, q, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim4, Dim8, Dim1, p, m, q, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim1, Dim8, Dim4, p, j, q, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim8, Dim4, i, j, k, l, m, n, o, p, j, q, m, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, k, l, n, o, j, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, j, q, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim1, Dim8, Dim4, p, j, q, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim4, Dim1, p, q, m, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim4, Dim1, i, j, k, l, m, n, o, p, q, m, j, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, k, l, n, o, j, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, m, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim4, Dim1, p, q, m, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim1, Dim4, p, q, j, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim1, Dim4, i, j, k, l, m, n, o, p, q, j, m, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, k, l, n, o, j, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, j, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim1, Dim4, p, q, j, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim7, Dim8, n, j, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim7, Dim8, i, j, k, l, m, n, o, n, j, p, q, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, Dim8, i, k, l, m, o, j, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, j, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim7, Dim8, n, j, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim7, Dim8, j, n, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim7, Dim8, i, j, k, l, m, n, o, j, n, p, q, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, Dim8, i, k, l, m, o, j, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(j, n, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim7, Dim8, j, n, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim7, Dim1, Dim8, n, p, j, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim1, Dim8, i, j, k, l, m, n, o, n, p, j, q, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, Dim8, i, k, l, m, o, j, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, p, j, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim7, Dim1, Dim8, n, p, j, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim7, Dim5, Dim8, j, p, n, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim5, Dim8, i, j, k, l, m, n, o, j, p, n, q, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, Dim8, i, k, l, m, o, j, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(j, p, n, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim7, Dim5, Dim8, j, p, n, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim7, Dim8, Dim1, n, p, q, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim8, Dim1, i, j, k, l, m, n, o, n, p, q, j, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, Dim8, i, k, l, m, o, j, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, p, q, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim7, Dim8, Dim1, n, p, q, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim7, Dim8, Dim5, j, p, q, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim8, Dim5, i, j, k, l, m, n, o, j, p, q, n, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, Dim8, i, k, l, m, o, j, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(j, p, q, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim7, Dim8, Dim5, j, p, q, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim5, Dim1, Dim8, p, n, j, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim1, Dim8, i, j, k, l, m, n, o, p, n, j, q, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, Dim8, i, k, l, m, o, j, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, n, j, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim5, Dim1, Dim8, p, n, j, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim1, Dim5, Dim8, p, j, n, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim5, Dim8, i, j, k, l, m, n, o, p, j, n, q, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, Dim8, i, k, l, m, o, j, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, j, n, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim1, Dim5, Dim8, p, j, n, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim5, Dim8, Dim1, p, n, q, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim8, Dim1, i, j, k, l, m, n, o, p, n, q, j, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, Dim8, i, k, l, m, o, j, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, n, q, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim5, Dim8, Dim1, p, n, q, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim1, Dim8, Dim5, p, j, q, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim8, Dim5, i, j, k, l, m, n, o, p, j, q, n, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, Dim8, i, k, l, m, o, j, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, j, q, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim1, Dim8, Dim5, p, j, q, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim5, Dim1, p, q, n, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim5, Dim1, i, j, k, l, m, n, o, p, q, n, j, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, Dim8, i, k, l, m, o, j, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, q, n, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim5, Dim1, p, q, n, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim1, Dim5, p, q, j, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim1, Dim5, i, j, k, l, m, n, o, p, q, j, n, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, Dim8, i, k, l, m, o, j, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, q, j, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim1, Dim5, p, q, j, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim7, Dim8, o, j, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim1, Dim7, Dim8, i, j, k, l, m, n, o, o, j, p, q, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, i, k, l, m, n, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, i, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(o, j, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim7, Dim8, o, j, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim7, Dim8, j, o, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, j, o, p, q, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, i, k, l, m, n, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, i, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(j, o, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim7, Dim8, j, o, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim1, Dim8, o, p, j, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim1, Dim8, i, j, k, l, m, n, o, o, p, j, q, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, i, k, l, m, n, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, i, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(o, p, j, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim1, Dim8, o, p, j, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, o, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim7, Dim6, Dim8, j, p, o, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim6, Dim8, i, j, k, l, m, n, o, j, p, o, q, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, i, k, l, m, n, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, i, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(j, p, o, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim7, Dim6, Dim8, j, p, o, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim1, o, p, q, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim8, Dim1, i, j, k, l, m, n, o, o, p, q, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, i, k, l, m, n, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, i, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(o, p, q, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim1, o, p, q, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, q, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim1, Dim7, Dim8, Dim6, j, p, q, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim8, Dim6, i, j, k, l, m, n, o, j, p, q, o, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, i, k, l, m, n, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, i, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(j, p, q, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim7, Dim8, Dim6, j, p, q, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim6, Dim1, Dim8, p, o, j, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim1, Dim8, i, j, k, l, m, n, o, p, o, j, q, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, i, k, l, m, n, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, i, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, o, j, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim6, Dim1, Dim8, p, o, j, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, o, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim1, Dim6, Dim8, p, j, o, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim6, Dim8, i, j, k, l, m, n, o, p, j, o, q, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, i, k, l, m, n, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, i, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, j, o, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim1, Dim6, Dim8, p, j, o, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim6, Dim8, Dim1, p, o, q, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim8, Dim1, i, j, k, l, m, n, o, p, o, q, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, i, k, l, m, n, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, i, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, o, q, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim6, Dim8, Dim1, p, o, q, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, q, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim1, Dim8, Dim6, p, j, q, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim8, Dim6, i, j, k, l, m, n, o, p, j, q, o, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, i, k, l, m, n, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, i, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, j, q, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim1, Dim8, Dim6, p, j, q, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim6, Dim1, p, q, o, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim6, Dim1, i, j, k, l, m, n, o, p, q, o, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, i, k, l, m, n, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, i, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, q, o, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim6, Dim1, p, q, o, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim1, Dim6, p, q, j, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim1, Dim6, i, j, k, l, m, n, o, p, q, j, o, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, i, k, l, m, n, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, Dim8, i, k, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, q, j, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim1, Dim6, p, q, j, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim7, Dim8, l, k, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim7, Dim8, i, j, k, l, m, n, o, l, k, p, q, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, m, n, o, k, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, k, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim7, Dim8, l, k, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim7, Dim8, k, l, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, k, l, m, n, o, k, l, p, q, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, m, n, o, k, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, l, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim7, Dim8, k, l, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim7, Dim2, Dim8, l, p, k, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim2, Dim8, i, j, k, l, m, n, o, l, p, k, q, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, m, n, o, k, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, p, k, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim7, Dim2, Dim8, l, p, k, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim7, Dim3, Dim8, k, p, l, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim3, Dim8, i, j, k, l, m, n, o, k, p, l, q, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, m, n, o, k, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, p, l, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim7, Dim3, Dim8, k, p, l, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim7, Dim8, Dim2, l, p, q, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim8, Dim2, i, j, k, l, m, n, o, l, p, q, k, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, m, n, o, k, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, p, q, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim7, Dim8, Dim2, l, p, q, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim7, Dim8, Dim3, k, p, q, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim8, Dim3, i, j, k, l, m, n, o, k, p, q, l, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, m, n, o, k, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, p, q, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim7, Dim8, Dim3, k, p, q, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim3, Dim2, Dim8, p, l, k, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim2, Dim8, i, j, k, l, m, n, o, p, l, k, q, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, m, n, o, k, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, l, k, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim3, Dim2, Dim8, p, l, k, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim2, Dim3, Dim8, p, k, l, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim3, Dim8, i, j, k, l, m, n, o, p, k, l, q, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, m, n, o, k, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, k, l, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim2, Dim3, Dim8, p, k, l, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim3, Dim8, Dim2, p, l, q, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim8, Dim2, i, j, k, l, m, n, o, p, l, q, k, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, m, n, o, k, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, l, q, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim3, Dim8, Dim2, p, l, q, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim2, Dim8, Dim3, p, k, q, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim8, Dim3, i, j, k, l, m, n, o, p, k, q, l, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, m, n, o, k, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, k, q, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim2, Dim8, Dim3, p, k, q, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim3, Dim2, p, q, l, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim3, Dim2, i, j, k, l, m, n, o, p, q, l, k, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, m, n, o, k, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, l, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim3, Dim2, p, q, l, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim2, Dim3, p, q, k, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim2, Dim3, i, j, k, l, m, n, o, p, q, k, l, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, m, n, o, k, l, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, k, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim2, Dim3, p, q, k, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim7, Dim8, m, k, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim7, Dim8, i, j, k, l, m, n, o, m, k, p, q, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, l, n, o, k, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, k, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim7, Dim8, m, k, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim7, Dim8, k, m, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, k, l, m, n, o, k, m, p, q, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, l, n, o, k, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(k, m, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim7, Dim8, k, m, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim7, Dim2, Dim8, m, p, k, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim2, Dim8, i, j, k, l, m, n, o, m, p, k, q, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, l, n, o, k, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, p, k, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim7, Dim2, Dim8, m, p, k, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim7, Dim4, Dim8, k, p, m, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim4, Dim8, i, j, k, l, m, n, o, k, p, m, q, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, l, n, o, k, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(k, p, m, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim7, Dim4, Dim8, k, p, m, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim7, Dim8, Dim2, m, p, q, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim8, Dim2, i, j, k, l, m, n, o, m, p, q, k, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, l, n, o, k, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, p, q, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim7, Dim8, Dim2, m, p, q, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim7, Dim8, Dim4, k, p, q, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim8, Dim4, i, j, k, l, m, n, o, k, p, q, m, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, l, n, o, k, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(k, p, q, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim7, Dim8, Dim4, k, p, q, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim4, Dim2, Dim8, p, m, k, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim2, Dim8, i, j, k, l, m, n, o, p, m, k, q, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, l, n, o, k, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, m, k, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim4, Dim2, Dim8, p, m, k, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim2, Dim4, Dim8, p, k, m, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim4, Dim8, i, j, k, l, m, n, o, p, k, m, q, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, l, n, o, k, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, k, m, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim2, Dim4, Dim8, p, k, m, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim4, Dim8, Dim2, p, m, q, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim8, Dim2, i, j, k, l, m, n, o, p, m, q, k, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, l, n, o, k, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, m, q, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim4, Dim8, Dim2, p, m, q, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim2, Dim8, Dim4, p, k, q, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim8, Dim4, i, j, k, l, m, n, o, p, k, q, m, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, l, n, o, k, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, k, q, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim2, Dim8, Dim4, p, k, q, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim4, Dim2, p, q, m, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim4, Dim2, i, j, k, l, m, n, o, p, q, m, k, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, l, n, o, k, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, m, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim4, Dim2, p, q, m, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim2, Dim4, p, q, k, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim2, Dim4, i, j, k, l, m, n, o, p, q, k, m, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, l, n, o, k, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, k, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim2, Dim4, p, q, k, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim7, Dim8, n, k, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim7, Dim8, i, j, k, l, m, n, o, n, k, p, q, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, l, m, o, k, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, k, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim7, Dim8, n, k, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim7, Dim8, k, n, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, k, l, m, n, o, k, n, p, q, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, l, m, o, k, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(k, n, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim7, Dim8, k, n, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim7, Dim2, Dim8, n, p, k, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim2, Dim8, i, j, k, l, m, n, o, n, p, k, q, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, l, m, o, k, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, p, k, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim7, Dim2, Dim8, n, p, k, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim7, Dim5, Dim8, k, p, n, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim5, Dim8, i, j, k, l, m, n, o, k, p, n, q, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, l, m, o, k, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(k, p, n, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim7, Dim5, Dim8, k, p, n, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim7, Dim8, Dim2, n, p, q, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim8, Dim2, i, j, k, l, m, n, o, n, p, q, k, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, l, m, o, k, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, p, q, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim7, Dim8, Dim2, n, p, q, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim7, Dim8, Dim5, k, p, q, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim8, Dim5, i, j, k, l, m, n, o, k, p, q, n, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, l, m, o, k, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(k, p, q, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim7, Dim8, Dim5, k, p, q, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim5, Dim2, Dim8, p, n, k, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim2, Dim8, i, j, k, l, m, n, o, p, n, k, q, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, l, m, o, k, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, n, k, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim5, Dim2, Dim8, p, n, k, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim2, Dim5, Dim8, p, k, n, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim5, Dim8, i, j, k, l, m, n, o, p, k, n, q, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, l, m, o, k, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, k, n, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim2, Dim5, Dim8, p, k, n, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim5, Dim8, Dim2, p, n, q, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim8, Dim2, i, j, k, l, m, n, o, p, n, q, k, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, l, m, o, k, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, n, q, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim5, Dim8, Dim2, p, n, q, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim2, Dim8, Dim5, p, k, q, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim8, Dim5, i, j, k, l, m, n, o, p, k, q, n, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, l, m, o, k, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, k, q, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim2, Dim8, Dim5, p, k, q, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim5, Dim2, p, q, n, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim5, Dim2, i, j, k, l, m, n, o, p, q, n, k, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, l, m, o, k, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, q, n, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim5, Dim2, p, q, n, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim2, Dim5, p, q, k, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim2, Dim5, i, j, k, l, m, n, o, p, q, k, n, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, l, m, o, k, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(p, q, k, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim2, Dim5, p, q, k, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim7, Dim8, o, k, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim2, Dim7, Dim8, i, j, k, l, m, n, o, o, k, p, q, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, l, m, n, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, Dim8, i, j, l, m, n, p, q>(TensorExpr(a, b));
}

//B(o, k, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim7, Dim8, o, k, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim7, Dim8, k, o, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, k, o, p, q, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, l, m, n, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, Dim8, i, j, l, m, n, p, q>(TensorExpr(a, b));
}

//B(k, o, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim7, Dim8, k, o, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim2, Dim8, o, p, k, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim2, Dim8, i, j, k, l, m, n, o, o, p, k, q, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, l, m, n, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, Dim8, i, j, l, m, n, p, q>(TensorExpr(a, b));
}

//B(o, p, k, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim2, Dim8, o, p, k, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, o, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim7, Dim6, Dim8, k, p, o, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim6, Dim8, i, j, k, l, m, n, o, k, p, o, q, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, l, m, n, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, Dim8, i, j, l, m, n, p, q>(TensorExpr(a, b));
}

//B(k, p, o, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim7, Dim6, Dim8, k, p, o, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim2, o, p, q, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim8, Dim2, i, j, k, l, m, n, o, o, p, q, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, l, m, n, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, Dim8, i, j, l, m, n, p, q>(TensorExpr(a, b));
}

//B(o, p, q, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim2, o, p, q, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, q, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim2, Dim7, Dim8, Dim6, k, p, q, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim8, Dim6, i, j, k, l, m, n, o, k, p, q, o, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, l, m, n, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, Dim8, i, j, l, m, n, p, q>(TensorExpr(a, b));
}

//B(k, p, q, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim7, Dim8, Dim6, k, p, q, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim6, Dim2, Dim8, p, o, k, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim2, Dim8, i, j, k, l, m, n, o, p, o, k, q, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, l, m, n, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, Dim8, i, j, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, o, k, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim6, Dim2, Dim8, p, o, k, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, o, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim2, Dim6, Dim8, p, k, o, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim6, Dim8, i, j, k, l, m, n, o, p, k, o, q, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, l, m, n, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, Dim8, i, j, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, k, o, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim2, Dim6, Dim8, p, k, o, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim6, Dim8, Dim2, p, o, q, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim8, Dim2, i, j, k, l, m, n, o, p, o, q, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, l, m, n, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, Dim8, i, j, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, o, q, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim6, Dim8, Dim2, p, o, q, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, q, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim2, Dim8, Dim6, p, k, q, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim8, Dim6, i, j, k, l, m, n, o, p, k, q, o, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, l, m, n, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, Dim8, i, j, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, k, q, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim2, Dim8, Dim6, p, k, q, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim6, Dim2, p, q, o, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim6, Dim2, i, j, k, l, m, n, o, p, q, o, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, l, m, n, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, Dim8, i, j, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, q, o, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim6, Dim2, p, q, o, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim2, Dim6, p, q, k, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim2, Dim6, i, j, k, l, m, n, o, p, q, k, o, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, l, m, n, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, Dim8, i, j, l, m, n, p, q>(TensorExpr(a, b));
}

//B(p, q, k, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim2, Dim6, p, q, k, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim7, Dim8, m, l, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim7, Dim8, i, j, k, l, m, n, o, m, l, p, q, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, n, o, l, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(m, l, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim7, Dim8, m, l, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim7, Dim8, l, m, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, n, o, l, m, p, q, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, n, o, l, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(l, m, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim7, Dim8, l, m, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim7, Dim3, Dim8, m, p, l, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim3, Dim8, i, j, k, l, m, n, o, m, p, l, q, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, n, o, l, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(m, p, l, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim7, Dim3, Dim8, m, p, l, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim7, Dim4, Dim8, l, p, m, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim4, Dim8, i, j, k, l, m, n, o, l, p, m, q, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, n, o, l, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(l, p, m, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim7, Dim4, Dim8, l, p, m, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim7, Dim8, Dim3, m, p, q, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim8, Dim3, i, j, k, l, m, n, o, m, p, q, l, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, n, o, l, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(m, p, q, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim7, Dim8, Dim3, m, p, q, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim7, Dim8, Dim4, l, p, q, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim8, Dim4, i, j, k, l, m, n, o, l, p, q, m, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, n, o, l, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(l, p, q, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim7, Dim8, Dim4, l, p, q, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim4, Dim3, Dim8, p, m, l, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim3, Dim8, i, j, k, l, m, n, o, p, m, l, q, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, n, o, l, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(p, m, l, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim4, Dim3, Dim8, p, m, l, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim3, Dim4, Dim8, p, l, m, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim4, Dim8, i, j, k, l, m, n, o, p, l, m, q, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, n, o, l, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(p, l, m, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim3, Dim4, Dim8, p, l, m, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim4, Dim8, Dim3, p, m, q, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim8, Dim3, i, j, k, l, m, n, o, p, m, q, l, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, n, o, l, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(p, m, q, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim4, Dim8, Dim3, p, m, q, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim3, Dim8, Dim4, p, l, q, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim8, Dim4, i, j, k, l, m, n, o, p, l, q, m, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, n, o, l, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(p, l, q, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim3, Dim8, Dim4, p, l, q, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim4, Dim3, p, q, m, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim4, Dim3, i, j, k, l, m, n, o, p, q, m, l, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, n, o, l, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, m, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim4, Dim3, p, q, m, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim3, Dim4, p, q, l, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim3, Dim4, i, j, k, l, m, n, o, p, q, l, m, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, n, o, l, m, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(p, q, l, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim3, Dim4, p, q, l, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim7, Dim8, n, l, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim7, Dim8, i, j, k, l, m, n, o, n, l, p, q, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, m, o, l, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(n, l, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim7, Dim8, n, l, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim7, Dim8, l, n, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, l, m, n, o, l, n, p, q, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, m, o, l, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(l, n, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim7, Dim8, l, n, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim7, Dim3, Dim8, n, p, l, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim3, Dim8, i, j, k, l, m, n, o, n, p, l, q, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, m, o, l, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(n, p, l, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim7, Dim3, Dim8, n, p, l, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim7, Dim5, Dim8, l, p, n, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim5, Dim8, i, j, k, l, m, n, o, l, p, n, q, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, m, o, l, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(l, p, n, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim7, Dim5, Dim8, l, p, n, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim7, Dim8, Dim3, n, p, q, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim8, Dim3, i, j, k, l, m, n, o, n, p, q, l, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, m, o, l, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(n, p, q, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim7, Dim8, Dim3, n, p, q, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim7, Dim8, Dim5, l, p, q, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim8, Dim5, i, j, k, l, m, n, o, l, p, q, n, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, m, o, l, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(l, p, q, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim7, Dim8, Dim5, l, p, q, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim5, Dim3, Dim8, p, n, l, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim3, Dim8, i, j, k, l, m, n, o, p, n, l, q, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, m, o, l, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(p, n, l, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim5, Dim3, Dim8, p, n, l, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim3, Dim5, Dim8, p, l, n, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim5, Dim8, i, j, k, l, m, n, o, p, l, n, q, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, m, o, l, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(p, l, n, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim3, Dim5, Dim8, p, l, n, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim5, Dim8, Dim3, p, n, q, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim8, Dim3, i, j, k, l, m, n, o, p, n, q, l, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, m, o, l, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(p, n, q, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim5, Dim8, Dim3, p, n, q, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim3, Dim8, Dim5, p, l, q, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim8, Dim5, i, j, k, l, m, n, o, p, l, q, n, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, m, o, l, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(p, l, q, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim3, Dim8, Dim5, p, l, q, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim5, Dim3, p, q, n, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim5, Dim3, i, j, k, l, m, n, o, p, q, n, l, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, m, o, l, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(p, q, n, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim5, Dim3, p, q, n, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim3, Dim5, p, q, l, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim3, Dim5, i, j, k, l, m, n, o, p, q, l, n, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, m, o, l, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(p, q, l, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim3, Dim5, p, q, l, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim7, Dim8, o, l, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim3, Dim7, Dim8, i, j, k, l, m, n, o, o, l, p, q, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, m, n, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, Dim8, i, j, k, m, n, p, q>(TensorExpr(a, b));
}

//B(o, l, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim7, Dim8, o, l, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim7, Dim8, l, o, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, l, o, p, q, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, m, n, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, Dim8, i, j, k, m, n, p, q>(TensorExpr(a, b));
}

//B(l, o, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim7, Dim8, l, o, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim3, Dim8, o, p, l, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim3, Dim8, i, j, k, l, m, n, o, o, p, l, q, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, m, n, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, Dim8, i, j, k, m, n, p, q>(TensorExpr(a, b));
}

//B(o, p, l, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim3, Dim8, o, p, l, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, o, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim7, Dim6, Dim8, l, p, o, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim6, Dim8, i, j, k, l, m, n, o, l, p, o, q, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, m, n, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, Dim8, i, j, k, m, n, p, q>(TensorExpr(a, b));
}

//B(l, p, o, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim7, Dim6, Dim8, l, p, o, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim3, o, p, q, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim8, Dim3, i, j, k, l, m, n, o, o, p, q, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, m, n, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, Dim8, i, j, k, m, n, p, q>(TensorExpr(a, b));
}

//B(o, p, q, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim3, o, p, q, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, q, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim3, Dim7, Dim8, Dim6, l, p, q, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim8, Dim6, i, j, k, l, m, n, o, l, p, q, o, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, m, n, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, Dim8, i, j, k, m, n, p, q>(TensorExpr(a, b));
}

//B(l, p, q, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim7, Dim8, Dim6, l, p, q, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim6, Dim3, Dim8, p, o, l, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim3, Dim8, i, j, k, l, m, n, o, p, o, l, q, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, m, n, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, Dim8, i, j, k, m, n, p, q>(TensorExpr(a, b));
}

//B(p, o, l, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim6, Dim3, Dim8, p, o, l, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, o, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim3, Dim6, Dim8, p, l, o, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim6, Dim8, i, j, k, l, m, n, o, p, l, o, q, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, m, n, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, Dim8, i, j, k, m, n, p, q>(TensorExpr(a, b));
}

//B(p, l, o, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim3, Dim6, Dim8, p, l, o, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim6, Dim8, Dim3, p, o, q, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim8, Dim3, i, j, k, l, m, n, o, p, o, q, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, m, n, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, Dim8, i, j, k, m, n, p, q>(TensorExpr(a, b));
}

//B(p, o, q, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim6, Dim8, Dim3, p, o, q, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, q, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim3, Dim8, Dim6, p, l, q, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim8, Dim6, i, j, k, l, m, n, o, p, l, q, o, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, m, n, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, Dim8, i, j, k, m, n, p, q>(TensorExpr(a, b));
}

//B(p, l, q, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim3, Dim8, Dim6, p, l, q, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim6, Dim3, p, q, o, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim6, Dim3, i, j, k, l, m, n, o, p, q, o, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, m, n, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, Dim8, i, j, k, m, n, p, q>(TensorExpr(a, b));
}

//B(p, q, o, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim6, Dim3, p, q, o, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim3, Dim6, p, q, l, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim3, Dim6, i, j, k, l, m, n, o, p, q, l, o, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, m, n, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, Dim8, i, j, k, m, n, p, q>(TensorExpr(a, b));
}

//B(p, q, l, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim3, Dim6, p, q, l, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim7, Dim8, n, m, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim7, Dim8, i, j, k, l, m, n, o, n, m, p, q, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, o, m, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(n, m, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim7, Dim8, n, m, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim7, Dim8, m, n, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, m, n, o, m, n, p, q, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, o, m, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(m, n, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim7, Dim8, m, n, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim7, Dim4, Dim8, n, p, m, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim4, Dim8, i, j, k, l, m, n, o, n, p, m, q, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, o, m, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(n, p, m, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim7, Dim4, Dim8, n, p, m, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim7, Dim5, Dim8, m, p, n, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim5, Dim8, i, j, k, l, m, n, o, m, p, n, q, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, o, m, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(m, p, n, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim7, Dim5, Dim8, m, p, n, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim7, Dim8, Dim4, n, p, q, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim8, Dim4, i, j, k, l, m, n, o, n, p, q, m, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, o, m, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(n, p, q, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim7, Dim8, Dim4, n, p, q, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim7, Dim8, Dim5, m, p, q, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim8, Dim5, i, j, k, l, m, n, o, m, p, q, n, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, o, m, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(m, p, q, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim7, Dim8, Dim5, m, p, q, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim5, Dim4, Dim8, p, n, m, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim4, Dim8, i, j, k, l, m, n, o, p, n, m, q, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, o, m, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(p, n, m, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim5, Dim4, Dim8, p, n, m, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim4, Dim5, Dim8, p, m, n, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim5, Dim8, i, j, k, l, m, n, o, p, m, n, q, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, o, m, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(p, m, n, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim4, Dim5, Dim8, p, m, n, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim5, Dim8, Dim4, p, n, q, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim8, Dim4, i, j, k, l, m, n, o, p, n, q, m, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, o, m, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(p, n, q, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim5, Dim8, Dim4, p, n, q, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim4, Dim8, Dim5, p, m, q, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim8, Dim5, i, j, k, l, m, n, o, p, m, q, n, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, o, m, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(p, m, q, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim4, Dim8, Dim5, p, m, q, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim5, Dim4, p, q, n, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim5, Dim4, i, j, k, l, m, n, o, p, q, n, m, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, o, m, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(p, q, n, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim5, Dim4, p, q, n, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim4, Dim5, p, q, m, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim4, Dim5, i, j, k, l, m, n, o, p, q, m, n, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, o, m, n, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(p, q, m, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim4, Dim5, p, q, m, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim7, Dim8, o, m, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim4, Dim7, Dim8, i, j, k, l, m, n, o, o, m, p, q, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, n, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, Dim8, i, j, k, l, n, p, q>(TensorExpr(a, b));
}

//B(o, m, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim7, Dim8, o, m, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim7, Dim8, m, o, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, m, o, p, q, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, n, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, Dim8, i, j, k, l, n, p, q>(TensorExpr(a, b));
}

//B(m, o, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim7, Dim8, m, o, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim4, Dim8, o, p, m, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim4, Dim8, i, j, k, l, m, n, o, o, p, m, q, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, n, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, Dim8, i, j, k, l, n, p, q>(TensorExpr(a, b));
}

//B(o, p, m, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim4, Dim8, o, p, m, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, o, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim7, Dim6, Dim8, m, p, o, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim6, Dim8, i, j, k, l, m, n, o, m, p, o, q, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, n, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, Dim8, i, j, k, l, n, p, q>(TensorExpr(a, b));
}

//B(m, p, o, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim7, Dim6, Dim8, m, p, o, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim4, o, p, q, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim8, Dim4, i, j, k, l, m, n, o, o, p, q, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, n, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, Dim8, i, j, k, l, n, p, q>(TensorExpr(a, b));
}

//B(o, p, q, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim4, o, p, q, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, q, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim4, Dim7, Dim8, Dim6, m, p, q, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim8, Dim6, i, j, k, l, m, n, o, m, p, q, o, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, n, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, Dim8, i, j, k, l, n, p, q>(TensorExpr(a, b));
}

//B(m, p, q, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim7, Dim8, Dim6, m, p, q, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim6, Dim4, Dim8, p, o, m, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim4, Dim8, i, j, k, l, m, n, o, p, o, m, q, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, n, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, Dim8, i, j, k, l, n, p, q>(TensorExpr(a, b));
}

//B(p, o, m, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim6, Dim4, Dim8, p, o, m, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, o, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim4, Dim6, Dim8, p, m, o, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim6, Dim8, i, j, k, l, m, n, o, p, m, o, q, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, n, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, Dim8, i, j, k, l, n, p, q>(TensorExpr(a, b));
}

//B(p, m, o, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim4, Dim6, Dim8, p, m, o, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim6, Dim8, Dim4, p, o, q, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim8, Dim4, i, j, k, l, m, n, o, p, o, q, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, n, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, Dim8, i, j, k, l, n, p, q>(TensorExpr(a, b));
}

//B(p, o, q, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim6, Dim8, Dim4, p, o, q, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, q, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim4, Dim8, Dim6, p, m, q, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim8, Dim6, i, j, k, l, m, n, o, p, m, q, o, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, n, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, Dim8, i, j, k, l, n, p, q>(TensorExpr(a, b));
}

//B(p, m, q, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim4, Dim8, Dim6, p, m, q, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim6, Dim4, p, q, o, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim6, Dim4, i, j, k, l, m, n, o, p, q, o, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, n, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, Dim8, i, j, k, l, n, p, q>(TensorExpr(a, b));
}

//B(p, q, o, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim6, Dim4, p, q, o, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim4, Dim6, p, q, m, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim4, Dim6, i, j, k, l, m, n, o, p, q, m, o, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, n, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, Dim8, i, j, k, l, n, p, q>(TensorExpr(a, b));
}

//B(p, q, m, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim4, Dim6, p, q, m, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim7, Dim8, o, n, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim5, Dim7, Dim8, i, j, k, l, m, n, o, o, n, p, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, p, q>(TensorExpr(a, b));
}

//B(o, n, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim7, Dim8, o, n, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim8, n, o, p, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, n, o, p, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, p, q>(TensorExpr(a, b));
}

//B(n, o, p, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim8, n, o, p, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim5, Dim8, o, p, n, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim5, Dim8, i, j, k, l, m, n, o, o, p, n, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, p, q>(TensorExpr(a, b));
}

//B(o, p, n, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim5, Dim8, o, p, n, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, o, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim7, Dim6, Dim8, n, p, o, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim6, Dim8, i, j, k, l, m, n, o, n, p, o, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, p, q>(TensorExpr(a, b));
}

//B(n, p, o, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim7, Dim6, Dim8, n, p, o, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim5, o, p, q, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim8, Dim5, i, j, k, l, m, n, o, o, p, q, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, p, q>(TensorExpr(a, b));
}

//B(o, p, q, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim5, o, p, q, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, q, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim5, Dim7, Dim8, Dim6, n, p, q, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim8, Dim6, i, j, k, l, m, n, o, n, p, q, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, p, q>(TensorExpr(a, b));
}

//B(n, p, q, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim7, Dim8, Dim6, n, p, q, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim6, Dim5, Dim8, p, o, n, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim5, Dim8, i, j, k, l, m, n, o, p, o, n, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, p, q>(TensorExpr(a, b));
}

//B(p, o, n, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim6, Dim5, Dim8, p, o, n, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, o, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim5, Dim6, Dim8, p, n, o, q> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim6, Dim8, i, j, k, l, m, n, o, p, n, o, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, p, q>(TensorExpr(a, b));
}

//B(p, n, o, q, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim5, Dim6, Dim8, p, n, o, q> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim6, Dim8, Dim5, p, o, q, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim8, Dim5, i, j, k, l, m, n, o, p, o, q, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, p, q>(TensorExpr(a, b));
}

//B(p, o, q, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim6, Dim8, Dim5, p, o, q, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, q, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim5, Dim8, Dim6, p, n, q, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim8, Dim6, i, j, k, l, m, n, o, p, n, q, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, p, q>(TensorExpr(a, b));
}

//B(p, n, q, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim5, Dim8, Dim6, p, n, q, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim6, Dim5, p, q, o, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim6, Dim5, i, j, k, l, m, n, o, p, q, o, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, p, q>(TensorExpr(a, b));
}

//B(p, q, o, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim6, Dim5, p, q, o, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, q, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor4_Expr<B, U, Dim7, Dim8, Dim5, Dim6, p, q, n, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim5, Dim6, i, j, k, l, m, n, o, p, q, n, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, p, q>(TensorExpr(a, b));
}

//B(p, q, n, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor4_Expr<B, U, Dim7, Dim8, Dim5, Dim6, p, q, n, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

} //namespace FTensor

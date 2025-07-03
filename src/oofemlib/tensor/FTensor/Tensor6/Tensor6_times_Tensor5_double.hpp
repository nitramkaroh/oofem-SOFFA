#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, n, )*B(o, p, q, r, s, ) double contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, int Dim9, int Dim10, char i, char j, char k, char l, char m, char n, char o, char p, char q, char r, char s, int DimA, int DimB, int DimC, int DimD, int DimX, int DimY, int DimZ, int DimW, int DimV, char a, char b, char c, char d, char x, char y, char z, char w, char v>
class Tensor6_times_Tensor5_double
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim9, Dim10, o, p, q, r, s> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const int N7) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        {
          // Permutation is where the indices get checked.
          result += Permutation6<DimA, DimB, DimC, DimD, DimX, DimY, a, b, c, d, x, y>().eval(
            iterA, N1, N2, N3, N4, xx, yy) * Permutation5<DimX, DimY, DimZ, DimW, DimV, x, y, z, w, v>().eval(iterB, xx, yy, N5, N6, N7);
        }
    return result;
  }

  Tensor6_times_Tensor5_double(
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &iter_a,
    const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim9, Dim10, o, p, q, r, s> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, n, )*B(j, i, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim0, Dim6, Dim7, Dim8, j, i, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, Dim7, Dim8, i, j, k, l, m, n, j, i, o, p, q, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, i, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim0, Dim6, Dim7, Dim8, j, i, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim1, Dim6, Dim7, Dim8, i, j, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, i, j, k, l, m, n, i, j, o, p, q, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, j, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim1, Dim6, Dim7, Dim8, i, j, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim0, Dim7, Dim8, j, o, i, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, Dim7, Dim8, i, j, k, l, m, n, j, o, i, p, q, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, o, i, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim0, Dim7, Dim8, j, o, i, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim1, Dim7, Dim8, i, o, j, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, Dim7, Dim8, i, j, k, l, m, n, i, o, j, p, q, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, o, j, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim1, Dim7, Dim8, i, o, j, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim0, Dim8, j, o, p, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim0, Dim8, i, j, k, l, m, n, j, o, p, i, q, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, o, p, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim0, Dim8, j, o, p, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim1, Dim8, i, o, p, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim1, Dim8, i, j, k, l, m, n, i, o, p, j, q, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, o, p, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim1, Dim8, i, o, p, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim8, Dim0, j, o, p, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, Dim0, i, j, k, l, m, n, j, o, p, q, i, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, o, p, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim8, Dim0, j, o, p, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim8, Dim1, i, o, p, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, Dim1, i, j, k, l, m, n, i, o, p, q, j, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, o, p, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim8, Dim1, i, o, p, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim0, Dim7, Dim8, o, j, i, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim7, Dim8, i, j, k, l, m, n, o, j, i, p, q, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, j, i, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim0, Dim7, Dim8, o, j, i, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim1, Dim7, Dim8, o, i, j, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, Dim8, i, j, k, l, m, n, o, i, j, p, q, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, i, j, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim1, Dim7, Dim8, o, i, j, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim0, Dim8, o, j, p, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim0, Dim8, i, j, k, l, m, n, o, j, p, i, q, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, j, p, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim0, Dim8, o, j, p, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim1, Dim8, o, i, p, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim1, Dim8, i, j, k, l, m, n, o, i, p, j, q, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, i, p, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim1, Dim8, o, i, p, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim8, Dim0, o, j, p, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim8, Dim0, i, j, k, l, m, n, o, j, p, q, i, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, j, p, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim8, Dim0, o, j, p, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim8, Dim1, o, i, p, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim8, Dim1, i, j, k, l, m, n, o, i, p, q, j, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, i, p, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim8, Dim1, o, i, p, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim0, Dim8, o, p, j, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim0, Dim8, i, j, k, l, m, n, o, p, j, i, q, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, j, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim0, Dim8, o, p, j, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim1, Dim8, o, p, i, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim1, Dim8, i, j, k, l, m, n, o, p, i, j, q, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, i, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim1, Dim8, o, p, i, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim8, Dim0, o, p, j, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim8, Dim0, i, j, k, l, m, n, o, p, j, q, i, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, j, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim8, Dim0, o, p, j, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim8, Dim1, o, p, i, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim8, Dim1, i, j, k, l, m, n, o, p, i, q, j, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, i, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim8, Dim1, o, p, i, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim1, Dim0, o, p, q, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim1, Dim0, i, j, k, l, m, n, o, p, q, j, i, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim1, Dim0, o, p, q, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim0, Dim1, o, p, q, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim0, Dim1, i, j, k, l, m, n, o, p, q, i, j, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, Dim8, k, l, m, n, i, j, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, k, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim0, Dim1, o, p, q, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim0, Dim6, Dim7, Dim8, k, i, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, Dim7, Dim8, i, j, k, l, m, n, k, i, o, p, q, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, i, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim0, Dim6, Dim7, Dim8, k, i, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim2, Dim6, Dim7, Dim8, i, k, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, i, j, k, l, m, n, i, k, o, p, q, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, k, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim2, Dim6, Dim7, Dim8, i, k, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim0, Dim7, Dim8, k, o, i, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, Dim7, Dim8, i, j, k, l, m, n, k, o, i, p, q, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, o, i, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim0, Dim7, Dim8, k, o, i, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim2, Dim7, Dim8, i, o, k, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, Dim7, Dim8, i, j, k, l, m, n, i, o, k, p, q, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, o, k, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim2, Dim7, Dim8, i, o, k, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim0, Dim8, k, o, p, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim0, Dim8, i, j, k, l, m, n, k, o, p, i, q, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, o, p, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim0, Dim8, k, o, p, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim2, Dim8, i, o, p, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim2, Dim8, i, j, k, l, m, n, i, o, p, k, q, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, o, p, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim2, Dim8, i, o, p, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim8, Dim0, k, o, p, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, Dim0, i, j, k, l, m, n, k, o, p, q, i, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, o, p, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim8, Dim0, k, o, p, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim8, Dim2, i, o, p, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, Dim2, i, j, k, l, m, n, i, o, p, q, k, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, o, p, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim8, Dim2, i, o, p, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim0, Dim7, Dim8, o, k, i, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim7, Dim8, i, j, k, l, m, n, o, k, i, p, q, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, k, i, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim0, Dim7, Dim8, o, k, i, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim2, Dim7, Dim8, o, i, k, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, Dim8, i, j, k, l, m, n, o, i, k, p, q, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, i, k, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim2, Dim7, Dim8, o, i, k, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim0, Dim8, o, k, p, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim0, Dim8, i, j, k, l, m, n, o, k, p, i, q, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, k, p, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim0, Dim8, o, k, p, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim2, Dim8, o, i, p, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim2, Dim8, i, j, k, l, m, n, o, i, p, k, q, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, i, p, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim2, Dim8, o, i, p, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim8, Dim0, o, k, p, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim8, Dim0, i, j, k, l, m, n, o, k, p, q, i, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, k, p, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim8, Dim0, o, k, p, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim8, Dim2, o, i, p, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim8, Dim2, i, j, k, l, m, n, o, i, p, q, k, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, i, p, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim8, Dim2, o, i, p, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim0, Dim8, o, p, k, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim0, Dim8, i, j, k, l, m, n, o, p, k, i, q, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, k, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim0, Dim8, o, p, k, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim2, Dim8, o, p, i, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim2, Dim8, i, j, k, l, m, n, o, p, i, k, q, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, i, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim2, Dim8, o, p, i, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim8, Dim0, o, p, k, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim8, Dim0, i, j, k, l, m, n, o, p, k, q, i, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, k, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim8, Dim0, o, p, k, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim8, Dim2, o, p, i, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim8, Dim2, i, j, k, l, m, n, o, p, i, q, k, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, i, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim8, Dim2, o, p, i, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim2, Dim0, o, p, q, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim2, Dim0, i, j, k, l, m, n, o, p, q, k, i, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim2, Dim0, o, p, q, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim0, Dim2, o, p, q, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim0, Dim2, i, j, k, l, m, n, o, p, q, i, k, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, Dim8, j, l, m, n, i, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, j, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim0, Dim2, o, p, q, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim0, Dim6, Dim7, Dim8, l, i, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim6, Dim7, Dim8, i, j, k, l, m, n, l, i, o, p, q, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, i, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim0, Dim6, Dim7, Dim8, l, i, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim3, Dim6, Dim7, Dim8, i, l, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, i, j, k, l, m, n, i, l, o, p, q, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, l, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim3, Dim6, Dim7, Dim8, i, l, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim0, Dim7, Dim8, l, o, i, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim0, Dim7, Dim8, i, j, k, l, m, n, l, o, i, p, q, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, o, i, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim0, Dim7, Dim8, l, o, i, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim3, Dim7, Dim8, i, o, l, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim3, Dim7, Dim8, i, j, k, l, m, n, i, o, l, p, q, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, o, l, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim3, Dim7, Dim8, i, o, l, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim0, Dim8, l, o, p, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim0, Dim8, i, j, k, l, m, n, l, o, p, i, q, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, o, p, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim0, Dim8, l, o, p, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim3, Dim8, i, o, p, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim3, Dim8, i, j, k, l, m, n, i, o, p, l, q, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, o, p, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim3, Dim8, i, o, p, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim8, Dim0, l, o, p, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, Dim0, i, j, k, l, m, n, l, o, p, q, i, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, o, p, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim8, Dim0, l, o, p, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim8, Dim3, i, o, p, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, Dim3, i, j, k, l, m, n, i, o, p, q, l, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(i, o, p, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim8, Dim3, i, o, p, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim0, Dim7, Dim8, o, l, i, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim7, Dim8, i, j, k, l, m, n, o, l, i, p, q, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, l, i, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim0, Dim7, Dim8, o, l, i, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim3, Dim7, Dim8, o, i, l, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, Dim8, i, j, k, l, m, n, o, i, l, p, q, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, i, l, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim3, Dim7, Dim8, o, i, l, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim0, Dim8, o, l, p, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim0, Dim8, i, j, k, l, m, n, o, l, p, i, q, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, l, p, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim0, Dim8, o, l, p, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim3, Dim8, o, i, p, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim3, Dim8, i, j, k, l, m, n, o, i, p, l, q, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, i, p, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim3, Dim8, o, i, p, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim8, Dim0, o, l, p, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim8, Dim0, i, j, k, l, m, n, o, l, p, q, i, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, l, p, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim8, Dim0, o, l, p, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim8, Dim3, o, i, p, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim8, Dim3, i, j, k, l, m, n, o, i, p, q, l, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, i, p, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim8, Dim3, o, i, p, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim0, Dim8, o, p, l, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim0, Dim8, i, j, k, l, m, n, o, p, l, i, q, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, l, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim0, Dim8, o, p, l, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim3, Dim8, o, p, i, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim3, Dim8, i, j, k, l, m, n, o, p, i, l, q, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, i, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim3, Dim8, o, p, i, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim8, Dim0, o, p, l, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim8, Dim0, i, j, k, l, m, n, o, p, l, q, i, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, l, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim8, Dim0, o, p, l, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim8, Dim3, o, p, i, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim8, Dim3, i, j, k, l, m, n, o, p, i, q, l, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, i, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim8, Dim3, o, p, i, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim3, Dim0, o, p, q, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim3, Dim0, i, j, k, l, m, n, o, p, q, l, i, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim3, Dim0, o, p, q, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim0, Dim3, o, p, q, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim0, Dim3, i, j, k, l, m, n, o, p, q, i, l, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, Dim8, j, k, m, n, i, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, j, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim0, Dim3, o, p, q, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim0, Dim6, Dim7, Dim8, m, i, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim6, Dim7, Dim8, i, j, k, l, m, n, m, i, o, p, q, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, i, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim0, Dim6, Dim7, Dim8, m, i, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim4, Dim6, Dim7, Dim8, i, m, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, i, j, k, l, m, n, i, m, o, p, q, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(i, m, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim4, Dim6, Dim7, Dim8, i, m, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim0, Dim7, Dim8, m, o, i, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim0, Dim7, Dim8, i, j, k, l, m, n, m, o, i, p, q, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, o, i, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim0, Dim7, Dim8, m, o, i, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim4, Dim7, Dim8, i, o, m, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim4, Dim7, Dim8, i, j, k, l, m, n, i, o, m, p, q, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(i, o, m, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim4, Dim7, Dim8, i, o, m, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim0, Dim8, m, o, p, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim0, Dim8, i, j, k, l, m, n, m, o, p, i, q, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, o, p, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim0, Dim8, m, o, p, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim4, Dim8, i, o, p, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim4, Dim8, i, j, k, l, m, n, i, o, p, m, q, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(i, o, p, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim4, Dim8, i, o, p, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim8, Dim0, m, o, p, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim8, Dim0, i, j, k, l, m, n, m, o, p, q, i, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, o, p, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim8, Dim0, m, o, p, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim8, Dim4, i, o, p, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, Dim4, i, j, k, l, m, n, i, o, p, q, m, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(i, o, p, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim8, Dim4, i, o, p, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim0, Dim7, Dim8, o, m, i, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim7, Dim8, i, j, k, l, m, n, o, m, i, p, q, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, m, i, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim0, Dim7, Dim8, o, m, i, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim4, Dim7, Dim8, o, i, m, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim7, Dim8, i, j, k, l, m, n, o, i, m, p, q, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, i, m, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim4, Dim7, Dim8, o, i, m, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim0, Dim8, o, m, p, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim0, Dim8, i, j, k, l, m, n, o, m, p, i, q, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, m, p, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim0, Dim8, o, m, p, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim4, Dim8, o, i, p, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim4, Dim8, i, j, k, l, m, n, o, i, p, m, q, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, i, p, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim4, Dim8, o, i, p, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim8, Dim0, o, m, p, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim8, Dim0, i, j, k, l, m, n, o, m, p, q, i, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, m, p, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim8, Dim0, o, m, p, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim8, Dim4, o, i, p, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim8, Dim4, i, j, k, l, m, n, o, i, p, q, m, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, i, p, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim8, Dim4, o, i, p, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim0, Dim8, o, p, m, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim0, Dim8, i, j, k, l, m, n, o, p, m, i, q, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, m, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim0, Dim8, o, p, m, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim4, Dim8, o, p, i, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim4, Dim8, i, j, k, l, m, n, o, p, i, m, q, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, i, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim4, Dim8, o, p, i, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim8, Dim0, o, p, m, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim8, Dim0, i, j, k, l, m, n, o, p, m, q, i, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, m, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim8, Dim0, o, p, m, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim8, Dim4, o, p, i, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim8, Dim4, i, j, k, l, m, n, o, p, i, q, m, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, i, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim8, Dim4, o, p, i, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim4, Dim0, o, p, q, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim4, Dim0, i, j, k, l, m, n, o, p, q, m, i, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim4, Dim0, o, p, q, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim0, Dim4, o, p, q, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim0, Dim4, i, j, k, l, m, n, o, p, q, i, m, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, Dim8, j, k, l, n, i, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, j, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim0, Dim4, o, p, q, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim7, Dim8, n, i, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim6, Dim7, Dim8, i, j, k, l, m, n, n, i, o, p, q, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, i, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim7, Dim8, n, i, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim7, Dim8, i, n, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, i, n, o, p, q, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(i, n, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim7, Dim8, i, n, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim7, Dim8, n, o, i, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim0, Dim7, Dim8, i, j, k, l, m, n, n, o, i, p, q, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, o, i, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim7, Dim8, n, o, i, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim5, Dim7, Dim8, i, o, n, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim5, Dim7, Dim8, i, j, k, l, m, n, i, o, n, p, q, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(i, o, n, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim5, Dim7, Dim8, i, o, n, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim0, Dim8, n, o, p, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim0, Dim8, i, j, k, l, m, n, n, o, p, i, q, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, o, p, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim0, Dim8, n, o, p, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim5, Dim8, i, o, p, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim5, Dim8, i, j, k, l, m, n, i, o, p, n, q, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(i, o, p, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim5, Dim8, i, o, p, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim0, n, o, p, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim8, Dim0, i, j, k, l, m, n, n, o, p, q, i, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, o, p, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim0, n, o, p, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim8, Dim5, i, o, p, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim8, Dim5, i, j, k, l, m, n, i, o, p, q, n, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(i, o, p, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim0, Dim6, Dim7, Dim8, Dim5, i, o, p, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, i, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim0, Dim7, Dim8, o, n, i, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim7, Dim8, i, j, k, l, m, n, o, n, i, p, q, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, n, i, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim0, Dim7, Dim8, o, n, i, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim5, Dim7, Dim8, o, i, n, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim7, Dim8, i, j, k, l, m, n, o, i, n, p, q, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, i, n, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim5, Dim7, Dim8, o, i, n, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim0, Dim8, o, n, p, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim0, Dim8, i, j, k, l, m, n, o, n, p, i, q, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, n, p, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim0, Dim8, o, n, p, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim5, Dim8, o, i, p, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim5, Dim8, i, j, k, l, m, n, o, i, p, n, q, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, i, p, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim5, Dim8, o, i, p, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim8, Dim0, o, n, p, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim8, Dim0, i, j, k, l, m, n, o, n, p, q, i, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, n, p, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim8, Dim0, o, n, p, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim8, Dim5, o, i, p, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim8, Dim5, i, j, k, l, m, n, o, i, p, q, n, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, i, p, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim0, Dim7, Dim8, Dim5, o, i, p, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, i, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim0, Dim8, o, p, n, i, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim0, Dim8, i, j, k, l, m, n, o, p, n, i, q, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, n, i, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim0, Dim8, o, p, n, i, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim5, Dim8, o, p, i, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim5, Dim8, i, j, k, l, m, n, o, p, i, n, q, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, i, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim5, Dim8, o, p, i, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, q, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim8, Dim0, o, p, n, q, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim8, Dim0, i, j, k, l, m, n, o, p, n, q, i, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, n, q, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim8, Dim0, o, p, n, q, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim8, Dim5, o, p, i, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim8, Dim5, i, j, k, l, m, n, o, p, i, q, n, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, i, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim0, Dim8, Dim5, o, p, i, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim5, Dim0, o, p, q, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim5, Dim0, i, j, k, l, m, n, o, p, q, n, i, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim5, Dim0, o, p, q, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim0, Dim5, o, p, q, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim0, Dim5, i, j, k, l, m, n, o, p, q, i, n, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim8, j, k, l, m, i, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, j, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim0, Dim5, o, p, q, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim1, Dim6, Dim7, Dim8, k, j, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, Dim7, Dim8, i, j, k, l, m, n, k, j, o, p, q, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, j, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim1, Dim6, Dim7, Dim8, k, j, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim2, Dim6, Dim7, Dim8, j, k, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, j, k, l, m, n, j, k, o, p, q, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, k, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim2, Dim6, Dim7, Dim8, j, k, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim1, Dim7, Dim8, k, o, j, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, Dim7, Dim8, i, j, k, l, m, n, k, o, j, p, q, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, o, j, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim1, Dim7, Dim8, k, o, j, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim2, Dim7, Dim8, j, o, k, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, Dim7, Dim8, i, j, k, l, m, n, j, o, k, p, q, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, o, k, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim2, Dim7, Dim8, j, o, k, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim1, Dim8, k, o, p, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim1, Dim8, i, j, k, l, m, n, k, o, p, j, q, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, o, p, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim1, Dim8, k, o, p, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim2, Dim8, j, o, p, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim2, Dim8, i, j, k, l, m, n, j, o, p, k, q, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, o, p, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim2, Dim8, j, o, p, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim8, Dim1, k, o, p, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, Dim1, i, j, k, l, m, n, k, o, p, q, j, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, o, p, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim8, Dim1, k, o, p, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim8, Dim2, j, o, p, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, Dim2, i, j, k, l, m, n, j, o, p, q, k, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, o, p, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim8, Dim2, j, o, p, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim1, Dim7, Dim8, o, k, j, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim7, Dim8, i, j, k, l, m, n, o, k, j, p, q, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, k, j, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim1, Dim7, Dim8, o, k, j, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim2, Dim7, Dim8, o, j, k, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, Dim8, i, j, k, l, m, n, o, j, k, p, q, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, j, k, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim2, Dim7, Dim8, o, j, k, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim1, Dim8, o, k, p, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim1, Dim8, i, j, k, l, m, n, o, k, p, j, q, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, k, p, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim1, Dim8, o, k, p, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim2, Dim8, o, j, p, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim2, Dim8, i, j, k, l, m, n, o, j, p, k, q, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, j, p, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim2, Dim8, o, j, p, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim8, Dim1, o, k, p, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim8, Dim1, i, j, k, l, m, n, o, k, p, q, j, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, k, p, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim8, Dim1, o, k, p, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim8, Dim2, o, j, p, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim8, Dim2, i, j, k, l, m, n, o, j, p, q, k, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, j, p, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim8, Dim2, o, j, p, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim1, Dim8, o, p, k, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim1, Dim8, i, j, k, l, m, n, o, p, k, j, q, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, k, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim1, Dim8, o, p, k, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim2, Dim8, o, p, j, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim2, Dim8, i, j, k, l, m, n, o, p, j, k, q, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, j, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim2, Dim8, o, p, j, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim8, Dim1, o, p, k, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim8, Dim1, i, j, k, l, m, n, o, p, k, q, j, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, k, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim8, Dim1, o, p, k, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim8, Dim2, o, p, j, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim8, Dim2, i, j, k, l, m, n, o, p, j, q, k, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, j, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim8, Dim2, o, p, j, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim2, Dim1, o, p, q, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim2, Dim1, i, j, k, l, m, n, o, p, q, k, j, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim2, Dim1, o, p, q, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim1, Dim2, o, p, q, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim1, Dim2, i, j, k, l, m, n, o, p, q, j, k, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, Dim8, i, l, m, n, j, k, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, l, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim1, Dim2, o, p, q, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim1, Dim6, Dim7, Dim8, l, j, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim6, Dim7, Dim8, i, j, k, l, m, n, l, j, o, p, q, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, j, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim1, Dim6, Dim7, Dim8, l, j, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim3, Dim6, Dim7, Dim8, j, l, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, j, k, l, m, n, j, l, o, p, q, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, l, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim3, Dim6, Dim7, Dim8, j, l, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim1, Dim7, Dim8, l, o, j, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim1, Dim7, Dim8, i, j, k, l, m, n, l, o, j, p, q, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, o, j, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim1, Dim7, Dim8, l, o, j, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim3, Dim7, Dim8, j, o, l, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim3, Dim7, Dim8, i, j, k, l, m, n, j, o, l, p, q, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, o, l, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim3, Dim7, Dim8, j, o, l, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim1, Dim8, l, o, p, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim1, Dim8, i, j, k, l, m, n, l, o, p, j, q, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, o, p, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim1, Dim8, l, o, p, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim3, Dim8, j, o, p, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim3, Dim8, i, j, k, l, m, n, j, o, p, l, q, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, o, p, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim3, Dim8, j, o, p, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim8, Dim1, l, o, p, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, Dim1, i, j, k, l, m, n, l, o, p, q, j, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, o, p, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim8, Dim1, l, o, p, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim8, Dim3, j, o, p, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, Dim3, i, j, k, l, m, n, j, o, p, q, l, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(j, o, p, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim8, Dim3, j, o, p, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim1, Dim7, Dim8, o, l, j, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim7, Dim8, i, j, k, l, m, n, o, l, j, p, q, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, l, j, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim1, Dim7, Dim8, o, l, j, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim3, Dim7, Dim8, o, j, l, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, Dim8, i, j, k, l, m, n, o, j, l, p, q, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, j, l, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim3, Dim7, Dim8, o, j, l, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim1, Dim8, o, l, p, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim1, Dim8, i, j, k, l, m, n, o, l, p, j, q, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, l, p, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim1, Dim8, o, l, p, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim3, Dim8, o, j, p, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim3, Dim8, i, j, k, l, m, n, o, j, p, l, q, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, j, p, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim3, Dim8, o, j, p, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim8, Dim1, o, l, p, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim8, Dim1, i, j, k, l, m, n, o, l, p, q, j, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, l, p, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim8, Dim1, o, l, p, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim8, Dim3, o, j, p, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim8, Dim3, i, j, k, l, m, n, o, j, p, q, l, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, j, p, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim8, Dim3, o, j, p, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim1, Dim8, o, p, l, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim1, Dim8, i, j, k, l, m, n, o, p, l, j, q, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, l, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim1, Dim8, o, p, l, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim3, Dim8, o, p, j, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim3, Dim8, i, j, k, l, m, n, o, p, j, l, q, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, j, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim3, Dim8, o, p, j, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim8, Dim1, o, p, l, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim8, Dim1, i, j, k, l, m, n, o, p, l, q, j, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, l, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim8, Dim1, o, p, l, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim8, Dim3, o, p, j, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim8, Dim3, i, j, k, l, m, n, o, p, j, q, l, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, j, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim8, Dim3, o, p, j, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim3, Dim1, o, p, q, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim3, Dim1, i, j, k, l, m, n, o, p, q, l, j, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim3, Dim1, o, p, q, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim1, Dim3, o, p, q, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim1, Dim3, i, j, k, l, m, n, o, p, q, j, l, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, Dim8, i, k, m, n, j, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim8, i, k, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim1, Dim3, o, p, q, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim1, Dim6, Dim7, Dim8, m, j, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim6, Dim7, Dim8, i, j, k, l, m, n, m, j, o, p, q, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, j, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim1, Dim6, Dim7, Dim8, m, j, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim4, Dim6, Dim7, Dim8, j, m, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, j, k, l, m, n, j, m, o, p, q, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(j, m, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim4, Dim6, Dim7, Dim8, j, m, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim1, Dim7, Dim8, m, o, j, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim1, Dim7, Dim8, i, j, k, l, m, n, m, o, j, p, q, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, o, j, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim1, Dim7, Dim8, m, o, j, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim4, Dim7, Dim8, j, o, m, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim4, Dim7, Dim8, i, j, k, l, m, n, j, o, m, p, q, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(j, o, m, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim4, Dim7, Dim8, j, o, m, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim1, Dim8, m, o, p, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim1, Dim8, i, j, k, l, m, n, m, o, p, j, q, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, o, p, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim1, Dim8, m, o, p, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim4, Dim8, j, o, p, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim4, Dim8, i, j, k, l, m, n, j, o, p, m, q, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(j, o, p, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim4, Dim8, j, o, p, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim8, Dim1, m, o, p, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim8, Dim1, i, j, k, l, m, n, m, o, p, q, j, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, o, p, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim8, Dim1, m, o, p, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim8, Dim4, j, o, p, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, Dim4, i, j, k, l, m, n, j, o, p, q, m, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(j, o, p, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim8, Dim4, j, o, p, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim1, Dim7, Dim8, o, m, j, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim7, Dim8, i, j, k, l, m, n, o, m, j, p, q, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, m, j, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim1, Dim7, Dim8, o, m, j, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim4, Dim7, Dim8, o, j, m, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim7, Dim8, i, j, k, l, m, n, o, j, m, p, q, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, j, m, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim4, Dim7, Dim8, o, j, m, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim1, Dim8, o, m, p, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim1, Dim8, i, j, k, l, m, n, o, m, p, j, q, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, m, p, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim1, Dim8, o, m, p, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim4, Dim8, o, j, p, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim4, Dim8, i, j, k, l, m, n, o, j, p, m, q, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, j, p, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim4, Dim8, o, j, p, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim8, Dim1, o, m, p, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim8, Dim1, i, j, k, l, m, n, o, m, p, q, j, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, m, p, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim8, Dim1, o, m, p, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim8, Dim4, o, j, p, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim8, Dim4, i, j, k, l, m, n, o, j, p, q, m, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, j, p, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim8, Dim4, o, j, p, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim1, Dim8, o, p, m, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim1, Dim8, i, j, k, l, m, n, o, p, m, j, q, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, m, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim1, Dim8, o, p, m, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim4, Dim8, o, p, j, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim4, Dim8, i, j, k, l, m, n, o, p, j, m, q, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, j, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim4, Dim8, o, p, j, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim8, Dim1, o, p, m, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim8, Dim1, i, j, k, l, m, n, o, p, m, q, j, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, m, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim8, Dim1, o, p, m, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim8, Dim4, o, p, j, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim8, Dim4, i, j, k, l, m, n, o, p, j, q, m, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, j, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim8, Dim4, o, p, j, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim4, Dim1, o, p, q, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim4, Dim1, i, j, k, l, m, n, o, p, q, m, j, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim4, Dim1, o, p, q, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim1, Dim4, o, p, q, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim1, Dim4, i, j, k, l, m, n, o, p, q, j, m, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, Dim8, i, k, l, n, j, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim8, i, k, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim1, Dim4, o, p, q, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim7, Dim8, n, j, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim6, Dim7, Dim8, i, j, k, l, m, n, n, j, o, p, q, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, j, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim7, Dim8, n, j, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim7, Dim8, j, n, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, j, n, o, p, q, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(j, n, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim7, Dim8, j, n, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim7, Dim8, n, o, j, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim1, Dim7, Dim8, i, j, k, l, m, n, n, o, j, p, q, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, o, j, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim7, Dim8, n, o, j, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim5, Dim7, Dim8, j, o, n, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim5, Dim7, Dim8, i, j, k, l, m, n, j, o, n, p, q, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(j, o, n, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim5, Dim7, Dim8, j, o, n, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim1, Dim8, n, o, p, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim1, Dim8, i, j, k, l, m, n, n, o, p, j, q, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, o, p, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim1, Dim8, n, o, p, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim5, Dim8, j, o, p, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim5, Dim8, i, j, k, l, m, n, j, o, p, n, q, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(j, o, p, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim5, Dim8, j, o, p, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim1, n, o, p, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim8, Dim1, i, j, k, l, m, n, n, o, p, q, j, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, o, p, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim1, n, o, p, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim8, Dim5, j, o, p, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim8, Dim5, i, j, k, l, m, n, j, o, p, q, n, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(j, o, p, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim1, Dim6, Dim7, Dim8, Dim5, j, o, p, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, j, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim1, Dim7, Dim8, o, n, j, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim7, Dim8, i, j, k, l, m, n, o, n, j, p, q, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, n, j, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim1, Dim7, Dim8, o, n, j, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim5, Dim7, Dim8, o, j, n, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim7, Dim8, i, j, k, l, m, n, o, j, n, p, q, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, j, n, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim5, Dim7, Dim8, o, j, n, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim1, Dim8, o, n, p, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim1, Dim8, i, j, k, l, m, n, o, n, p, j, q, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, n, p, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim1, Dim8, o, n, p, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim5, Dim8, o, j, p, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim5, Dim8, i, j, k, l, m, n, o, j, p, n, q, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, j, p, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim5, Dim8, o, j, p, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim8, Dim1, o, n, p, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim8, Dim1, i, j, k, l, m, n, o, n, p, q, j, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, n, p, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim8, Dim1, o, n, p, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim8, Dim5, o, j, p, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim8, Dim5, i, j, k, l, m, n, o, j, p, q, n, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, j, p, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim1, Dim7, Dim8, Dim5, o, j, p, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, j, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim1, Dim8, o, p, n, j, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim1, Dim8, i, j, k, l, m, n, o, p, n, j, q, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, n, j, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim1, Dim8, o, p, n, j, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim5, Dim8, o, p, j, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim5, Dim8, i, j, k, l, m, n, o, p, j, n, q, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, j, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim5, Dim8, o, p, j, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, q, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim8, Dim1, o, p, n, q, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim8, Dim1, i, j, k, l, m, n, o, p, n, q, j, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, n, q, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim8, Dim1, o, p, n, q, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim8, Dim5, o, p, j, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim8, Dim5, i, j, k, l, m, n, o, p, j, q, n, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, j, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim1, Dim8, Dim5, o, p, j, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim5, Dim1, o, p, q, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim5, Dim1, i, j, k, l, m, n, o, p, q, n, j, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim5, Dim1, o, p, q, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim1, Dim5, o, p, q, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim1, Dim5, i, j, k, l, m, n, o, p, q, j, n, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim8, i, k, l, m, j, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim8, i, k, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim1, Dim5, o, p, q, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim2, Dim6, Dim7, Dim8, l, k, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim6, Dim7, Dim8, i, j, k, l, m, n, l, k, o, p, q, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, k, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim2, Dim6, Dim7, Dim8, l, k, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim3, Dim6, Dim7, Dim8, k, l, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, m, n, k, l, o, p, q, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, l, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim3, Dim6, Dim7, Dim8, k, l, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim2, Dim7, Dim8, l, o, k, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim2, Dim7, Dim8, i, j, k, l, m, n, l, o, k, p, q, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, o, k, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim2, Dim7, Dim8, l, o, k, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim3, Dim7, Dim8, k, o, l, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim3, Dim7, Dim8, i, j, k, l, m, n, k, o, l, p, q, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, o, l, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim3, Dim7, Dim8, k, o, l, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim2, Dim8, l, o, p, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim2, Dim8, i, j, k, l, m, n, l, o, p, k, q, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, o, p, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim2, Dim8, l, o, p, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim3, Dim8, k, o, p, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim3, Dim8, i, j, k, l, m, n, k, o, p, l, q, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, o, p, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim3, Dim8, k, o, p, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim8, Dim2, l, o, p, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, Dim2, i, j, k, l, m, n, l, o, p, q, k, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(l, o, p, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim8, Dim2, l, o, p, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim8, Dim3, k, o, p, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, Dim3, i, j, k, l, m, n, k, o, p, q, l, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(k, o, p, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim8, Dim3, k, o, p, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim2, Dim7, Dim8, o, l, k, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim7, Dim8, i, j, k, l, m, n, o, l, k, p, q, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, l, k, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim2, Dim7, Dim8, o, l, k, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim3, Dim7, Dim8, o, k, l, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, Dim8, i, j, k, l, m, n, o, k, l, p, q, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, k, l, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim3, Dim7, Dim8, o, k, l, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim2, Dim8, o, l, p, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim2, Dim8, i, j, k, l, m, n, o, l, p, k, q, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, l, p, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim2, Dim8, o, l, p, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim3, Dim8, o, k, p, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim3, Dim8, i, j, k, l, m, n, o, k, p, l, q, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, k, p, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim3, Dim8, o, k, p, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim8, Dim2, o, l, p, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim8, Dim2, i, j, k, l, m, n, o, l, p, q, k, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, l, p, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim8, Dim2, o, l, p, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim8, Dim3, o, k, p, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim8, Dim3, i, j, k, l, m, n, o, k, p, q, l, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, k, p, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim8, Dim3, o, k, p, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim2, Dim8, o, p, l, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim2, Dim8, i, j, k, l, m, n, o, p, l, k, q, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, l, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim2, Dim8, o, p, l, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim3, Dim8, o, p, k, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim3, Dim8, i, j, k, l, m, n, o, p, k, l, q, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, k, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim3, Dim8, o, p, k, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim8, Dim2, o, p, l, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim8, Dim2, i, j, k, l, m, n, o, p, l, q, k, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, l, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim8, Dim2, o, p, l, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim8, Dim3, o, p, k, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim8, Dim3, i, j, k, l, m, n, o, p, k, q, l, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, k, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim8, Dim3, o, p, k, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim3, Dim2, o, p, q, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim3, Dim2, i, j, k, l, m, n, o, p, q, l, k, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim3, Dim2, o, p, q, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim2, Dim3, o, p, q, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim2, Dim3, i, j, k, l, m, n, o, p, q, k, l, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, m, n, k, l, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, m, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim2, Dim3, o, p, q, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim2, Dim6, Dim7, Dim8, m, k, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim6, Dim7, Dim8, i, j, k, l, m, n, m, k, o, p, q, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, k, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim2, Dim6, Dim7, Dim8, m, k, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim4, Dim6, Dim7, Dim8, k, m, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, l, m, n, k, m, o, p, q, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(k, m, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim4, Dim6, Dim7, Dim8, k, m, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim2, Dim7, Dim8, m, o, k, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim2, Dim7, Dim8, i, j, k, l, m, n, m, o, k, p, q, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, o, k, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim2, Dim7, Dim8, m, o, k, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim4, Dim7, Dim8, k, o, m, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim4, Dim7, Dim8, i, j, k, l, m, n, k, o, m, p, q, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(k, o, m, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim4, Dim7, Dim8, k, o, m, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim2, Dim8, m, o, p, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim2, Dim8, i, j, k, l, m, n, m, o, p, k, q, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, o, p, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim2, Dim8, m, o, p, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim4, Dim8, k, o, p, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim4, Dim8, i, j, k, l, m, n, k, o, p, m, q, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(k, o, p, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim4, Dim8, k, o, p, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim8, Dim2, m, o, p, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim8, Dim2, i, j, k, l, m, n, m, o, p, q, k, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(m, o, p, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim8, Dim2, m, o, p, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim8, Dim4, k, o, p, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, Dim4, i, j, k, l, m, n, k, o, p, q, m, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(k, o, p, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim8, Dim4, k, o, p, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim2, Dim7, Dim8, o, m, k, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim7, Dim8, i, j, k, l, m, n, o, m, k, p, q, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, m, k, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim2, Dim7, Dim8, o, m, k, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim4, Dim7, Dim8, o, k, m, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim7, Dim8, i, j, k, l, m, n, o, k, m, p, q, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, k, m, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim4, Dim7, Dim8, o, k, m, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim2, Dim8, o, m, p, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim2, Dim8, i, j, k, l, m, n, o, m, p, k, q, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, m, p, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim2, Dim8, o, m, p, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim4, Dim8, o, k, p, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim4, Dim8, i, j, k, l, m, n, o, k, p, m, q, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, k, p, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim4, Dim8, o, k, p, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim8, Dim2, o, m, p, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim8, Dim2, i, j, k, l, m, n, o, m, p, q, k, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, m, p, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim8, Dim2, o, m, p, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim8, Dim4, o, k, p, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim8, Dim4, i, j, k, l, m, n, o, k, p, q, m, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, k, p, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim8, Dim4, o, k, p, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim2, Dim8, o, p, m, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim2, Dim8, i, j, k, l, m, n, o, p, m, k, q, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, m, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim2, Dim8, o, p, m, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim4, Dim8, o, p, k, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim4, Dim8, i, j, k, l, m, n, o, p, k, m, q, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, k, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim4, Dim8, o, p, k, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim8, Dim2, o, p, m, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim8, Dim2, i, j, k, l, m, n, o, p, m, q, k, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, m, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim8, Dim2, o, p, m, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim8, Dim4, o, p, k, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim8, Dim4, i, j, k, l, m, n, o, p, k, q, m, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, k, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim8, Dim4, o, p, k, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim4, Dim2, o, p, q, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim4, Dim2, i, j, k, l, m, n, o, p, q, m, k, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim4, Dim2, o, p, q, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim2, Dim4, o, p, q, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim2, Dim4, i, j, k, l, m, n, o, p, q, k, m, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, l, n, k, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, l, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim2, Dim4, o, p, q, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim7, Dim8, n, k, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim6, Dim7, Dim8, i, j, k, l, m, n, n, k, o, p, q, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, k, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim7, Dim8, n, k, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim7, Dim8, k, n, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, k, n, o, p, q, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(k, n, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim7, Dim8, k, n, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim7, Dim8, n, o, k, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim2, Dim7, Dim8, i, j, k, l, m, n, n, o, k, p, q, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, o, k, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim7, Dim8, n, o, k, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim5, Dim7, Dim8, k, o, n, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim5, Dim7, Dim8, i, j, k, l, m, n, k, o, n, p, q, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(k, o, n, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim5, Dim7, Dim8, k, o, n, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim2, Dim8, n, o, p, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim2, Dim8, i, j, k, l, m, n, n, o, p, k, q, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, o, p, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim2, Dim8, n, o, p, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim5, Dim8, k, o, p, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim5, Dim8, i, j, k, l, m, n, k, o, p, n, q, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(k, o, p, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim5, Dim8, k, o, p, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim2, n, o, p, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim8, Dim2, i, j, k, l, m, n, n, o, p, q, k, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(n, o, p, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim2, n, o, p, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim8, Dim5, k, o, p, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim8, Dim5, i, j, k, l, m, n, k, o, p, q, n, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(k, o, p, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim2, Dim6, Dim7, Dim8, Dim5, k, o, p, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, k, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim2, Dim7, Dim8, o, n, k, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim7, Dim8, i, j, k, l, m, n, o, n, k, p, q, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, n, k, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim2, Dim7, Dim8, o, n, k, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim5, Dim7, Dim8, o, k, n, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim7, Dim8, i, j, k, l, m, n, o, k, n, p, q, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, k, n, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim5, Dim7, Dim8, o, k, n, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim2, Dim8, o, n, p, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim2, Dim8, i, j, k, l, m, n, o, n, p, k, q, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, n, p, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim2, Dim8, o, n, p, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim5, Dim8, o, k, p, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim5, Dim8, i, j, k, l, m, n, o, k, p, n, q, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, k, p, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim5, Dim8, o, k, p, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim8, Dim2, o, n, p, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim8, Dim2, i, j, k, l, m, n, o, n, p, q, k, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, n, p, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim8, Dim2, o, n, p, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim8, Dim5, o, k, p, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim8, Dim5, i, j, k, l, m, n, o, k, p, q, n, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, k, p, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim2, Dim7, Dim8, Dim5, o, k, p, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, k, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim2, Dim8, o, p, n, k, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim2, Dim8, i, j, k, l, m, n, o, p, n, k, q, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, n, k, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim2, Dim8, o, p, n, k, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim5, Dim8, o, p, k, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim5, Dim8, i, j, k, l, m, n, o, p, k, n, q, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, k, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim5, Dim8, o, p, k, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, q, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim8, Dim2, o, p, n, q, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim8, Dim2, i, j, k, l, m, n, o, p, n, q, k, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, n, q, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim8, Dim2, o, p, n, q, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim8, Dim5, o, p, k, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim8, Dim5, i, j, k, l, m, n, o, p, k, q, n, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, k, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim2, Dim8, Dim5, o, p, k, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim5, Dim2, o, p, q, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim5, Dim2, i, j, k, l, m, n, o, p, q, n, k, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim5, Dim2, o, p, q, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim2, Dim5, o, p, q, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim2, Dim5, i, j, k, l, m, n, o, p, q, k, n, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, l, m, k, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, l, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim2, Dim5, o, p, q, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim3, Dim6, Dim7, Dim8, m, l, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim6, Dim7, Dim8, i, j, k, l, m, n, m, l, o, p, q, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(m, l, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim3, Dim6, Dim7, Dim8, m, l, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim4, Dim6, Dim7, Dim8, l, m, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, l, m, n, l, m, o, p, q, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(l, m, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim4, Dim6, Dim7, Dim8, l, m, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim3, Dim7, Dim8, m, o, l, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim3, Dim7, Dim8, i, j, k, l, m, n, m, o, l, p, q, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(m, o, l, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim3, Dim7, Dim8, m, o, l, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim4, Dim7, Dim8, l, o, m, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim4, Dim7, Dim8, i, j, k, l, m, n, l, o, m, p, q, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(l, o, m, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim4, Dim7, Dim8, l, o, m, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim3, Dim8, m, o, p, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim3, Dim8, i, j, k, l, m, n, m, o, p, l, q, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(m, o, p, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim3, Dim8, m, o, p, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim4, Dim8, l, o, p, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim4, Dim8, i, j, k, l, m, n, l, o, p, m, q, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(l, o, p, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim4, Dim8, l, o, p, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim8, Dim3, m, o, p, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim8, Dim3, i, j, k, l, m, n, m, o, p, q, l, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(m, o, p, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim8, Dim3, m, o, p, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim8, Dim4, l, o, p, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, Dim4, i, j, k, l, m, n, l, o, p, q, m, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(l, o, p, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim8, Dim4, l, o, p, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim3, Dim7, Dim8, o, m, l, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim7, Dim8, i, j, k, l, m, n, o, m, l, p, q, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(o, m, l, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim3, Dim7, Dim8, o, m, l, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim4, Dim7, Dim8, o, l, m, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim7, Dim8, i, j, k, l, m, n, o, l, m, p, q, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(o, l, m, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim4, Dim7, Dim8, o, l, m, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim3, Dim8, o, m, p, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim3, Dim8, i, j, k, l, m, n, o, m, p, l, q, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(o, m, p, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim3, Dim8, o, m, p, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim4, Dim8, o, l, p, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim4, Dim8, i, j, k, l, m, n, o, l, p, m, q, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(o, l, p, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim4, Dim8, o, l, p, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim8, Dim3, o, m, p, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim8, Dim3, i, j, k, l, m, n, o, m, p, q, l, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(o, m, p, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim8, Dim3, o, m, p, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim8, Dim4, o, l, p, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim8, Dim4, i, j, k, l, m, n, o, l, p, q, m, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(o, l, p, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim8, Dim4, o, l, p, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim3, Dim8, o, p, m, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim3, Dim8, i, j, k, l, m, n, o, p, m, l, q, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, m, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim3, Dim8, o, p, m, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim4, Dim8, o, p, l, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim4, Dim8, i, j, k, l, m, n, o, p, l, m, q, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, l, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim4, Dim8, o, p, l, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim8, Dim3, o, p, m, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim8, Dim3, i, j, k, l, m, n, o, p, m, q, l, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, m, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim8, Dim3, o, p, m, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim8, Dim4, o, p, l, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim8, Dim4, i, j, k, l, m, n, o, p, l, q, m, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, l, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim8, Dim4, o, p, l, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim4, Dim3, o, p, q, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim4, Dim3, i, j, k, l, m, n, o, p, q, m, l, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim4, Dim3, o, p, q, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim3, Dim4, o, p, q, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim3, Dim4, i, j, k, l, m, n, o, p, q, l, m, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, Dim8, i, j, k, n, l, m, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim8, i, j, k, n, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim3, Dim4, o, p, q, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim7, Dim8, n, l, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim6, Dim7, Dim8, i, j, k, l, m, n, n, l, o, p, q, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(n, l, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim7, Dim8, n, l, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim7, Dim8, l, n, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, l, n, o, p, q, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(l, n, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim7, Dim8, l, n, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim7, Dim8, n, o, l, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim3, Dim7, Dim8, i, j, k, l, m, n, n, o, l, p, q, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(n, o, l, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim7, Dim8, n, o, l, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim5, Dim7, Dim8, l, o, n, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim5, Dim7, Dim8, i, j, k, l, m, n, l, o, n, p, q, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(l, o, n, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim5, Dim7, Dim8, l, o, n, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim3, Dim8, n, o, p, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim3, Dim8, i, j, k, l, m, n, n, o, p, l, q, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(n, o, p, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim3, Dim8, n, o, p, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim5, Dim8, l, o, p, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim5, Dim8, i, j, k, l, m, n, l, o, p, n, q, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(l, o, p, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim5, Dim8, l, o, p, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim3, n, o, p, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim8, Dim3, i, j, k, l, m, n, n, o, p, q, l, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(n, o, p, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim3, n, o, p, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim8, Dim5, l, o, p, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim8, Dim5, i, j, k, l, m, n, l, o, p, q, n, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(l, o, p, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim3, Dim6, Dim7, Dim8, Dim5, l, o, p, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, l, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim3, Dim7, Dim8, o, n, l, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim7, Dim8, i, j, k, l, m, n, o, n, l, p, q, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(o, n, l, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim3, Dim7, Dim8, o, n, l, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim5, Dim7, Dim8, o, l, n, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim7, Dim8, i, j, k, l, m, n, o, l, n, p, q, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(o, l, n, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim5, Dim7, Dim8, o, l, n, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim3, Dim8, o, n, p, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim3, Dim8, i, j, k, l, m, n, o, n, p, l, q, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(o, n, p, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim3, Dim8, o, n, p, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim5, Dim8, o, l, p, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim5, Dim8, i, j, k, l, m, n, o, l, p, n, q, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(o, l, p, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim5, Dim8, o, l, p, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim8, Dim3, o, n, p, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim8, Dim3, i, j, k, l, m, n, o, n, p, q, l, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(o, n, p, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim8, Dim3, o, n, p, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim8, Dim5, o, l, p, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim8, Dim5, i, j, k, l, m, n, o, l, p, q, n, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(o, l, p, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim3, Dim7, Dim8, Dim5, o, l, p, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, l, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim3, Dim8, o, p, n, l, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim3, Dim8, i, j, k, l, m, n, o, p, n, l, q, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, n, l, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim3, Dim8, o, p, n, l, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim5, Dim8, o, p, l, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim5, Dim8, i, j, k, l, m, n, o, p, l, n, q, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, l, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim5, Dim8, o, p, l, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, q, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim8, Dim3, o, p, n, q, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim8, Dim3, i, j, k, l, m, n, o, p, n, q, l, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, n, q, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim8, Dim3, o, p, n, q, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim8, Dim5, o, p, l, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim8, Dim5, i, j, k, l, m, n, o, p, l, q, n, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, l, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim3, Dim8, Dim5, o, p, l, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim5, Dim3, o, p, q, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim5, Dim3, i, j, k, l, m, n, o, p, q, n, l, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim5, Dim3, o, p, q, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim3, Dim5, o, p, q, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim3, Dim5, i, j, k, l, m, n, o, p, q, l, n, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, Dim8, i, j, k, m, l, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim8, i, j, k, m, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim3, Dim5, o, p, q, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim7, Dim8, n, m, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim6, Dim7, Dim8, i, j, k, l, m, n, n, m, o, p, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(n, m, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim7, Dim8, n, m, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, o, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim7, Dim8, m, n, o, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, m, n, o, p, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(m, n, o, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim7, Dim8, m, n, o, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim7, Dim8, n, o, m, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim4, Dim7, Dim8, i, j, k, l, m, n, n, o, m, p, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(n, o, m, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim7, Dim8, n, o, m, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim5, Dim7, Dim8, m, o, n, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim5, Dim7, Dim8, i, j, k, l, m, n, m, o, n, p, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(m, o, n, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim5, Dim7, Dim8, m, o, n, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim4, Dim8, n, o, p, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim4, Dim8, i, j, k, l, m, n, n, o, p, m, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(n, o, p, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim4, Dim8, n, o, p, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim5, Dim8, m, o, p, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim5, Dim8, i, j, k, l, m, n, m, o, p, n, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(m, o, p, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim5, Dim8, m, o, p, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim4, n, o, p, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim8, Dim4, i, j, k, l, m, n, n, o, p, q, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(n, o, p, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim4, n, o, p, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim8, Dim5, m, o, p, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim8, Dim5, i, j, k, l, m, n, m, o, p, q, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(m, o, p, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim4, Dim6, Dim7, Dim8, Dim5, m, o, p, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, m, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim4, Dim7, Dim8, o, n, m, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim7, Dim8, i, j, k, l, m, n, o, n, m, p, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(o, n, m, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim4, Dim7, Dim8, o, n, m, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, n, p, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim5, Dim7, Dim8, o, m, n, p, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim7, Dim8, i, j, k, l, m, n, o, m, n, p, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(o, m, n, p, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim5, Dim7, Dim8, o, m, n, p, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim4, Dim8, o, n, p, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim4, Dim8, i, j, k, l, m, n, o, n, p, m, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(o, n, p, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim4, Dim8, o, n, p, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim5, Dim8, o, m, p, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim5, Dim8, i, j, k, l, m, n, o, m, p, n, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(o, m, p, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim5, Dim8, o, m, p, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim8, Dim4, o, n, p, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim8, Dim4, i, j, k, l, m, n, o, n, p, q, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(o, n, p, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim5, Dim7, Dim8, Dim4, o, n, p, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim8, Dim5, o, m, p, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim8, Dim5, i, j, k, l, m, n, o, m, p, q, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(o, m, p, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim4, Dim7, Dim8, Dim5, o, m, p, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, m, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim4, Dim8, o, p, n, m, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim4, Dim8, i, j, k, l, m, n, o, p, n, m, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(o, p, n, m, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim4, Dim8, o, p, n, m, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, n, q, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim5, Dim8, o, p, m, n, q> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim5, Dim8, i, j, k, l, m, n, o, p, m, n, q, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(o, p, m, n, q, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim5, Dim8, o, p, m, n, q> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, q, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim8, Dim4, o, p, n, q, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim8, Dim4, i, j, k, l, m, n, o, p, n, q, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(o, p, n, q, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim5, Dim8, Dim4, o, p, n, q, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, q, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim8, Dim5, o, p, m, q, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim8, Dim5, i, j, k, l, m, n, o, p, m, q, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(o, p, m, q, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim4, Dim8, Dim5, o, p, m, q, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim5, Dim4, o, p, q, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim5, Dim4, i, j, k, l, m, n, o, p, q, n, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim5, Dim4, o, p, q, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, q, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim4, Dim5, o, p, q, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, Dim4, Dim5, i, j, k, l, m, n, o, p, q, m, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim8, i, j, k, l, m, n, o, p, q>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim8, i, j, k, l, o, p, q>(TensorExpr(a, b));
}

//B(o, p, q, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q>
auto operator*(const Tensor5_Expr<B, U, Dim6, Dim7, Dim8, Dim4, Dim5, o, p, q, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

} //namespace FTensor

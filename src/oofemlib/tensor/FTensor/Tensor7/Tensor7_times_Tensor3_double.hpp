#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, n, o, )*B(p, q, r, ) double contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, int Dim9, char i, char j, char k, char l, char m, char n, char o, char p, char q, char r, int DimA, int DimB, int DimC, int DimD, int DimE, int DimX, int DimY, int DimZ, char a, char b, char c, char d, char e, char x, char y, char z>
class Tensor7_times_Tensor3_double
{
  Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
  Tensor3_Expr<B, U, Dim7, Dim8, Dim9, p, q, r> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        {
          // Permutation is where the indices get checked.
          result += Permutation7<DimA, DimB, DimC, DimD, DimE, DimX, DimY, a, b, c, d, e, x, y>().eval(
            iterA, N1, N2, N3, N4, N5, xx, yy) * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx, yy, N6);
        }
    return result;
  }

  Tensor7_times_Tensor3_double(
    const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &iter_a,
    const Tensor3_Expr<B, U, Dim7, Dim8, Dim9, p, q, r> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, n, o, )*B(j, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim7, j, i, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim7, i, j, k, l, m, n, o, j, i, p, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, k, l, m, n, o, i, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, i, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim7, j, i, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim7, i, j, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, i, j, k, l, m, n, o, i, j, p, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, k, l, m, n, o, i, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, j, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim7, i, j, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim7, Dim0, j, p, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim0, i, j, k, l, m, n, o, j, p, i, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, k, l, m, n, o, i, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, p, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim7, Dim0, j, p, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim7, Dim1, i, p, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim1, i, j, k, l, m, n, o, i, p, j, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, k, l, m, n, o, i, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, p, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim7, Dim1, i, p, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim1, Dim0, p, j, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim0, i, j, k, l, m, n, o, p, j, i, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, k, l, m, n, o, i, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(p, j, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim1, Dim0, p, j, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim0, Dim1, p, i, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim1, i, j, k, l, m, n, o, p, i, j, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, k, l, m, n, o, i, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(p, i, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim0, Dim1, p, i, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim7, k, i, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim7, i, j, k, l, m, n, o, k, i, p, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, j, l, m, n, o, i, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, i, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim7, k, i, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim7, i, k, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, i, j, k, l, m, n, o, i, k, p, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, j, l, m, n, o, i, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, k, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim7, i, k, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim7, Dim0, k, p, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim0, i, j, k, l, m, n, o, k, p, i, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, j, l, m, n, o, i, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, p, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim7, Dim0, k, p, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim7, Dim2, i, p, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim2, i, j, k, l, m, n, o, i, p, k, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, j, l, m, n, o, i, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, p, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim7, Dim2, i, p, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim2, Dim0, p, k, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim0, i, j, k, l, m, n, o, p, k, i, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, j, l, m, n, o, i, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(p, k, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim2, Dim0, p, k, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim0, Dim2, p, i, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim2, i, j, k, l, m, n, o, p, i, k, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, j, l, m, n, o, i, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(p, i, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim0, Dim2, p, i, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim7, l, i, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim7, i, j, k, l, m, n, o, l, i, p, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, j, k, m, n, o, i, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, i, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim7, l, i, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim7, i, l, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, i, j, k, l, m, n, o, i, l, p, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, j, k, m, n, o, i, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(i, l, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim7, i, l, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim7, Dim0, l, p, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim0, i, j, k, l, m, n, o, l, p, i, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, j, k, m, n, o, i, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, p, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim7, Dim0, l, p, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim7, Dim3, i, p, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim3, i, j, k, l, m, n, o, i, p, l, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, j, k, m, n, o, i, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(i, p, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim7, Dim3, i, p, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim3, Dim0, p, l, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim0, i, j, k, l, m, n, o, p, l, i, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, j, k, m, n, o, i, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(p, l, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim3, Dim0, p, l, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim0, Dim3, p, i, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim3, i, j, k, l, m, n, o, p, i, l, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, j, k, m, n, o, i, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(p, i, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim0, Dim3, p, i, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim7, m, i, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim7, i, j, k, l, m, n, o, m, i, p, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, j, k, l, n, o, i, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, i, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim7, m, i, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim7, i, m, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim7, i, j, k, l, m, n, o, i, m, p, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, j, k, l, n, o, i, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(i, m, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim7, i, m, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim7, Dim0, m, p, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim0, i, j, k, l, m, n, o, m, p, i, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, j, k, l, n, o, i, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, p, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim7, Dim0, m, p, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim7, Dim4, i, p, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim4, i, j, k, l, m, n, o, i, p, m, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, j, k, l, n, o, i, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(i, p, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim7, Dim4, i, p, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim4, Dim0, p, m, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim0, i, j, k, l, m, n, o, p, m, i, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, j, k, l, n, o, i, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(p, m, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim4, Dim0, p, m, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim0, Dim4, p, i, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim4, i, j, k, l, m, n, o, p, i, m, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, j, k, l, n, o, i, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(p, i, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim0, Dim4, p, i, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim7, n, i, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim7, i, j, k, l, m, n, o, n, i, p, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, j, k, l, m, o, i, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(n, i, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim7, n, i, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim7, i, n, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim7, i, j, k, l, m, n, o, i, n, p, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, j, k, l, m, o, i, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(i, n, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim7, i, n, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim7, Dim0, n, p, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim0, i, j, k, l, m, n, o, n, p, i, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, j, k, l, m, o, i, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(n, p, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim7, Dim0, n, p, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim7, Dim5, i, p, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim5, i, j, k, l, m, n, o, i, p, n, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, j, k, l, m, o, i, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(i, p, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim7, Dim5, i, p, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim5, Dim0, p, n, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim0, i, j, k, l, m, n, o, p, n, i, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, j, k, l, m, o, i, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(p, n, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim5, Dim0, p, n, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim0, Dim5, p, i, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim5, i, j, k, l, m, n, o, p, i, n, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, j, k, l, m, o, i, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(p, i, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim0, Dim5, p, i, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim0, Dim7, o, i, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim0, Dim7, i, j, k, l, m, n, o, o, i, p, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, j, k, l, m, n, i, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, j, k, l, m, n, p>(TensorExpr(a, b));
}

//B(o, i, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim0, Dim7, o, i, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim6, Dim7, i, o, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim6, Dim7, i, j, k, l, m, n, o, i, o, p, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, j, k, l, m, n, i, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, j, k, l, m, n, p>(TensorExpr(a, b));
}

//B(i, o, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim6, Dim7, i, o, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim7, Dim0, o, p, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim0, i, j, k, l, m, n, o, o, p, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, j, k, l, m, n, i, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, j, k, l, m, n, p>(TensorExpr(a, b));
}

//B(o, p, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim7, Dim0, o, p, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(i, p, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim0, Dim7, Dim6, i, p, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim6, i, j, k, l, m, n, o, i, p, o, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, j, k, l, m, n, i, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, j, k, l, m, n, p>(TensorExpr(a, b));
}

//B(i, p, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim7, Dim6, i, p, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim6, Dim0, p, o, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim0, i, j, k, l, m, n, o, p, o, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, j, k, l, m, n, i, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, j, k, l, m, n, p>(TensorExpr(a, b));
}

//B(p, o, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim6, Dim0, p, o, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim0, Dim6, p, i, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim6, i, j, k, l, m, n, o, p, i, o, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, j, k, l, m, n, i, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, j, k, l, m, n, p>(TensorExpr(a, b));
}

//B(p, i, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim0, Dim6, p, i, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim7, k, j, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim7, i, j, k, l, m, n, o, k, j, p, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, i, l, m, n, o, j, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, j, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim7, k, j, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim7, j, k, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, i, j, k, l, m, n, o, j, k, p, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, i, l, m, n, o, j, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, k, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim7, j, k, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim7, Dim1, k, p, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim1, i, j, k, l, m, n, o, k, p, j, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, i, l, m, n, o, j, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, p, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim7, Dim1, k, p, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim7, Dim2, j, p, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim2, i, j, k, l, m, n, o, j, p, k, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, i, l, m, n, o, j, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, p, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim7, Dim2, j, p, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim2, Dim1, p, k, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim1, i, j, k, l, m, n, o, p, k, j, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, i, l, m, n, o, j, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(p, k, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim2, Dim1, p, k, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim1, Dim2, p, j, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim2, i, j, k, l, m, n, o, p, j, k, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, i, l, m, n, o, j, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(p, j, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim1, Dim2, p, j, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim7, l, j, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim7, i, j, k, l, m, n, o, l, j, p, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, i, k, m, n, o, j, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, j, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim7, l, j, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim7, j, l, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, i, j, k, l, m, n, o, j, l, p, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, i, k, m, n, o, j, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(j, l, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim7, j, l, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim7, Dim1, l, p, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim1, i, j, k, l, m, n, o, l, p, j, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, i, k, m, n, o, j, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, p, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim7, Dim1, l, p, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim7, Dim3, j, p, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim3, i, j, k, l, m, n, o, j, p, l, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, i, k, m, n, o, j, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(j, p, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim7, Dim3, j, p, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim3, Dim1, p, l, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim1, i, j, k, l, m, n, o, p, l, j, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, i, k, m, n, o, j, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(p, l, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim3, Dim1, p, l, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim1, Dim3, p, j, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim3, i, j, k, l, m, n, o, p, j, l, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, i, k, m, n, o, j, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(p, j, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim1, Dim3, p, j, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim7, m, j, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim7, i, j, k, l, m, n, o, m, j, p, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, i, k, l, n, o, j, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, j, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim7, m, j, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim7, j, m, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim7, i, j, k, l, m, n, o, j, m, p, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, i, k, l, n, o, j, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(j, m, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim7, j, m, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim7, Dim1, m, p, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim1, i, j, k, l, m, n, o, m, p, j, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, i, k, l, n, o, j, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, p, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim7, Dim1, m, p, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim7, Dim4, j, p, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim4, i, j, k, l, m, n, o, j, p, m, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, i, k, l, n, o, j, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(j, p, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim7, Dim4, j, p, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim4, Dim1, p, m, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim1, i, j, k, l, m, n, o, p, m, j, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, i, k, l, n, o, j, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(p, m, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim4, Dim1, p, m, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim1, Dim4, p, j, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim4, i, j, k, l, m, n, o, p, j, m, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, i, k, l, n, o, j, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(p, j, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim1, Dim4, p, j, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim7, n, j, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim7, i, j, k, l, m, n, o, n, j, p, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, i, k, l, m, o, j, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(n, j, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim7, n, j, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim7, j, n, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim7, i, j, k, l, m, n, o, j, n, p, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, i, k, l, m, o, j, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(j, n, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim7, j, n, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim7, Dim1, n, p, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim1, i, j, k, l, m, n, o, n, p, j, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, i, k, l, m, o, j, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(n, p, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim7, Dim1, n, p, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim7, Dim5, j, p, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim5, i, j, k, l, m, n, o, j, p, n, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, i, k, l, m, o, j, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(j, p, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim7, Dim5, j, p, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim5, Dim1, p, n, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim1, i, j, k, l, m, n, o, p, n, j, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, i, k, l, m, o, j, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(p, n, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim5, Dim1, p, n, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim1, Dim5, p, j, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim5, i, j, k, l, m, n, o, p, j, n, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, i, k, l, m, o, j, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(p, j, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim1, Dim5, p, j, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim1, Dim7, o, j, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim1, Dim7, i, j, k, l, m, n, o, o, j, p, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, i, k, l, m, n, j, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, i, k, l, m, n, p>(TensorExpr(a, b));
}

//B(o, j, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim1, Dim7, o, j, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim6, Dim7, j, o, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim6, Dim7, i, j, k, l, m, n, o, j, o, p, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, i, k, l, m, n, j, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, i, k, l, m, n, p>(TensorExpr(a, b));
}

//B(j, o, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim6, Dim7, j, o, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim7, Dim1, o, p, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim1, i, j, k, l, m, n, o, o, p, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, i, k, l, m, n, j, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, i, k, l, m, n, p>(TensorExpr(a, b));
}

//B(o, p, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim7, Dim1, o, p, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim1, Dim7, Dim6, j, p, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim6, i, j, k, l, m, n, o, j, p, o, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, i, k, l, m, n, j, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, i, k, l, m, n, p>(TensorExpr(a, b));
}

//B(j, p, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim7, Dim6, j, p, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim6, Dim1, p, o, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim1, i, j, k, l, m, n, o, p, o, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, i, k, l, m, n, j, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, i, k, l, m, n, p>(TensorExpr(a, b));
}

//B(p, o, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim6, Dim1, p, o, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim1, Dim6, p, j, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim6, i, j, k, l, m, n, o, p, j, o, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, i, k, l, m, n, j, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim7, i, k, l, m, n, p>(TensorExpr(a, b));
}

//B(p, j, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim1, Dim6, p, j, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim7, l, k, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim7, i, j, k, l, m, n, o, l, k, p, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, i, j, m, n, o, k, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(l, k, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim7, l, k, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim7, k, l, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, i, j, k, l, m, n, o, k, l, p, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, i, j, m, n, o, k, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(k, l, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim7, k, l, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim7, Dim2, l, p, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim2, i, j, k, l, m, n, o, l, p, k, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, i, j, m, n, o, k, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(l, p, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim7, Dim2, l, p, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim7, Dim3, k, p, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim3, i, j, k, l, m, n, o, k, p, l, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, i, j, m, n, o, k, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(k, p, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim7, Dim3, k, p, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim3, Dim2, p, l, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim2, i, j, k, l, m, n, o, p, l, k, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, i, j, m, n, o, k, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(p, l, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim3, Dim2, p, l, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim2, Dim3, p, k, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim3, i, j, k, l, m, n, o, p, k, l, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, i, j, m, n, o, k, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(p, k, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim2, Dim3, p, k, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim7, m, k, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim7, i, j, k, l, m, n, o, m, k, p, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, i, j, l, n, o, k, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(m, k, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim7, m, k, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim7, k, m, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim7, i, j, k, l, m, n, o, k, m, p, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, i, j, l, n, o, k, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(k, m, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim7, k, m, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim7, Dim2, m, p, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim2, i, j, k, l, m, n, o, m, p, k, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, i, j, l, n, o, k, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(m, p, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim7, Dim2, m, p, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim7, Dim4, k, p, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim4, i, j, k, l, m, n, o, k, p, m, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, i, j, l, n, o, k, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(k, p, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim7, Dim4, k, p, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim4, Dim2, p, m, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim2, i, j, k, l, m, n, o, p, m, k, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, i, j, l, n, o, k, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(p, m, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim4, Dim2, p, m, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim2, Dim4, p, k, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim4, i, j, k, l, m, n, o, p, k, m, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, i, j, l, n, o, k, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(p, k, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim2, Dim4, p, k, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim7, n, k, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim7, i, j, k, l, m, n, o, n, k, p, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, i, j, l, m, o, k, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(n, k, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim7, n, k, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim7, k, n, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim7, i, j, k, l, m, n, o, k, n, p, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, i, j, l, m, o, k, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(k, n, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim7, k, n, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim7, Dim2, n, p, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim2, i, j, k, l, m, n, o, n, p, k, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, i, j, l, m, o, k, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(n, p, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim7, Dim2, n, p, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim7, Dim5, k, p, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim5, i, j, k, l, m, n, o, k, p, n, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, i, j, l, m, o, k, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(k, p, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim7, Dim5, k, p, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim5, Dim2, p, n, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim2, i, j, k, l, m, n, o, p, n, k, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, i, j, l, m, o, k, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(p, n, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim5, Dim2, p, n, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim2, Dim5, p, k, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim5, i, j, k, l, m, n, o, p, k, n, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, i, j, l, m, o, k, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(p, k, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim2, Dim5, p, k, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim2, Dim7, o, k, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim2, Dim7, i, j, k, l, m, n, o, o, k, p, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, i, j, l, m, n, k, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, i, j, l, m, n, p>(TensorExpr(a, b));
}

//B(o, k, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim2, Dim7, o, k, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim6, Dim7, k, o, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim6, Dim7, i, j, k, l, m, n, o, k, o, p, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, i, j, l, m, n, k, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, i, j, l, m, n, p>(TensorExpr(a, b));
}

//B(k, o, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim6, Dim7, k, o, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim7, Dim2, o, p, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim2, i, j, k, l, m, n, o, o, p, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, i, j, l, m, n, k, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, i, j, l, m, n, p>(TensorExpr(a, b));
}

//B(o, p, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim7, Dim2, o, p, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim2, Dim7, Dim6, k, p, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim6, i, j, k, l, m, n, o, k, p, o, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, i, j, l, m, n, k, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, i, j, l, m, n, p>(TensorExpr(a, b));
}

//B(k, p, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim7, Dim6, k, p, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim6, Dim2, p, o, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim2, i, j, k, l, m, n, o, p, o, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, i, j, l, m, n, k, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, i, j, l, m, n, p>(TensorExpr(a, b));
}

//B(p, o, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim6, Dim2, p, o, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim2, Dim6, p, k, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim6, i, j, k, l, m, n, o, p, k, o, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, i, j, l, m, n, k, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim7, i, j, l, m, n, p>(TensorExpr(a, b));
}

//B(p, k, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim2, Dim6, p, k, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim7, m, l, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim7, i, j, k, l, m, n, o, m, l, p, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, i, j, k, n, o, l, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(m, l, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim7, m, l, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim7, l, m, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim7, i, j, k, l, m, n, o, l, m, p, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, i, j, k, n, o, l, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(l, m, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim7, l, m, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim7, Dim3, m, p, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim3, i, j, k, l, m, n, o, m, p, l, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, i, j, k, n, o, l, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(m, p, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim7, Dim3, m, p, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim7, Dim4, l, p, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim4, i, j, k, l, m, n, o, l, p, m, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, i, j, k, n, o, l, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(l, p, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim7, Dim4, l, p, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim4, Dim3, p, m, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim3, i, j, k, l, m, n, o, p, m, l, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, i, j, k, n, o, l, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(p, m, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim4, Dim3, p, m, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim3, Dim4, p, l, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim4, i, j, k, l, m, n, o, p, l, m, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, i, j, k, n, o, l, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(p, l, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim3, Dim4, p, l, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim7, n, l, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim7, i, j, k, l, m, n, o, n, l, p, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, i, j, k, m, o, l, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(n, l, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim7, n, l, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim7, l, n, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim7, i, j, k, l, m, n, o, l, n, p, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, i, j, k, m, o, l, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(l, n, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim7, l, n, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim7, Dim3, n, p, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim3, i, j, k, l, m, n, o, n, p, l, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, i, j, k, m, o, l, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(n, p, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim7, Dim3, n, p, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim7, Dim5, l, p, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim5, i, j, k, l, m, n, o, l, p, n, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, i, j, k, m, o, l, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(l, p, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim7, Dim5, l, p, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim5, Dim3, p, n, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim3, i, j, k, l, m, n, o, p, n, l, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, i, j, k, m, o, l, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(p, n, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim5, Dim3, p, n, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim3, Dim5, p, l, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim5, i, j, k, l, m, n, o, p, l, n, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, i, j, k, m, o, l, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(p, l, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim3, Dim5, p, l, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim3, Dim7, o, l, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim3, Dim7, i, j, k, l, m, n, o, o, l, p, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, i, j, k, m, n, l, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, i, j, k, m, n, p>(TensorExpr(a, b));
}

//B(o, l, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim3, Dim7, o, l, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim6, Dim7, l, o, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim6, Dim7, i, j, k, l, m, n, o, l, o, p, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, i, j, k, m, n, l, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, i, j, k, m, n, p>(TensorExpr(a, b));
}

//B(l, o, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim6, Dim7, l, o, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim7, Dim3, o, p, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim3, i, j, k, l, m, n, o, o, p, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, i, j, k, m, n, l, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, i, j, k, m, n, p>(TensorExpr(a, b));
}

//B(o, p, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim7, Dim3, o, p, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim3, Dim7, Dim6, l, p, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim6, i, j, k, l, m, n, o, l, p, o, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, i, j, k, m, n, l, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, i, j, k, m, n, p>(TensorExpr(a, b));
}

//B(l, p, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim7, Dim6, l, p, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim6, Dim3, p, o, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim3, i, j, k, l, m, n, o, p, o, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, i, j, k, m, n, l, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, i, j, k, m, n, p>(TensorExpr(a, b));
}

//B(p, o, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim6, Dim3, p, o, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim3, Dim6, p, l, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim6, i, j, k, l, m, n, o, p, l, o, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, i, j, k, m, n, l, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim7, i, j, k, m, n, p>(TensorExpr(a, b));
}

//B(p, l, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim3, Dim6, p, l, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim7, n, m, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim7, i, j, k, l, m, n, o, n, m, p, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, i, j, k, l, o, m, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(n, m, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim7, n, m, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim7, m, n, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim7, i, j, k, l, m, n, o, m, n, p, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, i, j, k, l, o, m, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(m, n, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim7, m, n, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim7, Dim4, n, p, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim4, i, j, k, l, m, n, o, n, p, m, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, i, j, k, l, o, m, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(n, p, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim7, Dim4, n, p, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim7, Dim5, m, p, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim5, i, j, k, l, m, n, o, m, p, n, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, i, j, k, l, o, m, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(m, p, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim7, Dim5, m, p, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim5, Dim4, p, n, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim4, i, j, k, l, m, n, o, p, n, m, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, i, j, k, l, o, m, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(p, n, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim5, Dim4, p, n, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim4, Dim5, p, m, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim5, i, j, k, l, m, n, o, p, m, n, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, i, j, k, l, o, m, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(p, m, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim4, Dim5, p, m, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim4, Dim7, o, m, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim4, Dim7, i, j, k, l, m, n, o, o, m, p, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, i, j, k, l, n, m, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, i, j, k, l, n, p>(TensorExpr(a, b));
}

//B(o, m, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim4, Dim7, o, m, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim6, Dim7, m, o, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim6, Dim7, i, j, k, l, m, n, o, m, o, p, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, i, j, k, l, n, m, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, i, j, k, l, n, p>(TensorExpr(a, b));
}

//B(m, o, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim6, Dim7, m, o, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim7, Dim4, o, p, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim4, i, j, k, l, m, n, o, o, p, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, i, j, k, l, n, m, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, i, j, k, l, n, p>(TensorExpr(a, b));
}

//B(o, p, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim7, Dim4, o, p, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim4, Dim7, Dim6, m, p, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim6, i, j, k, l, m, n, o, m, p, o, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, i, j, k, l, n, m, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, i, j, k, l, n, p>(TensorExpr(a, b));
}

//B(m, p, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim7, Dim6, m, p, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim6, Dim4, p, o, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim4, i, j, k, l, m, n, o, p, o, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, i, j, k, l, n, m, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, i, j, k, l, n, p>(TensorExpr(a, b));
}

//B(p, o, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim6, Dim4, p, o, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim4, Dim6, p, m, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim6, i, j, k, l, m, n, o, p, m, o, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, i, j, k, l, n, m, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim7, i, j, k, l, n, p>(TensorExpr(a, b));
}

//B(p, m, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim4, Dim6, p, m, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim5, Dim7, o, n, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim5, Dim7, i, j, k, l, m, n, o, o, n, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, i, j, k, l, m, p>(TensorExpr(a, b));
}

//B(o, n, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim5, Dim7, o, n, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim7, n, o, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, n, o, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, i, j, k, l, m, p>(TensorExpr(a, b));
}

//B(n, o, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim7, n, o, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim6, Dim7, Dim5, o, p, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, Dim5, i, j, k, l, m, n, o, o, p, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, i, j, k, l, m, p>(TensorExpr(a, b));
}

//B(o, p, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim7, Dim5, o, p, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim5, Dim7, Dim6, n, p, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim6, i, j, k, l, m, n, o, n, p, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, i, j, k, l, m, p>(TensorExpr(a, b));
}

//B(n, p, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim7, Dim6, n, p, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim6, Dim5, p, o, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, Dim5, i, j, k, l, m, n, o, p, o, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, i, j, k, l, m, p>(TensorExpr(a, b));
}

//B(p, o, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim6, Dim5, p, o, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor3_Expr<B, U, Dim7, Dim5, Dim6, p, n, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim6, i, j, k, l, m, n, o, p, n, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim7, i, j, k, l, m, p>(TensorExpr(a, b));
}

//B(p, n, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim7, Dim5, Dim6, p, n, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

} //namespace FTensor

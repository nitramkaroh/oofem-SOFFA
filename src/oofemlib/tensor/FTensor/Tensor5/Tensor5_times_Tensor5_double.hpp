#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, )*B(n, o, p, q, r, ) double contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, int Dim9, char i, char j, char k, char l, char m, char n, char o, char p, char q, char r, int DimA, int DimB, int DimC, int DimX, int DimY, int DimZ, int DimW, int DimV, char a, char b, char c, char x, char y, char z, char w, char v>
class Tensor5_times_Tensor5_double
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim9, n, o, p, q, r> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        {
          // Permutation is where the indices get checked.
          result += Permutation5<DimA, DimB, DimC, DimX, DimY, a, b, c, x, y>().eval(
            iterA, N1, N2, N3, xx, yy) * Permutation5<DimX, DimY, DimZ, DimW, DimV, x, y, z, w, v>().eval(iterB, xx, yy, N4, N5, N6);
        }
    return result;
  }

  Tensor5_times_Tensor5_double(
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &iter_a,
    const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim8, Dim9, n, o, p, q, r> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, )*B(j, i, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim0, Dim5, Dim6, Dim7, j, i, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, Dim6, Dim7, i, j, k, l, m, j, i, n, o, p, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, j, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim1, Dim5, Dim6, Dim7, i, j, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, i, j, k, l, m, i, j, n, o, p, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim0, Dim6, Dim7, j, n, i, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, Dim6, Dim7, i, j, k, l, m, j, n, i, o, p, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim1, Dim6, Dim7, i, n, j, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, Dim6, Dim7, i, j, k, l, m, i, n, j, o, p, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim0, Dim7, j, n, o, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim0, Dim7, i, j, k, l, m, j, n, o, i, p, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim1, Dim7, i, n, o, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim1, Dim7, i, j, k, l, m, i, n, o, j, p, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim7, Dim0, j, n, o, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim0, i, j, k, l, m, j, n, o, p, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim7, Dim1, i, n, o, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim1, i, j, k, l, m, i, n, o, p, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim0, Dim6, Dim7, n, j, i, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, Dim7, i, j, k, l, m, n, j, i, o, p, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim1, Dim6, Dim7, n, i, j, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, i, j, k, l, m, n, i, j, o, p, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim0, Dim7, n, j, o, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, Dim7, i, j, k, l, m, n, j, o, i, p, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim1, Dim7, n, i, o, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, Dim7, i, j, k, l, m, n, i, o, j, p, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim7, Dim0, n, j, o, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim0, i, j, k, l, m, n, j, o, p, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim7, Dim1, n, i, o, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim1, i, j, k, l, m, n, i, o, p, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim0, Dim7, n, o, j, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim7, i, j, k, l, m, n, o, j, i, p, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim1, Dim7, n, o, i, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, i, j, k, l, m, n, o, i, j, p, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim7, Dim0, n, o, j, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim0, i, j, k, l, m, n, o, j, p, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim7, Dim1, n, o, i, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim1, i, j, k, l, m, n, o, i, p, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim1, Dim0, n, o, p, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim0, i, j, k, l, m, n, o, p, j, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim0, Dim1, n, o, p, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim1, i, j, k, l, m, n, o, p, i, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim6, Dim7, k, l, m, i, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, i, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim0, Dim5, Dim6, Dim7, k, i, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, Dim6, Dim7, i, j, k, l, m, k, i, n, o, p, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, k, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim2, Dim5, Dim6, Dim7, i, k, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, i, j, k, l, m, i, k, n, o, p, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim0, Dim6, Dim7, k, n, i, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, Dim6, Dim7, i, j, k, l, m, k, n, i, o, p, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim2, Dim6, Dim7, i, n, k, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, Dim6, Dim7, i, j, k, l, m, i, n, k, o, p, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim0, Dim7, k, n, o, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim0, Dim7, i, j, k, l, m, k, n, o, i, p, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim2, Dim7, i, n, o, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim2, Dim7, i, j, k, l, m, i, n, o, k, p, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim7, Dim0, k, n, o, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim0, i, j, k, l, m, k, n, o, p, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim7, Dim2, i, n, o, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim2, i, j, k, l, m, i, n, o, p, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim0, Dim6, Dim7, n, k, i, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, Dim7, i, j, k, l, m, n, k, i, o, p, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim2, Dim6, Dim7, n, i, k, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, i, j, k, l, m, n, i, k, o, p, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim0, Dim7, n, k, o, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, Dim7, i, j, k, l, m, n, k, o, i, p, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim2, Dim7, n, i, o, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, Dim7, i, j, k, l, m, n, i, o, k, p, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim7, Dim0, n, k, o, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim0, i, j, k, l, m, n, k, o, p, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim7, Dim2, n, i, o, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim2, i, j, k, l, m, n, i, o, p, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim0, Dim7, n, o, k, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim7, i, j, k, l, m, n, o, k, i, p, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim2, Dim7, n, o, i, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, i, j, k, l, m, n, o, i, k, p, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim7, Dim0, n, o, k, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim0, i, j, k, l, m, n, o, k, p, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim7, Dim2, n, o, i, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim2, i, j, k, l, m, n, o, i, p, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim2, Dim0, n, o, p, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim0, i, j, k, l, m, n, o, p, k, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim0, Dim2, n, o, p, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim2, i, j, k, l, m, n, o, p, i, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, Dim6, Dim7, j, l, m, i, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, i, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim0, Dim5, Dim6, Dim7, l, i, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, Dim6, Dim7, i, j, k, l, m, l, i, n, o, p, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, l, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim3, Dim5, Dim6, Dim7, i, l, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, i, j, k, l, m, i, l, n, o, p, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim0, Dim6, Dim7, l, n, i, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, Dim6, Dim7, i, j, k, l, m, l, n, i, o, p, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim3, Dim6, Dim7, i, n, l, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, Dim6, Dim7, i, j, k, l, m, i, n, l, o, p, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim0, Dim7, l, n, o, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim0, Dim7, i, j, k, l, m, l, n, o, i, p, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim3, Dim7, i, n, o, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim3, Dim7, i, j, k, l, m, i, n, o, l, p, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim7, Dim0, l, n, o, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim7, Dim0, i, j, k, l, m, l, n, o, p, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim7, Dim3, i, n, o, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim3, i, j, k, l, m, i, n, o, p, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim0, Dim6, Dim7, n, l, i, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim6, Dim7, i, j, k, l, m, n, l, i, o, p, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim3, Dim6, Dim7, n, i, l, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, i, j, k, l, m, n, i, l, o, p, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim0, Dim7, n, l, o, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim0, Dim7, i, j, k, l, m, n, l, o, i, p, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim3, Dim7, n, i, o, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim3, Dim7, i, j, k, l, m, n, i, o, l, p, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim7, Dim0, n, l, o, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim0, i, j, k, l, m, n, l, o, p, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim7, Dim3, n, i, o, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim3, i, j, k, l, m, n, i, o, p, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim0, Dim7, n, o, l, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim7, i, j, k, l, m, n, o, l, i, p, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim3, Dim7, n, o, i, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, i, j, k, l, m, n, o, i, l, p, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim7, Dim0, n, o, l, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim0, i, j, k, l, m, n, o, l, p, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim7, Dim3, n, o, i, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim3, i, j, k, l, m, n, o, i, p, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim3, Dim0, n, o, p, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim0, i, j, k, l, m, n, o, p, l, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim0, Dim3, n, o, p, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim3, i, j, k, l, m, n, o, p, i, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, Dim6, Dim7, j, k, m, i, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, i, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim0, Dim5, Dim6, Dim7, m, i, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim5, Dim6, Dim7, i, j, k, l, m, m, i, n, o, p, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, m, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim4, Dim5, Dim6, Dim7, i, m, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, i, m, n, o, p, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim0, Dim6, Dim7, m, n, i, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim0, Dim6, Dim7, i, j, k, l, m, m, n, i, o, p, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim4, Dim6, Dim7, i, n, m, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim4, Dim6, Dim7, i, j, k, l, m, i, n, m, o, p, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim0, Dim7, m, n, o, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim0, Dim7, i, j, k, l, m, m, n, o, i, p, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim4, Dim7, i, n, o, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim4, Dim7, i, j, k, l, m, i, n, o, m, p, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim7, Dim0, m, n, o, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim7, Dim0, i, j, k, l, m, m, n, o, p, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(i, n, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim0, Dim5, Dim6, Dim7, Dim4, i, n, o, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, Dim4, i, j, k, l, m, i, n, o, p, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim0, Dim6, Dim7, n, m, i, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim6, Dim7, i, j, k, l, m, n, m, i, o, p, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim4, Dim6, Dim7, n, i, m, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim6, Dim7, i, j, k, l, m, n, i, m, o, p, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim0, Dim7, n, m, o, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim0, Dim7, i, j, k, l, m, n, m, o, i, p, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim4, Dim7, n, i, o, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim4, Dim7, i, j, k, l, m, n, i, o, m, p, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim7, Dim0, n, m, o, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim0, i, j, k, l, m, n, m, o, p, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, i, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim0, Dim6, Dim7, Dim4, n, i, o, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim4, i, j, k, l, m, n, i, o, p, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim0, Dim7, n, o, m, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim7, i, j, k, l, m, n, o, m, i, p, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim4, Dim7, n, o, i, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim7, i, j, k, l, m, n, o, i, m, p, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim7, Dim0, n, o, m, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim0, i, j, k, l, m, n, o, m, p, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, i, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim0, Dim7, Dim4, n, o, i, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim4, i, j, k, l, m, n, o, i, p, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim4, Dim0, n, o, p, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim0, i, j, k, l, m, n, o, p, m, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim0, Dim4, n, o, p, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim4, i, j, k, l, m, n, o, p, i, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, Dim7, j, k, l, i, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, j, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim1, Dim5, Dim6, Dim7, k, j, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, Dim6, Dim7, i, j, k, l, m, k, j, n, o, p, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, k, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim2, Dim5, Dim6, Dim7, j, k, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, l, m, j, k, n, o, p, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim1, Dim6, Dim7, k, n, j, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, Dim6, Dim7, i, j, k, l, m, k, n, j, o, p, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim2, Dim6, Dim7, j, n, k, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, Dim6, Dim7, i, j, k, l, m, j, n, k, o, p, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim1, Dim7, k, n, o, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim1, Dim7, i, j, k, l, m, k, n, o, j, p, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim2, Dim7, j, n, o, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim2, Dim7, i, j, k, l, m, j, n, o, k, p, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim7, Dim1, k, n, o, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim1, i, j, k, l, m, k, n, o, p, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim7, Dim2, j, n, o, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim2, i, j, k, l, m, j, n, o, p, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim1, Dim6, Dim7, n, k, j, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, Dim7, i, j, k, l, m, n, k, j, o, p, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim2, Dim6, Dim7, n, j, k, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, i, j, k, l, m, n, j, k, o, p, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim1, Dim7, n, k, o, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, Dim7, i, j, k, l, m, n, k, o, j, p, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim2, Dim7, n, j, o, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, Dim7, i, j, k, l, m, n, j, o, k, p, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim7, Dim1, n, k, o, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim1, i, j, k, l, m, n, k, o, p, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim7, Dim2, n, j, o, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim2, i, j, k, l, m, n, j, o, p, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim1, Dim7, n, o, k, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim7, i, j, k, l, m, n, o, k, j, p, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim2, Dim7, n, o, j, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, i, j, k, l, m, n, o, j, k, p, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim7, Dim1, n, o, k, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim1, i, j, k, l, m, n, o, k, p, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim7, Dim2, n, o, j, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim2, i, j, k, l, m, n, o, j, p, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim2, Dim1, n, o, p, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim1, i, j, k, l, m, n, o, p, k, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim1, Dim2, n, o, p, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim2, i, j, k, l, m, n, o, p, j, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, Dim6, Dim7, i, l, m, j, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, j, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim1, Dim5, Dim6, Dim7, l, j, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, Dim6, Dim7, i, j, k, l, m, l, j, n, o, p, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, l, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim3, Dim5, Dim6, Dim7, j, l, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, k, l, m, j, l, n, o, p, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim1, Dim6, Dim7, l, n, j, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, Dim6, Dim7, i, j, k, l, m, l, n, j, o, p, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim3, Dim6, Dim7, j, n, l, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, Dim6, Dim7, i, j, k, l, m, j, n, l, o, p, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim1, Dim7, l, n, o, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim1, Dim7, i, j, k, l, m, l, n, o, j, p, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim3, Dim7, j, n, o, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim3, Dim7, i, j, k, l, m, j, n, o, l, p, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim7, Dim1, l, n, o, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim7, Dim1, i, j, k, l, m, l, n, o, p, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim7, Dim3, j, n, o, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim3, i, j, k, l, m, j, n, o, p, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim1, Dim6, Dim7, n, l, j, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim6, Dim7, i, j, k, l, m, n, l, j, o, p, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim3, Dim6, Dim7, n, j, l, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, i, j, k, l, m, n, j, l, o, p, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim1, Dim7, n, l, o, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim1, Dim7, i, j, k, l, m, n, l, o, j, p, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim3, Dim7, n, j, o, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim3, Dim7, i, j, k, l, m, n, j, o, l, p, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim7, Dim1, n, l, o, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim1, i, j, k, l, m, n, l, o, p, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim7, Dim3, n, j, o, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim3, i, j, k, l, m, n, j, o, p, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim1, Dim7, n, o, l, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim7, i, j, k, l, m, n, o, l, j, p, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim3, Dim7, n, o, j, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, i, j, k, l, m, n, o, j, l, p, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim7, Dim1, n, o, l, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim1, i, j, k, l, m, n, o, l, p, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim7, Dim3, n, o, j, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim3, i, j, k, l, m, n, o, j, p, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim3, Dim1, n, o, p, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim1, i, j, k, l, m, n, o, p, l, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim1, Dim3, n, o, p, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim3, i, j, k, l, m, n, o, p, j, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, Dim6, Dim7, i, k, m, j, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, j, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim1, Dim5, Dim6, Dim7, m, j, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim5, Dim6, Dim7, i, j, k, l, m, m, j, n, o, p, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, m, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim4, Dim5, Dim6, Dim7, j, m, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, j, m, n, o, p, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim1, Dim6, Dim7, m, n, j, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim1, Dim6, Dim7, i, j, k, l, m, m, n, j, o, p, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim4, Dim6, Dim7, j, n, m, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim4, Dim6, Dim7, i, j, k, l, m, j, n, m, o, p, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim1, Dim7, m, n, o, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim1, Dim7, i, j, k, l, m, m, n, o, j, p, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim4, Dim7, j, n, o, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim4, Dim7, i, j, k, l, m, j, n, o, m, p, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim7, Dim1, m, n, o, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim7, Dim1, i, j, k, l, m, m, n, o, p, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(j, n, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim1, Dim5, Dim6, Dim7, Dim4, j, n, o, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, Dim4, i, j, k, l, m, j, n, o, p, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim1, Dim6, Dim7, n, m, j, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim6, Dim7, i, j, k, l, m, n, m, j, o, p, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim4, Dim6, Dim7, n, j, m, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim6, Dim7, i, j, k, l, m, n, j, m, o, p, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim1, Dim7, n, m, o, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim1, Dim7, i, j, k, l, m, n, m, o, j, p, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim4, Dim7, n, j, o, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim4, Dim7, i, j, k, l, m, n, j, o, m, p, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim7, Dim1, n, m, o, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim1, i, j, k, l, m, n, m, o, p, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, j, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim1, Dim6, Dim7, Dim4, n, j, o, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim4, i, j, k, l, m, n, j, o, p, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim1, Dim7, n, o, m, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim7, i, j, k, l, m, n, o, m, j, p, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim4, Dim7, n, o, j, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim7, i, j, k, l, m, n, o, j, m, p, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim7, Dim1, n, o, m, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim1, i, j, k, l, m, n, o, m, p, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, j, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim1, Dim7, Dim4, n, o, j, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim4, i, j, k, l, m, n, o, j, p, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim4, Dim1, n, o, p, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim1, i, j, k, l, m, n, o, p, m, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim1, Dim4, n, o, p, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim4, i, j, k, l, m, n, o, p, j, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, Dim7, i, k, l, j, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, k, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim2, Dim5, Dim6, Dim7, l, k, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, Dim6, Dim7, i, j, k, l, m, l, k, n, o, p, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, l, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim3, Dim5, Dim6, Dim7, k, l, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, k, l, m, k, l, n, o, p, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim2, Dim6, Dim7, l, n, k, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, Dim6, Dim7, i, j, k, l, m, l, n, k, o, p, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim3, Dim6, Dim7, k, n, l, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, Dim6, Dim7, i, j, k, l, m, k, n, l, o, p, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim2, Dim7, l, n, o, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim2, Dim7, i, j, k, l, m, l, n, o, k, p, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim3, Dim7, k, n, o, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim3, Dim7, i, j, k, l, m, k, n, o, l, p, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim7, Dim2, l, n, o, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim7, Dim2, i, j, k, l, m, l, n, o, p, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim7, Dim3, k, n, o, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim3, i, j, k, l, m, k, n, o, p, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim2, Dim6, Dim7, n, l, k, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim6, Dim7, i, j, k, l, m, n, l, k, o, p, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim3, Dim6, Dim7, n, k, l, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, i, j, k, l, m, n, k, l, o, p, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim2, Dim7, n, l, o, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim2, Dim7, i, j, k, l, m, n, l, o, k, p, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim3, Dim7, n, k, o, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim3, Dim7, i, j, k, l, m, n, k, o, l, p, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim7, Dim2, n, l, o, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim2, i, j, k, l, m, n, l, o, p, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim7, Dim3, n, k, o, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim3, i, j, k, l, m, n, k, o, p, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim2, Dim7, n, o, l, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim7, i, j, k, l, m, n, o, l, k, p, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim3, Dim7, n, o, k, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, i, j, k, l, m, n, o, k, l, p, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim7, Dim2, n, o, l, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim2, i, j, k, l, m, n, o, l, p, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim7, Dim3, n, o, k, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim3, i, j, k, l, m, n, o, k, p, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim3, Dim2, n, o, p, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim2, i, j, k, l, m, n, o, p, l, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim2, Dim3, n, o, p, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim3, i, j, k, l, m, n, o, p, k, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, m, k, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, k, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim2, Dim5, Dim6, Dim7, m, k, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim5, Dim6, Dim7, i, j, k, l, m, m, k, n, o, p, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, m, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim4, Dim5, Dim6, Dim7, k, m, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, k, m, n, o, p, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim2, Dim6, Dim7, m, n, k, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim2, Dim6, Dim7, i, j, k, l, m, m, n, k, o, p, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim4, Dim6, Dim7, k, n, m, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim4, Dim6, Dim7, i, j, k, l, m, k, n, m, o, p, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim2, Dim7, m, n, o, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim2, Dim7, i, j, k, l, m, m, n, o, k, p, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim4, Dim7, k, n, o, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim4, Dim7, i, j, k, l, m, k, n, o, m, p, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim7, Dim2, m, n, o, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim7, Dim2, i, j, k, l, m, m, n, o, p, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(k, n, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim2, Dim5, Dim6, Dim7, Dim4, k, n, o, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, Dim4, i, j, k, l, m, k, n, o, p, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim2, Dim6, Dim7, n, m, k, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim6, Dim7, i, j, k, l, m, n, m, k, o, p, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim4, Dim6, Dim7, n, k, m, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim6, Dim7, i, j, k, l, m, n, k, m, o, p, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim2, Dim7, n, m, o, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim2, Dim7, i, j, k, l, m, n, m, o, k, p, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim4, Dim7, n, k, o, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim4, Dim7, i, j, k, l, m, n, k, o, m, p, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim7, Dim2, n, m, o, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim2, i, j, k, l, m, n, m, o, p, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, k, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim2, Dim6, Dim7, Dim4, n, k, o, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim4, i, j, k, l, m, n, k, o, p, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim2, Dim7, n, o, m, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim7, i, j, k, l, m, n, o, m, k, p, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim4, Dim7, n, o, k, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim7, i, j, k, l, m, n, o, k, m, p, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim7, Dim2, n, o, m, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim2, i, j, k, l, m, n, o, m, p, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, k, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim2, Dim7, Dim4, n, o, k, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim4, i, j, k, l, m, n, o, k, p, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim4, Dim2, n, o, p, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim2, i, j, k, l, m, n, o, p, m, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim2, Dim4, n, o, p, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim4, i, j, k, l, m, n, o, p, k, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, l, k, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, l, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim3, Dim5, Dim6, Dim7, m, l, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim5, Dim6, Dim7, i, j, k, l, m, m, l, n, o, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, m, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim4, Dim5, Dim6, Dim7, l, m, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, l, m, n, o, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim3, Dim6, Dim7, m, n, l, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim3, Dim6, Dim7, i, j, k, l, m, m, n, l, o, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim4, Dim6, Dim7, l, n, m, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim4, Dim6, Dim7, i, j, k, l, m, l, n, m, o, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim3, Dim7, m, n, o, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim3, Dim7, i, j, k, l, m, m, n, o, l, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim4, Dim7, l, n, o, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim4, Dim7, i, j, k, l, m, l, n, o, m, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(m, n, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim4, Dim5, Dim6, Dim7, Dim3, m, n, o, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim7, Dim3, i, j, k, l, m, m, n, o, p, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(l, n, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim3, Dim5, Dim6, Dim7, Dim4, l, n, o, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim7, Dim4, i, j, k, l, m, l, n, o, p, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim3, Dim6, Dim7, n, m, l, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim6, Dim7, i, j, k, l, m, n, m, l, o, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim4, Dim6, Dim7, n, l, m, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim6, Dim7, i, j, k, l, m, n, l, m, o, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim3, Dim7, n, m, o, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim3, Dim7, i, j, k, l, m, n, m, o, l, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim4, Dim7, n, l, o, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim4, Dim7, i, j, k, l, m, n, l, o, m, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, m, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim4, Dim6, Dim7, Dim3, n, m, o, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim3, i, j, k, l, m, n, m, o, p, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, l, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim3, Dim6, Dim7, Dim4, n, l, o, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim4, i, j, k, l, m, n, l, o, p, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim3, Dim7, n, o, m, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim7, i, j, k, l, m, n, o, m, l, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim4, Dim7, n, o, l, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim7, i, j, k, l, m, n, o, l, m, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, m, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim4, Dim7, Dim3, n, o, m, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim3, i, j, k, l, m, n, o, m, p, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, l, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim3, Dim7, Dim4, n, o, l, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim4, i, j, k, l, m, n, o, l, p, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim4, Dim3, n, o, p, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim3, i, j, k, l, m, n, o, p, m, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//A(i, j, k, l, m, )*B(n, o, p, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor5_Expr<B, U, Dim5, Dim6, Dim7, Dim3, Dim4, n, o, p, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor5_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim4, i, j, k, l, m, n, o, p, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

} //namespace FTensor

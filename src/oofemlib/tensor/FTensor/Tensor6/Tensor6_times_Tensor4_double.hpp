#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, n, )*B(o, p, q, r, ) double contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, int Dim9, char i, char j, char k, char l, char m, char n, char o, char p, char q, char r, int DimA, int DimB, int DimC, int DimD, int DimX, int DimY, int DimZ, int DimW, char a, char b, char c, char d, char x, char y, char z, char w>
class Tensor6_times_Tensor4_double
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim9, o, p, q, r> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX, ++xx)
      for(int yy = 0; yy < DimY, ++yy)
        {
          // Permutation is where the indices get checked.
          result += Permutation6<DimA, DimB, DimC, DimD, DimX, DimY, a, b, c, d, x, y>().eval(
            iterA, N1, N2, N3, N4, xx, yy) * Permutation4<DimX, DimY, DimZ, DimW, x, y, z, w>().eval(iterB, xx, yy, N5, N6);
        }
    return result;
  }

  Tensor6_times_Tensor4_double(
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &iter_a,
    const Tensor4_Expr<B, U, Dim6, Dim7, Dim8, Dim9, o, p, q, r> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, n, )*B(j, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim6, Dim7, j, i, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, Dim7, i, j, k, l, m, n, j, i, o, p, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim6, Dim7, k, l, m, n, j, i, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, i, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim6, Dim7, j, i, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim6, Dim7, i, j, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, i, j, k, l, m, n, i, j, o, p, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim6, Dim7, k, l, m, n, i, j, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, j, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim6, Dim7, i, j, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim0, Dim7, j, o, i, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, Dim7, i, j, k, l, m, n, j, o, i, p, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim0, Dim7, k, l, m, n, j, o, i, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, o, i, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim0, Dim7, j, o, i, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim1, Dim7, i, o, j, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, Dim7, i, j, k, l, m, n, i, o, j, p, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim1, Dim7, k, l, m, n, i, o, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, o, j, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim1, Dim7, i, o, j, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim7, Dim0, j, o, p, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim0, i, j, k, l, m, n, j, o, p, i, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim0, k, l, m, n, j, o, p, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, o, p, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim7, Dim0, j, o, p, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim7, Dim1, i, o, p, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim1, i, j, k, l, m, n, i, o, p, j, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim1, k, l, m, n, i, o, p, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, o, p, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim7, Dim1, i, o, p, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim0, Dim7, o, j, i, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim7, i, j, k, l, m, n, o, j, i, p, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim0, Dim7, k, l, m, n, o, j, i, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, j, i, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim0, Dim7, o, j, i, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim1, Dim7, o, i, j, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, i, j, k, l, m, n, o, i, j, p, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim1, Dim7, k, l, m, n, o, i, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, i, j, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim1, Dim7, o, i, j, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim7, Dim0, o, j, p, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim0, i, j, k, l, m, n, o, j, p, i, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim0, k, l, m, n, o, j, p, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, j, p, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim7, Dim0, o, j, p, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim7, Dim1, o, i, p, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim1, i, j, k, l, m, n, o, i, p, j, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim1, k, l, m, n, o, i, p, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, i, p, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim7, Dim1, o, i, p, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim1, Dim0, o, p, j, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim0, i, j, k, l, m, n, o, p, j, i, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim0, k, l, m, n, o, p, j, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, j, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim1, Dim0, o, p, j, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim0, Dim1, o, p, i, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim1, i, j, k, l, m, n, o, p, i, j, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim1, k, l, m, n, o, p, i, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, i, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim0, Dim1, o, p, i, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim6, Dim7, k, i, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, Dim7, i, j, k, l, m, n, k, i, o, p, Dim1, Dim3, Dim4, Dim5, Dim2, Dim0, Dim6, Dim7, j, l, m, n, k, i, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, i, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim6, Dim7, k, i, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim6, Dim7, i, k, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, i, j, k, l, m, n, i, k, o, p, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, Dim6, Dim7, j, l, m, n, i, k, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, k, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim6, Dim7, i, k, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim0, Dim7, k, o, i, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, Dim7, i, j, k, l, m, n, k, o, i, p, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim0, Dim7, j, l, m, n, k, o, i, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, o, i, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim0, Dim7, k, o, i, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim2, Dim7, i, o, k, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, Dim7, i, j, k, l, m, n, i, o, k, p, Dim1, Dim3, Dim4, Dim5, Dim0, Dim6, Dim2, Dim7, j, l, m, n, i, o, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, o, k, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim2, Dim7, i, o, k, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim7, Dim0, k, o, p, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim0, i, j, k, l, m, n, k, o, p, i, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim0, j, l, m, n, k, o, p, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, o, p, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim7, Dim0, k, o, p, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim7, Dim2, i, o, p, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim2, i, j, k, l, m, n, i, o, p, k, Dim1, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim2, j, l, m, n, i, o, p, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, o, p, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim7, Dim2, i, o, p, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim0, Dim7, o, k, i, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim7, i, j, k, l, m, n, o, k, i, p, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, Dim0, Dim7, j, l, m, n, o, k, i, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, k, i, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim0, Dim7, o, k, i, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim2, Dim7, o, i, k, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, i, j, k, l, m, n, o, i, k, p, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim2, Dim7, j, l, m, n, o, i, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, i, k, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim2, Dim7, o, i, k, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim7, Dim0, o, k, p, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim0, i, j, k, l, m, n, o, k, p, i, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim0, j, l, m, n, o, k, p, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, k, p, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim7, Dim0, o, k, p, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim7, Dim2, o, i, p, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim2, i, j, k, l, m, n, o, i, p, k, Dim1, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim2, j, l, m, n, o, i, p, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, i, p, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim7, Dim2, o, i, p, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim2, Dim0, o, p, k, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim0, i, j, k, l, m, n, o, p, k, i, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim0, j, l, m, n, o, p, k, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, k, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim2, Dim0, o, p, k, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim0, Dim2, o, p, i, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim2, i, j, k, l, m, n, o, p, i, k, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim2, j, l, m, n, o, p, i, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, i, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim0, Dim2, o, p, i, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim6, Dim7, l, i, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim6, Dim7, i, j, k, l, m, n, l, i, o, p, Dim1, Dim2, Dim4, Dim5, Dim3, Dim0, Dim6, Dim7, j, k, m, n, l, i, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, i, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim6, Dim7, l, i, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim6, Dim7, i, l, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, i, j, k, l, m, n, i, l, o, p, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, Dim6, Dim7, j, k, m, n, i, l, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(i, l, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim6, Dim7, i, l, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim0, Dim7, l, o, i, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim0, Dim7, i, j, k, l, m, n, l, o, i, p, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim0, Dim7, j, k, m, n, l, o, i, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, o, i, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim0, Dim7, l, o, i, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim3, Dim7, i, o, l, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim3, Dim7, i, j, k, l, m, n, i, o, l, p, Dim1, Dim2, Dim4, Dim5, Dim0, Dim6, Dim3, Dim7, j, k, m, n, i, o, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(i, o, l, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim3, Dim7, i, o, l, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim7, Dim0, l, o, p, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim0, i, j, k, l, m, n, l, o, p, i, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim0, j, k, m, n, l, o, p, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, o, p, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim7, Dim0, l, o, p, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim7, Dim3, i, o, p, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim3, i, j, k, l, m, n, i, o, p, l, Dim1, Dim2, Dim4, Dim5, Dim0, Dim6, Dim7, Dim3, j, k, m, n, i, o, p, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(i, o, p, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim7, Dim3, i, o, p, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim0, Dim7, o, l, i, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim0, Dim7, i, j, k, l, m, n, o, l, i, p, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, Dim0, Dim7, j, k, m, n, o, l, i, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, l, i, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim0, Dim7, o, l, i, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim3, Dim7, o, i, l, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, i, j, k, l, m, n, o, i, l, p, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim3, Dim7, j, k, m, n, o, i, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, i, l, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim3, Dim7, o, i, l, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim7, Dim0, o, l, p, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim0, i, j, k, l, m, n, o, l, p, i, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, Dim7, Dim0, j, k, m, n, o, l, p, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, l, p, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim7, Dim0, o, l, p, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim7, Dim3, o, i, p, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim3, i, j, k, l, m, n, o, i, p, l, Dim1, Dim2, Dim4, Dim5, Dim6, Dim0, Dim7, Dim3, j, k, m, n, o, i, p, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, i, p, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim7, Dim3, o, i, p, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim3, Dim0, o, p, l, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim0, i, j, k, l, m, n, o, p, l, i, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim3, Dim0, j, k, m, n, o, p, l, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, l, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim3, Dim0, o, p, l, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim0, Dim3, o, p, i, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim3, i, j, k, l, m, n, o, p, i, l, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim0, Dim3, j, k, m, n, o, p, i, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, i, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim0, Dim3, o, p, i, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim6, Dim7, m, i, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim6, Dim7, i, j, k, l, m, n, m, i, o, p, Dim1, Dim2, Dim3, Dim5, Dim4, Dim0, Dim6, Dim7, j, k, l, n, m, i, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, i, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim6, Dim7, m, i, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim6, Dim7, i, m, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim6, Dim7, i, j, k, l, m, n, i, m, o, p, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, Dim6, Dim7, j, k, l, n, i, m, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(i, m, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim6, Dim7, i, m, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim0, Dim7, m, o, i, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim0, Dim7, i, j, k, l, m, n, m, o, i, p, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim0, Dim7, j, k, l, n, m, o, i, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, o, i, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim0, Dim7, m, o, i, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim4, Dim7, i, o, m, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim4, Dim7, i, j, k, l, m, n, i, o, m, p, Dim1, Dim2, Dim3, Dim5, Dim0, Dim6, Dim4, Dim7, j, k, l, n, i, o, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(i, o, m, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim4, Dim7, i, o, m, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim7, Dim0, m, o, p, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim0, i, j, k, l, m, n, m, o, p, i, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim0, j, k, l, n, m, o, p, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, o, p, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim7, Dim0, m, o, p, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim7, Dim4, i, o, p, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim4, i, j, k, l, m, n, i, o, p, m, Dim1, Dim2, Dim3, Dim5, Dim0, Dim6, Dim7, Dim4, j, k, l, n, i, o, p, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(i, o, p, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim7, Dim4, i, o, p, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim0, Dim7, o, m, i, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim0, Dim7, i, j, k, l, m, n, o, m, i, p, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, Dim0, Dim7, j, k, l, n, o, m, i, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, m, i, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim0, Dim7, o, m, i, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim4, Dim7, o, i, m, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim4, Dim7, i, j, k, l, m, n, o, i, m, p, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim4, Dim7, j, k, l, n, o, i, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, i, m, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim4, Dim7, o, i, m, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim7, Dim0, o, m, p, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim0, i, j, k, l, m, n, o, m, p, i, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, Dim7, Dim0, j, k, l, n, o, m, p, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, m, p, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim7, Dim0, o, m, p, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim7, Dim4, o, i, p, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim4, i, j, k, l, m, n, o, i, p, m, Dim1, Dim2, Dim3, Dim5, Dim6, Dim0, Dim7, Dim4, j, k, l, n, o, i, p, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, i, p, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim7, Dim4, o, i, p, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim4, Dim0, o, p, m, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim0, i, j, k, l, m, n, o, p, m, i, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim4, Dim0, j, k, l, n, o, p, m, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, p, m, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim4, Dim0, o, p, m, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim0, Dim4, o, p, i, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim4, i, j, k, l, m, n, o, p, i, m, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim0, Dim4, j, k, l, n, o, p, i, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, p, i, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim0, Dim4, o, p, i, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim7, n, i, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, Dim6, Dim7, i, j, k, l, m, n, n, i, o, p, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, j, k, l, m, n, i, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(n, i, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim7, n, i, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim7, i, n, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, Dim6, Dim7, i, j, k, l, m, n, i, n, o, p, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, j, k, l, m, i, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(i, n, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim7, i, n, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim7, n, o, i, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim0, Dim7, i, j, k, l, m, n, n, o, i, p, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, j, k, l, m, n, o, i, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(n, o, i, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim7, n, o, i, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim5, Dim7, i, o, n, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim5, Dim7, i, j, k, l, m, n, i, o, n, p, Dim1, Dim2, Dim3, Dim4, Dim0, Dim6, Dim5, Dim7, j, k, l, m, i, o, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(i, o, n, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim5, Dim7, i, o, n, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim0, n, o, p, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim0, i, j, k, l, m, n, n, o, p, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, j, k, l, m, n, o, p, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(n, o, p, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim0, n, o, p, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(i, o, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim0, Dim6, Dim7, Dim5, i, o, p, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, Dim5, i, j, k, l, m, n, i, o, p, n, Dim1, Dim2, Dim3, Dim4, Dim0, Dim6, Dim7, Dim5, j, k, l, m, i, o, p, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(i, o, p, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim6, Dim7, Dim5, i, o, p, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim0, Dim7, o, n, i, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim0, Dim7, i, j, k, l, m, n, o, n, i, p, Dim1, Dim2, Dim3, Dim4, Dim6, Dim5, Dim0, Dim7, j, k, l, m, o, n, i, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, n, i, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim0, Dim7, o, n, i, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim5, Dim7, o, i, n, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim5, Dim7, i, j, k, l, m, n, o, i, n, p, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim5, Dim7, j, k, l, m, o, i, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, i, n, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim5, Dim7, o, i, n, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim7, Dim0, o, n, p, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim0, i, j, k, l, m, n, o, n, p, i, Dim1, Dim2, Dim3, Dim4, Dim6, Dim5, Dim7, Dim0, j, k, l, m, o, n, p, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, n, p, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim7, Dim0, o, n, p, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim0, Dim7, Dim5, o, i, p, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, Dim5, i, j, k, l, m, n, o, i, p, n, Dim1, Dim2, Dim3, Dim4, Dim6, Dim0, Dim7, Dim5, j, k, l, m, o, i, p, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, i, p, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim0, Dim7, Dim5, o, i, p, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim5, Dim0, o, p, n, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim0, i, j, k, l, m, n, o, p, n, i, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim5, Dim0, j, k, l, m, o, p, n, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, p, n, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim5, Dim0, o, p, n, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim0, Dim5, o, p, i, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, Dim5, i, j, k, l, m, n, o, p, i, n, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim0, Dim5, j, k, l, m, o, p, i, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, p, i, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim0, Dim5, o, p, i, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim6, Dim7, k, j, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, Dim7, i, j, k, l, m, n, k, j, o, p, Dim0, Dim3, Dim4, Dim5, Dim2, Dim1, Dim6, Dim7, i, l, m, n, k, j, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, j, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim6, Dim7, k, j, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim6, Dim7, j, k, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, i, j, k, l, m, n, j, k, o, p, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, Dim6, Dim7, i, l, m, n, j, k, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, k, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim6, Dim7, j, k, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim1, Dim7, k, o, j, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, Dim7, i, j, k, l, m, n, k, o, j, p, Dim0, Dim3, Dim4, Dim5, Dim2, Dim6, Dim1, Dim7, i, l, m, n, k, o, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, o, j, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim1, Dim7, k, o, j, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim2, Dim7, j, o, k, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, Dim7, i, j, k, l, m, n, j, o, k, p, Dim0, Dim3, Dim4, Dim5, Dim1, Dim6, Dim2, Dim7, i, l, m, n, j, o, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, o, k, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim2, Dim7, j, o, k, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim7, Dim1, k, o, p, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim1, i, j, k, l, m, n, k, o, p, j, Dim0, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim1, i, l, m, n, k, o, p, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, o, p, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim7, Dim1, k, o, p, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim7, Dim2, j, o, p, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim2, i, j, k, l, m, n, j, o, p, k, Dim0, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim2, i, l, m, n, j, o, p, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, o, p, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim7, Dim2, j, o, p, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim1, Dim7, o, k, j, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim7, i, j, k, l, m, n, o, k, j, p, Dim0, Dim3, Dim4, Dim5, Dim6, Dim2, Dim1, Dim7, i, l, m, n, o, k, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, k, j, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim1, Dim7, o, k, j, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim2, Dim7, o, j, k, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, i, j, k, l, m, n, o, j, k, p, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim2, Dim7, i, l, m, n, o, j, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, j, k, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim2, Dim7, o, j, k, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim7, Dim1, o, k, p, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim1, i, j, k, l, m, n, o, k, p, j, Dim0, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim1, i, l, m, n, o, k, p, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, k, p, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim7, Dim1, o, k, p, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim7, Dim2, o, j, p, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim2, i, j, k, l, m, n, o, j, p, k, Dim0, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim2, i, l, m, n, o, j, p, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, j, p, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim7, Dim2, o, j, p, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim2, Dim1, o, p, k, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim1, i, j, k, l, m, n, o, p, k, j, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim1, i, l, m, n, o, p, k, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, k, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim2, Dim1, o, p, k, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim1, Dim2, o, p, j, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim2, i, j, k, l, m, n, o, p, j, k, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim2, i, l, m, n, o, p, j, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, Dim6, Dim7, i, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, j, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim1, Dim2, o, p, j, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim6, Dim7, l, j, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim6, Dim7, i, j, k, l, m, n, l, j, o, p, Dim0, Dim2, Dim4, Dim5, Dim3, Dim1, Dim6, Dim7, i, k, m, n, l, j, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, j, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim6, Dim7, l, j, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim6, Dim7, j, l, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, i, j, k, l, m, n, j, l, o, p, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, Dim6, Dim7, i, k, m, n, j, l, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(j, l, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim6, Dim7, j, l, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim1, Dim7, l, o, j, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim1, Dim7, i, j, k, l, m, n, l, o, j, p, Dim0, Dim2, Dim4, Dim5, Dim3, Dim6, Dim1, Dim7, i, k, m, n, l, o, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, o, j, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim1, Dim7, l, o, j, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim3, Dim7, j, o, l, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim3, Dim7, i, j, k, l, m, n, j, o, l, p, Dim0, Dim2, Dim4, Dim5, Dim1, Dim6, Dim3, Dim7, i, k, m, n, j, o, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(j, o, l, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim3, Dim7, j, o, l, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim7, Dim1, l, o, p, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim1, i, j, k, l, m, n, l, o, p, j, Dim0, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, Dim1, i, k, m, n, l, o, p, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, o, p, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim7, Dim1, l, o, p, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim7, Dim3, j, o, p, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim3, i, j, k, l, m, n, j, o, p, l, Dim0, Dim2, Dim4, Dim5, Dim1, Dim6, Dim7, Dim3, i, k, m, n, j, o, p, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(j, o, p, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim7, Dim3, j, o, p, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim1, Dim7, o, l, j, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim1, Dim7, i, j, k, l, m, n, o, l, j, p, Dim0, Dim2, Dim4, Dim5, Dim6, Dim3, Dim1, Dim7, i, k, m, n, o, l, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, l, j, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim1, Dim7, o, l, j, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim3, Dim7, o, j, l, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, i, j, k, l, m, n, o, j, l, p, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim3, Dim7, i, k, m, n, o, j, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, j, l, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim3, Dim7, o, j, l, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim7, Dim1, o, l, p, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim1, i, j, k, l, m, n, o, l, p, j, Dim0, Dim2, Dim4, Dim5, Dim6, Dim3, Dim7, Dim1, i, k, m, n, o, l, p, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, l, p, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim7, Dim1, o, l, p, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim7, Dim3, o, j, p, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim3, i, j, k, l, m, n, o, j, p, l, Dim0, Dim2, Dim4, Dim5, Dim6, Dim1, Dim7, Dim3, i, k, m, n, o, j, p, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, j, p, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim7, Dim3, o, j, p, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim3, Dim1, o, p, l, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim1, i, j, k, l, m, n, o, p, l, j, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim3, Dim1, i, k, m, n, o, p, l, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, l, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim3, Dim1, o, p, l, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim1, Dim3, o, p, j, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim3, i, j, k, l, m, n, o, p, j, l, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, Dim1, Dim3, i, k, m, n, o, p, j, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, Dim6, Dim7, i, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, j, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim1, Dim3, o, p, j, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim6, Dim7, m, j, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim6, Dim7, i, j, k, l, m, n, m, j, o, p, Dim0, Dim2, Dim3, Dim5, Dim4, Dim1, Dim6, Dim7, i, k, l, n, m, j, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, j, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim6, Dim7, m, j, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim6, Dim7, j, m, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim6, Dim7, i, j, k, l, m, n, j, m, o, p, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, Dim6, Dim7, i, k, l, n, j, m, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(j, m, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim6, Dim7, j, m, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim1, Dim7, m, o, j, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim1, Dim7, i, j, k, l, m, n, m, o, j, p, Dim0, Dim2, Dim3, Dim5, Dim4, Dim6, Dim1, Dim7, i, k, l, n, m, o, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, o, j, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim1, Dim7, m, o, j, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim4, Dim7, j, o, m, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim4, Dim7, i, j, k, l, m, n, j, o, m, p, Dim0, Dim2, Dim3, Dim5, Dim1, Dim6, Dim4, Dim7, i, k, l, n, j, o, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(j, o, m, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim4, Dim7, j, o, m, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim7, Dim1, m, o, p, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim1, i, j, k, l, m, n, m, o, p, j, Dim0, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, Dim1, i, k, l, n, m, o, p, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, o, p, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim7, Dim1, m, o, p, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim7, Dim4, j, o, p, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim4, i, j, k, l, m, n, j, o, p, m, Dim0, Dim2, Dim3, Dim5, Dim1, Dim6, Dim7, Dim4, i, k, l, n, j, o, p, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(j, o, p, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim7, Dim4, j, o, p, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim1, Dim7, o, m, j, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim1, Dim7, i, j, k, l, m, n, o, m, j, p, Dim0, Dim2, Dim3, Dim5, Dim6, Dim4, Dim1, Dim7, i, k, l, n, o, m, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, m, j, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim1, Dim7, o, m, j, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim4, Dim7, o, j, m, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim4, Dim7, i, j, k, l, m, n, o, j, m, p, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim4, Dim7, i, k, l, n, o, j, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, j, m, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim4, Dim7, o, j, m, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim7, Dim1, o, m, p, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim1, i, j, k, l, m, n, o, m, p, j, Dim0, Dim2, Dim3, Dim5, Dim6, Dim4, Dim7, Dim1, i, k, l, n, o, m, p, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, m, p, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim7, Dim1, o, m, p, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim7, Dim4, o, j, p, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim4, i, j, k, l, m, n, o, j, p, m, Dim0, Dim2, Dim3, Dim5, Dim6, Dim1, Dim7, Dim4, i, k, l, n, o, j, p, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, j, p, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim7, Dim4, o, j, p, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim4, Dim1, o, p, m, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim1, i, j, k, l, m, n, o, p, m, j, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim4, Dim1, i, k, l, n, o, p, m, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, p, m, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim4, Dim1, o, p, m, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim1, Dim4, o, p, j, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim4, i, j, k, l, m, n, o, p, j, m, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, Dim1, Dim4, i, k, l, n, o, p, j, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, Dim6, Dim7, i, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, p, j, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim1, Dim4, o, p, j, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim7, n, j, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, Dim6, Dim7, i, j, k, l, m, n, n, j, o, p, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, i, k, l, m, n, j, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(n, j, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim7, n, j, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim7, j, n, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, Dim6, Dim7, i, j, k, l, m, n, j, n, o, p, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, i, k, l, m, j, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(j, n, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim7, j, n, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim7, n, o, j, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim1, Dim7, i, j, k, l, m, n, n, o, j, p, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, i, k, l, m, n, o, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(n, o, j, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim7, n, o, j, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim5, Dim7, j, o, n, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim5, Dim7, i, j, k, l, m, n, j, o, n, p, Dim0, Dim2, Dim3, Dim4, Dim1, Dim6, Dim5, Dim7, i, k, l, m, j, o, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(j, o, n, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim5, Dim7, j, o, n, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim1, n, o, p, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim1, i, j, k, l, m, n, n, o, p, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, i, k, l, m, n, o, p, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(n, o, p, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim1, n, o, p, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim1, Dim6, Dim7, Dim5, j, o, p, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, Dim5, i, j, k, l, m, n, j, o, p, n, Dim0, Dim2, Dim3, Dim4, Dim1, Dim6, Dim7, Dim5, i, k, l, m, j, o, p, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(j, o, p, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim6, Dim7, Dim5, j, o, p, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim1, Dim7, o, n, j, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim1, Dim7, i, j, k, l, m, n, o, n, j, p, Dim0, Dim2, Dim3, Dim4, Dim6, Dim5, Dim1, Dim7, i, k, l, m, o, n, j, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, n, j, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim1, Dim7, o, n, j, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim5, Dim7, o, j, n, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim5, Dim7, i, j, k, l, m, n, o, j, n, p, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim5, Dim7, i, k, l, m, o, j, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, j, n, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim5, Dim7, o, j, n, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim7, Dim1, o, n, p, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim1, i, j, k, l, m, n, o, n, p, j, Dim0, Dim2, Dim3, Dim4, Dim6, Dim5, Dim7, Dim1, i, k, l, m, o, n, p, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, n, p, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim7, Dim1, o, n, p, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim1, Dim7, Dim5, o, j, p, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, Dim5, i, j, k, l, m, n, o, j, p, n, Dim0, Dim2, Dim3, Dim4, Dim6, Dim1, Dim7, Dim5, i, k, l, m, o, j, p, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, j, p, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim1, Dim7, Dim5, o, j, p, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim5, Dim1, o, p, n, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim1, i, j, k, l, m, n, o, p, n, j, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim5, Dim1, i, k, l, m, o, p, n, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, p, n, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim5, Dim1, o, p, n, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim1, Dim5, o, p, j, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, Dim5, i, j, k, l, m, n, o, p, j, n, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, Dim1, Dim5, i, k, l, m, o, p, j, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim6, Dim7, i, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, p, j, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim1, Dim5, o, p, j, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim6, Dim7, l, k, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim6, Dim7, i, j, k, l, m, n, l, k, o, p, Dim0, Dim1, Dim4, Dim5, Dim3, Dim2, Dim6, Dim7, i, j, m, n, l, k, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(l, k, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim6, Dim7, l, k, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim6, Dim7, k, l, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, i, j, k, l, m, n, k, l, o, p, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, Dim6, Dim7, i, j, m, n, k, l, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(k, l, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim6, Dim7, k, l, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim2, Dim7, l, o, k, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim2, Dim7, i, j, k, l, m, n, l, o, k, p, Dim0, Dim1, Dim4, Dim5, Dim3, Dim6, Dim2, Dim7, i, j, m, n, l, o, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(l, o, k, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim2, Dim7, l, o, k, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim3, Dim7, k, o, l, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim3, Dim7, i, j, k, l, m, n, k, o, l, p, Dim0, Dim1, Dim4, Dim5, Dim2, Dim6, Dim3, Dim7, i, j, m, n, k, o, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(k, o, l, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim3, Dim7, k, o, l, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim7, Dim2, l, o, p, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim2, i, j, k, l, m, n, l, o, p, k, Dim0, Dim1, Dim4, Dim5, Dim3, Dim6, Dim7, Dim2, i, j, m, n, l, o, p, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(l, o, p, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim7, Dim2, l, o, p, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim7, Dim3, k, o, p, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim3, i, j, k, l, m, n, k, o, p, l, Dim0, Dim1, Dim4, Dim5, Dim2, Dim6, Dim7, Dim3, i, j, m, n, k, o, p, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(k, o, p, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim7, Dim3, k, o, p, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim2, Dim7, o, l, k, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim2, Dim7, i, j, k, l, m, n, o, l, k, p, Dim0, Dim1, Dim4, Dim5, Dim6, Dim3, Dim2, Dim7, i, j, m, n, o, l, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(o, l, k, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim2, Dim7, o, l, k, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim3, Dim7, o, k, l, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, i, j, k, l, m, n, o, k, l, p, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim3, Dim7, i, j, m, n, o, k, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(o, k, l, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim3, Dim7, o, k, l, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim7, Dim2, o, l, p, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim2, i, j, k, l, m, n, o, l, p, k, Dim0, Dim1, Dim4, Dim5, Dim6, Dim3, Dim7, Dim2, i, j, m, n, o, l, p, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(o, l, p, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim7, Dim2, o, l, p, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim7, Dim3, o, k, p, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim3, i, j, k, l, m, n, o, k, p, l, Dim0, Dim1, Dim4, Dim5, Dim6, Dim2, Dim7, Dim3, i, j, m, n, o, k, p, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(o, k, p, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim7, Dim3, o, k, p, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim3, Dim2, o, p, l, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim2, i, j, k, l, m, n, o, p, l, k, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim3, Dim2, i, j, m, n, o, p, l, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, l, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim3, Dim2, o, p, l, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim2, Dim3, o, p, k, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim3, i, j, k, l, m, n, o, p, k, l, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, Dim2, Dim3, i, j, m, n, o, p, k, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, Dim6, Dim7, i, j, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, k, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim2, Dim3, o, p, k, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim6, Dim7, m, k, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim6, Dim7, i, j, k, l, m, n, m, k, o, p, Dim0, Dim1, Dim3, Dim5, Dim4, Dim2, Dim6, Dim7, i, j, l, n, m, k, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(m, k, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim6, Dim7, m, k, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim6, Dim7, k, m, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim6, Dim7, i, j, k, l, m, n, k, m, o, p, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, Dim6, Dim7, i, j, l, n, k, m, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(k, m, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim6, Dim7, k, m, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim2, Dim7, m, o, k, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim2, Dim7, i, j, k, l, m, n, m, o, k, p, Dim0, Dim1, Dim3, Dim5, Dim4, Dim6, Dim2, Dim7, i, j, l, n, m, o, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(m, o, k, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim2, Dim7, m, o, k, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim4, Dim7, k, o, m, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim4, Dim7, i, j, k, l, m, n, k, o, m, p, Dim0, Dim1, Dim3, Dim5, Dim2, Dim6, Dim4, Dim7, i, j, l, n, k, o, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(k, o, m, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim4, Dim7, k, o, m, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim7, Dim2, m, o, p, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim2, i, j, k, l, m, n, m, o, p, k, Dim0, Dim1, Dim3, Dim5, Dim4, Dim6, Dim7, Dim2, i, j, l, n, m, o, p, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(m, o, p, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim7, Dim2, m, o, p, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim7, Dim4, k, o, p, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim4, i, j, k, l, m, n, k, o, p, m, Dim0, Dim1, Dim3, Dim5, Dim2, Dim6, Dim7, Dim4, i, j, l, n, k, o, p, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(k, o, p, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim7, Dim4, k, o, p, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim2, Dim7, o, m, k, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim2, Dim7, i, j, k, l, m, n, o, m, k, p, Dim0, Dim1, Dim3, Dim5, Dim6, Dim4, Dim2, Dim7, i, j, l, n, o, m, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(o, m, k, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim2, Dim7, o, m, k, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim4, Dim7, o, k, m, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim4, Dim7, i, j, k, l, m, n, o, k, m, p, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim4, Dim7, i, j, l, n, o, k, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(o, k, m, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim4, Dim7, o, k, m, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim7, Dim2, o, m, p, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim2, i, j, k, l, m, n, o, m, p, k, Dim0, Dim1, Dim3, Dim5, Dim6, Dim4, Dim7, Dim2, i, j, l, n, o, m, p, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(o, m, p, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim7, Dim2, o, m, p, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim7, Dim4, o, k, p, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim4, i, j, k, l, m, n, o, k, p, m, Dim0, Dim1, Dim3, Dim5, Dim6, Dim2, Dim7, Dim4, i, j, l, n, o, k, p, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(o, k, p, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim7, Dim4, o, k, p, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim4, Dim2, o, p, m, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim2, i, j, k, l, m, n, o, p, m, k, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim4, Dim2, i, j, l, n, o, p, m, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(o, p, m, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim4, Dim2, o, p, m, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim2, Dim4, o, p, k, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim4, i, j, k, l, m, n, o, p, k, m, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, Dim2, Dim4, i, j, l, n, o, p, k, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, Dim6, Dim7, i, j, l, n, o, p>(TensorExpr(a, b));
}

//B(o, p, k, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim2, Dim4, o, p, k, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim7, n, k, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, Dim6, Dim7, i, j, k, l, m, n, n, k, o, p, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, i, j, l, m, n, k, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(n, k, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim7, n, k, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim7, k, n, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, Dim6, Dim7, i, j, k, l, m, n, k, n, o, p, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, i, j, l, m, k, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(k, n, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim7, k, n, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim7, n, o, k, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim2, Dim7, i, j, k, l, m, n, n, o, k, p, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, i, j, l, m, n, o, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(n, o, k, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim7, n, o, k, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim5, Dim7, k, o, n, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim5, Dim7, i, j, k, l, m, n, k, o, n, p, Dim0, Dim1, Dim3, Dim4, Dim2, Dim6, Dim5, Dim7, i, j, l, m, k, o, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(k, o, n, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim5, Dim7, k, o, n, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim2, n, o, p, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim2, i, j, k, l, m, n, n, o, p, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, i, j, l, m, n, o, p, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(n, o, p, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim2, n, o, p, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim2, Dim6, Dim7, Dim5, k, o, p, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, Dim5, i, j, k, l, m, n, k, o, p, n, Dim0, Dim1, Dim3, Dim4, Dim2, Dim6, Dim7, Dim5, i, j, l, m, k, o, p, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(k, o, p, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim6, Dim7, Dim5, k, o, p, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim2, Dim7, o, n, k, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim2, Dim7, i, j, k, l, m, n, o, n, k, p, Dim0, Dim1, Dim3, Dim4, Dim6, Dim5, Dim2, Dim7, i, j, l, m, o, n, k, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(o, n, k, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim2, Dim7, o, n, k, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim5, Dim7, o, k, n, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim5, Dim7, i, j, k, l, m, n, o, k, n, p, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim5, Dim7, i, j, l, m, o, k, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(o, k, n, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim5, Dim7, o, k, n, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim7, Dim2, o, n, p, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim2, i, j, k, l, m, n, o, n, p, k, Dim0, Dim1, Dim3, Dim4, Dim6, Dim5, Dim7, Dim2, i, j, l, m, o, n, p, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(o, n, p, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim7, Dim2, o, n, p, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim2, Dim7, Dim5, o, k, p, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, Dim5, i, j, k, l, m, n, o, k, p, n, Dim0, Dim1, Dim3, Dim4, Dim6, Dim2, Dim7, Dim5, i, j, l, m, o, k, p, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(o, k, p, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim2, Dim7, Dim5, o, k, p, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim5, Dim2, o, p, n, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim2, i, j, k, l, m, n, o, p, n, k, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim5, Dim2, i, j, l, m, o, p, n, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(o, p, n, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim5, Dim2, o, p, n, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim2, Dim5, o, p, k, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, Dim5, i, j, k, l, m, n, o, p, k, n, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, Dim2, Dim5, i, j, l, m, o, p, k, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim6, Dim7, i, j, l, m, o, p>(TensorExpr(a, b));
}

//B(o, p, k, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim2, Dim5, o, p, k, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim6, Dim7, m, l, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim6, Dim7, i, j, k, l, m, n, m, l, o, p, Dim0, Dim1, Dim2, Dim5, Dim4, Dim3, Dim6, Dim7, i, j, k, n, m, l, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(m, l, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim6, Dim7, m, l, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim6, Dim7, l, m, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim6, Dim7, i, j, k, l, m, n, l, m, o, p, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, Dim6, Dim7, i, j, k, n, l, m, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(l, m, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim6, Dim7, l, m, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim3, Dim7, m, o, l, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim3, Dim7, i, j, k, l, m, n, m, o, l, p, Dim0, Dim1, Dim2, Dim5, Dim4, Dim6, Dim3, Dim7, i, j, k, n, m, o, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(m, o, l, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim3, Dim7, m, o, l, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim4, Dim7, l, o, m, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim4, Dim7, i, j, k, l, m, n, l, o, m, p, Dim0, Dim1, Dim2, Dim5, Dim3, Dim6, Dim4, Dim7, i, j, k, n, l, o, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(l, o, m, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim4, Dim7, l, o, m, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim7, Dim3, m, o, p, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim3, i, j, k, l, m, n, m, o, p, l, Dim0, Dim1, Dim2, Dim5, Dim4, Dim6, Dim7, Dim3, i, j, k, n, m, o, p, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(m, o, p, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim7, Dim3, m, o, p, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim7, Dim4, l, o, p, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim4, i, j, k, l, m, n, l, o, p, m, Dim0, Dim1, Dim2, Dim5, Dim3, Dim6, Dim7, Dim4, i, j, k, n, l, o, p, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(l, o, p, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim7, Dim4, l, o, p, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim3, Dim7, o, m, l, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim3, Dim7, i, j, k, l, m, n, o, m, l, p, Dim0, Dim1, Dim2, Dim5, Dim6, Dim4, Dim3, Dim7, i, j, k, n, o, m, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(o, m, l, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim3, Dim7, o, m, l, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim4, Dim7, o, l, m, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim4, Dim7, i, j, k, l, m, n, o, l, m, p, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim4, Dim7, i, j, k, n, o, l, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(o, l, m, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim4, Dim7, o, l, m, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim7, Dim3, o, m, p, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim3, i, j, k, l, m, n, o, m, p, l, Dim0, Dim1, Dim2, Dim5, Dim6, Dim4, Dim7, Dim3, i, j, k, n, o, m, p, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(o, m, p, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim7, Dim3, o, m, p, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim7, Dim4, o, l, p, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim4, i, j, k, l, m, n, o, l, p, m, Dim0, Dim1, Dim2, Dim5, Dim6, Dim3, Dim7, Dim4, i, j, k, n, o, l, p, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(o, l, p, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim7, Dim4, o, l, p, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim4, Dim3, o, p, m, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim3, i, j, k, l, m, n, o, p, m, l, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim4, Dim3, i, j, k, n, o, p, m, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(o, p, m, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim4, Dim3, o, p, m, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim3, Dim4, o, p, l, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim4, i, j, k, l, m, n, o, p, l, m, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, Dim3, Dim4, i, j, k, n, o, p, l, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, Dim6, Dim7, i, j, k, n, o, p>(TensorExpr(a, b));
}

//B(o, p, l, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim3, Dim4, o, p, l, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim7, n, l, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, Dim6, Dim7, i, j, k, l, m, n, n, l, o, p, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, i, j, k, m, n, l, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(n, l, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim7, n, l, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim7, l, n, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, Dim6, Dim7, i, j, k, l, m, n, l, n, o, p, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, i, j, k, m, l, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(l, n, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim7, l, n, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim7, n, o, l, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim3, Dim7, i, j, k, l, m, n, n, o, l, p, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, Dim7, i, j, k, m, n, o, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(n, o, l, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim7, n, o, l, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim5, Dim7, l, o, n, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim5, Dim7, i, j, k, l, m, n, l, o, n, p, Dim0, Dim1, Dim2, Dim4, Dim3, Dim6, Dim5, Dim7, i, j, k, m, l, o, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(l, o, n, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim5, Dim7, l, o, n, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim3, n, o, p, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim3, i, j, k, l, m, n, n, o, p, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim3, i, j, k, m, n, o, p, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(n, o, p, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim3, n, o, p, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim3, Dim6, Dim7, Dim5, l, o, p, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, Dim5, i, j, k, l, m, n, l, o, p, n, Dim0, Dim1, Dim2, Dim4, Dim3, Dim6, Dim7, Dim5, i, j, k, m, l, o, p, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(l, o, p, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim6, Dim7, Dim5, l, o, p, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim3, Dim7, o, n, l, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim3, Dim7, i, j, k, l, m, n, o, n, l, p, Dim0, Dim1, Dim2, Dim4, Dim6, Dim5, Dim3, Dim7, i, j, k, m, o, n, l, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(o, n, l, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim3, Dim7, o, n, l, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim5, Dim7, o, l, n, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim5, Dim7, i, j, k, l, m, n, o, l, n, p, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim5, Dim7, i, j, k, m, o, l, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(o, l, n, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim5, Dim7, o, l, n, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim7, Dim3, o, n, p, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim3, i, j, k, l, m, n, o, n, p, l, Dim0, Dim1, Dim2, Dim4, Dim6, Dim5, Dim7, Dim3, i, j, k, m, o, n, p, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(o, n, p, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim7, Dim3, o, n, p, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim3, Dim7, Dim5, o, l, p, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, Dim5, i, j, k, l, m, n, o, l, p, n, Dim0, Dim1, Dim2, Dim4, Dim6, Dim3, Dim7, Dim5, i, j, k, m, o, l, p, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(o, l, p, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim3, Dim7, Dim5, o, l, p, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim5, Dim3, o, p, n, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim3, i, j, k, l, m, n, o, p, n, l, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim5, Dim3, i, j, k, m, o, p, n, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(o, p, n, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim5, Dim3, o, p, n, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim3, Dim5, o, p, l, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, Dim5, i, j, k, l, m, n, o, p, l, n, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, Dim3, Dim5, i, j, k, m, o, p, l, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim6, Dim7, i, j, k, m, o, p>(TensorExpr(a, b));
}

//B(o, p, l, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim3, Dim5, o, p, l, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim7, n, m, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, Dim6, Dim7, i, j, k, l, m, n, n, m, o, p, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, i, j, k, l, n, m, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(n, m, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim7, n, m, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim7, m, n, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, m, n, o, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(m, n, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim7, m, n, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim7, n, o, m, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim4, Dim7, i, j, k, l, m, n, n, o, m, p, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, Dim7, i, j, k, l, n, o, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(n, o, m, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim7, n, o, m, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim5, Dim7, m, o, n, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim5, Dim7, i, j, k, l, m, n, m, o, n, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim5, Dim7, i, j, k, l, m, o, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(m, o, n, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim5, Dim7, m, o, n, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim4, n, o, p, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, Dim4, i, j, k, l, m, n, n, o, p, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim4, i, j, k, l, n, o, p, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(n, o, p, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim4, n, o, p, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim4, Dim6, Dim7, Dim5, m, o, p, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, Dim5, i, j, k, l, m, n, m, o, p, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim5, i, j, k, l, m, o, p, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(m, o, p, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim6, Dim7, Dim5, m, o, p, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim4, Dim7, o, n, m, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim4, Dim7, i, j, k, l, m, n, o, n, m, p, Dim0, Dim1, Dim2, Dim3, Dim6, Dim5, Dim4, Dim7, i, j, k, l, o, n, m, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(o, n, m, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim4, Dim7, o, n, m, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim5, Dim7, o, m, n, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim5, Dim7, i, j, k, l, m, n, o, m, n, p, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim5, Dim7, i, j, k, l, o, m, n, p>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(o, m, n, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim5, Dim7, o, m, n, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim5, Dim7, Dim4, o, n, p, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, Dim4, i, j, k, l, m, n, o, n, p, m, Dim0, Dim1, Dim2, Dim3, Dim6, Dim5, Dim7, Dim4, i, j, k, l, o, n, p, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(o, n, p, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim5, Dim7, Dim4, o, n, p, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim4, Dim7, Dim5, o, m, p, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, Dim5, i, j, k, l, m, n, o, m, p, n, Dim0, Dim1, Dim2, Dim3, Dim6, Dim4, Dim7, Dim5, i, j, k, l, o, m, p, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(o, m, p, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim4, Dim7, Dim5, o, m, p, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim5, Dim4, o, p, n, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, Dim4, i, j, k, l, m, n, o, p, n, m, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim5, Dim4, i, j, k, l, o, p, n, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(o, p, n, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim5, Dim4, o, p, n, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor4_Expr<B, U, Dim6, Dim7, Dim4, Dim5, o, p, m, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, Dim5, i, j, k, l, m, n, o, p, m, n, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, Dim4, Dim5, i, j, k, l, o, p, m, n>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim6, Dim7, i, j, k, l, o, p>(TensorExpr(a, b));
}

//B(o, p, m, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim6, Dim7, Dim4, Dim5, o, p, m, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

} //namespace FTensor

#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, n, )*B(o, p, q, ) single contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q, int DimA, int DimB, int DimC, int DimD, int DimE, int DimX, int DimY, int DimZ, char a, char b, char c, char d, char e, char x, char y, char z>
class Tensor6_times_Tensor3_single
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor3_Expr<B, U, Dim6, Dim7, Dim8, o, p, q> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const int N7) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX, ++xx)
      {
        // Permutation is where the indices get checked.
        result += Permutation6<DimA, DimB, DimC, DimD, DimE, DimX, a, b, c, d, e, x>().eval(
          iterA, N1, N2, N3, N4, N5, xx) * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx, N6, N7);
      }
    return result;
  }

  Tensor6_times_Tensor3_single(
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &iter_a,
    const Tensor3_Expr<B, U, Dim6, Dim7, Dim8, o, p, q> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, n, )*B(i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim0, Dim6, Dim7, i, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, i, j, k, l, m, n, i, o, p, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, j, k, l, m, n, i, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, j, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim6, Dim7, i, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim0, Dim7, o, i, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, i, j, k, l, m, n, o, i, p, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, j, k, l, m, n, o, i, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, j, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, i, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim0, Dim7, o, i, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim7, Dim0, o, p, i> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, i, j, k, l, m, n, o, p, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, j, k, l, m, n, o, p, i>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, j, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, i, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim7, Dim0, o, p, i> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim1, Dim6, Dim7, j, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, i, j, k, l, m, n, j, o, p, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, i, k, l, m, n, j, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim6, Dim7, j, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim1, Dim7, o, j, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, i, j, k, l, m, n, o, j, p, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, i, k, l, m, n, o, j, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, j, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim1, Dim7, o, j, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim7, Dim1, o, p, j> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, i, j, k, l, m, n, o, p, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, i, k, l, m, n, o, p, j>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, j, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim7, Dim1, o, p, j> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim2, Dim6, Dim7, k, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, i, j, k, l, m, n, k, o, p, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, i, j, l, m, n, k, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim6, Dim7, k, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim2, Dim7, o, k, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, i, j, k, l, m, n, o, k, p, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, i, j, l, m, n, o, k, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, k, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim2, Dim7, o, k, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim7, Dim2, o, p, k> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, i, j, k, l, m, n, o, p, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, i, j, l, m, n, o, p, k>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, k, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim7, Dim2, o, p, k> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim3, Dim6, Dim7, l, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, i, j, k, l, m, n, l, o, p, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, i, j, k, m, n, l, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim6, Dim7, l, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim3, Dim7, o, l, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, i, j, k, l, m, n, o, l, p, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, Dim7, i, j, k, m, n, o, l, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, l, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim3, Dim7, o, l, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim7, Dim3, o, p, l> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, i, j, k, l, m, n, o, p, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim3, i, j, k, m, n, o, p, l>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(o, p, l, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim7, Dim3, o, p, l> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim4, Dim6, Dim7, m, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, i, j, k, l, m, n, m, o, p, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, i, j, k, l, n, m, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim6, Dim7, m, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim4, Dim7, o, m, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, i, j, k, l, m, n, o, m, p, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, Dim7, i, j, k, l, n, o, m, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, m, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim4, Dim7, o, m, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim7, Dim4, o, p, m> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, i, j, k, l, m, n, o, p, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim4, i, j, k, l, n, o, p, m>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(o, p, m, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim7, Dim4, o, p, m> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim7, n, o, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim6, Dim7, i, j, k, l, m, n, n, o, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, i, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(n, o, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim7, n, o, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim5, Dim7, o, n, p> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, i, j, k, l, m, n, o, n, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim5, Dim7, i, j, k, l, m, o, n, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, i, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, n, p, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim5, Dim7, o, n, p> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, )*B(o, p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
               const Tensor3_Expr<B, U, Dim6, Dim7, Dim5, o, p, n> &b)
{
  using TensorExpr
    = Tensor6_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, i, j, k, l, m, n, o, p, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, Dim5, i, j, k, l, m, o, p, n>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, i, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(o, p, n, )*A(i, j, k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor3_Expr<B, U, Dim6, Dim7, Dim5, o, p, n> &b,
               const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a)
{
  return a * b;
}

} //namespace FTensor

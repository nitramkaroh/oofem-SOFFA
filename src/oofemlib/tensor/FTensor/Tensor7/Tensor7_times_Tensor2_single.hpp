#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, n, o, )*B(p, q, ) single contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q, int DimA, int DimB, int DimC, int DimD, int DimE, int DimF, int DimX, int DimY, char a, char b, char c, char d, char e, char f, char x, char y>
class Tensor7_times_Tensor2_single
{
  Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
  Tensor2_Expr<B, U, Dim7, Dim8, p, q> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const int N7) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      {
        // Permutation is where the indices get checked.
        result += Permutation7<DimA, DimB, DimC, DimD, DimE, DimF, DimX, a, b, c, d, e, f, x>().eval(
          iterA, N1, N2, N3, N4, N5, N6, xx) * Permutation2<DimX, DimY, x, y>().eval(iterB, xx, N7);
      }
    return result;
  }

  Tensor7_times_Tensor2_single(
    const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &iter_a,
    const Tensor2_Expr<B, U, Dim7, Dim8, p, q> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, n, o, )*B(i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim0, Dim7, i, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, i, j, k, l, m, n, o, i, p, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, j, k, l, m, n, o, i, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, j, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim0, Dim7, i, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim7, Dim0, p, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, i, j, k, l, m, n, o, p, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, j, k, l, m, n, o, i, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, j, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(p, i, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim7, Dim0, p, i> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim1, Dim7, j, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, i, j, k, l, m, n, o, j, p, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, i, k, l, m, n, o, j, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim1, Dim7, j, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim7, Dim1, p, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, i, j, k, l, m, n, o, p, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, i, k, l, m, n, o, j, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(p, j, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim7, Dim1, p, j> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim2, Dim7, k, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, i, j, k, l, m, n, o, k, p, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, i, j, l, m, n, o, k, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim2, Dim7, k, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim7, Dim2, p, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, i, j, k, l, m, n, o, p, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, i, j, l, m, n, o, k, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(p, k, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim7, Dim2, p, k> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim3, Dim7, l, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, i, j, k, l, m, n, o, l, p, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, Dim7, i, j, k, m, n, o, l, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim3, Dim7, l, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim7, Dim3, p, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, i, j, k, l, m, n, o, p, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, Dim7, i, j, k, m, n, o, l, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(p, l, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim7, Dim3, p, l> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim4, Dim7, m, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, i, j, k, l, m, n, o, m, p, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, Dim7, i, j, k, l, n, o, m, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim4, Dim7, m, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim7, Dim4, p, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, i, j, k, l, m, n, o, p, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, Dim7, i, j, k, l, n, o, m, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(p, m, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim7, Dim4, p, m> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(n, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim5, Dim7, n, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim5, Dim7, i, j, k, l, m, n, o, n, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim5, Dim7, i, j, k, l, m, o, n, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, i, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(n, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim5, Dim7, n, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim7, Dim5, p, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim5, i, j, k, l, m, n, o, p, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim5, Dim7, i, j, k, l, m, o, n, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, Dim7, i, j, k, l, m, o, p>(TensorExpr(a, b));
}

//B(p, n, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim7, Dim5, p, n> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim6, Dim7, o, p> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim6, Dim7, i, j, k, l, m, n, o, o, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, i, j, k, l, m, n, p>(TensorExpr(a, b));
}

//B(o, p, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim6, Dim7, o, p> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

//A(i, j, k, l, m, n, o, )*B(p, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
               const Tensor2_Expr<B, U, Dim7, Dim6, p, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim6, i, j, k, l, m, n, o, p, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim7, i, j, k, l, m, n, p>(TensorExpr(a, b));
}

//B(p, o, )*A(i, j, k, l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor2_Expr<B, U, Dim7, Dim6, p, o> &b,
               const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  return a * b;
}

} //namespace FTensor

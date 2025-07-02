#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, )*B(n, o, p, q, ) single contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q, int DimA, int DimB, int DimC, int DimD, int DimX, int DimY, int DimZ, int DimW, char a, char b, char c, char d, char x, char y, char z, char w>
class Tensor5_times_Tensor4_single
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim8, n, o, p, q> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const int N7) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      {
        // Permutation is where the indices get checked.
        result += Permutation5<DimA, DimB, DimC, DimD, DimX, a, b, c, d, x>().eval(
          iterA, N1, N2, N3, N4, xx) * Permutation4<DimX, DimY, DimZ, DimW, x, y, z, w>().eval(iterB, xx, N5, N6, N7);
      }
    return result;
  }

  Tensor5_times_Tensor4_single(
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &iter_a,
    const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim8, n, o, p, q> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, )*B(i, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim7, i, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, i, j, k, l, m, i, n, o, p, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, Dim7, j, k, l, m, i, n, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, j, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(i, n, o, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim6, Dim7, i, n, o, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim7, n, i, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, i, j, k, l, m, n, i, o, p, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, Dim7, j, k, l, m, n, i, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, j, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(n, i, o, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim6, Dim7, n, i, o, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, i, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim7, n, o, i, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, i, j, k, l, m, n, o, i, p, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, Dim7, j, k, l, m, n, o, i, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, j, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(n, o, i, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim0, Dim7, n, o, i, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, p, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim0, n, o, p, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, i, j, k, l, m, n, o, p, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim0, j, k, l, m, n, o, p, i>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, j, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(n, o, p, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim0, n, o, p, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim7, j, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, i, j, k, l, m, j, n, o, p, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, Dim7, i, k, l, m, j, n, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(j, n, o, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim6, Dim7, j, n, o, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim7, n, j, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, i, j, k, l, m, n, j, o, p, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, Dim7, i, k, l, m, n, j, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(n, j, o, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim6, Dim7, n, j, o, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, j, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim7, n, o, j, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, i, j, k, l, m, n, o, j, p, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, Dim7, i, k, l, m, n, o, j, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(n, o, j, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim1, Dim7, n, o, j, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, p, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim1, n, o, p, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, i, j, k, l, m, n, o, p, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim1, i, k, l, m, n, o, p, j>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, k, l, m, n, o, p>(TensorExpr(a, b));
}

//B(n, o, p, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim1, n, o, p, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim7, k, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, i, j, k, l, m, k, n, o, p, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, Dim7, i, j, l, m, k, n, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(k, n, o, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim6, Dim7, k, n, o, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim7, n, k, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, i, j, k, l, m, n, k, o, p, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, Dim7, i, j, l, m, n, k, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(n, k, o, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim6, Dim7, n, k, o, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, k, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim7, n, o, k, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, i, j, k, l, m, n, o, k, p, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, Dim7, i, j, l, m, n, o, k, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(n, o, k, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim2, Dim7, n, o, k, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, p, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim2, n, o, p, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, i, j, k, l, m, n, o, p, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, Dim2, i, j, l, m, n, o, p, k>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, l, m, n, o, p>(TensorExpr(a, b));
}

//B(n, o, p, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim2, n, o, p, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim7, l, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, Dim7, i, j, k, l, m, l, n, o, p, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, Dim7, i, j, k, m, l, n, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(l, n, o, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim6, Dim7, l, n, o, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim7, n, l, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, Dim7, i, j, k, l, m, n, l, o, p, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, Dim7, i, j, k, m, n, l, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(n, l, o, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim6, Dim7, n, l, o, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, l, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim7, n, o, l, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, Dim7, i, j, k, l, m, n, o, l, p, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, Dim7, i, j, k, m, n, o, l, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(n, o, l, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim3, Dim7, n, o, l, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, p, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim3, n, o, p, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim3, i, j, k, l, m, n, o, p, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, Dim3, i, j, k, m, n, o, p, l>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim7, i, j, k, m, n, o, p>(TensorExpr(a, b));
}

//B(n, o, p, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim3, n, o, p, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim7, m, n, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, m, n, o, p, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, i, j, k, l, m, n, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(m, n, o, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim7, m, n, o, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, o, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim7, n, m, o, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, Dim7, i, j, k, l, m, n, m, o, p, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, Dim7, i, j, k, l, n, m, o, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(n, m, o, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim6, Dim7, n, m, o, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, m, p, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim7, n, o, m, p> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, Dim7, i, j, k, l, m, n, o, m, p, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, Dim7, i, j, k, l, n, o, m, p>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(n, o, m, p, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim4, Dim7, n, o, m, p> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, p, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim4, n, o, p, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim7, Dim4, i, j, k, l, m, n, o, p, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, Dim4, i, j, k, l, n, o, p, m>;
  return Tensor7_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim7, i, j, k, l, n, o, p>(TensorExpr(a, b));
}

//B(n, o, p, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim4, n, o, p, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

} //namespace FTensor

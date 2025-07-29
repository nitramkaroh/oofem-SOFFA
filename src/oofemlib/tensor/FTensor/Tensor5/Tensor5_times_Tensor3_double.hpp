#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, )*B(n, o, p, ) double contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p, int DimA, int DimB, int DimC, int DimX, int DimY, int DimZ, char a, char b, char c, char x, char y, char z>
class Tensor5_times_Tensor3_double
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor3_Expr<B, U, Dim5, Dim6, Dim7, n, o, p> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        {
          // Permutation is where the indices get checked.
          result += Permutation5<DimA, DimB, DimC, DimX, DimY, a, b, c, x, y>().eval(
            iterA, N1, N2, N3, xx, yy) * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx, yy, N4);
        }
    return result;
  }

  Tensor5_times_Tensor3_double(
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &iter_a,
    const Tensor3_Expr<B, U, Dim5, Dim6, Dim7, n, o, p> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, )*B(j, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim5, j, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, i, j, k, l, m, j, i, n, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//B(j, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim5, j, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim5, i, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, i, j, k, l, m, i, j, n, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//B(i, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim5, i, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim0, j, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, i, j, k, l, m, j, n, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//B(j, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim0, j, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim1, i, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, i, j, k, l, m, i, n, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//B(i, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim1, i, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim0, n, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, i, j, k, l, m, n, j, i, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//B(n, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim0, n, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim1, n, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, i, j, k, l, m, n, i, j, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//B(n, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim1, n, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim5, k, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, i, j, k, l, m, k, i, n, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//B(k, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim5, k, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim5, i, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, i, j, k, l, m, i, k, n, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//B(i, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim5, i, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim0, k, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, i, j, k, l, m, k, n, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//B(k, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim0, k, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim2, i, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, i, j, k, l, m, i, n, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//B(i, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim2, i, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim0, n, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, i, j, k, l, m, n, k, i, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//B(n, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim0, n, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim2, n, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, i, j, k, l, m, n, i, k, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//B(n, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim2, n, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim5, l, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, i, j, k, l, m, l, i, n, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//B(l, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim5, l, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim5, i, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, i, j, k, l, m, i, l, n, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//B(i, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim5, i, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim0, l, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, i, j, k, l, m, l, n, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//B(l, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim0, l, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim3, i, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, i, j, k, l, m, i, n, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//B(i, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim3, i, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim0, n, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, i, j, k, l, m, n, l, i, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//B(n, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim0, n, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim3, n, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, i, j, k, l, m, n, i, l, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//B(n, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim3, n, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim5, m, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim5, i, j, k, l, m, m, i, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>(TensorExpr(a, b));
}

//B(m, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim5, m, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim5, i, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim5, i, j, k, l, m, i, m, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>(TensorExpr(a, b));
}

//B(i, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim5, i, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim0, m, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim0, i, j, k, l, m, m, n, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>(TensorExpr(a, b));
}

//B(m, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim0, m, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim4, i, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim4, i, j, k, l, m, i, n, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>(TensorExpr(a, b));
}

//B(i, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim4, i, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim0, n, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, i, j, k, l, m, n, m, i, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>(TensorExpr(a, b));
}

//B(n, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim0, n, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim4, n, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, i, j, k, l, m, n, i, m, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>(TensorExpr(a, b));
}

//B(n, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim4, n, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim5, k, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, i, j, k, l, m, k, j, n, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//B(k, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim5, k, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim5, j, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, i, j, k, l, m, j, k, n, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//B(j, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim5, j, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim1, k, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, i, j, k, l, m, k, n, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//B(k, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim1, k, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim2, j, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, i, j, k, l, m, j, n, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//B(j, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim2, j, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim1, n, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, i, j, k, l, m, n, k, j, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//B(n, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim1, n, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim2, n, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, i, j, k, l, m, n, j, k, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//B(n, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim2, n, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim5, l, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, i, j, k, l, m, l, j, n, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//B(l, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim5, l, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim5, j, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, i, j, k, l, m, j, l, n, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//B(j, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim5, j, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim1, l, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, i, j, k, l, m, l, n, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//B(l, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim1, l, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim3, j, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, i, j, k, l, m, j, n, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//B(j, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim3, j, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim1, n, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, i, j, k, l, m, n, l, j, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//B(n, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim1, n, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim3, n, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, i, j, k, l, m, n, j, l, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//B(n, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim3, n, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim5, m, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim5, i, j, k, l, m, m, j, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>(TensorExpr(a, b));
}

//B(m, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim5, m, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim5, j, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim5, i, j, k, l, m, j, m, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>(TensorExpr(a, b));
}

//B(j, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim5, j, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim1, m, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim1, i, j, k, l, m, m, n, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>(TensorExpr(a, b));
}

//B(m, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim1, m, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim4, j, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim4, i, j, k, l, m, j, n, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>(TensorExpr(a, b));
}

//B(j, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim4, j, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim1, n, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, i, j, k, l, m, n, m, j, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>(TensorExpr(a, b));
}

//B(n, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim1, n, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim4, n, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, i, j, k, l, m, n, j, m, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>(TensorExpr(a, b));
}

//B(n, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim4, n, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim5, l, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, i, j, k, l, m, l, k, n, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//B(l, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim5, l, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim5, k, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, i, j, k, l, m, k, l, n, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//B(k, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim5, k, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim2, l, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, i, j, k, l, m, l, n, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//B(l, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim2, l, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim3, k, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, i, j, k, l, m, k, n, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//B(k, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim3, k, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim2, n, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, i, j, k, l, m, n, l, k, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//B(n, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim2, n, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim3, n, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, i, j, k, l, m, n, k, l, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//B(n, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim3, n, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim5, m, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim5, i, j, k, l, m, m, k, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>(TensorExpr(a, b));
}

//B(m, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim5, m, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim5, k, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim5, i, j, k, l, m, k, m, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>(TensorExpr(a, b));
}

//B(k, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim5, k, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim2, m, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim2, i, j, k, l, m, m, n, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>(TensorExpr(a, b));
}

//B(m, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim2, m, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim4, k, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim4, i, j, k, l, m, k, n, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>(TensorExpr(a, b));
}

//B(k, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim4, k, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim2, n, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, i, j, k, l, m, n, m, k, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>(TensorExpr(a, b));
}

//B(n, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim2, n, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim4, n, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, i, j, k, l, m, n, k, m, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>(TensorExpr(a, b));
}

//B(n, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim4, n, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim5, m, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim5, i, j, k, l, m, m, l, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>(TensorExpr(a, b));
}

//B(m, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim5, m, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim5, i, j, k, l, m, l, m, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>(TensorExpr(a, b));
}

//B(l, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim3, m, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim3, i, j, k, l, m, m, n, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>(TensorExpr(a, b));
}

//B(m, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim3, m, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim4, l, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim4, i, j, k, l, m, l, n, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>(TensorExpr(a, b));
}

//B(l, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim4, l, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim3, n, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, i, j, k, l, m, n, m, l, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>(TensorExpr(a, b));
}

//B(n, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim3, n, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim4, n, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, i, j, k, l, m, n, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>(TensorExpr(a, b));
}

//B(n, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim4, n, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

} //namespace FTensor

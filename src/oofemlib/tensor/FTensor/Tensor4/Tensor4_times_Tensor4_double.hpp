#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, )*B(m, n, o, p, ) double contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p, int DimA, int DimB, int DimX, int DimY, int DimZ, int DimW, char a, char b, char x, char y, char z, char w>
class Tensor4_times_Tensor4_double
{
  Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
  Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim7, m, n, o, p> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX, ++xx)
      for(int yy = 0; yy < DimY, ++yy)
        {
          // Permutation is where the indices get checked.
          result += Permutation4<DimA, DimB, DimX, DimY, a, b, x, y>().eval(
            iterA, N1, N2, xx, yy) * Permutation4<DimX, DimY, DimZ, DimW, x, y, z, w>().eval(iterB, xx, yy, N3, N4);
        }
    return result;
  }

  Tensor4_times_Tensor4_double(
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &iter_a,
    const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim7, m, n, o, p> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, )*B(j, i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim5, j, i, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim0, Dim4, Dim5, i, j, k, l, j, i, m, n, Dim2, Dim3, Dim1, Dim0, Dim4, Dim5, k, l, j, i, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim5, i, j, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, i, j, k, l, i, j, m, n, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim5, j, m, i, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim0, Dim5, i, j, k, l, j, m, i, n, Dim2, Dim3, Dim1, Dim4, Dim0, Dim5, k, l, j, m, i, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim5, i, m, j, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim1, Dim5, i, j, k, l, i, m, j, n, Dim2, Dim3, Dim0, Dim4, Dim1, Dim5, k, l, i, m, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim0, j, m, n, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim5, Dim0, i, j, k, l, j, m, n, i, Dim2, Dim3, Dim1, Dim4, Dim5, Dim0, k, l, j, m, n, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim1, i, m, n, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim1, i, j, k, l, i, m, n, j, Dim2, Dim3, Dim0, Dim4, Dim5, Dim1, k, l, i, m, n, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim5, m, j, i, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, i, j, k, l, m, j, i, n, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, k, l, m, j, i, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim5, m, i, j, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, i, j, k, l, m, i, j, n, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim0, m, j, n, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, i, j, k, l, m, j, n, i, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, k, l, m, j, n, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim1, m, i, n, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, i, j, k, l, m, i, n, j, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, k, l, m, i, n, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim0, m, n, j, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, i, j, k, l, m, n, j, i, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, k, l, m, n, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim1, m, n, i, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, i, j, k, l, m, n, i, j, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, k, l, m, n, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim5, k, i, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim0, Dim4, Dim5, i, j, k, l, k, i, m, n, Dim1, Dim3, Dim2, Dim0, Dim4, Dim5, j, l, k, i, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim5, i, k, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim2, Dim4, Dim5, i, j, k, l, i, k, m, n, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim5, k, m, i, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim0, Dim5, i, j, k, l, k, m, i, n, Dim1, Dim3, Dim2, Dim4, Dim0, Dim5, j, l, k, m, i, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim5, i, m, k, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim2, Dim5, i, j, k, l, i, m, k, n, Dim1, Dim3, Dim0, Dim4, Dim2, Dim5, j, l, i, m, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim0, k, m, n, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim5, Dim0, i, j, k, l, k, m, n, i, Dim1, Dim3, Dim2, Dim4, Dim5, Dim0, j, l, k, m, n, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim2, i, m, n, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim2, i, j, k, l, i, m, n, k, Dim1, Dim3, Dim0, Dim4, Dim5, Dim2, j, l, i, m, n, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim5, m, k, i, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, i, j, k, l, m, k, i, n, Dim1, Dim3, Dim4, Dim2, Dim0, Dim5, j, l, m, k, i, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim5, m, i, k, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, i, j, k, l, m, i, k, n, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim0, m, k, n, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, i, j, k, l, m, k, n, i, Dim1, Dim3, Dim4, Dim2, Dim5, Dim0, j, l, m, k, n, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim2, m, i, n, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, i, j, k, l, m, i, n, k, Dim1, Dim3, Dim4, Dim0, Dim5, Dim2, j, l, m, i, n, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim0, m, n, k, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, i, j, k, l, m, n, k, i, Dim1, Dim3, Dim4, Dim5, Dim2, Dim0, j, l, m, n, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim2, m, n, i, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, i, j, k, l, m, n, i, k, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, j, l, m, n, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim5, l, i, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim0, Dim4, Dim5, i, j, k, l, l, i, m, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim5, i, l, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim3, Dim4, Dim5, i, j, k, l, i, l, m, n, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim5, l, m, i, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim0, Dim5, i, j, k, l, l, m, i, n, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, j, k, l, m, i, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim5, i, m, l, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim3, Dim5, i, j, k, l, i, m, l, n, Dim1, Dim2, Dim0, Dim4, Dim3, Dim5, j, k, i, m, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim0, l, m, n, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim5, Dim0, i, j, k, l, l, m, n, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, j, k, l, m, n, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim3, i, m, n, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim3, i, j, k, l, i, m, n, l, Dim1, Dim2, Dim0, Dim4, Dim5, Dim3, j, k, i, m, n, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim5, m, l, i, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, i, j, k, l, m, l, i, n, Dim1, Dim2, Dim4, Dim3, Dim0, Dim5, j, k, m, l, i, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim5, m, i, l, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, i, j, k, l, m, i, l, n, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim0, m, l, n, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, i, j, k, l, m, l, n, i, Dim1, Dim2, Dim4, Dim3, Dim5, Dim0, j, k, m, l, n, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim3, m, i, n, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, i, j, k, l, m, i, n, l, Dim1, Dim2, Dim4, Dim0, Dim5, Dim3, j, k, m, i, n, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim0, m, n, l, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, i, j, k, l, m, n, l, i, Dim1, Dim2, Dim4, Dim5, Dim3, Dim0, j, k, m, n, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim3, m, n, i, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, i, j, k, l, m, n, i, l, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, j, k, m, n, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim5, k, j, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim1, Dim4, Dim5, i, j, k, l, k, j, m, n, Dim0, Dim3, Dim2, Dim1, Dim4, Dim5, i, l, k, j, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim5, j, k, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim2, Dim4, Dim5, i, j, k, l, j, k, m, n, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim5, k, m, j, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim1, Dim5, i, j, k, l, k, m, j, n, Dim0, Dim3, Dim2, Dim4, Dim1, Dim5, i, l, k, m, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim5, j, m, k, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim2, Dim5, i, j, k, l, j, m, k, n, Dim0, Dim3, Dim1, Dim4, Dim2, Dim5, i, l, j, m, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim1, k, m, n, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim5, Dim1, i, j, k, l, k, m, n, j, Dim0, Dim3, Dim2, Dim4, Dim5, Dim1, i, l, k, m, n, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim2, j, m, n, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim5, Dim2, i, j, k, l, j, m, n, k, Dim0, Dim3, Dim1, Dim4, Dim5, Dim2, i, l, j, m, n, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim5, m, k, j, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, i, j, k, l, m, k, j, n, Dim0, Dim3, Dim4, Dim2, Dim1, Dim5, i, l, m, k, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim5, m, j, k, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, i, j, k, l, m, j, k, n, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim1, m, k, n, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, i, j, k, l, m, k, n, j, Dim0, Dim3, Dim4, Dim2, Dim5, Dim1, i, l, m, k, n, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim2, m, j, n, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, i, j, k, l, m, j, n, k, Dim0, Dim3, Dim4, Dim1, Dim5, Dim2, i, l, m, j, n, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim1, m, n, k, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, i, j, k, l, m, n, k, j, Dim0, Dim3, Dim4, Dim5, Dim2, Dim1, i, l, m, n, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim2, m, n, j, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, i, j, k, l, m, n, j, k, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, i, l, m, n, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim5, l, j, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim1, Dim4, Dim5, i, j, k, l, l, j, m, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim5, j, l, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim3, Dim4, Dim5, i, j, k, l, j, l, m, n, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim5, l, m, j, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim1, Dim5, i, j, k, l, l, m, j, n, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, i, k, l, m, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim5, j, m, l, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim3, Dim5, i, j, k, l, j, m, l, n, Dim0, Dim2, Dim1, Dim4, Dim3, Dim5, i, k, j, m, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim1, l, m, n, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim5, Dim1, i, j, k, l, l, m, n, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, i, k, l, m, n, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim3, j, m, n, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim5, Dim3, i, j, k, l, j, m, n, l, Dim0, Dim2, Dim1, Dim4, Dim5, Dim3, i, k, j, m, n, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim5, m, l, j, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, i, j, k, l, m, l, j, n, Dim0, Dim2, Dim4, Dim3, Dim1, Dim5, i, k, m, l, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim5, m, j, l, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, i, j, k, l, m, j, l, n, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim1, m, l, n, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, i, j, k, l, m, l, n, j, Dim0, Dim2, Dim4, Dim3, Dim5, Dim1, i, k, m, l, n, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim3, m, j, n, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, i, j, k, l, m, j, n, l, Dim0, Dim2, Dim4, Dim1, Dim5, Dim3, i, k, m, j, n, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim1, m, n, l, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, i, j, k, l, m, n, l, j, Dim0, Dim2, Dim4, Dim5, Dim3, Dim1, i, k, m, n, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim3, m, n, j, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, i, j, k, l, m, n, j, l, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, i, k, m, n, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim5, l, k, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim2, Dim4, Dim5, i, j, k, l, l, k, m, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim5, k, l, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim3, Dim4, Dim5, i, j, k, l, k, l, m, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim5, l, m, k, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim2, Dim5, i, j, k, l, l, m, k, n, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, i, j, l, m, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim5, k, m, l, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim3, Dim5, i, j, k, l, k, m, l, n, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, i, j, k, m, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim2, l, m, n, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim5, Dim2, i, j, k, l, l, m, n, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, i, j, l, m, n, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim3, k, m, n, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim5, Dim3, i, j, k, l, k, m, n, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, i, j, k, m, n, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim5, m, l, k, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, i, j, k, l, m, l, k, n, Dim0, Dim1, Dim4, Dim3, Dim2, Dim5, i, j, m, l, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim5, m, k, l, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, i, j, k, l, m, k, l, n, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim2, m, l, n, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, i, j, k, l, m, l, n, k, Dim0, Dim1, Dim4, Dim3, Dim5, Dim2, i, j, m, l, n, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim3, m, k, n, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, i, j, k, l, m, k, n, l, Dim0, Dim1, Dim4, Dim2, Dim5, Dim3, i, j, m, k, n, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim2, m, n, l, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, i, j, k, l, m, n, l, k, Dim0, Dim1, Dim4, Dim5, Dim3, Dim2, i, j, m, n, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim3, m, n, k, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, i, j, k, l, m, n, k, l, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, i, j, m, n, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>(TensorExpr(a, b));
}

} //namespace FTensor

#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, )*B(m, n, o, ) double contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o, int DimA, int DimB, int DimX, int DimY, int DimZ, char a, char b, char x, char y, char z>
class Tensor4_times_Tensor3_double
{
  Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
  Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      for(int yy = 0; yy < DimY; ++yy)
        {
          // Permutation is where the indices get checked.
          result += Permutation4<DimA, DimB, DimX, DimY, a, b, x, y>().eval(
            iterA, N1, N2, xx, yy) * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx, yy, N3);
        }
    return result;
  }

  Tensor4_times_Tensor3_double(
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &iter_a,
    const Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, )*B(j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim4, j, i, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim0, Dim4, i, j, k, l, j, i, m, Dim2, Dim3, Dim1, Dim0, Dim4, k, l, j, i, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, k, l, m>(TensorExpr(a, b));
}

//B(j, i, m, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim4, j, i, m> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim4, i, j, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim1, Dim4, i, j, k, l, i, j, m, Dim2, Dim3, Dim0, Dim1, Dim4, k, l, i, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, k, l, m>(TensorExpr(a, b));
}

//B(i, j, m, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim4, i, j, m> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim0, j, m, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim0, i, j, k, l, j, m, i, Dim2, Dim3, Dim1, Dim4, Dim0, k, l, j, m, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, k, l, m>(TensorExpr(a, b));
}

//B(j, m, i, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim0, j, m, i> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim1, i, m, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim1, i, j, k, l, i, m, j, Dim2, Dim3, Dim0, Dim4, Dim1, k, l, i, m, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, k, l, m>(TensorExpr(a, b));
}

//B(i, m, j, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim1, i, m, j> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim0, m, j, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, i, j, k, l, m, j, i, Dim2, Dim3, Dim4, Dim1, Dim0, k, l, m, j, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, k, l, m>(TensorExpr(a, b));
}

//B(m, j, i, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim0, m, j, i> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim1, m, i, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, i, j, k, l, m, i, j, Dim2, Dim3, Dim4, Dim0, Dim1, k, l, m, i, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim4, k, l, m>(TensorExpr(a, b));
}

//B(m, i, j, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim1, m, i, j> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim4, k, i, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim0, Dim4, i, j, k, l, k, i, m, Dim1, Dim3, Dim2, Dim0, Dim4, j, l, k, i, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, j, l, m>(TensorExpr(a, b));
}

//B(k, i, m, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim4, k, i, m> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim4, i, k, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim2, Dim4, i, j, k, l, i, k, m, Dim1, Dim3, Dim0, Dim2, Dim4, j, l, i, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, j, l, m>(TensorExpr(a, b));
}

//B(i, k, m, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim4, i, k, m> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim0, k, m, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim0, i, j, k, l, k, m, i, Dim1, Dim3, Dim2, Dim4, Dim0, j, l, k, m, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, j, l, m>(TensorExpr(a, b));
}

//B(k, m, i, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim0, k, m, i> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim2, i, m, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim2, i, j, k, l, i, m, k, Dim1, Dim3, Dim0, Dim4, Dim2, j, l, i, m, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, j, l, m>(TensorExpr(a, b));
}

//B(i, m, k, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim2, i, m, k> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim0, m, k, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, i, j, k, l, m, k, i, Dim1, Dim3, Dim4, Dim2, Dim0, j, l, m, k, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, j, l, m>(TensorExpr(a, b));
}

//B(m, k, i, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim0, m, k, i> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim2, m, i, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, i, j, k, l, m, i, k, Dim1, Dim3, Dim4, Dim0, Dim2, j, l, m, i, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim4, j, l, m>(TensorExpr(a, b));
}

//B(m, i, k, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim2, m, i, k> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim4, l, i, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim0, Dim4, i, j, k, l, l, i, m, Dim1, Dim2, Dim3, Dim0, Dim4, j, k, l, i, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, j, k, m>(TensorExpr(a, b));
}

//B(l, i, m, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim4, l, i, m> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim4, i, l, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim3, Dim4, i, j, k, l, i, l, m, Dim1, Dim2, Dim0, Dim3, Dim4, j, k, i, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, j, k, m>(TensorExpr(a, b));
}

//B(i, l, m, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim4, i, l, m> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim0, l, m, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim0, i, j, k, l, l, m, i, Dim1, Dim2, Dim3, Dim4, Dim0, j, k, l, m, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, j, k, m>(TensorExpr(a, b));
}

//B(l, m, i, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim0, l, m, i> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim3, i, m, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim3, i, j, k, l, i, m, l, Dim1, Dim2, Dim0, Dim4, Dim3, j, k, i, m, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, j, k, m>(TensorExpr(a, b));
}

//B(i, m, l, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim3, i, m, l> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim0, m, l, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, i, j, k, l, m, l, i, Dim1, Dim2, Dim4, Dim3, Dim0, j, k, m, l, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, j, k, m>(TensorExpr(a, b));
}

//B(m, l, i, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim0, m, l, i> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim3, m, i, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, i, j, k, l, m, i, l, Dim1, Dim2, Dim4, Dim0, Dim3, j, k, m, i, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim4, j, k, m>(TensorExpr(a, b));
}

//B(m, i, l, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim3, m, i, l> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim4, k, j, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim1, Dim4, i, j, k, l, k, j, m, Dim0, Dim3, Dim2, Dim1, Dim4, i, l, k, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, i, l, m>(TensorExpr(a, b));
}

//B(k, j, m, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim4, k, j, m> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim4, j, k, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim2, Dim4, i, j, k, l, j, k, m, Dim0, Dim3, Dim1, Dim2, Dim4, i, l, j, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, i, l, m>(TensorExpr(a, b));
}

//B(j, k, m, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim4, j, k, m> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim1, k, m, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim1, i, j, k, l, k, m, j, Dim0, Dim3, Dim2, Dim4, Dim1, i, l, k, m, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, i, l, m>(TensorExpr(a, b));
}

//B(k, m, j, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim1, k, m, j> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim2, j, m, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim2, i, j, k, l, j, m, k, Dim0, Dim3, Dim1, Dim4, Dim2, i, l, j, m, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, i, l, m>(TensorExpr(a, b));
}

//B(j, m, k, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim2, j, m, k> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim1, m, k, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, i, j, k, l, m, k, j, Dim0, Dim3, Dim4, Dim2, Dim1, i, l, m, k, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, i, l, m>(TensorExpr(a, b));
}

//B(m, k, j, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim1, m, k, j> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim2, m, j, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, i, j, k, l, m, j, k, Dim0, Dim3, Dim4, Dim1, Dim2, i, l, m, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim4, i, l, m>(TensorExpr(a, b));
}

//B(m, j, k, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim2, m, j, k> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim4, l, j, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim1, Dim4, i, j, k, l, l, j, m, Dim0, Dim2, Dim3, Dim1, Dim4, i, k, l, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, i, k, m>(TensorExpr(a, b));
}

//B(l, j, m, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim4, l, j, m> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim4, j, l, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim3, Dim4, i, j, k, l, j, l, m, Dim0, Dim2, Dim1, Dim3, Dim4, i, k, j, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, i, k, m>(TensorExpr(a, b));
}

//B(j, l, m, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim4, j, l, m> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim1, l, m, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim1, i, j, k, l, l, m, j, Dim0, Dim2, Dim3, Dim4, Dim1, i, k, l, m, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, i, k, m>(TensorExpr(a, b));
}

//B(l, m, j, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim1, l, m, j> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim3, j, m, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim3, i, j, k, l, j, m, l, Dim0, Dim2, Dim1, Dim4, Dim3, i, k, j, m, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, i, k, m>(TensorExpr(a, b));
}

//B(j, m, l, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim3, j, m, l> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim1, m, l, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, i, j, k, l, m, l, j, Dim0, Dim2, Dim4, Dim3, Dim1, i, k, m, l, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, i, k, m>(TensorExpr(a, b));
}

//B(m, l, j, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim1, m, l, j> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim3, m, j, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, i, j, k, l, m, j, l, Dim0, Dim2, Dim4, Dim1, Dim3, i, k, m, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim4, i, k, m>(TensorExpr(a, b));
}

//B(m, j, l, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim3, m, j, l> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim4, l, k, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim2, Dim4, i, j, k, l, l, k, m, Dim0, Dim1, Dim3, Dim2, Dim4, i, j, l, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, i, j, m>(TensorExpr(a, b));
}

//B(l, k, m, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim4, l, k, m> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim3, Dim4, i, j, k, l, k, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, i, j, m>(TensorExpr(a, b));
}

//B(k, l, m, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim2, l, m, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim2, i, j, k, l, l, m, k, Dim0, Dim1, Dim3, Dim4, Dim2, i, j, l, m, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, i, j, m>(TensorExpr(a, b));
}

//B(l, m, k, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim2, l, m, k> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim3, k, m, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim3, i, j, k, l, k, m, l, Dim0, Dim1, Dim2, Dim4, Dim3, i, j, k, m, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, i, j, m>(TensorExpr(a, b));
}

//B(k, m, l, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim3, k, m, l> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim2, m, l, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, i, j, k, l, m, l, k, Dim0, Dim1, Dim4, Dim3, Dim2, i, j, m, l, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, i, j, m>(TensorExpr(a, b));
}

//B(m, l, k, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim2, m, l, k> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim3, m, k, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, i, j, k, l, m, k, l, Dim0, Dim1, Dim4, Dim2, Dim3, i, j, m, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim4, i, j, m>(TensorExpr(a, b));
}

//B(m, k, l, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim3, m, k, l> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

} //namespace FTensor

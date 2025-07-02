#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, )*B(m, n, o, p, ) triple contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p, int DimA, int DimX, int DimY, int DimZ, int DimW, char a, char x, char y, char z, char w>
class Tensor4_times_Tensor4_triple
{
  Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
  Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim7, m, n, o, p> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX, ++xx)
      for(int yy = 0; yy < DimY, ++yy)
        for(int zz = 0; zz < DimZ, ++zz)
          {
            // Permutation is where the indices get checked.
            result += Permutation4<DimA, DimX, DimY, DimZ, a, x, y, z>().eval(
              iterA, N1, xx, yy, zz) * Permutation4<DimX, DimY, DimZ, DimW, x, y, z, w>().eval(iterB, xx, yy, zz, N2);
          }
    return result;
  }

  Tensor4_times_Tensor4_triple(
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &iter_a,
    const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim7, m, n, o, p> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, )*B(k, j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim4, k, j, i, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim1, Dim0, Dim4, i, j, k, l, k, j, i, m, Dim3, Dim2, Dim1, Dim0, Dim4, l, k, j, i, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim4, j, k, i, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim2, Dim0, Dim4, i, j, k, l, j, k, i, m, Dim3, Dim1, Dim2, Dim0, Dim4, l, j, k, i, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim4, k, i, j, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim0, Dim1, Dim4, i, j, k, l, k, i, j, m, Dim3, Dim2, Dim0, Dim1, Dim4, l, k, i, j, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim4, j, i, k, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim0, Dim2, Dim4, i, j, k, l, j, i, k, m, Dim3, Dim1, Dim0, Dim2, Dim4, l, j, i, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim4, i, k, j, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim2, Dim1, Dim4, i, j, k, l, i, k, j, m, Dim3, Dim0, Dim2, Dim1, Dim4, l, i, k, j, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim4, i, j, k, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim1, Dim2, Dim4, i, j, k, l, i, j, k, m, Dim3, Dim0, Dim1, Dim2, Dim4, l, i, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim0, k, j, m, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim1, Dim4, Dim0, i, j, k, l, k, j, m, i, Dim3, Dim2, Dim1, Dim4, Dim0, l, k, j, m, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim0, j, k, m, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim2, Dim4, Dim0, i, j, k, l, j, k, m, i, Dim3, Dim1, Dim2, Dim4, Dim0, l, j, k, m, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim1, k, i, m, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim0, Dim4, Dim1, i, j, k, l, k, i, m, j, Dim3, Dim2, Dim0, Dim4, Dim1, l, k, i, m, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim2, j, i, m, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim0, Dim4, Dim2, i, j, k, l, j, i, m, k, Dim3, Dim1, Dim0, Dim4, Dim2, l, j, i, m, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim1, i, k, m, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim2, Dim4, Dim1, i, j, k, l, i, k, m, j, Dim3, Dim0, Dim2, Dim4, Dim1, l, i, k, m, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim2, i, j, m, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim1, Dim4, Dim2, i, j, k, l, i, j, m, k, Dim3, Dim0, Dim1, Dim4, Dim2, l, i, j, m, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim0, k, m, j, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim1, Dim0, i, j, k, l, k, m, j, i, Dim3, Dim2, Dim4, Dim1, Dim0, l, k, m, j, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim0, j, m, k, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim2, Dim0, i, j, k, l, j, m, k, i, Dim3, Dim1, Dim4, Dim2, Dim0, l, j, m, k, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim1, k, m, i, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim0, Dim1, i, j, k, l, k, m, i, j, Dim3, Dim2, Dim4, Dim0, Dim1, l, k, m, i, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim2, j, m, i, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim0, Dim2, i, j, k, l, j, m, i, k, Dim3, Dim1, Dim4, Dim0, Dim2, l, j, m, i, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim1, i, m, k, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim2, Dim1, i, j, k, l, i, m, k, j, Dim3, Dim0, Dim4, Dim2, Dim1, l, i, m, k, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim2, i, m, j, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim1, Dim2, i, j, k, l, i, m, j, k, Dim3, Dim0, Dim4, Dim1, Dim2, l, i, m, j, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim0, m, k, j, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim0, i, j, k, l, m, k, j, i, Dim3, Dim4, Dim2, Dim1, Dim0, l, m, k, j, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim0, m, j, k, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim0, i, j, k, l, m, j, k, i, Dim3, Dim4, Dim1, Dim2, Dim0, l, m, j, k, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim1, m, k, i, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim1, i, j, k, l, m, k, i, j, Dim3, Dim4, Dim2, Dim0, Dim1, l, m, k, i, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim2, m, j, i, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim2, i, j, k, l, m, j, i, k, Dim3, Dim4, Dim1, Dim0, Dim2, l, m, j, i, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim1, m, i, k, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim1, i, j, k, l, m, i, k, j, Dim3, Dim4, Dim0, Dim2, Dim1, l, m, i, k, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim2, m, i, j, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim2, i, j, k, l, m, i, j, k, Dim3, Dim4, Dim0, Dim1, Dim2, l, m, i, j, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, l, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim4, l, j, i, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim1, Dim0, Dim4, i, j, k, l, l, j, i, m, Dim2, Dim3, Dim1, Dim0, Dim4, k, l, j, i, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim4, j, l, i, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim3, Dim0, Dim4, i, j, k, l, j, l, i, m, Dim2, Dim1, Dim3, Dim0, Dim4, k, j, l, i, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim4, l, i, j, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim0, Dim1, Dim4, i, j, k, l, l, i, j, m, Dim2, Dim3, Dim0, Dim1, Dim4, k, l, i, j, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim4, j, i, l, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim0, Dim3, Dim4, i, j, k, l, j, i, l, m, Dim2, Dim1, Dim0, Dim3, Dim4, k, j, i, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim4, i, l, j, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim3, Dim1, Dim4, i, j, k, l, i, l, j, m, Dim2, Dim0, Dim3, Dim1, Dim4, k, i, l, j, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim4, i, j, l, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim1, Dim3, Dim4, i, j, k, l, i, j, l, m, Dim2, Dim0, Dim1, Dim3, Dim4, k, i, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim0, l, j, m, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim1, Dim4, Dim0, i, j, k, l, l, j, m, i, Dim2, Dim3, Dim1, Dim4, Dim0, k, l, j, m, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim0, j, l, m, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim3, Dim4, Dim0, i, j, k, l, j, l, m, i, Dim2, Dim1, Dim3, Dim4, Dim0, k, j, l, m, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim1, l, i, m, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim0, Dim4, Dim1, i, j, k, l, l, i, m, j, Dim2, Dim3, Dim0, Dim4, Dim1, k, l, i, m, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim3, j, i, m, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim0, Dim4, Dim3, i, j, k, l, j, i, m, l, Dim2, Dim1, Dim0, Dim4, Dim3, k, j, i, m, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim1, i, l, m, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim3, Dim4, Dim1, i, j, k, l, i, l, m, j, Dim2, Dim0, Dim3, Dim4, Dim1, k, i, l, m, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim3, i, j, m, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim1, Dim4, Dim3, i, j, k, l, i, j, m, l, Dim2, Dim0, Dim1, Dim4, Dim3, k, i, j, m, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim0, l, m, j, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim1, Dim0, i, j, k, l, l, m, j, i, Dim2, Dim3, Dim4, Dim1, Dim0, k, l, m, j, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim0, j, m, l, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim3, Dim0, i, j, k, l, j, m, l, i, Dim2, Dim1, Dim4, Dim3, Dim0, k, j, m, l, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim1, l, m, i, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim0, Dim1, i, j, k, l, l, m, i, j, Dim2, Dim3, Dim4, Dim0, Dim1, k, l, m, i, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim3, j, m, i, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim0, Dim3, i, j, k, l, j, m, i, l, Dim2, Dim1, Dim4, Dim0, Dim3, k, j, m, i, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim1, i, m, l, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim3, Dim1, i, j, k, l, i, m, l, j, Dim2, Dim0, Dim4, Dim3, Dim1, k, i, m, l, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim3, i, m, j, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim1, Dim3, i, j, k, l, i, m, j, l, Dim2, Dim0, Dim4, Dim1, Dim3, k, i, m, j, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim0, m, l, j, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim0, i, j, k, l, m, l, j, i, Dim2, Dim4, Dim3, Dim1, Dim0, k, m, l, j, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim0, m, j, l, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim0, i, j, k, l, m, j, l, i, Dim2, Dim4, Dim1, Dim3, Dim0, k, m, j, l, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim1, m, l, i, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim1, i, j, k, l, m, l, i, j, Dim2, Dim4, Dim3, Dim0, Dim1, k, m, l, i, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim3, m, j, i, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim3, i, j, k, l, m, j, i, l, Dim2, Dim4, Dim1, Dim0, Dim3, k, m, j, i, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim1, m, i, l, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim1, i, j, k, l, m, i, l, j, Dim2, Dim4, Dim0, Dim3, Dim1, k, m, i, l, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim3, m, i, j, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim3, i, j, k, l, m, i, j, l, Dim2, Dim4, Dim0, Dim1, Dim3, k, m, i, j, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, k, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim4, l, k, i, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim2, Dim0, Dim4, i, j, k, l, l, k, i, m, Dim1, Dim3, Dim2, Dim0, Dim4, j, l, k, i, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim4, k, l, i, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim3, Dim0, Dim4, i, j, k, l, k, l, i, m, Dim1, Dim2, Dim3, Dim0, Dim4, j, k, l, i, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim4, l, i, k, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim0, Dim2, Dim4, i, j, k, l, l, i, k, m, Dim1, Dim3, Dim0, Dim2, Dim4, j, l, i, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim4, k, i, l, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim0, Dim3, Dim4, i, j, k, l, k, i, l, m, Dim1, Dim2, Dim0, Dim3, Dim4, j, k, i, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim4, i, l, k, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim3, Dim2, Dim4, i, j, k, l, i, l, k, m, Dim1, Dim0, Dim3, Dim2, Dim4, j, i, l, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim4, i, k, l, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim2, Dim3, Dim4, i, j, k, l, i, k, l, m, Dim1, Dim0, Dim2, Dim3, Dim4, j, i, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim0, l, k, m, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim2, Dim4, Dim0, i, j, k, l, l, k, m, i, Dim1, Dim3, Dim2, Dim4, Dim0, j, l, k, m, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim0, k, l, m, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim3, Dim4, Dim0, i, j, k, l, k, l, m, i, Dim1, Dim2, Dim3, Dim4, Dim0, j, k, l, m, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim2, l, i, m, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim0, Dim4, Dim2, i, j, k, l, l, i, m, k, Dim1, Dim3, Dim0, Dim4, Dim2, j, l, i, m, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim3, k, i, m, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim0, Dim4, Dim3, i, j, k, l, k, i, m, l, Dim1, Dim2, Dim0, Dim4, Dim3, j, k, i, m, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim2, i, l, m, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim3, Dim4, Dim2, i, j, k, l, i, l, m, k, Dim1, Dim0, Dim3, Dim4, Dim2, j, i, l, m, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim3, i, k, m, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim2, Dim4, Dim3, i, j, k, l, i, k, m, l, Dim1, Dim0, Dim2, Dim4, Dim3, j, i, k, m, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim0, l, m, k, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim2, Dim0, i, j, k, l, l, m, k, i, Dim1, Dim3, Dim4, Dim2, Dim0, j, l, m, k, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim0, k, m, l, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim3, Dim0, i, j, k, l, k, m, l, i, Dim1, Dim2, Dim4, Dim3, Dim0, j, k, m, l, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim2, l, m, i, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim0, Dim2, i, j, k, l, l, m, i, k, Dim1, Dim3, Dim4, Dim0, Dim2, j, l, m, i, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim3, k, m, i, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim0, Dim3, i, j, k, l, k, m, i, l, Dim1, Dim2, Dim4, Dim0, Dim3, j, k, m, i, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim2, i, m, l, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim3, Dim2, i, j, k, l, i, m, l, k, Dim1, Dim0, Dim4, Dim3, Dim2, j, i, m, l, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(i, m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim3, i, m, k, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim2, Dim3, i, j, k, l, i, m, k, l, Dim1, Dim0, Dim4, Dim2, Dim3, j, i, m, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim0, m, l, k, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim0, i, j, k, l, m, l, k, i, Dim1, Dim4, Dim3, Dim2, Dim0, j, m, l, k, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim0, m, k, l, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim0, i, j, k, l, m, k, l, i, Dim1, Dim4, Dim2, Dim3, Dim0, j, m, k, l, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim2, m, l, i, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim2, i, j, k, l, m, l, i, k, Dim1, Dim4, Dim3, Dim0, Dim2, j, m, l, i, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim3, m, k, i, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim3, i, j, k, l, m, k, i, l, Dim1, Dim4, Dim2, Dim0, Dim3, j, m, k, i, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim2, m, i, l, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim2, i, j, k, l, m, i, l, k, Dim1, Dim4, Dim0, Dim3, Dim2, j, m, i, l, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim3, m, i, k, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim3, i, j, k, l, m, i, k, l, Dim1, Dim4, Dim0, Dim2, Dim3, j, m, i, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, j, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim4, l, k, j, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim2, Dim1, Dim4, i, j, k, l, l, k, j, m, Dim0, Dim3, Dim2, Dim1, Dim4, i, l, k, j, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim4, k, l, j, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim3, Dim1, Dim4, i, j, k, l, k, l, j, m, Dim0, Dim2, Dim3, Dim1, Dim4, i, k, l, j, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim4, l, j, k, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim1, Dim2, Dim4, i, j, k, l, l, j, k, m, Dim0, Dim3, Dim1, Dim2, Dim4, i, l, j, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim4, k, j, l, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim1, Dim3, Dim4, i, j, k, l, k, j, l, m, Dim0, Dim2, Dim1, Dim3, Dim4, i, k, j, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim4, j, l, k, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim3, Dim2, Dim4, i, j, k, l, j, l, k, m, Dim0, Dim1, Dim3, Dim2, Dim4, i, j, l, k, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim4, j, k, l, m> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim2, Dim3, Dim4, i, j, k, l, j, k, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim1, l, k, m, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim2, Dim4, Dim1, i, j, k, l, l, k, m, j, Dim0, Dim3, Dim2, Dim4, Dim1, i, l, k, m, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim1, k, l, m, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim3, Dim4, Dim1, i, j, k, l, k, l, m, j, Dim0, Dim2, Dim3, Dim4, Dim1, i, k, l, m, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim2, l, j, m, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim1, Dim4, Dim2, i, j, k, l, l, j, m, k, Dim0, Dim3, Dim1, Dim4, Dim2, i, l, j, m, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim3, k, j, m, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim1, Dim4, Dim3, i, j, k, l, k, j, m, l, Dim0, Dim2, Dim1, Dim4, Dim3, i, k, j, m, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim2, j, l, m, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim3, Dim4, Dim2, i, j, k, l, j, l, m, k, Dim0, Dim1, Dim3, Dim4, Dim2, i, j, l, m, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim3, j, k, m, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim2, Dim4, Dim3, i, j, k, l, j, k, m, l, Dim0, Dim1, Dim2, Dim4, Dim3, i, j, k, m, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim1, l, m, k, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim2, Dim1, i, j, k, l, l, m, k, j, Dim0, Dim3, Dim4, Dim2, Dim1, i, l, m, k, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim1, k, m, l, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim3, Dim1, i, j, k, l, k, m, l, j, Dim0, Dim2, Dim4, Dim3, Dim1, i, k, m, l, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim2, l, m, j, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim1, Dim2, i, j, k, l, l, m, j, k, Dim0, Dim3, Dim4, Dim1, Dim2, i, l, m, j, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim3, k, m, j, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim1, Dim3, i, j, k, l, k, m, j, l, Dim0, Dim2, Dim4, Dim1, Dim3, i, k, m, j, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim2, j, m, l, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim3, Dim2, i, j, k, l, j, m, l, k, Dim0, Dim1, Dim4, Dim3, Dim2, i, j, m, l, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim3, j, m, k, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim2, Dim3, i, j, k, l, j, m, k, l, Dim0, Dim1, Dim4, Dim2, Dim3, i, j, m, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim1, m, l, k, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim1, i, j, k, l, m, l, k, j, Dim0, Dim4, Dim3, Dim2, Dim1, i, m, l, k, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim1, m, k, l, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim1, i, j, k, l, m, k, l, j, Dim0, Dim4, Dim2, Dim3, Dim1, i, m, k, l, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim2, m, l, j, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim2, i, j, k, l, m, l, j, k, Dim0, Dim4, Dim3, Dim1, Dim2, i, m, l, j, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim3, m, k, j, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim3, i, j, k, l, m, k, j, l, Dim0, Dim4, Dim2, Dim1, Dim3, i, m, k, j, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim2, m, j, l, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim2, i, j, k, l, m, j, l, k, Dim0, Dim4, Dim1, Dim3, Dim2, i, m, j, l, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim3, m, j, k, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim3, i, j, k, l, m, j, k, l, Dim0, Dim4, Dim1, Dim2, Dim3, i, m, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, i, m>(TensorExpr(a, b));
}

} //namespace FTensor

#pragma once

namespace FTensor
{
//Universal A(i, j, k, )*B(l, m, n, ) double contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n, int DimA, int DimX, int DimY, int DimZ, char a, char x, char y, char z>
class Tensor3_times_Tensor3_double
{
  Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
  Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX, ++xx)
      for(int yy = 0; yy < DimY, ++yy)
        {
          // Permutation is where the indices get checked.
          result += Permutation3<DimA, DimX, DimY, a, x, y>().eval(
            iterA, N1, xx, yy) * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx, yy, N2);
        }
    return result;
  }

  Tensor3_times_Tensor3_double(
    const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &iter_a,
    const Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, )*B(j, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim1, Dim0, Dim3, j, i, l> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim1, Dim0, Dim3, i, j, k, j, i, l, Dim2, Dim1, Dim0, Dim3, k, j, i, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(i, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim0, Dim1, Dim3, i, j, l> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim0, Dim1, Dim3, i, j, k, i, j, l, Dim2, Dim0, Dim1, Dim3, k, i, j, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(j, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim0, j, l, i> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim1, Dim3, Dim0, i, j, k, j, l, i, Dim2, Dim1, Dim3, Dim0, k, j, l, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(i, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim1, i, l, j> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim0, Dim3, Dim1, i, j, k, i, l, j, Dim2, Dim0, Dim3, Dim1, k, i, l, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(l, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim0, l, j, i> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim0, i, j, k, l, j, i, Dim2, Dim3, Dim1, Dim0, k, l, j, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(l, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim1, l, i, j> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim1, i, j, k, l, i, j, Dim2, Dim3, Dim0, Dim1, k, l, i, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, k, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(k, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim2, Dim0, Dim3, k, i, l> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim2, Dim0, Dim3, i, j, k, k, i, l, Dim1, Dim2, Dim0, Dim3, j, k, i, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(i, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim0, Dim2, Dim3, i, k, l> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim0, Dim2, Dim3, i, j, k, i, k, l, Dim1, Dim0, Dim2, Dim3, j, i, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(k, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim0, k, l, i> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim2, Dim3, Dim0, i, j, k, k, l, i, Dim1, Dim2, Dim3, Dim0, j, k, l, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(i, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim0, Dim3, Dim2, i, l, k> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim0, Dim3, Dim2, i, j, k, i, l, k, Dim1, Dim0, Dim3, Dim2, j, i, l, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(l, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim0, l, k, i> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim0, i, j, k, l, k, i, Dim1, Dim3, Dim2, Dim0, j, l, k, i>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(l, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim3, Dim0, Dim2, l, i, k> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim2, i, j, k, l, i, k, Dim1, Dim3, Dim0, Dim2, j, l, i, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, j, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(k, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim2, Dim1, Dim3, k, j, l> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim2, Dim1, Dim3, i, j, k, k, j, l, Dim0, Dim2, Dim1, Dim3, i, k, j, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim1, Dim2, Dim3, j, k, l> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim1, Dim2, Dim3, i, j, k, j, k, l, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(k, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim2, Dim3, Dim1, k, l, j> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim2, Dim3, Dim1, i, j, k, k, l, j, Dim0, Dim2, Dim3, Dim1, i, k, l, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(j, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim1, Dim3, Dim2, j, l, k> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim1, Dim3, Dim2, i, j, k, j, l, k, Dim0, Dim1, Dim3, Dim2, i, j, l, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(l, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim3, Dim2, Dim1, l, k, j> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim1, i, j, k, l, k, j, Dim0, Dim3, Dim2, Dim1, i, l, k, j>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

//A(i, j, k, )*B(l, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
               const Tensor3_Expr<B, U, Dim3, Dim1, Dim2, l, j, k> &b)
{
  using TensorExpr
    = Tensor3_times_Tensor3_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim2, i, j, k, l, j, k, Dim0, Dim3, Dim1, Dim2, i, l, j, k>;
  return Tensor2_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, i, l>(TensorExpr(a, b));
}

} //namespace FTensor

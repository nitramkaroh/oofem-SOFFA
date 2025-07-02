#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, )*B(n, o, p, q, ) triple contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, int Dim8, char i, char j, char k, char l, char m, char n, char o, char p, char q, int DimA, int DimB, int DimX, int DimY, int DimZ, int DimW, char a, char b, char x, char y, char z, char w>
class Tensor5_times_Tensor4_triple
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim8, n, o, p, q> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX, ++xx)
      for(int yy = 0; yy < DimY, ++yy)
        for(int zz = 0; zz < DimZ, ++zz)
          {
            // Permutation is where the indices get checked.
            result += Permutation5<DimA, DimB, DimX, DimY, DimZ, a, b, x, y, z>().eval(
              iterA, N1, N2, xx, yy, zz) * Permutation4<DimX, DimY, DimZ, DimW, x, y, z, w>().eval(iterB, xx, yy, zz, N3);
          }
    return result;
  }

  Tensor5_times_Tensor4_triple(
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &iter_a,
    const Tensor4_Expr<B, U, Dim5, Dim6, Dim7, Dim8, n, o, p, q> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, )*B(k, j, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim5, k, j, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim0, Dim5, i, j, k, l, m, k, j, i, n, Dim3, Dim4, Dim2, Dim1, Dim0, Dim5, l, m, k, j, i, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(k, j, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim5, k, j, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim5, j, k, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim0, Dim5, i, j, k, l, m, j, k, i, n, Dim3, Dim4, Dim1, Dim2, Dim0, Dim5, l, m, j, k, i, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(j, k, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim5, j, k, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim5, k, i, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim1, Dim5, i, j, k, l, m, k, i, j, n, Dim3, Dim4, Dim2, Dim0, Dim1, Dim5, l, m, k, i, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(k, i, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim5, k, i, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, i, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim5, j, i, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim2, Dim5, i, j, k, l, m, j, i, k, n, Dim3, Dim4, Dim1, Dim0, Dim2, Dim5, l, m, j, i, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(j, i, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim5, j, i, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim5, i, k, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim1, Dim5, i, j, k, l, m, i, k, j, n, Dim3, Dim4, Dim0, Dim2, Dim1, Dim5, l, m, i, k, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(i, k, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim5, i, k, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim5, i, j, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, i, j, k, l, m, i, j, k, n, Dim3, Dim4, Dim0, Dim1, Dim2, Dim5, l, m, i, j, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(i, j, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim5, i, j, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim0, k, j, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, Dim0, i, j, k, l, m, k, j, n, i, Dim3, Dim4, Dim2, Dim1, Dim5, Dim0, l, m, k, j, n, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(k, j, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim0, k, j, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim0, j, k, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, Dim0, i, j, k, l, m, j, k, n, i, Dim3, Dim4, Dim1, Dim2, Dim5, Dim0, l, m, j, k, n, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(j, k, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim0, j, k, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim1, k, i, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, Dim1, i, j, k, l, m, k, i, n, j, Dim3, Dim4, Dim2, Dim0, Dim5, Dim1, l, m, k, i, n, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(k, i, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim1, k, i, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, i, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim2, j, i, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, Dim2, i, j, k, l, m, j, i, n, k, Dim3, Dim4, Dim1, Dim0, Dim5, Dim2, l, m, j, i, n, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(j, i, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim2, j, i, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim1, i, k, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, Dim1, i, j, k, l, m, i, k, n, j, Dim3, Dim4, Dim0, Dim2, Dim5, Dim1, l, m, i, k, n, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(i, k, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim1, i, k, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim2, i, j, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim2, i, j, k, l, m, i, j, n, k, Dim3, Dim4, Dim0, Dim1, Dim5, Dim2, l, m, i, j, n, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(i, j, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim2, i, j, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim0, k, n, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, Dim0, i, j, k, l, m, k, n, j, i, Dim3, Dim4, Dim2, Dim5, Dim1, Dim0, l, m, k, n, j, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(k, n, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim0, k, n, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim0, j, n, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, Dim0, i, j, k, l, m, j, n, k, i, Dim3, Dim4, Dim1, Dim5, Dim2, Dim0, l, m, j, n, k, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(j, n, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim0, j, n, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim1, k, n, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, Dim1, i, j, k, l, m, k, n, i, j, Dim3, Dim4, Dim2, Dim5, Dim0, Dim1, l, m, k, n, i, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(k, n, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim1, k, n, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim2, j, n, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, Dim2, i, j, k, l, m, j, n, i, k, Dim3, Dim4, Dim1, Dim5, Dim0, Dim2, l, m, j, n, i, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(j, n, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim2, j, n, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim1, i, n, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, Dim1, i, j, k, l, m, i, n, k, j, Dim3, Dim4, Dim0, Dim5, Dim2, Dim1, l, m, i, n, k, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(i, n, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim1, i, n, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim2, i, n, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, Dim2, i, j, k, l, m, i, n, j, k, Dim3, Dim4, Dim0, Dim5, Dim1, Dim2, l, m, i, n, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(i, n, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim2, i, n, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim0, n, k, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim0, i, j, k, l, m, n, k, j, i, Dim3, Dim4, Dim5, Dim2, Dim1, Dim0, l, m, n, k, j, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(n, k, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim0, n, k, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim0, n, j, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim0, i, j, k, l, m, n, j, k, i, Dim3, Dim4, Dim5, Dim1, Dim2, Dim0, l, m, n, j, k, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(n, j, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim0, n, j, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim1, n, k, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim1, i, j, k, l, m, n, k, i, j, Dim3, Dim4, Dim5, Dim2, Dim0, Dim1, l, m, n, k, i, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(n, k, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim1, n, k, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim2, n, j, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim2, i, j, k, l, m, n, j, i, k, Dim3, Dim4, Dim5, Dim1, Dim0, Dim2, l, m, n, j, i, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(n, j, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim2, n, j, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim1, n, i, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim1, i, j, k, l, m, n, i, k, j, Dim3, Dim4, Dim5, Dim0, Dim2, Dim1, l, m, n, i, k, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(n, i, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim1, n, i, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim2, n, i, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, i, j, k, l, m, n, i, j, k, Dim3, Dim4, Dim5, Dim0, Dim1, Dim2, l, m, n, i, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim3, Dim4, Dim5, l, m, n>(TensorExpr(a, b));
}

//B(n, i, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim2, n, i, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim5, l, j, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim0, Dim5, i, j, k, l, m, l, j, i, n, Dim2, Dim4, Dim3, Dim1, Dim0, Dim5, k, m, l, j, i, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(l, j, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim5, l, j, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim5, j, l, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim0, Dim5, i, j, k, l, m, j, l, i, n, Dim2, Dim4, Dim1, Dim3, Dim0, Dim5, k, m, j, l, i, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(j, l, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim5, j, l, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim5, l, i, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim1, Dim5, i, j, k, l, m, l, i, j, n, Dim2, Dim4, Dim3, Dim0, Dim1, Dim5, k, m, l, i, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(l, i, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim5, l, i, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, i, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim5, j, i, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim3, Dim5, i, j, k, l, m, j, i, l, n, Dim2, Dim4, Dim1, Dim0, Dim3, Dim5, k, m, j, i, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(j, i, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim5, j, i, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim5, i, l, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim1, Dim5, i, j, k, l, m, i, l, j, n, Dim2, Dim4, Dim0, Dim3, Dim1, Dim5, k, m, i, l, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(i, l, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim5, i, l, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim5, i, j, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim3, Dim5, i, j, k, l, m, i, j, l, n, Dim2, Dim4, Dim0, Dim1, Dim3, Dim5, k, m, i, j, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(i, j, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim5, i, j, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim0, l, j, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, Dim0, i, j, k, l, m, l, j, n, i, Dim2, Dim4, Dim3, Dim1, Dim5, Dim0, k, m, l, j, n, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(l, j, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim0, l, j, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim0, j, l, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, Dim0, i, j, k, l, m, j, l, n, i, Dim2, Dim4, Dim1, Dim3, Dim5, Dim0, k, m, j, l, n, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(j, l, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim0, j, l, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim1, l, i, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, Dim1, i, j, k, l, m, l, i, n, j, Dim2, Dim4, Dim3, Dim0, Dim5, Dim1, k, m, l, i, n, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(l, i, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim1, l, i, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, i, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim3, j, i, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, Dim3, i, j, k, l, m, j, i, n, l, Dim2, Dim4, Dim1, Dim0, Dim5, Dim3, k, m, j, i, n, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(j, i, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim3, j, i, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim1, i, l, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, Dim1, i, j, k, l, m, i, l, n, j, Dim2, Dim4, Dim0, Dim3, Dim5, Dim1, k, m, i, l, n, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(i, l, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim1, i, l, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim3, i, j, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim3, i, j, k, l, m, i, j, n, l, Dim2, Dim4, Dim0, Dim1, Dim5, Dim3, k, m, i, j, n, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(i, j, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim3, i, j, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim0, l, n, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, Dim0, i, j, k, l, m, l, n, j, i, Dim2, Dim4, Dim3, Dim5, Dim1, Dim0, k, m, l, n, j, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(l, n, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim0, l, n, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim0, j, n, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, Dim0, i, j, k, l, m, j, n, l, i, Dim2, Dim4, Dim1, Dim5, Dim3, Dim0, k, m, j, n, l, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(j, n, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim0, j, n, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim1, l, n, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, Dim1, i, j, k, l, m, l, n, i, j, Dim2, Dim4, Dim3, Dim5, Dim0, Dim1, k, m, l, n, i, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(l, n, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim1, l, n, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim3, j, n, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, Dim3, i, j, k, l, m, j, n, i, l, Dim2, Dim4, Dim1, Dim5, Dim0, Dim3, k, m, j, n, i, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(j, n, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim3, j, n, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim1, i, n, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, Dim1, i, j, k, l, m, i, n, l, j, Dim2, Dim4, Dim0, Dim5, Dim3, Dim1, k, m, i, n, l, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(i, n, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim1, i, n, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim3, i, n, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, Dim3, i, j, k, l, m, i, n, j, l, Dim2, Dim4, Dim0, Dim5, Dim1, Dim3, k, m, i, n, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(i, n, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim3, i, n, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim0, n, l, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim0, i, j, k, l, m, n, l, j, i, Dim2, Dim4, Dim5, Dim3, Dim1, Dim0, k, m, n, l, j, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(n, l, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim0, n, l, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim0, n, j, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim0, i, j, k, l, m, n, j, l, i, Dim2, Dim4, Dim5, Dim1, Dim3, Dim0, k, m, n, j, l, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(n, j, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim0, n, j, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim1, n, l, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim1, i, j, k, l, m, n, l, i, j, Dim2, Dim4, Dim5, Dim3, Dim0, Dim1, k, m, n, l, i, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(n, l, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim1, n, l, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim3, n, j, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim3, i, j, k, l, m, n, j, i, l, Dim2, Dim4, Dim5, Dim1, Dim0, Dim3, k, m, n, j, i, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(n, j, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim3, n, j, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim1, n, i, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim1, i, j, k, l, m, n, i, l, j, Dim2, Dim4, Dim5, Dim0, Dim3, Dim1, k, m, n, i, l, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(n, i, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim1, n, i, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim3, n, i, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim3, i, j, k, l, m, n, i, j, l, Dim2, Dim4, Dim5, Dim0, Dim1, Dim3, k, m, n, i, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim4, Dim5, k, m, n>(TensorExpr(a, b));
}

//B(n, i, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim3, n, i, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim5, m, j, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim0, Dim5, i, j, k, l, m, m, j, i, n, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, k, l, m, j, i, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(m, j, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim0, Dim5, m, j, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim5, j, m, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim0, Dim5, i, j, k, l, m, j, m, i, n, Dim2, Dim3, Dim1, Dim4, Dim0, Dim5, k, l, j, m, i, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(j, m, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim0, Dim5, j, m, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim5, m, i, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim1, Dim5, i, j, k, l, m, m, i, j, n, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, k, l, m, i, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(m, i, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim1, Dim5, m, i, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim5, j, i, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim4, Dim5, i, j, k, l, m, j, i, m, n, Dim2, Dim3, Dim1, Dim0, Dim4, Dim5, k, l, j, i, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(j, i, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim4, Dim5, j, i, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim5, i, m, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim1, Dim5, i, j, k, l, m, i, m, j, n, Dim2, Dim3, Dim0, Dim4, Dim1, Dim5, k, l, i, m, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(i, m, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim1, Dim5, i, m, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim5, i, j, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim4, Dim5, i, j, k, l, m, i, j, m, n, Dim2, Dim3, Dim0, Dim1, Dim4, Dim5, k, l, i, j, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(i, j, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim4, Dim5, i, j, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim0, m, j, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim5, Dim0, i, j, k, l, m, m, j, n, i, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, k, l, m, j, n, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(m, j, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim0, m, j, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim0, j, m, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim5, Dim0, i, j, k, l, m, j, m, n, i, Dim2, Dim3, Dim1, Dim4, Dim5, Dim0, k, l, j, m, n, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(j, m, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim0, j, m, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim1, m, i, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim5, Dim1, i, j, k, l, m, m, i, n, j, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, k, l, m, i, n, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(m, i, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim1, m, i, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, i, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim4, j, i, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim5, Dim4, i, j, k, l, m, j, i, n, m, Dim2, Dim3, Dim1, Dim0, Dim5, Dim4, k, l, j, i, n, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(j, i, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim0, Dim5, Dim4, j, i, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim1, i, m, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim5, Dim1, i, j, k, l, m, i, m, n, j, Dim2, Dim3, Dim0, Dim4, Dim5, Dim1, k, l, i, m, n, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(i, m, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim1, i, m, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, j, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim4, i, j, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim5, Dim4, i, j, k, l, m, i, j, n, m, Dim2, Dim3, Dim0, Dim1, Dim5, Dim4, k, l, i, j, n, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(i, j, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim1, Dim5, Dim4, i, j, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim0, m, n, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim1, Dim0, i, j, k, l, m, m, n, j, i, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, k, l, m, n, j, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(m, n, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim0, m, n, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim0, j, n, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim4, Dim0, i, j, k, l, m, j, n, m, i, Dim2, Dim3, Dim1, Dim5, Dim4, Dim0, k, l, j, n, m, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(j, n, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim0, j, n, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim1, m, n, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim0, Dim1, i, j, k, l, m, m, n, i, j, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, k, l, m, n, i, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(m, n, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim1, m, n, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim4, j, n, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim0, Dim4, i, j, k, l, m, j, n, i, m, Dim2, Dim3, Dim1, Dim5, Dim0, Dim4, k, l, j, n, i, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(j, n, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim0, Dim4, j, n, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim1, i, n, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim4, Dim1, i, j, k, l, m, i, n, m, j, Dim2, Dim3, Dim0, Dim5, Dim4, Dim1, k, l, i, n, m, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(i, n, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim1, i, n, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim4, i, n, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim1, Dim4, i, j, k, l, m, i, n, j, m, Dim2, Dim3, Dim0, Dim5, Dim1, Dim4, k, l, i, n, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(i, n, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim1, Dim4, i, n, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, j, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim0, n, m, j, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim0, i, j, k, l, m, n, m, j, i, Dim2, Dim3, Dim5, Dim4, Dim1, Dim0, k, l, n, m, j, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(n, m, j, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim0, n, m, j, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim0, n, j, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim0, i, j, k, l, m, n, j, m, i, Dim2, Dim3, Dim5, Dim1, Dim4, Dim0, k, l, n, j, m, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(n, j, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim0, n, j, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, i, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim1, n, m, i, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim1, i, j, k, l, m, n, m, i, j, Dim2, Dim3, Dim5, Dim4, Dim0, Dim1, k, l, n, m, i, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(n, m, i, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim1, n, m, i, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim4, n, j, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, Dim4, i, j, k, l, m, n, j, i, m, Dim2, Dim3, Dim5, Dim1, Dim0, Dim4, k, l, n, j, i, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(n, j, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim0, Dim4, n, j, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim1, n, i, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim1, i, j, k, l, m, n, i, m, j, Dim2, Dim3, Dim5, Dim0, Dim4, Dim1, k, l, n, i, m, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(n, i, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim1, n, i, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim4, n, i, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, Dim4, i, j, k, l, m, n, i, j, m, Dim2, Dim3, Dim5, Dim0, Dim1, Dim4, k, l, n, i, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim2, Dim3, Dim5, k, l, n>(TensorExpr(a, b));
}

//B(n, i, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim1, Dim4, n, i, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim5, l, k, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim0, Dim5, i, j, k, l, m, l, k, i, n, Dim1, Dim4, Dim3, Dim2, Dim0, Dim5, j, m, l, k, i, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(l, k, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim5, l, k, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim5, k, l, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim0, Dim5, i, j, k, l, m, k, l, i, n, Dim1, Dim4, Dim2, Dim3, Dim0, Dim5, j, m, k, l, i, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(k, l, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim5, k, l, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim5, l, i, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim2, Dim5, i, j, k, l, m, l, i, k, n, Dim1, Dim4, Dim3, Dim0, Dim2, Dim5, j, m, l, i, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(l, i, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim5, l, i, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim5, k, i, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim3, Dim5, i, j, k, l, m, k, i, l, n, Dim1, Dim4, Dim2, Dim0, Dim3, Dim5, j, m, k, i, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(k, i, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim5, k, i, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim5, i, l, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim2, Dim5, i, j, k, l, m, i, l, k, n, Dim1, Dim4, Dim0, Dim3, Dim2, Dim5, j, m, i, l, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(i, l, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim5, i, l, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim5, i, k, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim3, Dim5, i, j, k, l, m, i, k, l, n, Dim1, Dim4, Dim0, Dim2, Dim3, Dim5, j, m, i, k, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(i, k, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim5, i, k, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim0, l, k, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, Dim0, i, j, k, l, m, l, k, n, i, Dim1, Dim4, Dim3, Dim2, Dim5, Dim0, j, m, l, k, n, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(l, k, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim0, l, k, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim0, k, l, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, Dim0, i, j, k, l, m, k, l, n, i, Dim1, Dim4, Dim2, Dim3, Dim5, Dim0, j, m, k, l, n, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(k, l, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim0, k, l, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim2, l, i, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, Dim2, i, j, k, l, m, l, i, n, k, Dim1, Dim4, Dim3, Dim0, Dim5, Dim2, j, m, l, i, n, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(l, i, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim2, l, i, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim3, k, i, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, Dim3, i, j, k, l, m, k, i, n, l, Dim1, Dim4, Dim2, Dim0, Dim5, Dim3, j, m, k, i, n, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(k, i, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim3, k, i, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim2, i, l, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, Dim2, i, j, k, l, m, i, l, n, k, Dim1, Dim4, Dim0, Dim3, Dim5, Dim2, j, m, i, l, n, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(i, l, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim2, i, l, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim3, i, k, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, Dim3, i, j, k, l, m, i, k, n, l, Dim1, Dim4, Dim0, Dim2, Dim5, Dim3, j, m, i, k, n, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(i, k, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim3, i, k, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim0, l, n, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, Dim0, i, j, k, l, m, l, n, k, i, Dim1, Dim4, Dim3, Dim5, Dim2, Dim0, j, m, l, n, k, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(l, n, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim0, l, n, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim0, k, n, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, Dim0, i, j, k, l, m, k, n, l, i, Dim1, Dim4, Dim2, Dim5, Dim3, Dim0, j, m, k, n, l, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(k, n, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim0, k, n, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim2, l, n, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, Dim2, i, j, k, l, m, l, n, i, k, Dim1, Dim4, Dim3, Dim5, Dim0, Dim2, j, m, l, n, i, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(l, n, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim2, l, n, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim3, k, n, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, Dim3, i, j, k, l, m, k, n, i, l, Dim1, Dim4, Dim2, Dim5, Dim0, Dim3, j, m, k, n, i, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(k, n, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim3, k, n, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim2, i, n, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, Dim2, i, j, k, l, m, i, n, l, k, Dim1, Dim4, Dim0, Dim5, Dim3, Dim2, j, m, i, n, l, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(i, n, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim2, i, n, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim3, i, n, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, Dim3, i, j, k, l, m, i, n, k, l, Dim1, Dim4, Dim0, Dim5, Dim2, Dim3, j, m, i, n, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(i, n, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim3, i, n, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim0, n, l, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim0, i, j, k, l, m, n, l, k, i, Dim1, Dim4, Dim5, Dim3, Dim2, Dim0, j, m, n, l, k, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(n, l, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim0, n, l, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim0, n, k, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim0, i, j, k, l, m, n, k, l, i, Dim1, Dim4, Dim5, Dim2, Dim3, Dim0, j, m, n, k, l, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(n, k, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim0, n, k, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim2, n, l, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim2, i, j, k, l, m, n, l, i, k, Dim1, Dim4, Dim5, Dim3, Dim0, Dim2, j, m, n, l, i, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(n, l, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim2, n, l, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim3, n, k, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim3, i, j, k, l, m, n, k, i, l, Dim1, Dim4, Dim5, Dim2, Dim0, Dim3, j, m, n, k, i, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(n, k, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim3, n, k, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim2, n, i, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim2, i, j, k, l, m, n, i, l, k, Dim1, Dim4, Dim5, Dim0, Dim3, Dim2, j, m, n, i, l, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(n, i, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim2, n, i, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim3, n, i, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim3, i, j, k, l, m, n, i, k, l, Dim1, Dim4, Dim5, Dim0, Dim2, Dim3, j, m, n, i, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim4, Dim5, j, m, n>(TensorExpr(a, b));
}

//B(n, i, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim3, n, i, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim5, m, k, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim0, Dim5, i, j, k, l, m, m, k, i, n, Dim1, Dim3, Dim4, Dim2, Dim0, Dim5, j, l, m, k, i, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(m, k, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim0, Dim5, m, k, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim5, k, m, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim0, Dim5, i, j, k, l, m, k, m, i, n, Dim1, Dim3, Dim2, Dim4, Dim0, Dim5, j, l, k, m, i, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(k, m, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim0, Dim5, k, m, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim5, m, i, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim2, Dim5, i, j, k, l, m, m, i, k, n, Dim1, Dim3, Dim4, Dim0, Dim2, Dim5, j, l, m, i, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(m, i, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim2, Dim5, m, i, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim5, k, i, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim4, Dim5, i, j, k, l, m, k, i, m, n, Dim1, Dim3, Dim2, Dim0, Dim4, Dim5, j, l, k, i, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(k, i, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim4, Dim5, k, i, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim5, i, m, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim2, Dim5, i, j, k, l, m, i, m, k, n, Dim1, Dim3, Dim0, Dim4, Dim2, Dim5, j, l, i, m, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(i, m, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim2, Dim5, i, m, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim5, i, k, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim4, Dim5, i, j, k, l, m, i, k, m, n, Dim1, Dim3, Dim0, Dim2, Dim4, Dim5, j, l, i, k, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(i, k, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim4, Dim5, i, k, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim0, m, k, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim5, Dim0, i, j, k, l, m, m, k, n, i, Dim1, Dim3, Dim4, Dim2, Dim5, Dim0, j, l, m, k, n, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(m, k, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim0, m, k, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim0, k, m, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim5, Dim0, i, j, k, l, m, k, m, n, i, Dim1, Dim3, Dim2, Dim4, Dim5, Dim0, j, l, k, m, n, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(k, m, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim0, k, m, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim2, m, i, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim5, Dim2, i, j, k, l, m, m, i, n, k, Dim1, Dim3, Dim4, Dim0, Dim5, Dim2, j, l, m, i, n, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(m, i, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim2, m, i, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, i, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim4, k, i, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim5, Dim4, i, j, k, l, m, k, i, n, m, Dim1, Dim3, Dim2, Dim0, Dim5, Dim4, j, l, k, i, n, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(k, i, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim0, Dim5, Dim4, k, i, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim2, i, m, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim5, Dim2, i, j, k, l, m, i, m, n, k, Dim1, Dim3, Dim0, Dim4, Dim5, Dim2, j, l, i, m, n, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(i, m, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim2, i, m, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, k, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim4, i, k, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim5, Dim4, i, j, k, l, m, i, k, n, m, Dim1, Dim3, Dim0, Dim2, Dim5, Dim4, j, l, i, k, n, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(i, k, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim2, Dim5, Dim4, i, k, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim0, m, n, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim2, Dim0, i, j, k, l, m, m, n, k, i, Dim1, Dim3, Dim4, Dim5, Dim2, Dim0, j, l, m, n, k, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(m, n, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim0, m, n, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim0, k, n, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim4, Dim0, i, j, k, l, m, k, n, m, i, Dim1, Dim3, Dim2, Dim5, Dim4, Dim0, j, l, k, n, m, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(k, n, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim0, k, n, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim2, m, n, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim0, Dim2, i, j, k, l, m, m, n, i, k, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, j, l, m, n, i, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(m, n, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim2, m, n, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim4, k, n, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim0, Dim4, i, j, k, l, m, k, n, i, m, Dim1, Dim3, Dim2, Dim5, Dim0, Dim4, j, l, k, n, i, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(k, n, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim0, Dim4, k, n, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim2, i, n, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim4, Dim2, i, j, k, l, m, i, n, m, k, Dim1, Dim3, Dim0, Dim5, Dim4, Dim2, j, l, i, n, m, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(i, n, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim2, i, n, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim4, i, n, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim2, Dim4, i, j, k, l, m, i, n, k, m, Dim1, Dim3, Dim0, Dim5, Dim2, Dim4, j, l, i, n, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(i, n, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim2, Dim4, i, n, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, k, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim0, n, m, k, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim0, i, j, k, l, m, n, m, k, i, Dim1, Dim3, Dim5, Dim4, Dim2, Dim0, j, l, n, m, k, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(n, m, k, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim0, n, m, k, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim0, n, k, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim0, i, j, k, l, m, n, k, m, i, Dim1, Dim3, Dim5, Dim2, Dim4, Dim0, j, l, n, k, m, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(n, k, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim0, n, k, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, i, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim2, n, m, i, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim2, i, j, k, l, m, n, m, i, k, Dim1, Dim3, Dim5, Dim4, Dim0, Dim2, j, l, n, m, i, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(n, m, i, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim2, n, m, i, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim4, n, k, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, Dim4, i, j, k, l, m, n, k, i, m, Dim1, Dim3, Dim5, Dim2, Dim0, Dim4, j, l, n, k, i, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(n, k, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim0, Dim4, n, k, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim2, n, i, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim2, i, j, k, l, m, n, i, m, k, Dim1, Dim3, Dim5, Dim0, Dim4, Dim2, j, l, n, i, m, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(n, i, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim2, n, i, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim4, n, i, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, Dim4, i, j, k, l, m, n, i, k, m, Dim1, Dim3, Dim5, Dim0, Dim2, Dim4, j, l, n, i, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim3, Dim5, j, l, n>(TensorExpr(a, b));
}

//B(n, i, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim2, Dim4, n, i, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim5, m, l, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim0, Dim5, i, j, k, l, m, m, l, i, n, Dim1, Dim2, Dim4, Dim3, Dim0, Dim5, j, k, m, l, i, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(m, l, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim0, Dim5, m, l, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim5, l, m, i, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim0, Dim5, i, j, k, l, m, l, m, i, n, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, j, k, l, m, i, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(l, m, i, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim0, Dim5, l, m, i, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim5, m, i, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim3, Dim5, i, j, k, l, m, m, i, l, n, Dim1, Dim2, Dim4, Dim0, Dim3, Dim5, j, k, m, i, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(m, i, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim3, Dim5, m, i, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim5, l, i, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim4, Dim5, i, j, k, l, m, l, i, m, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(l, i, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim4, Dim5, l, i, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim5, i, m, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim3, Dim5, i, j, k, l, m, i, m, l, n, Dim1, Dim2, Dim0, Dim4, Dim3, Dim5, j, k, i, m, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(i, m, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim3, Dim5, i, m, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim5, i, l, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim4, Dim5, i, j, k, l, m, i, l, m, n, Dim1, Dim2, Dim0, Dim3, Dim4, Dim5, j, k, i, l, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(i, l, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim4, Dim5, i, l, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim0, m, l, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim5, Dim0, i, j, k, l, m, m, l, n, i, Dim1, Dim2, Dim4, Dim3, Dim5, Dim0, j, k, m, l, n, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(m, l, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim0, m, l, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim0, l, m, n, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim5, Dim0, i, j, k, l, m, l, m, n, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, j, k, l, m, n, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(l, m, n, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim0, l, m, n, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, i, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim3, m, i, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim5, Dim3, i, j, k, l, m, m, i, n, l, Dim1, Dim2, Dim4, Dim0, Dim5, Dim3, j, k, m, i, n, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(m, i, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim3, m, i, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, i, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim4, l, i, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim5, Dim4, i, j, k, l, m, l, i, n, m, Dim1, Dim2, Dim3, Dim0, Dim5, Dim4, j, k, l, i, n, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(l, i, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim0, Dim5, Dim4, l, i, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim3, i, m, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim5, Dim3, i, j, k, l, m, i, m, n, l, Dim1, Dim2, Dim0, Dim4, Dim5, Dim3, j, k, i, m, n, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(i, m, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim3, i, m, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, l, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim4, i, l, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim5, Dim4, i, j, k, l, m, i, l, n, m, Dim1, Dim2, Dim0, Dim3, Dim5, Dim4, j, k, i, l, n, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(i, l, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim3, Dim5, Dim4, i, l, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim0, m, n, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim3, Dim0, i, j, k, l, m, m, n, l, i, Dim1, Dim2, Dim4, Dim5, Dim3, Dim0, j, k, m, n, l, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(m, n, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim0, m, n, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim0, l, n, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim4, Dim0, i, j, k, l, m, l, n, m, i, Dim1, Dim2, Dim3, Dim5, Dim4, Dim0, j, k, l, n, m, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(l, n, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim0, l, n, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim3, m, n, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim0, Dim3, i, j, k, l, m, m, n, i, l, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, j, k, m, n, i, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(m, n, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim3, m, n, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim4, l, n, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim0, Dim4, i, j, k, l, m, l, n, i, m, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, j, k, l, n, i, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(l, n, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim0, Dim4, l, n, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim3, i, n, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim4, Dim3, i, j, k, l, m, i, n, m, l, Dim1, Dim2, Dim0, Dim5, Dim4, Dim3, j, k, i, n, m, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(i, n, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim4, Dim3, i, n, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(i, n, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim4, i, n, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim3, Dim4, i, j, k, l, m, i, n, l, m, Dim1, Dim2, Dim0, Dim5, Dim3, Dim4, j, k, i, n, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(i, n, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim0, Dim5, Dim3, Dim4, i, n, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, l, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim0, n, m, l, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim0, i, j, k, l, m, n, m, l, i, Dim1, Dim2, Dim5, Dim4, Dim3, Dim0, j, k, n, m, l, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(n, m, l, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim0, n, m, l, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, m, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim0, n, l, m, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim0, i, j, k, l, m, n, l, m, i, Dim1, Dim2, Dim5, Dim3, Dim4, Dim0, j, k, n, l, m, i>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(n, l, m, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim0, n, l, m, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, i, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim3, n, m, i, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, Dim3, i, j, k, l, m, n, m, i, l, Dim1, Dim2, Dim5, Dim4, Dim0, Dim3, j, k, n, m, i, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(n, m, i, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim0, Dim3, n, m, i, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, i, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim4, n, l, i, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, Dim4, i, j, k, l, m, n, l, i, m, Dim1, Dim2, Dim5, Dim3, Dim0, Dim4, j, k, n, l, i, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(n, l, i, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim0, Dim4, n, l, i, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim3, n, i, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, Dim3, i, j, k, l, m, n, i, m, l, Dim1, Dim2, Dim5, Dim0, Dim4, Dim3, j, k, n, i, m, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(n, i, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim4, Dim3, n, i, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim4, n, i, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, Dim4, i, j, k, l, m, n, i, l, m, Dim1, Dim2, Dim5, Dim0, Dim3, Dim4, j, k, n, i, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim5, j, k, n>(TensorExpr(a, b));
}

//B(n, i, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim0, Dim3, Dim4, n, i, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim5, l, k, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim1, Dim5, i, j, k, l, m, l, k, j, n, Dim0, Dim4, Dim3, Dim2, Dim1, Dim5, i, m, l, k, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(l, k, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim5, l, k, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim5, k, l, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim1, Dim5, i, j, k, l, m, k, l, j, n, Dim0, Dim4, Dim2, Dim3, Dim1, Dim5, i, m, k, l, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(k, l, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim5, k, l, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim5, l, j, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim2, Dim5, i, j, k, l, m, l, j, k, n, Dim0, Dim4, Dim3, Dim1, Dim2, Dim5, i, m, l, j, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(l, j, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim5, l, j, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim5, k, j, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim3, Dim5, i, j, k, l, m, k, j, l, n, Dim0, Dim4, Dim2, Dim1, Dim3, Dim5, i, m, k, j, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(k, j, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim5, k, j, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim5, j, l, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim2, Dim5, i, j, k, l, m, j, l, k, n, Dim0, Dim4, Dim1, Dim3, Dim2, Dim5, i, m, j, l, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(j, l, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim5, j, l, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim5, j, k, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim3, Dim5, i, j, k, l, m, j, k, l, n, Dim0, Dim4, Dim1, Dim2, Dim3, Dim5, i, m, j, k, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(j, k, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim5, j, k, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim1, l, k, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, Dim1, i, j, k, l, m, l, k, n, j, Dim0, Dim4, Dim3, Dim2, Dim5, Dim1, i, m, l, k, n, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(l, k, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim1, l, k, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim1, k, l, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, Dim1, i, j, k, l, m, k, l, n, j, Dim0, Dim4, Dim2, Dim3, Dim5, Dim1, i, m, k, l, n, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(k, l, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim1, k, l, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim2, l, j, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, Dim2, i, j, k, l, m, l, j, n, k, Dim0, Dim4, Dim3, Dim1, Dim5, Dim2, i, m, l, j, n, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(l, j, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim2, l, j, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim3, k, j, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, Dim3, i, j, k, l, m, k, j, n, l, Dim0, Dim4, Dim2, Dim1, Dim5, Dim3, i, m, k, j, n, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(k, j, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim3, k, j, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim2, j, l, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, Dim2, i, j, k, l, m, j, l, n, k, Dim0, Dim4, Dim1, Dim3, Dim5, Dim2, i, m, j, l, n, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(j, l, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim2, j, l, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim3, j, k, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, Dim3, i, j, k, l, m, j, k, n, l, Dim0, Dim4, Dim1, Dim2, Dim5, Dim3, i, m, j, k, n, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(j, k, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim3, j, k, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim1, l, n, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, Dim1, i, j, k, l, m, l, n, k, j, Dim0, Dim4, Dim3, Dim5, Dim2, Dim1, i, m, l, n, k, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(l, n, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim1, l, n, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim1, k, n, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, Dim1, i, j, k, l, m, k, n, l, j, Dim0, Dim4, Dim2, Dim5, Dim3, Dim1, i, m, k, n, l, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(k, n, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim1, k, n, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim2, l, n, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, Dim2, i, j, k, l, m, l, n, j, k, Dim0, Dim4, Dim3, Dim5, Dim1, Dim2, i, m, l, n, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(l, n, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim2, l, n, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim3, k, n, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, Dim3, i, j, k, l, m, k, n, j, l, Dim0, Dim4, Dim2, Dim5, Dim1, Dim3, i, m, k, n, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(k, n, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim3, k, n, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim2, j, n, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, Dim2, i, j, k, l, m, j, n, l, k, Dim0, Dim4, Dim1, Dim5, Dim3, Dim2, i, m, j, n, l, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(j, n, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim2, j, n, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim3, j, n, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, Dim3, i, j, k, l, m, j, n, k, l, Dim0, Dim4, Dim1, Dim5, Dim2, Dim3, i, m, j, n, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(j, n, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim3, j, n, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim1, n, l, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim1, i, j, k, l, m, n, l, k, j, Dim0, Dim4, Dim5, Dim3, Dim2, Dim1, i, m, n, l, k, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(n, l, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim1, n, l, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim1, n, k, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim1, i, j, k, l, m, n, k, l, j, Dim0, Dim4, Dim5, Dim2, Dim3, Dim1, i, m, n, k, l, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(n, k, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim1, n, k, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim2, n, l, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim2, i, j, k, l, m, n, l, j, k, Dim0, Dim4, Dim5, Dim3, Dim1, Dim2, i, m, n, l, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(n, l, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim2, n, l, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim3, n, k, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim3, i, j, k, l, m, n, k, j, l, Dim0, Dim4, Dim5, Dim2, Dim1, Dim3, i, m, n, k, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(n, k, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim3, n, k, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim2, n, j, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim2, i, j, k, l, m, n, j, l, k, Dim0, Dim4, Dim5, Dim1, Dim3, Dim2, i, m, n, j, l, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(n, j, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim2, n, j, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim3, n, j, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim3, i, j, k, l, m, n, j, k, l, Dim0, Dim4, Dim5, Dim1, Dim2, Dim3, i, m, n, j, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim4, Dim5, i, m, n>(TensorExpr(a, b));
}

//B(n, j, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim3, n, j, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim5, m, k, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim1, Dim5, i, j, k, l, m, m, k, j, n, Dim0, Dim3, Dim4, Dim2, Dim1, Dim5, i, l, m, k, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(m, k, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim1, Dim5, m, k, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim5, k, m, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim1, Dim5, i, j, k, l, m, k, m, j, n, Dim0, Dim3, Dim2, Dim4, Dim1, Dim5, i, l, k, m, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(k, m, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim1, Dim5, k, m, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim5, m, j, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim2, Dim5, i, j, k, l, m, m, j, k, n, Dim0, Dim3, Dim4, Dim1, Dim2, Dim5, i, l, m, j, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(m, j, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim2, Dim5, m, j, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim5, k, j, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim4, Dim5, i, j, k, l, m, k, j, m, n, Dim0, Dim3, Dim2, Dim1, Dim4, Dim5, i, l, k, j, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(k, j, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim4, Dim5, k, j, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim5, j, m, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim2, Dim5, i, j, k, l, m, j, m, k, n, Dim0, Dim3, Dim1, Dim4, Dim2, Dim5, i, l, j, m, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(j, m, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim2, Dim5, j, m, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim5, j, k, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim4, Dim5, i, j, k, l, m, j, k, m, n, Dim0, Dim3, Dim1, Dim2, Dim4, Dim5, i, l, j, k, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(j, k, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim4, Dim5, j, k, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim1, m, k, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim5, Dim1, i, j, k, l, m, m, k, n, j, Dim0, Dim3, Dim4, Dim2, Dim5, Dim1, i, l, m, k, n, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(m, k, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim1, m, k, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim1, k, m, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim5, Dim1, i, j, k, l, m, k, m, n, j, Dim0, Dim3, Dim2, Dim4, Dim5, Dim1, i, l, k, m, n, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(k, m, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim1, k, m, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim2, m, j, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim5, Dim2, i, j, k, l, m, m, j, n, k, Dim0, Dim3, Dim4, Dim1, Dim5, Dim2, i, l, m, j, n, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(m, j, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim2, m, j, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, j, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim4, k, j, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim5, Dim4, i, j, k, l, m, k, j, n, m, Dim0, Dim3, Dim2, Dim1, Dim5, Dim4, i, l, k, j, n, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(k, j, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim1, Dim5, Dim4, k, j, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim2, j, m, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim5, Dim2, i, j, k, l, m, j, m, n, k, Dim0, Dim3, Dim1, Dim4, Dim5, Dim2, i, l, j, m, n, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(j, m, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim2, j, m, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, k, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim4, j, k, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim5, Dim4, i, j, k, l, m, j, k, n, m, Dim0, Dim3, Dim1, Dim2, Dim5, Dim4, i, l, j, k, n, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(j, k, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim2, Dim5, Dim4, j, k, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim1, m, n, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim2, Dim1, i, j, k, l, m, m, n, k, j, Dim0, Dim3, Dim4, Dim5, Dim2, Dim1, i, l, m, n, k, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(m, n, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim1, m, n, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim1, k, n, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim4, Dim1, i, j, k, l, m, k, n, m, j, Dim0, Dim3, Dim2, Dim5, Dim4, Dim1, i, l, k, n, m, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(k, n, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim1, k, n, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim2, m, n, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim1, Dim2, i, j, k, l, m, m, n, j, k, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, i, l, m, n, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(m, n, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim2, m, n, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim4, k, n, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim1, Dim4, i, j, k, l, m, k, n, j, m, Dim0, Dim3, Dim2, Dim5, Dim1, Dim4, i, l, k, n, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(k, n, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim1, Dim4, k, n, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim2, j, n, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim4, Dim2, i, j, k, l, m, j, n, m, k, Dim0, Dim3, Dim1, Dim5, Dim4, Dim2, i, l, j, n, m, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(j, n, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim2, j, n, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim4, j, n, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim2, Dim4, i, j, k, l, m, j, n, k, m, Dim0, Dim3, Dim1, Dim5, Dim2, Dim4, i, l, j, n, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(j, n, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim2, Dim4, j, n, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, k, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim1, n, m, k, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim1, i, j, k, l, m, n, m, k, j, Dim0, Dim3, Dim5, Dim4, Dim2, Dim1, i, l, n, m, k, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(n, m, k, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim1, n, m, k, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim1, n, k, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim1, i, j, k, l, m, n, k, m, j, Dim0, Dim3, Dim5, Dim2, Dim4, Dim1, i, l, n, k, m, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(n, k, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim1, n, k, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, j, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim2, n, m, j, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim2, i, j, k, l, m, n, m, j, k, Dim0, Dim3, Dim5, Dim4, Dim1, Dim2, i, l, n, m, j, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(n, m, j, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim2, n, m, j, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim4, n, k, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, Dim4, i, j, k, l, m, n, k, j, m, Dim0, Dim3, Dim5, Dim2, Dim1, Dim4, i, l, n, k, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(n, k, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim1, Dim4, n, k, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim2, n, j, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim2, i, j, k, l, m, n, j, m, k, Dim0, Dim3, Dim5, Dim1, Dim4, Dim2, i, l, n, j, m, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(n, j, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim2, n, j, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim4, n, j, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, Dim4, i, j, k, l, m, n, j, k, m, Dim0, Dim3, Dim5, Dim1, Dim2, Dim4, i, l, n, j, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim3, Dim5, i, l, n>(TensorExpr(a, b));
}

//B(n, j, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim2, Dim4, n, j, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim5, m, l, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim1, Dim5, i, j, k, l, m, m, l, j, n, Dim0, Dim2, Dim4, Dim3, Dim1, Dim5, i, k, m, l, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(m, l, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim1, Dim5, m, l, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim5, l, m, j, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim1, Dim5, i, j, k, l, m, l, m, j, n, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, i, k, l, m, j, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(l, m, j, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim1, Dim5, l, m, j, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim5, m, j, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim3, Dim5, i, j, k, l, m, m, j, l, n, Dim0, Dim2, Dim4, Dim1, Dim3, Dim5, i, k, m, j, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(m, j, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim3, Dim5, m, j, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim5, l, j, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim4, Dim5, i, j, k, l, m, l, j, m, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(l, j, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim4, Dim5, l, j, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim5, j, m, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim3, Dim5, i, j, k, l, m, j, m, l, n, Dim0, Dim2, Dim1, Dim4, Dim3, Dim5, i, k, j, m, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(j, m, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim3, Dim5, j, m, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim5, j, l, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim4, Dim5, i, j, k, l, m, j, l, m, n, Dim0, Dim2, Dim1, Dim3, Dim4, Dim5, i, k, j, l, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(j, l, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim4, Dim5, j, l, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim1, m, l, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim5, Dim1, i, j, k, l, m, m, l, n, j, Dim0, Dim2, Dim4, Dim3, Dim5, Dim1, i, k, m, l, n, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(m, l, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim1, m, l, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim1, l, m, n, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim5, Dim1, i, j, k, l, m, l, m, n, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, i, k, l, m, n, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(l, m, n, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim1, l, m, n, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, j, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim3, m, j, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim5, Dim3, i, j, k, l, m, m, j, n, l, Dim0, Dim2, Dim4, Dim1, Dim5, Dim3, i, k, m, j, n, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(m, j, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim3, m, j, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, j, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim4, l, j, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim5, Dim4, i, j, k, l, m, l, j, n, m, Dim0, Dim2, Dim3, Dim1, Dim5, Dim4, i, k, l, j, n, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(l, j, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim1, Dim5, Dim4, l, j, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim3, j, m, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim5, Dim3, i, j, k, l, m, j, m, n, l, Dim0, Dim2, Dim1, Dim4, Dim5, Dim3, i, k, j, m, n, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(j, m, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim3, j, m, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, l, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim4, j, l, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim5, Dim4, i, j, k, l, m, j, l, n, m, Dim0, Dim2, Dim1, Dim3, Dim5, Dim4, i, k, j, l, n, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(j, l, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim3, Dim5, Dim4, j, l, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim1, m, n, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim3, Dim1, i, j, k, l, m, m, n, l, j, Dim0, Dim2, Dim4, Dim5, Dim3, Dim1, i, k, m, n, l, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(m, n, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim1, m, n, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim1, l, n, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim4, Dim1, i, j, k, l, m, l, n, m, j, Dim0, Dim2, Dim3, Dim5, Dim4, Dim1, i, k, l, n, m, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(l, n, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim1, l, n, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim3, m, n, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim1, Dim3, i, j, k, l, m, m, n, j, l, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, i, k, m, n, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(m, n, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim3, m, n, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim4, l, n, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim1, Dim4, i, j, k, l, m, l, n, j, m, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, i, k, l, n, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(l, n, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim1, Dim4, l, n, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim3, j, n, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim4, Dim3, i, j, k, l, m, j, n, m, l, Dim0, Dim2, Dim1, Dim5, Dim4, Dim3, i, k, j, n, m, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(j, n, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim4, Dim3, j, n, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim4, j, n, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim3, Dim4, i, j, k, l, m, j, n, l, m, Dim0, Dim2, Dim1, Dim5, Dim3, Dim4, i, k, j, n, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(j, n, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim1, Dim5, Dim3, Dim4, j, n, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, l, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim1, n, m, l, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim1, i, j, k, l, m, n, m, l, j, Dim0, Dim2, Dim5, Dim4, Dim3, Dim1, i, k, n, m, l, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(n, m, l, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim1, n, m, l, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, m, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim1, n, l, m, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim1, i, j, k, l, m, n, l, m, j, Dim0, Dim2, Dim5, Dim3, Dim4, Dim1, i, k, n, l, m, j>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(n, l, m, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim1, n, l, m, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, j, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim3, n, m, j, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, Dim3, i, j, k, l, m, n, m, j, l, Dim0, Dim2, Dim5, Dim4, Dim1, Dim3, i, k, n, m, j, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(n, m, j, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim1, Dim3, n, m, j, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, j, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim4, n, l, j, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, Dim4, i, j, k, l, m, n, l, j, m, Dim0, Dim2, Dim5, Dim3, Dim1, Dim4, i, k, n, l, j, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(n, l, j, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim1, Dim4, n, l, j, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim3, n, j, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, Dim3, i, j, k, l, m, n, j, m, l, Dim0, Dim2, Dim5, Dim1, Dim4, Dim3, i, k, n, j, m, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(n, j, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim4, Dim3, n, j, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim4, n, j, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, Dim4, i, j, k, l, m, n, j, l, m, Dim0, Dim2, Dim5, Dim1, Dim3, Dim4, i, k, n, j, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim5, i, k, n>(TensorExpr(a, b));
}

//B(n, j, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim1, Dim3, Dim4, n, j, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim5, m, l, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim2, Dim5, i, j, k, l, m, m, l, k, n, Dim0, Dim1, Dim4, Dim3, Dim2, Dim5, i, j, m, l, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(m, l, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim2, Dim5, m, l, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim5, l, m, k, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim2, Dim5, i, j, k, l, m, l, m, k, n, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, i, j, l, m, k, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(l, m, k, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim2, Dim5, l, m, k, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim5, m, k, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim3, Dim5, i, j, k, l, m, m, k, l, n, Dim0, Dim1, Dim4, Dim2, Dim3, Dim5, i, j, m, k, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(m, k, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim3, Dim5, m, k, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim5, l, k, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim4, Dim5, i, j, k, l, m, l, k, m, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(l, k, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim4, Dim5, l, k, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim5, k, m, l, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim3, Dim5, i, j, k, l, m, k, m, l, n, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, i, j, k, m, l, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(k, m, l, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim3, Dim5, k, m, l, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim5, k, l, m, n> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, k, l, m, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(k, l, m, n, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim4, Dim5, k, l, m, n> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, l, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim2, m, l, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim5, Dim2, i, j, k, l, m, m, l, n, k, Dim0, Dim1, Dim4, Dim3, Dim5, Dim2, i, j, m, l, n, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(m, l, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim2, m, l, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim2, l, m, n, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim5, Dim2, i, j, k, l, m, l, m, n, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, i, j, l, m, n, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(l, m, n, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim2, l, m, n, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, k, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim3, m, k, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim5, Dim3, i, j, k, l, m, m, k, n, l, Dim0, Dim1, Dim4, Dim2, Dim5, Dim3, i, j, m, k, n, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(m, k, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim3, m, k, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, k, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim4, l, k, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim5, Dim4, i, j, k, l, m, l, k, n, m, Dim0, Dim1, Dim3, Dim2, Dim5, Dim4, i, j, l, k, n, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(l, k, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim2, Dim5, Dim4, l, k, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim3, k, m, n, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim5, Dim3, i, j, k, l, m, k, m, n, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, i, j, k, m, n, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(k, m, n, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim3, k, m, n, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, l, n, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim4, k, l, n, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim5, Dim4, i, j, k, l, m, k, l, n, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, i, j, k, l, n, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(k, l, n, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim3, Dim5, Dim4, k, l, n, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim2, m, n, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim3, Dim2, i, j, k, l, m, m, n, l, k, Dim0, Dim1, Dim4, Dim5, Dim3, Dim2, i, j, m, n, l, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(m, n, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim2, m, n, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim2, l, n, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim4, Dim2, i, j, k, l, m, l, n, m, k, Dim0, Dim1, Dim3, Dim5, Dim4, Dim2, i, j, l, n, m, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(l, n, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim4, Dim2, l, n, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim3, m, n, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim2, Dim3, i, j, k, l, m, m, n, k, l, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, i, j, m, n, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(m, n, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim3, m, n, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim4, l, n, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim2, Dim4, i, j, k, l, m, l, n, k, m, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, i, j, l, n, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(l, n, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim3, Dim5, Dim2, Dim4, l, n, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim3, k, n, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim4, Dim3, i, j, k, l, m, k, n, m, l, Dim0, Dim1, Dim2, Dim5, Dim4, Dim3, i, j, k, n, m, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(k, n, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim4, Dim3, k, n, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim4, k, n, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim3, Dim4, i, j, k, l, m, k, n, l, m, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, i, j, k, n, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(k, n, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim2, Dim5, Dim3, Dim4, k, n, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, l, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim2, n, m, l, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, Dim2, i, j, k, l, m, n, m, l, k, Dim0, Dim1, Dim5, Dim4, Dim3, Dim2, i, j, n, m, l, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(n, m, l, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim3, Dim2, n, m, l, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, m, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim2, n, l, m, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, Dim2, i, j, k, l, m, n, l, m, k, Dim0, Dim1, Dim5, Dim3, Dim4, Dim2, i, j, n, l, m, k>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(n, l, m, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim4, Dim2, n, l, m, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim3, n, m, k, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, Dim3, i, j, k, l, m, n, m, k, l, Dim0, Dim1, Dim5, Dim4, Dim2, Dim3, i, j, n, m, k, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(n, m, k, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim4, Dim2, Dim3, n, m, k, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, k, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim4, n, l, k, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, Dim4, i, j, k, l, m, n, l, k, m, Dim0, Dim1, Dim5, Dim3, Dim2, Dim4, i, j, n, l, k, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(n, l, k, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim3, Dim2, Dim4, n, l, k, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, m, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim3, n, k, m, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, Dim3, i, j, k, l, m, n, k, m, l, Dim0, Dim1, Dim5, Dim2, Dim4, Dim3, i, j, n, k, m, l>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(n, k, m, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim4, Dim3, n, k, m, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim4, n, k, l, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor4_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, Dim4, i, j, k, l, m, n, k, l, m, Dim0, Dim1, Dim5, Dim2, Dim3, Dim4, i, j, n, k, l, m>;
  return Tensor3_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim5, i, j, n>(TensorExpr(a, b));
}

//B(n, k, l, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<B, U, Dim5, Dim2, Dim3, Dim4, n, k, l, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

} //namespace FTensor

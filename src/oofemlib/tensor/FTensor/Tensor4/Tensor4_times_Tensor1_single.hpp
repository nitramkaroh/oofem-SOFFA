#pragma once

namespace FTensor
{

  // A(i,j,k,l) * B(m,n) double contraction
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l,
            char m, int DimA, int DimB, int DimC, int DimX, char a,
            char b, char c, char x>
  class Tensor4_times_Tensor1_double
  {
    Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    Tensor1_Expr<B, U, Dim4, m> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2, const int N3) const
    {
      typename promote<T, U>::V result(0);
      for(int xx = 0; xx < DimX; ++xx)
	{
            // Permutation is where the indices get checked.
            result += Permutation4<DimA, DimB, DimC, DimX, a, b, c, x>().eval(
                        iterA, N1, N2, N3, xx)
	              * iterB(xx);
	}
      return result;
    }

    Tensor4_times_Tensor1_double(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &iter_a,
      const Tensor1_Expr<B, U, Dim4, m> &iter_b)
        : iterA(iter_a), iterB(iter_b)
    {}
  };


  
  // A(i,j,k,l)*B(l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim3, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor1_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3,
				     i, j, k, l, l, Dim0, Dim1, Dim2,
                                     Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, i,
                        j, k>(TensorExpr(a, b));
  }

  // B(l) * A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor1_Expr<B, U, Dim3, l> &b,
	    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  

   // A(i,j,k,l)*B(k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim2, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor1_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2,
				     i, j, k, l, k, Dim0, Dim1, Dim2,
                                     Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, i,
                        j, l>(TensorExpr(a, b));
  }


  // B(l) * A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor1_Expr<B, U, Dim2, k> &b,
	    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }


   // A(i,j,k,l)*B(j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim1, j> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor1_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1,
				     i, j, k, l, j, Dim0, Dim1, Dim2,
                                     Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, i,
                        k, l>(TensorExpr(a, b));
  }



  // B(l) * A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor1_Expr<B, U, Dim1, j> &b,
	    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }



  
   // A(i,j,k,l)*B(i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim0, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor1_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0,
				     i, j, k, l, i, Dim0, Dim1, Dim2,
                                     Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, j,
                        k, l>(TensorExpr(a, b));
  }


  // B(i) * A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor1_Expr<B, U, Dim0, i> &b,
	    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }
  
  
  
}

#pragma once

namespace FTensor
{
  // A(i,j,k,l)*B(l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_times_Tensor1_single4
  {
    Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    Tensor1_Expr<B, U, Dim3, l> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<Current_Dim> &) const
    {
      return iterA(N1, N2, N3, Current_Dim - 1) * iterB(Current_Dim - 1)
	+ eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, N2, N3, 0) * iterB(0);
    }

  public:
    Tensor4_times_Tensor1_single4(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor1_Expr<B, U, Dim3, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim3>());
    }
    
  };

 template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3,
            char i, char j, char k, char l>
 Tensor3_Expr<Tensor4_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3,  i, j, k,l>, typename promote<T, U>::V, Dim0, Dim1, Dim2, i, j, k>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim3, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, i, j, k>(TensorExpr(a, b));
  }
  //////////////////////////////////////////////////////////////////////////////////////
  // B(l) * A(i,j,k,l)
 template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3,
            char i, char j, char k, char l>
 Tensor3_Expr<Tensor4_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3,  i, j, k,l>, typename promote<T, U>::V, Dim0, Dim1, Dim2, i, j, k>
 operator*(const Tensor1_Expr<B, U, Dim3, l> &b,
	    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    using TensorExpr
      = Tensor4_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, i, j, k>(TensorExpr(a, b));
  }
  
  /////////////////////////////////////////////////////////////////////////////////////
  // A(i,j,k,l)*B(k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_times_Tensor1_single3
  {
    Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    Tensor1_Expr<B, U, Dim2, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<Current_Dim> &) const
    {
      return iterA(N1, N2, Current_Dim - 1,  N3) * iterB(Current_Dim - 1)
	+ eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, N2, 0, N3) * iterB(0);
    }

  public:
    Tensor4_times_Tensor1_single3(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor1_Expr<B, U, Dim2, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim2>());
    }
    
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3,
            char i, char j, char k, char l>
  Tensor3_Expr<Tensor4_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3,  i, j, k, l>, typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j, l>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim2, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j, l>(TensorExpr(a, b));
  }
  // B(k) * A(i,j,k,l)
    template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3,
            char i, char j, char k, char l>
  Tensor3_Expr<Tensor4_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3,  i, j, k, l>, typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j, l>
  operator*(const Tensor1_Expr<B, U, Dim2, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    using TensorExpr
      = Tensor4_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j, l>(TensorExpr(a, b));
  }
  //////////////////////////////////////////////////////////////////////////////////////
    // A(i,j,k,l)*B(j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_times_Tensor1_single2
  {
    Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    Tensor1_Expr<B, U, Dim1, j> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1, N2,  N3) * iterB(Current_Dim - 1)
	+ eval(N1,  N2,  N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, 0, N2, N3) * iterB(0);
    }

  public:
    Tensor4_times_Tensor1_single2(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor1_Expr<B, U, Dim1, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim1>());
    }
    
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3,
            char i, char j, char k, char l>
  Tensor3_Expr<Tensor4_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3,  i, j, k, l>, typename promote<T, U>::V, Dim0, Dim1, Dim2, i, k, l>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim1, j> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, i, k, l>(TensorExpr(a, b));
  }
  // A(i,j,k,l)*B(j)
    template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3,
            char i, char j, char k, char l>
  Tensor3_Expr<Tensor4_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3,  i, j, k, l>, typename promote<T, U>::V, Dim0, Dim1, Dim2, i, k, l>
  operator*(const Tensor1_Expr<B, U, Dim1, j> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    using TensorExpr
      = Tensor4_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, i, k, l>(TensorExpr(a, b));
  }
    //////////////////////////////////////////////////////////////////////////////////////
    // A(i,j,k,l)*B(i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor4_times_Tensor1_single1
  {
    Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    Tensor1_Expr<B, U, Dim0, i> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<Current_Dim> &) const
    {
      return iterA(Current_Dim - 1, N1, N2,  N3) * iterB(Current_Dim - 1)
	+ eval(N1,  N2,  N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(0, N1, N2, N3) * iterB(0);
    }

  public:
    Tensor4_times_Tensor1_single1(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
      const Tensor1_Expr<B, U, Dim0, i> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim0>());
    }
    
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3,
            char i, char j, char k, char l>
  Tensor3_Expr<Tensor4_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3,  i, j, k, l>, typename promote<T, U>::V, Dim0, Dim1, Dim2, j, k, l>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim0, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, j, k, l>(TensorExpr(a, b));
  }
  // *B(i) *A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3,
            char i, char j, char k, char l>
  Tensor3_Expr<Tensor4_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3,  i, j, k, l>, typename promote<T, U>::V, Dim0, Dim1, Dim2, j, k, l>
  operator*(const Tensor1_Expr<B, U, Dim0, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    using TensorExpr
      = Tensor4_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, j, k, l>(TensorExpr(a, b));
  }
  //////////////////////////////////////////////

  
  
}

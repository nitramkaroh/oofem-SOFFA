#pragma once

namespace FTensor
{
////////////////////////// A(i, j, k, l, m, n, o, )*B(i, ) ////////////////////////// 
// Base class:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
class Tensor7_times_Tensor1_single1
{
  Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
  Tensor1_Expr<B, U, Dim0, i> iterB;

  template<int Current_Dim>
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<Current_Dim> &) const
  {
    return iterA(Current_Dim - 1, N1, N2, N3, N4, N5, N6) * iterB(Current_Dim - 1) + eval(N1, N2, N3, N4, N5, N6, Number<Current_Dim - 1>());
  }
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<1> &) const
  {
    return iterA(0, N1, N2, N3, N4, N5, N6) * iterB(0);
  }

public:
  Tensor7_times_Tensor1_single1(
    const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
    const Tensor1_Expr<B, U, Dim0, i> &b)
       : iterA(a), iterB(b)
  {}

  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    return eval(N1, N2, N3, N4, N5, N6, Number<Dim0>());
  }

};

// operator for A*b:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, j, k, l, m, n, o>
operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
          const Tensor1_Expr<B, U, Dim0, i> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, j, k, l, m, n, o>(TensorExpr(a, b));
}

// operator for b*A:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, j, k, l, m, n, o>
operator*(const Tensor1_Expr<B, U, Dim0, i> &b,
          const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, j, k, l, m, n, o>(TensorExpr(a, b));
}

////////////////////////// A(i, j, k, l, m, n, o, )*B(j, ) ////////////////////////// 
// Base class:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
class Tensor7_times_Tensor1_single2
{
  Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
  Tensor1_Expr<B, U, Dim1, j> iterB;

  template<int Current_Dim>
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<Current_Dim> &) const
  {
    return iterA(N1, Current_Dim - 1, N2, N3, N4, N5, N6) * iterB(Current_Dim - 1) + eval(N1, N2, N3, N4, N5, N6, Number<Current_Dim - 1>());
  }
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<1> &) const
  {
    return iterA(N1, 0, N2, N3, N4, N5, N6) * iterB(0);
  }

public:
  Tensor7_times_Tensor1_single2(
    const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
    const Tensor1_Expr<B, U, Dim1, j> &b)
       : iterA(a), iterB(b)
  {}

  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    return eval(N1, N2, N3, N4, N5, N6, Number<Dim1>());
  }

};

// operator for A*b:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, i, k, l, m, n, o>
operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
          const Tensor1_Expr<B, U, Dim1, j> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, i, k, l, m, n, o>(TensorExpr(a, b));
}

// operator for b*A:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, i, k, l, m, n, o>
operator*(const Tensor1_Expr<B, U, Dim1, j> &b,
          const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, i, k, l, m, n, o>(TensorExpr(a, b));
}

////////////////////////// A(i, j, k, l, m, n, o, )*B(k, ) ////////////////////////// 
// Base class:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
class Tensor7_times_Tensor1_single3
{
  Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
  Tensor1_Expr<B, U, Dim2, k> iterB;

  template<int Current_Dim>
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<Current_Dim> &) const
  {
    return iterA(N1, N2, Current_Dim - 1, N3, N4, N5, N6) * iterB(Current_Dim - 1) + eval(N1, N2, N3, N4, N5, N6, Number<Current_Dim - 1>());
  }
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<1> &) const
  {
    return iterA(N1, N2, 0, N3, N4, N5, N6) * iterB(0);
  }

public:
  Tensor7_times_Tensor1_single3(
    const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
    const Tensor1_Expr<B, U, Dim2, k> &b)
       : iterA(a), iterB(b)
  {}

  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    return eval(N1, N2, N3, N4, N5, N6, Number<Dim2>());
  }

};

// operator for A*b:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, i, j, l, m, n, o>
operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
          const Tensor1_Expr<B, U, Dim2, k> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, i, j, l, m, n, o>(TensorExpr(a, b));
}

// operator for b*A:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, i, j, l, m, n, o>
operator*(const Tensor1_Expr<B, U, Dim2, k> &b,
          const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, i, j, l, m, n, o>(TensorExpr(a, b));
}

////////////////////////// A(i, j, k, l, m, n, o, )*B(l, ) ////////////////////////// 
// Base class:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
class Tensor7_times_Tensor1_single4
{
  Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
  Tensor1_Expr<B, U, Dim3, l> iterB;

  template<int Current_Dim>
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<Current_Dim> &) const
  {
    return iterA(N1, N2, N3, Current_Dim - 1, N4, N5, N6) * iterB(Current_Dim - 1) + eval(N1, N2, N3, N4, N5, N6, Number<Current_Dim - 1>());
  }
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<1> &) const
  {
    return iterA(N1, N2, N3, 0, N4, N5, N6) * iterB(0);
  }

public:
  Tensor7_times_Tensor1_single4(
    const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
    const Tensor1_Expr<B, U, Dim3, l> &b)
       : iterA(a), iterB(b)
  {}

  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    return eval(N1, N2, N3, N4, N5, N6, Number<Dim3>());
  }

};

// operator for A*b:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, i, j, k, m, n, o>
operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
          const Tensor1_Expr<B, U, Dim3, l> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, i, j, k, m, n, o>(TensorExpr(a, b));
}

// operator for b*A:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, i, j, k, m, n, o>
operator*(const Tensor1_Expr<B, U, Dim3, l> &b,
          const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, i, j, k, m, n, o>(TensorExpr(a, b));
}

////////////////////////// A(i, j, k, l, m, n, o, )*B(m, ) ////////////////////////// 
// Base class:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
class Tensor7_times_Tensor1_single5
{
  Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
  Tensor1_Expr<B, U, Dim4, m> iterB;

  template<int Current_Dim>
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<Current_Dim> &) const
  {
    return iterA(N1, N2, N3, N4, Current_Dim - 1, N5, N6) * iterB(Current_Dim - 1) + eval(N1, N2, N3, N4, N5, N6, Number<Current_Dim - 1>());
  }
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<1> &) const
  {
    return iterA(N1, N2, N3, N4, 0, N5, N6) * iterB(0);
  }

public:
  Tensor7_times_Tensor1_single5(
    const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
    const Tensor1_Expr<B, U, Dim4, m> &b)
       : iterA(a), iterB(b)
  {}

  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    return eval(N1, N2, N3, N4, N5, N6, Number<Dim4>());
  }

};

// operator for A*b:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single5<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, i, j, k, l, n, o>
operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
          const Tensor1_Expr<B, U, Dim4, m> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single5<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, i, j, k, l, n, o>(TensorExpr(a, b));
}

// operator for b*A:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single5<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, i, j, k, l, n, o>
operator*(const Tensor1_Expr<B, U, Dim4, m> &b,
          const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single5<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, i, j, k, l, n, o>(TensorExpr(a, b));
}

////////////////////////// A(i, j, k, l, m, n, o, )*B(n, ) ////////////////////////// 
// Base class:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
class Tensor7_times_Tensor1_single6
{
  Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
  Tensor1_Expr<B, U, Dim5, n> iterB;

  template<int Current_Dim>
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<Current_Dim> &) const
  {
    return iterA(N1, N2, N3, N4, N5, Current_Dim - 1, N6) * iterB(Current_Dim - 1) + eval(N1, N2, N3, N4, N5, N6, Number<Current_Dim - 1>());
  }
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<1> &) const
  {
    return iterA(N1, N2, N3, N4, N5, 0, N6) * iterB(0);
  }

public:
  Tensor7_times_Tensor1_single6(
    const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
    const Tensor1_Expr<B, U, Dim5, n> &b)
       : iterA(a), iterB(b)
  {}

  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    return eval(N1, N2, N3, N4, N5, N6, Number<Dim5>());
  }

};

// operator for A*b:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single6<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, i, j, k, l, m, o>
operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
          const Tensor1_Expr<B, U, Dim5, n> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single6<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, i, j, k, l, m, o>(TensorExpr(a, b));
}

// operator for b*A:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single6<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, i, j, k, l, m, o>
operator*(const Tensor1_Expr<B, U, Dim5, n> &b,
          const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single6<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim6, i, j, k, l, m, o>(TensorExpr(a, b));
}

////////////////////////// A(i, j, k, l, m, n, o, )*B(o, ) ////////////////////////// 
// Base class:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
class Tensor7_times_Tensor1_single7
{
  Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> iterA;
  Tensor1_Expr<B, U, Dim6, o> iterB;

  template<int Current_Dim>
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<Current_Dim> &) const
  {
    return iterA(N1, N2, N3, N4, N5, N6, Current_Dim - 1) * iterB(Current_Dim - 1) + eval(N1, N2, N3, N4, N5, N6, Number<Current_Dim - 1>());
  }
  typename promote<T, U>::V
  eval(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const Number<1> &) const
  {
    return iterA(N1, N2, N3, N4, N5, N6, 0) * iterB(0);
  }

public:
  Tensor7_times_Tensor1_single7(
    const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
    const Tensor1_Expr<B, U, Dim6, o> &b)
       : iterA(a), iterB(b)
  {}

  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    return eval(N1, N2, N3, N4, N5, N6, Number<Dim6>());
  }

};

// operator for A*b:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single7<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>
operator*(const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a,
          const Tensor1_Expr<B, U, Dim6, o> &b)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single7<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>(TensorExpr(a, b));
}

// operator for b*A:
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
Tensor6_Expr<Tensor7_times_Tensor1_single7<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>
operator*(const Tensor1_Expr<B, U, Dim6, o> &b,
          const Tensor7_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o> &a)
{
  using TensorExpr
    = Tensor7_times_Tensor1_single7<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>(TensorExpr(a, b));
}

} //namespace FTensor

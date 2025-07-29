// Subtracts Tensor7 - Tensor7 -> Tensor7

#pragma once

namespace FTensor
{
  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim2_0, int Dim3_0, int Dim4_0, int Dim5_0, int Dim6_0, int Dim0_1, int Dim1_1, int Dim2_1,
            int Dim3_1, int Dim4_1, int Dim5_1, int Dim6_1, char i0, char j0, char k0, char l0, char m0, char n0, char o0, char i1, char j1,
            char k1, char l1, char m1, char n1, char o1>
  class Tensor7_minus_Tensor7
  {
    Tensor7_Expr<A, T, Dim0_0, Dim1_0, Dim2_0, Dim3_0, Dim4_0, Dim5_0, Dim6_0, i0, j0, k0, l0, m0, n0, o0> iterA;
    Tensor7_Expr<B, U, Dim0_1, Dim1_1, Dim2_1, Dim3_1, Dim4_1, Dim5_1, Dim6_1, i1, j1, k1, l1, m1, n1, o1> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6, const int N7) const
    {
      return iterA(N1, N2, N3, N4, N5, N6, N7) - permute(iterA, iterB, N1, N2, N3, N4, N5, N6, N7);
    }

    Tensor7_minus_Tensor7(const Tensor7_Expr<A, T, Dim0_0, Dim1_0, Dim2_0,
                                             Dim3_0, Dim4_0, Dim5_0, Dim6_0, i0, j0, k0, l0, m0, n0, o0> &a,
                          const Tensor7_Expr<B, U, Dim0_1, Dim1_1, Dim2_1,
                                             Dim3_1, Dim4_1, Dim5_1, Dim6_1, i1, j1, k1, l1, m1, n1, o1> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim2_0, int Dim3_0, int Dim4_0, int Dim5_0, int Dim6_0, int Dim0_1, int Dim1_1, int Dim2_1,
            int Dim3_1, int Dim4_1, int Dim5_1, int Dim6_1, char i0, char j0, char k0, char l0, char m0, char n0, char o0, char i1, char j1,
            char k1, char l1, char m1, char n1, char o1>
  auto operator-(
    const Tensor7_Expr<A, T, Dim0_0, Dim1_0, Dim2_0, Dim3_0, Dim4_0, Dim5_0, Dim6_0, i0, j0, k0, l0, m0, n0, o0> &a,
    const Tensor7_Expr<B, U, Dim0_1, Dim1_1, Dim2_1, Dim3_1, Dim4_1, Dim5_1, Dim6_1, i1, j1, k1, l1, m1, n1, o1> &b)
  {
    using TensorExpr
      = Tensor7_minus_Tensor7<A, B, T, U, Dim0_0, Dim1_0, Dim2_0, Dim3_0, Dim4_0, Dim5_0, Dim6_0,
                              Dim0_1, Dim1_1, Dim2_1, Dim3_1, Dim4_1, Dim5_1, Dim6_1, i0, j0, k0, l0, m0, n0, o0,
                              i1, j1, k1, l1, m1, n1, o1>;
    return Tensor7_Expr<TensorExpr, typename promote<T, U>::V, Dim0_0, Dim1_0,
                        Dim2_0, Dim3_0, Dim4_0, Dim5_0, Dim6_0, i0, j0, k0, l0, m0, n0, o0>(TensorExpr(a, b));
  }
}

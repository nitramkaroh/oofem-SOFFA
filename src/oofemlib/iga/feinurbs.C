/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "feinurbs.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "iga.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace oofem {

void
NURBSInterpolation :: initializeFrom(InputRecord &ir)
{
    BSplineInterpolation :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, weights, _IFT_NURBSInterpolation_weights);

    int numberOfWeights = 1;
    for ( int i = 0; i < nsd; i++ ) {
        numberOfWeights *= numberOfControlPoints [ i ];
    }
    if ( weights.giveSize() != numberOfWeights ) {
        throw ValueInputException(ir, _IFT_NURBSInterpolation_weights, "the number of weights must equal the number of control points");
    }
    // the rational basis functions are divided by sum(N*w)
    for ( double w : weights ) {
        if ( !std::isfinite(w) || w <= 0. ) {
            throw ValueInputException(ir, _IFT_NURBSInterpolation_weights, "weights must be finite and positive");
        }
    }
}


void NURBSInterpolation :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    const FEIIGAElementGeometryWrapper &gw = static_cast< const FEIIGAElementGeometryWrapper& >(cellgeo);
    IntArray span(nsd);
    double sum = 0.0, val;
    int count, c = 1;
    std :: vector< FloatArray > N(nsd);

    if ( gw.knotSpan ) {
        span = * gw.knotSpan;
    } else {
        for ( int i = 0; i < nsd; i++ ) {
            span[i] = this->findSpan(numberOfControlPoints [ i ] - 1, degree [ i ], lcoords[i], knotVector [ i ]);
        }
    }

    for ( int i = 0; i < nsd; i++ ) {
        this->basisFuns(N [ i ], span[i], lcoords[i], degree [ i ], knotVector [ i ]);
    }

    count = giveNumberOfKnotSpanBasisFunctions(span);
    answer.resize(count);

    if ( nsd == 1 ) {
        int uind = span[0] - degree [ 0 ];
        int ind = uind + 1;
        for ( int k = 0; k <= degree [ 0 ]; k++ ) {
            answer.at(c++) = val = N [ 0 ][k] * weights.at(ind + k);       // Nu*w
            sum += val;
        }
    } else if ( nsd == 2 ) {
        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int l = 0; l <= degree [ 1 ]; l++ ) {
            for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                answer.at(c++) = val = N [ 0 ][k] * N [ 1 ][l] * weights.at(ind + k); // Nu*Nv*w
                sum += val;
            }

            ind += numberOfControlPoints [ 0 ];
        }
    } else if ( nsd == 3 ) {
        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int tind = span[2] - degree [ 2 ];
        int ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int m = 0; m <= degree [ 2 ]; m++ ) {
            int indx = ind;
            for ( int l = 0; l <= degree [ 1 ]; l++ ) {
                for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                    answer.at(c++) = val = N [ 0 ][k] * N [ 1 ][l] * N [ 2 ][m] * weights.at(ind + k);     // Nu*Nv*Nt*w
                    sum += val;
                }

                ind += numberOfControlPoints [ 0 ];
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];
        }
    } else {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }

    answer.times(1./sum);
}


double NURBSInterpolation :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    const FEIIGAElementGeometryWrapper &gw = static_cast< const FEIIGAElementGeometryWrapper& >(cellgeo);
    FloatMatrix jacobian(nsd, nsd);
    IntArray span(nsd);
    double Jacob = 0.;
    int count;
    std :: vector< FloatArray > N(nsd);
    std :: vector< FloatMatrix > ders(nsd);

    if ( gw.knotSpan ) {
        span = * gw.knotSpan;
    } else {
        for ( int i = 0; i < nsd; i++ ) {
            span[i] = this->findSpan(numberOfControlPoints [ i ] - 1, degree [ i ], lcoords[i], knotVector [ i ]);
        }
    }

    for ( int i = 0; i < nsd; i++ ) {
        this->dersBasisFuns(1, lcoords[i], span[i], degree [ i ], knotVector [ i ], ders [ i ]);
    }

    count = giveNumberOfKnotSpanBasisFunctions(span);
    answer.resize(count, nsd);

    std :: vector< FloatArray > Aders(nsd);
    FloatArray wders;          // 0th and 1st derivatives in w direction on BSpline

    for ( int i = 0; i < nsd; i++ ) {
        Aders [ i ].resize(nsd + 1);
        Aders [ i ].zero();
    }

    wders.resize(nsd + 1);
    wders.zero();

    if ( nsd == 1 ) {
        // calculate values and derivatives of nonrational Bspline curve with weights at first (Aders, wders)
        int uind = span[0] - degree [ 0 ];
        int ind = uind + 1;
        for ( int k = 0; k <= degree [ 0 ]; k++ ) {
            const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
            double w = this->weights.at(ind + k);

            Aders [ 0 ][0] += ders [ 0 ](0, k) * vertexCoords[0] * w;   // xw=sum(Nu*x*w)
            wders[0]    += ders [ 0 ](0, k) * w;                               // w=sum(Nu*w)

            Aders [ 0 ][1] += ders [ 0 ](1, k) * vertexCoords[0] * w;   // dxw/du=sum(dNu/du*x*w)
            wders[1]    += ders [ 0 ](1, k) * w;                               // dw/du=sum(dNu/du*w)
        }

        double weight = wders[0];

        // calculation of jacobian matrix according to Eq 4.7
        jacobian(0, 0) = ( Aders [ 0 ][1] - wders[1] * Aders [ 0 ][0] / weight ) / weight; // dx/du

        Jacob = jacobian.giveDeterminant();

        //calculation of derivatives of NURBS basis functions with respect to local parameters is not covered by NURBS book
        double product = Jacob * weight * weight;
        int cnt = 0;
        ind = uind + 1;
        for ( int k = 0; k <= degree [ 0 ]; k++ ) {
            double w = this->weights.at(ind + k);
            // [dNu/du*w*sum(Nu*w) - Nu*w*sum(dNu/du*w)] / [J*sum(Nu*w)^2]
            answer(cnt, 0) = ( ders [ 0 ](1, k) * w * weight - ders [ 0 ](0, k) * w * wders[1] ) / product;
            cnt++;
        }
    } else if ( nsd == 2 ) {
        FloatArray tmp1(nsd + 1), tmp2(nsd + 1);    // allow for weight

        // calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int l = 0; l <= degree [ 1 ]; l++ ) {
            tmp1.zero();
            tmp2.zero();
            for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
                double w = this->weights.at(ind + k);

                tmp1[0] += ders [ 0 ](0, k) * vertexCoords[0] * w; // sum(Nu*x*w)
                tmp1[1] += ders [ 0 ](0, k) * vertexCoords[1] * w; // sum(Nu*y*w)
                tmp1[2] += ders [ 0 ](0, k) * w;           // sum(Nu*w)

                tmp2[0] += ders [ 0 ](1, k) * vertexCoords[0] * w; // sum(dNu/du*x*w)
                tmp2[1] += ders [ 0 ](1, k) * vertexCoords[1] * w; // sum(dNu/du*y*w)
                tmp2[2] += ders [ 0 ](1, k) * w;           // sum(dNu/du*w)
            }

            ind += numberOfControlPoints [ 0 ];

            Aders [ 0 ][0] += ders [ 1 ](0, l) * tmp1[0]; // xw=sum(Nv*sum(Nu*x*w)
            Aders [ 1 ][0] += ders [ 1 ](0, l) * tmp1[1]; // yw=sum(Nv*sum(Nu*y*w)
            wders[0]    += ders [ 1 ](0, l) * tmp1[2]; // w=sum(Nv*sum(Nu*w)

            Aders [ 0 ][1] += ders [ 1 ](0, l) * tmp2[0]; // dxw/du=sum(Nv*sum(dNu/du*x*w)
            Aders [ 1 ][1] += ders [ 1 ](0, l) * tmp2[1]; // dyw/du=sum(Nv*sum(dNu/du*y*w)
            wders[1]    += ders [ 1 ](0, l) * tmp2[2];        // dw/du=sum(Nv*sum(dNu/du*w)

            Aders [ 0 ][2] += ders [ 1 ](1, l) * tmp1[0]; // dxw/dv=sum(dNv/dv*sum(Nu*x*w)
            Aders [ 1 ][2] += ders [ 1 ](1, l) * tmp1[1]; // dyw/dv=sum(dNv/dv*sum(Nu*y*w)
            wders[2]    += ders [ 1 ](1, l) * tmp1[2]; // dw/dv=sum(dNv/dv*sum(Nu*w)
        }

        double weight = wders[0];

        // calculation of jacobian matrix according to Eq 4.19
        tmp1[0] = Aders [ 0 ][0] / weight;
        tmp1[1] = Aders [ 1 ][0] / weight;
        jacobian(0, 0) = ( Aders [ 0 ][1] - wders[1] * tmp1[0] ) / weight; // dx/du
        jacobian(0, 1) = ( Aders [ 1 ][1] - wders[1] * tmp1[1] ) / weight; // dy/du
        jacobian(1, 0) = ( Aders [ 0 ][2] - wders[2] * tmp1[0] ) / weight; // dx/dv
        jacobian(1, 1) = ( Aders [ 1 ][2] - wders[2] * tmp1[1] ) / weight; // dy/dv

        Jacob = jacobian.giveDeterminant();

        //calculation of derivatives of NURBS basis functions with respect to local parameters is not covered by NURBS book
        double product = Jacob * weight * weight;
        int cnt = 0;
        ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int l = 0; l <= degree [ 1 ]; l++ ) {
            for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                double w = this->weights.at(ind + k);
                // dNu/du*Nv*w*sum(Nu*Nv*w) - Nu*Nv*w*sum(dNu/du*Nv*w)
                tmp1[0] = ders [ 0 ](1, k) * ders [ 1 ](0, l) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, l) * w * wders[1];
                // Nu*dNv/dv*w*sum(Nu*Nv*w) - Nu*Nv*w*sum(Nu*dNv/dv*w)
                tmp1[1] = ders [ 0 ](0, k) * ders [ 1 ](1, l) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, l) * w * wders[2];

                answer(cnt, 0) = ( +jacobian(1, 1) * tmp1[0] - jacobian(0, 1) * tmp1[1] ) / product;
                answer(cnt, 1) = ( -jacobian(1, 0) * tmp1[0] + jacobian(0, 0) * tmp1[1] ) / product;
                cnt++;
            }

            ind += numberOfControlPoints [ 0 ];
        }
    } else if ( nsd == 3 ) {
        FloatArray tmp1(nsd + 1), tmp2(nsd + 1);    // allow for weight
        FloatArray temp1(nsd + 1), temp2(nsd + 1), temp3(nsd + 1); // allow for weight

        // calculate values and derivatives of nonrational Bspline solid with weights at first (Aders, wders)
        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int tind = span[2] - degree [ 2 ];
        int ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int m = 0; m <= degree [ 2 ]; m++ ) {
            temp1.zero();
            temp2.zero();
            temp3.zero();
            int indx = ind;
            for ( int l = 0; l <= degree [ 1 ]; l++ ) {
                tmp1.zero();
                tmp2.zero();
                for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                    const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
                    double w = this->weights.at(ind + k);

                    tmp1[0] += ders [ 0 ](0, k) * vertexCoords[0] * w;              // sum(Nu*x*w)
                    tmp1[1] += ders [ 0 ](0, k) * vertexCoords[1] * w;              // sum(Nu*y*w)
                    tmp1[2] += ders [ 0 ](0, k) * vertexCoords[2] * w;              // sum(Nu*z*w)
                    tmp1[3] += ders [ 0 ](0, k) * w;                                       // sum(Nu*w)

                    tmp2[0] += ders [ 0 ](1, k) * vertexCoords[0] * w;              // sum(dNu/du*x*w)
                    tmp2[1] += ders [ 0 ](1, k) * vertexCoords[1] * w;              // sum(dNu/du*y*w)
                    tmp2[2] += ders [ 0 ](1, k) * vertexCoords[2] * w;              // sum(dNu/du*z*w)
                    tmp2[3] += ders [ 0 ](1, k) * w;                                       // sum(dNu/du*w)
                }

                ind += numberOfControlPoints [ 0 ];

                temp1[0] += ders [ 1 ](0, l) * tmp1[0];            // sum(Nv*sum(Nu*x*w))
                temp1[1] += ders [ 1 ](0, l) * tmp1[1];            // sum(Nv*sum(Nu*y*w))
                temp1[2] += ders [ 1 ](0, l) * tmp1[2];            // sum(Nv*sum(Nu*z*w))
                temp1[3] += ders [ 1 ](0, l) * tmp1[3];            // sum(Nv*sum(Nu*w))

                temp2[0] += ders [ 1 ](0, l) * tmp2[0];            // sum(Nv*sum(dNu/du*x*w))
                temp2[1] += ders [ 1 ](0, l) * tmp2[1];            // sum(Nv*sum(dNu/du*y*w))
                temp2[2] += ders [ 1 ](0, l) * tmp2[2];            // sum(Nv*sum(dNu/du*z*w))
                temp2[3] += ders [ 1 ](0, l) * tmp2[3];            // sum(Nv*sum(dNu/du*w))

                temp3[0] += ders [ 1 ](1, l) * tmp1[0];            // sum(dNv/dv*sum(Nu*x*w))
                temp3[1] += ders [ 1 ](1, l) * tmp1[1];            // sum(dNv/dv*sum(Nu*y*w))
                temp3[2] += ders [ 1 ](1, l) * tmp1[2];            // sum(dNv/dv*sum(Nu*z*w))
                temp3[3] += ders [ 1 ](1, l) * tmp1[3];            // sum(dNv/dv*sum(Nu*w))
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];

            Aders [ 0 ][0] += ders [ 2 ](0, m) * temp1[0];     // x=sum(Nt*sum(Nv*sum(Nu*x*w)))
            Aders [ 1 ][0] += ders [ 2 ](0, m) * temp1[1];     // y=sum(Nt*sum(Nv*sum(Nu*y*w)))
            Aders [ 2 ][0] += ders [ 2 ](0, m) * temp1[2];     // z=sum(Nt*sum(Nv*sum(Nu*z*w)))
            wders[0]    += ders [ 2 ](0, m) * temp1[3];        // w=sum(Nt*sum(Nv*sum(Nu*w)))

            Aders [ 0 ][1] += ders [ 2 ](0, m) * temp2[0];     // dx/du=sum(Nt*sum(Nv*sum(dNu/du*x*w)))
            Aders [ 1 ][1] += ders [ 2 ](0, m) * temp2[1];     // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y*w)))
            Aders [ 2 ][1] += ders [ 2 ](0, m) * temp2[2];     // dz/du=sum(Nt*sum(Nv*sum(dNu/du*z*w)))
            wders[1]    += ders [ 2 ](0, m) * temp2[3];        // dw/du=sum(Nt*sum(Nv*sum(dNu/du*w)))

            Aders [ 0 ][2] += ders [ 2 ](0, m) * temp3[0];     // dx/dv=sum(Nt*sum(dNv/dv*sum(Nu*x*w)))
            Aders [ 1 ][2] += ders [ 2 ](0, m) * temp3[1];     // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y*w)))
            Aders [ 2 ][2] += ders [ 2 ](0, m) * temp3[2];     // dz/dv=sum(Nt*sum(dNv/dv*sum(Nu*z*w)))
            wders[2]    += ders [ 2 ](0, m) * temp3[3];        // dw/dv=sum(Nt*sum(dNv/dv*sum(Nu*w)))

            Aders [ 0 ][3] += ders [ 2 ](1, m) * temp1[0];     // dx/dt=sum(dNt/dt*sum(Nv*sum(Nu*x*w)))
            Aders [ 1 ][3] += ders [ 2 ](1, m) * temp1[1];     // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y*w)))
            Aders [ 2 ][3] += ders [ 2 ](1, m) * temp1[2];     // dz/dt=sum(dNt/dt*sum(Nv*sum(Nu*z*w)))
            wders[3]    += ders [ 2 ](1, m) * temp1[3];        // dw/dt=sum(dNt/dt*sum(Nv*sum(Nu*w)))
        }

        double weight = wders[0];

        // calculation of jacobian matrix
        tmp1[0] = Aders [ 0 ][0] / weight;
        tmp1[1] = Aders [ 1 ][0] / weight;
        tmp1[2] = Aders [ 2 ][0] / weight;
        jacobian(0, 0) = ( Aders [ 0 ][1] - wders[1] * tmp1[0] ) / weight; // dx/du
        jacobian(0, 1) = ( Aders [ 1 ][1] - wders[1] * tmp1[1] ) / weight; // dy/du
        jacobian(0, 2) = ( Aders [ 2 ][1] - wders[1] * tmp1[2] ) / weight; // dz/du
        jacobian(1, 0) = ( Aders [ 0 ][2] - wders[2] * tmp1[0] ) / weight; // dx/dv
        jacobian(1, 1) = ( Aders [ 1 ][2] - wders[2] * tmp1[1] ) / weight; // dy/dv
        jacobian(1, 2) = ( Aders [ 2 ][2] - wders[2] * tmp1[2] ) / weight; // dz/dv
        jacobian(2, 0) = ( Aders [ 0 ][3] - wders[3] * tmp1[0] ) / weight; // dx/dt
        jacobian(2, 1) = ( Aders [ 1 ][3] - wders[3] * tmp1[1] ) / weight; // dy/dt
        jacobian(2, 2) = ( Aders [ 2 ][3] - wders[3] * tmp1[2] ) / weight; // dz/dt

        Jacob = jacobian.giveDeterminant();

        //calculation of derivatives of NURBS basis functions with respect to local parameters is not covered by NURBS book
        double product = Jacob * weight * weight;
        int cnt = 0;
        ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int m = 0; m <= degree [ 2 ]; m++ ) {
            int indx = ind;
            for ( int l = 0; l <= degree [ 1 ]; l++ ) {
                for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                    double w = this->weights.at(ind + k);
                    // dNu/du*Nv*Nt*w*sum(Nu*Nv*Nt*w) - Nu*Nv*Nt*w*sum(dNu/du*Nv*Nt*w)
                    tmp1[0] = ders [ 0 ](1, k) * ders [ 1 ](0, l) * ders [ 2 ](0, m) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, l) * ders [ 2 ](0, m) * w * wders[1];
                    // Nu*dNv/dv*Nt*w*sum(Nu*Nv*Nt*w) - Nu*Nv*Nt*w*sum(Nu*dNv/dv*Nt*w)
                    tmp1[1] = ders [ 0 ](0, k) * ders [ 1 ](1, l) * ders [ 2 ](0, m) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, l) * ders [ 2 ](0, m) * w * wders[2];
                    // Nu*Nv*dNt/dt*w*sum(Nu*Nv*Nt*w) - Nu*Nv*Nt*w*sum(Nu*Nv*dNt/dt*w)
                    tmp1[2] = ders [ 0 ](0, k) * ders [ 1 ](0, l) * ders [ 2 ](1, m) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, l) * ders [ 2 ](0, m) * w * wders[3];

                    answer(cnt, 0) = ( ( jacobian(1, 1) * jacobian(2, 2) - jacobian(1, 2) * jacobian(2, 1) ) * tmp1[0] +
                                      ( jacobian(0, 2) * jacobian(2, 1) - jacobian(0, 1) * jacobian(2, 2) ) * tmp1[1] +
                                      ( jacobian(0, 1) * jacobian(1, 2) - jacobian(0, 2) * jacobian(1, 1) ) * tmp1[2] ) / product;                                                      // dN/dx
                    answer(cnt, 1) = ( ( jacobian(1, 2) * jacobian(2, 0) - jacobian(1, 0) * jacobian(2, 2) ) * tmp1[0] +
                                      ( jacobian(0, 0) * jacobian(2, 2) - jacobian(0, 2) * jacobian(2, 0) ) * tmp1[1] +
                                      ( jacobian(0, 2) * jacobian(1, 0) - jacobian(0, 0) * jacobian(1, 2) ) * tmp1[2] ) / product;                                                      // dN/dy
                    answer(cnt, 2) = ( ( jacobian(1, 0) * jacobian(2, 1) - jacobian(1, 1) * jacobian(2, 0) ) * tmp1[0] +
                                      ( jacobian(0, 1) * jacobian(2, 0) - jacobian(0, 0) * jacobian(2, 1) ) * tmp1[1] +
                                      ( jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0) ) * tmp1[2] ) / product;                                                      // dN/dz
                    cnt++;
                }

                ind += numberOfControlPoints [ 0 ];
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];
        }
    } else {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }


    return Jacob;
}


void NURBSInterpolation :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    /* Based on SurfacePoint A4.3 implementation*/
    const FEIIGAElementGeometryWrapper &gw = static_cast< const FEIIGAElementGeometryWrapper& >(cellgeo);
    IntArray span(nsd);
    double weight = 0.0;
    std :: vector< FloatArray > N(nsd);

    if ( gw.knotSpan ) {
        span = * gw.knotSpan;
    } else {
        for ( int i = 0; i < nsd; i++ ) {
            span[i] = this->findSpan(numberOfControlPoints [ i ] - 1, degree [ i ], lcoords[i], knotVector [ i ]);
        }
    }

    for ( int i = 0; i < nsd; i++ ) {
        this->basisFuns(N [ i ], span[i], lcoords[i], degree [ i ], knotVector [ i ]);
    }

    answer.resize(nsd);
    answer.zero();

    if ( nsd == 1 ) {
        int uind = span[0] - degree [ 0 ];
        int ind = uind + 1;
        for ( int k = 0; k <= degree [ 0 ]; k++ ) {
            const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
            double w = this->weights.at(ind + k);
            answer[0] += N [ 0 ][k] * vertexCoords[0] * w;       // xw=sum(Nu*x*w)
            weight    += N [ 0 ][k] * w;                                // w=sum(Nu*w)
        }
    } else if ( nsd == 2 ) {
        FloatArray tmp(nsd + 1);       // allow for weight

        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int l = 0; l <= degree [ 1 ]; l++ ) {
            tmp.zero();
            for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
                double w = this->weights.at(ind + k);
                tmp[0] += N [ 0 ][k] * vertexCoords[0] * w; // sum(Nu*x*w)
                tmp[1] += N [ 0 ][k] * vertexCoords[1] * w; // sum(Nu*y*w)
                tmp[2] += N [ 0 ][k] * w;            // sum(Nu*w)
            }

            ind += numberOfControlPoints [ 0 ];
            answer[0] += N [ 1 ][l] * tmp[0]; // xw=sum(Nv*Nu*x*w)
            answer[1] += N [ 1 ][l] * tmp[1]; // yw=sum(Nv*Nu*y*w)
            weight    += N [ 1 ][l] * tmp[2]; // w=sum(Nv*Nu*w)
        }
    } else if ( nsd == 3 ) {
        FloatArray tmp(nsd + 1), temp(nsd + 1);    // allow for weight

        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int tind = span[2] - degree [ 2 ];
        int ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int m = 0; m <= degree [ 2 ]; m++ ) {
            temp.zero();
            int indx = ind;
            for ( int l = 0; l <= degree [ 1 ]; l++ ) {
                tmp.zero();
                for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                    const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
                    double w = this->weights.at(ind + k);
                    tmp[0] += N [ 0 ][k] * vertexCoords[0] * w;               // sum(Nu*x*w)
                    tmp[1] += N [ 0 ][k] * vertexCoords[1] * w;               // sum(Nu*y*w)
                    tmp[2] += N [ 0 ][k] * vertexCoords[2] * w;               // sum(Nu*z*w)
                    tmp[3] += N [ 0 ][k] * w;                                        // sum(Nu*w)
                }

                ind += numberOfControlPoints [ 0 ];

                temp[0] += N [ 1 ][l] * tmp[0];             // sum(Nv*Nu*x*w)
                temp[1] += N [ 1 ][l] * tmp[1];             // sum(Nv*Nu*y*w)
                temp[2] += N [ 1 ][l] * tmp[2];             // sum(Nv*Nu*z*w)
                temp[3] += N [ 1 ][l] * tmp[3];             // sum(Nv*Nu*w)
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];

            answer[0] += N [ 2 ][m] * temp[0]; // xw=sum(Nv*Nu*Nt*x*w)
            answer[1] += N [ 2 ][m] * temp[1]; // yw=sum(Nv*Nu*Nt*y*w)
            answer[2] += N [ 2 ][m] * temp[2]; // zw=sum(Nv*Nu*Nt*z*w)
            weight    += N [ 2 ][m] * temp[3]; // w=sum(Nv*Nu*Nt*w)
        }
    } else {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }

    answer.times(1.0 / weight);
}


void NURBSInterpolation :: giveJacobianMatrixAt(FloatMatrix &jacobian, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    //
    // Based on Algorithm A4.4 (p. 137) for d=1
    //
    const FEIIGAElementGeometryWrapper &gw = static_cast< const FEIIGAElementGeometryWrapper& >(cellgeo);
    IntArray span(nsd);
    std :: vector< FloatMatrix > ders(nsd);
    jacobian.resize(nsd, nsd);

    if ( gw.knotSpan ) {
        span = * gw.knotSpan;
    } else {
        for ( int i = 0; i < nsd; i++ ) {
            span[i] = this->findSpan(numberOfControlPoints [ i ] - 1, degree [ i ], lcoords[i], knotVector [ i ]);
        }
    }

    for ( int i = 0; i < nsd; i++ ) {
        this->dersBasisFuns(1, lcoords[i], span[i], degree [ i ], knotVector [ i ], ders [ i ]);
    }

    std :: vector< FloatArray > Aders(nsd);
    FloatArray wders;          // 0th and 1st derivatives in w direction on BSpline

    for ( int i = 0; i < nsd; i++ ) {
        Aders [ i ].resize(nsd + 1);
        Aders [ i ].zero();
    }

    wders.resize(nsd + 1);
    wders.zero();

    if ( nsd == 1 ) {
        // calculate values and derivatives of nonrational Bspline curve with weights at first (Aders, wders)
        int uind = span[0] - degree [ 0 ];
        int ind = uind + 1;
        for ( int k = 0; k <= degree [ 0 ]; k++ ) {
            const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
            double w = this->weights.at(ind + k);

            Aders [ 0 ][0] += ders [ 0 ](0, k) * vertexCoords[0] * w;   // xw=sum(Nu*x*w)
            wders[0]    += ders [ 0 ](0, k) * w;                               // w=sum(Nu*w)

            Aders [ 0 ][1] += ders [ 0 ](1, k) * vertexCoords[0] * w;   // dxw/du=sum(dNu/du*x*w)
            wders[1]    += ders [ 0 ](1, k) * w;                               // dw/du=sum(dNu/du*w)
        }

        double weight = wders[0];

        // calculation of jacobian matrix according to Eq 4.7
        jacobian(0, 0) = ( Aders [ 0 ][1] - wders[1] * Aders [ 0 ][0] / weight ) / weight; // dx/du
    } else if ( nsd == 2 ) {
        FloatArray tmp1(nsd + 1), tmp2(nsd + 1);    // allow for weight

        // calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int l = 0; l <= degree [ 1 ]; l++ ) {
            tmp1.zero();
            tmp2.zero();
            for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
                double w = this->weights.at(ind + k);

                tmp1[0] += ders [ 0 ](0, k) * vertexCoords[0] * w; // sum(Nu*x*w)
                tmp1[1] += ders [ 0 ](0, k) * vertexCoords[1] * w; // sum(Nu*y*w)
                tmp1[2] += ders [ 0 ](0, k) * w;           // sum(Nu*w)

                tmp2[0] += ders [ 0 ](1, k) * vertexCoords[0] * w; // sum(dNu/du*x*w)
                tmp2[1] += ders [ 0 ](1, k) * vertexCoords[1] * w; // sum(dNu/du*y*w)
                tmp2[2] += ders [ 0 ](1, k) * w;           // sum(dNu/du*w)
            }

            ind += numberOfControlPoints [ 0 ];

            Aders [ 0 ][0] += ders [ 1 ](0, l) * tmp1[0]; // xw=sum(Nv*sum(Nu*x*w)
            Aders [ 1 ][0] += ders [ 1 ](0, l) * tmp1[1]; // yw=sum(Nv*sum(Nu*y*w)
            wders[0]    += ders [ 1 ](0, l) * tmp1[2]; // w=sum(Nv*sum(Nu*w)

            Aders [ 0 ][1] += ders [ 1 ](0, l) * tmp2[0]; // dxw/du=sum(Nv*sum(dNu/du*x*w)
            Aders [ 1 ][1] += ders [ 1 ](0, l) * tmp2[1]; // dyw/du=sum(Nv*sum(dNu/du*y*w)
            wders[1]    += ders [ 1 ](0, l) * tmp2[2];        // dw/du=sum(Nv*sum(dNu/du*w)

            Aders [ 0 ][2] += ders [ 1 ](1, l) * tmp1[0]; // dxw/dv=sum(dNv/dv*sum(Nu*x*w)
            Aders [ 1 ][2] += ders [ 1 ](1, l) * tmp1[1]; // dyw/dv=sum(dNv/dv*sum(Nu*y*w)
            wders[2]    += ders [ 1 ](1, l) * tmp1[2]; // dw/dv=sum(dNv/dv*sum(Nu*w)
        }

        double weight = wders[0];

        // calculation of jacobian matrix according to Eq 4.19
        tmp1[0] = Aders [ 0 ][0] / weight;
        tmp1[1] = Aders [ 1 ][0] / weight;
        jacobian(0, 0) = ( Aders [ 0 ][1] - wders[1] * tmp1[0] ) / weight; // dx/du
        jacobian(0, 1) = ( Aders [ 1 ][1] - wders[1] * tmp1[1] ) / weight; // dy/du
        jacobian(1, 0) = ( Aders [ 0 ][2] - wders[2] * tmp1[0] ) / weight; // dx/dv
        jacobian(1, 1) = ( Aders [ 1 ][2] - wders[2] * tmp1[1] ) / weight; // dy/dv
    } else if ( nsd == 3 ) {
        FloatArray tmp1(nsd + 1), tmp2(nsd + 1);    // allow for weight
        FloatArray temp1(nsd + 1), temp2(nsd + 1), temp3(nsd + 1); // allow for weight

        // calculate values and derivatives of nonrational Bspline solid with weights at first (Aders, wders)
        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int tind = span[2] - degree [ 2 ];
        int ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int m = 0; m <= degree [ 2 ]; m++ ) {
            temp1.zero();
            temp2.zero();
            temp3.zero();
            int indx = ind;
            for ( int l = 0; l <= degree [ 1 ]; l++ ) {
                tmp1.zero();
                tmp2.zero();
                for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                    const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
                    double w = this->weights.at(ind + k);

                    tmp1[0] += ders [ 0 ](0, k) * vertexCoords[0] * w;              // sum(Nu*x*w)
                    tmp1[1] += ders [ 0 ](0, k) * vertexCoords[1] * w;              // sum(Nu*y*w)
                    tmp1[2] += ders [ 0 ](0, k) * vertexCoords[2] * w;              // sum(Nu*z*w)
                    tmp1[3] += ders [ 0 ](0, k) * w;                                       // sum(Nu*w)

                    tmp2[0] += ders [ 0 ](1, k) * vertexCoords[0] * w;              // sum(dNu/du*x*w)
                    tmp2[1] += ders [ 0 ](1, k) * vertexCoords[1] * w;              // sum(dNu/du*y*w)
                    tmp2[2] += ders [ 0 ](1, k) * vertexCoords[2] * w;              // sum(dNu/du*z*w)
                    tmp2[3] += ders [ 0 ](1, k) * w;                                       // sum(dNu/du*w)
                }

                ind += numberOfControlPoints [ 0 ];

                temp1[0] += ders [ 1 ](0, l) * tmp1[0];            // sum(Nv*sum(Nu*x*w))
                temp1[1] += ders [ 1 ](0, l) * tmp1[1];            // sum(Nv*sum(Nu*y*w))
                temp1[2] += ders [ 1 ](0, l) * tmp1[2];            // sum(Nv*sum(Nu*z*w))
                temp1[3] += ders [ 1 ](0, l) * tmp1[3];            // sum(Nv*sum(Nu*w))

                temp2[0] += ders [ 1 ](0, l) * tmp2[0];            // sum(Nv*sum(dNu/du*x*w))
                temp2[1] += ders [ 1 ](0, l) * tmp2[1];            // sum(Nv*sum(dNu/du*y*w))
                temp2[2] += ders [ 1 ](0, l) * tmp2[2];            // sum(Nv*sum(dNu/du*z*w))
                temp2[3] += ders [ 1 ](0, l) * tmp2[3];            // sum(Nv*sum(dNu/du*w))

                temp3[0] += ders [ 1 ](1, l) * tmp1[0];            // sum(dNv/dv*sum(Nu*x*w))
                temp3[1] += ders [ 1 ](1, l) * tmp1[1];            // sum(dNv/dv*sum(Nu*y*w))
                temp3[2] += ders [ 1 ](1, l) * tmp1[2];            // sum(dNv/dv*sum(Nu*z*w))
                temp3[3] += ders [ 1 ](1, l) * tmp1[3];            // sum(dNv/dv*sum(Nu*w))
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];

            Aders [ 0 ][0] += ders [ 2 ](0, m) * temp1[0];     // x=sum(Nt*sum(Nv*sum(Nu*x*w)))
            Aders [ 1 ][0] += ders [ 2 ](0, m) * temp1[1];     // y=sum(Nt*sum(Nv*sum(Nu*y*w)))
            Aders [ 2 ][0] += ders [ 2 ](0, m) * temp1[2];     // z=sum(Nt*sum(Nv*sum(Nu*z*w)))
            wders[0]    += ders [ 2 ](0, m) * temp1[3];        // w=sum(Nt*sum(Nv*sum(Nu*w)))

            Aders [ 0 ][1] += ders [ 2 ](0, m) * temp2[0];     // dx/du=sum(Nt*sum(Nv*sum(dNu/du*x*w)))
            Aders [ 1 ][1] += ders [ 2 ](0, m) * temp2[1];     // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y*w)))
            Aders [ 2 ][1] += ders [ 2 ](0, m) * temp2[2];     // dz/du=sum(Nt*sum(Nv*sum(dNu/du*z*w)))
            wders[1]    += ders [ 2 ](0, m) * temp2[3];        // dw/du=sum(Nt*sum(Nv*sum(dNu/du*w)))

            Aders [ 0 ][2] += ders [ 2 ](0, m) * temp3[0];     // dx/dv=sum(Nt*sum(dNv/dv*sum(Nu*x*w)))
            Aders [ 1 ][2] += ders [ 2 ](0, m) * temp3[1];     // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y*w)))
            Aders [ 2 ][2] += ders [ 2 ](0, m) * temp3[2];     // dz/dv=sum(Nt*sum(dNv/dv*sum(Nu*z*w)))
            wders[2]    += ders [ 2 ](0, m) * temp3[3];        // dw/dv=sum(Nt*sum(dNv/dv*sum(Nu*w)))

            Aders [ 0 ][3] += ders [ 2 ](1, m) * temp1[0];     // dx/dt=sum(dNt/dt*sum(Nv*sum(Nu*x*w)))
            Aders [ 1 ][3] += ders [ 2 ](1, m) * temp1[1];     // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y*w)))
            Aders [ 2 ][3] += ders [ 2 ](1, m) * temp1[2];     // dz/dt=sum(dNt/dt*sum(Nv*sum(Nu*z*w)))
            wders[3]    += ders [ 2 ](1, m) * temp1[3];        // dw/dt=sum(dNt/dt*sum(Nv*sum(Nu*w)))
        }

        double weight = wders[0];

        // calculation of jacobian matrix
        tmp1[0] = Aders [ 0 ][0] / weight;
        tmp1[1] = Aders [ 1 ][0] / weight;
        tmp1[2] = Aders [ 2 ][0] / weight;
        jacobian(0, 0) = ( Aders [ 0 ][1] - wders[1] * tmp1[0] ) / weight; // dx/du
        jacobian(0, 1) = ( Aders [ 1 ][1] - wders[1] * tmp1[1] ) / weight; // dy/du
        jacobian(0, 2) = ( Aders [ 2 ][1] - wders[1] * tmp1[2] ) / weight; // dz/du
        jacobian(1, 0) = ( Aders [ 0 ][2] - wders[2] * tmp1[0] ) / weight; // dx/dv
        jacobian(1, 1) = ( Aders [ 1 ][2] - wders[2] * tmp1[1] ) / weight; // dy/dv
        jacobian(1, 2) = ( Aders [ 2 ][2] - wders[2] * tmp1[2] ) / weight; // dz/dv
        jacobian(2, 0) = ( Aders [ 0 ][3] - wders[3] * tmp1[0] ) / weight; // dx/dt
        jacobian(2, 1) = ( Aders [ 1 ][3] - wders[3] * tmp1[1] ) / weight; // dy/dt
        jacobian(2, 2) = ( Aders [ 2 ][3] - wders[3] * tmp1[2] ) / weight; // dz/dt
    } else {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }

}


///////////////////////
// NURBS3dSurfaceInterpolation

void NURBS3dSurfaceInterpolation::initializeFrom( InputRecord &ir )
{
    BSpline3dSurfaceInterpolation::initializeFrom( ir );

    IR_GIVE_FIELD( ir, weights, _IFT_NURBSInterpolation_weights );
    if ( weights.giveSize() != numberOfControlPoints[0] * numberOfControlPoints[1] ) {
        throw ValueInputException( ir, _IFT_NURBSInterpolation_weights, "the number of weights must equal the number of control points" );
    }
    // the rational basis functions are divided by sum(N*w)
    double wmax = 0.;
    for ( double w : weights ) {
        if ( !std::isfinite( w ) || w <= 0. ) {
            throw ValueInputException( ir, _IFT_NURBSInterpolation_weights, "weights must be finite and positive" );
        }
        wmax = std::max( wmax, w );
    }
    // Only the ratios of the weights matter. Divided by the largest one they lie in (0, 1], so the weighted sums of the
    // B-spline derivatives in evalBasis stay in their range (a common factor like 1e307 would overflow them). A division,
    // since 1/wmax may itself overflow; exact when the largest weight is 1. A ratio below the normal range of double
    // would lose its precision or vanish.
    for ( double &w : weights ) {
        w /= wmax;
        if ( !( w >= std::numeric_limits< double >::min() ) ) {
            throw ValueInputException( ir, _IFT_NURBSInterpolation_weights,
                "the ratio of the smallest to the largest weight is below the range of double" );
        }
    }
}


void NURBS3dSurfaceInterpolation::evalBasis( int nder, const FloatArray &lcoords, const FEICellGeometry &cellgeo, FloatMatrix &answer,
    IntArray &cps ) const
{
    BSpline3dSurfaceInterpolation::evalBasis( nder, lcoords, cellgeo, answer, cps );

    // W = sum w N and its derivatives, in the columns of answer (W > 0: positive weights, B-splines >= 0 summing to 1)
    int count = answer.giveNumberOfRows();
    double W[6] = { 0., 0., 0., 0., 0., 0. };
    for ( int a = 0; a < count; a++ ) {
        double w = this->weights.at( cps[a] );
        for ( int c = 0; c < 6; c++ ) {
            W[c] += w * answer( a, c );
        }
    }

    // R = w N / W and its derivatives from R W = w N: columns N, N_,xi, N_,eta, N_,xixi, N_,etaeta, N_,xieta
    for ( int a = 0; a < count; a++ ) {
        double w  = this->weights.at( cps[a] );
        double R  = w * answer( a, 0 ) / W[0];
        double Ru = ( w * answer( a, 1 ) - R * W[1] ) / W[0];
        double Rv = ( w * answer( a, 2 ) - R * W[2] ) / W[0];
        if ( nder >= 2 ) {
            answer( a, 3 ) = ( w * answer( a, 3 ) - 2. * Ru * W[1] - R * W[3] ) / W[0];
            answer( a, 4 ) = ( w * answer( a, 4 ) - 2. * Rv * W[2] - R * W[4] ) / W[0];
            answer( a, 5 ) = ( w * answer( a, 5 ) - Ru * W[2] - Rv * W[1] - R * W[5] ) / W[0];
        }
        answer( a, 0 ) = R;
        if ( nder >= 1 ) {
            answer( a, 1 ) = Ru;
            answer( a, 2 ) = Rv;
        }
    }
}


void NURBS3dSurfaceInterpolation::evalN( FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const
{
    FloatMatrix basis;
    IntArray cps;
    this->evalBasis( 0, lcoords, cellgeo, basis, cps );
    answer.resize( basis.giveNumberOfRows() );
    for ( int a = 0; a < basis.giveNumberOfRows(); a++ ) {
        answer[a] = basis( a, 0 );
    }
}


void NURBS3dSurfaceInterpolation::local2global( FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const
{
    FloatMatrix basis;
    IntArray cps;
    this->evalBasis( 0, lcoords, cellgeo, basis, cps );
    answer.resize( 3 );
    answer.zero();
    for ( int a = 0; a < basis.giveNumberOfRows(); a++ ) {
        const auto &X = giveControlPointCoordinates( cellgeo, cps[a] );
        for ( int i = 0; i < 3; i++ ) {
            answer[i] += basis( a, 0 ) * X[i];
        }
    }
}


void NURBS3dSurfaceInterpolation::evalSecondDerivatives( FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const
{
    IntArray cps;
    this->evalBasis( 2, lcoords, cellgeo, answer, cps );
}

} // end namespace oofem

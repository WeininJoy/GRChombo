/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXSTATICVORTEX_HPP_
#define COMPLEXSTATICVORTEX_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the initial conditions for a static vortex
// starting with cylindrical coords
// ds^2 = -dt^2 + exp(-4*pi*G*C^2*r^2n)*(dr^2 + r^2 d\psi^2) + dz^2
// \phi = C*r^n * exp(i*n*\psi)
// Use coordinate transform to transfer g, K, shift in cylindrical coords
// to h, A, shift in Cartesian coords 

class ComplexStaticVortex
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double Amp = 0.1;                      //!<< The amplitude (C) of the static vortex
        int n = 1;                             //!<< The winding number (n) of the static vortex
        std::array<double, CH_SPACEDIM> center; //!< The center of the vortex
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    const params_t m_params;
    const double m_dx;

    ComplexStaticVortex(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx) {}

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // set up vars for the metric and extrinsic curvature, shift and lapse in
        // cylindrical coords
        Tensor<2, data_t> cylindrical_g;
        Tensor<2, data_t> cylindrical_K;
        Tensor<1, data_t> cylindrical_shift;
        data_t vortex_lapse;
        
        // get position and set vars

        Vars<data_t> metric_vars;
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        // Compute the components in cylindrical coords
        compute_metric_vortex(cylindrical_g, cylindrical_K, cylindrical_shift, vortex_lapse, coords);
        
        // work out where we are on the grid
        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;

        using namespace CoordinateTransformations;
        // Convert spherical components to cartesian components using coordinate
        // transforms
        vars.h = spherical_to_cartesian_LL(spherical_g, x, y, z);
        vars.A = spherical_to_cartesian_LL(spherical_K, x, y, z);
        vars.shift = spherical_to_cartesian_U(spherical_shift, x, y, z);

        using namespace TensorAlgebra;
        // Convert to BSSN vars
        data_t deth = compute_determinant(vars.h);
        auto h_UU = compute_inverse_sym(vars.h);
        vars.chi = pow(deth, -1. / 3.);

        // transform extrinsic curvature into A and TrK - note h is still non
        // conformal version which is what we need here
        vars.K = compute_trace(vars.A, h_UU);
        make_trace_free(vars.A, vars.h, h_UU);

        // Make conformal
        FOR(i, j)
        {
            vars.h[i][j] *= vars.chi;
            vars.A[i][j] *= vars.chi;
        }

        // use a pre collapsed lapse, could also use analytic one
        // vars.lapse = kerr_lapse;
        vars.lapse = pow(vars.chi, 0.5);

        // Populate the variables on the grid
        // NB We stil need to set Gamma^i which is NON ZERO
        // but we do this via a separate class/compute function
        // as we need the gradients of the metric which are not yet available
        current_cell.store_vars(vars);
    }

};

#endif /* COMPLEXSTATICVORTEX_HPP_ */

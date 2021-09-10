/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_chi, // the conformal factor of the metric, not evolved
    c_rho, // the energy density of the SF
    //c_cir_x,  //circulation (x direction)
    //c_cir_y,  //circulation (y direction)
    //c_cir_z,  //circulation (z direction)

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi", "rho"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */

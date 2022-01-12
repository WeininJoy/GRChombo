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
    c_cir,  // circulation
    c_gamma_00, // gamma[0][0]
    c_K_00, // K_tensor[0][0]
    c_K, // K
    c_shift_0, // shift[0]
    c_d1_gamma_000, // d1_gamma[0][0][0]

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi", "rho", "cir", "gamma_00", "K_00" ,"K", "shift_0", "d1_gamma_000"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */

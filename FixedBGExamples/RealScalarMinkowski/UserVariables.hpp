/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "DiagnosticVariables.hpp"
#include <array>
#include <string>

// assign an enum to each variable
enum
{
    c_phi, // matter field added
    c_Pi,  //(minus) conjugate momentum
    c_cir_x,  //circulation (x direction)
    c_cir_y,  //circulation (y direction)
    c_cir_z,  //circulation (z direction)

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {"phi_Re",
                                                                 "Pi_Re",
                                                                 "cir_x",
                                                                 "cir_y",
                                                                 "cir_z"};
}

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */

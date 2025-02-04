/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "DiagnosticVariables.hpp"
#include "CCZ4UserVariables.hpp"

// assign an enum to each variable
enum
{
    c_phi_Re = NUM_CCZ4_VARS, // matter field added
    c_Pi_Re,  //(minus) conjugate momentum
    c_phi_Im, // matter field added
    c_Pi_Im,  //(minus) conjugate momentum

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_names = {"phi_Re","Pi_Re", "phi_Im", "Pi_Im"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */

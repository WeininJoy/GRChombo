/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCALARPOTENTIAL_HPP_
#define SCALARPOTENTIAL_HPP_

#include "simd.hpp"

class ScalarPotential
{
  public:
    struct params_t
    {
        double scalar_mass;
    };

  private:
    const params_t m_params;

  public:
    //! The constructor
    ScalarPotential(const params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi_Re, data_t &dVdphi_Im,
                           const vars_t<data_t> &vars) const
    {
        const double m = m_params.scalar_mass;
        // The potential value at phi
        V_of_phi = 0.5 * m * m * pow(vars.phi_Re,2) * pow(vars.phi_Im,2);

        // The potential gradient at phi
        dVdphi_Re = m * m * vars.phi_Re;
        dVdphi_Im = m * m * vars.phi_Im;
    }
};

#endif /* SCALARPOTENTIAL_HPP_ */

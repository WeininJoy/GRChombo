/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "BoostedIsotropicBHFixedBG2.hpp"
#include "ExcisionDiagnostics.hpp"
#include "ExcisionEvolution.hpp"
#include "FixedBGDensity.hpp"
#include "FixedBGEvolution.hpp"
#include "FixedBGScalarField.hpp"
#include "InitialConditions.hpp"
#include "ScalarPotential.hpp"
#include "RealScalarCirculation.hpp"
#include "CirculationExtraction.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero to avoid undefined values,
    // then initial conditions for fields
    SetValue set_zero(0.0);
    BoostedIsotropicBHFixedBG boosted_bh(m_p.bg_params,
                                         m_dx); // just calculates chi
    InitialConditions<BoostedIsotropicBHFixedBG> set_phi(
        m_p.field_amplitude, m_p.potential_params.scalar_mass, m_p.center,
        boosted_bh, m_dx);
    auto compute_pack = make_compute_pack(set_zero, boosted_bh);
    BoxLoops::loop(compute_pack, m_state_diagnostics, m_state_diagnostics,
                   SKIP_GHOST_CELLS);
    BoxLoops::loop(set_phi, m_state_new, m_state_new, FILL_GHOST_CELLS);

    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<ScalarFieldWithPotential, BoostedIsotropicBHFixedBG>(
            m_dx, m_p.center, boosted_bh),
        m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd());
}

// Things to do before outputting a plot file
void ScalarFieldLevel::prePlotLevel() {

    ScalarPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoostedIsotropicBHFixedBG boosted_bh(m_p.bg_params, m_dx);
    RealScalarCirculation<ScalarFieldWithPotential, BoostedIsotropicBHFixedBG>
        circulation(scalar_field, boosted_bh, potential, m_dx, m_p.center);  //Calculate circulation

    BoxLoops::loop(circulation, m_state_new, m_state_diagnostics, SKIP_GHOST_CELLS); //pass m_state_diagnostics as the output to have the circ vars as diagnostic vars 
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    // We don't want undefined values floating around in the constraints so
    // zero these
    ScalarPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoostedIsotropicBHFixedBG boosted_bh(m_p.bg_params, m_dx);
    FixedBGEvolution<ScalarFieldWithPotential, BoostedIsotropicBHFixedBG>
        my_matter(scalar_field, boosted_bh, m_p.sigma, m_dx, m_p.center);

    BoxLoops::loop(my_matter, a_soln, a_rhs, SKIP_GHOST_CELLS);

    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<ScalarFieldWithPotential, BoostedIsotropicBHFixedBG>(
            m_dx, m_p.center, boosted_bh),
        a_soln, a_rhs, SKIP_GHOST_CELLS, disable_simd());
}

void ScalarFieldLevel::specificPostTimeStep()
{
    // At any level, but after the coarsest timestep
    int min_level = 0;
    bool calculate_quantities = at_level_timestep_multiple(min_level);

    if (calculate_quantities)
    {
        fillAllGhosts();
        ScalarPotential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        BoostedIsotropicBHFixedBG boosted_bh(m_p.bg_params, m_dx);
        FixedBGDensity<ScalarFieldWithPotential, BoostedIsotropicBHFixedBG>
            density(scalar_field, boosted_bh, m_dx, m_p.center);
        BoxLoops::loop(density, m_state_new, m_state_diagnostics,
                       SKIP_GHOST_CELLS);
        // excise within horizon and outside extraction radius
        BoxLoops::loop(
            ExcisionDiagnostics<ScalarFieldWithPotential,
                                BoostedIsotropicBHFixedBG>(
                m_dx, m_p.center, boosted_bh, m_p.inner_r, m_p.outer_r),
            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
            disable_simd());
    }

    // write out the integral after each coarse timestep on rank 0
    if (m_level == 0)
    {
        bool first_step = (m_time == m_dt);
        // integrate the densities and sources and write to a file
        AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
        double rho_sum = amr_reductions.sum(c_rho);

        SmallDataIO integral_file("VolumeIntegrals", m_dt, m_time,
                                  m_restart_time, SmallDataIO::APPEND,
                                  first_step);
        // remove any duplicate data if this is post restart
        integral_file.remove_duplicate_time_data();
        std::vector<double> data_for_writing = {rho_sum};
        // write data
        if (first_step)
        {
            integral_file.write_header_line({"rho"});
        }
        integral_file.write_time_data_line(data_for_writing);

        // set up an interpolator
        // pass the boundary params so that we can use symmetries if
        // applicable
        AMRInterpolator<Lagrange<4>> interpolator(
            m_gr_amr, m_p.origin, m_p.dx, m_p.boundary_params,
            m_p.verbosity);

        // this should fill all ghosts including the boundary ones according
        // to the conditions set in params.txt
        interpolator.refresh();

        // set up the query and execute it
        int num_points = 1000;
        CirculationExtraction extraction(c_cir, num_points, m_p.L, m_p.center,
                                        m_dt, m_time);
        extraction.execute_query(&interpolator, "outputs_1000");
    }
}

void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.regrid_length,
                                              m_p.center),
                   current_state, tagging_criterion, disable_simd());
}

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
// #include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
// #include "SetValue.hpp"
// #include "SmallDataIO.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4RHS.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "InitialConditions.hpp"
#include "ComplexStaticVortex.hpp"
#include "ScalarPotential.hpp"
// #include "ComplexScalarCirculation.hpp"
// #include "CirculationExtraction.hpp"
#include "ComplexScalarField.hpp"
#include "SetValue.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                        PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
                       m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, EXCLUDE_GHOST_CELLS,
                       disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero to avoid undefined values,
    // then initial conditions for complex scalar field 
    BoxLoops::loop(
        make_compute_pack(SetValue(0.), ComplexStaticVortex(m_p.vortex_params, m_dx),
                          InitialConditions(m_p.initial_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
    
    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

#ifdef CH_USE_HDF5
// Things to do before outputting a plot file
void ScalarFieldLevel::prePlotLevel() {

    fillAllGhosts();
    ScalarPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(
        MatterConstraints<ScalarFieldWithPotential>(
            scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);


    /*
    ScalarPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    ComplexStaticVortex static_vor(m_p.vortex_params, m_dx);
    ComplexScalarCirculation<ScalarFieldWithPotential, ComplexStaticVortex>
        circulation(scalar_field, static_vor, potential, m_dx, m_p.center);  //Calculate circulation

    // pass m_state_diagnostics as the output to have the circ vars as diagnostic vars 
    BoxLoops::loop(circulation, m_state_new, m_state_diagnostics, SKIP_GHOST_CELLS); 
    */
}
#endif

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    // We don't want undefined values floating around in the constraints so
    // zero these

    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),  PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
                        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ComplexScalarField
    ScalarPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

void ScalarFieldLevel::specificPostTimeStep()
{
    /*
    // At any level, but after the coarsest timestep
    int min_level = 0;
    bool calculate_quantities = at_level_timestep_multiple(min_level);

    
    if (calculate_quantities)
    {
        fillAllGhosts();
        ScalarPotential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        ComplexStaticVortex static_vor(m_p.vortex_params, m_dx);
        FixedBGDensity<ScalarFieldWithPotential, ComplexStaticVortex>
            density(scalar_field, static_vor, m_dx, m_p.center);
        BoxLoops::loop(density, m_state_new, m_state_diagnostics,
                       SKIP_GHOST_CELLS);
    }
    */

    /*
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
    */
        /*
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
    */
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void ScalarFieldLevel::preTagCells()
{
    // we don't need any ghosts filled for the fixed grids tagging criterion
    // used here so don't fill any
}

void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.regrid_length,
                                              m_p.center),
                   current_state, tagging_criterion);
}

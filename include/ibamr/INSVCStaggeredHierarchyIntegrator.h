// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_INSVCStaggeredHierarchyIntegrator
#define included_IBAMR_INSVCStaggeredHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/SideDataSynchronization.h"
#include "ibtk/ibtk_enums.h"

#include "CellVariable.h"
#include "EdgeVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyEdgeDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchyNodeDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "NodeVariable.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace IBAMR
{
class BrinkmanPenalizationStrategy;
class ConvectiveOperator;
} // namespace IBAMR
namespace IBTK
{
class PoissonSolver;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchLevel;
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class BasePatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
} // namespace mesh
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSVCStaggeredHierarchyIntegrator provides an abstract interface for time
 * integrator for a staggered-grid, incompressible Navier-Stokes solver on an AMR grid hierarchy
 * with variable coefficients.
 *
 * An optional re-scaling factor c can be specified to minimize the loss of floating point precision for poorly
 * scaling linear systems. The scaling acts on the momentum part of the saddle-point system, yielding
 * \f$ c A + G(c x_p) = c b_u\f$
 * in which the viscous block, the pressure degrees of freedom, and the velocity RHS have been scaled.
 *
 * Scaling \f$ c \f$ is chosen such that
 * \f$ c(\frac{\rho}{dt} - \frac{\mu}{dx^2}) \sim \frac{1}{dx} \f$.
 * The above scaling is chosen from the incompressiblity operator which scales
 * as \f$ (\nabla \cdot) \sim \frac{1}{dx} \f$.
 * Assuming \f$ dt \sim dx \f$ and \f$ 1/dx = N \f$, we have
 * \f$ c \sim \frac{1}{\rho - \mu N} \f$. Here \f$ N \f$ is the number of cells for
 * a unit length of the physical domain.
 *
 * Different levels of patch hierarchy can have different scaling because
 * of the difference in the grid spacing \f$ dx \f$. Therefore, an array of scale
 * factors is read from the input file (corresponding to different levels).
 * If the scale array does not contain values for all the levels in the hierarchy,
 * it is filled by the most finest scaling factor provided by the user (for the missing
 * finer levels).
 */

class INSVCStaggeredHierarchyIntegrator : public INSHierarchyIntegrator
{
public:
    /*!
     * The constructor for class INSVCStaggeredHierarchyIntegrator sets some
     * default values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    INSVCStaggeredHierarchyIntegrator(std::string object_name,
                                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                      bool register_for_restart = true);

    /*!
     * The destructor for class INSVCStaggeredHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~INSVCStaggeredHierarchyIntegrator();

    /*!
     * Get the convective operator being used by this solver class.
     *
     * If the time integrator is configured to solve the time-dependent
     * (creeping) Stokes equations, then the returned pointer will be NULL.
     *
     * If the convective operator has not already been constructed, and if the
     * time integrator is not configured to solve the time-dependent (creeping)
     * Stokes equations, then this function will initialize the default type of
     * convective operator, which may be set in the class input database.
     */
    SAMRAI::tbox::Pointer<ConvectiveOperator> getConvectiveOperator() override;

    /*!
     * Get the subdomain solver for the velocity subsystem.  Such solvers can be
     * useful in constructing block preconditioners.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getVelocitySubdomainSolver() override;

    /*!
     * Get the subdomain solver for the pressure subsystem.  Such solvers can be
     * useful in constructing block preconditioners.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getPressureSubdomainSolver() override;

    /*!
     * Register a solver for the time-dependent incompressible Stokes equations.
     */
    void setStokesSolver(SAMRAI::tbox::Pointer<StaggeredStokesSolver> stokes_solver);

    /*!
     * Get the solver for the time-dependent incompressible Stokes equations
     * used by this solver class.
     */
    SAMRAI::tbox::Pointer<StaggeredStokesSolver> getStokesSolver();

    /*!
     * Indicate that the Stokes solver should be (re-)initialized before the
     * next time step.
     */
    void setStokesSolverNeedsInit();

    /*!
     * Virtual method to initialize the variables, basic communications
     * algorithms, solvers, and other data structures used by a concrete time
     * integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    virtual void
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Virtual method to initialize the AMR patch hierarchy and data defined on the hierarchy at
     * the start of a computation.  If the computation is begun from a restart
     * file, the patch hierarchy and patch data are read from the hierarchy
     * database.  Otherwise, the patch hierarchy and patch data are initialized
     * by the gridding algorithm associated with the integrator object.
     *
     * The implementation of this function assumes that the hierarchy exists
     * upon entry to the function, but that it contains no patch levels.  On
     * return from this function, the state of the integrator object will be
     * such that it is possible to step through time via the advanceHierarchy()
     * function.
     */
    virtual void
    initializePatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                             SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Virtual method to prepare to advance the data from current_time to new_time.
     */
    virtual void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Virtual method to clean up data following call(s) to integrateHierarchy().
     */
    virtual void postprocessIntegrateHierarchy(double current_time,
                                               double new_time,
                                               bool skip_synchronize_new_state_data,
                                               int num_cycles = 1) override;

    /*!
     * Explicitly remove nullspace components from a solution vector.
     */
    void removeNullSpace(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol_vec);

    /*!
     * Register a variable mass density variable with the hierarchy integrator.
     */
    void registerMassDensityVariable(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var);

    /*!
     * Get the mass density variable registered with the hierarchy integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getMassDensityVariable() const;

    /*!
     * Register a variable viscosity variable with the hierarchy integrator.
     */
    void registerViscosityVariable(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var);

    /*!
     * Get the viscosity variable registered with the hierarchy integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getViscosityVariable() const;

    /*!
     * Set the interpolation type used for material properties rho
     */
    void setDensityVCInterpolationType(const IBTK::VCInterpType vc_interp_type);

    /*!
     * Set the interpolation type used for material properties mu
     */
    void setViscosityVCInterpolationType(const IBTK::VCInterpType vc_interp_type);

    /*!
     * \brief Function to reset fluid density or viscosity if they are
     * maintained by this integrator.
     */
    using ResetFluidPropertiesFcnPtr = void (*)(int property_idx,
                                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > property_var,
                                                SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                int cycle_num,
                                                double time,
                                                double current_time,
                                                double new_time,
                                                void* ctx);

    /*!
     * \brief Register function to reset fluid density.
     */
    void registerResetFluidDensityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Register function to reset fluid viscosity.
     */
    void registerResetFluidViscosityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Register BrinkmanPenalizationStrategy objects to add Brinkman penalization term
     * in the momentum equation.
     */
    virtual void
    registerBrinkmanPenalizationStrategy(SAMRAI::tbox::Pointer<IBAMR::BrinkmanPenalizationStrategy> brinkman_force);

    /*!
     * \brief Supply initial conditions for the density field, if maintained by the fluid integrator.
     */
    void registerMassDensityInitialConditions(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> rho_init_fcn);

    /*!
     * \brief Supply initial conditions for the viscosity field, if maintained by the fluid integrator.
     */
    void registerViscosityInitialConditions(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> mu_init_fcn);

    /*
     * \brief Pure virtual method to supply boundary conditions for the density field, if maintained by the fluid
     * integrator.
     */
    virtual void registerMassDensityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* rho_bc_coef) = 0;

    /*
     * \brief Supply boundary conditions for the density field, if maintained by the fluid integrator.
     *
     * \note The boundary conditions set here will be overwritten if viscosity if being advected.
     */
    void registerViscosityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* mu_bc_coef);

    /*
     * \brief Set the transported viscosity variable if it is being maintained by the advection-diffusion integrator.
     *
     * \note The variable set here MUST be registered and maintained by the advection-diffusion integrator.
     *
     * \note If multiple advection diffusion integrators are registered, you can specify which advection diffusion
     * integrator is used to evolve the viscosity.
     */
    void
    setTransportedViscosityVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > mu_adv_diff_var,
                                    unsigned int adv_diff_idx = 0);

    /*!
     * \brief Get the transported viscosity variable that is being manintained by an advection-diffusion integrator
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > getTransportedViscosityVariable() const;

    /*!
     * \brief Get the side-centered density patch data index, which will always be the newest one used in the linear
     * operator i.e. rho_sc in rho_sc*u^{n+1} term.
     *
     * \note These patch data will not be deallocated at the end of the time step, so they can be used for various
     * applications
     */
    inline int getLinearOperatorRhoPatchDataIndex() const
    {
        return d_rho_linear_op_idx;
    }

    /*!
     * \brief Get whether or not density is constant
     */
    inline bool rhoIsConstant() const
    {
        return d_rho_is_const;
    }

    /*!
     * \brief Get whether or not viscosity is constant
     */
    inline bool muIsConstant() const
    {
        return d_mu_is_const;
    }

    /*!
     * \brief Get the cell-centered viscosity patch data index, which will always be the newest one used in the linear
     * operator.
     *
     * \note These patch data will not be deallocated at the end of the time step, so they can be used for various
     * applications
     */
    inline int getLinearOperatorMuPatchDataIndex() const
    {
        return d_mu_linear_op_idx;
    }

    /*!
     * \brief Get the interpolated viscosity patch data index, which will always be the newest one used in the linear
     * operator.
     *
     * \note These patch data will not be deallocated at the end of the time step, so they can be used for various
     * applications
     */
    inline int getInterpolatedLinearOperatorMuPatchDataIndex() const
    {
        return d_mu_interp_linear_op_idx;
    }

    /*!
     * \brief Get the scaling factor used for A, p and u_rhs.
     */
    inline SAMRAI::tbox::Array<double> getScalingFactor() const
    {
        return d_A_scale;
    }

    /*!
     * \brief Get the viscosity boundary conditions
     */
    inline SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getViscosityBoundaryConditions() const
    {
        return d_mu_bc_coef;
    }

    /*!
     * \brief Get the Brinkman penalization objects registered with this class.
     */
    const std::vector<SAMRAI::tbox::Pointer<IBAMR::BrinkmanPenalizationStrategy> >&
    getBrinkmanPenalizationStrategy() const
    {
        return d_brinkman_force;
    } // getBrinkmanPenalizationStrategy

    /*!
     * \brief Get "old" velocity variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getOldVelocityVariable() const
    {
        return d_U_old_var;
    } // getOldVelocityVariable

protected:
    /*!
     * L1 norm of the discrete divergence of the fluid velocity before regridding.
     */
    double d_div_U_norm_1_pre = 0.0;

    /*!
     * L2 norm of the discrete divergence of the fluid velocity before regridding.
     */
    double d_div_U_norm_2_pre = 0.0;

    /*!
     * L-infinity norm of the discrete divergence of the fluid velocity before regridding.
     */
    double d_div_U_norm_oo_pre = 0.0;

    /*!
     * L1 norm of the discrete divergence of the fluid velocity after regridding.
     */
    double d_div_U_norm_1_post = 0.0;

    /*!
     * L2 norm of the discrete divergence of the fluid velocity after regridding.
     */
    double d_div_U_norm_2_post = 0.0;

    /*!
     * L-infinity norm of the discrete divergence of the fluid velocity after regridding.
     */
    double d_div_U_norm_oo_post = 0.0;

    /*!
     * Whether we need to perform a regrid projection when (re-)initializing composite hierarchy data.
     */
    bool d_do_regrid_projection = false;

    /*!
     * Determine the largest stable timestep on an individual patch.
     */
    double getStableTimestep(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const override;

    /*!
     * Prepare the current hierarchy for regridding. Here we calculate the divergence.
     */
    void regridHierarchyBeginSpecialized() override;

    /*!
     * Update the current hierarchy data after regridding. Here we recalculate
     * the divergence and, if it has grown by a factor more than
     * d_regrid_max_div_growth_factor, we then project the velocity field onto
     * a divergence-free set of grid functions.
     */
    void regridHierarchyEndSpecialized() override;

    /*!
     * Perform data initialization after the entire hierarchy has been constructed.
     */
    void initializeCompositeHierarchyDataSpecialized(double init_data_time, bool initial_time) override;

    /*!
     * Virtual method to initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     */
    virtual void
    initializeLevelDataSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                   int level_number,
                                   double init_data_time,
                                   bool can_be_refined,
                                   bool initial_time,
                                   SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
                                   bool allocate_data) override;

    /*!
     * Virtual method to reset cached hierarchy dependent data.
     */
    virtual void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level) override;

    /*!
     * Virtual method to set integer tags to "one" in cells where refinement of the given level
     * should occur according to the magnitude of the fluid vorticity.
     */
    virtual void
    applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                     int level_number,
                                     double error_data_time,
                                     int tag_index,
                                     bool initial_time,
                                     bool uses_richardson_extrapolation_too) override;

    /*!
     * Virtual method to prepare variables for plotting.
     */
    virtual void setupPlotDataSpecialized() override;

    /*!
     * Copy data from a side-centered variable to a face-centered variable.
     */
    void copySideToFace(const int U_fc_idx,
                        const int U_sc_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM, double> > d_hier_fc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> > d_hier_sc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyNodeDataOpsReal<NDIM, double> > d_hier_nc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyEdgeDataOpsReal<NDIM, double> > d_hier_ec_data_ops;

    /*
     * Boundary condition and data synchronization operators.
     */
    SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> d_bc_helper;
    SAMRAI::tbox::Pointer<IBTK::SideDataSynchronization> d_side_synch_op;

    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_rho_bdry_bc_fill_op, d_mu_bdry_bc_fill_op;

    /*!
     * Double precision values are (optional) factors used to rescale the
     * density and viscosity for plotting.
     *
     * Boolean values indicates whether to output various quantities for
     * visualization.
     */
    double d_rho_scale = 1.0, d_mu_scale = 1.0;
    bool d_output_rho = false, d_output_mu = false;

    /*
     * Hierarchy operators and solvers.
     */
    int d_coarsest_reset_ln, d_finest_reset_ln;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_U_scratch_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_U_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_U_adv_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_N_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_P_scratch_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_P_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_sol_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_rhs_vec;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > > d_nul_vecs;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > > d_U_nul_vecs;
    bool d_vectors_need_init, d_explicitly_remove_nullspace = false;

    std::string d_stokes_solver_type = StaggeredStokesSolverManager::UNDEFINED,
                d_stokes_precond_type = StaggeredStokesSolverManager::UNDEFINED;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_stokes_solver_db, d_stokes_precond_db;
    SAMRAI::tbox::Pointer<StaggeredStokesSolver> d_stokes_solver;
    bool d_stokes_solver_needs_init;

    /*!
     * Fluid solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_old_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_U_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_P_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_F_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_F_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Q_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_N_old_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Omega_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Div_U_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Omega_Norm_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_regrid_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_src_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_indicator_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_F_div_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_EE_var;

    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_rho_var, d_mu_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_pressure_D_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_pressure_rhs_D_var;
#if (NDIM == 2)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > d_velocity_D_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > d_velocity_rhs_D_var;
#elif (NDIM == 3)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > d_velocity_D_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > d_velocity_rhs_D_var;
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_velocity_D_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_velocity_C_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_velocity_L_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_velocity_rhs_C_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_N_full_var;

    std::string d_N_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_N_refine_type = "CONSERVATIVE_LINEAR_REFINE";

    std::string d_mu_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_mu_refine_type = "CONSERVATIVE_LINEAR_REFINE";
    std::string d_mu_bdry_extrap_type = "CONSTANT";

    std::string d_rho_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_rho_refine_type = "CONSERVATIVE_LINEAR_REFINE";
    std::string d_rho_bdry_extrap_type = "CONSTANT";

    /*!
     * Interpolated material property variables.
     */
#if (NDIM == 2)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > d_mu_interp_var;
#elif (NDIM == 3)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > d_mu_interp_var;
#endif

    /*!
     * Functions resetting rho and mu if they are maintained by this integrator.
     */
    std::vector<ResetFluidPropertiesFcnPtr> d_reset_rho_fcns, d_reset_mu_fcns;
    std::vector<void*> d_reset_rho_fcns_ctx, d_reset_mu_fcns_ctx;

    /*!
     * Brinkman force strategy objects registered with this integrator.
     */
    std::vector<SAMRAI::tbox::Pointer<IBAMR::BrinkmanPenalizationStrategy> > d_brinkman_force;

    /*!
     * Temporary storage variables that contain intermediate quantities
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_temp_sc_var;
    int d_temp_sc_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_temp_cc_var;
    int d_temp_cc_idx;

    /*
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_U_current_idx, d_U_new_idx, d_U_scratch_idx;
    int d_U_old_current_idx, d_U_old_new_idx, d_U_old_scratch_idx;
    int d_P_current_idx, d_P_new_idx, d_P_scratch_idx;
    int d_F_current_idx, d_F_new_idx, d_F_scratch_idx;
    int d_Q_current_idx, d_Q_new_idx, d_Q_scratch_idx;
    int d_N_old_current_idx, d_N_old_new_idx, d_N_old_scratch_idx;
    int d_mu_current_idx, d_mu_new_idx, d_mu_scratch_idx;

    /*
     * Patch data descriptor indices for all "plot" variables managed by the
     * integrator.
     *
     * Plot variables have one context: current.
     */
    int d_U_cc_idx, d_F_cc_idx, d_Omega_idx, d_Div_U_idx, d_EE_idx;

    /*
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    int d_Omega_Norm_idx, d_U_regrid_idx, d_U_src_idx, d_indicator_idx, d_F_div_idx;
    int d_velocity_C_idx, d_velocity_L_idx, d_velocity_D_idx, d_velocity_D_cc_idx, d_pressure_D_idx;
    int d_velocity_rhs_C_idx, d_velocity_rhs_D_idx, d_pressure_rhs_D_idx;
    int d_mu_interp_idx;
    int d_N_full_idx;

    /*
     * Persistent patch data indices for the density and viscosity used in the linear operators
     */
    int d_mu_linear_op_idx, d_mu_interp_linear_op_idx, d_rho_linear_op_idx;

    /*
     * Variables to indicate if either rho or mu is constant.
     */
    bool d_rho_is_const = false, d_mu_is_const = false;

    /*
     * Variable to indicate the type of interpolation to be done for rho and mu.
     */
    IBTK::VCInterpType d_rho_vc_interp_type, d_mu_vc_interp_type;

    /*
     * Variable to indicate the scaling factor used for A, p and u_rhs.
     */
    SAMRAI::tbox::Array<double> d_A_scale;

    /*
     * Variable to set how often the preconditioner is reinitialized.
     */
    int d_precond_reinit_interval = 1;

    /*
     * Objects to set initial condition for density and viscosity when they are maintained by the fluid integrator.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_rho_init_fcn, d_mu_init_fcn;

    /*
     * Boundary condition objects for viscosity, which is provided by an appropriate advection-diffusion
     * integrator
     * or set by the fluid integrator.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_mu_bc_coef = nullptr;

    /*
     * Variable to keep track of a transported viscosity variable maintained by an advection-diffusion integrator
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_mu_adv_diff_var;

    /*
     * Index to track which advection diffusion integrator maintains the viscosity.
     */
    unsigned int d_mu_adv_diff_idx = 0;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSVCStaggeredHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSVCStaggeredHierarchyIntegrator(const INSVCStaggeredHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSVCStaggeredHierarchyIntegrator& operator=(const INSVCStaggeredHierarchyIntegrator& that) = delete;

    /*!
     * Preprocess the operators and solvers used by the hierarchy integrator.
     */
    void preprocessOperatorsAndSolvers(double current_time, double new_time);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_INSVCStaggeredHierarchyIntegrator

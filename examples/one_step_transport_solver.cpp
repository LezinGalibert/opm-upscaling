#include <config.h>
#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <dune/common/version.hh>

#include <opm/porsol/euler/EulerUpstreamImplicit.hpp>
#include <opm/porsol/common/SimulatorTraits.hpp>

#include <opm/common/utility/parameters/ParameterGroup.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/upscaling/UpscalerBase.hpp>



// ------- Typedefs -------
namespace Opm {
    template <class IsotropyPolicy>
    struct Implicit
    {
        template <class GridInterface, class BoundaryConditions>
        struct TransportSolver
        {
            enum { Dimension = GridInterface::Dimension };
            typedef typename IsotropyPolicy::template ResProp<Dimension>::Type RP;

            typedef EulerUpstreamImplicit<GridInterface,
                                          RP,
                                          BoundaryConditions> Type;

        };
    };

    typedef SimulatorTraits<Isotropic, Implicit> UpscalingTraitsBasicImplicit;
}

using namespace Opm;

typedef UpscalerBase<UpscalingTraitsBasicImplicit> Super;
typedef BasicBoundaryConditions<true, true> BCs;


typedef Dune::CpGrid GridType;
	enum { Dimension = GridType::dimension };
typedef GridInterfaceEuler<GridType> GridInterface


typedef UpscalingTraitsBasicImplicit::TransportSolver<GridInterface, Super::BCs>::Type TransportSolver;
typedef UpscalingTraitsBasicImplicit::FlowSolver<GridInterface, BCs>::Type FlowSolver;

typedef UpscalingTraitsBasicImplicit::ResProp<Dimension>::Type ResProp;


int main (int argc, char** argv){
    // Initialize.
	Opm::ParameterGroup param(argc, argv);

	int num_sats = param.getDefault("num_sats", 4);
	double min_sat = param.getDefault("min_sat", 0.2);
	double max_sat = param.getDefault("max_sat", 0.8);

	std::vector<double> saturations
	saturations.resize(num_sats);
	for (int i = 0; i < num_sats; ++i) {
		double factor = num_sats == 1 ? 0 : double(i)/double(num_sats - 1);
		saturations[i] = (1.0 - factor)*min_sat + factor*max_sat;
	}

    int flow_direction = param.getDefault("flow_direction", 0);
    bool start_from_cl = param.getDefault("start_from_cl", true);
    double stepsize= Opm::unit::convert::from(param.getDefault("init_stepsize", init_stepsize), Opm::unit::day);

    Opm::SparseVector<double> injection(num_cells);

    float pdrop = param.getDefault("pdrop", 1e2);
    Dune::FieldVector<double, 3> gravity(0.0);

    Dune::CpGrid grid; //init in **setupGridAndProps**
    GridInterface ginterf; //init in **ginterf_.init**
    BoundaryConditionType bctype = param.get<int>("boundary_condition_type");
    bctype_ = static_cast<BoundaryConditionType>(bct); //init like in **UpscalerBase<Traits>::initImpl**
    bool twodim_hack = param.getDefault("2d_hack", twodim_hack); //init like in **UpscalerBase<Traits>::initImpl**
    BCs bcond; //init in **setupGridAndProps**

    ResProp res_prop; //init in **setupGridAndProps**
    std::vector<double> src(num_cells, 0.0);
    int residual_tolerance = param.getDefault("residual_tolerance", residual_tolerance); //init like in **UpscalerBase<Traits>::initImpl**
    int linsolver_verbosity = param.getDefault("linsolver_verbosity", linsolver_verbosity); //init like in **UpscalerBase<Traits>::initImpl**
    int linsolver_type = param.getDefault("linsolver_type", linsolver_type); //init like in **UpscalerBase<Traits>::initImpl**



    // Set up initial saturation profile.
    std::vector<double> initial_saturation;
    // Initializes initial saturations uniformly
    initial_saturation.resize(num_cells, saturations[i]);


    // take out grid/res_prop initializers?
    setupGridAndProps(temp_param, grid_, res_prop_);
    ginterf_.init(grid_);

    // Set up boundary conditions.
    setupUpscalingConditions(ginterf, bctype, flow_direction, pdrop, saturations[i], twodim_hack, bcond);

    transport_solver_.initObj(ginterf, res_prop, bcond);

    // Run pressure solver.
    flow_solver_.solve(res_prop, initial_saturation, >bcond, src, residual_tolerance, linsolver_verbosity, linsolver_type);

    bool converged = transport_solver.transportSolve(initial_saturation, stepsize, gravity, flow_solver_.getSolution(), injection);

    return 1;

}

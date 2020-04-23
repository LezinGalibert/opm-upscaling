// how to init saturations? CL or uniform?
// how to import parameters
// how to init ginterf & res_prop
// how to implement loop: loss function, stop criteria, what to do in case of non convergence
//

std::vector<double> saturations = INPUT FROM PARAMETERFILE
bool init_from_cl
float pdrop
int flow_direction

// ------- Typedefs -------
typedef typename Traits::template TransportSolver<GridInterface, typename Super::BCs>::Type TransportSolver;
typedef typename Traits::template FlowSolver<GridInterface, BCs>::Type FlowSolver;



// Set up initial saturation profile.
std::vector<double> saturation = initial_saturation;

if init_from_cl {
    // Initializes initial saturations at capillary limit (recommended? speeds up solver)
	upscaler.setToCapillaryLimit(saturations[i], initial_saturation);
}
else {
	// Initializes initial saturations uniformly
	initial_saturation.resize(num_cells, saturations[i]);
}



// Set up boundary conditions.
setupUpscalingConditions(this->ginterf_, this->bctype_, flow_direction,
		pdrop, saturations[i], this->twodim_hack_, this->bcond_);

transport_solver_.initObj(this->ginterf_, this->res_prop_, this->bcond_);



// Run pressure solver.
this->flow_solver_.solve(this->res_prop_, saturation, this->bcond_, src,
                         this->residual_tolerance_, this->linsolver_verbosity_, this->linsolver_type_);

//Need for following line??
double max_mod = this->flow_solver_.postProcessFluxes();
std::cout << "Max mod = " << max_mod << std::endl;



while ((!stationary) && (it_count < max_it_)) { // && transport_cost < max_transport_cost_)
    // Run transport solver.
    std::cout << "Running transport step " << it_count << " with stepsize "
              << stepsize/Opm::unit::year << " years." << std::endl;
    bool converged = transport_solver_.transportSolve(saturation, stepsize, gravity,
                                                      this->flow_solver_.getSolution(), injection);

    // Run pressure solver.
    if (converged) {
        init_saturation = saturation;

        CALCULATE STOP CRITERIA:
		 - saturation distribution?
		 - flow rate?
		 - other physical metrics?
		 - which norms: L1, L2, inf, other?

    }
    else (!converged){
    	DO SOMETHING WITH CRITERIA OR TIMESTEP --> reloop
    }
















    }

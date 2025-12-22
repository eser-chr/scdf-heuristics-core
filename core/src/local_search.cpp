#include <iostream>
#include "solvers.hpp"
#include "structures.hpp"


Solution LS::local_search(
    const Instance &I,
    const Solution &initial_sol,
    const Neighborhood::NeighborhoodFactory &neigh_factory,
    StepFunction::Func step_function,
    StoppingCriterion &criterion,
    int* iteration_ptr)
{

    Solution sol = initial_sol; // copy
    double f = utils::objective(I, sol);
    size_t iteration = 0;
    static std::mt19937 rng(std::random_device{}());

    criterion.reset();

    while (!criterion(iteration, f))
    {
        // build neighborhood for current sol.
        // Use of ptr bcs Neighborhood is abstract and need to rebuild many times
        // and thus not use of reference.

        auto neigh = neigh_factory(I, sol);
        auto mov = step_function(*neigh, rng); // best/first move in THIS neighborhood

        if (!mov.has_value())
            break; // local optimum w.r.t. this neighborhood

        sol = neigh->apply(*mov);
        f = utils::objective(I, sol);
        ++iteration;
    }
    if (iteration_ptr != nullptr){
        *iteration_ptr = iteration;
    }

    return sol;
}

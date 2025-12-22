#include <vector>
#include <memory>
#include <functional>
#include "neighborhoods.hpp"
#include "solvers.hpp"
#include "structures.hpp"


Solution VND::vnd(
    const Instance &I,
    const Solution &initial_sol,
    const Neighborhood::NeighborhoodFactories &neighborhood_factories,
    StepFunction::Func step_function,
    StoppingCriterion &stopping_criterion,
    int *iteration_ptr)
{
    Solution sol = initial_sol;
    double f = utils::objective(I, sol);

    size_t i = 0;
    size_t K = (size_t)neighborhood_factories.size();

    while (i < K)
    {
        Neighborhood::NeighborhoodFactory single_neigh = neighborhood_factories[i];
        Solution new_sol = LS::local_search(
            I,
            sol,
            single_neigh,
            step_function,
            stopping_criterion);

        double f_new = utils::objective(I, new_sol);
        i++;
        if (f_new < f)
        {
            sol = new_sol;
            f = f_new;
            i = 0;
        }
    }

    if (iteration_ptr != nullptr)
        *iteration_ptr = i;

    return sol;
}

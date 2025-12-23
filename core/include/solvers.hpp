#pragma once
#include <vector>
#include "structures.hpp"
#include "neighborhoods.hpp"
#include "step_function.hpp"
#include "stopping_criteria.hpp"

namespace DC // Deterministic Construction
{

    Solution construction(
        const Instance &I);
};
namespace RC // Random Construction
{
    Solution construction(
        const Instance &I,
        double lamda);
};

namespace BS
{
    struct BeamState
    {
        double score;
        std::vector<int> route;
        int cargo;
        std::vector<int> active;
        std::vector<int> remaining;
    };

    Solution beam_search(const Instance &I, double a, int beam_width = 5);

};
namespace LS
{
    Solution local_search(
        const Instance &I,
        const Solution &initial_sol,
        const Neighborhood::NeighborhoodFactory &neigh_factories,
        StepFunction::Func step_function,
        StoppingCriterion &criterion,
        int *iteration_ptr = nullptr);
};

namespace VND
{

    Solution vnd(
        const Instance &I,
        const Solution &initial_sol,
        const Neighborhood::NeighborhoodFactories &neighborhood_factories,
        StepFunction::Func step_function,
        StoppingCriterion &stopping_criterion,
        int *iteration_ptr = nullptr);

};
namespace GRASP // Replace with the real randomized constructor
{

    Solution randomized_constructor_simple(
        const Instance &I,
        double a,
        double alpha);

    Solution grasp(
        const Instance &I,
        std::function<Solution(const Instance &)> randomized_constructor,
        const Neighborhood::NeighborhoodFactories &neighborhoods,
        StepFunction::Func step_function,
        StoppingCriterion &stopping_outer,
        StoppingCriterion &stopping_local,
        int *iteration_ptr = nullptr);

    Solution meta_grasp( // Use of VND
        const Instance &I,
        std::function<Solution(const Instance &)> randomized_constructor,
        const Neighborhood::NeighborhoodFactories &neighborhoods,
        StepFunction::Func step_function,
        StoppingCriterion &stopping_outer,
        StoppingCriterion &stopping_local,
        int *iteration_ptr = nullptr);

};

namespace SA
{
    Solution simulated_annealing(
        const Instance &I,
        const Solution &initial_sol,
        const Neighborhood::NeighborhoodFactories &neighborhood_factories,
        double T_start,
        double T_end,
        double cooling,
        StepFunction::Func step_function,
        StoppingCriterion &stopping_criterion,
        int *iteration_ptr = nullptr);
};

namespace GA{
    Solution genetic_algorithm(Instance const& I);
}
namespace AC{};
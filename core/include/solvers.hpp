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
        int cargo;
        double score;
        std::vector<int> route;
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

namespace LN
{
    // The idea is to encode the initial solution to requests. Then remove requests by
    // finding the most heavy one for each route. This is done by rebuilding the route using
    // beam search. e.g rA: 1,2,3,4,5 -> create_rotue(1,2,3,4). Try for all combo and find the heaviest request.

    // Then we remove C requests.
    // Next we need to add x requests in the encoding till we reach gamma deliveries.
    // For each track and possible request we calculate how much larger the route would be.
    // We take the ideal request for each track and build a first list.

    // A->5
    // B->5
    // C->3

    // Proposed additions = [(A,5), (B, 5), (C, 3), ...]
    // If that list has more than x DIFFERENT requests we are good. We sort the list and add them progressively.
    // Possible implementation as an unordered_set (double requests are removed automatically)
    // then we pick random routes/request to add to the solution.

    // If it has less (bcs some tracks propose the same request) then we add only as many not duplicates and then repeat the
    // process until x requests have been fulfilled.

    // Tracks x Possible requests ->Beam Searches e.g 5x2 = 10 per first iter, less in the next ones.
    struct BestSolution
    {
        double objective;
        Solution sol;
    };
    struct RequestPair
    {
        int request_removed;
        int vehicle;
        double delta;
        std::vector<int> rest_requests;
    };
    RequestPair find_heaviest_request_in_route(Instance const &I, Encoding const &encoding, int vehicle, int beam_width);
    RequestPair find_best_request_to_add(Instance const &I, Encoding const &encoding, int vehicle, int beam_width);
    Encoding apply_removal(Instance const &I, Encoding const &encoding, std::vector<RequestPair> const &to_be_removed);
    Encoding apply_addition(Instance const &I, Encoding const &encoding, std::vector<RequestPair> const &to_be_removed);
    Encoding remove_requests(Instance const &I, Encoding const &encoding, int k, int beam_width);
    Encoding append_requests(Instance const &I, Encoding const &encoding, int k, int beam_width);
    Solution large_neighborhood(Instance const &I, Solution const &sol, int k, size_t iters, int bw_remove, int bw_append);

};

namespace GA
{
    struct BestSolution
    {
        double objective;
        Solution sol;
    };
    BestSolution get_best_solution(Instance const &I, std::vector<Encoding> const& encodings);
    std::vector<Encoding> generate_initial_population(Instance const &I, int k1);
    std::vector<Encoding> reproduce(Instance const &I, std::vector<Encoding> const& parents);
    // If a combination results in a request being delivered by  strictly two vehicles this function resolves that.
    // It uses a uniform distribution.
    // It is meant to be used only inside the plus operator
    std::vector<int> select_indices_next_generation(Instance const &I, std::vector<Encoding> const &population, int k1);
    void mutate(Instance const &I, std::vector<Encoding> &population, int k2);
   
    /**
     * @param k1 Number of population to be reproduced
     * @param k2 Number of mutated requests. 0= no mutation
     * @param iters Number of iterations-generations
     * @param beam_width Beam width of the beam search used to create the route
     */
    Solution genetic_algorithm(Instance const &I, int k1, int k2, int iters, int beam_width);
};
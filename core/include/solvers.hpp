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

namespace LN
{
    struct ModifyRequestInfo
    {
        int vehicle_id; // vehicle id as in the solution. i.e vehicle_id = k -> sol.route[k] the corresponding route.
        int request_id; // request id as in the Instance
        int p_idx;      // pickup index in the route of that vehicle
        int d_idx;      // delivery index in the route of that vehicle
        double delta;   // difference in distance

        // ModifyRequestInfo(int vehicle_id, int request_id, int p_idx, int d_idx, double delta) : vehicle_id(vehicle_id), request_id(request_id), p_idx(p_idx), d_idx(d_idx), delta(delta) {}
    };
    struct ModifyRequestInfoFull
    {
        int vehicle_id;         // vehicle id as in the solution. i.e vehicle_id = k -> sol.route[k] the corresponding route.
        int request_id;         // request id as in the Instance
        std::vector<int> route; // The proposed route.
        double delta;           // difference in distance

        // ModifyRequestInfo(int vehicle_id, int request_id, int p_idx, int d_idx, double delta) : vehicle_id(vehicle_id), request_id(request_id), p_idx(p_idx), d_idx(d_idx), delta(delta) {}
    };

    enum class Modification
    {
        Remove = 0,
        Append = 1
    };
    // std::vector<int> apply_modification(Modification modification, ModifyRequestInfo const &info);
    std::vector<int> apply_modification(Instance const &I, Solution const &sol, Modification modification, ModifyRequestInfo const &info);
    double calc_delta_of_request(Instance const &I, Solution const &sol, Modification modification, ModifyRequestInfo const &info);

    Solution large_neighborhood(Instance const &I, Solution const &sol, int k);
};

namespace GA
{
    class Encoding
    {
        using dna_t = std::vector<std::vector<bool>>;
        // Check if you can use a single vector with offset
        dna_t dna;

    public:
        Encoding() = default;
        Encoding(dna_t &&dna);
        Encoding(Instance const &I, Solution const &sol);

        // Decoding process. Uses greedy algo from DC construction to
        // reconstruct a solution from an encoding
        Solution to_sol(Instance const &I) const;
        bool is_encoding_correct(Instance const &I) const;
        Encoding operator+(Encoding const &other) const;
        Encoding add(Instance const& I, Encoding const &other) const;
        // void _resolve_degeneracies(dna_t & child);
        void set_vehicle_for_request(int vehicle, int request);
    };

    std::vector<Encoding> generate_initial_population(Instance const &I, int k1);
    std::vector<Encoding> reproduce(Instance const &I, std::vector<Encoding> parents);
    // If a combination results in a request being delivered by  strictly two vehicles this function resolves that.
    // It uses a uniform distribution.
    // It is meant to be used only inside the plus operator
    std::vector<int> select_indices_next_generation(Instance const &I, std::vector<Encoding> const &population, int k1);
    void mutate(Instance const &I, std::vector<Encoding> &population, int k2);
    /**
     * @param k1 Number of population to be reproduced
     * @param k2 Number of mutated requests. 0= no mutation
     * @param iters Number of iterations-generations
     */
    Solution genetic_algorithm(Instance const &I, int k1, int k2, int iters);

};
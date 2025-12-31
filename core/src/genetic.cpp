#include <cassert>
#include <iostream>
#include <limits>
#include "structures.hpp"
#include "solvers.hpp"

std::vector<Encoding> GA::generate_initial_population(Instance const &I, int k1)
{
    std::vector<Encoding> tortn;
    tortn.reserve(k1);

    if (k1 < 3)
    {
        std::cout << " I need a population of more than 3 to operate. \n";
        std::abort();
    }

    size_t size_of_dc = (size_t)k1 / 3;
    size_t counter = 0;

    for (; counter < size_of_dc; counter++)
    {
        // Solution tmp = BS::beam_search(I, 0.9, 5);

        Solution tmp = DC::construction(I);
        Encoding encoding(I, tmp);
        tortn.push_back(std::move(encoding));
    }
    for (; counter < k1; counter++)
    {
        Solution tmp = BS::beam_search(I, 0.9, 5);
        Encoding encoding(I, tmp);
        tortn.push_back(std::move(encoding));
    }
    return tortn;
}

std::vector<Encoding> GA::reproduce(Instance const &I, std::vector<Encoding> const &parents)
{
    std::vector<Encoding> tortn;
    tortn.reserve((parents.size() * (parents.size() - 1)) / 2);

    for (size_t i = 0; i < parents.size(); i++)
    {
        for (size_t j = i + 1; j < parents.size(); j++)
        {

            Encoding new_encoding = parents[i].add(I, parents[j]);
            // assert(new_encoding.total_num_of_requests()==I.gamma);
            if (new_encoding.total_num_of_requests() != I.gamma)
            {
                std::cout << new_encoding.total_num_of_requests() << I.gamma << std::endl;
            }
            tortn.emplace_back(std::move(new_encoding));
            // tortn.emplace_back(parents[i] + parents[j]);
        }
    }
    return tortn;
}

void GA::mutate(Instance const &I, std::vector<Encoding> &population, int k2)
{
    if (k2 == 0)
        return;

    static thread_local std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> request_dist(0, I.n - 1);
    std::uniform_int_distribution<int> vehicle_dist(0, I.nK - 1);

    for (auto &encoding : population)
    {
        for (int i = 0; i < k2; i++)
        {
            int request_to_mutate = request_dist(gen);
            int vehicle_to = vehicle_dist(gen);
            encoding.set_vehicle_for_request(vehicle_to, request_to_mutate);
        }
    }
}
std::vector<int> GA::select_indices_next_generation(Instance const &I, std::vector<Encoding> const &population, int k1)
{
    // std::vector<double> objectives(population.size());
    std::vector<double> objectives;
    objectives.reserve(population.size());

    for(auto const& enc: population){
        Solution tmp = enc.to_sol(I);
        objectives.push_back(utils::objective(I, tmp));
    }
    auto indices = numerical::argsort(objectives);
    return std::vector<int>(indices.begin(), indices.begin()+k1);
}

GA::BestSolution GA::get_best_solution(Instance const &I, std::vector<Encoding> const& encodings)
{
    Solution sol;
    double objective = std::numeric_limits<double>::infinity();
    BestSolution best_sol;

    for (size_t i = 0; i < encodings.size(); i++)
    {
        Solution tmp = encodings[i].to_sol(I);
        double obj = utils::objective(I, tmp);
        if (obj < objective)
        {
            sol = tmp;
            objective = obj;
        }
    }
    best_sol.objective = objective;
    best_sol.sol = sol;
    return best_sol;
}

Solution GA::genetic_algorithm(Instance const &I, int k1, int k2, int iters, int beam_width)
{
    assert(beam_width > 0);
    assert(iters > 0);
    auto population = generate_initial_population(I, k1);
    BestSolution best_sol = get_best_solution(I, population);
    for (size_t iter = 0; iter < iters; iter++)
    {
        assert(population.size() == k1);
        auto offsprings = reproduce(I, population);
        assert(offsprings.size() >= k1);
        mutate(I, offsprings, k2);
        auto indices_of_survivors = select_indices_next_generation(I, offsprings, k1);
        std::vector<Encoding> new_population;
        new_population.reserve(k1);
        for (size_t i = 0; i < k1; i++)
        {
            new_population.push_back(std::move(offsprings[indices_of_survivors[i]]));
        }
        BestSolution new_best_sol = get_best_solution(I, new_population);

        if (new_best_sol.objective < best_sol.objective)
        {
            best_sol = std::move(new_best_sol);
        }

        std::swap(population, new_population); // Delete the older data
    }
    assert(best_sol.sol.is_solution_feasible(I));
    return best_sol.sol;
}

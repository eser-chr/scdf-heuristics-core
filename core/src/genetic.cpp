#include <cassert>
#include <iostream>
#include "structures.hpp"
#include "solvers.hpp"

// Decleration. COnsider put it in a include file.
// We use it also in large neighborhood. Not good for readability. Error prone.
std::vector<int> build_route_greedy(
    const Instance &I,
    const std::vector<int> &reqs);

GA::Encoding::Encoding(Instance const &I, Solution const &sol)
{

    dna.reserve(I.nK);
    for (size_t j = 0; j < I.nK; j++)
    {
        dna.emplace_back(I.n, false);
    }

    for (size_t j = 0; j < sol.routes.size(); j++)
    {
        auto const route = sol.routes[j];
        for (auto node : route)
        {
            assert(node != 0);
            if (node <= I.n)
            {
                dna[j][node - 1] = true;
            }
        }
    }
}

void GA::Encoding::set_vehicle_for_request(int vehicle, int request)
{
    int rows = dna.size();
    int cols = dna[0].size();

    for (size_t i = 0; i < rows; i++)
    {
        dna[i][request] = (vehicle == (int)i);
    }
}

bool GA::Encoding::is_encoding_correct(Instance const &I) const
{
    size_t rows = dna.size();
    size_t cols = dna[0].size();

    for (auto col : dna)
    {
        if (col.size() != cols)
            return false;
    }

    for (size_t j = 0; j < cols; j++)
    {
        // Here we check if a request has more or equal to two routes.
        // 0, 1 are fine.

        bool flag = true;
        for (size_t i = 0; i < rows; i++)
        {
            // If flag was already false i.e a previous route had this request
            // and  also the current route has this then something is wrong
            if (!flag && dna[i][j])
                return false;
            flag = !dna[i][j];
        }
    }

    if (cols != I.n)
    {
        std::cerr << " cols of encoding do not match requests in Instance" << std::endl;
        return false;
    }

    return true;
}

Solution GA::Encoding::to_sol(Instance const &I) const
{
    Solution new_sol;
    for (auto const &encoding_route : dna)
    {
        assert(encoding_route.size() == I.n);
        std::vector<int> route_requests;
        for (size_t i = 0; i < encoding_route.size(); ++i)
        {
            if (encoding_route[i])
                route_requests.push_back(i);
        }

        auto route = build_route_greedy(I, route_requests);
        new_sol.routes.push_back(route);
    }

    new_sol.compute_cached_values_from_routes(I);
    assert(new_sol.is_solution_feasible(I));
    return new_sol;
}

GA::Encoding::Encoding(dna_t &&other_dna) : dna(std::move(other_dna)) {}

GA::Encoding GA::Encoding::operator+(Encoding const &other) const
{
    static thread_local std::mt19937 rng(std::random_device{}());
    // std::random_device rd;
    // std::mt19937 rng(rd());
    std::uniform_int_distribution<int> dist(0, 1);

    dna_t offspring = dna;
    size_t rows = dna.size();
    size_t cols = dna[0].size();

    for (size_t col = 0; col < cols; col++)
    {
        size_t r = dist(rng); // flip coin
        // std::array<size_t, 2> proposed_rows;
        bool is_set = false;
        int row_set = -1;
        for (size_t row = 0; row < rows; row++)
        {
            // Possible improvement if we test and. If both are true then we can break early.
            bool tmp = dna[row][col] || other.dna[row][col];
            if (!is_set && tmp)
            {
                is_set = tmp;
                row_set = row;
                offspring[row][col] = tmp;
            }
            if (is_set && tmp)
            {
                if (r == 1)
                {
                    // The second offer wins. Remove the first one and
                    // add here
                    offspring[row_set][col] = false;
                    row_set = row;
                    offspring[row_set][col] = true;
                }
                break; // Break since there is not other true value for this request
            }
            // For the rest values offspring[row][col] should be false as in the first dna.
        }
    }

    return {std::move(offspring)};
}

GA::Encoding GA::Encoding::add(Instance const &I, Encoding const &other) const
{
    static thread_local std::mt19937 rng(std::random_device{}());

    size_t rows = dna.size();
    size_t cols = dna[0].size();

    // Step 1: Categorize requests
    std::vector<int> both_parents;    // Delivered by both parents
    std::vector<int> one_parent_only; // Delivered by exactly one parent

    for (size_t col = 0; col < cols; col++)
    {
        bool in_this = false;
        bool in_other = false;
        int row_this = -1, row_other = -1;

        // Check if request is delivered in this parent
        for (size_t row = 0; row < rows; row++)
        {
            if (dna[row][col])
            {
                in_this = true;
                row_this = row;
                break;
            }
        }

        // Check if request is delivered in other parent
        for (size_t row = 0; row < rows; row++)
        {
            if (other.dna[row][col])
            {
                in_other = true;
                row_other = row;
                break;
            }
        }

        if (in_this && in_other)
        {
            both_parents.push_back(col);
        }
        else if (in_this || in_other)
        {
            one_parent_only.push_back(col);
        }
    }

    // Step 2: Determine how many from one_parent_only to include
    // Assuming you have access to gamma (see note below)
    int gamma = I.gamma /* I.gamma - you need to pass this */;
    int needed = gamma - (int)both_parents.size();
    needed = std::max(0, std::min(needed, (int)one_parent_only.size()));

    // Step 3: Randomly select 'needed' requests from one_parent_only
    std::shuffle(one_parent_only.begin(), one_parent_only.end(), rng);

    // Step 4: Build offspring - initialize with all false
    dna_t offspring(rows, std::vector<bool>(cols, false));

    // Add requests from both parents (inherit vehicle assignment randomly)
    for (int col : both_parents)
    {
        int row_this = -1, row_other = -1;

        for (size_t row = 0; row < rows; row++)
        {
            if (dna[row][col])
                row_this = row;
            if (other.dna[row][col])
                row_other = row;
        }

        // Randomly choose which parent's vehicle assignment to use
        std::uniform_int_distribution<int> coin(0, 1);
        int chosen_row = (coin(rng) == 0) ? row_this : row_other;
        offspring[chosen_row][col] = true;
    }

    // Add randomly selected requests from one_parent_only
    for (int i = 0; i < needed; i++)
    {
        int col = one_parent_only[i];

        // Find which parent has this request and which vehicle
        for (size_t row = 0; row < rows; row++)
        {
            if (dna[row][col] || other.dna[row][col])
            {
                offspring[row][col] = true;
                break;
            }
        }
    }

    return {std::move(offspring)};
}

std::vector<GA::Encoding> GA::generate_initial_population(Instance const &I, int k1)
{
    std::vector<Encoding> tortn;
    tortn.reserve(k1);

    for (size_t i = 0; i < k1; i++)
    {
        Solution tmp = RC::construction(I, 0.9);
        Encoding encoding(I, tmp);
        tortn.push_back(std::move(encoding));
    }
    return tortn;
}

std::vector<GA::Encoding> GA::reproduce(Instance const &I, std::vector<Encoding> parents)
{
    std::vector<Encoding> tortn;
    tortn.reserve((parents.size() * (parents.size() - 1)) / 2);

    for (size_t i = 0; i < parents.size(); i++)
    {
        for (size_t j = i + 1; j < parents.size(); j++)
        {
            // Encoding new_encoding =
            tortn.emplace_back(parents[i].add(I,parents[j]));
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
    std::vector<double> objectives(population.size());
    for (size_t i = 0; i < population.size(); i++)
    {
        Solution tmp = population[i].to_sol(I);
        objectives[i] = -utils::objective(I, tmp); // We use minus so we can argsort later from min to max
    }
    auto indices = numerical::argsort(objectives);
    std::vector<int> tortn(k1);
    assert(indices.size() >= k1);

    // tortn.f
    // (indices.begin(), indices.begin()+k1);
    for (size_t i = 0; i < k1; i++)
    {
        tortn[i] = indices[i];
    }

    return tortn;
}

Solution GA::genetic_algorithm(Instance const &I, int k1, int k2, int iters)
{
    auto population = generate_initial_population(I, k1);
    for (size_t iter = 0; iter < iters; iter++)
    {

        assert(population.size() == k1);
        auto offsprings = reproduce(I, population);
        assert(offsprings.size() >= k1);
        mutate(I, offsprings, k2);
        auto indices_of_survivors = select_indices_next_generation(I, offsprings, k1);
        std::vector<Encoding> new_population(k1);
        for (size_t i = 0; i < k1; i++)
        {
            new_population[i] = offsprings[indices_of_survivors[i]];
        }
        std::swap(population, new_population); // Delete the older data
    }

    auto argsort_indices = select_indices_next_generation(I, population, 1);
    Solution tortn = population[argsort_indices[0]].to_sol(I); // Use the best index!
    tortn.compute_cached_values_from_routes(I);
    assert(tortn.is_solution_feasible(I));
    return tortn;
}

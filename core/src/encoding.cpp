#include "solvers.hpp"
#include <cassert>
#include <iostream>
#include <functional>
#include <algorithm>

std::vector<int> create_track_route(Instance const &I, int beam_width, std::vector<int> const &requests);

Encoding::Encoding(Instance const &I, Solution const &sol)
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

void Encoding::set_vehicle_for_request(int vehicle, int request)
{
    int rows = dna.size();
    int cols = dna[0].size();

    for (size_t i = 0; i < rows; i++)
    {
        dna[i][request] = (vehicle == (int)i);
    }
}

bool Encoding::is_encoding_correct(Instance const &I) const
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
Solution Encoding::to_sol(Instance const &I, int beam_width) const
{
    if(!cached_solution){
        cached_solution = std::move(_compute_solution(I, beam_width));
    }

    return *cached_solution;
}

Solution Encoding::_compute_solution(Instance const &I, int beam_width) const{


    assert(beam_width>0);
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

        auto route = create_track_route(I, beam_width, route_requests);
        new_sol.routes.push_back(route);
    }

    new_sol.compute_cached_values_from_routes(I);
    // assert(new_sol.is_solution_feasible(I));
    return new_sol;
}

int Encoding::total_num_of_requests() const
{
    size_t rows = dna.size();
    size_t cols = dna[0].size();

    int counter = 0;
    for (size_t col = 0; col < cols; col++)
    {

        for (size_t row = 0; row < rows; row++)
        {
            if (dna[row][col])
            {
                counter++;
                break;
            }
        }
    }
    return counter;
}

Encoding::Encoding(dna_t &&other_dna) : dna(std::move(other_dna)) {}

Encoding Encoding::operator+(Encoding const &other) const
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

Encoding Encoding::add(Instance const &I, Encoding const &other) const
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

    int needed = I.gamma - (int)both_parents.size();
    std::shuffle(one_parent_only.begin(), one_parent_only.end(), rng);
    std::shuffle(both_parents.begin(), both_parents.end(), rng);

    if (needed < 0)
    {

        both_parents = std::vector<int>(both_parents.begin(), both_parents.begin() + both_parents.size() + needed);
    }

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

std::vector<int> Encoding::get_requests_of_route(int route) const
{
    size_t rows = dna.size();
    assert(route < rows);
    size_t cols = dna[0].size();

    std::vector<int> requests;
    requests.reserve(cols);
    for (size_t col = 0; col < cols; col++)
    {
        if (dna[route][col])
            requests.push_back((int)col);
    }

    return requests;
}

Encoding::dna_t const &Encoding::get_dna() const
{
    return dna;
}
Encoding::dna_t Encoding::get_dna_copy() const
{
    return dna;
}

std::vector<int> Encoding::get_non_delivered_requests() const{
    size_t rows = dna.size();
    size_t cols = dna[0].size();

    std::vector<int> non_delivered;

    for(size_t col =0; col<cols;++col){
        bool is_delivered = false;
        for(size_t row = 0; row<rows; row++){
            if(dna[row][col]){
                is_delivered=true;
                break;
            }
        }

        if(is_delivered)
            continue;
        
        non_delivered.push_back((int)col);


    }

    return non_delivered;
}

int Encoding::get_num_vehicles() const{
    return dna.size();
}
int Encoding::get_num_requests() const{
    return dna[0].size();
}
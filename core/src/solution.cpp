#include <fstream>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "structures.hpp"

void Solution::write_solution(const std::string &path, const std::string &instance_name) const
{
    std::cout << "nK = " << routes.size() << std::endl;
    std::ofstream f(path);
    if (!f)
    {
        throw std::runtime_error("Could not open file for writing: " + path);
    }

    f << instance_name << "\n";

    for (const auto &route : routes)
    {
        if (route.empty())
        {
            f << "\n"; // blank line for empty route
            std::cout << "One empty route" << std::endl;
        }
        else
        {
            for (size_t i = 0; i < route.size(); ++i)
            {
                f << route[i];
                if (i + 1 < route.size())
                    f << " ";
            }
            f << "\n";
        }
    }
}

bool Solution::is_solution_feasible(Instance const &I)
{
    // Are gamma requests fullfiled?
    // Are pickups delivered ?
    // Is capacity ever larger than allowed?
    int full_req = 0;
    int total_size_of_nodes = 0;
    for (auto const &route : routes)
    {
        int capacity = 0;
        total_size_of_nodes += route.size();

        for (auto node : route)
        {
            assert(node != 0);
            capacity += I.load_change[node];
            if (capacity > I.C)
            {
                std::cerr << "While asserting solution. Capacity found to be incorrect" << std::endl;
                return false;
            }

            if (node > 0 && node < I.n + 1)
            {
                bool is_deliver_node_in_route =
                    std::find(route.begin(), route.end(), I.n + node) != route.end();

                if (!is_deliver_node_in_route)
                {
                    std::cerr << "While asserting solution. Deliver was not inside the route " << node << std::endl;
                    return false;
                }
                full_req++;
            }

            if (node >= I.n + 1 && node < 2 * I.n + 1)
            {
                bool is_pickup_inside_the_route =
                    std::find(route.begin(), route.end(), node - I.n) != route.end();

                if (!is_pickup_inside_the_route)
                {
                    std::cerr << "While asserting solution. Pickup was not inside the route" << node << std::endl;
                    return false;
                }
            }
        }
    }

    assert(total_size_of_nodes % 2 == 0);
    assert(total_size_of_nodes / 2 == full_req);

    if (full_req < I.gamma)
    {
        std::cerr << "While asserting solution. Not enough requests delivered " << full_req << "  -  " << I.gamma << std::endl;
        return false;
    }
    if (full_req > I.gamma)
        std::cout << "While asserting solution. Not an error, but more than gamma requests delivered" << std::endl;

    return true;
}

void Solution::compute_cached_values_from_routes(Instance const &I)
{

    assert(I.nK == routes.size());
    assert(routes_distances.size() == 0 || routes_distances.size() == routes.size());

    if (routes_distances.size() == 0)
        routes_distances.resize(I.nK, 0.0);

    total_distance = 0.0;
    sum_of_squares = 0.0;
    for (size_t i = 0; i < routes.size(); i++)
    {
        double tmp = utils::calc_route_distance(I, routes[i]);
        routes_distances[i] = tmp;
        total_distance+=tmp;
        sum_of_squares+= tmp*tmp;
    }
    fairness = I.fairness;
}
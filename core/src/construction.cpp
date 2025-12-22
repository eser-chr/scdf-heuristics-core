/*
The idea of this construction heristic is to clusterize the requests into regions and then add them into routes.
Another possible idea is to use a convex hull for clusters. Maybe in Python.
*/

#include <vector>
#include <algorithm>
#include <random>
#include <limits>
#include <iostream>
#include "structures.hpp"
#include "solvers.hpp"

std::vector<int> build_route_greedy(
    const Instance &I,
    const std::vector<int> &reqs)
{
    std::vector<int> unpicked = reqs;
    std::vector<int> active;
    std::vector<int> route;
    int cargo = 0;
    int last = 0; // depot

    auto can_i_pickup_request = [&](int r)
    { return cargo + I.demands[r] <= I.C; };

    while (!unpicked.empty() || !active.empty())
    {
        int best_node = -1, best_req = -1;
        bool pick = false;
        double best_d{std::numeric_limits<double>::infinity()};

        // pickups
        for (int r : unpicked)
        {
            if (!can_i_pickup_request(r))
                continue;
            int pn = 1 + r;
            double d = I.dist[last][pn]; //distance from depot pseudo heuristic
            if (d < best_d)
            {
                best_d = d;
                best_node = pn;
                best_req = r;
                pick = true;
            }
        }

        // deliveries
        for (int r : active)
        {
            int dn = 1 + I.n + r;
            double d = I.dist[last][dn];
            if (d < best_d)
            {
                best_d = d;
                best_node = dn;
                best_req = r;
                pick = false;
            }
        }

        // Finish by delivering all of them
        if (best_node == -1)
        {
            double bd{std::numeric_limits<double>::infinity()};
            for (int r : active)
            {
                int dn = 1 + I.n + r;
                double d = I.dist[last][dn];
                if (d < bd)
                {
                    bd = d;
                    best_node = dn;
                    best_req = r;
                    pick = false;
                }
            }
        }

        route.push_back(best_node);

        // cargo logisitcs
        if (pick)
        {
            cargo += I.demands[best_req];
            active.push_back(best_req);
            unpicked.erase(std::remove(unpicked.begin(), unpicked.end(), best_req),
                           unpicked.end());
        }
        else
        {
            cargo -= I.demands[best_req];
            active.erase(std::remove(active.begin(), active.end(), best_req),
                         active.end());
        }

        last = best_node;
    }

    return route;
}
// Returns the first gamma requests, sorted by a simple cost heuristic.
auto select_gamma_requests(Instance const &I)
{
    std::vector<double> cost(I.n);

    for (size_t i = 0; i < I.n; ++i)
        cost[i] = I.demands[i] * I.dist[1 + i][1 + i + I.n];

    auto argsort = numerical::argsort(cost);

    std::vector<int> to_rtn;
    to_rtn.reserve(I.gamma);
    for (int k = 0; k < I.gamma; ++k)
        to_rtn.push_back(static_cast<int>(argsort[k]));

    return to_rtn;
}

Solution DC::construction(
    const Instance &I)
{

    auto indices_of_requests_to_serve = select_gamma_requests(I);

    if (indices_of_requests_to_serve.size() != I.gamma)
    {
        throw std::runtime_error("Assertion failed");
    }

    std::vector<int> assign = clusters::balanced_kmeans(I, indices_of_requests_to_serve, 20, 20);
    gt::Matrix<int> per_track(I.nK); // per track requests-responibilities


    for (size_t i = 0; i < indices_of_requests_to_serve.size(); i++)
    {
        int req = indices_of_requests_to_serve[i]; // actual request ID
        int k = assign[i];                         // cluster index
        per_track[k].push_back(req);
    }

    gt::Matrix<int> routes;
    routes.reserve(I.nK);
    for (int k = 0; k < I.nK; k++)
    {
        auto route = build_route_greedy(I, per_track[k]);
        routes.push_back(std::move(route));
    }

    Solution sol;
    sol.routes = std::move(routes);
    auto all_distances = utils::all_route_distances(I, sol);
    std::vector<double> sq_distances;
    sq_distances.resize(all_distances.size());
    std::transform(all_distances.begin(), all_distances.end(), sq_distances.begin(), [](auto val)
                   { return val * val; });
    sol.total_distance = std::accumulate(all_distances.begin(), all_distances.end(), 0.0);

    sol.sum_of_squares = std::accumulate(sq_distances.begin(), sq_distances.end(), 0.0);
    sol.routes_distances = utils::all_route_distances(I, sol);
    return sol;
}

#include <cassert>
#include <vector>
#include <iostream>
#include "solvers.hpp"
#include "structures.hpp"

std::vector<int> build_route_greedy(
    const Instance &I,
    const std::vector<int> &reqs);

std::vector<int> LN::apply_modification(Instance const &I, Solution const &sol, Modification modification, ModifyRequestInfo const &info)
{
    auto const &old_route = sol.routes[info.vehicle_id];
    std::vector<int> to_rtn;

    if (modification == Modification::Append)
    {
        to_rtn.reserve(old_route.size() + 2);
        int i = 0;
        for (; i < info.p_idx; i++)
            to_rtn.emplace_back(old_route[i]);
        to_rtn.emplace_back(info.request_id + 1);
        for (; i < info.d_idx; i++)
            to_rtn.emplace_back(old_route[i]);
        to_rtn.emplace_back(info.request_id + I.n + 1);
        for (; i < old_route.size(); i++)
            to_rtn.emplace_back(old_route[i]);
        return to_rtn;
    }
    else
    {
        to_rtn.reserve(old_route.size() - 2);
        for (size_t i = 0; i < old_route.size(); i++)
        {
            if (i == info.p_idx || i == info.d_idx)
                continue;
            to_rtn.emplace_back(old_route[i]);
        }
        return to_rtn;
    }
}

// Can be accelerated possibly if required.
double LN::calc_delta_of_request(Instance const &I, Solution const &sol, Modification modification, ModifyRequestInfo const &info)
{
    auto new_route = apply_modification(I, sol, modification, info);
    double old_distance = utils::calc_route_distance(I, sol, info.vehicle_id);
    double new_distance = utils::calc_route_distance(I, new_route);
    return new_distance - old_distance;
}

std::vector<int> find_requests_in_route(std::vector<int> const &route, int N)
{
    std::vector<int> requests;
    requests.reserve(route.size() / 2);

    for (int node : route)
    {
        assert(node != 0);
        if (node < N + 1)
        {
            requests.emplace_back(node-1);
        }
    }

    return requests;
}

std::pair<int, int> find_pd_indices_of_request(std::vector<int> const &route, int N, int request)
{
    int pidx = -1;
    int didx = -1;
    for (size_t i = 0; i < route.size(); i++)
    {
        if (route[i] == request + 1)
            pidx = i;
        if (route[i] == request + N + 1)
            didx = i;
    }

    assert(pidx != -1);
    assert(didx != -1);

    return {pidx, didx};
}

std::vector<LN::ModifyRequestInfo> build_remove_request_list(Instance const &I, Solution const &sol)
{
    std::vector<LN::ModifyRequestInfo> tortn;
    for (auto i = 0; i < sol.routes.size(); ++i)
    {
        auto const &route = sol.routes[i];
        auto const requests = find_requests_in_route(route, I.n);
        for (auto request : requests)
        {
            auto [pidx, didx] = find_pd_indices_of_request(route, I.n, request);
            LN::ModifyRequestInfo mri{i, request, pidx, didx};
            double delta = LN::calc_delta_of_request(I, sol, LN::Modification::Remove, mri);
            mri.delta = delta;
            tortn.push_back(mri);
        }
    }
    return tortn;
}

std::vector<LN::ModifyRequestInfoFull> build_append_request_list(Instance const &I, Solution const &partial_sol)
{
    std::vector<bool> requests_fulfilled(I.n, false);

    for (auto const &route : partial_sol.routes)
    {
        auto requests = find_requests_in_route(route, I.n);
        for (auto request : requests)
            requests_fulfilled[request] = true;
    }

    std::vector<LN::ModifyRequestInfoFull> to_rtn;
    for (size_t request_id = 1; request_id < I.n + 1; request_id++)
    {
        if (requests_fulfilled[request_id])
            continue;

        // for every vehicle find the best way to put this request. Possible rebuild of path.
        for (int vehicle_id = 0; vehicle_id < I.nK; vehicle_id++)
        {
            auto requests = find_requests_in_route(partial_sol.routes[vehicle_id], I.n);
            requests.push_back(request_id);
            auto new_suggested_route = build_route_greedy(I, requests);

            double old_distance = utils::calc_route_distance(I, partial_sol.routes[vehicle_id]);
            double new_distance = utils::calc_route_distance(I, new_suggested_route);
            LN::ModifyRequestInfoFull tmp;
            tmp.request_id = request_id;
            tmp.vehicle_id = vehicle_id;
            tmp.route = std::move(new_suggested_route);
            tmp.delta = new_distance - old_distance;

            to_rtn.push_back(std::move(tmp));
        }
    }

    return to_rtn;
}

// std::vector<LN::ModifyRequestInfoFull> build_append_request_list(
//     Instance const &I, 
//     Solution const &partial_sol,
//     int max_proposals_per_request)  // ✅ Limit to top 5 vehicles per request
// {
//     std::vector<bool> requests_fulfilled(I.n, false);

//     for (auto const &route : partial_sol.routes)
//     {
//         auto requests = find_requests_in_route(route, I.n);
//         for (auto request : requests)
//             requests_fulfilled[request] = true;
//     }

//     std::vector<LN::ModifyRequestInfoFull> to_rtn;
    
//     for (size_t request_id = 0; request_id < I.n; request_id++)
//     {
//         if (requests_fulfilled[request_id])
//             continue;

//         // Calculate cost to add this request to each vehicle
//         std::vector<std::pair<int, double>> vehicle_costs;
//         vehicle_costs.reserve(I.nK);
        
//         for (int vehicle_id = 0; vehicle_id < I.nK; vehicle_id++)
//         {
//             auto requests = find_requests_in_route(partial_sol.routes[vehicle_id], I.n);
//             requests.push_back(request_id);
//             auto new_route = build_route_greedy(I, requests);

//             double old_distance = utils::calc_route_distance(I, partial_sol.routes[vehicle_id]);
//             double new_distance = utils::calc_route_distance(I, new_route);
            
//             vehicle_costs.push_back({vehicle_id, new_distance - old_distance});
//         }
        
//         // ✅ Only keep top N best vehicles for this request
//         std::partial_sort(vehicle_costs.begin(), 
//                          vehicle_costs.begin() + std::min(max_proposals_per_request, static_cast<int>(vehicle_costs.size())),
//                          vehicle_costs.end(),
//                          [](auto const &a, auto const &b) { return a.second < b.second; });
        
//         int limit = std::min(max_proposals_per_request, static_cast<int>(vehicle_costs.size()));
//         for (int i = 0; i < limit; i++)
//         {
//             int vehicle_id = vehicle_costs[i].first;
            
//             auto requests = find_requests_in_route(partial_sol.routes[vehicle_id], I.n);
//             requests.push_back(request_id);
//             auto new_route = build_route_greedy(I, requests);
            
//             LN::ModifyRequestInfoFull tmp;
//             tmp.request_id = request_id;
//             tmp.vehicle_id = vehicle_id;
//             tmp.route = std::move(new_route);
//             tmp.delta = vehicle_costs[i].second;
            
//             to_rtn.push_back(std::move(tmp));
//         }
//     }

//     return to_rtn;
// }


Solution LN::large_neighborhood(Instance const &I, Solution const &sol, int k)
{
    Solution partial_sol = sol;
    { // I put this in the block so memory can be freed from deleting the remove nodes_list
        // Adjust so it picks from different routes.

        auto remove_request_list = build_remove_request_list(I, sol);
        assert(remove_request_list.size() >= k);
        std::partial_sort(remove_request_list.begin(), remove_request_list.begin() + k, remove_request_list.end(), [](ModifyRequestInfo const &a, ModifyRequestInfo const &b)
                          { return a.delta > b.delta; });

        for (auto i = 0; i < k; i++)
        {
            auto new_route = apply_modification(I, partial_sol, Modification::Remove, remove_request_list[i]);
            partial_sol.routes[remove_request_list[i].vehicle_id] = new_route;
        }
    }

    for (int i = 0; i < k; i++)
    {
        auto proposals = build_append_request_list(I, partial_sol);
        if (proposals.empty())
        {
            std::cerr << " Proposals are empty why" << std::endl;
            std::abort();
        }
        // std::(proposals.begin(), proposals.end(), );
        auto min_it = std::min_element(proposals.begin(), proposals.end(), [](ModifyRequestInfoFull const &a, ModifyRequestInfoFull const &b)
                                       { return a.delta < b.delta; });
        auto best_proposal = *min_it;
        partial_sol.routes[best_proposal.vehicle_id] = best_proposal.route;

        // apply_modification(I, partial_sol, Modification::Append, proposals[0]);
    }

    assert(partial_sol.is_solution_feasible(I));

    // recalculate stuff

    return partial_sol;
};

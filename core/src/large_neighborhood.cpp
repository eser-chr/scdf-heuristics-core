#include <cassert>
#include <vector>
#include <iostream>
#include <algorithm>
#include "solvers.hpp"
#include "structures.hpp"

std::vector<int> create_track_route(Instance const &I, int beam_width, std::vector<int> const &requests);

LN::RequestPair LN::find_heaviest_request_in_route(Instance const &I, Encoding const &encoding, int vehicle, int beam_width)
{
    assert(vehicle >= 0 && vehicle < encoding.get_num_vehicles());
    double best_delta = std::numeric_limits<double>::infinity();
    int best_request = -1; // Heaviest
    std::vector<int> rest_requests;

    auto requests = encoding.get_requests_of_route(vehicle);

    if (requests.empty() || requests.size() == 1)
    {
        // Either this track delivers 1 or zero requests. There is no
        // point in removing that so we return an empty removal request.
        // See below in remove requests that this requestpair is ignored.
        // Return sentinel or throw exception
        RequestPair tortn{-1, vehicle, 0.0, {}};
        // tortn.delta = 0.0;
        // tortn.vehicle = vehicle;
        // tortn.request_removed = -1;
        // tortn.rest_requests = {};
        return tortn;
    }

    auto route = create_track_route(I, beam_width, requests);
    if (route.empty())
    {
        std::cerr << "ERROR: create_track_route returned empty route for " << requests.size() << " requests" << std::endl;
        std::abort();
    }

    double original_distance = utils::calc_route_distance(I, route);
    for (auto const &request : requests)
    {
        std::vector<int> new_requests;
        new_requests.reserve(requests.size() - 1);

        for (auto sec_request : requests)
        {
            if (sec_request == request)
                continue;
            new_requests.push_back(sec_request);
        }

        auto new_route = create_track_route(I, beam_width, new_requests);
        double new_distance = utils::calc_route_distance(I, new_route);
        double delta = new_distance - original_distance;
        if (delta < best_delta)
        {
            best_delta = delta;
            best_request = request;
            rest_requests = new_requests;
        }
    }

    assert(best_request != -1);

    RequestPair tortn{best_request, vehicle, best_delta, std::move(rest_requests)};
    // tortn.delta = best_delta;
    // tortn.vehicle = vehicle;
    // tortn.request_removed = best_request;
    // tortn.rest_requests = rest_requests;
    return tortn;
}

Encoding LN::apply_removal(Instance const &I, Encoding const &encoding, std::vector<LN::RequestPair> const &to_be_removed)
{
    auto new_dna = encoding.get_dna_copy();
    if (to_be_removed.size() == 0)
    {
        return {std::move(new_dna)};
    }
    auto num_vehicles = encoding.get_num_vehicles();
    auto num_requests = encoding.get_num_requests();
    assert(num_vehicles == I.nK);
    assert(num_requests == I.n);

    for (auto const &element : to_be_removed)
    {
        auto vehicle = element.vehicle;
        auto request = element.request_removed;
        assert(vehicle < new_dna.size());
        assert(request < new_dna[0].size());


        assert(new_dna[vehicle][request]);
        new_dna[vehicle][request] = false;
    }
    return Encoding{std::move(new_dna)};
}

Encoding LN::apply_addition(Instance const &I, Encoding const &encoding, std::vector<LN::RequestPair> const &to_be_removed)
{
    auto new_dna = encoding.get_dna_copy();

    for (auto const &element : to_be_removed)
    {
        auto vehicle = element.vehicle;
        auto request = element.request_removed;
        assert(!new_dna[vehicle][request]);
        new_dna[vehicle][request] = true;
    }

    return {std::move(new_dna)};
}

Encoding LN::remove_requests(Instance const &I, Encoding const &encoding, int k, int beam_width)
{
    auto const &dna = encoding.get_dna();
    Encoding new_encoding = encoding;
    assert(dna.size() == I.nK);
    int removed = 0;
    while (removed < k)
    {
        std::vector<RequestPair> to_be_removed;
        for (int vehicle = 0; vehicle < I.nK; vehicle++)
        {
            auto request_pair = find_heaviest_request_in_route(I, new_encoding, vehicle, beam_width);
            if (request_pair.request_removed == -1)
                continue;
            to_be_removed.push_back(request_pair);
        }

        assert(new_encoding.get_num_requests() > to_be_removed.size());

        if (k - removed < to_be_removed.size())
        {
            // If I have more valid remove options than I need to remove.
            int how_many_to_remove = k - removed;

            std::partial_sort(to_be_removed.begin(), to_be_removed.begin() + how_many_to_remove,
                              to_be_removed.end(), [](RequestPair const &a, RequestPair const &b)
                              { return a.delta < b.delta; });

            to_be_removed = std::vector<RequestPair>(to_be_removed.begin(), to_be_removed.begin() + how_many_to_remove);
        }

        new_encoding = apply_removal(I, new_encoding, to_be_removed);
        removed += to_be_removed.size();
    }
    return new_encoding;
}

LN::RequestPair LN::find_best_request_to_add(Instance const &I, Encoding const &encoding, int vehicle, int beam_width)
{

    double best_delta = std::numeric_limits<double>::infinity();
    int best_request = -1;           // Heaviest
    std::vector<int> total_requests; // Later to store in the solution

    auto const &dna = encoding.get_dna();
    auto non_delivered_requests = encoding.get_non_delivered_requests();
    auto delivered_requests = encoding.get_requests_of_route(vehicle);
    auto route = create_track_route(I, beam_width, delivered_requests);
    double original_distance = utils::calc_route_distance(I, route);

    for (auto request : non_delivered_requests)
    {
        auto new_delivered_requests = delivered_requests;
        new_delivered_requests.push_back(request);
        auto new_route = create_track_route(I, beam_width, new_delivered_requests);
        double new_distance = utils::calc_route_distance(I, new_route);

        double delta = new_distance - original_distance;
        // assert(delta >= 0);

        if (delta < best_delta)
        {
            best_delta = delta;
            best_request = request;
            total_requests = std::move(new_delivered_requests);
        }
    }
    RequestPair tortn{ best_request, vehicle, best_delta, std::move(total_requests)};
    // tortn.delta = best_delta;
    // tortn.vehicle = vehicle;
    // tortn.request_removed = best_request;
    // tortn.rest_requests = total_requests;

    return tortn;
}

Encoding LN::append_requests(Instance const &I, Encoding const &encoding, int k, int beam_width)
{
    auto const &dna = encoding.get_dna();
    Encoding new_encoding = encoding;

    assert(dna.size() == I.nK);
    int addition_iterations = k / I.nK;
    int additions_in_last_iter = k % I.nK;
    for (size_t iter = 0; iter < addition_iterations; iter++)
    {

        std::vector<RequestPair> to_be_appended;
        for (int vehicle = 0; vehicle < I.nK; vehicle++)
        {
            to_be_appended.push_back(find_best_request_to_add(I, new_encoding, vehicle, beam_width));
        }

        new_encoding = apply_addition(I, new_encoding, to_be_appended);
    }

    // Final for modulo requests
    std::vector<RequestPair> to_be_appended;
    for (int vehicle = 0; vehicle < I.nK; vehicle++)
    {
        to_be_appended.push_back(find_best_request_to_add(I, new_encoding, vehicle, beam_width));
    }
    std::partial_sort(to_be_appended.begin(),
                      to_be_appended.begin() + additions_in_last_iter,
                      to_be_appended.end(),
                      [](RequestPair const &a, RequestPair const &b)
                      { return a.delta < b.delta; });

    // Truncate to only those elements
    auto to_be_appended_final = std::vector<RequestPair>(to_be_appended.begin(), to_be_appended.begin() + additions_in_last_iter);

    new_encoding = apply_addition(I, new_encoding, to_be_appended_final);

    return new_encoding;
}

Solution LN::large_neighborhood(Instance const &I, Solution const &sol, int k, size_t iters, int bw_remove, int bw_append)
{
    Encoding encoding(I, sol);
    auto new_encoding = encoding;

    double best_objective = utils::objective(I, sol);
    Solution best_sol = sol;

    for (size_t iter = 0; iter < iters; ++iter)
    {
        new_encoding = remove_requests(I, new_encoding, k, bw_remove);
        new_encoding = append_requests(I, new_encoding, k, bw_append);

        Solution tmp = new_encoding.to_sol(I);
        double tmp_objective = utils::objective(I, tmp);

        if (tmp_objective < best_objective)
        {
            best_objective = tmp_objective;
            best_sol = tmp;
        }
    }
    assert(best_sol.is_solution_feasible(I));
    return best_sol;
}

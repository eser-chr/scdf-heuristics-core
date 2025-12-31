#include <vector>
#include <algorithm>
#include <limits>
#include <cassert>
#include <iostream>
#include "solvers.hpp"
#include "structures.hpp"

std::vector<int> create_simple_sequential_route(Instance const &I, std::vector<int> const &requests)
{

    std::vector<int> route;
    route.reserve(2 * requests.size());

    // Simple strategy: pick up all, then deliver all
    for (int req : requests)
        route.push_back(1 + req); // pickup

    for (int req : requests)
        route.push_back(1 + I.n + req); // delivery

    return route;
}

std::vector<int> create_track_route(Instance const &I, int beam_width, std::vector<int> const &requests)
{

    if (requests.empty())
    {
        std::cerr << " In create_track_route -> requests were empty !!!" << std::endl;
        return {};
    }

    std::vector<BS::BeamState> beam_states{BS::BeamState{
        0, 0.0, {}, {}, requests}};

    size_t max_steps = 4 * requests.size();

    for (size_t step = 0; step < max_steps; ++step)
    {
        std::vector<BS::BeamState> new_beam;

        for (const auto &st : beam_states)
        {
            if (st.active.empty() && st.remaining.empty())
            {
                new_beam.push_back(st);
                continue;
            }

            // Pickups
            for (int req : st.remaining)
            {
                if (st.cargo + I.demands[req] <= I.C)
                {
                    int p = 1 + req;

                    if (p < 0 || p >= I.dist.size())
                    {
                        std::cerr << "ERROR: Invalid pickup node " << p << " for request " << req << std::endl;
                        continue;
                    }

                    std::vector<int> new_route = st.route;
                    new_route.push_back(p);

                    std::vector<int> new_active = st.active;
                    new_active.push_back(req);

                    std::vector<int> new_remaining;
                    new_remaining.reserve(st.remaining.size() - 1);
                    for (int r : st.remaining)
                    {
                        if (r != req)
                            new_remaining.push_back(r);
                    }

                    int last = st.route.empty() ? 0 : st.route.back();

                    if (last < 0 || last >= I.dist.size() || p >= I.dist[last].size())
                    {
                        std::cerr << "ERROR: Invalid dist access [" << last << "][" << p << "]" << std::endl;
                        continue;
                    }

                    double new_score = st.score + I.dist[last][p];

                    new_beam.push_back(BS::BeamState{
                        st.cargo + I.demands[req],
                        new_score,
                        std::move(new_route),
                        std::move(new_active),
                        std::move(new_remaining)});
                }
            }

            // Deliveries
            for (int req : st.active)
            {
                int d = 1 + I.n + req;

                if (d < 0 || d >= I.dist.size())
                {
                    std::cerr << "ERROR: Invalid delivery node " << d << " for request " << req << std::endl;
                    continue;
                }

                std::vector<int> new_route = st.route;
                new_route.push_back(d);
                int new_cargo = st.cargo - I.demands[req];

                std::vector<int> new_active;
                new_active.reserve(st.active.size() - 1);
                for (int r : st.active)
                {
                    if (r != req)
                        new_active.push_back(r);
                }

                std::vector<int> new_remaining = st.remaining;

                int last = st.route.empty() ? 0 : st.route.back();

                if (last < 0 || last >= I.dist.size() || d >= I.dist[last].size())
                {
                    std::cerr << "ERROR: Invalid dist access [" << last << "][" << d << "]" << std::endl;
                    continue;
                }

                double new_score = st.score + I.dist[last][d];

                new_beam.push_back(BS::BeamState{
                    new_cargo,
                    new_score,
                    std::move(new_route),
                    std::move(new_active),
                    std::move(new_remaining)});
            }
        }

        if (new_beam.empty())
        {
            std::cerr << "WARNING: Beam became empty at step " << step << std::endl;
            break;
        }

        std::sort(new_beam.begin(), new_beam.end(),
                  [](const BS::BeamState &a, const BS::BeamState &b)
                  {
                      return a.score < b.score;
                  });

        if ((int)new_beam.size() > beam_width)
            new_beam.resize(beam_width);

        beam_states = std::move(new_beam);

        bool all_complete = true;
        for (const auto &st : beam_states)
        {
            if (!st.active.empty() || !st.remaining.empty())
            {
                all_complete = false;
                break;
            }
        }

        if (all_complete)
            break;
    }

    double best_score = std::numeric_limits<double>::infinity();
    std::vector<int> best_route;

    for (size_t idx = 0; idx < beam_states.size(); ++idx)
    {
        const auto &st = beam_states[idx];

        if (!st.active.empty() || !st.remaining.empty())
        {
            std::cout << "[SELECT] Skipping incomplete state " << idx << std::endl;
            continue;
        }

        std::vector<int> final_route = st.route;
        for (size_t i = 0; i < final_route.size(); ++i)
        {
            int node = final_route[i];
            if (node < 0 || node >= I.dist.size())
            {
                std::cerr << "ERROR: Invalid node " << node << " at position " << i
                          << " in route of size " << final_route.size() << std::endl;
                return {};
            }
        }

        double d = utils::calc_route_distance(I, final_route);

        if (d < best_score)
        {
            best_score = d;
            best_route = std::move(final_route);
        }
    }

    if (best_route.empty())
    {
        std::cerr << "ERROR: No complete route found for " << requests.size() << " requests" << std::endl;
        return create_simple_sequential_route(I, requests);
    }

    assert(best_route.size() == 2 * requests.size());

    return best_route;
}

Solution BS::beam_search(const Instance &I, double a, int beam_width)
{
    std::vector<double> costs = utils::calc_my_metric(I, a);
    auto perm = numerical::argsort(costs);
    std::vector<int> to_deliver_requests(perm.begin(), perm.begin() + I.gamma);

    // per_track_requests
    gt::Matrix<int> per_track_requests(I.nK);
    for (int t = 0; t < I.nK; ++t)
    {
        for (size_t idx = t; idx < to_deliver_requests.size(); idx += I.nK)
        {
            per_track_requests[t].push_back(to_deliver_requests[idx]);
        }
    }

    gt::Matrix<int> routes;
    routes.reserve(I.nK);

    for (int track = 0; track < I.nK; ++track)
        routes.push_back(std::move(create_track_route(I, beam_width, per_track_requests[track])));

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
    sol.fairness = I.fairness;
    return sol;
}

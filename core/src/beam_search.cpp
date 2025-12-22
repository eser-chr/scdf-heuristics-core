#include <vector>
#include <algorithm>
#include <limits>
#include <cassert>
#include "solvers.hpp"
#include "structures.hpp"

// void flush_deliveries(
//     const Instance &I,
//     std::vector<int> &route,       // will be extended
//     const std::vector<int> &active // remaining pickups (to deliver)
// )
// {
//     std::vector<int> remaining = active;

//     int last = route.empty() ? 0 : route.back();

//     while (!remaining.empty())
//     {
//         int best_r = -1;
//         double best_d = std::numeric_limits<double>::infinity();

//         for (int r : remaining)
//         {
//             int deliver_node = 1 + I.n + r;
//             double d = I.dist[last][deliver_node];
//             if (d < best_d)
//             {
//                 best_d = d;
//                 best_r = r;
//             }
//         }

//         // Append best delivery
//         int deliver_node = 1 + I.n + best_r;
//         route.push_back(deliver_node);

//         // Update last
//         last = deliver_node;

//         // Remove delivered request
//         remaining.erase(
//             std::remove(remaining.begin(), remaining.end(), best_r),
//             remaining.end());
//     }
// }

std::vector<int> create_track_route(Instance const &I, int beam_width, std::vector<int> const &requests)
{
    // std::vector<int> remaining = requests;
    std::vector<BS::BeamState> beam_states{BS::BeamState{
        0.0,
        {},       // route
        0,        // cargo
        {},       // active
        requests // remaining
    }};


    for (size_t step = 0; step < 2*requests.size(); ++step)
    {
        std::vector<BS::BeamState> new_beam;

        for (const auto &st : beam_states) // Append possible paths
        {
            // Remaining to be picked up
            for (int req : st.remaining)
            {
                if (st.cargo + I.demands[req] <= I.C)
                {
                    int p = 1 + req;
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
                    double new_score = st.score + I.dist[last][p];

                    new_beam.push_back(BS::BeamState{
                        new_score,
                        std::move(new_route),
                        st.cargo + I.demands[req],
                        std::move(new_active),
                        std::move(new_remaining)});
                }
            }

            // Active and need to be delivered
            for (int req : st.active)
            {
                int d = 1 + I.n + req;
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

                std::vector<int> new_remaining = st.remaining; // unchanged

                int last = st.route.empty() ? 0 : st.route.back();
                double new_score = st.score + I.dist[last][d];

                new_beam.push_back(BS::BeamState{
                    new_score,
                    std::move(new_route),
                    new_cargo,
                    std::move(new_active),
                    std::move(new_remaining)});
            }
        }

        if (new_beam.empty())
            break;

        std::sort(new_beam.begin(), new_beam.end(),
                  [](const BS::BeamState &a, const BS::BeamState &b)
                  {
                      return a.score < b.score;
                  });

        if ((int)new_beam.size() > beam_width)
            new_beam.resize(beam_width);

        beam_states = std::move(new_beam);
    }

    // finalise each candidate by dropping all active
    double best_score = std::numeric_limits<double>::infinity();
    std::vector<int> best_route;

    for (const auto &st : beam_states)
    {
        std::vector<int> final_route = st.route;
        // flush_deliveries(I, final_route, st.active);
        int d = utils::calc_route_distance(I, final_route);
        if (d < best_score)
        {
            best_score = d;
            best_route = std::move(final_route);
        }
    }

    assert(best_route.size() / 2 == requests.size());
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
    return sol;
}

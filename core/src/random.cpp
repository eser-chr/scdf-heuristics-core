#include <vector>
#include <algorithm>
#include <random>
#include <limits>
#include "structures.hpp"
#include "solvers.hpp"

struct Candidate
{
    int req;      // request id
    int node;     // actual node index in the graph
    double dist;  // distance from last
    bool is_pick; // true = pickup, false = delivery
};

std::vector<Candidate> collect_candidates(
    const Instance &I,
    const std::vector<int> &unpicked,
    const std::vector<int> &active,
    int last,
    int cargo)
{
    std::vector<Candidate> C;
    C.reserve(unpicked.size() + active.size());

    // pickups
    for (int r : unpicked)
    {
        if (cargo + I.demands[r] <= I.C)
        {
            int pn = 1 + r;
            C.push_back({r, pn, I.dist[last][pn], true});
        }
    }

    // deliveries
    for (int r : active)
    {
        int dn = 1 + I.n + r;
        C.push_back({r, dn, I.dist[last][dn], false});
    }

    return C;
}

// ---------------------------------------------------------------
// 2. Pick the best candidate: greedy or exponential random
// ---------------------------------------------------------------
int choose_candidate_index(
    const std::vector<Candidate> &C,
    bool greedy,
    double lambda_exp)
{

    static thread_local std::mt19937 rng(std::random_device{}());
    std::vector<int> idx(C.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(),
              [&](int a, int b)
              { return C[a].dist < C[b].dist; });

    if (greedy)
        return idx[0];

    // exponential softmin
    std::vector<double> weights;
    weights.reserve(idx.size());

    for (int id : idx)
        weights.push_back(std::exp(-lambda_exp * C[id].dist));

    std::discrete_distribution<int> dist(weights.begin(), weights.end());
    int pos = dist(rng); // position in sorted list
    return idx[pos];     // convert to original index
}

// ---------------------------------------------------------------
// 3. Deliver all active when no regular candidate exists
// ---------------------------------------------------------------
Candidate nearest_delivery_fallback(
    const Instance &I,
    const std::vector<int> &active,
    int last)
{
    int best_r = -1;
    int best_node = -1;
    double best_d = 1e18;

    for (int r : active)
    {
        int dn = 1 + I.n + r;
        double d = I.dist[last][dn];
        if (d < best_d)
        {
            best_d = d;
            best_r = r;
            best_node = dn;
        }
    }

    return {best_r, best_node, best_d, false};
}

// ---------------------------------------------------------------
// 4. Main route construction
// ---------------------------------------------------------------
std::vector<int> build_route(
    const Instance &I,
    const std::vector<int> &reqs,
    bool greedy,
    double lambda_exp)
{
    std::vector<int> unpicked = reqs;
    std::vector<int> active;
    std::vector<int> route;

    int cargo = 0;
    int last = 0; // depot

    while (!unpicked.empty() || !active.empty())
    {
        auto C = collect_candidates(I, unpicked, active, last, cargo);

        Candidate choice{-1, -1, 0.0, true};

        if (C.empty())
        {
            choice = nearest_delivery_fallback(I, active, last);
        }
        else
        {
            int ci = choose_candidate_index(C, greedy, lambda_exp);
            choice = C[ci];
        }

        route.push_back(choice.node);

        if (choice.is_pick)
        {
            cargo += I.demands[choice.req];
            active.push_back(choice.req);
            unpicked.erase(std::remove(unpicked.begin(), unpicked.end(), choice.req),
                           unpicked.end());
        }
        else
        {
            cargo -= I.demands[choice.req];
            active.erase(std::remove(active.begin(), active.end(), choice.req),
                         active.end());
        }

        last = choice.node;
    }

    return route;
}



auto select_gamma_requests_random(Instance const &I)
{
    std::vector<double> cost(I.n);

    for (size_t i = 0; i < I.n; ++i)
    {
        // auto pickup_coords = I.coords[1 + i];
        // auto delivery_coords = I.coords[1 + i + I.n];
        cost[i] = I.demands[i] * I.dist[1+i][1+i+I.n];
    }

    auto argsort = numerical::argsort(cost);

    std::vector<int> to_rtn;
    to_rtn.reserve(I.gamma);
    for (int k = 0; k < I.gamma; ++k)
        to_rtn.push_back(static_cast<int>(argsort[k]));

    return to_rtn;
}

Solution RC::construction(
    const Instance &I,
    double lamda)
{
    auto indices_of_requests_to_serve = select_gamma_requests_random(I);

    std::vector<int> assign = clusters::balanced_kmeans(I, indices_of_requests_to_serve, 20, 20);
    gt::Matrix<int> per_track(I.nK); // per track requests-responibilities
    // for (int r = 0; r < I.n; r++)
    //     per_track[assign[r]].push_back(r);

    for (int i = 0; i < (int)indices_of_requests_to_serve.size(); i++)
    {
        int req = indices_of_requests_to_serve[i]; // actual request ID
        int k = assign[i];                         // cluster index
        per_track[k].push_back(req);
    }

    gt::Matrix<int> routes;
    routes.reserve(I.nK);
    for (int k = 0; k < I.nK; k++)
    {
        auto route = build_route(I, per_track[k], false, lamda);
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
    sol.fairness = I.fairness;
    return sol;
}

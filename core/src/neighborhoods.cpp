#include <algorithm>
#include <numeric>
#include <cassert>
#include "neighborhoods.hpp"
#include "structures.hpp"

#include <iostream>
// Creates a vector of info for each request in a route.
std::vector<PickupDeliveryInfo>
pickup_delivery_positions(const Instance &I, const std::vector<int> &route)
{
    struct Pos
    {
        int p_idx = -1;
        int d_idx = -1;
    };
    std::unordered_map<int, Pos> pos;

    for (int idx = 0; idx < (int)route.size(); ++idx)
    {
        int node = route[idx];
        int req = I.request_of_node[node];
        if (req >= 0)
        {
            if (I.load_change[node] > 0)
                pos[req].p_idx = idx;
            else
                pos[req].d_idx = idx;
        }
    }

    std::vector<PickupDeliveryInfo> out;
    out.reserve(pos.size());
    for (const auto &kv : pos)
    {
        int req = kv.first;
        int p_idx = kv.second.p_idx;
        int d_idx = kv.second.d_idx;
        if (p_idx >= 0 && d_idx >= 0)
        {
            int pickup_node = 2 * req + 1;
            int delivery_node = 2 * req + 2;
            out.push_back(PickupDeliveryInfo{p_idx, d_idx, req, pickup_node, delivery_node});
        }
    }
    return out;
}

std::vector<GenericMove> IntraRouteNeighborhood::generate() const
{
    std::vector<GenericMove> to_rtn;
    for (int r = 0; r < (int)sol.routes.size(); ++r)
    {
        const auto &route = sol.routes[r];
        int m = (int)route.size();
        for (int k = 0; k < m; ++k)
        {
            for (int l = k + 1; l < m; ++l) // No for l<k+1 because of double counting.
            {
                to_rtn.push_back(GenericMove{type, {r, k, l}});
            }
        }
    }
    return to_rtn;
}

std::optional<GenericMove>
IntraRouteNeighborhood::generate_random(std::mt19937 &rng) const
{
    if (sol.routes.empty())
        return std::nullopt;

    std::uniform_int_distribution<int> route_dist(0, sol.routes.size() - 1);

    for (int tries = 0; tries < this->MAX_TRIES_RANDOM; tries++)
    {
        int r = route_dist(rng);
        int m = sol.routes[r].size();

        // Must have at least 2 positions
        if (m < 2)
            continue;

        std::uniform_int_distribution<int> idx_dist(0, m - 1);

        // 2) Pick random indices k < l
        int k = idx_dist(rng);
        int l = idx_dist(rng);
        if (k == l)
            continue;
        if (k > l)
            std::swap(k, l);

        GenericMove mvs{type, {r, k, l}};
        if (is_valid(mvs))
            return mvs;
    }

    return std::nullopt;
}

bool IntraRouteNeighborhood::is_valid(const GenericMove &mov) const
{
    assert(mov.type == type);

    int r = mov.data[0];
    int k = mov.data[1];
    int l = mov.data[2];

    auto route = sol.routes[r];
    std::swap(route[k], route[l]);

    auto cargo = utils::calc_route_cargo(I, route);
    for (int c : cargo)
    {
        if (c >= I.C)
            return false;
    }
    return true;
}

double IntraRouteNeighborhood::calc_delta(const GenericMove &mov) const
{
    assert(mov.type == type);
    int r = mov.data[0];
    int k = mov.data[1];
    int l = mov.data[2];

    // if (k == l)
    //     return 0.0;
    if (k >= l)
    {
        std::cerr << "calc delta k>=l failed.assertion";
        std::abort();
    }

    const auto &route = sol.routes[r];
    const auto &dist = I.dist;

    int x = route[k];
    int y = route[l];

    int A = (k > 0) ? route[k - 1] : 0;
    int B = (k + 1 < (int)route.size()) ? route[k + 1] : 0;
    int C = (l > 0) ? route[l - 1] : 0;
    int D = (l + 1 < (int)route.size()) ? route[l + 1] : 0;

    int delta_d;

    if (l == k + 1)
    {
        delta_d =
            dist[A][y] + dist[y][x] + dist[x][D] -
            (dist[A][x] + dist[x][y] + dist[y][D]);
    }
    else
    {
        delta_d =
            dist[A][y] + dist[y][B] + dist[C][x] + dist[x][D] -
            (dist[A][x] + dist[x][B] + dist[C][y] + dist[y][D]);
    }

    double d_old = utils::calc_route_distance(I, route);
    double d_new = d_old + delta_d;

    double S_old = sol.total_distance;
    double Q_old = sol.sum_of_squares;

    double S_new = S_old - d_old + d_new;
    double Q_new = Q_old - d_old * d_old + d_new * d_new;

    double J_old = (S_old * S_old) / (I.nK * Q_old);
    double J_new = (S_new * S_new) / (I.nK * Q_new);

    return delta_d + I.rho * (J_old - J_new);
}

Solution IntraRouteNeighborhood::apply(const GenericMove &mov) const
{
    assert(mov.type == type);

    int r = mov.data[0];
    int k = mov.data[1];
    int l = mov.data[2];

    Solution new_sol = sol;
    std::swap(new_sol.routes[r][k], new_sol.routes[r][l]);

    auto all_distances = utils::all_route_distances(I, new_sol);
    std::vector<double> sq_distances;
    sq_distances.resize(all_distances.size());
    std::transform(all_distances.begin(), all_distances.end(), sq_distances.begin(), [](auto val)
                   { return val * val; });
    new_sol.total_distance = std::accumulate(all_distances.begin(), all_distances.end(), 0.0);

    new_sol.sum_of_squares = std::accumulate(sq_distances.begin(), sq_distances.end(), 0.0);
    return new_sol;
}

// =====================================================================
// 2. PairRelocateNeighborhood
// =====================================================================

// std::vector<GenericMove> PairRelocateNeighborhood::generate() const
// {
//     std::vector<GenericMove> to_rtn;
//     int R = (int)sol.routes.size();

//     for (int r_from = 0; r_from < R; ++r_from)
//     {
//         const auto &routeA = sol.routes[r_from];
//         auto infos = pickup_delivery_positions(I, routeA);

//         for (const auto &info : infos)
//         {
//             for (int r_to = 0; r_to < R; ++r_to)
//             {
//                 if (r_to == r_from)
//                     continue;

//                 int lenB = (int)sol.routes[r_to].size();
//                 for (int p_new = 0; p_new <= lenB; ++p_new)
//                 {
//                     for (int d_new = p_new + 1; d_new <= lenB + 1; ++d_new)
//                     {
//                         to_rtn.push_back(GenericMove{
//                             type,
//                             {r_from, info.p_idx, info.d_idx,
//                              r_to, p_new, d_new,
//                              info.req, info.pickup_node, info.delivery_node}});
//                     }
//                 }
//             }
//         }
//     }
//     return to_rtn;
// }

// std::optional<GenericMove>
// PairRelocateNeighborhood::generate_random(std::mt19937 &rng) const
// {
//     int R = (int)sol.routes.size();
//     if (R < 2)
//         return std::nullopt;

//     // 1) Pick r_from randomly
//     std::uniform_int_distribution<int> route_dist(0, R - 1);

//     for (int tries = 0; tries < this->MAX_TRIES_RANDOM; tries++)
//     {
//         int r_from = route_dist(rng);
//         const auto &routeA = sol.routes[r_from];

//         // Need at least one pickup-delivery pair
//         auto infos = pickup_delivery_positions(I, routeA);
//         if (infos.empty())
//             continue;

//         // 2) Choose a random pair (pickup_idx, delivery_idx)
//         std::uniform_int_distribution<int> info_dist(0, (int)infos.size() - 1);
//         const auto &info = infos[info_dist(rng)];

//         // 3) Pick r_to randomly, different from r_from
//         int r_to = r_from;
//         if (R > 1)
//         {
//             do
//             {
//                 r_to = route_dist(rng);
//             } while (r_to == r_from);
//         }

//         const auto &routeB = sol.routes[r_to];
//         int lenB = (int)routeB.size();
//         if (lenB < 0)
//             continue;

//         // 4) Random insertion positions for pickup and delivery
//         std::uniform_int_distribution<int> pos_dist(0, lenB + 1); // inclusive

//         int p_new = pos_dist(rng);
//         std::uniform_int_distribution<int> d_dist(p_new + 1, lenB + 2);
//         int d_new = d_dist(rng);

//         // 5) Build the candidate move
//         GenericMove m{
//             2,
//             {r_from, info.p_idx, info.d_idx,
//              r_to, p_new, d_new,
//              info.req, info.pickup_node, info.delivery_node}};

//         // 6) Validate it
//         if (is_valid(m))
//             return m;
//     }

//     return std::nullopt;
// }

// bool PairRelocateNeighborhood::is_valid(const GenericMove &mov) const
// {
//     int r_from = mov.data[0];
//     int p_old = mov.data[1];
//     int d_old = mov.data[2];
//     int r_to = mov.data[3];
//     int p_new = mov.data[4];
//     int d_new = mov.data[5];
//     int pnode = mov.data[7];
//     int dnode = mov.data[8];

//     auto routeA = sol.routes[r_from];
//     auto routeB = sol.routes[r_to];

//     // FIX 1: erase in correct order depending on indices
//     if (p_old < d_old)
//     {
//         routeA.erase(routeA.begin() + d_old);
//         routeA.erase(routeA.begin() + p_old);
//     }
//     else
//     {
//         routeA.erase(routeA.begin() + p_old);
//         routeA.erase(routeA.begin() + d_old);
//     }

//     if (p_new < 0 || p_new > (int)routeB.size())
//         return false;
//     if (d_new < 0 || d_new > (int)routeB.size() + 1)
//         return false;

//     routeB.insert(routeB.begin() + p_new, pnode);
//     // after inserting pickup, delivery index is in [p_new+1, size]
//     if (d_new > (int)routeB.size())
//         return false;
//     routeB.insert(routeB.begin() + d_new, dnode);

//     for (int c : utils::calc_route_cargo(I, routeA))
//         if (c > I.C)
//             return false;
//     for (int c : utils::calc_route_cargo(I, routeB))
//         if (c > I.C)
//             return false;

//     return true;
// }

// Solution PairRelocateNeighborhood::apply(const GenericMove &mov) const
// {
//     int r_from = mov.data[0];
//     int p_old = mov.data[1];
//     int d_old = mov.data[2];
//     int r_to = mov.data[3];
//     int p_new = mov.data[4];
//     int d_new = mov.data[5];
//     int pnode = mov.data[7];
//     int dnode = mov.data[8];

//     Solution new_sol = sol;

//     auto &routeA = new_sol.routes[r_from];
//     auto &routeB = new_sol.routes[r_to];

//     if (p_old < d_old)
//     {
//         routeA.erase(routeA.begin() + d_old);
//         routeA.erase(routeA.begin() + p_old);
//     }
//     else
//     {
//         routeA.erase(routeA.begin() + p_old);
//         routeA.erase(routeA.begin() + d_old);
//     }

//     routeB.insert(routeB.begin() + p_new, pnode);
//     routeB.insert(routeB.begin() + d_new, dnode);

//     auto all_distances = utils::all_route_distances(I, new_sol);
//     std::vector<double> sq_distances;
//     sq_distances.resize(all_distances.size());
//     std::transform(all_distances.begin(), all_distances.end(), sq_distances.begin(), [](auto val)
//                    { return val * val; });
//     new_sol.total_distance = std::accumulate(all_distances.begin(), all_distances.end(), 0.0);

//     new_sol.sum_of_squares = std::accumulate(sq_distances.begin(), sq_distances.end(), 0.0);

//     return new_sol;
// }
// double PairRelocateNeighborhood::calc_delta(const GenericMove &mov) const
// {
//     int r_from = mov.data[0];
//     int p_old  = mov.data[1];
//     int d_old  = mov.data[2];
//     int r_to   = mov.data[3];
//     int p_new  = mov.data[4];
//     int d_new  = mov.data[5];
//     int pnode  = mov.data[7];
//     int dnode  = mov.data[8];

//     const auto &routeA_orig = sol.routes[r_from];
//     const auto &routeB_orig = sol.routes[r_to];

//     double dA_old = utils::route_distance(I, routeA_orig);
//     double dB_old = utils::route_distance(I, routeB_orig);

//     auto routeA = routeA_orig;
//     auto routeB = routeB_orig;

//     if (p_old < d_old)
//     {
//         routeA.erase(routeA.begin() + d_old);
//         routeA.erase(routeA.begin() + p_old);
//     }
//     else
//     {
//         routeA.erase(routeA.begin() + p_old);
//         routeA.erase(routeA.begin() + d_old);
//     }

//     routeB.insert(routeB.begin() + p_new, pnode);
//     if (d_new > (int)routeB.size())
//         return std::numeric_limits<double>::infinity();
//     routeB.insert(routeB.begin() + d_new, dnode);

//     double dA_new = utils::route_distance(I, routeA);
//     double dB_new = utils::route_distance(I, routeB);

//     double S_old = sol.total_distance;
//     double Q_old = sol.sum_of_squares;

//     double S_new = S_old - dA_old - dB_old + dA_new + dB_new;
//     double Q_new = Q_old
//                  - dA_old * dA_old + dA_new * dA_new
//                  - dB_old * dB_old + dB_new * dB_new;

//     double J_old = (S_old * S_old) / (I.nK * Q_old);
//     double J_new = (S_new * S_new) / (I.nK * Q_new);

//     double delta_d = (dA_new - dA_old) + (dB_new - dB_old);
//     return delta_d + I.rho * (J_old - J_new);
// }

std::vector<int> get_request_indices(Instance const &I, std::vector<int> const &route)
{
    std::vector<int> to_rtn;
    to_rtn.reserve(route.size() / 2);

    for (int node : route)
    {
        if (node < I.n + 1)
        {
            to_rtn.push_back(node-1);
        }
    }
    return to_rtn;
}

std::vector<GenericMove> RequestMove::generate() const
{
    std::vector<GenericMove> to_rtn;
    for (int from = 0; I.nK; from++)
    {
        auto request_indices = get_request_indices(I, sol.routes[from]);
        for (int to = from + 1; I.nK; to++)
        {
            for (auto request : request_indices)
                to_rtn.push_back(GenericMove{type, {from, to, request}});
        }
    }
    return to_rtn;
}

std::optional<GenericMove> RequestMove::generate_random(std::mt19937 &rng) const
{
    assert(sol.routes.size()==I.nK);
    std::uniform_int_distribution<int> routes_dist(0, I.nK - 1);
    for (size_t i = 0; i < this->MAX_TRIES_RANDOM; ++i)
    {
        int from = routes_dist(rng);
        int to = routes_dist(rng);

        if (sol.routes[from].size() < 2)
            continue;
        if (from == to)
            continue;

        auto request_indices = get_request_indices(I, sol.routes[from]);
        std::uniform_int_distribution<int> request_dist(0, request_indices.size() - 1);
        int req = request_indices[request_dist(rng)];

        return GenericMove{type, {from, to, req}};
    }
    return std::nullopt;
}

bool RequestMove::is_valid(GenericMove const &mov) const
{
    assert(mov.type == type);

    return true;
}

double RequestMove::calc_delta(GenericMove const &move) const
{
    assert(move.type == type);
    int from = move.data[0];
    int to = move.data[1];
    int request = move.data[2];

    double d_old_from = sol.routes_distances[from];
    double d_old_to = sol.routes_distances[to];

    // auto new_from = sol.routes[from];
    // auto new_to = sol.routes[to];

    // Remove old delivery and pick up.
    //  Add new delivery and pick up.
    std::vector<int> new_from, new_to;
    new_from.reserve(sol.routes[from].size() - 2);
    new_to.reserve(sol.routes[to].size() + 2);

    for (auto node : sol.routes[from])
    {
        if (node != request + 1 || node != request + I.n + 1)
            new_from.push_back(node);
    }

    for (auto node : sol.routes[to])
        new_to.push_back(node);

    new_to.push_back(request + 1);
    new_to.push_back(request + 1 + I.n);

    double d_new_from = utils::calc_route_distance(I, new_from);
    double d_new_to = utils::calc_route_distance(I, new_to);

    //         double d_old = utils::route_distance(I, route);
    //     double d_new = d_old + delta_d;

    //     double S_old = sol.total_distance;
    //     double Q_old = sol.sum_of_squares;

    //     double S_new = S_old - d_old + d_new;
    //     double Q_new = Q_old - d_old * d_old + d_new * d_new;

    //     double J_old = (S_old * S_old) / (I.nK * Q_old);
    //     double J_new = (S_new * S_new) / (I.nK * Q_new);

    //     return delta_d + I.rho * (J_old - J_new);
    // }

    double delta_d = d_new_from - d_old_from + d_new_to - d_old_to;

    // Needs to implement the delta in the fairness function.
    return delta_d;

}

Solution RequestMove::apply(GenericMove const &move) const
{
    assert(move.type == type);
    int from = move.data[0];
    int to = move.data[1];
    int request = move.data[2];

    std::vector<int> new_from, new_to;
    new_from.reserve(sol.routes[from].size() - 2);
    new_to.reserve(sol.routes[to].size() + 2);

    for (auto node : sol.routes[from])
    {
        if (node != request + 1 || node != request + I.n + 1)
            new_from.push_back(node);
    }

    for (auto node : sol.routes[to])
        new_to.push_back(node);

    new_to.push_back(request + 1);
    new_to.push_back(request + 1 + I.n);

    Solution new_sol = sol;
    new_sol.routes[from] = std::move(new_from);
    new_sol.routes[to] = std::move(new_to);

    double d_new_from = utils::calc_route_distance(I, new_from);
    double d_new_to = utils::calc_route_distance(I, new_to);

    

    new_sol.routes_distances[from] = d_new_from;
    new_sol.routes_distances[to] = d_new_to;
    new_sol.total_distance = std::accumulate(new_sol.routes_distances.begin(), new_sol.routes_distances.end(), 0.0);
    new_sol.sum_of_squares = std::accumulate(new_sol.routes_distances.begin(), new_sol.routes_distances.end(), 0.0, [](double sum, double x)
                                             { return sum + x * x; });

    return sol ;
}

// =====================================================================
// 3. TwoOptNeighborhood
// =====================================================================

std::vector<GenericMove> TwoOptNeighborhood::generate() const
{
    std::vector<GenericMove> to_rtn;
    for (int r = 0; r < (int)sol.routes.size(); ++r)
    {
        const auto &route = sol.routes[r];
        int m = (int)route.size();
        for (int i = 0; i < m - 2; ++i)
        {
            for (int j = i + 2; j < m; ++j)
            {
                to_rtn.push_back(GenericMove{3, {r, i, j}});
            }
        }
    }
    return to_rtn;
}

std::optional<GenericMove>
TwoOptNeighborhood::generate_random(std::mt19937 &rng) const
{
    if (sol.routes.empty())
        return std::nullopt;

    std::uniform_int_distribution<int> route_dist(0, sol.routes.size() - 1);
    int rid = route_dist(rng);

    const auto &r = sol.routes[rid];
    if (r.size() < 4)
        return std::nullopt;

    std::uniform_int_distribution<int> idx_dist(0, r.size() - 1);

    for (int t = 0; t < this->MAX_TRIES_RANDOM; t++)
    { // few random attempts
        int i = idx_dist(rng);
        int j = idx_dist(rng);
        if (i == j)
            continue;
        if (i > j)
            std::swap(i, j);

        GenericMove m{3, {rid, i, j}};

        if (is_valid(m))
            return m;
    }

    return std::nullopt;
}

bool TwoOptNeighborhood::is_valid(const GenericMove &mov) const
{
    int r = mov.data[0];
    int i = mov.data[1];
    int j = mov.data[2];

    auto route = sol.routes[r];
    std::reverse(route.begin() + i, route.begin() + j + 1);

    for (int c : utils::calc_route_cargo(I, route))
        if (c > I.C)
            return false;

    int nreq = I.n;
    std::vector<int> p(nreq, -1), d(nreq, -1);

    for (int idx = 0; idx < (int)route.size(); ++idx)
    {
        int node = route[idx];
        int req = I.request_of_node[node];
        if (req < 0)
            continue;

        if (I.load_change[node] > 0)
        {
            if (p[req] == -1)
                p[req] = idx;
        }
        else
        {
            if (d[req] == -1)
                d[req] = idx;
        }
    }

    for (int req = 0; req < nreq; ++req)
    {
        if (p[req] != -1 && d[req] != -1 && p[req] > d[req])
            return false;
    }

    return true;
}

double TwoOptNeighborhood::calc_delta(const GenericMove &mov) const
{
    int r = mov.data[0];
    int i = mov.data[1];
    int j = mov.data[2];

    const auto &route = sol.routes[r];
    const auto &dist = I.dist;

    int A = (i > 0) ? route[i - 1] : 0;
    int x = route[i];
    int y = route[j];
    int B = (j + 1 < (int)route.size()) ? route[j + 1] : 0;

    int removed = dist[A][x] + dist[y][B];
    int added = dist[A][y] + dist[x][B];

    double delta_d = added - removed;

    double d_old = utils::calc_route_distance(I, route);
    double d_new = d_old + delta_d;

    double S_old = sol.total_distance;
    double Q_old = sol.sum_of_squares;

    double S_new = S_old - d_old + d_new;
    double Q_new = Q_old - d_old * d_old + d_new * d_new;

    double J_old = (S_old * S_old) / (I.nK * Q_old);
    double J_new = (S_new * S_new) / (I.nK * Q_new);

    return delta_d + I.rho * (J_old - J_new);
}

Solution TwoOptNeighborhood::apply(const GenericMove &mov) const
{
    int r = mov.data[0];
    int i = mov.data[1];
    int j = mov.data[2];

    Solution new_sol = sol;
    std::reverse(new_sol.routes[r].begin() + i, new_sol.routes[r].begin() + j + 1);

    auto all_distances = utils::all_route_distances(I, new_sol);
    std::vector<double> sq_distances;
    sq_distances.resize(all_distances.size());
    std::transform(all_distances.begin(), all_distances.end(), sq_distances.begin(), [](auto val)
                   { return val * val; });
    new_sol.total_distance = std::accumulate(all_distances.begin(), all_distances.end(), 0.0);

    new_sol.sum_of_squares = std::accumulate(sq_distances.begin(), sq_distances.end(), 0.0);

    return new_sol;
}

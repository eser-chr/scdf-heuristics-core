#include <algorithm>
#include <iostream>
#include <sstream>
#include <cassert>
#include <unordered_set>
#include "structures.hpp"
namespace utils
{
    // Returns the exscan of cargo
    std::vector<int> calc_route_cargo(
        const Instance &inst,
        const std::vector<int> &route)
    {
        std::vector<int> out(route.size());
        int load = 0;

        for (size_t t = 0; t < route.size(); t++)
        {
            int node = route[t];
            load += inst.load_change[node];
            out[t] = load;
        }
        return out;
    }

    double calc_route_distance(const Instance &inst, const std::vector<int> &route)
    {
        if (route.empty())
            return 0;

        double d = inst.dist[0][route[0]];
        for (size_t i = 0; i + 1 < route.size(); i++)
            d += inst.dist[route[i]][route[i + 1]];
        d += inst.dist[route.back()][0];

        return d;
    }

    double calc_route_distance(Instance const &I, Solution const &sol, int route_idx)
    {
        return calc_route_distance(I, sol.routes[route_idx]);
    }

    std::vector<double> all_route_distances(
        const Instance &inst,
        const Solution &sol)
    {
        std::vector<double> res(sol.routes.size());
        for (size_t i = 0; i < sol.routes.size(); i++)
            res[i] = (double)calc_route_distance(inst, sol.routes[i]);
        return res;
    }

    double jain_fairness(Instance const &I, const std::vector<double> &dists)
    {
        if (dists.empty())
            throw std::runtime_error("dist has length 0!");

        double sum = 0.0;
        double sq_sum = 0.0;

        for (double x : dists)
        {
            sum += x;
            sq_sum += x * x;
        }

        double num = sum * sum;
        double den = I.nK * sq_sum;

        if (den == 0.0)
        {
            std::stringstream ss;
            ss << "Division by zero in Jain fairness. "
               << "len=" << dists.size() << ", sq_sum=" << sq_sum;
            throw std::runtime_error(ss.str());
        }

        return num / den;
    }

    double max_min_fairness(Instance const &I, std::vector<double> const &dists)
    {
        if (dists.empty())
            throw std::runtime_error("dist has length 0!");

        auto [min_it, max_it] = std::minmax_element(dists.begin(), dists.end());
        double min = *min_it;
        double max = *max_it;

        return min / max;
    }

    double gini_cefficient_nominator(Instance const &I, std::vector<double> const &dists)
    {
        double nominator = 0.0;

        for (size_t i = 0; i < dists.size(); ++i)
        {
            for (size_t j = i + 1; j < dists.size(); ++j)
            { // Avoid double & self counting
                nominator += std::abs(dists[i] - dists[j]);
            }
        }
        return nominator;
    }

    double gini_cefficient(Instance const &I, std::vector<double> const &dists)
    {
        double denominator = 0.0;
        double nominator = 0.0;

        for (size_t i = 0; i < dists.size(); ++i)
        {
            auto const &d_i = dists[i];
            denominator += d_i;

            for (size_t j = i + 1; j < dists.size(); ++j)
            { // Avoid double & self counting
                nominator += std::abs(d_i - dists[j]);
            }
        }

        // Because of not double counting we do not need the 2 in the denominator
        return 1 - (nominator / denominator);
    }

    bool is_route_feasible(const Instance &inst,
                           const std::vector<int> &route)
    {
        int load = 0;
        std::unordered_set<int> picked;
        std::unordered_set<int> dropped;

        for (int node : route)
        {

            int req = inst.request_of_node[node];
            if (req < 0)
                return false; // invalid node inside route

            load += inst.load_change[node];
            if (load > inst.C || load < 0)
                return false; // capacity violation

            if (inst.load_change[node] > 0)
                picked.insert(req);
            else
            {
                if (!picked.count(req))
                    return false; // drop before pickup
                dropped.insert(req);
            }
        }

        // must have at least gamma drops
        return (int)dropped.size() >= inst.gamma;
    }

    // double objective(const Instance &inst, const Solution &sol)
    // {
    //     auto dists = all_route_distances(inst, sol);

    //     double sum_dist = std::accumulate(dists.begin(), dists.end(), 0.0);
    //     double fairness = jain_fairness(inst, dists);

    //     return sum_dist + inst.rho * (1.0 - fairness);
    // }
    double objective(Instance const &I, Solution const &sol)
    {
        assert(I.fairness == sol.fairness);
        auto dists = all_route_distances(I, sol);

        double sum_dist = std::accumulate(dists.begin(), dists.end(), 0.0);
        double fairness;
        if (sol.fairness == "jain")
            fairness = jain_fairness(I, dists);
        if (sol.fairness == "gini")
            fairness = gini_cefficient(I, dists);
        if (sol.fairness == "maxmin")
            fairness = max_min_fairness(I, dists);

        // double fairness = fairness_func(I, dists);

        return sum_dist + I.rho * (1.0 - fairness);
    }

    std::vector<double> calc_my_metric(const Instance &I, double a)
    {
        int n = I.n;

        // Solo trip from depot pickup delivery depot.
        std::vector<double> solo(n);
        for (int req = 0; req < n; ++req)
        {
            int p = 1 + req;
            int d = 1 + n + req;
            int s = I.dist[0][p] + I.dist[p][d] + I.dist[d][0];
            solo[req] = static_cast<double>(s);
        }

        double max_dist = 0.0;
        for (double v : solo)
            max_dist = std::max(max_dist, v);

        int max_dem_int = 0;
        for (int c : I.demands)
            max_dem_int = std::max(max_dem_int, c);

        if (max_dist == 0.0)
            max_dist = 1.0;
        if (max_dem_int == 0)
            max_dem_int = 1;

        double max_dem = static_cast<double>(max_dem_int);

        std::vector<double> costs(n);
        for (int i = 0; i < n; ++i)
        {
            double dist_norm = solo[i] / max_dist;
            double dem_norm = static_cast<double>(I.demands[i]) / max_dem;
            costs[i] = a * dist_norm + (1.0 - a) * dem_norm;
        }
        return costs;
    }

} // namespace utils

namespace numerical
{
    template <typename T>
    std::vector<int> argsort(const std::vector<T> &org)
    {
        std::vector<int> perm(org.size());
        std::iota(perm.begin(), perm.end(), 0);
        std::sort(perm.begin(), perm.end(),
                  [&](int i, int j)
                  { return org[i] < org[j]; });
        return perm;
    }

    template <typename T>
    T select_uniformly(const std::vector<T> &org, std::mt19937 &rng)
    {
        std::uniform_int_distribution<int> idx(0, org.size());
        size_t i = idx(rng);
        return org[i];
    }

    double calc_distance_between_nodes(gt::Coords const &p1, gt::Coords const &p2)
    {
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        return std::sqrt(dx * dx + dy * dy);
    }

};
template std::vector<int> numerical::argsort<double>(const std::vector<double> &);
template std::vector<int> numerical::argsort<int>(const std::vector<int> &);
template double numerical::select_uniformly(const std::vector<double> &org, std::mt19937 &rng);
template int numerical::select_uniformly(const std::vector<int> &org, std::mt19937 &rng);
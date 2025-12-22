#include <vector>
#include <algorithm>
#include <random>
#include <functional>
#include "solvers.hpp"
#include "structures.hpp"


Solution GRASP::randomized_constructor_simple(
    const Instance &I,
    double a,
    double alpha)
{
    const int n  = I.n;
    const int nK = I.nK;
    const int max_tries = 100;

    // heuristic costs
    std::vector<double> costs = utils::calc_my_metric(I, a);

    // argsort(costs)
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(),
              [&](int i, int j) { return costs[i] < costs[j]; });

    Solution sol;
    sol.routes.assign(nK, std::vector<int>{});

    int served = 0;
    std::vector<bool> used(n, false);

    static thread_local std::mt19937 rng(std::random_device{}());

    while (served < I.gamma)
    {
        // remaining requests (not yet used)
        std::vector<int> remaining;
        remaining.reserve(n);
        for (int r : perm)
            if (!used[r])
                remaining.push_back(r);

        if (remaining.empty())
            break;

        int k = std::max(1, int(alpha * remaining.size()));
        if (k > (int)remaining.size())
            k = (int)remaining.size();

        std::uniform_int_distribution<int> pick_rcl(0, k - 1);
        int req = remaining[pick_rcl(rng)];
        used[req] = true;

        int pickup = 1 + req;
        int drop   = 1 + I.n + req;
        int dem    = I.demands[req];

        bool inserted = false;

        for (int t = 0; t < max_tries; ++t)
        {
            std::uniform_int_distribution<int> pick_route(0, nK - 1);
            int vk = pick_route(rng);
            auto &route = sol.routes[vk];
            int m = (int)route.size();

            // handle empty route explicitly
            if (m == 0)
            {
                if (dem <= I.C) {
                    route.push_back(pickup);
                    route.push_back(drop);
                    inserted = true;
                    served++;
                    break;
                }
                // cannot serve on empty vehicle (capacity too small), try another vehicle/attempt
                continue;
            }

            // non-empty route: choose ip, jp safely
            // ip in [0, m-1], jp in [ip+1, m]
            std::uniform_int_distribution<int> pick_ip(0, m - 1);
            int ip = pick_ip(rng);

            std::uniform_int_distribution<int> pick_jp(ip + 1, m);
            int jp = pick_jp(rng);

            // build candidate route:
            std::vector<int> new_r;
            new_r.reserve(m + 2);

            // [0, ip)
            new_r.insert(new_r.end(), route.begin(), route.begin() + ip);
            // pickup
            new_r.push_back(pickup);
            // [ip, jp)
            new_r.insert(new_r.end(), route.begin() + ip, route.begin() + jp);
            // drop
            new_r.push_back(drop);
            // [jp, m)
            new_r.insert(new_r.end(), route.begin() + jp, route.end());

            auto cargo = utils::calc_route_cargo(I, new_r);
            bool ok = true;
            for (int c : cargo)
            {
                if (c < 0 || c > I.C)
                {
                    ok = false;
                    break;
                }
            }

            if (ok)
            {
                route = std::move(new_r);
                inserted = true;
                served++;
                break;
            }
        }

        // failed to insert this request in max_tries attempts → skip it
        if (!inserted)
            continue;
    }

    return sol;
}


Solution GRASP::grasp(
    const Instance &I,
    std::function<Solution(const Instance &)> randomized_constructor,
    const Neighborhood::NeighborhoodFactories &neighborhoods,
    StepFunction::Func step_function,
    StoppingCriterion &stopping_outer,
    StoppingCriterion &stopping_local,
    int *iteration_ptr)
{
    Solution best_sol; // final result
    double best_f = std::numeric_limits<double>::infinity();

    int step = 0;

    // reset local-search stop criterion every restart
    stopping_local.reset();
    stopping_outer.reset();

    while (!stopping_outer(step, best_f))
    {
        // === construct a new initial solution ===
        Solution sol0 = randomized_constructor(I);

        stopping_local.reset();
        Solution sol1 = LS::local_search(
            I,
            sol0,
            neighborhoods[step % neighborhoods.size()], // rotate neighborhood
            step_function,
            stopping_local);

        double f1 = utils::objective(I, sol1);

        if (f1 < best_f)
        {
            best_f = f1;
            best_sol = sol1;
        }

        step++;
    }
    if(iteration_ptr!=nullptr){
        *iteration_ptr = step;
    }

    return best_sol;
}

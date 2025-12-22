#include <chrono>
#include <iostream>
#include <optional>
#include "neighborhoods.hpp"
#include "solvers.hpp"
#include "structures.hpp"

#include "stopping_criteria.hpp"
#include "step_function.hpp"

bool do_i_accept(double delta, double T, std::mt19937 &rng)
{
    static std::uniform_real_distribution<double> uni(0.0, 1.0);

    if (delta < 0.0)
    {
        return true;
    }
    else
    {
        double p = std::exp(-delta / T);
        if (uni(rng) < p)
            return true;
    }
    return false;
}

Solution SA::simulated_annealing(
    const Instance &I,
    const Solution &initial_sol,
    const Neighborhood::NeighborhoodFactories &neighborhood_factories,
    double T_start,
    double T_end,
    double cooling,
    StepFunction::Func step_function, // Random move preferebaly
    StoppingCriterion &stopping_criterion,
    int *iteration_ptr )
{

    static std::mt19937 rng(std::random_device{}());
    Solution sol = initial_sol;
    Solution best_sol = sol;
    double f = utils::objective(I, sol);
    double best_f = f;
    double T = T_start;
    size_t i = 0;

    while (!stopping_criterion(i, best_f))
    {
        T = std::max(T, T_end);
        std::uniform_int_distribution<int> neigh_dist(0, (int)neighborhood_factories.size() - 1);
        int idx_neigh = neigh_dist(rng);
        auto neigh = neighborhood_factories[idx_neigh](I, sol);
        auto const &mov = step_function(*(neigh), rng);

        if (!mov.has_value())
        {
            break;
        }
        auto actual_move = mov.value();
        double delta = neigh->calc_delta(actual_move);
        bool accept = do_i_accept(delta, T, rng);

        if (accept)
        {
            sol = neigh->apply(actual_move);
            f += delta;

            if (f < best_f)
            {
                best_f = f;
                best_sol = sol;
            }
        }

        T *= cooling;
        i++;
    }
    if(iteration_ptr!= nullptr){
        *iteration_ptr = i;
    }

    return best_sol;
}

// Solution SA::simulated_annealing(
//     const Instance &I,
//     const Solution &initial_sol,
//     const Neighborhood::NeighborhoodFactories &neighborhood_factories,
//     double T_start,
//     double T_end,
//     double cooling,
//     StepFunction::Func step_function, // random move
//     StoppingCriterion &stopping_criterion)
// {
//     using clk = std::chrono::high_resolution_clock;

//     static std::mt19937 rng(std::random_device{}());
//     Solution sol = initial_sol;
//     Solution best_sol = sol;

//     double f = utils::objective(I, sol);
//     double best_f = f;
//     double T = T_start;
//     size_t iter = 0;

//     // Timing accumulators
//     double time_neigh = 0.0;
//     double time_step = 0.0;
//     double time_accept = 0.0;
//     auto t_start_total = clk::now();

//     while (!stopping_criterion(iter, best_f))
//     {
//         T = std::max(T, T_end);

//         // Pick neighborhood
//         std::uniform_int_distribution<int> neigh_dist(
//             0, (int)neighborhood_factories.size() - 1);
//         int idx_neigh = neigh_dist(rng);

//         // NEIGHBORHOOD GENERATION TIME
//         auto t0 = clk::now();
//         auto neigh = neighborhood_factories[idx_neigh](I, sol);
//         auto t1 = clk::now();
//         time_neigh += std::chrono::duration<double, std::milli>(t1 - t0).count();

//         // MOVE SELECTION TIME
//         t0 = clk::now();
//         auto const &mov = step_function(*(neigh));
//         t1 = clk::now();
//         time_step += std::chrono::duration<double, std::milli>(t1 - t0).count();

//         if (!mov.has_value())
//             break;

//         auto actual_move = mov.value();

//         // ACCEPTANCE / APPLY TIME
//         t0 = clk::now();
//         double delta = neigh->calc_delta(actual_move);
//         bool accept = do_i_accept(delta, T, rng);
//         t1 = clk::now();
//         time_accept += std::chrono::duration<double, std::milli>(t1 - t0).count();

//         if (accept)
//         {
//             sol = neigh->apply(actual_move);
//             f += delta;

//             if (f < best_f)
//             {
//                 best_f = f;
//                 best_sol = sol;
//             }
//         }

//         T *= cooling;
//         iter++;
//     }

//     auto t_end_total = clk::now();
//     double total_ms = std::chrono::duration<double, std::milli>(t_end_total - t_start_total).count();

//     // ---------------------------
//     // Final timing output
//     // ---------------------------
//     std::cout << "\n--- SA timing summary ---\n"
//               << "Total time:                " << total_ms    << " ms\n"
//               << "Neighborhood creation:     " << time_neigh  << " ms\n"
//               << "Move selection:            " << time_step   << " ms\n"
//               << "Acceptance & eval:         " << time_accept << " ms\n"
//               << "Iterations:                " << iter        << "\n"
//               << "Initial f:                 " << utils::objective(I, initial_sol) << "\n"
//               << "Best f:                    " << best_f      << "\n"
//               << "---------------------------\n";

//     return best_sol;
// }

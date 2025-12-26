
#include <fstream>
#include <filesystem>
#include <iostream>
#include <vector>
#include <map>
#include "solvers.hpp"
#include "structures.hpp"
#include "path_utils.hpp"

struct RES
{
    int N;
    std::string method;
    double objective;
    double duration; // ms
    fs::path instance_path;
};

void write_csv_results(const fs::path &output_path, const std::vector<RES> &results)
{
    std::ofstream out(output_path);
    if (!out.is_open())
    {

        std::cerr << "Error: cannot open results file: "
                  << output_path << std::endl;
        return;
    }

    out << "path,N,method,duration,objective\n";

    for (const auto &r : results)
    {
        out << r.instance_path.string() << ","
            << r.N << ","
            << r.method << ","
            << r.duration << ","
            << r.objective << "\n";
    }
    out.close();
}

int main(int argc, char **argv)
{
    std::cout << "Enter MAIN A I compare dc rc " << std::endl;
    auto [base_instances, base_output, _] = parse_paths(argc, argv);

    std::vector<int> Ns{50, 100};

    std::vector<RES> all_res;

    std::vector<Neighborhood::NeighborhoodFactory> neighborhoods = {
        [](const Instance &I, const Solution &s)
        { return std::make_unique<IntraRouteNeighborhood>(I, s); },

        // [](const Instance &I, const Solution &s)
        // { return std::make_unique<PairRelocateNeighborhood>(I, s); },

        [](Instance const &I, Solution const &s)
        { return std::make_unique<RequestMove>(I, s); },

        [](const Instance &I, const Solution &s)
        { return std::make_unique<TwoOptNeighborhood>(I, s); }};

    for (auto N : Ns)
    {
        ImprovementThreshold stopping_criterion((double)N / 10);
        std::string N_str = std::to_string(N);
        std::cout << "Enter " << N_str << std::endl;
        std::filesystem::path subdir = base_instances / N_str / "test";
        auto instance_paths = get_some_instance_paths(subdir, 4);

        for (auto const &instance : instance_paths)
        {
            std::cout << " -" << std::endl;
            Solution dr_sol;
            Instance I(instance, "maxmin");
            { // Deterministic

                std::cout<<"dc";
                Timer t;
                dr_sol = DC::construction(I);
                double exec_time = t.get_time();
                if (!dr_sol.is_solution_feasible(I))
                {
                    std::cerr << "Solution is not feasible ERROR?";
                }
                double objective = utils::objective(I, dr_sol);
                
                RES res;
                res.duration = exec_time;
                res.method = "DC ";
                res.instance_path = instance;
                res.objective = objective;
                res.N = N;
                all_res.push_back(res);
            }
            { // Random
                std::cout<<" rc";
                
                Timer t;
                auto rc_sol = RC::construction(I, 0.2); // after tuning
                double exec_time = t.get_time();
                if (!rc_sol.is_solution_feasible(I))
                {
                    std::cerr << "Solution is not feasible ERROR?";
                }
                double objective = utils::objective(I, rc_sol);
                
                RES res;
                res.duration = exec_time;
                res.method = "RC";
                res.instance_path = instance;
                res.objective = objective;
                res.N = N;
                all_res.push_back(res);
            }
            { // Beam Search
                
                std::cout<<" bs";
                Timer t;
                auto bs_sol = BS::beam_search(I, 0.9, 10); // after tuning
                double exec_time = t.get_time();
                if (!bs_sol.is_solution_feasible(I))
                {
                    std::cerr << "Solution is not feasible ERROR?";
                }
                double objective = utils::objective(I, bs_sol);
                
                RES res;
                res.duration = exec_time;
                res.method = "BS";
                res.instance_path = instance;
                res.objective = objective;
                res.N = N;
                all_res.push_back(res);
            }
            
            // { // LS
                
            //     std::cout<<" ls";
            //     Timer t;
            //     int total_iter = 0;
            //     auto sol_ls = LS::local_search(I, dr_sol, neighborhoods[1], StepFunction::first_improvement, stopping_criterion);
            //     double exec_time = t.get_time();
            //     if (!sol_ls.is_solution_feasible(I))
            //     {
            //         std::cerr << "Solution is not feasible ERROR?";
            //     }
            //     double objective = utils::objective(I, sol_ls);
                
            //     RES res;
            //     res.duration = exec_time;
            //     res.method = "LS";
            //     res.instance_path = instance;
            //     res.objective = objective;
            //     res.N = N;
            //     all_res.push_back(res);
            // }
            // { // VND
                
            //     std::cout<<" vnd";
            //     Timer t;
            //     int total_iter = 0;
            //     auto sol_vnd = VND::vnd(I, dr_sol, neighborhoods, StepFunction::first_improvement, stopping_criterion);
            //     double exec_time = t.get_time();
            //     if (!sol_vnd.is_solution_feasible(I))
            //     {
            //         std::cerr << "Solution is not feasible ERROR?";
            //     }
            //     double objective = utils::objective(I, sol_vnd);
                
            //     RES res;
            //     res.duration = exec_time;
            //     res.method = "VND";
            //     res.instance_path = instance;
            //     res.objective = objective;
            //     res.N = N;
            //     all_res.push_back(res);
            // }
            // { // SA
                
            //     std::cout<<" sa";
            //     int total_iters = 0;
            //     Timer t;
            //     auto sol_sa = SA::simulated_annealing(I, dr_sol, neighborhoods,
            //                                           100.0, 1e-3, 0.97, StepFunction::first_improvement, stopping_criterion, &total_iters);
            //     if (total_iters == 0)
            //         std::cout << "I have done no move \n";

            //     double exec_time = t.get_time();
            //     double objective = utils::objective(I, sol_sa);

            //     RES res;
            //     res.duration = exec_time;
            //     res.method = "SA";
            //     res.instance_path = instance;
            //     res.objective = objective;
            //     res.N = N;
            //     all_res.push_back(res);
            // }

            // {
            //     std::cout<<" LN";
            //     Timer t;
            //     auto ln_sol = LN::large_neighborhood(I, dr_sol, 5); // after tuning
            //     double exec_time = t.get_time();
            //     if (!ln_sol.is_solution_feasible(I))
            //     {
            //         std::cerr << "Solution is not feasible ERROR?";
            //     }
            //     double objective = utils::objective(I, ln_sol);
                
            //     RES res;
            //     res.duration = exec_time;
            //     res.method = "LN";
            //     res.instance_path = instance;
            //     res.objective = objective;
            //     res.N = N;
            //     all_res.push_back(res);
            // }
            {
                std::cout<<" GA";
                Timer t;
                auto ga_sol = GA::genetic_algorithm(I, 10,1, 20); // after tuning
                double exec_time = t.get_time();
                if (!ga_sol.is_solution_feasible(I))
                {
                    std::cerr << "Solution is not feasible ERROR?";
                }
                double objective = utils::objective(I, ga_sol);
                
                RES res;
                res.duration = exec_time;
                res.method = "GA";
                res.instance_path = instance;
                res.objective = objective;
                res.N = N;
                all_res.push_back(res);
            }
        }
        std::cout << "\n";
    }

    write_csv_results(base_output / "dr_rc_vnd_timing.csv", all_res);
}

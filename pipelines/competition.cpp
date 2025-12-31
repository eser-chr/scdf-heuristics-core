#include <map>
#include <chrono>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <cassert>
#include "solvers.hpp"
#include "structures.hpp"
#include "path_utils.hpp"

struct IFile
{
    int N;
    std::string name;
    std::filesystem::path input;
};

struct RES
{
    double duration;
    double objective;
    std::string method;
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

    out << "method,duration,objective\n";

    for (const auto &r : results)
    {
        out << r.method << ","
            << r.duration << ","
            << r.objective
            << "\n";
    }
    out.close();
}

// class Timer
// {
//     std::chrono::high_resolution_clock::time_point start;

// public:
//     Timer() : start(std::chrono::high_resolution_clock::now()) {}
//     double get_time()
//     {
//         auto end = std::chrono::high_resolution_clock::now();
//         return std::chrono::duration<double, std::milli>(end - start).count();
//     }
// };

void deal_with_single_file(IFile const &ifile, std::map<std::string, bool> const &what_to_run, std::filesystem::path const &output)
{
    auto const &instance_path = ifile.input;
    auto const &instance_name = ifile.name;
    auto output_folder = output / instance_name;
    std::filesystem::create_directories(output_folder);

    std::vector<Neighborhood::NeighborhoodFactory> neighborhoods = {
        [](const Instance &I, const Solution &s)
        { return std::make_unique<IntraRouteNeighborhood>(I, s); },

        [](Instance const &I, Solution const &s)
        { return std::make_unique<RequestMove>(I, s); },

        [](const Instance &I, const Solution &s)
        { return std::make_unique<TwoOptNeighborhood>(I, s); }};

    try
    {
        Instance I(instance_path);

        Solution sol_drc = DC::construction(I);
        assert(sol_drc.is_solution_feasible(I));

        sol_drc.write_solution(output_folder / "dc.txt", instance_name);
        MaxIterations stopping(500);

        Solution sol_rc, sol_ls, sol_beam, sol_vnd, sol_sa, sol_grasp, sol_ga, sol_ln;
        std::vector<RES> results;
        results.reserve(10);

        if (what_to_run.at("RANDOM"))
        {
            Timer t;
            sol_rc = RC::construction(I, 0.2); // Param as in tuning
            double time = t.get_time();
            assert(sol_rc.is_solution_feasible(I));

            RES res{time, utils::objective(I, sol_rc), "RANDOM"};
            sol_rc.write_solution(output_folder / "rc.txt", instance_name);
            results.push_back(res);
        }
        // ---------------- LS ----------------
        if (what_to_run.at("LS"))
        {
            Timer t;
            sol_ls = LS::local_search(
                I,
                sol_drc,
                neighborhoods[0],
                StepFunction::first_improvement,
                stopping);
            double time = t.get_time();
            assert(sol_ls.is_solution_feasible(I));
            double obj = utils::objective(I, sol_ls);
            results.push_back(RES{time, obj, "LS"});
            sol_ls.write_solution(output_folder / "ls.txt", instance_name);
        }

        // ---------------- BS ----------------
        if (what_to_run.at("BS"))
        {
            Timer t;
            sol_beam = BS::beam_search(I, 1.0, 7);
            double time = t.get_time();
            assert(sol_beam.is_solution_feasible(I));

            double obj = utils::objective(I, sol_beam);
            results.push_back(RES{time, obj, "BS"});
            sol_beam.write_solution(output_folder / "bs.txt", instance_name);
        }

        // ---------------- VND ----------------
        if (what_to_run.at("VND"))
        {
            MaxIterations stopping_vnd(10000);
            Timer t;
            sol_vnd = VND::vnd(I,
                               sol_drc,
                               neighborhoods,
                               StepFunction::first_improvement,
                               stopping_vnd);
            double time = t.get_time();
            assert(sol_vnd.is_solution_feasible(I));

            double obj = utils::objective(I, sol_vnd);
            results.push_back(RES{time, obj, "VND"});
            sol_vnd.write_solution(output_folder / "vnd.txt", instance_name);
        }

        // ---------------- SA ----------------
        if (what_to_run.at("SA"))
        {
            MaxIterations stopping_sa(500);
            Timer t;

            sol_sa = SA::simulated_annealing(I,
                                             sol_drc,
                                             neighborhoods,
                                             1,
                                             0.1,
                                             0.995,
                                             StepFunction::random_step,
                                             stopping_sa);
            double time = t.get_time();
            assert(sol_sa.is_solution_feasible(I));

            double obj = utils::objective(I, sol_sa);
            results.push_back(RES{time, obj, "SA"});
            sol_sa.write_solution(output_folder / "sa.txt", instance_name);
        }

        // ---------------- GRASP ----------------
        if (what_to_run.at("GRASP"))
        {
            Timer t;

            auto constructor = [&](const Instance &I)
            {
                return GRASP::randomized_constructor_simple(I, 1.0, 0.5);
            };

            MaxIterations stopping_outer(100);
            MaxIterations stopping_local(2000);

            sol_grasp = GRASP::grasp(
                I,
                constructor,
                neighborhoods,
                StepFunction::first_improvement,
                stopping_outer,
                stopping_local);
            double time = t.get_time();
            assert(sol_grasp.is_solution_feasible(I));

            double obj = utils::objective(I, sol_grasp);
            results.push_back(RES{time, obj, "GRASP"});
            sol_grasp.write_solution(output_folder / "grasp.txt", instance_name);
        }
        // ---------------- LN ----------------
        if (what_to_run.at("LN"))
        {
            Timer t;
            sol_ln = LN::large_neighborhood(I, sol_vnd, 4, 20, 8, 8);
            double time = t.get_time();
            assert(sol_ln.is_solution_feasible(I));

            double obj = utils::objective(I, sol_ln);
            results.push_back(RES{time, obj, "LN"});
            sol_ln.write_solution(output_folder / "ln.txt", instance_name);
        }
        // ---------------- GA ----------------
        if (what_to_run.at("GA"))
        {
            Timer t;

            sol_ga = GA::genetic_algorithm(I, 10, 1, 20, 8);
            double time = t.get_time();
            assert(sol_ga.is_solution_feasible(I));

            double obj = utils::objective(I, sol_ga);
            results.push_back(RES{time, obj, "GA"});
            sol_ga.write_solution(output_folder / "ga.txt", instance_name);
        }

        write_csv_results(output_folder / "results.csv", results);
    }
    catch (const std::exception &e)
    {
        std::cerr << "ERROR processing " << instance_name << ": " << e.what() << "\n";
    }
    catch (...)
    {
        std::cerr << "ERROR: Unknown error processing " << instance_name << "\n";
    }
}

int main(int argc, char *argv[])
{
    auto [base_instances, output, _] = parse_paths(argc, argv);
    std::vector<IFile> files{
        // IFile{100, "instance61_nreq100_nveh2_gamma91", base_instances / "heuristics/instances/100/competition/instance61_nreq100_nveh2_gamma91.txt"},
        IFile{1000, "instance61_nreq1000_nveh20_gamma879", base_instances / "heuristics/instances/1000/competition/instance61_nreq1000_nveh20_gamma879.txt"},
        // IFile{2000, "instance61_nreq2000_nveh40_gamma1829", base_instances / "heuristics/instances/2000/competition/instance61_nreq2000_nveh40_gamma1829.txt"},
    };

    std::map<std::string, bool> what_to_run = {
        {"RANDOM", true},
        {"LS", false},
        {"BS", true},
        {"VND", true},
        {"SA", false},
        {"GRASP", false},
        {"LN", true},
        {"GA", true},

    };

    for (auto const &ifile : files)
    {
        deal_with_single_file(ifile, what_to_run, output);
    }
}

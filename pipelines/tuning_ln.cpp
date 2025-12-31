/**
 * Tuning script for large neighborhood !!!
 */

#include <fstream>
#include <filesystem>
#include <iostream>
#include <vector>
#include "solvers.hpp"
#include "structures.hpp"
#include "path_utils.hpp"

void write_binary(const std::vector<std::vector<double>> &data, const std::string &filename)
{
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open())
    {
        std::cerr << "Error: cannot open binary file: " << filename << "\n";
        return;
    }
    size_t rows = data.size();
    file.write(reinterpret_cast<const char *>(&rows), sizeof(rows));
    for (const auto &row : data)
    {
        size_t cols = row.size();
        file.write(reinterpret_cast<const char *>(&cols), sizeof(cols));
        file.write(reinterpret_cast<const char *>(row.data()), cols * sizeof(double));
    }
}

struct RES
{
    int N;
    int k;
    int bw1;
    int bw2;
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

    out << "path,N,k,bw1,bw2,duration,objective\n";

    for (const auto &r : results)
    {
        out << r.instance_path.string() << ","
            << r.N << ","
            << r.k << ","
            << r.bw1 << ","
            << r.bw2 << ","
            << r.duration << ","
            << r.objective << "\n";
    }

    out.close();
}

struct combo
{
    int k;
    int bw1;
    int bw2;
};

combo get_combo(size_t counter, std::vector<int> const &ks, 
                std::vector<int> const &bw1s, std::vector<int> const &bw2s)
{
    auto Nks = ks.size();
    auto Nbw1s = bw1s.size();
    auto Nbw2s = bw2s.size();

    size_t idx_bw2 = counter % Nbw2s;
    size_t idx_bw1 = (counter / Nbw2s) % Nbw1s;
    size_t idx_k = counter / (Nbw2s * Nbw1s);

    if (idx_k >= Nks) {
        std::cerr << "Error: counter out of range in get_combo\n";
        return combo{ks[0], bw1s[0], bw2s[0]};  // Or throw exception
    }

    return combo{ks[idx_k], bw1s[idx_bw1], bw2s[idx_bw2]};
}


int main(int argc, char **argv)
{
    std::cout << "Start" << std::endl;
    auto [base_instances, base_output, _] = parse_paths(argc, argv);

    std::vector<int> Ns{50};
    std::vector<int> ks{2};
    std::vector<int> bw1s{5};
    std::vector<int> bw2s{5, 10};
    std::vector<RES> all_res;
    std::vector<std::vector<double>> objectives_over_time;

    for (auto N : Ns)
    {
        // ImprovementThreshold stopping_criterion((double)N / 10);
        std::string N_str = std::to_string(N);
        std::cout << "Enter " << N_str << std::endl;
        std::filesystem::path subdir = base_instances / N_str / "test";
        auto instance_paths = get_some_instance_paths(subdir, 1);

        for (auto const &instance : instance_paths)
        {
            size_t total_size = ks.size() * bw1s.size() * bw2s.size();
            Instance I(instance, "jain");

            auto dr_sol = DC::construction(I);

            for (size_t counter = 0; counter < total_size; ++counter)
            {

                auto [k, bw1, bw2] = get_combo(counter, ks, bw1s, bw2s);
                std::vector<double> objective_over_time;

                Timer t;
                auto ln_sol = LN::large_neighborhood(I, dr_sol, k, 20, bw1, bw2, &objective_over_time);
                double exec_time = t.get_time();
                if (!ln_sol.is_solution_feasible(I))
                {
                    std::cerr << "Solution is not feasible ERROR? \n";
                }
                double objective = utils::objective(I, ln_sol);
                RES res;
                res.duration = exec_time;
                res.instance_path = instance;
                res.objective = objective;
                res.N = N;
                res.k = k;
                res.bw1 = bw1;
                res.bw2 = bw2;
                all_res.push_back(res);
                objectives_over_time.push_back(objective_over_time);
            }
        }
        std::cout << "\n";
    }

    write_csv_results(base_output / "tuning_ln.csv", all_res);
    write_binary(objectives_over_time, (base_output / "tuning_ln_objectives_over_time.bin").string());
}

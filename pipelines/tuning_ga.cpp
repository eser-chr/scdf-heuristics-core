/**
 * Tuning script for genetic algorithm!!!
 */

#include <fstream>
#include <filesystem>
#include <iostream>
#include <vector>
#include <map>
#include "solvers.hpp"
#include "structures.hpp"
#include "path_utils.hpp"

void write_binary(const std::vector<std::vector<double>> &data, const std::string &filename)
{
    std::ofstream file(filename, std::ios::binary);
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
    int k1;
    int k2;
    int beam_width;
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

    out << "path,N,k1,k2,beam_width,duration,objective\n";

    for (const auto &r : results)
    {
        out << r.instance_path.string() << ","
            << r.N << ","
            << r.k1 << ","
            << r.k2 << ","
            << r.beam_width << ","
            << r.duration << ","
            << r.objective << "\n";
    }

    out.close();
}

struct combo
{
    int k1;
    int k2;
    int beam_width;
};

combo get_combo(size_t counter, std::vector<int> const &k1s, std::vector<int> const &k2s, std::vector<int> const &beam_widths)
{
    // static size_t counter = 0;

    auto Nk1 = k1s.size();
    auto Nk2 = k2s.size();
    auto Nbeam_width = beam_widths.size();

    size_t idx_beam = counter % Nbeam_width;
    size_t idx_k2 = (counter / Nbeam_width) % Nk2;
    size_t idx_k1 = counter / (Nbeam_width * Nk2);

    // counter++;

    return combo{k1s[idx_k1], k2s[idx_k2], beam_widths[idx_beam]};
}

int main(int argc, char **argv)
{
    std::cout << "Start" << std::endl;
    auto [base_instances, base_output, _] = parse_paths(argc, argv);

    std::vector<int> Ns{50};
    std::vector<int> k1s{5};
    std::vector<int> k2s{0};
    std::vector<int> beam_widths{5,10};
    std::vector<RES> all_res;
    std::vector<std::vector<double>> objectives_over_time;

    for (auto N : Ns)
    {
        ImprovementThreshold stopping_criterion((double)N / 10);
        std::string N_str = std::to_string(N);
        std::cout << "Enter " << N_str << std::endl;
        std::filesystem::path subdir = base_instances / N_str / "test";
        auto instance_paths = get_some_instance_paths(subdir, 1);

        for (auto const &instance : instance_paths)
        {
            size_t total_size = k1s.size() * k2s.size() * beam_widths.size();
            std::cout << " -" << std::endl;
            Instance I(instance, "maxmin");

            for (size_t counter = 0; counter < total_size; ++counter)
            {

                auto [k1, k2, beam_width] = get_combo(counter, k1s, k2s, beam_widths);
                std::vector<double> objective_over_time;

                Timer t;
                auto ga_sol = GA::genetic_algorithm(I, k1, k2, 20, beam_width, &objective_over_time); // after tuning
                double exec_time = t.get_time();
                if (!ga_sol.is_solution_feasible(I))
                {
                    std::cerr << "Solution is not feasible ERROR?";
                }
                double objective = utils::objective(I, ga_sol);
                RES res;
                res.duration = exec_time;
                res.instance_path = instance;
                res.objective = objective;
                res.N = N;
                res.k1 = k1;
                res.k2 = k2;
                res.beam_width = beam_width;
                all_res.push_back(res);
                objectives_over_time.push_back(objective_over_time);
            }
        }
        std::cout << "\n";
    }

    write_csv_results(base_output / "tuning_ga.csv", all_res);
    // write_binary(objectives_over_time, base_output / "tuning_ga_objectives_over_time.csv");
    write_binary(objectives_over_time, (base_output/"tuning_ga_objectives_over_time.bin").string());

}

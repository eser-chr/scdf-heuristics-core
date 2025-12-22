
#include <fstream>
#include <filesystem>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
namespace fs = std::filesystem;

auto get_instance_paths(const fs::path &folder)
{
    if (!fs::exists(folder))
        throw std::runtime_error("folder Does not exist " + folder.string());

    std::vector<fs::path> instances{};
    for (const auto &entry : fs::directory_iterator(folder))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".txt")
        {
            // std::cout << entry.path() << std::endl;
            instances.push_back(entry.path());
        }
    }
    return instances;
}

auto get_some_instance_paths(const fs::path &folder, int num_of_instances)
{
    auto paths = get_instance_paths(folder);
    if (num_of_instances >= paths.size())
    {
        return paths;
    }

    std::mt19937 rng(std::random_device{}());
    std::vector<int> indices(paths.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    std::vector<fs::path> to_return;
    to_return.reserve(num_of_instances);
    for (size_t i = 0; i < num_of_instances; i++)
        to_return.push_back(paths[indices[i]]);

    return to_return;
}

struct ParsedPaths
{
    fs::path base_instances;
    fs::path base_output;
    int idx;
};

ParsedPaths parse_paths(int argc, char **argv)
{
    if (argc != 3 && argc !=4)
    {
        std::cerr << "Usage: ./run <instances_path> <output_path>\n";
        std::exit(1);
    }

    fs::path instances = argv[1];
    fs::path output = argv[2];
    int n{};
    if(argc==4){
        n = std::atoi(argv[3]);
    }



    if (!fs::exists(instances) || !fs::is_directory(instances))
    {
        std::cerr << "Error: instances_path is not a valid directory\n";
        std::exit(1);
    }

    if (!fs::exists(output))
    {
        fs::create_directories(output);
    }

    return ParsedPaths{instances, output, n};
}

class Timer
{
    std::chrono::high_resolution_clock::time_point start;

public:
    Timer() : start(std::chrono::high_resolution_clock::now()) {}
    double get_time()
    {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start).count();
    }
};
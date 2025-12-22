#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "structures.hpp"


void Instance::printme() const
{
    std::cout << name << " \n";
    std::cout << "---------" << "\n";
    std::cout << n << " " << nK << " " << C << " " << gamma << " " << rho << std::endl;
}

Instance::Instance(const std::string &path)
{
    load_from_file(path);
    if (!is_instance_correct())
    {
        throw std::runtime_error("Instance \"" + path + "\" failed validation.");
    }
}

bool Instance::is_instance_correct()
{
    auto fail = [&](const std::string &msg)
    {
        std::cerr << "[Instance check failed] " << msg << "\n";
        return false;
    };

    // --- basic sizes ---
    if (n <= 0)
        return fail("n <= 0");
    if (nK <= 0)
        return fail("nK <= 0");
    if (gamma <= 0)
        return fail("gamma <= 0");
    if (C <= 0)
        return fail("C <= 0");

    // --- structural consistency ---
    if ((int)coords.size() != 2*n+1)
        return fail("coords.size() != n (" + std::to_string(coords.size()) + " vs " + std::to_string(n) + ")");

    if ((int)demands.size() != n)
        return fail("demands.size() != n (" + std::to_string(demands.size()) + " vs " + std::to_string(n) + ")");

    if ((int)dist.size() != 2*n+1)
        return fail("dist.size() != n (" + std::to_string(dist.size()) + " vs " + std::to_string(n) + ")");

    for (int i = 0; i < n; i++) {
        if ((int)dist[i].size() != 2*n+1)
            return fail("dist[" + std::to_string(i) + "].size() != n");
    }

    // --- demand sanity ---
    for (int i = 0; i < n; i++) {
        if (demands[i] < 0)
            return fail("demands[" + std::to_string(i) + "] < 0");
        if (demands[i] > C)
            return fail("demands[" + std::to_string(i) + "] > C (" + std::to_string(demands[i]) + " > " + std::to_string(C) + ")");
    }

    // --- gamma & nK ---
    if (gamma > n)
        return fail("gamma > n (" + std::to_string(gamma) + " > " + std::to_string(n) + ")");

    if (nK > n)
        return fail("nK > n (" + std::to_string(nK) + " > " + std::to_string(n) + ")");

    // --- coordinate sanity ---
    for (int i = 0; i < n; i++) {
        if (!std::isfinite(coords[i].x))
            return fail("coords[" + std::to_string(i) + "].x is not finite");
        if (!std::isfinite(coords[i].y))
            return fail("coords[" + std::to_string(i) + "].y is not finite");
    }

    // --- distance matrix sanity ---
    for (int u = 0; u < n; u++) {

        if (dist[u][u] != 0)
            return fail("dist[" + std::to_string(u) + "][" + std::to_string(u) + "] != 0");

        for (int v = 0; v < n; v++) {
            if (dist[u][v] < 0)
                return fail("dist[" + std::to_string(u) + "][" + std::to_string(v) + "] < 0");

            if (dist[u][v] != dist[v][u])
                return fail("dist not symmetric at ("
                            + std::to_string(u) + "," + std::to_string(v) + ")");
        }
    }

    return true;
}

void Instance::load_from_file(const std::string &path)
{
    std::ifstream f(path);
    if (!f)
        throw std::runtime_error("Could not open instance file: " + path);

    name = path;

    // read all lines
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(f, line))
        lines.push_back(line);

    // header
    {
        std::stringstream ss(lines[0]);
        ss >> n >> nK >> C >> gamma >> rho;
    }

    // find markers
    size_t idx_dem = -1, idx_loc = -1;
    for (size_t i = 0; i < lines.size(); ++i)
    {
        if (lines[i].rfind("# demands", 0) == 0)
            idx_dem = i;
        if (lines[i].rfind("# request locations", 0) == 0)
            idx_loc = i;
    }
    if (idx_dem < 0 || idx_loc < 0)
        throw std::runtime_error("Bad instance file: missing markers");

    // parse demands
    {
        std::vector<int> tokens;
        for (int i = idx_dem + 1; i < idx_loc; i++)
        {
            std::stringstream ss(lines[i]);
            int x;
            while (ss >> x)
                tokens.push_back(x);
        }
        if ((int)tokens.size() != n)
            throw std::runtime_error("Bad file: wrong number of demands");
        demands = tokens;
    }

    f.close();



    // Build up the rest structure
    // number of nodes including depot
    int nV = 1 + 2 * n;

    coords.reserve(nV);
    for (int i = idx_loc + 1; i <= idx_loc + nV; ++i)
    {
        std::stringstream ss(lines[i]);
        double x, y;
        ss >> x >> y;
        coords.emplace_back(x, y);
    }

    // build distance matrix
    dist.assign(nV, std::vector<double>(nV, 0.0));
    for (int u = 0; u < nV; u++)
    {
        for (int v = u+1; v < nV; v++)
        {
            
            double dx = coords[u].x - coords[v].x;
            double dy = coords[u].y - coords[v].y;
            dist[u][v] = std::sqrt(dx * dx + dy * dy);
            dist[v][u] = dist[u][v];
        }
    }

    // request_of_node + load_change
    request_of_node.assign(nV, -1);
    load_change.assign(nV, 0);

    for (int i = 0; i < n; i++)
    {
        int p = 1 + i;
        int d = 1 + n + i;

        request_of_node[p] = i;
        request_of_node[d] = i;

        load_change[p] = +demands[i];
        load_change[d] = -demands[i];
    }
}
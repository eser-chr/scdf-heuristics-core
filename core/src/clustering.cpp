#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <limits>
#include <numeric>
#include "structures.hpp"

clusters::ClusterCenters::ClusterCenters(const Instance &I, std::vector<int> const &reqs)
    : centers(I.nK, {0.0, 0.0})
{
    std::vector<int> indices = reqs;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(indices.begin(), indices.end(), gen);

    int used = std::min(I.nK, (int)indices.size());
    for (int k = 0; k < used; ++k)
    {
        int r = indices[k];
        centers[k] = I.coords[1 + r];
    }

    // if fewer requests than vehicles, remaining centers stay at (0,0)
}

// update using only selected requests
void clusters::ClusterCenters::update_centers(const Instance &I,
                                              const std::vector<int> &reqs,
                                              const std::vector<int> &assign)
{
    centers.assign(I.nK, {0.0, 0.0});
    std::vector<int> count(I.nK, 0);

    for (std::size_t i = 0; i < reqs.size(); ++i)
    {
        int r = reqs[i];   // original request index
        int k = assign[i]; // cluster index

        if (k < 0 || k >= I.nK)
        {
            std::cerr << "Invalid cluster index k=" << k << " for req " << reqs[i] << "\n";
            std::abort();
        }

        centers[k].x += I.coords[1 + r].x;
        centers[k].y += I.coords[1 + r].y;
        count[k]++;
    }

    for (int k = 0; k < I.nK; ++k)
    {
        if (count[k] > 0)
        {
            centers[k].x /= count[k];
            centers[k].y /= count[k];
        }
    }
}

// assignment only over reqs
std::vector<int> clusters::balanced_assign(
    const Instance &I,
    const clusters::ClusterCenters &C,
    const std::vector<int> &reqs,
    double target_load)
{
    std::vector<int> assign(reqs.size(), -1);
    std::vector<double> load(I.nK, 0.0);

    for (std::size_t i = 0; i < reqs.size(); ++i)
    {
        int r = reqs[i];
        double best_score = std::numeric_limits<double>::infinity();
        int best_k = 0;

        for (int k = 0; k < I.nK; ++k)
        {
            // double distance = calc_dist2(I.coords[1 + r], C.centers[k]);
            double distance = numerical::calc_distance_between_nodes(I.coords[1 + r], C.centers[k]);
            double load_after = load[k] + I.demands[r];
            double load_dev = std::abs(load_after - target_load);
            double score = distance*distance + load_dev * load_dev;

            if (score < best_score)
            {
                best_score = score;
                best_k = k;
            }
        }

        assign[i] = best_k;
        load[best_k] += I.demands[r];
    }

    return assign;
}

// k-means over a subset of requests
std::vector<int> clusters::balanced_kmeans(
    const Instance &I,
    const std::vector<int> &reqs,
    int iters = 20,
    int restarts = 20)
{
    double total_dem = 0.0;
    for (int r : reqs)
        total_dem += I.demands[r];
    double target_load = total_dem / I.nK;

    std::vector<int> best_assign(reqs.size(), 0);
    double best_score = std::numeric_limits<double>::infinity();

    for (int s = 0; s < restarts; ++s)
    {
        ClusterCenters C(I, reqs);
        std::vector<int> assign(reqs.size(), 0);

        for (int it = 0; it < iters; ++it)
        {
            assign = balanced_assign(I, C, reqs, target_load);
            C.update_centers(I, reqs, assign);
        }

        double score = 0.0;
        for (std::size_t i = 0; i < reqs.size(); ++i)
        {
            int r = reqs[i];
            int k = assign[i];
            double distance = numerical::calc_distance_between_nodes(I.coords[1 + r], C.centers[k]);
            score += distance*distance;
        }

        if (score < best_score)
        {
            best_score = score;
            best_assign = assign;
        }
    }

    return best_assign; // size = reqs.size(), cluster index per selected request
}

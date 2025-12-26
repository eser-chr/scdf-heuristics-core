#pragma once
#include <string>
#include <vector>
#include <string>
#include <random>
#include <functional>
namespace gt // GeneralTypes
{
    // using precision_t = double;
    template <typename T>
    using Matrix = std::vector<std::vector<T>>;
    struct Coords
    {
        double x;
        double y;
    };
};

class Instance
{
    using Matrix = gt::Matrix<double>;

public:
    std::string name;
    int n;      // # requests
    int nK;     // # vehicles
    int C;      // vehicle capacity
    int gamma;  // min served requests
    double rho; // fairness weight
    std::string fairness;

    std::vector<int> demands;         // demands per request size n
    Matrix dist;                      // (1 + 2n) x (1 + 2n) distance over nodes, precompute
    std::vector<int> request_of_node; // size (1 + 2n)
    std::vector<int> load_change;     // how much does the load change if a vehicle passes throught that node
    std::vector<gt::Coords> coords;

    Instance(std::string const &path, std::string const & fairness = "jain");
    void printme() const;

private:
    void load_from_file(const std::string &path);
    bool is_instance_correct();
};

struct Solution
{
    gt::Matrix<int> routes; // size = nK
    std::vector<double> routes_distances;
    double total_distance;
    double sum_of_squares; // sum of dists[i]**2. For efficient delta eval
    std::string fairness;

    Solution() = default;
    void write_solution(const std::string &path, const std::string &instance_name) const;
    bool is_solution_feasible(Instance const& I);
    void compute_cached_values_from_routes(Instance const& I);
};

namespace utils
{

    std::vector<int> calc_route_cargo(
        const Instance &inst, const std::vector<int> &route);

    // int route_distance(
    //     const Instance &inst, const std::vector<int> &route);

    double calc_route_distance(Instance const& I, Solution const & sol, int route_idx);
    double calc_route_distance(Instance const& I, std::vector<int> const& route);

    std::vector<double> all_route_distances(
        const Instance &inst, const Solution &sol);

    double jain_fairness(Instance const &I,
                         const std::vector<double> &dists);

    double max_min_fairness(Instance const &I, std::vector<double> const &dists);

    double gini_cefficient_nominator(Instance const &I, std::vector<double> const &dists);
    double gini_cefficient(Instance const &I, std::vector<double> const &dists);

    bool is_route_feasible(
        const Instance &inst, const std::vector<int> &route);

    double objective(
        const Instance &inst, const Solution &sol);

    // overload if different fairness is required.
    // double objective(Instance const &I, Solution const &sol, std::function<double(Instance const &I, std::vector<double> const &dists)> fairness_func);

    std::vector<double> calc_my_metric(const Instance &I, double a);
} // namespace utils

namespace numerical
{
    template <typename T>
    std::vector<int> argsort(const std::vector<T> &org);

    template <typename T>
    T select_uniformly(const std::vector<T> &org, std::mt19937 &rng);

    double calc_distance_between_nodes(gt::Coords const&p1, gt::Coords const& p2);

};

namespace clusters
{

    class ClusterCenters
    {
    public:
        using centers_t = std::vector<gt::Coords>;
        centers_t centers;

        ClusterCenters(const Instance &I, std::vector<int> const &reqs);
        // Update via an assignment. Calcs the point average of that clusters requests(pickup only!)
        void update_centers(const Instance &I, const std::vector<int> &reqs, const std::vector<int> &assign);
    };

    std::vector<int> balanced_assign(
        const Instance &I,
        const ClusterCenters &C,
        const std::vector<int> &reqs,
        double target_load);

    std::vector<int> balanced_kmeans(
        const Instance &I,
        const std::vector<int> &reqs,
        int iters,
        int restarts);
};
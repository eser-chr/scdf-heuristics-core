#pragma once
#include <vector>
#include <functional>
#include <memory>
#include <optional>
#include <algorithm>
#include <numeric>
#include <algorithm>
#include "structures.hpp"

struct GenericMove
{
    int type;              // 1 = intra, 2 = pair-relocate, 3 = 2-opt
    std::vector<int> data; // payload
};

struct PickupDeliveryInfo
{
    int p_idx;
    int d_idx;
    int req;
    int pickup_node;
    int delivery_node;
};

std::vector<PickupDeliveryInfo>
pickup_delivery_positions(const Instance &I, const std::vector<int> &route);

class Neighborhood
{
protected:
    size_t MAX_TRIES_RANDOM = 100;

public:
    using NeighborhoodFactory =
        std::function<std::unique_ptr<Neighborhood>(const Instance &, const Solution &)>;
    using NeighborhoodFactories = std::vector<NeighborhoodFactory>;

    const Instance &I;
    const Solution &sol;
    double f;

    Neighborhood(const Instance &I_, const Solution &sol_)
        : I(I_), sol(sol_), f(utils::objective(I_, sol_)) {}

    virtual ~Neighborhood() = default;

    virtual std::vector<GenericMove> generate() const = 0;
    virtual std::optional<GenericMove> generate_random(std::mt19937 &rng) const = 0;
    virtual bool is_valid(const GenericMove &mov) const = 0;
    virtual double calc_delta(const GenericMove &mov) const = 0;
    virtual double calc_delta_jain(const GenericMove &mov) const = 0;
    virtual double calc_delta_maxmin(const GenericMove &mov) const = 0;
    virtual double calc_delta_gini(const GenericMove &mov) const = 0;
    virtual Solution apply(const GenericMove &mov) const = 0;
};

class IntraRouteNeighborhood : public Neighborhood
{
public:
    int const type = 1;
    std::string const name = "IntraRoute";
    IntraRouteNeighborhood(const Instance &I_, const Solution &sol_)
    : Neighborhood(I_, sol_) {}
    std::vector<GenericMove> generate() const override;
    std::optional<GenericMove> generate_random(std::mt19937 &rng) const;
    bool is_valid(const GenericMove &mov) const override;
    double calc_delta(const GenericMove &mov) const override;
    double calc_delta_jain(const GenericMove &mov) const override;
    double calc_delta_maxmin(const GenericMove &mov) const override;
    double calc_delta_gini(const GenericMove &mov) const override;
    Solution apply(const GenericMove &mov) const override;
    
};
/**
 * We swap in between routes requests. Since a request has to be delivered fully by a only
 * one vehicle then we move the whole request.
 */
// class PairRelocateNeighborhood : public Neighborhood
// {
//     public:
//     int const type = 2;
//     std::string const name = "PairLocate";
    
//     PairRelocateNeighborhood(const Instance &I_, const Solution &sol_)
//     : Neighborhood(I_, sol_) {}
//     std::vector<GenericMove> generate() const override;
//     std::optional<GenericMove> generate_random(std::mt19937 &rng) const;
//     bool is_valid(const GenericMove &mov) const override;
//     double calc_delta(const GenericMove &mov) const override;
//     Solution apply(const GenericMove &mov) const override;
    
// };

/**
 * Moves the request from one vehicle to another. Move is described by the triplet 
 * from , to , request id. 
 * The nodes of that request are placed at the end. Maybe good for escaping neighborhoods.
 */
class RequestMove : public Neighborhood
{
    public:
    int const type = 2;
    std::string const name = "PairLocate";
    
    RequestMove(const Instance &I_, const Solution &sol_)
    : Neighborhood(I_, sol_) {}
    std::vector<GenericMove> generate() const override;
    std::optional<GenericMove> generate_random(std::mt19937 &rng) const;
    bool is_valid(const GenericMove &mov) const override;
    double calc_delta(const GenericMove &mov) const override;
    double calc_delta_jain(const GenericMove &mov) const override;
    double calc_delta_maxmin(const GenericMove &mov) const override;
    double calc_delta_gini(const GenericMove &mov) const override;
    Solution apply(const GenericMove &mov) const override;
    
};




class TwoOptNeighborhood : public Neighborhood
{
    public:
    int const type = 3;
    std::string const name = "TwoOpt";
    
    TwoOptNeighborhood(const Instance &I_, const Solution &sol_)
    : Neighborhood(I_, sol_) {}
    std::vector<GenericMove> generate() const override;
    std::optional<GenericMove> generate_random(std::mt19937 &rng) const;
    bool is_valid(const GenericMove &mov) const override;
    double calc_delta(const GenericMove &mov) const override;
    double calc_delta_jain(const GenericMove &mov) const override;
    double calc_delta_maxmin(const GenericMove &mov) const override;
    double calc_delta_gini(const GenericMove &mov) const override;
    Solution apply(const GenericMove &mov) const override;

};
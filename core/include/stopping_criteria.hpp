#pragma once
#include <cmath>
#include <vector>
#include <memory>

class StoppingCriterion
{
public:
    virtual ~StoppingCriterion() = default;
    virtual void reset() = 0;
    virtual bool operator()(int iteration, double f) = 0;
};

class MaxIterations : public StoppingCriterion
{
    int max_iters;

public:
    explicit MaxIterations(int m)
        : max_iters(m) {}

    inline void reset() override {}
    inline bool operator()(int iteration, double) override
    {
        return iteration >= max_iters;
    }
};

class ObjectiveThreshold : public StoppingCriterion
{
    double threshold;

public:
    explicit ObjectiveThreshold(double thr)
        : threshold(thr) {}

    inline void reset() override {}

    inline bool operator()(int, double f) override
    {
        return f <= threshold;
    }
};

class ImprovementThreshold : public StoppingCriterion
{
    double eps;
    double last_f;
    bool first;

public:
    explicit ImprovementThreshold(double e)
        : eps(e), last_f(0.0), first(true) {}

    inline void reset() override
    {
        first = true;
    }

    inline bool operator()(int, double f) override
    {
        if (first)
        {
            first = false;
            last_f = f;
            return false;
        }
        double diff = std::fabs(last_f - f);
        last_f = f;
        return diff < eps;
    }
};

// MAYBE COMBO OF CRITERIA
// -----------------------


class AnyCriterion : public StoppingCriterion
{
    std::vector<std::shared_ptr<StoppingCriterion>> criteria;

public:
    AnyCriterion(std::initializer_list<std::shared_ptr<StoppingCriterion>> lst)
        : criteria(lst) {}

    void reset() override
    {
        for (auto &c : criteria)
            c->reset();
    }

    bool operator()(int iteration, double f) override
    {
        for (auto &c : criteria)
            if ((*c)(iteration, f))
                return true;
        return false;
    }
};

class AllCriterion : public StoppingCriterion
{
    std::vector<std::shared_ptr<StoppingCriterion>> criteria;

public:
    AllCriterion(std::initializer_list<std::shared_ptr<StoppingCriterion>> lst)
        : criteria(lst) {}

    void reset() override
    {
        for (auto &c : criteria)
            c->reset();
    }

    bool operator()(int iteration, double f) override
    {
        for (auto &c : criteria)
            if (!(*c)(iteration, f))
                return false;
        return true;
    }
};
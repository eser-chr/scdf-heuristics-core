#pragma once
#include <optional>
#include <functional>
#include "neighborhoods.hpp"
#include <random>

namespace StepFunction
{
    using Return_t = std::optional<GenericMove>;
    using Func = std::function<Return_t(const Neighborhood &, std::mt19937 &)>;


    inline Return_t first_improvement(const Neighborhood &N, std::mt19937 &rng)
    {
        const size_t _MAX_TRIES = 1000; // random tries before giving up

        for (size_t t = 0; t < _MAX_TRIES; t++)
        {
            auto mov = N.generate_random(rng);
            if (!mov.has_value())
                continue;

            if (!N.is_valid(*mov))
                continue;

            if (N.calc_delta(*mov) < 0)
                return mov;
        }


        return std::nullopt;
    }

    inline Return_t best_improvement(const Neighborhood &N, std::mt19937 &)
    {
        std::vector<GenericMove> moves =  N.generate();

        double best_delta = 0;
        Return_t best = std::nullopt;

        for (const auto &m : moves)
        {
            if (!N.is_valid(m))
                continue;
            double d = N.calc_delta(m);
            if (d < best_delta)
            {
                best_delta = d;
                best = m;
            }
        }
        return best;
    }

    inline Return_t random_step(const Neighborhood &N, std::mt19937 &rng)
    {
        // static thread_local std::mt19937 rng(std::random_device{}());
        return N.generate_random(rng);
    }

} // namespace StepFunction
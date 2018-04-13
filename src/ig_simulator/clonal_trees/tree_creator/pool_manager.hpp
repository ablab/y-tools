//
// Created by Andrew Bzikadze on 4/9/17.
//

#pragma once

#include <cstddef>
#include "cartesian_tree.hpp"
#include "ig_simulator_utils.hpp"
#include "ig_simulator_config.hpp"

namespace ig_simulator {

class AbstractPoolManager {
protected:
    Treap<> pool;
    mutable std::bernoulli_distribution ret_to_pool_distr;
    size_t max_index;

public:
    explicit AbstractPoolManager(double ret_prob):
            pool(),
            ret_to_pool_distr(check_probability(ret_prob)),
            max_index(1)
    {
        pool.Insert(0, 1);
    }

    AbstractPoolManager& operator=(const AbstractPoolManager&) = delete;
    AbstractPoolManager& operator=(AbstractPoolManager&&) = delete;

    static std::unique_ptr<AbstractPoolManager> CreatePoolManager(PoolManagerStrategy pool_manager_strategy, double ret_prob);

    size_t MaxIndex() const { return max_index; }
    void Erase(size_t index) {
        VERIFY(index < max_index);
        pool.Erase(index);
    }

    size_t Size() const { return pool.Size(); }
    virtual std::pair<size_t, bool> GetIndex(size_t n_insert) = 0;
};

using AbstractPoolManagerCPtr = std::unique_ptr<AbstractPoolManager>;


class UniformPoolManager final : public AbstractPoolManager {
public:
    explicit UniformPoolManager(double ret_prob):
        AbstractPoolManager(ret_prob)
    { }

    std::pair<size_t, bool> GetIndex(size_t n_insert) override;
};

class WideTreePoolManager final : public AbstractPoolManager {
public:
    explicit WideTreePoolManager(double ret_prob):
        AbstractPoolManager(ret_prob)
    { }

    std::pair<size_t, bool> GetIndex(size_t n_insert) override;
};

class DeepTreePoolManager final : public AbstractPoolManager {
public:
    explicit DeepTreePoolManager(double ret_prob):
        AbstractPoolManager(ret_prob)
    { }

    std::pair<size_t, bool> GetIndex(size_t n_insert) override;
};

} // End namespace ig_simulator

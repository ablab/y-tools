//
// Created by Andrew Bzikadze on 4/9/17.
//

#include <limits>
#include <ig_simulator_config.hpp>
#include "pool_manager.hpp"

namespace ig_simulator {

std::unique_ptr<AbstractPoolManager>
AbstractPoolManager::CreatePoolManager(const PoolManagerStrategy pool_manager_strategy, double ret_prob) {
    if (pool_manager_strategy == PoolManagerStrategy::UniformPoolManager) {
        return std::unique_ptr<AbstractPoolManager>(new UniformPoolManager(ret_prob));
    } else if (pool_manager_strategy == PoolManagerStrategy::DeepTreePoolManager) {
        return std::unique_ptr<AbstractPoolManager>(new DeepTreePoolManager(ret_prob));
    } else if (pool_manager_strategy == PoolManagerStrategy::WideTreePoolManager) {
        return std::unique_ptr<AbstractPoolManager>(new WideTreePoolManager(ret_prob));
    }
    VERIFY(false);
}

std::pair<size_t, bool> UniformPoolManager::GetIndex(size_t n_insert) {
    size_t raw_index = random_index(1, pool.Sum());
    size_t index, freq;
    std::tie(index, freq) = pool.LowerBound(raw_index);
    VERIFY(freq == 1);

    for (size_t i = 0; i < n_insert; ++i) {
        pool.Insert(max_index++, 1);
    }

    bool ret_to_pool = ret_to_pool_distr(MTSingleton::GetInstance());
    if (not ret_to_pool) {
        pool.Erase(index, freq);
    }
    return { index, ret_to_pool };
}

std::pair<size_t, bool> WideTreePoolManager::GetIndex(size_t n_insert) {
    size_t raw_index = random_index(1, pool.Sum());
    size_t index, freq;
    std::tie(index, freq) = pool.LowerBound(raw_index);

    for (size_t i = 0; i < n_insert; ++i) {
        pool.Insert(max_index++, 1);
    }

    bool ret_to_pool = ret_to_pool_distr(MTSingleton::GetInstance());
    if (ret_to_pool) {
        pool.SetFreq(index, freq, freq + 1);
    } else {
        pool.Erase(index, freq);
    }
    return { index, ret_to_pool };
}

std::pair<size_t, bool> DeepTreePoolManager::GetIndex(size_t n_insert) {
    size_t raw_index = random_index(1, pool.Sum());
    size_t index, freq;
    std::tie(index, freq) = pool.LowerBound(raw_index);

    size_t new_freq = freq + 1;
    if (freq < std::numeric_limits<unsigned int>::max()) {
        new_freq += static_cast<size_t>(static_cast<double>(new_freq) * 0.5);
    } else {
        new_freq += static_cast<size_t>(sqrt(static_cast<double>(new_freq)));
    }

    for (size_t i = 0; i < n_insert; ++i) {
        pool.Insert(max_index++, new_freq);
    }

    bool ret_to_pool = ret_to_pool_distr(MTSingleton::GetInstance());
    if (ret_to_pool) {
        pool.SetFreq(index, freq, new_freq);
    } else {
        pool.Erase(index, freq);
    }
    return { index, ret_to_pool };
}

} // End namespace ig_simulator

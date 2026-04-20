#pragma once

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "node_global.hpp"

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

// forward declaration
inline std::vector<int>
make_children_from_order_with_placeholder(
    const std::vector<int>& order,
    const std::unordered_map<int, int>* placeholder_nodes,
    const std::vector<int>* local_to_global);

// =====================================================
// Strong DSD fallback: heuristic Shannon split
// =====================================================
inline int strong_else_decompose(
    const std::string& mf,
    const std::vector<int>& order,
    int depth,
    int pivot_node,
    const std::vector<int>* local_to_global,
    const std::unordered_map<int, int>* placeholder_nodes,
    int (*strong_rec)(
        const std::string&,
        const std::vector<int>&,
        int,
        const std::vector<int>*,
        const std::unordered_map<int, int>*))
{
    (void)pivot_node;

    const int n = static_cast<int>(order.size());
    if (n == 0)
    {
        return new_node(mf, {});
    }

    if ((size_t(1) << n) != mf.size())
    {
        return new_node(mf, {});
    }

    for (char c : mf)
    {
        if (c != '0' && c != '1')
        {
            throw std::runtime_error("strong_else_decompose: mf contains non-binary char");
        }
    }

    kitty::dynamic_truth_table tt(n);
    kitty::create_from_binary_string(tt, mf);

    std::vector<uint32_t> support;
    support.reserve(n);
    for (uint32_t i = 0; i < static_cast<uint32_t>(n); ++i)
    {
        if (kitty::has_var(tt, i))
        {
            support.push_back(i);
        }
    }

    if (support.empty())
    {
        return new_node(mf, {});
    }

    uint32_t var_min = support.front();
    int best_score = std::numeric_limits<int>::max();
    for (auto var : support)
    {
        const auto co00 = kitty::cofactor0(tt, var);
        const auto co11 = kitty::cofactor1(tt, var);
        int s_co0 = 0;
        int s_co1 = 0;
        for (uint32_t s_sd = 0; s_sd < co00.num_vars(); ++s_sd)
        {
            if (kitty::has_var(co00, s_sd))
            {
                s_co0 += 1;
            }
        }
        for (uint32_t s_sd = 0; s_sd < co11.num_vars(); ++s_sd)
        {
            if (kitty::has_var(co11, s_sd))
            {
                s_co1 += 1;
            }
        }
        const int add = s_co0 + s_co1;
        if (add < best_score)
        {
            best_score = add;
            var_min = var;
        }
    }

    const uint32_t chosen_var = var_min;
    size_t order_index = static_cast<size_t>(n - 1u - var_min);
    std::vector<int> reordered_order = order;
    if (order_index < reordered_order.size())
    {
        const int pivot_var = reordered_order.at(order_index);
        reordered_order.erase(reordered_order.begin() + static_cast<long>(order_index));
        reordered_order.insert(reordered_order.begin(), pivot_var);
    }

    if (var_min != static_cast<uint32_t>(n - 1))
    {
        kitty::swap_inplace(tt, var_min, static_cast<uint32_t>(n - 1));
        var_min = static_cast<uint32_t>(n - 1);
        order_index = 0;
    }

    std::string indent(static_cast<size_t>(depth) * 2, ' ');
    std::cout << indent << "⚙️ strong_else_dec: pick var_id "
              << reordered_order.front()
              << " (kitty var " << chosen_var
              << ", order index " << order_index << ") as pivot (MSB)\n";
    std::cout << indent << "   order(MSB->LSB) = { ";
    for (int v : reordered_order)
    {
        std::cout << v << " ";
    }
    std::cout << "}\n";
    std::cout << indent << "   reordered tt = " << kitty::to_binary(tt) << "\n";

    const auto children = make_children_from_order_with_placeholder(
        reordered_order, placeholder_nodes, local_to_global);
    if (children.empty())
    {
        return new_node(mf, {});
    }
    const int pivot = children.front();

    auto build_cofactor_tt = [&](bool value) {
        const uint32_t new_vars = static_cast<uint32_t>(n - 1);
        kitty::dynamic_truth_table new_tt(new_vars);
        const uint64_t limit = 1ull << new_vars;

        for (uint64_t x = 0; x < limit; ++x)
        {
            uint64_t old = 0;
            for (uint32_t b = 0; b < new_vars; ++b)
            {
                if (x & (1ull << b))
                {
                    const uint32_t old_pos = b < var_min ? b : b + 1u;
                    old |= (1ull << old_pos);
                }
            }
            if (value)
            {
                old |= (1ull << var_min);
            }

            if (kitty::get_bit(tt, old))
            {
                kitty::set_bit(new_tt, x);
            }
        }

        std::vector<int> out_order = reordered_order;
        if (!out_order.empty())
        {
            out_order.erase(out_order.begin());
        }

        return std::pair<std::string, std::vector<int>>{
            kitty::to_binary(new_tt),
            std::move(out_order)};
    };

    auto [f_pos, pos_order] = build_cofactor_tt(true);
    auto [f_neg, neg_order] = build_cofactor_tt(false);

    const int pos_node =
        strong_rec(f_pos, pos_order, depth + 1, local_to_global, placeholder_nodes);
    const int neg_node =
        strong_rec(f_neg, neg_order, depth + 1, local_to_global, placeholder_nodes);

    const int pos_term = new_node("1000", {pivot, pos_node});
    const int neg_term = new_node("0010", {pivot, neg_node});

    return new_node("1110", {pos_term, neg_term});
}
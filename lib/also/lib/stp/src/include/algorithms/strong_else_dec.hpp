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
    const std::vector<int>& order,          // MSB -> LSB (你原约定)
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
    const int n = static_cast<int>(order.size());
    if (n == 0) return new_node(mf, {});
    if ((size_t(1) << n) != mf.size()) return new_node(mf, {});

    for (char c : mf)
    {
        if (c != '0' && c != '1')
            throw std::runtime_error("strong_else_decompose: mf contains non-binary char");
    }

    // ---------------------------
    // 1) 构建 tt，并建立 kitty<->语义映射
    // ---------------------------
    kitty::dynamic_truth_table tt(n);
    kitty::create_from_binary_string(tt, mf);

    // 约定：kitty var i 对应 "LSB->MSB" 变量位；
    // 而你的 order 是 MSB->LSB。
    // 所以：kitty i <-> order[n-1-i]
    std::vector<int> kitty_to_semantic(n);
    for (int i = 0; i < n; ++i)
        kitty_to_semantic[i] = order[n - 1 - i];

    auto recompute_support = [&](kitty::dynamic_truth_table const& t) {
        std::vector<uint32_t> supp;
        supp.reserve(t.num_vars());
        for (uint32_t i = 0; i < t.num_vars(); ++i)
            if (kitty::has_var(t, i)) supp.push_back(i);
        return supp;
    };

    std::vector<uint32_t> support = recompute_support(tt);
    if (support.empty())
        return new_node(mf, {});

    // ---------------------------
    // 2) 启发式选 var_min：min (supp(co0)+supp(co1))
    // ---------------------------
    uint32_t var_min = support.front();
    int best_score = std::numeric_limits<int>::max();

    for (auto var : support)
    {
        const auto co0 = kitty::cofactor0(tt, var);
        const auto co1 = kitty::cofactor1(tt, var);

        int s0 = 0, s1 = 0;
        for (uint32_t i = 0; i < co0.num_vars(); ++i) if (kitty::has_var(co0, i)) ++s0;
        for (uint32_t i = 0; i < co1.num_vars(); ++i) if (kitty::has_var(co1, i)) ++s1;

        const int score = s0 + s1;
        if (score < best_score)
        {
            best_score = score;
            var_min = var;
        }
    }

    const uint32_t chosen_var = var_min;

    // ---------------------------
    // 3) 让 chosen_var 变成 MSB (kitty 的最高位)，并同步更新映射
    // ---------------------------
    const uint32_t msb = static_cast<uint32_t>(n - 1);
    if (var_min != msb)
    {
        kitty::swap_inplace(tt, var_min, msb);
        std::swap(kitty_to_semantic[var_min], kitty_to_semantic[msb]);
        var_min = msb;
    }

    // pivot 的语义变量就是 kitty MSB 对应的语义变量
    const int pivot_var_semantic = kitty_to_semantic[var_min];

    // 用 pivot_var_semantic 把 order 重排：把 pivot 放到 MSB 位置（order.front）
    std::vector<int> reordered_order = order;
    {
        auto it = std::find(reordered_order.begin(), reordered_order.end(), pivot_var_semantic);
        if (it != reordered_order.end())
        {
            int v = *it;
            reordered_order.erase(it);
            reordered_order.insert(reordered_order.begin(), v);
        }
    }

    if (pivot_node < 0)
    {
        std::cout << "⚙️ strong_else_dec: pick var_id "
                  << pivot_var_semantic
                  << " (kitty var " << chosen_var
                  << " -> moved to MSB kitty " << var_min << ")\n";

        std::cout << "   order(MSB->LSB) = { ";
        for (int v : reordered_order) std::cout << v << " ";
        std::cout << "}\n";

        std::cout << "   reordered tt = " << kitty::to_binary(tt) << "\n";
    }

    // ---------------------------
    // 4) 构 children，pivot 是 children.front()
    // ---------------------------
    const auto children = make_children_from_order_with_placeholder(
        reordered_order, placeholder_nodes, local_to_global);
    if (children.empty()) return new_node(mf, {});
    const int pivot = children.front();

    // ---------------------------
    // 5) 构建 cofactor 的 (tt + order) ，并做“正确的缩维/重编号”
    // ---------------------------
    auto build_cofactor = [&](bool value) -> std::pair<std::string, std::vector<int>>
    {
        const uint32_t new_vars = static_cast<uint32_t>(n - 1); // 去掉 MSB 一维
        kitty::dynamic_truth_table new_tt(new_vars);

        // 新 tt 的 kitty index b (0..new_vars-1) 对应旧 tt 的同 index b
        // 因为我们把 pivot 放在 MSB(var_min=msb)，去掉的就是最高位。
        const uint64_t limit = 1ull << new_vars;
        for (uint64_t x = 0; x < limit; ++x)
        {
            uint64_t old = x;
            if (value) old |= (1ull << var_min); // 设旧 tt 的 MSB
            if (kitty::get_bit(tt, old))
                kitty::set_bit(new_tt, x);
        }

        // 去掉 pivot 后的 kitty->semantic 映射（先只是删掉 MSB）
        std::vector<int> map_after_drop;
        map_after_drop.reserve(new_vars);
        for (uint32_t i = 0; i < new_vars; ++i)
            map_after_drop.push_back(kitty_to_semantic[i]); // 0..new_vars-1 原样保留

        // 重新算 support，并把 new_tt 压缩到最小 support
        std::vector<uint32_t> new_support = recompute_support(new_tt);

        // 如果 support 没满，做 shrink，并同步 shrink 映射
        if (new_support.size() != new_vars)
        {
            kitty::dynamic_truth_table reduced_tt(static_cast<uint32_t>(new_support.size()));
            const uint64_t reduced_limit = 1ull << new_support.size();

            for (uint64_t x = 0; x < reduced_limit; ++x)
            {
                uint64_t old = 0;
                for (uint32_t b = 0; b < new_support.size(); ++b)
                {
                    if (x & (1ull << b))
                        old |= (1ull << new_support[b]); // 映射回 new_tt 的变量位
                }
                if (kitty::get_bit(new_tt, old))
                    kitty::set_bit(reduced_tt, x);
            }

            // ✅ 正确的 reduced_order：直接用“被保留的变量”的语义映射
            std::vector<int> reduced_order; // MSB->LSB
            reduced_order.reserve(new_support.size());

            // reduced_tt 的 kitty 变量顺序是 0..k-1（LSB->MSB）
            // 我们要输出 order 是 MSB->LSB
            for (int i = static_cast<int>(new_support.size()) - 1; i >= 0; --i)
            {
                const uint32_t kept_old_kitty = new_support[static_cast<uint32_t>(i)];
                reduced_order.push_back(map_after_drop[kept_old_kitty]);
            }

            return {kitty::to_binary(reduced_tt), std::move(reduced_order)};
        }
        else
        {
            // 没有缩维：order 就是 map_after_drop 的 MSB->LSB
            std::vector<int> out_order;
            out_order.reserve(new_vars);
            for (int i = static_cast<int>(new_vars) - 1; i >= 0; --i)
                out_order.push_back(map_after_drop[static_cast<uint32_t>(i)]);

            return {kitty::to_binary(new_tt), std::move(out_order)};
        }
    };

    auto [f_pos, pos_order] = build_cofactor(true);
    auto [f_neg, neg_order] = build_cofactor(false);

    const int pos_node =
        strong_rec(f_pos, pos_order, depth + 1, local_to_global, placeholder_nodes);
    const int neg_node =
        strong_rec(f_neg, neg_order, depth + 1, local_to_global, placeholder_nodes);

    // ITE / MUX 节点：你原来写的是 "10101100"（4输入 LUT 的 ITE）
    // children: {pivot, neg, pos} 这里按你原本的约定保留
    return new_node("10101100", {pivot, neg_node, pos_node});
}

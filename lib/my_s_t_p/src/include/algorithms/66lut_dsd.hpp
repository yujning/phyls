#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <unordered_set>
#include "node_global.hpp"

// =====================================================
// Debug switch
// =====================================================
inline bool LUT66_DSD_DEBUG_PRINT = false;

// =====================================================
// Pretty print: TT + order (must be paired)
// =====================================================
inline void print_tt_with_order_66(
    const std::string& title,
    const std::string& tt,
    const std::vector<int>& order,
    int depth = 0)
{
    if (!LUT66_DSD_DEBUG_PRINT) return;

    std::string indent((size_t)depth * 2, ' ');
    std::cout << indent << "📌 " << title << "\n";
    std::cout << indent << "   TT    = " << tt << "\n";
    std::cout << indent << "   order = { ";
    for (int v : order) std::cout << v << " ";
    std::cout << "}\n";
}

// =====================================================
// 66-LUT DSD result
// =====================================================
struct Lut66DsdResult {
    bool found = false;
    size_t L = 0;              // = 2^{|Mx|}
    std::string Mx;            // = block0 + block1, length 2*L
    std::string My;            // length 2^{|My|}

    // 0-based positions in current order
    std::vector<int> mx_pos;
    std::vector<int> my_pos;

    // var IDs in MSB->LSB order (by position)
    std::vector<int> mx_vars_msb2lsb;
    std::vector<int> my_vars_msb2lsb;

    // optional debug
    std::string block0;
    std::string block1;
    std::string reordered_tt;  // concatenated blocks by My assignment (My|Mx)
};

// =====================================================
// next_combination: comb is strictly increasing indices in [0..m-1]
// =====================================================
inline bool next_combination_66(std::vector<int>& comb, int m)
{
    int k = (int)comb.size();
    for (int i = k - 1; i >= 0; --i) {
        if (comb[i] < m - k + i) {
            ++comb[i];
            for (int j = i + 1; j < k; ++j) {
                comb[j] = comb[j - 1] + 1;
            }
            return true;
        }
    }
    return false;
}

// =====================================================
// Extract block (Mx subfunction) for a fixed My assignment
// - order positions: 0..n-1 (0 is MSB position)
// - mx_pos_sorted / my_pos_sorted must be sorted ascending by position (MSB->LSB)
// - my_assignment encoded in My's MSB->LSB order
// - mx_index encoded in Mx's MSB->LSB order
// =====================================================
inline std::string extract_block_for_mx_66(
    const std::string& mf,
    int n,
    const std::vector<int>& mx_pos_sorted,
    const std::vector<int>& my_pos_sorted,
    uint64_t my_assignment)
{
    int k = (int)mx_pos_sorted.size();
    int m = (int)my_pos_sorted.size();

    const size_t sub_size = 1ull << k;
    std::string sub(sub_size, '0');

    for (size_t mx_index = 0; mx_index < sub_size; ++mx_index)
    {
        uint64_t full_index = 0;

        // 1) place My bits (MSB->LSB)
        for (int t = 0; t < m; ++t)
        {
            int pos = my_pos_sorted[t];
            int bit = (int)((my_assignment >> (m - 1 - t)) & 1ull);
            full_index |= (uint64_t(bit) << (n - 1 - pos));
        }

        // 2) place Mx bits (MSB->LSB)
        for (int t = 0; t < k; ++t)
        {
            int pos = mx_pos_sorted[t];
            int bit = (int)((mx_index >> (k - 1 - t)) & 1ull);
            full_index |= (uint64_t(bit) << (n - 1 - pos));
        }

        sub[mx_index] = mf[(size_t)full_index];
    }

    return sub;
}

// =====================================================
// Helper: print current candidate split info
// =====================================================
inline void print_candidate_info_66(
    int depth,
    int k,
    int m,
    const std::vector<int>& mx_vars_msb2lsb,
    const std::vector<int>& my_vars_msb2lsb)
{
    if (!LUT66_DSD_DEBUG_PRINT) return;

    std::string indent((size_t)depth * 2, ' ');
     std::cout << indent << "🔎 66-DSD 尝试 |My|=" << m << " (<=6), |Mx|=" << k << " (<=5)\n";
    std::cout << indent << "   My={ ";
    for (int v : my_vars_msb2lsb) std::cout << v << " ";
    std::cout << "}  Mx={ ";
    for (int v : mx_vars_msb2lsb) std::cout << v << " ";
    std::cout << "}\n";
}

// =====================================================
// ★ 66-LUT Strong DSD (subset enumeration, single split)
// Constraint: |My| <= K, |Mx| <= K-1, and same split rule as strong_dsd.hpp (exactly 2 blocks)
// Search: iterate k within feasible range; enumerate Mx subset by var_id ascending.
// =====================================================
inline Lut66DsdResult run_66lut_dsd_by_mx_subset(
    const std::string& mf,
    const std::vector<int>& order,
    int depth_for_print = 0)
{
    Lut66DsdResult out;

    int n = (int)order.size();
    if (n <= 1) return out;
    if ((size_t(1) << n) != mf.size()) return out;

    // (var_id, pos) sorted by var_id (ascending) for stable enumeration
    struct VarPos { int var; int pos; };
    std::vector<VarPos> vp;
    vp.reserve(n);
    for (int pos = 0; pos < n; ++pos) vp.push_back({order[pos], pos});
    std::sort(vp.begin(), vp.end(), [](const VarPos& x, const VarPos& y){
        return x.var < y.var;
    });

    // K-LUT constraints:
    // m = n-k <= K  => k >= n-K
    // k <= K-1
    const int lut_k = LUT_PAIR_FANIN_LIMIT;
    const int mx_k = mx_lut_fanin_limit();
    int k_min = std::max(1, n - lut_k);
    int k_max = std::min(mx_k, n - 1);

    for (int k = k_max; k >= k_min; --k)
    {
        std::vector<int> comb(k);
        for (int i = 0; i < k; ++i) comb[i] = i;

        while (true)
        {
            // mx positions from comb
            std::vector<int> mx_pos;
            mx_pos.reserve(k);
            for (int idx : comb) mx_pos.push_back(vp[idx].pos);

            // my positions = rest
            std::vector<int> my_pos;
            my_pos.reserve(n - k);
            {
                std::vector<char> is_mx(n, 0);
                for (int p : mx_pos) is_mx[p] = 1;
                for (int p = 0; p < n; ++p) if (!is_mx[p]) my_pos.push_back(p);
            }

            // sort by position (MSB->LSB)
            std::sort(mx_pos.begin(), mx_pos.end());
            std::sort(my_pos.begin(), my_pos.end());

            // build vars list in MSB->LSB by position
            std::vector<int> mx_vars_msb2lsb;
            std::vector<int> my_vars_msb2lsb;
            mx_vars_msb2lsb.reserve(mx_pos.size());
            my_vars_msb2lsb.reserve(my_pos.size());
            for (int p : mx_pos) mx_vars_msb2lsb.push_back(order[p]);
            for (int p : my_pos) my_vars_msb2lsb.push_back(order[p]);

            int m = n - k; // |My|
            if (m > lut_k || k > mx_k) {
                if (!next_combination_66(comb, (int)vp.size())) break;
                continue;
            }

            print_candidate_info_66(depth_for_print, k, m, mx_vars_msb2lsb, my_vars_msb2lsb);

            uint64_t my_count = 1ull << m;
            uint64_t L = 1ull << k;

            if (LUT66_DSD_DEBUG_PRINT)
            {
                std::string reordered;
                reordered.reserve((size_t)my_count * (size_t)L);

                for (uint64_t y = 0; y < my_count; ++y) {
                    std::string block = extract_block_for_mx_66(mf, n, mx_pos, my_pos, y);
                    reordered += block;
                }

                std::vector<int> reordered_order;
                reordered_order.reserve(n);
                reordered_order.insert(reordered_order.end(), my_vars_msb2lsb.begin(), my_vars_msb2lsb.end());
                reordered_order.insert(reordered_order.end(), mx_vars_msb2lsb.begin(), mx_vars_msb2lsb.end());

                print_tt_with_order_66("候选 split 的重排 TT (My|Mx)", reordered, reordered_order, depth_for_print);
            }

            std::unordered_map<std::string, int> block_index;
            std::vector<std::string> blocks;
            blocks.reserve(2);

            std::string My;
            My.reserve((size_t)my_count);

            bool too_many = false;

            for (uint64_t y = 0; y < my_count; ++y)
            {
                std::string block = extract_block_for_mx_66(mf, n, mx_pos, my_pos, y);

                auto it = block_index.find(block);
                if (it == block_index.end())
                {
                    if (blocks.size() >= 2) { too_many = true; break; }
                    int id = (int)blocks.size();
                    block_index.emplace(block, id);
                    blocks.push_back(block);
                    My.push_back(id == 0 ? '1' : '0');
                }
                else
                {
                    My.push_back(it->second == 0 ? '1' : '0');
                }
            }

            if (!too_many && blocks.size() == 2)
            {
                out.found = true;
                out.L = (size_t)L;
                out.Mx = blocks[0] + blocks[1];
                out.My = My;

                out.mx_pos = mx_pos;
                out.my_pos = my_pos;
                out.mx_vars_msb2lsb = mx_vars_msb2lsb;
                out.my_vars_msb2lsb = my_vars_msb2lsb;

                out.block0 = blocks[0];
                out.block1 = blocks[1];

                // reordered_tt (My|Mx)
                std::string reordered;
                reordered.reserve((size_t)my_count * (size_t)L);
                for (uint64_t y = 0; y < my_count; ++y) {
                    std::string block = extract_block_for_mx_66(mf, n, mx_pos, my_pos, y);
                    reordered += block;
                }
                out.reordered_tt = reordered;

                if (LUT66_DSD_DEBUG_PRINT)
                {
                    std::string indent((size_t)depth_for_print * 2, ' ');
                    std::cout << indent << "✅ 命中 66-LUT Strong DSD split\n";
                    std::cout << indent << "   block0 = " << blocks[0] << "\n";
                    std::cout << indent << "   block1 = " << blocks[1] << "\n";
                    print_tt_with_order_66("当前 split 的 My", My, my_vars_msb2lsb, depth_for_print);
                }

                return out;
            }

            if (!next_combination_66(comb, (int)vp.size())) break;
        }
    }

    return out;
}

// =====================================================
// children (MSB->LSB) + placeholder support (same spirit as strong_dsd.hpp)
// =====================================================
inline std::vector<int> make_children_from_order_with_placeholder_66(
    const std::vector<int>& order_msb2lsb,
    const std::unordered_map<int, int>* placeholder_nodes)
{
    std::vector<int> children;
    children.reserve(order_msb2lsb.size());

    // 关键：按 MSB -> LSB 直接放
    for (int var_id : order_msb2lsb)
    {
        if (placeholder_nodes) {
            auto it = placeholder_nodes->find(var_id);
            if (it != placeholder_nodes->end()) {
                children.push_back(it->second);
                continue;
            }
        }

          auto global_it = PLACEHOLDER_BINDINGS.find(var_id);
        if (global_it != PLACEHOLDER_BINDINGS.end())
        {
            children.push_back(global_it->second);
            continue;
        }
        children.push_back(new_in_node(var_id));
    }

    return children;
}

// =====================================================
// ★ 入口：只做一次 66-LUT DSD，然后按 “MY -> MX(含MY占位)” 方式建 DAG
// 原变量默认顺序：n,n-1,...,1 (MSB->LSB)
// =====================================================
inline bool run_66lut_dsd_and_build_dag(const TT& root_tt)
{
    int n = (int)root_tt.order.size();
    if ((size_t(1) << n) != root_tt.f01.size()) return false;

    RESET_NODE_GLOBAL();
    int max_var_id = root_tt.order.empty() ? 0
                                           : *std::max_element(root_tt.order.begin(),
                                                               root_tt.order.end());
    ORIGINAL_VAR_COUNT = max_var_id;
    auto res = run_66lut_dsd_by_mx_subset(root_tt.f01, root_tt.order,
                                          /*depth_for_print=*/0);
    if (!res.found) {
        if (LUT66_DSD_DEBUG_PRINT) std::cout << "❌ No valid 66-LUT Strong DSD split\n";
        return false;
    }

    if (LUT66_DSD_DEBUG_PRINT)
    {
        std::cout << "✅ 66-LUT DSD FOUND\n";
        std::cout << "   |My|=" << res.my_vars_msb2lsb.size() << "  |Mx|=" << res.mx_vars_msb2lsb.size() << "\n";
        std::cout << "   My vars (MSB->LSB): { ";
        for (int v : res.my_vars_msb2lsb) std::cout << v << " ";
        std::cout << "}\n";
        std::cout << "   Mx vars (MSB->LSB): { ";
        for (int v : res.mx_vars_msb2lsb) std::cout << v << " ";
        std::cout << "}\n";
        std::cout << "   MY = " << res.My << "\n";
        std::cout << "   MX = " << res.Mx << "\n";
    }

   // ----- build MY node (<=K vars) -----
    {
        std::vector<int> order_my = res.my_vars_msb2lsb; // MSB->LSB
        auto children_my = make_children_from_order_with_placeholder_66(order_my, nullptr);
        int my_node = new_node(res.My, children_my);

   // ----- build MY node (<=K vars) -----
        int k = (int)res.mx_vars_msb2lsb.size();
        int my_local_id = allocate_placeholder_var_id(nullptr);


        // 先把 MY 的变量做成一个集合，便于过滤
        std::unordered_set<int> my_var_set(
            res.my_vars_msb2lsb.begin(),
            res.my_vars_msb2lsb.end()
        );

            // MX 的输入顺序：MY 在 MSB，其后是“非 MY 的原变量”；同时去重，防止占位符重复出现
            std::vector<int> order_mx;
            order_mx.reserve(1 + res.mx_vars_msb2lsb.size());

            std::unordered_set<int> seen_vars;
            order_mx.push_back(my_local_id);
            seen_vars.insert(my_local_id);

            for (int v : res.mx_vars_msb2lsb)
            {
                // ★ 关键：排除属于 MY 的变量；如果 MX 列表里出现重复，也要跳过
                if (my_var_set.count(v) == 0 && seen_vars.insert(v).second)
                    order_mx.push_back(v);
            }
        

        std::unordered_map<int, int> placeholder;
        placeholder[my_local_id] = my_node;
          register_placeholder_bindings(placeholder);

        auto children_mx =
            make_children_from_order_with_placeholder_66(order_mx, &placeholder);

        int root_id = new_node(res.Mx, children_mx);
        ROOT_NODE_ID = root_id;

    }

    return true;
}

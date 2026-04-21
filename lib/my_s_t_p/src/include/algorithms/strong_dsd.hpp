#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cstdint>

#include "strong_else_dec.hpp"
#include "node_global.hpp"

// =====================================================
// Debug switch
// =====================================================
inline bool STRONG_DSD_DEBUG_PRINT = false;

// =====================================================
// Pretty print: TT + order (must be paired)
// =====================================================
inline void print_tt_with_order(
    const std::string& title,
    const std::string& tt,
    const std::vector<int>& order,
    int depth = 0)
{
    if (!STRONG_DSD_DEBUG_PRINT) return;

    std::string indent((size_t)depth * 2, ' ');
    std::cout << indent << "📌 " << title << "\n";
    std::cout << indent << "   TT    = " << tt << "\n";
    std::cout << indent << "   order = { ";
    for (int v : order) std::cout << v << " ";
    std::cout << "}\n";
}


// =====================================================
// StrongDsdResult
// =====================================================
struct StrongDsdResult {
    bool found = false;
    size_t L = 0;          // = 2^{|Mx|}
    std::string Mx;        // = block0 + block1, length 2*L
    std::string My;        // length 2^{|My|}
};

// =====================================================
// Split meta
// =====================================================
struct StrongDsdSplit {
    bool found = false;
    StrongDsdResult dsd;

    // 0-based positions in current order
    std::vector<int> mx_pos;
    std::vector<int> my_pos;

    // var IDs in MSB->LSB order (by position)
    std::vector<int> mx_vars_msb2lsb;
    std::vector<int> my_vars_msb2lsb;

    // optional debug info
    std::string block0;
    std::string block1;

    std::string reordered_tt;
};

// =====================================================
// next_combination: comb is strictly increasing indices in [0..m-1]
// =====================================================
inline bool next_combination(std::vector<int>& comb, int m)
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

inline bool is_all_same(const std::string& s, char c)
{
    for (char x : s) if (x != c) return false;
    return true;
}

// =====================================================
// Extract block (Mx subfunction) for a fixed My assignment
// - order positions: 0..n-1 (0 is MSB position)
// - mx_pos_sorted / my_pos_sorted must be sorted ascending by position (MSB->LSB)
// - my_assignment encoded in My's MSB->LSB order
// - mx_index encoded in Mx's MSB->LSB order
// =====================================================
inline std::string extract_block_for_mx(
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
inline void print_candidate_info(
    int depth,
    int k,
    const std::vector<int>& mx_vars_msb2lsb,
    const std::vector<int>& my_vars_msb2lsb)
{
    if (!STRONG_DSD_DEBUG_PRINT) return;

    std::string indent((size_t)depth * 2, ' ');
    std::cout << indent << "🔎 尝试 |Mx|=" << k << " : Mx={ ";
    for (int v : mx_vars_msb2lsb) std::cout << v << " ";
    std::cout << "}  My={ ";
    for (int v : my_vars_msb2lsb) std::cout << v << " ";
    std::cout << "}\n";
}

// =====================================================
// ★ Subset-enumeration Strong DSD
// Search order: |Mx| = ceil(n/2) down to 1
// Combination order: by var_id ascending (a,b,c...)
// =====================================================
inline StrongDsdSplit run_strong_dsd_by_mx_subset(
    const std::string& mf,
    const std::vector<int>& order,
    int depth_for_print = 0)
{
    StrongDsdSplit out;

    int n = (int)order.size();
    if (n <= 1) return out;

    if ((size_t(1) << n) != mf.size()) return out;

    // =====================================================
    // (var_id, pos) sorted by var_id ascending
    // =====================================================
    struct VarPos { int var; int pos; };
    std::vector<VarPos> vp;
    vp.reserve(n);
    for (int pos = 0; pos < n; ++pos)
        vp.push_back({order[pos], pos});

    std::sort(vp.begin(), vp.end(),
        [](const VarPos& a, const VarPos& b) {
            return a.var < b.var;
        });

    int min_k = 1;
    int max_k = n - 2;

    // =====================================================
    // preferred_var = newest (largest id)
    // =====================================================
    int preferred_var = *std::max_element(order.begin(), order.end());
    int preferred_idx = -1;
    for (int i = 0; i < (int)vp.size(); ++i) {
        if (vp[i].var == preferred_var) {
            preferred_idx = i;
            break;
        }
    }

    // =====================================================
    // try_combination: core DSD test (完全不改你原逻辑)
    // =====================================================
    auto try_combination = [&](int k, const std::vector<int>& comb) -> bool {

        // ---- mx_pos / my_pos ----
        std::vector<int> mx_pos;
        for (int idx : comb) mx_pos.push_back(vp[idx].pos);

        std::vector<int> my_pos;
        {
            std::vector<char> is_mx(n, 0);
            for (int p : mx_pos) is_mx[p] = 1;
            for (int p = 0; p < n; ++p)
                if (!is_mx[p]) my_pos.push_back(p);
        }

        std::sort(mx_pos.begin(), mx_pos.end());
        std::sort(my_pos.begin(), my_pos.end());

        std::vector<int> mx_vars_msb2lsb, my_vars_msb2lsb;
        for (int p : mx_pos) mx_vars_msb2lsb.push_back(order[p]);
        for (int p : my_pos) my_vars_msb2lsb.push_back(order[p]);

        print_candidate_info(depth_for_print, k,
                             mx_vars_msb2lsb,
                             my_vars_msb2lsb);

        int m = n - k;
        uint64_t my_count = 1ull << m;
        uint64_t L = 1ull << k;

        if (STRONG_DSD_DEBUG_PRINT) {
            std::string reordered;
            for (uint64_t y = 0; y < my_count; ++y) {
                reordered += extract_block_for_mx(
                    mf, n, mx_pos, my_pos, y);
            }

            std::vector<int> reordered_order;
            reordered_order.insert(
                reordered_order.end(),
                my_vars_msb2lsb.begin(),
                my_vars_msb2lsb.end());
            reordered_order.insert(
                reordered_order.end(),
                mx_vars_msb2lsb.begin(),
                mx_vars_msb2lsb.end());

            print_tt_with_order(
                "候选 split 的重排 TT (My|Mx)",
                reordered,
                reordered_order,
                depth_for_print);
        }

        std::unordered_map<std::string, int> block_index;
        std::vector<std::string> blocks;
        std::string My;

        bool too_many = false;

        for (uint64_t y = 0; y < my_count; ++y) {
            std::string block = extract_block_for_mx(
                mf, n, mx_pos, my_pos, y);

            auto it = block_index.find(block);
            if (it == block_index.end()) {
                if (blocks.size() >= 2) {
                    too_many = true;
                    break;
                }
                int id = (int)blocks.size();
                block_index.emplace(block, id);
                blocks.push_back(block);
                My.push_back(id == 0 ? '1' : '0');
            } else {
                My.push_back(it->second == 0 ? '1' : '0');
            }
        }

        if (!too_many && blocks.size() == 2) {
            out.found = true;
            out.dsd.found = true;
            out.dsd.L = (size_t)L;
            out.dsd.Mx = blocks[0] + blocks[1];
            out.dsd.My = My;

            out.mx_pos = mx_pos;
            out.my_pos = my_pos;
            out.mx_vars_msb2lsb = mx_vars_msb2lsb;
            out.my_vars_msb2lsb = my_vars_msb2lsb;

            if (STRONG_DSD_DEBUG_PRINT) {
                std::string indent(depth_for_print * 2, ' ');
                std::cout << indent << "✅ 命中 Strong DSD split\n";
                print_tt_with_order("当前 split 的 My",
                                    My,
                                    my_vars_msb2lsb,
                                    depth_for_print);
            }
            return true;
        }

        return false;
    };

    // =====================================================
    // PASS 1: preferred_var ∈ Mx（所有 k）
    // =====================================================
    if (preferred_idx >= 0) {
        for (int k = max_k; k >= min_k; --k) {
            std::vector<int> comb(k);
            for (int i = 0; i < k; ++i) comb[i] = i;

            while (true) {
                if (std::find(comb.begin(), comb.end(), preferred_idx)
                    != comb.end()) {
                    if (try_combination(k, comb))
                        return out;
                }
                if (!next_combination(comb, (int)vp.size()))
                    break;
            }
        }
    }

    // =====================================================
    // PASS 2: preferred_var ∉ Mx（兜底）
    // =====================================================
    if (preferred_idx >= 0) {
        for (int k = max_k; k >= min_k; --k) {
            std::vector<int> comb(k);
            for (int i = 0; i < k; ++i) comb[i] = i;

            while (true) {
                if (std::find(comb.begin(), comb.end(), preferred_idx)
                    == comb.end()) {
                    if (try_combination(k, comb))
                        return out;
                }
                if (!next_combination(comb, (int)vp.size()))
                    break;
            }
        }
    }

    return out;
}

// =====================================================
// children
// =====================================================
inline std::vector<int>
make_children_from_order_with_placeholder(
    const std::vector<int>& order,
    const std::unordered_map<int, int>* placeholder_nodes,
    const std::vector<int>* local_to_global)
{
    std::vector<int> children;
    children.reserve(order.size());

    // 关键：按 MSB -> LSB 直接放
    for (int var_id : order) {

        if (placeholder_nodes) {
            auto ph = placeholder_nodes->find(var_id);
            if (ph != placeholder_nodes->end()) {
                children.push_back(ph->second);
                continue;
            }
        }

        
        auto global_it = PLACEHOLDER_BINDINGS.find(var_id);
        if (global_it != PLACEHOLDER_BINDINGS.end())
        {
            children.push_back(global_it->second);
            continue;
        }

        int global_var_id = var_id;
        if (local_to_global &&
            var_id >= 0 &&
            var_id < (int)local_to_global->size())
        {
            int mapped = (*local_to_global)[var_id];
            if (mapped != 0) global_var_id = mapped;
        }

        children.push_back(new_in_node(global_var_id));
    }

    return children;
}

inline int resolve_global_var_id(
    int local_id,
    const std::vector<int>* local_to_global)
{
    if (local_to_global && local_id >= 0 && local_id < (int)local_to_global->size()) {
        int mapped = (*local_to_global)[local_id];
        if (mapped != 0) return mapped;
    }
    return local_id;
}
inline int resolve_var_node_id_no_side_effect(
    int var_id,  // 这是当前 order 里的 var_id（可能是局部变量）
    const std::vector<int>* local_to_global,
    const std::unordered_map<int,int>* placeholder_nodes)
{
    // 1) placeholder_nodes（局部占位）
    if (placeholder_nodes) {
        auto it = placeholder_nodes->find(var_id);
        if (it != placeholder_nodes->end()) return it->second;
    }

    // 2) 全局占位绑定（如果你用过 PLACEHOLDER_BINDINGS）
    auto git = PLACEHOLDER_BINDINGS.find(var_id);
    if (git != PLACEHOLDER_BINDINGS.end()) return git->second;

    // 3) local_to_global 映射
    int global_var = var_id;
    if (local_to_global &&
        var_id >= 0 &&
        var_id < (int)local_to_global->size() &&
        (*local_to_global)[var_id] != 0)
    {
        global_var = (*local_to_global)[var_id];
    }

    // 4) 最终返回 input node
    return new_in_node(global_var);
}

// =====================================================
// ★ Recursive Strong DSD (subset enumeration)
// =====================================================
inline int build_strong_dsd_nodes_impl(
    const std::string& mf,
    const std::vector<int>& order,
    int depth,
    const std::vector<int>* local_to_global,
    const std::unordered_map<int, int>* placeholder_nodes)
{
    // =====================================================
    // Entry
    // =====================================================
    print_tt_with_order("进入 Strong DSD", mf, order, depth);

    const int n = static_cast<int>(order.size());

    // =====================================================
    // (1) Leaf: 真值表最小（≤ 2-input）
    // =====================================================
    if (mf.size() <= 4)
    {
        print_tt_with_order("⏹ Stop (size <= 4)", mf, order, depth);

        // 仅在这里才允许 exact（-e）
        if (ENABLE_ELSE_DEC && n >= 3)
        {
            std::string indent((size_t)depth * 2, ' ');
            std::cout << indent
                      << "⚠️ Leaf exact 2-LUT (n=" << n << ")\n";

            int pivot_node = -1;
            if (!order.empty())
            {
                auto ch = make_children_from_order_with_placeholder(
                    order, placeholder_nodes, local_to_global);
                if (!ch.empty()) pivot_node = ch.front();
            }

            return strong_else_decompose(
                mf,
                order,
                depth,
                pivot_node,
                local_to_global,
                placeholder_nodes,
                build_strong_dsd_nodes_impl);
        }

        // 不开 -e 或 n<=2：直接落地
        auto children = make_children_from_order_with_placeholder(
            order, placeholder_nodes, local_to_global);
        return new_node(mf, children);
    }

    // =====================================================
    // (2) 核心：先尝试 Strong DSD split
    // =====================================================
    StrongDsdSplit split = run_strong_dsd_by_mx_subset(mf, order, depth);

    if (!split.found)
    {
        std::string indent((size_t)depth * 2, ' ');
       // std::cout << indent << "❌ Strong DSD: no valid split\n";

        // =================================================
        // (3) Strong 失败：fallback
        // =================================================
        if (ENABLE_ELSE_DEC && n > 4)
        {
            // 只做一层 Shannon，然后回 strong 主线
            int pivot_node =
                make_children_from_order_with_placeholder(
                    order, placeholder_nodes, local_to_global
                )[0];

            //std::cout << indent
                     // << "⚠️ Fallback: Shannon ONE layer (n=" << n << ")\n";

            return strong_else_decompose(
                mf,
                order,
                depth,
                pivot_node,
                local_to_global,
                placeholder_nodes,
                build_strong_dsd_nodes_impl);
        }

        // n<=4 且 strong 失败：
        // -e 已在 leaf 处理
        auto children = make_children_from_order_with_placeholder(
            order, placeholder_nodes, local_to_global);
        return new_node(mf, children);
    }

    // =====================================================
    // (4) Strong split found：正常递归
    // =====================================================
    const auto& result = split.dsd;

    {
        std::string indent((size_t)depth * 2, ' ');
       // std::cout << indent << "✅ L = " << result.L << "\n";
        //std::cout << indent << "Mx = " << result.Mx << "\n";
        //std::cout << indent << "My = " << result.My << "\n";

        //std::cout << indent << "My 使用变量（MSB->LSB）：{ ";
        //for (int v : split.my_vars_msb2lsb) std::cout << v << " ";
        //std::cout << "}\n";

        //std::cout << indent << "Mx 使用变量（MSB->LSB）：{ ";
        //for (int v : split.mx_vars_msb2lsb) std::cout << v << " ";
        //std::cout << "}\n";
    }

    // -------- recurse My --------
    const std::vector<int>& order_my = split.my_vars_msb2lsb;
    print_tt_with_order("递归进入 My", result.My, order_my, depth);

    int my_id = build_strong_dsd_nodes_impl(
        result.My,
        order_my,
        depth + 1,
        local_to_global,
        placeholder_nodes);

    // -------- recurse Mx --------
    int k = (int)split.mx_vars_msb2lsb.size();
    int my_local_id = k + 1;

    std::vector<int> order_mx;
    order_mx.reserve(k + 1);
    order_mx.push_back(my_local_id);
    for (int i = k; i >= 1; --i) order_mx.push_back(i);

    std::unordered_map<int, int> placeholder_nodes_mx;
    std::vector<int> local_to_global_mx(my_local_id + 1, 0);

    placeholder_nodes_mx[my_local_id] = my_id;

    const std::unordered_map<int, int> empty_ph;
    const auto& parent_ph = placeholder_nodes ? *placeholder_nodes : empty_ph;

    for (int i = 0; i < k; ++i)
    {
        int old_id = split.mx_vars_msb2lsb[i];
        int new_id = k - i;

        auto it = parent_ph.find(old_id);
        if (it != parent_ph.end())
            placeholder_nodes_mx[new_id] = it->second;
        else
            local_to_global_mx[new_id] =
                resolve_global_var_id(old_id, local_to_global);
    }

    print_tt_with_order("递归进入 Mx", result.Mx, order_mx, depth);

    return build_strong_dsd_nodes_impl(
        result.Mx,
        order_mx,
        depth + 1,
        &local_to_global_mx,
        &placeholder_nodes_mx);
}


inline bool strong_is_non_2input_node(int node_id)
{
    if (node_id <= 0) return false;
    const DSDNode* nd = nullptr;
    for (const auto& cand : NODE_LIST)
    {
        if (cand.id == node_id)
        {
            nd = &cand;
            break;
        }
    }
    if (!nd) return false;
    if (nd->func == "in" || nd->func == "0" || nd->func == "1") return false;
    return nd->child.size() > 2;
}

inline void strong_replace_node_everywhere(int old_id, int new_id)
{
    if (old_id == new_id) return;

    for (auto& nd : NODE_LIST)
    {
        for (auto& c : nd.child)
        {
            if (c == old_id) c = new_id;
        }
    }

    if (ROOT_NODE_ID == old_id)
        ROOT_NODE_ID = new_id;
}

inline int strong_refine_non_2input_node(int node_id)
{
    const DSDNode* nd = nullptr;
    for (const auto& cand : NODE_LIST)
    {
        if (cand.id == node_id)
        {
            nd = &cand;
            break;
        }
    }
    if (!nd) return node_id;
    const int n = static_cast<int>(nd->child.size());

    std::vector<int> order;
    order.reserve(n);

    std::unordered_map<int, int> placeholder_nodes;
    placeholder_nodes.reserve(n);

    for (int child_id : nd->child)
    {
        const DSDNode* child = nullptr;
        for (const auto& cand : NODE_LIST)
        {
            if (cand.id == child_id)
            {
                child = &cand;
                break;
            }
        }
        if (child && child->func == "in")
        {
            order.push_back(child->var_id);
        }
        else
        {
            int ph_id = allocate_placeholder_var_id(&placeholder_nodes);
            placeholder_nodes[ph_id] = child_id;
            order.push_back(ph_id);
        }
    }

    const int pivot_node = nd->child.empty() ? -1 : nd->child.front();

    return strong_else_decompose(
        nd->func,
        order,
        /*depth=*/0,
        pivot_node,
        /*local_to_global=*/nullptr,
        &placeholder_nodes,
        build_strong_dsd_nodes_impl);
}

inline void strong_refine_all_non_2input_nodes()
{
    std::vector<int> targets;
    targets.reserve(NODE_LIST.size());

    for (const auto& nd : NODE_LIST)
    {
        if (strong_is_non_2input_node(nd.id))
            targets.push_back(nd.id);
    }

    if (targets.empty())
    {
        std::cout << "✅ Strong DSD: no non-2input nodes to refine\n";
        return;
    }

   // std::cout << "🔧 Strong DSD: refining " << targets.size()
            //  << " non-2input nodes\n";

    for (int node_id : targets)
    {
        if (!strong_is_non_2input_node(node_id)) continue;
        int new_root = strong_refine_non_2input_node(node_id);
        strong_replace_node_everywhere(node_id, new_root);
    }
}
inline bool is_need_post_decompose(const DSDNode& nd)
{
    // 基本节点不处理
    if (nd.func == "in" || nd.func == "0" || nd.func == "1")
        return false;

    // 只关心 >2-input
    return nd.child.size() > 2;
}

inline void post_decompose_all_large_nodes_fixpoint()
{
   // std::cout << "🔧 Post-decompose: start fixpoint refinement\n";

    bool changed = true;
    int round = 0;

    while (changed)
    {
        changed = false;
        ++round;

        //std::cout << "🔁 Post-decompose round " << round << "\n";

        // ⚠️ 每一轮都重新扫描整个 NODE_LIST
        for (size_t i = 0; i < NODE_LIST.size(); ++i)
        {
            const DSDNode& nd = NODE_LIST[i];

            if (!is_need_post_decompose(nd))
                continue;

            int old_id = nd.id;

            // =====================================================
            // 🔑 关键：跳过已经“脱网”的节点
            // =====================================================
            bool referenced = (ROOT_NODE_ID == old_id);

            if (!referenced)
            {
                for (const auto& n2 : NODE_LIST)
                {
                    for (int c : n2.child)
                    {
                        if (c == old_id)
                        {
                            referenced = true;
                            break;
                        }
                    }
                    if (referenced) break;
                }
            }

            if (!referenced)
                continue;
            // =====================================================

            //std::cout << "  🔍 Found >2-input node: id=" << old_id
                  //  << " fanin=" << nd.child.size()
                  // << " func=" << nd.func << "\n";

            int new_id = strong_refine_non_2input_node(old_id);

            if (new_id != old_id)
            {
                //std::cout << "  ✂️ Refined node " << old_id
                      //  << " -> " << new_id << "\n";

                strong_replace_node_everywhere(old_id, new_id);

                // =====================================================
                // 🔥 关键：从 NODE_LIST 中物理删除 old 节点
                // =====================================================
                for (size_t k = 0; k < NODE_LIST.size(); ++k)
                {
                    if (NODE_LIST[k].id == old_id)
                    {
                        NODE_LIST.erase(NODE_LIST.begin() + k);
                        break;
                    }
                }
                // =====================================================

                changed = true;
                break;
            }


        }
    }

    //std::cout << "✅ Post-decompose finished: no >2-input nodes left\n";
}

inline int build_strong_dsd_nodes(
    const std::string& mf,
    const std::vector<int>& order,
    int depth = 0)
{
    int root_id =
        build_strong_dsd_nodes_impl(mf, order, depth, nullptr, nullptr);

    ROOT_NODE_ID = root_id;

    // ⭐ Strong DSD 完全结束后，再做后处理
    if (ENABLE_ELSE_DEC)
    {
        post_decompose_all_large_nodes_fixpoint();
        root_id = ROOT_NODE_ID;
    }

    return root_id;
}

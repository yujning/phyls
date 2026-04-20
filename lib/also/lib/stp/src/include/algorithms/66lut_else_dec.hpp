#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "node_global.hpp"
#include "66lut_dsd.hpp"
#include "strong_dsd.hpp"

// =====================================================
// Helper: reuse 66-LUT DSD (|My|<=6, |Mx|<=5) to build a small chain.
// Supports placeholder nodes so we can embed the result into a larger DAG.
// =====================================================
inline int build_66lut_chain(
    const std::string& tt,
    const std::vector<int>& order,
    const std::unordered_map<int, int>* placeholder_nodes)
{
    if ((size_t(1) << order.size()) != tt.size())
        return -1;

    // Base case: direct 6-LUT
    if (order.size() <= 6)
    {
        auto children =
            make_children_from_order_with_placeholder_66(order, placeholder_nodes);
        return new_node(tt, children);
    }

    // Try a 66-LUT strong DSD split first.
    auto res = run_66lut_dsd_by_mx_subset(tt, order, /*depth_for_print=*/0);
    if (!res.found)
        return -1;

    // Build MY (always <= 6)
    int my_node = build_66lut_chain(res.My, res.my_vars_msb2lsb, nullptr);
    if (my_node < 0)
        return -1;

    // Build MX with MY placeholder at the MSB position
       int placeholder_id = allocate_placeholder_var_id(placeholder_nodes);


    std::unordered_set<int> my_var_set(
        res.my_vars_msb2lsb.begin(), res.my_vars_msb2lsb.end());

    std::vector<int> order_mx;
    order_mx.reserve(res.mx_vars_msb2lsb.size() + 1);
    order_mx.push_back(placeholder_id);
    for (int v : res.mx_vars_msb2lsb)
        if (my_var_set.count(v) == 0)
            order_mx.push_back(v);

    std::unordered_map<int, int> merged_placeholder;
    merged_placeholder[placeholder_id] = my_node;
    if (placeholder_nodes)
        merged_placeholder.insert(placeholder_nodes->begin(), placeholder_nodes->end());
    register_placeholder_bindings(merged_placeholder);
    auto children_mx =
        make_children_from_order_with_placeholder_66(order_mx, &merged_placeholder);
    return new_node(res.Mx, children_mx);
}

// =====================================================
// 66-LUT fallback that leverages "dsd -f -s" first-level split
// 1) Perform one strong DSD split (no |Mx|/|My| constraint).
// 2) If either side has <= 6 vars, stop there; otherwise reuse 66-LUT flow
//    to further decompose the larger MX subfunction.
// =====================================================
inline bool run_66lut_else_dec_and_build_dag(const TT& root_tt)
{
    const int n = (int)root_tt.order.size();
    if ((size_t(1) << n) != root_tt.f01.size())
        return false;

    RESET_NODE_GLOBAL();
    ORIGINAL_VAR_COUNT = n;

    auto split = run_strong_dsd_by_mx_subset(root_tt.f01, root_tt.order, 0);
    if (!split.found)
        return false;

    const auto& res = split.dsd;

    // Build MY (stop immediately if it already fits in a 6-LUT)
    int my_node = build_66lut_chain(res.My, split.my_vars_msb2lsb, nullptr);
    if (my_node < 0)
        return false;

    // Prepare MX input order with MY as MSB placeholder
       const int placeholder_id = allocate_placeholder_var_id(nullptr);
    std::unordered_set<int> my_var_set(
        split.my_vars_msb2lsb.begin(), split.my_vars_msb2lsb.end());

    std::vector<int> order_mx;
    order_mx.reserve(split.mx_vars_msb2lsb.size() + 1);
    order_mx.push_back(placeholder_id);
    for (int v : split.mx_vars_msb2lsb)
        if (my_var_set.count(v) == 0)
            order_mx.push_back(v);

    std::unordered_map<int, int> placeholder;
    placeholder[placeholder_id] = my_node;
        register_placeholder_bindings(placeholder);
    int root_node = -1;
    if ((int)order_mx.size() <= 6)
    {
        // MX also small enough: stop after first split.
        auto children =
            make_children_from_order_with_placeholder_66(order_mx, &placeholder);
        root_node = new_node(res.Mx, children);
    }
    else
    {
        // MX is still large: recursively decompose with the 66-LUT flow.
        root_node = build_66lut_chain(res.Mx, order_mx, &placeholder);
    }

    return root_node >= 0;
}
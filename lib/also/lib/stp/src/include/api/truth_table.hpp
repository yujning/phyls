#pragma once

#include <sstream>
#include <string>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>

#include <optional>
#include <numeric>

#include "../algorithms/bi_decomposition.hpp"



namespace stp
{
struct bidecomposition_nodes
{
    int root_id{0};
    std::vector<DSDNode> nodes;
    std::vector<int> variable_order;
};

/**
 * Run DSD decomposition directly from a binary truth table string.
 *
 * @param binary_truth_table Truth table encoded as a binary string whose
 *                           length is a power of two.
 * @return true when the input is valid and the run succeeds.
 */
inline bool run_dsd_from_binary(const std::string& binary_truth_table)
{
    return run_dsd_recursive(binary_truth_table);
}

/**
 * Run DSD decomposition from a hexadecimal truth table string.
 *
 * @param hex_truth_table Hexadecimal string representing a truth table.
 * @return true when the input is valid and the run succeeds.
 */
inline bool run_dsd_from_hex(const std::string& hex_truth_table)
{
    auto hex = hex_truth_table;
    if (hex.rfind("0x", 0) == 0 || hex.rfind("0X", 0) == 0)
    {
        hex = hex.substr(2);
    }

    const auto bit_count = static_cast<unsigned>(hex.size() * 4);

    unsigned num_vars = 0;
    while ((1u << num_vars) < bit_count) ++num_vars;
    if ((1u << num_vars) != bit_count)
    {
        return false;
    }

    kitty::dynamic_truth_table tt(num_vars);
    kitty::create_from_hex_string(tt, hex);

    std::ostringstream oss;
    kitty::print_binary(tt, oss);
    return run_dsd_recursive(oss.str());
}

/**
 * Run bi-decomposition directly from a binary truth table string.
 *
 * @param binary_truth_table Truth table encoded as a binary string whose
 *                           length is a power of two.
 * @return true when the input is valid and the run succeeds.
 */
inline bool run_bidecomposition_from_binary(const std::string& binary_truth_table)
{
    return run_bi_decomp_recursive(binary_truth_table);
}

/**
 * Run bi-decomposition directly from a kitty truth table.
 *
 * @param tt kitty dynamic truth table.
 * @return true when the input is valid and the run succeeds.
 */
inline bool run_bidecomposition(const kitty::dynamic_truth_table& tt)
{
    std::ostringstream oss;
    kitty::print_binary(tt, oss);
    return run_bi_decomp_recursive(oss.str());
}

inline std::optional<bidecomposition_nodes> capture_bidecomposition(const kitty::dynamic_truth_table& tt)
{
    std::ostringstream oss;
    kitty::print_binary(tt, oss);

    RESET_NODE_GLOBAL();
    const auto prev_output = BD_MINIMAL_OUTPUT;
    BD_MINIMAL_OUTPUT = true;
    ENABLE_ELSE_DEC = true;

    const auto num_vars = tt.num_vars();
    ORIGINAL_VAR_COUNT = static_cast<int>(num_vars);

    TT root;
    root.f01 = oss.str();
    root.order.resize(num_vars);
    std::iota(root.order.begin(), root.order.end(), 1);

    for (int v = 1; v <= ORIGINAL_VAR_COUNT; ++v)
    {
        new_in_node(v);
    }

    auto root_shrunk = shrink_to_support(root);
    const auto root_id = bi_decomp_recursive(root_shrunk, 0);

    bidecomposition_nodes result{ root_id, NODE_LIST, FINAL_VAR_ORDER };

    BD_MINIMAL_OUTPUT = prev_output;

    if (root_id <= 0)
    {
        return std::nullopt;
    }

    return result;
}
/**
 * Run DSD directly from a kitty truth table.
 *
 * @param tt kitty dynamic truth table.
 * @return true when the input is valid and the run succeeds.
 */
inline bool run_dsd(const kitty::dynamic_truth_table& tt)
{
    std::ostringstream oss;
    kitty::print_binary(tt, oss);
    return run_dsd_recursive(oss.str());
}

} // namespace stp
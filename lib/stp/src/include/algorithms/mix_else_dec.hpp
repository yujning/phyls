#pragma once

// ⭐ 防止 std::is_trivial 污染 percy
#ifdef is_trivial
#undef is_trivial
#endif

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "node_global.hpp"
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>

#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/networks/klut.hpp>

inline int mix_resolve_var_node_id(
    int var_id,
    const std::vector<int>* local_to_global,
    const std::unordered_map<int, int>* placeholder_nodes)
{
  if (placeholder_nodes)
  {
    auto it = placeholder_nodes->find(var_id);
    if (it != placeholder_nodes->end())
      return it->second;
  }

  auto git = PLACEHOLDER_BINDINGS.find(var_id);
  if (git != PLACEHOLDER_BINDINGS.end())
    return git->second;

  int global_var = var_id;
  if (local_to_global &&
      var_id >= 0 &&
      var_id < (int)local_to_global->size() &&
      (*local_to_global)[var_id] != 0)
  {
    global_var = (*local_to_global)[var_id];
  }

  return new_in_node(global_var);
}

inline std::vector<int> mix_make_children_kitty_order(
    const std::vector<int>& order_msb2lsb,
    const std::vector<int>* local_to_global,
    const std::unordered_map<int, int>* placeholder_nodes)
{
  std::vector<int> ch;
  ch.reserve(order_msb2lsb.size());

  for (auto it = order_msb2lsb.rbegin(); it != order_msb2lsb.rend(); ++it)
    ch.push_back(mix_resolve_var_node_id(*it, local_to_global, placeholder_nodes));

  return ch;
}

inline int mix_exact_refine_2lut(
    const std::string& f01,
    const std::vector<int>& order,
    int depth,
    const std::vector<int>* local_to_global,
    const std::unordered_map<int, int>* placeholder_nodes)
{
  const int n = static_cast<int>(order.size());
  std::string indent(static_cast<size_t>(depth) * 2, ' ');

  std::cout << indent << "⚠️ Mix EXACT refine (n=" << n << ")\n";
  std::cout << indent << "f=" << f01 << "\n";

  for (char c : f01)
    if (c != '0' && c != '1')
      throw std::runtime_error("mix_exact_refine_2lut: f01 contains non-binary char");

  kitty::dynamic_truth_table tt(n);
  kitty::create_from_binary_string(tt, f01);
  std::cout << indent << "[DEBUG] kitty hex = " << kitty::to_hex(tt) << "\n";

  mockturtle::klut_network klut;
  std::vector<mockturtle::klut_network::signal> pis;
  pis.reserve(n);
  for (int i = 0; i < n; ++i)
    pis.push_back(klut.create_pi());

  mockturtle::exact_resynthesis<mockturtle::klut_network> resyn(2);
  resyn(klut, tt, pis.begin(), pis.end(),
        [&](auto const& s) { klut.create_po(s); });

  std::cout << indent << "Exact 2-LUT count = " << klut.num_gates() << "\n";

  std::unordered_map<mockturtle::klut_network::node, int> node_map;
  auto orig_children = mix_make_children_kitty_order(
      order, local_to_global, placeholder_nodes);

  klut.foreach_pi([&](auto const& n_pi, auto index) {
    if (index < orig_children.size())
      node_map[n_pi] = orig_children[index];
  });

  klut.foreach_gate([&](auto const& n_gate) {
    std::vector<int> childs;
    childs.reserve(klut.fanin_size(n_gate));

    klut.foreach_fanin(n_gate, [&](auto const& f_handle) {
      const auto f_node = klut.get_node(f_handle);

      if (klut.is_constant(f_node))
      {
        const bool is_one = klut.is_complemented(f_handle) ^ klut.constant_value(f_node);
        childs.push_back(new_node(is_one ? "1" : "0", {}));
      }
      else
      {
        auto cid = node_map.at(f_node);
        if (klut.is_complemented(f_handle))
          cid = new_node("01", { cid });
        childs.push_back(cid);
      }
    });

    std::reverse(childs.begin(), childs.end());

    auto func = klut.node_function(n_gate);
    std::string func_bin = kitty::to_binary(func);

    node_map[n_gate] = new_node(func_bin, childs);
  });

  auto po_sig = klut.po_at(0);
  auto po_node = klut.get_node(po_sig);

  if (klut.is_constant(po_node))
  {
    bool val = klut.constant_value(po_node);
    if (klut.is_complemented(po_sig))
      val = !val;
    return new_node(val ? "1" : "0", {});
  }

  auto root_id = node_map.at(po_node);
  if (klut.is_complemented(po_sig))
    root_id = new_node("01", { root_id });

  return root_id;
}

inline int mix_else_decompose(
    const TT& f,
    int depth,
    const std::vector<int>* local_to_global,
    const std::unordered_map<int, int>* placeholder_nodes,
    const std::function<int(
        const TT&,
        int,
        const std::vector<int>*,
        const std::unordered_map<int, int>*)>& mix_rec)
{
  const uint32_t bits = static_cast<uint32_t>(f.f01.size());
  const uint32_t n = static_cast<uint32_t>(std::log2(bits));

  if ((1u << n) != bits)
    throw std::runtime_error("mix_else_decompose: f01 length is not power-of-two");

  if (n > 4u)
  {
    std::string indent(static_cast<size_t>(depth) * 2, ' ');
    std::cout << indent << "⚠️ Mix fallback: Shannon decomposition (n=" << n << ")\n";
    std::cout << indent << "f=" << f.f01 << "\n";

    auto orig_children =
        mix_make_children_kitty_order(f.order, local_to_global, placeholder_nodes);
    if (orig_children.empty())
      throw std::runtime_error("mix_else_decompose: no children for Shannon split");

    const auto pivot_node = orig_children.back();
    const auto pivot_var = f.order.empty() ? -1 : f.order.front();

    const auto half_bits = bits / 2u;
    TT f_pos{ f.f01.substr(0, half_bits), {} };
    TT f_neg{ f.f01.substr(half_bits), {} };

    if (!f.order.empty())
    {
      f_pos.order.assign(f.order.begin() + 1, f.order.end());
      f_neg.order.assign(f.order.begin() + 1, f.order.end());
    }

    std::cout << indent << "  split depth " << depth
              << " pivot var=" << pivot_var
              << " pos f=" << f_pos.f01 << "\n";
    std::cout << indent << "  split depth " << depth
              << " pivot var=" << pivot_var
              << " neg f=" << f_neg.f01 << "\n" << std::flush;

    const auto pos_node =
        mix_rec(f_pos, depth + 1, local_to_global, placeholder_nodes);
    const auto neg_node =
        mix_rec(f_neg, depth + 1, local_to_global, placeholder_nodes);

    const auto pos_term = new_node("1000", { pivot_node, pos_node });
    const auto neg_term = new_node("0010", { pivot_node, neg_node });

    return new_node("1110", { pos_term, neg_term });
  }

  return mix_exact_refine_2lut(f.f01, f.order, depth, local_to_global, placeholder_nodes);
}
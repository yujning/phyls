#pragma once

// ⭐ 防止 std::is_trivial 污染 percy
#ifdef is_trivial
#undef is_trivial
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>


#include "node_global.hpp"
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
// 前向声明（定义在 stp_dsd.hpp 中）
struct TT;
static int build_small_tree(const TT& t);
static int dsd_factor(const TT& f, int depth);

inline int dsd_else_decompose(
    const TT& f,
    int depth)
{
  const uint32_t bits = static_cast<uint32_t>(f.f01.size());
  const uint32_t n = static_cast<uint32_t>(std::log2(bits));

  if ((1u << n) != bits)
    throw std::runtime_error("dsd_else_decompose: f01 length is not power-of-two");

  auto orig_children = make_children_from_order(f);

  if (n > 4u)
  {
    std::cout << "⚠️ depth " << depth
              << ": fallback Shannon decomposition (n=" << n << ")\n";
    std::cout << "f=" << f.f01 << "\n";
    std::cout << "输入节点 = { ";
    for (size_t i = 0; i < orig_children.size(); ++i)
    {
      std::cout << orig_children[i];
      if (i + 1 < orig_children.size())
        std::cout << ",";
    }
    std::cout << " } 变量 = { ";
    for (size_t i = 0; i < f.order.size(); ++i)
    {
      std::cout << f.order[i];
      if (i + 1 < f.order.size())
        std::cout << ",";
    }
    std::cout << " }\n";

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

    std::cout << "  split depth " << depth
              << " pivot var=" << pivot_var
              << " pos f=" << f_pos.f01 << "\n";
    std::cout << "  split depth " << depth
              << " pivot var=" << pivot_var
              << " neg f=" << f_neg.f01 << "\n" << std::flush;

    const auto pos_node = dsd_factor(f_pos, depth + 1);
    const auto neg_node = dsd_factor(f_neg, depth + 1);

    const auto pos_term = new_node("1000", { pivot_node, pos_node });
    const auto neg_term = new_node("0010", { pivot_node, neg_node });

    return new_node("1110", { pos_term, neg_term });
  }
  
  for (char c : f.f01)
    if (c != '0' && c != '1')
      throw std::runtime_error("dsd_else_decompose: f01 contains non-binary char");

  std::cout << "⚠️ depth " << depth
            << ": EXACT 2-LUT refine (n=" << n << ")\n";
  std::cout << "f=" << f.f01 << "\n";

  kitty::dynamic_truth_table tt(n);
  kitty::create_from_binary_string(tt, f.f01);
  std::cout << "[DEBUG] kitty hex = " << kitty::to_hex(tt) << "\n";

  mockturtle::klut_network klut;
  std::vector<mockturtle::klut_network::signal> pis;
  pis.reserve(n);

  for (uint32_t i = 0; i < n; ++i)
    pis.push_back(klut.create_pi());

  if (pis.size() != n)
    throw std::runtime_error("dsd_else_decompose: PI creation failed");


mockturtle::exact_resynthesis<mockturtle::klut_network> resyn(2);
  resyn(klut, tt, pis.begin(), pis.end(),
        [&](auto const& s) {
          klut.create_po(s);
        });

  std::cout << "Exact 2-LUT count = " << klut.num_gates() << "\n";

  std::unordered_map<mockturtle::klut_network::node, int> node_map;

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
    if (klut.is_complemented(po_sig)) val = !val;
    return new_node(val ? "1" : "0", {});
  }

  auto root_id = node_map.at(po_node);
  if (klut.is_complemented(po_sig))
  root_id = new_node("01", { root_id });
  return root_id;
}
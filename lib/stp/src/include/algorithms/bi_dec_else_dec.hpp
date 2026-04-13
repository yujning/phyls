#pragma once

// ⭐ 防止 std::is_trivial 污染 percy
#ifdef is_trivial
#undef is_trivial
#endif

#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <algorithm>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>
#include <kitty/constructors.hpp>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include "stp_dsd.hpp"
#include "node_global.hpp"   // new_node / new_in_node

// 你已有的接口（外部实现）
int new_node(const std::string&, const std::vector<int>&);

// 如果你系统里有常量节点，建议你接上（可选）
// int const0_node();
// int const1_node();

// 如果你已经有 build_small_tree(f) 就保留，否则你自己实现一个 fallback
int build_small_tree(const TT& f);

static int bi_decomp_recursive(const TT& f, int depth);



inline int else_decompose(
    const TT& f,
    const std::vector<int>& orig_children,
    int depth
)
{
  // f.f01 是二进制字符串，长度应该是 2^n
  const uint32_t bits = static_cast<uint32_t>(f.f01.size());
  const uint32_t n = static_cast<uint32_t>(std::log2(bits));

  if ((1u << n) != bits)
    throw std::runtime_error("else_decompose: f01 length is not power-of-two");

  if (n > 4u)
  {
    std::cout << "⚠️ depth " << depth
              << ": Shannon decomposition (n=" << n << ")\n";
     std::cout << "f=" << f.f01 << "\n";

    if (orig_children.empty())
      throw std::runtime_error("else_decompose: no children for Shannon split");

    // 香农分解先以高位变量为轴（orig_children 最后一个变量）
    const auto pivot = orig_children.back();

    const auto half_bits = bits / 2u;
    TT f_pos{ f.f01.substr(0, half_bits), {} };
    TT f_neg{ f.f01.substr(half_bits), {} };

    if (!f.order.empty())
    {
      f_pos.order.assign(f.order.begin() + 1, f.order.end());
      f_neg.order.assign(f.order.begin() + 1, f.order.end());
    }

        std::cout << "  split depth " << depth << " pos f=" << f_pos.f01 << "\n";
    std::cout << "  split depth " << depth << " neg f=" << f_neg.f01 << "\n"
              << std::flush;

    const auto pos_node = bi_decomp_recursive(f_pos, depth + 1);
    const auto neg_node = bi_decomp_recursive(f_neg, depth + 1);

    const auto pos_term = new_node("1000", { pivot, pos_node });
    const auto neg_term = new_node("0010", { pivot, neg_node });

    return new_node("1110", { pos_term, neg_term });
  }

  // 可选：检查字符串只含 0/1
  for (char c : f.f01)
    if (c != '0' && c != '1')
      throw std::runtime_error("else_decompose: f01 contains non-binary char");

  std::cout << "⚠️ depth " << depth
            << ": EXACT 2-LUT refine (n=" << n << ")\n";
  std::cout << "f=" << f.f01 << "\n";

  // 1) binary string -> kitty TT
  kitty::dynamic_truth_table tt(n);
  kitty::create_from_binary_string(tt, f.f01);
  std::cout << "[DEBUG] kitty hex = " << kitty::to_hex(tt) << "\n";

  // 2) build klut and create EXACTLY n PIs
  mockturtle::klut_network klut;
  std::vector<mockturtle::klut_network::signal> pis;
  pis.reserve(n);

  for (uint32_t i = 0; i < n; ++i)
    pis.push_back(klut.create_pi());

  // 防御：必须保证输入数量对
  if (pis.size() != n)
    throw std::runtime_error("else_decompose: PI creation failed");

  // 3) exact 2-LUT
  mockturtle::exact_resynthesis<mockturtle::klut_network> resyn(2);
  resyn(klut, tt, pis.begin(), pis.end(),
        [&](auto const& s) {
          klut.create_po(s);
        });

  std::cout << "Exact 2-LUT count = " << klut.num_gates() << "\n";

  // 3) 把 klut 网络翻译成全局节点（遵循 bd 打印格式）
  // 3) 把 klut 网络翻译成全局节点（遵循 bd 打印格式）
  std::unordered_map<mockturtle::klut_network::node, int> node_map;

  // 输入节点：按原始顺序绑定到 orig_children
  klut.foreach_pi([&](auto const& n, auto index) {
    if (index < orig_children.size())
      node_map[n] = orig_children[index];
  });

  // 中间 gate
  klut.foreach_gate([&](auto const& n) {
    std::vector<int> childs;
    childs.reserve(klut.fanin_size(n));

    klut.foreach_fanin(n, [&](auto const& f) {
      const auto f_node = klut.get_node(f);

      if (klut.is_constant(f_node))
      {
        // 只在真的需要时才建常量节点，避免额外打印
        const bool is_one = klut.is_complemented(f) ^ klut.constant_value(f_node);
        childs.push_back(new_node(is_one ? "1" : "0", {}));
      }
      else
      {
        auto cid = node_map.at(f_node);
        if (klut.is_complemented(f))
          cid = new_node("01", { cid });
        childs.push_back(cid);
        
      }
    });

  std::reverse(childs.begin(), childs.end());

  auto func = klut.node_function(n);
  std::string func_bin = kitty::to_binary(func);

  node_map[n] = new_node(func_bin, childs);

  });

  // 输出（只取第一个 PO）
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
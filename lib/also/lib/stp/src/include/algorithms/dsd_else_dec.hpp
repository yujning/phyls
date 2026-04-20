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
#include <limits>

#include "node_global.hpp"
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/operations.hpp>


// 前向声明（定义在 stp_dsd.hpp 中）
struct TT;

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


  
  for (char c : f.f01)
    if (c != '0' && c != '1')
      throw std::runtime_error("dsd_else_decompose: f01 contains non-binary char");

  std::cout << "⚠️ depth " << depth
     << ": Shannon decomposition (n=" << n << ")\n";
  std::cout << "f=" << f.f01 << "\n";

  kitty::dynamic_truth_table tt(n);
  kitty::create_from_binary_string(tt, f.f01);

  std::vector<uint32_t> support;
  support.reserve(n);

  for (uint32_t i = 0; i < n; ++i)
   {
    if (kitty::has_var(tt, i))
      support.push_back(i);
  }

  uint32_t var_min = support.empty() ? 0u : support.front();
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
        s_co0 += 1;
    }
    for (uint32_t s_sd = 0; s_sd < co11.num_vars(); ++s_sd)
    {
      if (kitty::has_var(co11, s_sd))
        s_co1 += 1;
    }
    const int add = s_co0 + s_co1;
    if (add < best_score)
    {
      best_score = add;
      var_min = var;
    }
  }

   const auto order_index = static_cast<size_t>(n - 1u - var_min);
  const auto pivot_node = orig_children.at(var_min);
  const auto pivot_var = f.order.empty() ? -1 : f.order.at(order_index);

 auto build_cofactor_tt = [&](bool value) {
    const uint32_t new_vars = n - 1u;
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
        old |= (1ull << var_min);

    if (kitty::get_bit(tt, old))
        kitty::set_bit(new_tt, x);
    }
    TT out;
    out.f01.resize(limit);
    for (uint64_t i = 0; i < limit; ++i)
      out.f01[i] = kitty::get_bit(new_tt, i) ? '1' : '0';
   out.order = f.order;
    if (!out.order.empty())
      out.order.erase(out.order.begin() + static_cast<long>(order_index));
    return out;
  };

     TT f_pos = build_cofactor_tt(false);
  TT f_neg = build_cofactor_tt(true);

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
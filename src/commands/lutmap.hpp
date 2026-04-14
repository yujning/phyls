/* phyLS: powerful heightened yielded Logic Synthesis
 * Copyright (C) 2022 */

/**
 * @file lutmap.hpp
 *
 * @brief performs FPGA technology mapping of the network
 *
 * @author Homyoung
 * @since  2022/12/21
 */

#ifndef LUTMAP_HPP
#define LUTMAP_HPP

#include <mockturtle/generators/arithmetic.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/networks/klut.hpp>
#include "../core/lut_mapper.hpp"

#include <optional>
#include <unordered_map>
#include <numeric>
#include <tuple>
#include <map>

#include <mockturtle/views/topo_view.hpp>
struct DSDNode {
  int id;
  std::string func;
  std::vector<int> child;
  int var_id = -1;
};
#include "../../lib/my_s_t_p/src/include/algorithms/node_global.hpp"

inline int new_node(const std::string& func, const std::vector<int>& child) {
  auto key = std::make_tuple(func, child);
  if (NODE_HASH.count(key)) return NODE_HASH[key];
  int id = NODE_ID++;
  NODE_LIST.push_back({id, func, child, -1});
  NODE_HASH[key] = id;
  return id;
}

inline int new_in_node(int var_id) {
  if (INPUT_NODE_CACHE.count(var_id)) return INPUT_NODE_CACHE[var_id];
  int id = NODE_ID++;
  NODE_LIST.push_back({id, "in", {}, var_id});
  INPUT_NODE_CACHE[var_id] = id;
  return id;
}

inline std::vector<int> make_children_from_order(const TT& t) {
  std::vector<int> ch;
  ch.reserve(t.order.size());
  for (auto it = t.order.rbegin(); it != t.order.rend(); ++it) {
    ch.push_back(new_in_node(*it));
  }
  return ch;
}
#include "../../lib/my_s_t_p/src/include/algorithms/strong_dsd.hpp"


using namespace std;
using namespace mockturtle;

namespace alice {

class lutmap_command : public command {
 public:
  explicit lutmap_command(const environment::ptr& env)
      : command(env, "FPGA technology mapping of the network [default = AIG]") {
    add_flag("--mig, -m", "FPGA technology mapping for MIG");
    add_flag("--xag, -g", "FPGA technology mapping for XAG");
    add_flag("--xmg, -x", "FPGA technology mapping for XMG");
    add_flag("--klut, -k", "FPGA technology mapping for k-LUT");
    add_option("--cut_size, -s", cut_size,
               "Maximum number of leaves for a cut [default = 6]");
    add_option(
        "--cut_limit, -l", cut_limit,
        "the input Maximum number of cuts for a node name [default = 25]");
    add_option("--relax_required, -r", relax_required,
               "delay relaxation ratio (%) [default = 0]");
    add_flag("--area, -a", "toggles area-oriented mapping [default = false]");
    add_flag("--recompute_cuts, -c",
             "recompute cuts at each step [default = true]");
    add_flag("--edge, -e", "Use edge count reduction [default = true]");
    add_flag("--cost_function, -f",
             "LUT map with cost function [default = false]");
    add_flag("--dominated_cuts, -d",
             "Remove the cuts that are contained in others [default = true]");
    add_option("--output, -o", filename, "the bench filename");
    add_flag("--verbose, -v", "print the information");
      add_flag("--stp", "decompose mapped kLUT with strong_dsd+else_dec and remap");
  }

 protected:
 static kitty::dynamic_truth_table tt_from_binary(std::string const& bits) {
    uint32_t num_vars = 0u;
    while ((1u << num_vars) < bits.size()) {
      ++num_vars;
    }
    kitty::dynamic_truth_table tt(num_vars);
    for (uint32_t i = 0u; i < bits.size(); ++i) {
      if (bits[i] == '1') {
        kitty::set_bit(tt, i);
      }
    }
    return tt;
  }

  static std::string tt_to_binary01(kitty::dynamic_truth_table const& tt) {
    std::string bits;
    bits.reserve(1u << tt.num_vars());
    for (uint32_t i = 0u; i < (1u << tt.num_vars()); ++i) {
      bits.push_back(kitty::get_bit(tt, i) ? '1' : '0');
    }
    return bits;
  }

static std::optional<klut_network> stp_decompose_klut_network(
    klut_network const& ntk)
{
  klut_network decomposed;
  topo_view topo{ntk};
  std::unordered_map<uint32_t, klut_network::signal> old2new;

  old2new[topo.node_to_index(topo.get_constant(false))] =
      decomposed.get_constant(false);
  if (topo.get_constant(false) != topo.get_constant(true)) {
    old2new[topo.node_to_index(topo.get_constant(true))] =
        decomposed.get_constant(true);
  }

  topo.foreach_pi([&](auto const& n) {
    old2new[topo.node_to_index(n)] = decomposed.create_pi();
  });

  topo.foreach_gate([&](auto const& n) {
    std::vector<klut_network::signal> fanins;
    topo.foreach_fanin(n, [&](auto const& f) {
      auto idx = topo.node_to_index(topo.get_node(f));
      auto it = old2new.find(idx);
      if (it == old2new.end()) {
        throw std::runtime_error("stp_decompose_klut_network: missing mapped fanin");
      }
      auto sig = it->second;
      fanins.push_back(topo.is_complemented(f) ? !sig : sig);
    });

    auto const node_tt = topo.node_function(n);
    const uint32_t fanin_size = static_cast<uint32_t>(fanins.size());

    // 小节点不做 STP，直接保留
    if (fanin_size <= 2u) {
      old2new[topo.node_to_index(n)] = decomposed.create_node(fanins, node_tt);
      return;
    }

    try {
      RESET_NODE_GLOBAL();
      ENABLE_ELSE_DEC = true;
      STRONG_DSD_DEBUG_PRINT = false;

      std::vector<int> order(fanin_size);
      std::iota(order.begin(), order.end(), 1);

      const auto root_id =
          build_strong_dsd_nodes(tt_to_binary01(node_tt), order);

      if (root_id <= 0) {
        old2new[topo.node_to_index(n)] =
            decomposed.create_node(fanins, node_tt);
        return;
      }

      // 建一个 id -> DSDNode* 的索引表，方便递归查找
      std::unordered_map<int, const DSDNode*> id2node;
      id2node.reserve(NODE_LIST.size());
      for (auto const& dsd_node : NODE_LIST) {
        id2node.emplace(dsd_node.id, &dsd_node);
      }

      // 递归缓存：DSD node id -> decomposed signal
      std::unordered_map<int, klut_network::signal> memo;
      memo.reserve(NODE_LIST.size());

      std::function<klut_network::signal(int)> build_signal =
          [&](int node_id) -> klut_network::signal {
        auto it_memo = memo.find(node_id);
        if (it_memo != memo.end()) {
          return it_memo->second;
        }

        auto it_node = id2node.find(node_id);
        if (it_node == id2node.end()) {
          throw std::runtime_error("stp_decompose_klut_network: unknown DSD node id");
        }

        const auto& dsd_node = *(it_node->second);

        // 输入节点：映射回当前 LUT 的 fanins
        if (dsd_node.func == "in") {
          if (dsd_node.var_id <= 0 ||
              dsd_node.var_id > static_cast<int>(fanins.size())) {
            throw std::runtime_error("stp_decompose_klut_network: input var_id out of range");
          }
          auto sig = fanins[static_cast<uint32_t>(dsd_node.var_id - 1)];
          memo[node_id] = sig;
          return sig;
        }

        // 常量
        if (dsd_node.func == "0" || dsd_node.func == "1") {
          auto sig = decomposed.get_constant(dsd_node.func == "1");
          memo[node_id] = sig;
          return sig;
        }

        // 普通内部节点：先递归构建 children，再 create_node
        std::vector<klut_network::signal> children;
        children.reserve(dsd_node.child.size());
        for (auto child_id : dsd_node.child) {
          children.push_back(build_signal(child_id));
        }

        auto sig =
            decomposed.create_node(children, tt_from_binary(dsd_node.func));
        memo[node_id] = sig;
        return sig;
      };

      auto root_sig = build_signal(root_id);
      old2new[topo.node_to_index(n)] = root_sig;
    } catch (std::exception const& e) {
      std::cerr << "[warning] --stp node decomposition failed: "
                << e.what()
                << ", fallback to original node\n";
      old2new[topo.node_to_index(n)] = decomposed.create_node(fanins, node_tt);
    }
  });

  topo.foreach_po([&](auto const& f) {
    auto idx = topo.node_to_index(topo.get_node(f));
    auto it = old2new.find(idx);
    if (it == old2new.end()) {
      throw std::runtime_error("stp_decompose_klut_network: missing mapped PO driver");
    }
    auto sig = it->second;
    decomposed.create_po(topo.is_complemented(f) ? !sig : sig);
  });

  return decomposed;
}
  struct lut_custom_cost {
    std::pair<uint32_t, uint32_t> operator()(uint32_t num_leaves) const {
      if (num_leaves < 2u) return {0u, 0u};
      return {num_leaves, 1u}; /* area, delay */
    }

    std::pair<uint32_t, uint32_t> operator()(
        kitty::dynamic_truth_table const& tt) const {
      if (tt.num_vars() < 2u) return {0u, 0u};
      return {tt.num_vars(), 1u}; /* area, delay */
    }
  };

  void execute() {
    if (is_set("mig")) {
      if (store<mig_network>().size() == 0u)
        std::cerr << "Error: Empty MIG network\n";
      else {
        auto mig = store<mig_network>().current();
        mapping_view mapped_mig{mig};
        phyLS::lut_map_params ps;
        if (is_set("area")) ps.area_oriented_mapping = true;
        if (is_set("relax_required")) ps.relax_required = relax_required;
        if (is_set("cut_size")) ps.cut_enumeration_ps.cut_size = cut_size;
        if (is_set("cut_limit")) ps.cut_enumeration_ps.cut_limit = cut_limit;
        if (is_set("recompute_cuts")) ps.recompute_cuts = false;
        if (is_set("edge")) ps.edge_optimization = false;
        if (is_set("dominated_cuts")) ps.remove_dominated_cuts = false;
        cout << "Mapped MIG into " << cut_size << "-LUT : ";
        phyLS::lut_map(mapped_mig, ps);
        mapped_mig.clear_mapping();
      }
    } else if (is_set("xag")) {
      if (store<xag_network>().size() == 0u)
        std::cerr << "Error: Empty XAG network\n";
      else {
        auto xag = store<xag_network>().current();
        mapping_view mapped_xag{xag};
        phyLS::lut_map_params ps;
        if (is_set("area")) ps.area_oriented_mapping = true;
        if (is_set("relax_required")) ps.relax_required = relax_required;
        if (is_set("cut_size")) ps.cut_enumeration_ps.cut_size = cut_size;
        if (is_set("cut_limit")) ps.cut_enumeration_ps.cut_limit = cut_limit;
        if (is_set("recompute_cuts")) ps.recompute_cuts = false;
        if (is_set("edge")) ps.edge_optimization = false;
        if (is_set("dominated_cuts")) ps.remove_dominated_cuts = false;
        cout << "Mapped XAG into " << cut_size << "-LUT : ";
        phyLS::lut_map(mapped_xag, ps);
        mapped_xag.clear_mapping();
      }
    } else if (is_set("xmg")) {
      if (store<xmg_network>().size() == 0u)
        std::cerr << "Error: Empty XMG network\n";
      else {
        auto xmg = store<xmg_network>().current();
        mapping_view mapped_xmg{xmg};
        phyLS::lut_map_params ps;
        if (is_set("area")) ps.area_oriented_mapping = true;
        if (is_set("relax_required")) ps.relax_required = relax_required;
        if (is_set("cut_size")) ps.cut_enumeration_ps.cut_size = cut_size;
        if (is_set("cut_limit")) ps.cut_enumeration_ps.cut_limit = cut_limit;
        if (is_set("recompute_cuts")) ps.recompute_cuts = false;
        if (is_set("edge")) ps.edge_optimization = false;
        if (is_set("dominated_cuts")) ps.remove_dominated_cuts = false;
        cout << "Mapped XMG into " << cut_size << "-LUT : ";
        phyLS::lut_map(mapped_xmg, ps);
        mapped_xmg.clear_mapping();
      }
    } else if (is_set("klut")) {
      if (store<klut_network>().size() == 0u)
        std::cerr << "Error: Empty k-LUT network\n";
      else {
        auto klut = store<klut_network>().current();
        mapping_view mapped_klut{klut};
        phyLS::lut_map_params ps;
        if (is_set("area")) ps.area_oriented_mapping = true;
        if (is_set("relax_required")) ps.relax_required = relax_required;
        if (is_set("cut_size")) ps.cut_enumeration_ps.cut_size = cut_size;
        if (is_set("cut_limit")) ps.cut_enumeration_ps.cut_limit = cut_limit;
        if (is_set("recompute_cuts")) ps.recompute_cuts = false;
        if (is_set("edge")) ps.edge_optimization = false;
        if (is_set("dominated_cuts")) ps.remove_dominated_cuts = false;
        cout << "Mapped kLUT into " << cut_size << "-LUT : ";
        phyLS::lut_map(mapped_klut, ps);
        //mapped_klut.clear_mapping();
        if (is_set("stp")) {
          auto stp_klut = stp_decompose_klut_network(klut);
          if (!stp_klut) {
            std::cerr << "[warning] --stp decomposition failed, keep first mapping\n";
            if (is_set("output")) {
              write_bench(mapped_klut, filename);
            } else {
              mapped_klut.clear_mapping();
            }
          } else {
            mapping_view mapped_stp{*stp_klut};
            cout << "Re-mapped STP decomposed kLUT into " << cut_size
                << "-LUT : ";
            phyLS::lut_map(mapped_stp, ps);
            if (is_set("output")) {
              write_bench(mapped_stp, filename);
            } else {
              mapped_stp.clear_mapping();
            }
          }
        } else {
          mapped_klut.clear_mapping();
        }
      }
    } else {
      if (store<aig_network>().size() == 0u)
        std::cerr << "Error: Empty AIG network\n";
      else {
        auto aig = store<aig_network>().current();
        mapping_view<aig_network, true> mapped_aig{aig};
        phyLS::lut_map_params ps;
        if (is_set("area")) ps.area_oriented_mapping = true;
        if (is_set("relax_required")) ps.relax_required = relax_required;
        if (is_set("cut_size")) ps.cut_enumeration_ps.cut_size = cut_size;
        if (is_set("cut_limit")) ps.cut_enumeration_ps.cut_limit = cut_limit;
        if (is_set("recompute_cuts")) ps.recompute_cuts = false;
        if (is_set("edge")) ps.edge_optimization = false;
        if (is_set("dominated_cuts")) ps.remove_dominated_cuts = false;
        cout << "Mapped AIG into " << cut_size << "-LUT : ";
        if (is_set("cost_function"))
          phyLS::lut_map<decltype(mapped_aig), true, lut_custom_cost>(
              mapped_aig, ps);
        else
          phyLS::lut_map(mapped_aig, ps);
        if (is_set("output")) {
          write_bench(mapped_aig, filename);
        } else {
          mapped_aig.clear_mapping();
        }
      }
    }
  }

 private:
  uint32_t cut_size{6u};
  uint32_t cut_limit{8u};
  uint32_t relax_required{0u};
  std::string filename = "lut.bench";
};

ALICE_ADD_COMMAND(lutmap, "Mapping")

}  // namespace alice

#endif
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
#include "../core/flow_detail.hpp"
#include "../../lib/my_s_t_p/src/include/algorithms/node_global.hpp"

#include <optional>
#include <functional>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <map>
#include <mockturtle/algorithms/collapse_mapped.hpp>
struct DSDNode {
  int id;
  std::string func;
  std::vector<int> child;
  int var_id = -1;
};

#include "../../lib/my_s_t_p/src/include/algorithms/strong_dsd.hpp"
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
        add_option("--bench", bench_filename,
               "export mapped 6-LUT bench after mapping");
    add_flag("--verbose, -v", "print the information");
   add_flag("--stp", "decompose original kLUT with lut_resyn -d style STP decomposition and remap");
    add_flag("--dec", "decompose original kLUT with lut_resyn -l style decomposition and remap");
  }

 protected:
   struct strong_dsd_nodes {
    int root_id{0};
    std::vector<DSDNode> nodes;
  };

    static std::optional<strong_dsd_nodes> capture_strong_dsd(
      kitty::dynamic_truth_table const& tt) {
    std::ostringstream oss;
    kitty::print_binary(tt, oss);

    RESET_NODE_GLOBAL();

  const auto prev_output = STRONG_DSD_DEBUG_PRINT;
    const auto prev_else_dec = ENABLE_ELSE_DEC;

    STRONG_DSD_DEBUG_PRINT = false;
    ENABLE_ELSE_DEC = true;

    const auto num_vars = tt.num_vars();
    ORIGINAL_VAR_COUNT = static_cast<int>(num_vars);

   std::vector<int> order(num_vars);
    std::iota(order.begin(), order.end(), 1);

     for (int v = 1; v <= ORIGINAL_VAR_COUNT; ++v) {
      new_in_node(v);
    }

 const auto root_id = build_strong_dsd_nodes(oss.str(), order, 0);
    strong_dsd_nodes result{root_id, NODE_LIST};

    STRONG_DSD_DEBUG_PRINT = prev_output;
    ENABLE_ELSE_DEC = prev_else_dec;

    if (root_id <= 0) {
      return std::nullopt;
    }

   return result;
  }


      struct stp_dsd_lut_resynthesis_local {
    template <typename LeavesIterator, typename Fn>
    void operator()(klut_network& ntk,
                    kitty::dynamic_truth_table const& function,
                    LeavesIterator begin, LeavesIterator end, Fn&& fn) const {
      std::vector<klut_network::signal> children(begin, end);
      if (children.size() <= 2u) {
        fn(ntk.create_node(children, function));
        return;
      }

         auto decomposition = capture_strong_dsd(function);
      if (!decomposition) {
        return;
      }
      std::unordered_map<int, klut_network::signal> var_to_signal;
      const auto n = children.size();
      for (auto i = 0u; i < n; ++i) {
        var_to_signal.emplace(static_cast<int>(n - i), children[i]);
      }

      std::unordered_map<int, DSDNode> node_lookup;
      for (auto const& node : decomposition->nodes) {
        node_lookup.try_emplace(node.id, node);
      }

 std::unordered_map<int, klut_network::signal> cache;


std::function<std::optional<klut_network::signal>(int)> build =
          [&](int id) -> std::optional<klut_network::signal> {
        if (auto it = cache.find(id); it != cache.end()) {
          return it->second;
        }

 auto node_it = node_lookup.find(id);
        if (node_it == node_lookup.end()) {
          return std::nullopt;
        }

            const auto& node = node_it->second;
        if (node.func == "in") {
          auto it = var_to_signal.find(node.var_id);
          if (it == var_to_signal.end()) return std::nullopt;
          cache.emplace(id, it->second);
          return it->second;
        }

       if (node.func == "0" || node.func == "1") {
          auto sig = ntk.get_constant(node.func == "1");
          cache.emplace(id, sig);
          return sig;
        }

   std::vector<klut_network::signal> args;
        args.reserve(node.child.size());
        for (auto cid : node.child) {
          auto c = build(cid);
          if (!c) return std::nullopt;
          args.push_back(*c);
        }
        if (args.empty()) return std::nullopt;

        uint32_t num_vars = 0u;
        while ((1u << num_vars) < node.func.size()) {
          ++num_vars;
        }
        auto tt_local = kitty::dynamic_truth_table(num_vars);
        kitty::create_from_binary_string(tt_local, node.func);
        auto sig = ntk.create_node(args, tt_local);
        cache.emplace(id, sig);
        return sig;
      };

     if (auto root = build(decomposition->root_id)) {
        fn(*root);
      }
    }
  };

  struct klut_dec_resynthesis {
    template <typename LeavesIterator, typename Fn>
    void operator()(klut_network& ntk,
                    kitty::dynamic_truth_table const& function,
                    LeavesIterator begin, LeavesIterator end, Fn&& fn) const {
      mockturtle::decomposition_flow_params ps;
      const std::vector<klut_network::signal> leaves(begin, end);
      const auto f = mockturtle::dsd_detail(ntk, function, leaves, ps);
      fn(f);
    }
  };

  static std::optional<klut_network> dec_decompose_klut_network(
      klut_network const& ntk) {
    try {
      klut_dec_resynthesis resyn;
      return node_resynthesis<klut_network>(ntk, resyn);
    } catch (std::exception const& e) {
     std::cerr << "[warning] --dec node decomposition failed: " << e.what()
                << ", fallback to original node\n";
      return std::nullopt;
    }
      }

 static std::optional<klut_network> stp_decompose_klut_network(
      klut_network const& ntk) {
    try {
      stp_dsd_lut_resynthesis_local resyn;
      return node_resynthesis<klut_network>(ntk, resyn);
    } catch (std::exception const& e) {
      std::cerr << "[warning] --stp node decomposition failed: " << e.what()
                << ", fallback to original node\n";
      return std::nullopt;
    }
      }

    template <typename MappedNetwork>
    void write_mapped_klut_bench_if_needed(const MappedNetwork& mapped) const {
      if (bench_filename.empty()) return;
      const auto klut = *collapse_mapped_network<klut_network>(mapped);
      write_bench(klut, bench_filename);
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

        phyLS::lut_map_params ps;
        if (is_set("area")) ps.area_oriented_mapping = true;
        if (is_set("relax_required")) ps.relax_required = relax_required;
        if (is_set("cut_size")) ps.cut_enumeration_ps.cut_size = cut_size;
        if (is_set("cut_limit")) ps.cut_enumeration_ps.cut_limit = cut_limit;
        if (is_set("recompute_cuts")) ps.recompute_cuts = false;
        if (is_set("edge")) ps.edge_optimization = false;
        if (is_set("dominated_cuts")) ps.remove_dominated_cuts = false;
        if (!is_set("stp") && !is_set("dec")) {
          mapping_view mapped_klut{klut};
          cout << "Mapped kLUT into " << cut_size << "-LUT : ";
          phyLS::lut_map(mapped_klut, ps);
           write_mapped_klut_bench_if_needed(mapped_klut);
          mapped_klut.clear_mapping();
          return;
        }

        if (is_set("stp") && is_set("dec")) {
          std::cerr << "[warning] both --stp and --dec are set, --stp is used\n";
        }
        if (is_set("stp")) {
          auto stp_klut = stp_decompose_klut_network(klut);
          if (!stp_klut) {
          std::cerr << "[warning] --stp decomposition failed, map original kLUT\n";
            mapping_view mapped_klut{klut};
            cout << "Mapped original kLUT into " << cut_size << "-LUT : ";
            phyLS::lut_map(mapped_klut, ps);
            if (is_set("output")) write_bench(mapped_klut, filename);
           // else mapped_klut.clear_mapping();
                       write_mapped_klut_bench_if_needed(mapped_klut);
            mapped_klut.clear_mapping();
          } else {
            mapping_view mapped_stp{*stp_klut};
            cout << "Re-mapped STP decomposed kLUT into " << cut_size
                << "-LUT : ";
            phyLS::lut_map(mapped_stp, ps);
            if (is_set("output")) {
              write_bench(mapped_stp, filename);
         write_mapped_klut_bench_if_needed(mapped_stp);
            mapped_stp.clear_mapping();
            }
          }
       } else if (is_set("dec")) {
          auto dec_klut = dec_decompose_klut_network(klut);
          if (!dec_klut) {
            std::cerr << "[warning] --dec decomposition failed, map original kLUT\n";
            mapping_view mapped_klut{klut};
            cout << "Mapped original kLUT into " << cut_size << "-LUT : ";
            phyLS::lut_map(mapped_klut, ps);
            if (is_set("output")) write_bench(mapped_klut, filename);
            write_mapped_klut_bench_if_needed(mapped_klut);
            mapped_klut.clear_mapping();
          } else {
            mapping_view mapped_dec{*dec_klut};
            cout << "Re-mapped DEC decomposed kLUT into " << cut_size
                 << "-LUT : ";
            phyLS::lut_map(mapped_dec, ps);
            if (is_set("output")) {
              write_bench(mapped_dec, filename);
            write_mapped_klut_bench_if_needed(mapped_dec);
            mapped_dec.clear_mapping();
            }
          }
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
        write_mapped_klut_bench_if_needed(mapped_aig);
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
 std::string bench_filename{};
};

ALICE_ADD_COMMAND(lutmap, "Mapping")

}  // namespace alice

#endif
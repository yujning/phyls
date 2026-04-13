#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

// globals / node creation
#include "node_global.hpp"

// kitty + mockturtle exact
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>

// =====================================================
// Strong DSD fallback: Shannon / Exact (placeholder-aware)
// 语义：
//   - n <= 4 : exact 2-LUT（终止）
//   - n > 4  : Shannon ONE layer，然后回 Strong DSD
// =====================================================


// ---------- helper: resolve a var_id to a node_id (placeholder-aware) ----------
inline int strong_resolve_var_node_id(
    int var_id,
    const std::vector<int>* local_to_global,
    const std::unordered_map<int,int>* placeholder_nodes )
{
  // 1) local placeholder
  if ( placeholder_nodes )
  {
    auto it = placeholder_nodes->find( var_id );
    if ( it != placeholder_nodes->end() )
      return it->second;
  }

  // 2) global placeholder bindings
  auto git = PLACEHOLDER_BINDINGS.find( var_id );
  if ( git != PLACEHOLDER_BINDINGS.end() )
    return git->second;

  // 3) local_to_global mapping
  int global_var = var_id;
  if ( local_to_global &&
       var_id >= 0 &&
       var_id < (int)local_to_global->size() &&
       (*local_to_global)[var_id] != 0 )
  {
    global_var = (*local_to_global)[var_id];
  }

  // 4) input node
  return new_in_node( global_var );
}

// ---------- helper: build children in kitty var order (LSB->MSB) ----------
inline std::vector<int> strong_make_children_kitty_order(
    const std::vector<int>& order_msb2lsb,
    const std::vector<int>* local_to_global,
    const std::unordered_map<int,int>* placeholder_nodes )
{
  std::vector<int> ch;
  ch.reserve( order_msb2lsb.size() );

  // kitty var0 = LSB, so bind from order.back() to order.front()
  for ( auto it = order_msb2lsb.rbegin(); it != order_msb2lsb.rend(); ++it )
  {
    ch.push_back( strong_resolve_var_node_id( *it, local_to_global, placeholder_nodes ) );
  }
  return ch;
}

// =====================================================
// exact refine (2-LUT) but placeholder-aware
// =====================================================
inline int strong_exact_refine_2lut(
    const std::string& mf,
    const std::vector<int>& order,
    int depth,
    const std::vector<int>* local_to_global,
    const std::unordered_map<int,int>* placeholder_nodes )
{
  const int n = (int)order.size();

  std::string indent((size_t)depth * 2, ' ');
  std::cout << indent << "⚠️ Strong EXACT refine (n=" << n << ")\n";
  std::cout << indent << "f=" << mf << "\n";

  for ( char c : mf )
    if ( c != '0' && c != '1' )
      throw std::runtime_error("strong_exact_refine_2lut: mf contains non-binary char");

  kitty::dynamic_truth_table tt( n );
  kitty::create_from_binary_string( tt, mf );
  std::cout << indent << "[DEBUG] kitty hex = " << kitty::to_hex( tt ) << "\n";

  mockturtle::klut_network klut;
  std::vector<mockturtle::klut_network::signal> pis;
  pis.reserve( n );

  for ( int i = 0; i < n; ++i )
    pis.push_back( klut.create_pi() );

  mockturtle::exact_resynthesis<mockturtle::klut_network> resyn( 2 );
  resyn( klut, tt, pis.begin(), pis.end(),
        [&]( auto const& s ) { klut.create_po( s ); } );

  std::cout << indent << "Exact 2-LUT count = " << klut.num_gates() << "\n";

  // Map KLUT nodes to your node IDs
  std::unordered_map<mockturtle::klut_network::node, int> node_map;

  // IMPORTANT: bind PI index in kitty order (LSB->MSB)
  auto orig_children = strong_make_children_kitty_order(
      order, local_to_global, placeholder_nodes );

  klut.foreach_pi( [&]( auto const& n_pi, auto index ) {
    if ( index < orig_children.size() )
      node_map[n_pi] = orig_children[index];
  } );

  klut.foreach_gate( [&]( auto const& n_gate ) {
    std::vector<int> childs;
    childs.reserve( klut.fanin_size( n_gate ) );

    klut.foreach_fanin( n_gate, [&]( auto const& f_handle ) {
      const auto f_node = klut.get_node( f_handle );

      if ( klut.is_constant( f_node ) )
      {
        const bool is_one = klut.is_complemented( f_handle ) ^ klut.constant_value( f_node );
        childs.push_back( new_node( is_one ? "1" : "0", {} ) );
      }
      else
      {
        auto cid = node_map.at( f_node );
        if ( klut.is_complemented( f_handle ) )
          cid = new_node( "01", { cid } );
        childs.push_back( cid );
      }
    } );

    // keep consistent with dsd_else_decompose:
    // reverse fanins before creating LUT node
    std::reverse( childs.begin(), childs.end() );

    auto func = klut.node_function( n_gate );
    std::string func_bin = kitty::to_binary( func );

    node_map[n_gate] = new_node( func_bin, childs );
  } );

  auto po_sig  = klut.po_at( 0 );
  auto po_node = klut.get_node( po_sig );

  if ( klut.is_constant( po_node ) )
  {
    bool val = klut.constant_value( po_node );
    if ( klut.is_complemented( po_sig ) ) val = !val;
    return new_node( val ? "1" : "0", {} );
  }

  auto root_id = node_map.at( po_node );
  if ( klut.is_complemented( po_sig ) )
    root_id = new_node( "01", { root_id } );

  return root_id;
}

// =====================================================
// Strong DSD fallback入口
// =====================================================
inline int strong_else_decompose(
    const std::string& mf,
    const std::vector<int>& order,
    int depth,

    // pivot node id (already placeholder-aware by caller)
    int pivot_node,

    const std::vector<int>* local_to_global,
    const std::unordered_map<int, int>* placeholder_nodes,

    // callback to Strong recursion
    int (*strong_rec)(
        const std::string&,
        const std::vector<int>&,
        int,
        const std::vector<int>*,
        const std::unordered_map<int, int>*))
{
  const int n = (int)order.size();
  std::string indent((size_t)depth * 2, ' ');

  // ---------- n <= 4 : exact 2-LUT ----------
  if ( n <= 4 )
  {
    return strong_exact_refine_2lut(
        mf, order, depth, local_to_global, placeholder_nodes );
  }

  // ---------- n > 4 : Shannon ONE layer ----------
  std::cout << indent << "⚠️ Strong fallback: Shannon ONE layer (n=" << n << ")\n";

  const size_t half = mf.size() / 2;
  std::string f_pos = mf.substr( 0, half );
  std::string f_neg = mf.substr( half );

  // pivot is MSB => remove order[0]
  std::vector<int> child_order( order.begin() + 1, order.end() );

  const int pos_node =
      strong_rec( f_pos, child_order, depth + 1, local_to_global, placeholder_nodes );

  const int neg_node =
      strong_rec( f_neg, child_order, depth + 1, local_to_global, placeholder_nodes );

  const int pos_term = new_node( "1000", { pivot_node, pos_node } );
  const int neg_term = new_node( "0010", { pivot_node, neg_node } );

  return new_node( "1110", { pos_term, neg_term } );
}

/* also: Advanced Logic Synthesis and Optimization tool
 * Copyright (C) 2019- Ningbo University, Ningbo, China */

/**
 * @file stp_dsd_resynthesis.hpp
 *
 * @brief Resynthesize KLUTs using STP strong DSD into 2-LUT structures.
 */

#ifndef STP_DSD_RESYNTHESIS_HPP
#define STP_DSD_RESYNTHESIS_HPP

#include <api/truth_table.hpp>

#include <algorithms/strong_dsd.hpp>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/node_resynthesis/sop_factoring.hpp>
#include <algorithm>

#include <numeric>
#include <sstream>

#include <functional>
#include <iostream>
#include <optional>
#include <unordered_map>
#include <vector>

namespace also
{

struct strong_dsd_nodes
{
  int root_id{0};
  std::vector<DSDNode> nodes;
};

inline std::optional<strong_dsd_nodes> capture_strong_dsd( kitty::dynamic_truth_table const& tt )
{
  std::ostringstream oss;
  kitty::print_binary( tt, oss );

  RESET_NODE_GLOBAL();

  const auto prev_output = STRONG_DSD_DEBUG_PRINT;
  const auto prev_else_dec = ENABLE_ELSE_DEC;

  STRONG_DSD_DEBUG_PRINT = false;
  ENABLE_ELSE_DEC = true;

  const auto num_vars = tt.num_vars();
  ORIGINAL_VAR_COUNT = static_cast<int>( num_vars );

  std::vector<int> order( num_vars );
  std::iota( order.begin(), order.end(), 1 );

  for ( int v = 1; v <= ORIGINAL_VAR_COUNT; ++v )
  {
    new_in_node( v );
  }

  const auto root_id = build_strong_dsd_nodes( oss.str(), order, 0 );

  strong_dsd_nodes result{ root_id, NODE_LIST };

  STRONG_DSD_DEBUG_PRINT = prev_output;
  ENABLE_ELSE_DEC = prev_else_dec;

  if ( root_id <= 0 )
  {
    return std::nullopt;
  }

  return result;
}

template<class Ntk>
class stp_dsd_lut_resynthesis
{
public:
  template<typename LeavesIterator, typename Fn>
  void operator()( Ntk& ntk, kitty::dynamic_truth_table const& function,
                  LeavesIterator begin, LeavesIterator end, Fn&& fn ) const
  {
    std::vector<typename Ntk::signal> children( begin, end );
    //std::cout << "[stp_dsd] klut truth table = ";
   // kitty::print_binary( function, std::cout );
    //std::cout << "\n";
    if ( children.size() <= 2u )
    {
      fn( ntk.create_node( children, function ) );
      return;
    }

    auto decomposition = capture_strong_dsd( function );
    if ( !decomposition )
    {
           std::cout << "[stp_dsd] decomposition failed.\n";
      return;
    }

    // ⭐⭐⭐ 关键修改：反转映射，匹配 bd 的变量编号约定
    std::unordered_map<int, typename Ntk::signal> var_to_signal;
    const auto n = children.size();
    for ( auto i = 0u; i < n; ++i )
    {
      // children[0] 是最低位 → 对应 bd 中的变量 n
      // children[n-1] 是最高位 → 对应 bd 中的变量 1
      var_to_signal.emplace( static_cast<int>( n - i ), children[i] );
    }

    std::unordered_map<int, DSDNode> node_lookup;
    for ( auto const& node : decomposition->nodes )
    {
      node_lookup.try_emplace( node.id, node );
    }

    std::unordered_map<int, typename Ntk::signal> cache;

    std::function<std::optional<typename Ntk::signal>( int )> build =
        [&]( int id ) -> std::optional<typename Ntk::signal> {
      if ( auto it = cache.find( id ); it != cache.end() )
      {
        return it->second;
      }

      auto node_it = node_lookup.find( id );
      if ( node_it == node_lookup.end() )
      {
        return std::nullopt;
      }

      const auto& node = node_it->second;

      typename Ntk::signal result{};

      if ( node.func == "in" )
      {
        if ( auto it = var_to_signal.find( node.var_id ); it != var_to_signal.end() )
        {
          result = it->second;
        }
        else
        {
          return std::nullopt;
        }
      }
      else if ( node.func == "0" )
      {
        result = ntk.get_constant( false );
      }
      else if ( node.func == "1" )
      {
        result = ntk.get_constant( true );
      }
      else
      {
        std::vector<int> child_ids = node.child;

        // ⭐⭐⭐ 关键：反转子节点顺序，匹配 mockturtle 约定
        std::reverse( child_ids.begin(), child_ids.end() );

        std::vector<typename Ntk::signal> fanins;
        fanins.reserve( child_ids.size() );
        for ( auto child_id : child_ids )
        {
          auto child_sig = build( child_id );
          if ( !child_sig )
          {
            return std::nullopt;
          }
          fanins.push_back( *child_sig );
        }

        if ( node.func == "01" && fanins.size() == 1u )
        {
          result = ntk.create_not( fanins.front() );
        }
        else
        {
          kitty::dynamic_truth_table tt( node.child.size() );
          for ( auto i = 0u; i < node.func.size(); ++i )
          {
            if ( node.func[i] == '1' )
            {
              const auto bit_index = node.func.size() - 1 - i;
              kitty::set_bit( tt, bit_index );
            }
          }
          result = ntk.create_node( fanins, tt );
        }
      }

      cache.emplace( id, result );
      return result;
    };

    if ( auto root = build( decomposition->root_id ) )
    {
      fn( *root );
    }
  }
};

template<class Ntk, class BaseResyn = mockturtle::sop_factoring<Ntk>>
class stp_dsd_resynthesis
{
public:
  explicit stp_dsd_resynthesis( BaseResyn base_resyn = {} )
      : base_resyn_( std::move( base_resyn ) )
  {}

  template<typename LeavesIterator, typename Fn>
  void operator()( Ntk& ntk, kitty::dynamic_truth_table const& function,
                   LeavesIterator begin, LeavesIterator end, Fn&& fn ) const
  {
    std::vector<typename Ntk::signal> children( begin, end );
    if ( children.size() <= 2u )
    {
      base_resyn_( ntk, function, children.begin(), children.end(), fn );
      return;
    }

    auto decomposition = capture_strong_dsd( function );
    if ( !decomposition )
    {
      std::cout << "[stp_dsd_resynthesis] fallback to base resynthesis (" << children.size()
                << " inputs)\n";      
      base_resyn_( ntk, function, children.begin(), children.end(), fn );
      return;
    }
    std::cout << "[stp_dsd_resynthesis] using strong DSD decomposition (" << children.size()
              << " inputs)\n";

    std::unordered_map<int, typename Ntk::signal> var_to_signal;
    const auto n = children.size();
    for ( auto i = 0u; i < n; ++i )
    {
      var_to_signal.emplace( static_cast<int>( n - i ), children[i] );
    }

    std::unordered_map<int, DSDNode> node_lookup;
    for ( auto const& node : decomposition->nodes )
    {
      node_lookup.try_emplace( node.id, node );
    }

    std::unordered_map<int, typename Ntk::signal> cache;

    std::function<std::optional<typename Ntk::signal>( int )> build =
        [&]( int id ) -> std::optional<typename Ntk::signal> {
      if ( auto it = cache.find( id ); it != cache.end() )
      {
        return it->second;
      }

      auto node_it = node_lookup.find( id );
      if ( node_it == node_lookup.end() )
      {
        return std::nullopt;
      }

      const auto& node = node_it->second;
      typename Ntk::signal result{};

      if ( node.func == "in" )
      {
        if ( auto it = var_to_signal.find( node.var_id ); it != var_to_signal.end() )
        {
          result = it->second;
        }
        else
        {
          return std::nullopt;
        }
      }
      else if ( node.func == "0" )
      {
        result = ntk.get_constant( false );
      }
      else if ( node.func == "1" )
      {
        result = ntk.get_constant( true );
      }
      else if ( node.func == "01" && node.child.size() == 1u )
      {
        auto child_sig = build( node.child.front() );
        if ( !child_sig )
        {
          return std::nullopt;
        }
        result = !*child_sig;
      }
      else
      {
        std::vector<int> child_ids = node.child;
        std::reverse( child_ids.begin(), child_ids.end() );

        std::vector<typename Ntk::signal> fanins;
        fanins.reserve( child_ids.size() );
        for ( auto child_id : child_ids )
        {
          auto child_sig = build( child_id );
          if ( !child_sig )
          {
            return std::nullopt;
          }
          fanins.push_back( *child_sig );
        }

        kitty::dynamic_truth_table tt( node.child.size() );
        for ( auto i = 0u; i < node.func.size(); ++i )
        {
          if ( node.func[i] == '1' )
          {
            const auto bit_index = node.func.size() - 1 - i;
            kitty::set_bit( tt, bit_index );
          }
        }

        bool success = false;
        base_resyn_( ntk, tt, fanins.begin(), fanins.end(),
                     [&]( auto const& f ) {
                       result = f;
                       success = true;
                       return false;
                     } );
        if ( !success )
        {
          return std::nullopt;
        }
      }

      cache.emplace( id, result );
      return result;
    };

    if ( auto root = build( decomposition->root_id ) )
    {
      fn( *root );
    }
    else
    {
       std::cout << "[stp_dsd] decomposition failed, fallback to base resynthesis.\n";
      base_resyn_( ntk, function, children.begin(), children.end(), fn );
    }
  }

private:
  BaseResyn base_resyn_;
};

} // namespace also



#endif
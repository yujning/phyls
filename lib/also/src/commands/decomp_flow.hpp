/* also: Advanced Logic Synthesis and Optimization tool
 * Copyright (C) 2019- Ningbo University, Ningbo, China */

/**
 * @file decomp_flow.hpp
 *
 * @brief decompose a truth table into a LUT-network
 *
 * @author Chengyu Ma
 * @author Zhufei Chu
 * @since  0.1
 */
#ifndef DECOMP_FLOW_HPP
#define DECOMP_FLOW_HPP

#include "../core/flow_detail.hpp"
#include <alice/alice.hpp>
#include <mockturtle/mockturtle.hpp>

namespace alice
{
class dec_command : public command
{
public:
  explicit dec_command( const environment::ptr& env ) : command( env, "Decompose a truth table into a LUT network." )
  {
    add_option( "--tt_from_hex, -a", tt_in, "truth table from hex" );
    add_option( "--tt_from_binary, -b", tt_in, "truth table from binary" );
    add_flag( "--noM3aj, -m", "don't ues M3aj,default = true" );
    add_flag( "--noXor, -x", "don't ues Xor,default = true" );
    add_flag( "--noMux, -u", "don't ues Mux,default = true" );
    add_flag( "--non-dsd, -n", "stop decomposition when func is non-dsd, default = false" );
    add_flag( "--shortsupport, -s", "The decomposition stops when the number of support <= 4, default = false " );
    add_flag( "--xmg", "convert klut to xmg" );
    add_flag( "--aig", "convert klut to aig" );
    add_flag( "--stp_dsd, -d", "use strong dsd with strong_else_dec in stp" );
  }

protected:
void execute()
{
  mockturtle::klut_network ntk;
  std::vector<mockturtle::klut_network::signal> children;
  kitty::dynamic_truth_table remainder;
  mockturtle::decomposition_flow_params ps;

  /* 1. 读取真值表，只创建 PI */
  if ( is_set( "--tt_from_hex" ) )
  {
    mockturtle::read_hex( tt_in, remainder, var_num, ntk, children );
  }
  else if ( is_set( "--tt_from_binary" ) )
  {
    mockturtle::read_binary( tt_in, remainder, var_num, ntk, children );
  }

  /* 2. 设置 dsd_detail 参数（仅在 non-STP 情况下使用） */
  if ( is_set( "noM3aj" ) ) ps.allow_m3aj = false;
  if ( is_set( "noXor" ) )  ps.allow_xor  = false;
  if ( is_set( "noMux" ) )  ps.allow_mux  = false;
  if ( is_set( "shortsupport" ) ) ps.allow_s = true;
  if ( is_set( "non-dsd" ) )      ps.allow_n = true;

  /* 3. 构建 klut 网络（核心分支） */
  if ( is_set( "stp_dsd" ) )
  {
    also::stp_dsd_lut_resynthesis<klut_network> resyn;
    klut_network::signal root{};

    resyn(
      ntk,
      remainder,
      children.begin(),
      children.end(),
      [&]( klut_network::signal const& s ) {
        root = s;
      } );

    ntk.create_po( root );
  }
  else
  {
    ntk.create_po(
      mockturtle::dsd_detail( ntk, remainder, children, ps ) );
  }

  

  /* 4. 图类型转换（必须在网络建好之后） */
  if ( is_set( "xmg" ) )
  {
    //xmg_network xmg = convert_klut_to_graph<xmg_network>( ntk );
    // store<xmg_network>().extend();
    // store<xmg_network>().current() = cleanup_dangling( xmg );

      xmg_npn_resynthesis resyn;
      exact_library_params eps;
      exact_library<xmg_network> lib( resyn, eps );
      map_params ps;
      map_stats st;
      auto mapped = mockturtle::map( ntk, lib, ps, &st );
      xmg_network xmg = cleanup_dangling( mapped );
      store<xmg_network>().extend();
      store<xmg_network>().current() = cleanup_dangling( xmg );
    return;
  }

  if ( is_set( "aig" ) )
  {
    aig_network aig = convert_klut_to_graph<aig_network>( ntk );
    store<aig_network>().extend();
    store<aig_network>().current() = cleanup_dangling( aig );
    return;
  }

  /* 5. 默认存 klut */
  store<klut_network>().extend();
  store<klut_network>().current() = cleanup_dangling( ntk );
}

private:
  std::string tt_in;
  int var_num;
};
ALICE_ADD_COMMAND( dec, "Decomposition" );
} // namespace alice

#endif

#ifndef DSD_HPP
#define DSD_HPP

#include <iostream>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <chrono>
#include <sstream>

#include <alice/alice.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>

// algorithms
#include "../include/algorithms/stp_dsd.hpp"
#include "../include/algorithms/reorder.hpp"     // all_reorders(raw)
#include "../include/algorithms/strong_dsd.hpp"  // build_strong_dsd_nodes(...)
#include "../include/algorithms/mix_dsd.hpp"     // run_dsd_recursive_mix(...)

namespace alice
{
  class dsd_command : public command
  {
  public:
    explicit dsd_command( const environment::ptr& env )
      : command( env, "Run STP decomposition / raw reorder" )
    {
      add_option( "-f, --factor", hex_input,
                  "hexadecimal number (must map to 2^n bits)" );

      add_option( "-x, --raw", raw_input,
                  "raw truth-table with don't care (x), length must be 2^n" );

      add_flag( "-s, --strong",
                "use strong DSD: find first L=2^k with exactly two block types (ACD)" );

      add_flag( "-m, --mix",
                "prefer DSD (-f) per layer, fallback to strong DSD when needed" );

      add_flag( "-e, --else",
                "enable prime fallback (Shannon / exact) for selected algorithm" );
    }

  protected:
    void execute() override
    {
      using clk = std::chrono::high_resolution_clock;

      const bool use_raw      = is_set( "raw" );
      const bool use_hex      = is_set( "factor" );
      const bool use_strong   = is_set( "strong" );
      const bool use_mix      = is_set( "mix" );
      const bool use_else_dec = is_set( "else" );

      if ( use_raw && use_hex )
      {
        std::cout << "❌ Options -f and -x cannot be used together.\n";
        return;
      }

      if ( use_strong && use_mix )
      {
        std::cout << "❌ Options -s and -m cannot be used together.\n";
        return;
      }

      // ------------------------------------------------------------
      // RAW MODE (-x)
      // ------------------------------------------------------------
      if ( use_raw )
      {
        const std::string raw = raw_input;

        if ( !is_power_of_two( raw.size() ) )
        {
          std::cout << "❌ Error: length (" << raw.size() << ") is not power of 2\n";
          return;
        }

        if ( raw.size() == 4 )
        {
          std::cout << "⚠ 输入函数已是 2-LUT，无需分解。\n";
          return;
        }

        std::cout << "➡ Raw truth-table mode (-x)\n";
        std::cout << "Input = " << raw << "\n";

        // 全局开关：让算法内部决定 -e 怎么用（STP/Strong/Mix 各自处理）
        ENABLE_ELSE_DEC = use_else_dec;

        // ---------- Strong DSD ----------
        if ( use_strong )
        {
          if ( raw.find_first_not_of( "01" ) != std::string::npos )
          {
            std::cout << "❌ Strong DSD requires raw input of only 0/1.\n";
            return;
          }

          const auto t1 = clk::now();

          RESET_NODE_GLOBAL();
          ENABLE_ELSE_DEC = use_else_dec;
          ORIGINAL_VAR_COUNT = static_cast<int>( std::log2( raw.size() ) );

          std::vector<int> order;
          order.reserve( ORIGINAL_VAR_COUNT );
          for ( int i = ORIGINAL_VAR_COUNT; i >= 1; --i )
            order.push_back( i );

          for ( int v = 1; v <= ORIGINAL_VAR_COUNT; ++v )
            new_in_node( v );

          // ✅ 永远走 strong（-e 不再劫持到 run_dsd_recursive）
          build_strong_dsd_nodes( raw, order, 0 );

          const auto t2 = clk::now();
          const auto us = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
          std::cout << "⏱ Strong DSD time = " << us << " us\n";
          return;
        }

        // ---------- Mix DSD ----------
        if ( use_mix )
        {
          if ( raw.find_first_not_of( "01" ) != std::string::npos )
          {
            std::cout << "❌ Mixed DSD requires raw input of only 0/1.\n";
            return;
          }

          const auto t1 = clk::now();

          // ✅ 永远走 mix（-e 由 mix 内部根据 ENABLE_ELSE_DEC 决定）
          run_dsd_recursive_mix( raw );

          const auto t2 = clk::now();
          const auto us = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
          std::cout << "⏱ Mixed DSD time = " << us << " us\n";
          return;
        }

        // ---------- Default: reorder ----------
        const auto t1 = clk::now();
        all_reorders( raw );
        const auto t2 = clk::now();

        const auto us = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        std::cout << "⏱ RAW decomposition time = " << us << " us\n";
        return;
      }

      // ------------------------------------------------------------
      // HEX MODE (-f)
      // ------------------------------------------------------------
      if ( !use_hex )
      {
        std::cout << "❌ Please use -f <hex> or -x <raw>\n";
        return;
      }

      std::string hex = hex_input;
      if ( hex.rfind( "0x", 0 ) == 0 || hex.rfind( "0X", 0 ) == 0 )
        hex = hex.substr( 2 );

      const unsigned bit_count = static_cast<unsigned>( hex.size() * 4 );
      unsigned num_vars        = 0;
      while ( ( 1u << num_vars ) < bit_count ) num_vars++;

      if ( ( 1u << num_vars ) != bit_count )
      {
        std::cout << "❌ Hex length is not 2^n bits\n";
        return;
      }

      if ( bit_count == 4 )
      {
        std::cout << "⚠ 输入函数已是 2-LUT，无需分解。\n";
        return;
      }

      kitty::dynamic_truth_table tt( num_vars );
      kitty::create_from_hex_string( tt, hex );

      std::ostringstream oss;
      kitty::print_binary( tt, oss );
      const std::string binary = oss.str();

      std::cout << "📘 Hex " << hex << " => binary " << binary
                << " (len = " << binary.size() << " vars = " << num_vars << ")\n";

      // 全局开关：让算法内部决定 -e 怎么用
      ENABLE_ELSE_DEC = use_else_dec;

      // ---------- Strong DSD ----------
      if ( use_strong )
      {
        const auto t1 = clk::now();

        RESET_NODE_GLOBAL();
        ENABLE_ELSE_DEC = use_else_dec;
        ORIGINAL_VAR_COUNT = static_cast<int>( std::log2( binary.size() ) );

        std::vector<int> order;
        order.reserve( ORIGINAL_VAR_COUNT );
        for ( int i = ORIGINAL_VAR_COUNT; i >= 1; --i )
          order.push_back( i );

        for ( int v = 1; v <= ORIGINAL_VAR_COUNT; ++v )
          new_in_node( v );

        // ✅ 永远走 strong（-e 不再劫持到 run_dsd_recursive）
        build_strong_dsd_nodes( binary, order, 0 );

        const auto t2 = clk::now();
        const auto us = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        std::cout << "⏱ Strong DSD time = " << us << " us\n";
        return;
      }

      // ---------- Mix DSD ----------
      if ( use_mix )
      {
        const auto t1 = clk::now();

        // ✅ 永远走 mix（-e 由 mix 内部根据 ENABLE_ELSE_DEC 决定）
        run_dsd_recursive_mix( binary );

        const auto t2 = clk::now();
        const auto us = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        std::cout << "⏱ Mixed DSD time = " << us << " us\n";
        return;
      }

      // ---------- Default: STP DSD ----------
      const auto t1 = clk::now();

      // ✅ STP DSD 自己会用 enable_else_dec（以及/或 ENABLE_ELSE_DEC）
      run_dsd_recursive( binary, use_else_dec );

      const auto t2 = clk::now();
      const auto us = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
      std::cout << "⏱ DSD execution time = " << us << " us\n";
    }

  private:
    std::string hex_input{};
    std::string raw_input{};
  };

  ALICE_ADD_COMMAND( dsd, "STP" )
} // namespace alice

#endif

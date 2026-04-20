/* also: Advanced Logic Synthesis and Optimization tool
 * Copyright ...
 */

#ifndef BALANCE_HPP
#define BALANCE_HPP

#include <mockturtle/mockturtle.hpp>
#include <alice/alice.hpp>

namespace alice
{

class balance_command : public command
{
public:
  explicit balance_command( const environment::ptr& env )
      : command( env, "Cut-based balancing of KLUT networks" )
  {
    add_flag( "--klut", "Balance KLUT network" );
  }

protected:
  void execute()
  {
    if( is_set( "klut" ) )
    {
      assert( store<klut_network>().size() > 0 );
      auto ntk = store<klut_network>().current();

      mockturtle::lut_map_params lps;  
      // 不设置 cut_size，使用默认 cut_enumeration 参数

      auto new_ntk = mockturtle::sop_balancing( ntk, lps );

      store<klut_network>().extend();
      store<klut_network>().current() = new_ntk;

      std::cout << "Balanced KLUT created.\n";
    }
    else
    {
      std::cout << "You must specify --klut.\n";
    }
  }
};

ALICE_ADD_COMMAND( balance, "Rewriting" )

} // namespace alice

#endif

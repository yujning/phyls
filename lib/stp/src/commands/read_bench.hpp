#ifndef RB_HPP
#define RB_HPP

#include <iostream>
#include <alice/alice.hpp>

#include "../include/algorithms/bench_lut.hpp"
#include "../include/algorithms/lut_func_cache.hpp"   

namespace alice
{

class read_bench_command : public command
{
public:
    explicit read_bench_command( const environment::ptr& env )
        : command( env, "Read BENCH LUT netlist" )
    {
        add_option( "file", filename, "BENCH file" )->required();
    }

protected:
    void execute() override
    {
        BenchNetlist net = read_bench_lut( filename );

        BENCH_NETLIST = net;
        BENCH_LOADED  = true;
        BENCH_SOURCE  = filename;

        // ✅ 清空 LUT resynthesis cache
        LutFuncCache::clear();

        std::cout << "📥 BENCH parsed\n";
        std::cout << "  Inputs  : " << net.inputs.size() << "\n";
        std::cout << "  Outputs : " << net.outputs.size() << "\n";
        std::cout << "  LUTs    : " << net.luts.size() << "\n\n";

        for ( const auto& kv : net.luts )
        {
            const auto& name = kv.first;
            const auto& lut  = kv.second;

            std::cout << "🔹 " << name << "\n";
            std::cout << "   hex    = " << lut.hex << "\n";
            std::cout << "   fanins = ";
            for ( const auto& f : lut.fanins )
                std::cout << f << " ";
            std::cout << "\n\n";
        }
    }

private:
    std::string filename;
};

ALICE_ADD_COMMAND( read_bench, "IO" );

} // namespace alice

#endif // RB_HPP

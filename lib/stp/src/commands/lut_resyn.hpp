#ifndef LUT_RESYN_HPP
#define LUT_RESYN_HPP

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <alice/alice.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>

#include "../include/algorithms/bench_lut.hpp"
#include "../include/algorithms/bi_decomposition.hpp"
#include "../include/algorithms/node_global.hpp"
#include "../include/algorithms/truth_table.hpp"

#include "../include/algorithms/stp_dsd.hpp"
#include "../include/algorithms/strong_dsd.hpp"
#include "../include/algorithms/mix_dsd.hpp"

#include "../include/algorithms/66lut_bidec.hpp"
#include "../include/algorithms/66lut_dsd.hpp"
#include "../include/algorithms/66lut_else_dec.hpp"

#include "../include/algorithms/lut_func_cache.hpp" // cache

namespace alice
{

enum class resyn_strategy
{
    bi_dec,
    dsd,
    strong_dsd,
    mix_dsd
};

/*============================================================*
 * BI-DECOMPOSITION
 *============================================================*/
static int run_bi_decomp_for_resyn(
    const std::string& binary01,
    bool enable_else_dec,
    bool enable_dsd_mix_fallback )
{
    bool prev_minimal_output = BD_MINIMAL_OUTPUT;
    bool prev_enable_else    = ENABLE_ELSE_DEC;
    bool prev_dsd_mix        = BD_ENABLE_DSD_MIX_FALLBACK;

    RESET_NODE_GLOBAL();
    ENABLE_ELSE_DEC            = enable_else_dec;
    BD_ENABLE_DSD_MIX_FALLBACK = enable_dsd_mix_fallback;
    BD_MINIMAL_OUTPUT          = true;

    if (!is_power_of_two(binary01.size()))
        throw std::runtime_error("input length must be power of two");

    int n = static_cast<int>(std::log2(binary01.size()));
    ORIGINAL_VAR_COUNT = n;

    TT root;
    root.f01 = binary01;
    root.order.resize(n);
    for (int i = 0; i < n; ++i)
        root.order[i] = n - i;

    for (int v = 1; v <= n; ++v)
        new_in_node(v);

    TT root_shrunk = shrink_to_support(root);
    int root_id = bi_decomp_recursive(root_shrunk, 0);

    BD_MINIMAL_OUTPUT          = prev_minimal_output;
    ENABLE_ELSE_DEC            = prev_enable_else;
    BD_ENABLE_DSD_MIX_FALLBACK = prev_dsd_mix;

    return root_id;
}

/*============================================================*
 * STRONG DSD (入口 reset -> set，避免 -e 粘住/被 reset 覆盖)
 *============================================================*/
static int run_strong_dsd_for_resyn(
    const std::string& binary01,
    bool enable_else_dec )
{
    if (!is_power_of_two(binary01.size()))
        throw std::runtime_error("input length must be power of two");

    int n = static_cast<int>(std::log2(binary01.size()));

    // ✅ 入口统一语义：reset → set
    RESET_NODE_GLOBAL();
    ENABLE_ELSE_DEC     = enable_else_dec;
    ORIGINAL_VAR_COUNT  = n;

    TT root;
    root.f01 = binary01;
    root.order.resize(n);
    for (int i = 0; i < n; ++i)
        root.order[i] = n - i;

    for (int v = 1; v <= n; ++v)
        new_in_node(v);

    TT root_shrunk = shrink_to_support(root);
    return build_strong_dsd_nodes(root_shrunk.f01, root_shrunk.order, 0);
}

/*============================================================*
 * DSD (入口 reset -> set，避免 -e 粘住/被 reset 覆盖)
 * 说明：run_dsd_recursive 内部可能依赖 ENABLE_ELSE_DEC，
 *       所以这里显式 reset + set，且不再硬编码 true。
 *============================================================*/
static int run_dsd_for_resyn(
    const std::string& binary01,
    bool enable_else_dec )
{
    if (!is_power_of_two(binary01.size()))
        throw std::runtime_error("input length must be power of two");

    int n = static_cast<int>(std::log2(binary01.size()));

    RESET_NODE_GLOBAL();
    ENABLE_ELSE_DEC    = enable_else_dec;
    ORIGINAL_VAR_COUNT = n;

    TT root;
    root.f01 = binary01;
    root.order.resize(n);
    for (int i = 0; i < n; ++i)
        root.order[i] = n - i;

    for (int v = 1; v <= n; ++v)
        new_in_node(v);

    TT root_shrunk = shrink_to_support(root);

    // ✅ 不要 run_dsd_recursive(binary01, true) 这种硬编码
    return run_dsd_recursive(root_shrunk.f01, enable_else_dec);
}

/*============================================================*
 * MIX DSD (入口 reset -> set)
 *============================================================*/
static int run_mix_dsd_for_resyn(
    const std::string& binary01,
    bool enable_else_dec )
{
    if (!is_power_of_two(binary01.size()))
        throw std::runtime_error("input length must be power of two");

    int n = static_cast<int>(std::log2(binary01.size()));

    RESET_NODE_GLOBAL();
    ENABLE_ELSE_DEC    = enable_else_dec;
    ORIGINAL_VAR_COUNT = n;

    TT root;
    root.f01 = binary01;
    root.order.resize(n);
    for (int i = 0; i < n; ++i)
        root.order[i] = n - i;

    for (int v = 1; v <= n; ++v)
        new_in_node(v);

    TT root_shrunk = shrink_to_support(root);

    if (enable_else_dec)
        return run_dsd_recursive(root_shrunk.f01, true);

    return run_dsd_recursive_mix(root_shrunk.f01);
}

/*============================================================*
 * LUT66
 *============================================================*/
static bool run_lut66_for_resyn(
    const std::string& binary01,
    int nvars,
    int& root_id,
    bool only_lut66 )
{
    TT root;
    root.f01 = binary01;
    root.order.resize(nvars);
    for (int i = 0; i < nvars; ++i)
        root.order[i] = nvars - i;

    TT root_shrunk = shrink_to_support(root);
    const unsigned shrunk_vars =
        static_cast<unsigned>(root_shrunk.order.size());

    const int max_var_id = root_shrunk.order.empty()
        ? 0
        : *std::max_element(root_shrunk.order.begin(),
                            root_shrunk.order.end());

    if (shrunk_vars <= 6)
    {
        RESET_NODE_GLOBAL();
        ORIGINAL_VAR_COUNT = max_var_id;

        std::vector<int> sorted_vars = root_shrunk.order;
        std::sort(sorted_vars.begin(), sorted_vars.end());
        sorted_vars.erase(
            std::unique(sorted_vars.begin(), sorted_vars.end()),
            sorted_vars.end());

        for (int var_id : sorted_vars)
            new_in_node(var_id);

        std::vector<int> children;
        for (int var_id : root_shrunk.order)
            children.push_back(new_in_node(var_id));

        new_node(root_shrunk.f01, children);

        root_id = NODE_LIST.back().id;
        return true;
    }

    bool success = run_66lut_dsd_and_build_dag(root_shrunk);
    if (!success)
        success = run_strong_bi_dec_and_build_dag(root_shrunk);
    if (!success && !only_lut66)
        success = run_66lut_else_dec_and_build_dag(root_shrunk);

    if (success)
        root_id = NODE_LIST.back().id;

    return success;
}

/*============================================================*
 * LUT RESYN COMMAND
 *============================================================*/
class lut_resyn_command : public command
{
public:
    explicit lut_resyn_command(const environment::ptr& env)
        : command(env, "Resynthesize LUT netlist to 2-LUT")
    {
        add_option("file", input_file, "BENCH file");
        add_option("-o,--output", output_file, "output BENCH file")->required();

        add_flag("-b,--bi_dec", use_bi_dec, "use bi-decomposition (default)");
        add_flag("-d,--dsd", use_dsd, "use DSD-based decomposition");
        add_flag("-s,--strong", use_strong_dsd, "use strong DSD (requires -d)");
        add_flag("-m,--mix", use_mix_dsd, "use mixed DSD (requires -d)");
        add_flag("-e,--else_dec", use_else_dec,
                 "exact synthesis + Shannon fallback");

        add_flag("--dm,--dsd_mix", use_dsd_mix_fallback,
                 "mixed DSD fallback when BD cannot proceed");

        add_flag("--lut66", use_lut66,
                 "per-LUT 66-LUT decomposition");

        add_flag("--only", use_lut66_only,
                 "66-LUT only, no fallback");
    }

protected:
    void execute() override
    {
        // ✅ 命令级：每次 lut_resyn 都重新算（不复用上次会话 cache）
        LutFuncCache::clear();

        if (output_file.empty())
        {
            std::cout << "❌ Output file missing\n";
            return;
        }

        BenchNetlist net;
        if (input_file.empty())
        {
            if (!BENCH_LOADED)
            {
                std::cout << "❌ No BENCH loaded\n";
                return;
            }
            net = BENCH_NETLIST;
        }
        else
        {
            net = read_bench_lut(input_file);
            BENCH_NETLIST = net;
            BENCH_LOADED  = true;
            BENCH_SOURCE  = input_file;
        }

        resyn_strategy strategy = resyn_strategy::bi_dec;
        if (use_dsd)
        {
            strategy = resyn_strategy::dsd;
            if (use_strong_dsd) strategy = resyn_strategy::strong_dsd;
            if (use_mix_dsd)    strategy = resyn_strategy::mix_dsd;
        }
        else if (use_bi_dec)
        {
            strategy = resyn_strategy::bi_dec;
        }

        std::ofstream fout(output_file);
        if (!fout)
        {
            std::cout << "❌ Cannot open " << output_file << "\n";
            return;
        }

        for (const auto& in : net.inputs)
            fout << "INPUT(" << in << ")\n";
        for (const auto& out : net.outputs)
            fout << "OUTPUT(" << out << ")\n";
        
        for (const auto& line : net.passthrough_lines)
            fout << line << "\n";
        fout << "\n";

        std::vector<std::string> lut_names;
        lut_names.reserve(net.luts.size());
        for (const auto& kv : net.luts)
            lut_names.push_back(kv.first);
        std::sort(lut_names.begin(), lut_names.end());

        int unique_id = 0;

                auto emit_original_lut = [&]( const std::string& lut_name, const LutNode& lut )
        {
            fout << lut_name << " = LUT " << lut.hex << " (";
            for (size_t i = 0; i < lut.fanins.size(); ++i)
            {
                if (i) fout << ", ";
                fout << lut.fanins[i];
            }
            fout << ")\n\n";
        };
        for (const auto& name : lut_names)
        {
            const auto& lut = net.luts.at(name);
                        try
            {

            if (lut.fanins.size() <= 2)
            {
                fout << name << " = LUT " << lut.hex << " (";
                for (size_t i = 0; i < lut.fanins.size(); ++i)
                {
                    if (i) fout << ", ";
                    fout << lut.fanins[i];
                }
                fout << ")\n";
                continue;
            }

            // cache key includes mode; but we still clear cache per command,
            // so no cross-command reuse.
            uint32_t mode = 0;
            mode |= static_cast<uint32_t>(strategy) & 0xFF;
            if (use_else_dec)         mode |= (1u << 8);
            if (use_dsd_mix_fallback) mode |= (1u << 9);
            if (use_lut66)            mode |= (1u << 10);
            if (use_lut66_only)       mode |= (1u << 11);

            std::string binary01 = hex_to_binary(lut.hex);
            LutFuncKey key{
                static_cast<uint32_t>(lut.fanins.size()),
                binary01,
                mode
            };

            if (!LutFuncCache::has(key))
            {
                int root_id = 0;

                if (use_lut66)
                {
                    bool success = run_lut66_for_resyn(
                        binary01,
                        static_cast<int>(lut.fanins.size()),
                        root_id,
                        use_lut66_only
                    );

                    if (!success && use_lut66_only)
                    {
                        emit_original_lut(name, lut);
                        continue;
                    }

                    if (!success && use_else_dec)
                    {
                        root_id = run_bi_decomp_for_resyn(binary01, true, true);
                        success = true;
                    }

                    if (!success)
                    {
           emit_original_lut(name, lut);
                        continue;
                    }
                }
                else
                {
                    // ✅ 彻底禁止硬编码 true 的调用路径
                    switch (strategy)
                    {
                    case resyn_strategy::bi_dec:
                        root_id = run_bi_decomp_for_resyn(
                            binary01, use_else_dec, use_dsd_mix_fallback);
                        break;

                    case resyn_strategy::dsd:
                        root_id = run_dsd_for_resyn(binary01, use_else_dec);
                        break;

                    case resyn_strategy::strong_dsd:
                        root_id = run_strong_dsd_for_resyn(binary01, use_else_dec);
                        break;

                    case resyn_strategy::mix_dsd:
                        root_id = run_mix_dsd_for_resyn(binary01, use_else_dec);
                        break;
                    }
                }

                CachedResyn entry;
                entry.nodes.reserve(NODE_LIST.size());
                for (const auto& n : NODE_LIST)
                {
                    CachedLutNode c;
                    c.id     = n.id;
                    c.func   = n.func;
                    c.var_id = n.var_id;
                    c.child  = n.child;
                    entry.nodes.push_back(std::move(c));
                }
                entry.root_id = root_id;

                LutFuncCache::insert(key, std::move(entry));
            }

            const auto& cached = LutFuncCache::get(key);

            std::map<int, std::string> name_of;

            // inputs
            for (const auto& node : cached.nodes)
            {
                if (node.func == "in")
                               {
                    if (node.var_id <= 0 ||
                        static_cast<size_t>(node.var_id) > lut.fanins.size())
                    {
                        throw std::runtime_error(
                            "internal variable index out of LUT fanin range");
                    }
                    name_of[node.id] = lut.fanins[node.var_id - 1];
                               }
            }

            // internal node names
            for (const auto& node : cached.nodes)
            {
                if (node.func == "in") continue;

                if (node.id == cached.root_id)
                    name_of[node.id] = name;
                else
                    name_of[node.id] = name + "_d" + std::to_string(++unique_id);
            }

            // emit
            for (const auto& node : cached.nodes)
            {
                if (node.func == "in") continue;

                fout << name_of[node.id] << " = LUT 0x"
                     << bin_to_hex(node.func) << " (";

                auto rev = node.child;
                std::reverse(rev.begin(), rev.end());

                for (size_t i = 0; i < rev.size(); ++i)
                {
                    if (i) fout << ", ";
                    fout << name_of[rev[i]];
                }
                fout << ")\n";
            }

            fout << "\n";
                        }
            catch (const std::exception& e)
            {
                std::cout << "⚠️  Skip resynthesis for LUT `" << name
                          << "`: " << e.what() << "\n";
                emit_original_lut(name, lut);
            }
        }

        std::cout << "✅ LUT resynthesis written to "
                  << output_file << "\n";

        use_bi_dec = false;
        use_dsd = false;
        use_strong_dsd = false;
        use_mix_dsd = false;
        use_else_dec = false;
        use_dsd_mix_fallback = false;
        use_lut66 = false;
        use_lut66_only = false;
    }

private:
    std::string input_file;
    std::string output_file;

    bool use_bi_dec = false;
    bool use_dsd = false;
    bool use_strong_dsd = false;
    bool use_mix_dsd = false;
    bool use_else_dec = false;
    bool use_dsd_mix_fallback = false;
    bool use_lut66 = false;
    bool use_lut66_only = false;
};

ALICE_ADD_COMMAND(lut_resyn, "STP")

} // namespace alice

#endif // LUT_RESYN_HPP

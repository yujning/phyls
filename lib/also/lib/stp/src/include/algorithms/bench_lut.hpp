#ifndef BENCH_LUT_HPP
#define BENCH_LUT_HPP

#include <string>
#include <vector>
#include <unordered_map>

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <cctype>


struct LutNode
{
    std::string name;                    // 节点名
    std::string hex;                     // 真值表（0x...）
    std::vector<std::string> fanins;     // 变量顺序（非常重要）
};

// ---------------- BENCH 网络 ----------------
struct BenchNetlist
{
    std::vector<std::string> inputs;
    std::vector<std::string> outputs;
    std::unordered_map<std::string, LutNode> luts;
    std::vector<std::string> passthrough_lines; //保留非lut的字段
};

// ---------------- BENCH 缓存 ----------------
inline BenchNetlist BENCH_NETLIST{};
inline bool BENCH_LOADED = false;
inline std::string BENCH_SOURCE{};


inline std::vector<std::string> tokenize( const std::string& line )
{
    std::vector<std::string> tokens;
    std::string cur;

    auto flush = [&]()
    {
        if (!cur.empty()) {
            tokens.push_back(cur);
            cur.clear();
        }
    };

    for (char c : line)
    {
        if (isspace(c) || c == ',' || c == '(' || c == ')' || c == '=')
            flush();
        else
            cur += c;
    }
    flush();
    return tokens;
}
inline BenchNetlist read_bench_lut( const std::string& filename )
{
    BenchNetlist net;
    std::ifstream fin(filename);

    if (!fin)
        throw std::runtime_error("cannot open file");

    std::string line;
    int lineno = 0;

    while (std::getline(fin, line))
    {
        lineno++;

        // 去注释
        auto pos = line.find('#');
        if (pos != std::string::npos)
            line = line.substr(0, pos);

        if (line.empty())
            continue;

        auto t = tokenize(line);
        if (t.empty())
            continue;

        // ---------- INPUT ----------
        if (t[0] == "INPUT")
        {
            net.inputs.push_back(t[1]);
        }
        // ---------- OUTPUT ----------
        else if (t[0] == "OUTPUT")
        {
            net.outputs.push_back(t[1]);
        }
        // ---------- LUT ----------
        else if (t.size() >= 4 && t[1] == "LUT")
        {
            LutNode node;
            node.name = t[0];
            node.hex  = t[2];

            if (node.hex.rfind("0x", 0) != 0)
            {
                throw std::runtime_error(
                    "line " + std::to_string(lineno) +
                    ": LUT hex must start with 0x"
                );
            }

            for (size_t i = 3; i < t.size(); ++i)
                node.fanins.push_back(t[i]);

            net.luts[node.name] = node;
        }
        else
        {
            // ⭐ 非 LUT 语句（如 vdd / gnd / assign），原样保留
            net.passthrough_lines.push_back(line);
        }

    }

    return net;
}

#endif
#ifndef WRITE_BENCH_HPP
#define WRITE_BENCH_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <alice/alice.hpp>
#include "../include/algorithms/truth_table.hpp"

// 来自 DSD 的节点 & 变量顺序
extern std::vector<DSDNode> NODE_LIST;
extern inline std::vector<int> FINAL_VAR_ORDER;
extern inline int ORIGINAL_VAR_COUNT;
extern inline int ROOT_NODE_ID;
namespace alice
{

inline bool is_binary_constant_func(const std::string& f)
{
    if (f.empty()) return false;
    return std::all_of(f.begin(), f.end(),
                       [&](char c){ return c == f[0]; });
}

// ========================================================
// 变量编号 → 字母
// varN → 'a', var(N-1) → 'b', …
// ========================================================
static std::string varname_from_id(int v)
{
    // int maxv = ORIGINAL_VAR_COUNT;
    // return std::string(1, char('a' + (maxv - v)));
        return std::string(1, char('a' + (v - 1)));
}

// ========================================================
// 安全二进制转 hex（任意长度，不会崩）
// ========================================================

// ========================================================
// write_bench
// ========================================================
class write_bench_command : public command
{
public:
    explicit write_bench_command(const environment::ptr& env)
        : command(env, "Write DSD result as LUT BENCH file")
    {
        add_option("file", filename,
            "output benchmark filename")->required();
    }

protected:
    void execute() override
    {
        if (filename.empty())
        {
            std::cout << "❌ No output file provided\n";
            return;
        }

        std::ofstream fout(filename);
        if (!fout)
        {
            std::cout << "❌ Cannot open " << filename << "\n";
            return;
        }

        if (NODE_LIST.empty())
        {
            std::cout << "❌ NODE_LIST is empty! (DSD not run?)\n";
            return;
        }

        if (ORIGINAL_VAR_COUNT == 0)
        {
            std::cout << "❌ ORIGINAL_VAR_COUNT is 0! (DSD not run?)\n";
            return;
        }

        // ====================================================
        // 1) 输入变量
        // ====================================================
        for (int v = 1; v <= ORIGINAL_VAR_COUNT; v++)
            fout << "INPUT(" << varname_from_id(v) << ")\n";

        fout << "OUTPUT(F0)\n\n";

        // ====================================================
        // 2) 节点命名
        // ====================================================
        std::map<int,std::string> name_of;
         
        
        for (auto &n : NODE_LIST)
        {
            if (n.func == "in")
                name_of[n.id] = varname_from_id(n.var_id);
        }

        

        for (auto &n : NODE_LIST)
        {
            if (n.func != "in")
                name_of[n.id] = "new_n" + std::to_string(n.id);
        }

        int root_id = ROOT_NODE_ID != 0 ? ROOT_NODE_ID : NODE_LIST.back().id;
        name_of[root_id] = "F0";

        // ====================================================
        // 3) 输出 LUT（🔥 全部 child 顺序反转！）
        // ====================================================
for (auto &n : NODE_LIST)
{
    if (n.func == "in")
        continue;

    // 🔥🔥🔥 第一优先级：binary 常量（不看 child 数）
    if (!n.func.empty() &&
        std::all_of(n.func.begin(), n.func.end(),
                    [&](char c){ return c == n.func[0]; }))
    {
        std::cout << "[DEBUG] const node " << n.id << " func=" << n.func << "\n";

        if (n.func[0] == '0')
            fout << name_of[n.id] << " = gnd\n";
        else
            fout << name_of[n.id] << " = vdd\n";
        continue;
    }

    // ===== 普通 LUT =====
    fout << name_of[n.id] << " = LUT 0x"
         << bin_to_hex(n.func) << " (";

    std::vector<int> rev = n.child;
    std::reverse(rev.begin(), rev.end());

    for (size_t i = 0; i < rev.size(); i++)
    {
        fout << name_of[rev[i]];
        if (i + 1 < rev.size())
            fout << ", ";
    }

    fout << ")\n";
}


        // ====================================================
        // Done
        // ====================================================
        std::cout << "✅ BENCH written to " << filename << "\n\n";

        std::cout << "📋 变量映射（最低位→'a'）：\n";
        //for (int v = ORIGINAL_VAR_COUNT; v >= 1; v--)
        for (int v = 1; v <= ORIGINAL_VAR_COUNT; v++)
            std::cout << "   变量" << v << " → '" 
                      << varname_from_id(v) << "'\n";
    }

private:
    std::string filename{};
};

ALICE_ADD_COMMAND(write_bench, "STP")

} // namespace alice

#endif

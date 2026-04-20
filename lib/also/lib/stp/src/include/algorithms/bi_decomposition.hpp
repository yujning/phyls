#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <set>
#include <kitty/kitty.hpp>

#include "reorder.hpp"
#include "node_global.hpp"
#include "bi_dec_else_dec.hpp"
#include "stp_dsd.hpp"
#include "mix_dsd.hpp"

using std::string;
using std::vector;
using std::set;

int new_node(const std::string&, const std::vector<int>&);
inline bool BD_MINIMAL_OUTPUT = false;
// When true, bi-decomposition search only considers cases with k2 == 0.
inline bool BD_ONLY_K2_EQ_0 = false;
// When true, abort BD recursion on failure to enable external fallback (e.g., DSD -m).
inline bool BD_ENABLE_DSD_MIX_FALLBACK = false;


// =====================================================
// BiDecompResult
// =====================================================
struct BiDecompResult {
    int k1, k2, k3;
    vector<int> Gamma;   // 原始变量编号
    vector<int> Theta;   // 原始变量编号
    vector<int> Lambda;  // 原始变量编号
    string F01;
    TT phi_tt;
    TT psi_tt;
};



// 打印结构矩阵（调试用，保持简单）
static void print_structure_matrix(
        int k1, int k2, int k3,
        const vector<vector<string>> &blk,
        const vector<int> &phi_bits,
        const vector<int> &psi_bits)
{
    int R = 1 << k1;
    int C = 1 << k2;

    std::cout << "Mf blocks:\n";
    for (int r = 0; r < R; r++)
    {
        std::cout << "  Row " << r << " : ";
        for (int c = 0; c < C; c++)
            std::cout << blk[r][c] << " ";
        std::cout << "\n";
    }

    std::cout << "\nφ table bits:\n";
    int idx = 0;
    for (int r = 0; r < R; r++)
    {
        std::cout << "  row " << r << " : ";
        for (int c = 0; c < C; c++)
            std::cout << phi_bits[idx++] << " ";
        std::cout << "\n";
    }

    std::cout << "\nψ bits:\n  ";
    for (int b : psi_bits) std::cout << b;
    std::cout << "\n\n";
}



static string apply_variable_reordering_swap(
    const string& f01,
    int n,
    const vector<int>& Gamma_indices,   // 1-based positions
    const vector<int>& Theta_indices,   // 1-based positions
    const vector<int>& Lambda_indices,  // 1-based positions
    int k1, int k2, int k3)
{
    // ---------- 1. 构造 new_order（原始变量编号，MSB->LSB） ----------
    vector<int> new_order;
    new_order.reserve(n);

    // 原始变量顺序是：in.order = 位置1..n
    // 这里 indices 是“位置”，不是变量编号
    for (int p : Gamma_indices)  new_order.push_back(p);
    for (int p : Theta_indices)  new_order.push_back(p);
    for (int p : Lambda_indices) new_order.push_back(p);

    // new_order 现在是【原始位置编号】
    // 我们需要把它映射为【变量编号顺序】
    // 即：变量编号 = 原 TT 的 order[pos-1]

    // ⚠️ 关键：f01 对应的原始变量顺序是 1..n（MSB->LSB）
    vector<int> old_order(n);
    for (int i = 0; i < n; ++i)
        old_order[i] = i + 1;

    vector<int> new_var_order(n);
    for (int i = 0; i < n; ++i)
        new_var_order[i] = new_order[i];

    // ---------- 2. 索引映射重排 ----------
    string out;
    out.resize(f01.size());

    // 建立 old_pos：变量编号 → bit 位置
    std::vector<int> old_pos(n + 1);
    for (int i = 0; i < n; ++i)
        old_pos[old_order[i]] = i;

    for (size_t new_idx = 0; new_idx < f01.size(); ++new_idx)
    {
        uint64_t old_idx = 0;

        for (int i = 0; i < n; ++i)
        {
            int bit = (new_idx >> (n - 1 - i)) & 1;
            int var = new_var_order[i];
            int pos = old_pos[var];

            old_idx |= (uint64_t(bit) << (n - 1 - pos));
        }

        out[new_idx] = f01[old_idx];
    }

    return out;
}


// =====================================================
// ⭐ k2=0 且 k3=1 特殊情况：直接从块序列提取全局 u,这里是将not节点融入顶层
// =====================================================
static vector<BiDecompResult>
handle_k2_eq_0_k3_eq_1_special(const TT& in, int k1, int k3)
{
    vector<BiDecompResult> results;

    const string &f01 = in.f01;
    int n = k1 + k3;  // k2=0
    
    if ((int)in.order.size() != n) return results;
    
    // ⭐ 只处理 k3=1 的情况
    if (k3 != 1) return results;

    int R = 1 << k1;  // 行数
    int B = 2;        // 块长度固定为 2

    std::cout << "\n🔷 特殊情况：k2=0, k3=1 (块长度=2)\n";

    // ========== 1. 提取所有块 ==========
    vector<string> blocks(R);
    for (int r = 0; r < R; ++r)
    {
        string block;
        for (int l = 0; l < B; ++l)
        {
            int idx = (r << k3) | l;
            block.push_back(f01[idx]);
        }
        blocks[r] = block;
    }

    std::cout << "📦 块序列：";
    for (const string &b : blocks) std::cout << b << " ";
    std::cout << "\n";

    // ========== 2. 统计不同的块类型 ==========
    set<string> unique_blocks(blocks.begin(), blocks.end());
    
    std::cout << "📊 不同的块类型：";
    for (const string &b : unique_blocks) std::cout << b << " ";
    std::cout << " (共 " << unique_blocks.size() << " 种)\n";

    // ========== 3. 检查是否可分解 ==========
    if (unique_blocks.size() > 2)
    {
        std::cout << "❌ 块类型超过 2 种，不可分解\n";
        return results;
    }

    if (unique_blocks.empty())
    {
        std::cout << "❌ 无有效块，不可分解\n";
        return results;
    }

    // ========== 4. 确定全局 u ==========
    vector<string> u_list(unique_blocks.begin(), unique_blocks.end());
    
    string global_u1 = u_list[0];
    string global_u2 = (u_list.size() == 2) ? u_list[1] : u_list[0];

    std::cout << "✅ 全局 u1 = " << global_u1 << ", u2 = " << global_u2 << "\n";

    // ========== 5. 构造 F ==========
    string F01 = global_u1 + global_u2;
    std::cout << "📌 F = " << F01 << "\n";

    // ========== 6. 强制 Mψ = [10...0] (恒等向量) ==========
    string Mpsi_fixed;
    Mpsi_fixed.push_back('1');
    for (int i = 1; i < B; ++i)
        Mpsi_fixed.push_back('0');

    std::cout << "📌 强制 Mψ = [" << Mpsi_fixed << "] (恒等向量)\n";

    // ========== 7. 构造 φ：根据块匹配 u1 或 u2 ==========
    // 定义 u 作用规则
    auto mul_u = [&](const string& u, const string& P) -> string {
        if (u == "10") return P;
        if (u == "01") {
            string out = P;
            for (char &c : out) c = (c=='0' ? '1' : '0');
            return out;
        }
        if (u == "11") return string(P.size(), '1');
        if (u == "00") return string(P.size(), '0');
        return P;
    };

    string g1 = mul_u(global_u1, Mpsi_fixed);
    string g2 = mul_u(global_u2, Mpsi_fixed);

    std::cout << "📌 u1·Mψ = " << g1 << "\n";
    std::cout << "📌 u2·Mψ = " << g2 << "\n";

    vector<int> phi_bits(R);
    bool valid = true;

    for (int r = 0; r < R; ++r)
    {
        if (blocks[r] == g1)
            phi_bits[r] = 1;
        else if (blocks[r] == g2)
            phi_bits[r] = 0;
        else
        {
            std::cout << "❌ 块 " << blocks[r] << " 无法匹配 u1·Mψ 或 u2·Mψ\n";
            valid = false;
            break;
        }
    }

    if (!valid)
        return results;

    std::cout << "✅ φ 构造成功：";
    for (int b : phi_bits) std::cout << b;
    std::cout << "\n";

    // ========== 8. 构造结果 ==========
    vector<int> Gamma, Lambda;
    for (int i = 0; i < k1; ++i)
        Gamma.push_back(in.order[i]);
    for (int i = 0; i < k3; ++i)
        Lambda.push_back(in.order[k1 + i]);

    BiDecompResult Rst;
    Rst.k1 = k1;
    Rst.k2 = 0;
    Rst.k3 = k3;
    Rst.Gamma = Gamma;
    Rst.Theta = {};  // 空
    Rst.Lambda = Lambda;
    Rst.F01 = F01;

    // φ(Γ) - 只依赖 Γ
    Rst.phi_tt.f01.resize(R);
    for (int i = 0; i < R; ++i)
        Rst.phi_tt.f01[i] = phi_bits[i] ? '1' : '0';
    Rst.phi_tt.order = Gamma;

    // ψ(Λ) - 只依赖 Λ
    Rst.psi_tt.f01 = Mpsi_fixed;
    Rst.psi_tt.order = Lambda;

    std::cout << "\n✅ k2=0 分解成功！\n";
    std::cout << "   Γ = { ";
    for (int v : Gamma) std::cout << v << " ";
    std::cout << "}\n";
    std::cout << "   Λ = { ";
    for (int v : Lambda) std::cout << v << " ";
    std::cout << "}\n";
    std::cout << "   φ = " << Rst.phi_tt.f01 << "\n";
    std::cout << "   ψ = " << Rst.psi_tt.f01 << "\n\n";

    results.push_back(Rst);
    return results;
}


// =====================================================
// 针对给定 k1,k2,k3，在当前 TT (in) 上尝试一次分解
//   注意：in.order 里存的是“原始变量编号”，顺序是 Γ,Θ,Λ
// =====================================================
static vector<BiDecompResult>
enumerate_one_case(const TT& in, int k1, int k2, int k3)
{
    vector<BiDecompResult> results;

    const string &f01 = in.f01;
    if (f01.empty()) return results;

    int n = (int)std::log2((double)f01.size());
    if ((int)in.order.size() != n) return results;
    if (k1 + k2 + k3 != n) return results;

    // 对称剪枝：Γ 和 Λ 等大的时候用首变量比较
    // if (k1 == k3)
    // {
    //     int gamma_first  = in.order[0];
    //     int lambda_first = in.order[k1 + k2];

    //     if (gamma_first > lambda_first)
    //         return results;
    // }

    int R = 1 << k1;
    int C = 1 << k2;
    int B = 1 << k3;

    if (k2 == 0 && k3 == 1)
        return handle_k2_eq_0_k3_eq_1_special(in, k1, k3);


    auto is_const_block2 = [&](const string& b){
        if (b.empty()) return false;
        for (char c : b)
            if (c != b[0]) return false;
        return true;
    };

    auto is_complement2 = [&](const string& a, const string& b){
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i)
            if (a[i] == b[i]) return false;
        return true;
    };


    const string& f01_used = f01;

    // 1) 构造块矩阵 Mf: blk[r][c]
   // 1) 构造块矩阵 Mf: blk[r][c]
    vector<vector<string>> blk(R, vector<string>(C));
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < C; ++c)
        {
            string s;
            for (int l = 0; l < B; ++l)
            {
                int idx = (r << (k2+k3)) | (c << k3) | l;
                s.push_back(f01_used[idx]);
            }
            blk[r][c] = s;
        }

    // 2) 每列做"模式分类"
    struct ColType
    {
        bool ok;
        vector<string> const_blocks;
        vector<string> nonconst_blocks;
    };

    vector<ColType> cols(C);

    for (int c = 0; c < C; ++c)
    {
        ColType ct{true, {}, {}};

        auto push_unique = [&](vector<string>& v, const string& s){
            for (auto &x : v) if (x == s) return;
            v.push_back(s);
        };

        for (int r = 0; r < R; ++r)
        {
            const string &b = blk[r][c];
            if (is_const_block2(b))
                push_unique(ct.const_blocks, b);
            else
                push_unique(ct.nonconst_blocks, b);
        }

        int nc = (int)ct.nonconst_blocks.size();
        int cc = (int)ct.const_blocks.size();

        // 合法情况判断
        if (nc == 0 && cc == 1)
        {
            // 纯常数列（单一常数）✓
        }
        else if (nc == 0 && cc == 2)
        {
            // 纯常数列（两种常数）✓
        }
        else if (nc == 1 && cc == 0)
        {
            // 纯非常数列 ✓
        }
        else if (nc == 1 && cc == 1)
        {
            // 一种常数 + 一种非常数 ✓
        }
        else if (nc == 2 && cc == 0 && is_complement2(ct.nonconst_blocks[0], ct.nonconst_blocks[1]))
        {
            // 互补非常数对 ✓
        }
        else
        {
            // 其他情况不合法
            ct.ok = false;
        }

        if (!ct.ok) return results;
        cols[c] = ct;
    }

    // 3) 收集必需的 u 类型（优先处理强约束列）
    set<string> required_u;

    // 第一遍：处理强约束列（互补对、混合列、纯非常数列）
    for (int c = 0; c < C; ++c)
    {
        const auto &ct = cols[c];
        int nc = (int)ct.nonconst_blocks.size();
        int cc = (int)ct.const_blocks.size();

        if (nc == 2 && cc == 0)
        {
            // 互补对：强制需要 10 和 01
            required_u.insert("10");
            required_u.insert("01");
        }
        else if (nc == 1 && cc == 1)
        {
            // 混合列：需要 非常数u + 对应常数u
            required_u.insert("10");  // 默认用10处理非常数块

            string const_val = ct.const_blocks[0];
            if (const_val == string(B, '1'))
                required_u.insert("11");
            else if (const_val == string(B, '0'))
                required_u.insert("00");
        }
        else if (nc == 1 && cc == 0)
        {
            // 纯非常数列
            required_u.insert("10");
        }
    }

    // 第二遍：处理纯常数列
    for (int c = 0; c < C; ++c)
    {
        const auto &ct = cols[c];
        int nc = (int)ct.nonconst_blocks.size();
        int cc = (int)ct.const_blocks.size();

        if (nc == 0 && cc == 1)
        {
            // 单一常数块的纯常数列
            if (required_u.size() < 2)
            {
                string const_val = ct.const_blocks[0];
                if (const_val == string(B, '1'))
                    required_u.insert("11");
                else if (const_val == string(B, '0'))
                    required_u.insert("00");
            }
        }
        else if (nc == 0 && cc == 2)
        {
            // 两种常数块的纯常数列
            if (required_u.size() < 2)
            {
                for (const string& b : ct.const_blocks)
                {
                    if (b == string(B, '0'))
                        required_u.insert("00");
                    else if (b == string(B, '1'))
                        required_u.insert("11");

                    if (required_u.size() >= 2) break;
                }
            }
        }
    }
    

    
    vector<string> u_types(required_u.begin(), required_u.end());

    // 4) 检查是否可分解
    if (u_types.size() > 2)
    {
        std::cout << "  ⚠️  需要 " << u_types.size() << " 种 u，不可分解（";
        for (const string &u : u_types) std::cout << u << " ";
        std::cout << "）\n";
        return results;
    }

    if (u_types.empty())
    {
        return results;
    }

    // 5) 确定全局 u1 和 u2
    string global_u1, global_u2;
    if (u_types.size() == 2)
    {
        global_u1 = u_types[0];
        global_u2 = u_types[1];
    }
    else if (u_types.size() == 1)
    {
        global_u1 = u_types[0];
        global_u2 = global_u1;  // 只有1种u，重复使用
    }

    // 6) 构造 F
    string F01;
    if (u_types.size() == 2)
        F01 = global_u1 + global_u2;
    else
        F01 = global_u1 + global_u1;

    // =================================================================
    // 7) ★★ 用严格 u·Mψ=W 规则重算 φ-table & ψ-table（替换你原来的 heuristic）
    // =================================================================

    // u 作用在列向量 P 上：支持 10,01,11,00 四种
    auto mul_u = [&](const string& u, const string& P)->string {
        if (u == "10") return P;                      // 恒等
        if (u == "01") {                              // 取反
            string out = P;
            for (char &c : out) c = (c=='0' ? '1' : '0');
            return out;
        }
        if (u == "11") return string(P.size(), '1');  // 常 1
        if (u == "00") return string(P.size(), '0');  // 常 0
        return P;
    };

    // 按列存 Mψ，每列一个长度为 B 的字符串
    vector<string> Mpsi_cols(C);
    // φ 是 R×C 的 0/1 表
    vector<vector<int>> phi_mat(R, vector<int>(C, 0));

    for (int c = 0; c < C; ++c)
    {
        // 这一列所有块
        vector<string> col_blocks(R);
        for (int r = 0; r < R; ++r)
            col_blocks[r] = blk[r][c];

        // 候选 P：这一列出现过的所有块去重
        vector<string> candidates;
        for (auto &b : col_blocks)
        {
            if (std::find(candidates.begin(), candidates.end(), b) == candidates.end())
                candidates.push_back(b);
        }

        bool solved = false;
        string chosen_P;

        // 尝试每一个候选 P 作为 Mψ 列向量
        for (auto &P : candidates)
        {
            bool ok = true;

            for (auto &W : col_blocks)
            {
                string g1 = mul_u(global_u1, P);
                string g2 = mul_u(global_u2, P);

                if (W != g1 && W != g2)
                {
                    ok = false;
                    break;
                }
            }

            if (ok)
            {
                solved = true;
                chosen_P = P;
                break;
            }
        }

        if (!solved)
        {
            std::cout << "  ⚠️  列 " << c 
                      << " 无法找到统一的 Mψ 使得所有块都来自 {u1,u2}·Mψ，判定该 (k1,k2,k3) 不可分解\n";
            return results;
        }

        // 选定这一列的 Mψ
        Mpsi_cols[c] = chosen_P;

        // 根据 chosen_P 填这一列的 φ：W == u1·P → φ=1，否则 φ=0（即用 u2）
        string g1 = mul_u(global_u1, chosen_P);
        string g2 = mul_u(global_u2, chosen_P);

        for (int r = 0; r < R; ++r)
        {
            const string &W = col_blocks[r];

            if (W == g1)
                phi_mat[r][c] = 1;
            else
                phi_mat[r][c] = 0;   // solvable 保证 W==g1 或 W==g2
        }
    }

    // 展平成你原来使用的一维 phi_bits / psi_bits
    vector<int> phi_bits(R*C);
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < C; ++c)
            phi_bits[r*C + c] = phi_mat[r][c];

    vector<int> psi_bits(C*B);
    for (int c = 0; c < C; ++c)
        for (int l = 0; l < B; ++l)
            psi_bits[c*B + l] = (Mpsi_cols[c][l] == '1' ? 1 : 0);

    // =================================================================
    // 8) 后面全是你原来的逻辑：切 Γ/Θ/Λ，构造 TT，打印，一行未动
    // =================================================================

    // 9) 从 in.order 切出 Γ,Θ,Λ（这里 in.order 已经是“原始编号的 [Γ,Θ,Λ] 顺序”）
    vector<int> Gamma, Theta, Lambda;
    for (int i = 0; i < k1; ++i) Gamma.push_back(in.order[i]);
    for (int i = 0; i < k2; ++i) Theta.push_back(in.order[k1 + i]);
    for (int i = 0; i < k3; ++i) Lambda.push_back(in.order[k1 + k2 + i]);

    BiDecompResult Rst;
    Rst.k1 = k1;
    Rst.k2 = k2;
    Rst.k3 = k3;
    Rst.Gamma = Gamma;
    Rst.Theta = Theta;
    Rst.Lambda = Lambda;
    Rst.F01 = F01;

    Rst.phi_tt.f01.resize(R*C);
    for (int i = 0; i < R*C; ++i)
        Rst.phi_tt.f01[i] = phi_bits[i] ? '1' : '0';

    // ★★ 这里是关键：phi_tt.order 里放“原始变量编号”，顺序为 Γ,Θ
    Rst.phi_tt.order.clear();
    for (int v : Gamma) Rst.phi_tt.order.push_back(v);
    for (int v : Theta) Rst.phi_tt.order.push_back(v);

    Rst.psi_tt.f01.resize(C*B);
    for (int i = 0; i < C*B; ++i)
        Rst.psi_tt.f01[i] = psi_bits[i] ? '1' : '0';

    // ★★ 同理：psi_tt.order = Θ,Λ（原始编号）
    Rst.psi_tt.order.clear();
    for (int v : Theta)  Rst.psi_tt.order.push_back(v);
    for (int v : Lambda) Rst.psi_tt.order.push_back(v);

    std::cout << "\n===== Matrix Form (Theorem 4.2) =====\n";
    std::cout << "k1=" << k1 << "  k2=" << k2 << "  k3=" << k3 << "\n";
    std::cout << "F = " << Rst.F01 << "\n";
    std::cout << "u_types: ";
    for (const string &u : u_types) std::cout << u << " ";
    std::cout << "\n";
    std::cout << "global_u1 = " << global_u1 << ", global_u2 = " << global_u2 << "\n\n";
    print_structure_matrix(k1, k2, k3, blk, phi_bits, psi_bits);

    results.push_back(Rst);
    return results;
}



static bool
find_first_bi_decomposition(const TT& in, BiDecompResult& out)
{
    const string &f01 = in.f01;
    if (f01.empty()) return false;

    int n = (int)std::log2((double)f01.size());
    if ((int)in.order.size() != n) return false;

    // 枚举 k2 和 k3 的大小
    const int k2_begin = BD_ONLY_K2_EQ_0 ? 0 : 0;
    const int k2_end   = BD_ONLY_K2_EQ_0 ? 0 : (n - 2);
    for (int k2 = k2_begin; k2 <= k2_end; ++k2)
    {
        int max_k3 = (n - k2) / 2;

        // 修改：k3 从 max_k3 递减到 1（或0，看您需求）
        for (int k3 = max_k3; k3 >= 1; --k3) 
        {
            int k1 = n - k2 - k3;
            if (k1 <= 0) continue;

            //std::cout << "\n========== 尝试 k1=" << k1 << ", k2=" << k2 << ", k3=" << k3 << " ==========\n";

            if (!BD_MINIMAL_OUTPUT)
            std::cout << "\n========== 尝试 k1=" << k1 << ", k2=" << k2 << ", k3=" << k3 << " ==========\n";

            // 先试试不重排的情况（变量已经是 [Γ,Θ,Λ] 顺序）
            auto sub = enumerate_one_case(in, k1, k2, k3);
            if (!sub.empty())
            {
                out = sub[0];
                //std::cout << "✓ 不需重排即可分解！\n";
                                if (!BD_MINIMAL_OUTPUT)
                    std::cout << "✓ 不需重排即可分解！\n";
                return true;
            }

            // 枚举 Θ 的所有 C(n, k2) 种组合（按位置，1-based）
            vector<bool> theta_mask(n, false);
            std::fill(theta_mask.begin(), theta_mask.begin() + k2, true);

            do {
                vector<int> Theta_pos;  // Θ的位置
                for (int i = 0; i < n; ++i)
                    if (theta_mask[i])
                        Theta_pos.push_back(i + 1);

                // 剩余位置用于分配 Γ 和 Λ
                vector<int> remaining_pos;
                for (int i = 0; i < n; ++i)
                    if (!theta_mask[i])
                        remaining_pos.push_back(i + 1);

                // 枚举 Λ 的所有 C(n-k2, k3) 种组合
                vector<bool> lambda_mask(remaining_pos.size(), false);
                std::fill(lambda_mask.begin(), lambda_mask.begin() + k3, true);

                do {
                    vector<int> Lambda_pos;  // Λ的位置
                    vector<int> Gamma_pos;   // Γ的位置

                    for (size_t i = 0; i < remaining_pos.size(); ++i)
                    {
                        if (lambda_mask[i])
                            Lambda_pos.push_back(remaining_pos[i]);
                        else
                            Gamma_pos.push_back(remaining_pos[i]);
                    }

                    // ⭐ 对称性剪枝：当 k1 == k3 时，要求 Γ 的首位置 < Λ 的首位置
                    if (k1 == k3 && Gamma_pos[0] > Lambda_pos[0])
                        continue;

                    // 打印当前尝试
                    // std::cout << "  尝试位置：Γ={";
                    // for (int p : Gamma_pos) std::cout << p << " ";
                    // std::cout << "}, Θ={";
                    // for (int p : Theta_pos) std::cout << p << " ";
                    // std::cout << "}, Λ={";
                    // for (int p : Lambda_pos) std::cout << p << " ";
                    // std::cout << "} → 变量 Γ={";
                    // for (int p : Gamma_pos) std::cout << in.order[p-1] << " ";
                    // std::cout << "}, Θ={";
                    // for (int p : Theta_pos) std::cout << in.order[p-1] << " ";
                    // std::cout << "}, Λ={";
                    // for (int p : Lambda_pos) std::cout << in.order[p-1] << " ";
                    // std::cout << "}\n";
                                        if (!BD_MINIMAL_OUTPUT)
                    {
                        // 打印当前尝试
                        std::cout << "  尝试位置：Γ={";
                        for (int p : Gamma_pos) std::cout << p << " ";
                        std::cout << "}, Θ={";
                        for (int p : Theta_pos) std::cout << p << " ";
                        std::cout << "}, Λ={";
                        for (int p : Lambda_pos) std::cout << p << " ";
                        std::cout << "} → 变量 Γ={";
                        for (int p : Gamma_pos) std::cout << in.order[p-1] << " ";
                        std::cout << "}, Θ={";
                        for (int p : Theta_pos) std::cout << in.order[p-1] << " ";
                        std::cout << "}, Λ={";
                        for (int p : Lambda_pos) std::cout << in.order[p-1] << " ";
                        std::cout << "}\n";
                    }

                    // ⭐ 重排真值表：按 [Γ, Θ, Λ] 的位置顺序
                    string reordered_f01 = apply_variable_reordering_swap(
                        f01, n,
                        Gamma_pos, Theta_pos, Lambda_pos,
                        k1, k2, k3
                    );


                        //std::cout << "📌 重排后的 f01（二进制） = " << reordered_f01 << "\n";
                        if (!BD_MINIMAL_OUTPUT)
                            std::cout << "📌 重排后的 f01（二进制） = " << reordered_f01 << "\n";

                    // 构造重排后的 TT，order 保存原始变量编号
                    TT reordered_tt;
                    reordered_tt.f01 = reordered_f01;
                    reordered_tt.order.clear();

                    // 按 [Γ, Θ, Λ] 顺序记录原始变量编号
                    for (int pos : Gamma_pos)
                        reordered_tt.order.push_back(in.order[pos - 1]);
                    for (int pos : Theta_pos)
                        reordered_tt.order.push_back(in.order[pos - 1]);
                    for (int pos : Lambda_pos)
                        reordered_tt.order.push_back(in.order[pos - 1]);

                    // 在重排后的真值表上尝试分解
                    sub = enumerate_one_case(reordered_tt, k1, k2, k3);

                    if (!sub.empty())
                    {
                        out = sub[0];
                        //std::cout << "    ✓ 找到分解！\n";
                                                if (!BD_MINIMAL_OUTPUT)
                            std::cout << "    ✓ 找到分解！\n";
                        return true;
                    }

                } while (std::prev_permutation(lambda_mask.begin(), lambda_mask.end()));

            } while (std::prev_permutation(theta_mask.begin(), theta_mask.end()));
        }
    }

    std::cout << "❌ 遍历所有 (k1,k2,k3) 和变量分组，未找到有效分解\n";
    return false;
}


// =====================================================
// 递归双分解（参照 DSD 的编号和递归方式）
// =====================================================
static int bi_decomp_recursive(const TT& f, int depth = 0)
{
    int len = f.f01.size();
    int nv  = f.order.size();

    // 基本情况：2输入或更少，直接建小树
    if (len <= 4)
        return build_small_tree(f);

            // 如果开启 DSD 混合模式，优先尝试 DSD -m；失败时再回到 BD
     // 如果开启 DSD 混合模式，优先尝试 DSD -m；失败时再回到 BD
    if (BD_ENABLE_DSD_MIX_FALLBACK)
    {
        int max_var = 0;
        for (int v : f.order)
            max_var = std::max(max_var, v);

        std::vector<int> local_to_global(max_var + 1, 0);
        for (int v : f.order)
            local_to_global[v] = v;

        auto mix_try = dsd_factor_mix_impl(f, depth, &local_to_global, nullptr, false);

        if (mix_try.decomposed && mix_try.fully_success)
        {
            std::cout << "✅ 深度 " << depth << "：DSD -m 完全成功，保持 DSD 分支\n";
            return mix_try.node_id;
        }

        std::cout << "ℹ️ 深度 " << depth << "：DSD -m 未完全成功，回退到 BD\n";

    
    }

    // 尝试找到第一个双分解
    BiDecompResult result;
    bool found = find_first_bi_decomposition(f, result);

    if (!found)
    {
        if (ENABLE_ELSE_DEC)
        {
            std::cout << "⚠️ 深度 " << depth
                      << "：无法双分解 → 触发 else_dec 回退 (n=" << nv << ")\n";

            auto orig_children = make_children_from_order(f);
            return else_decompose(f, orig_children, depth);
        }

        if (BD_ENABLE_DSD_MIX_FALLBACK)
        {
            std::cout << "⚠️ 深度 " << depth << "：无法双分解 → 触发 DSD -m 回退\n";
            int max_var = 0;
            for (int v : f.order)
                max_var = std::max(max_var, v);

            std::vector<int> local_to_global(max_var + 1, 0);
            for (int v : f.order)
                local_to_global[v] = v;


            auto mix = dsd_factor_mix_impl(f, depth, &local_to_global, nullptr, true);

            if (mix.fully_success && mix.node_id >= 0)
            {
                std::cout << "✅ 深度 " << depth << "：DSD -m fallback 完全成功\n";
                return mix.node_id;
            }

            std::cout << "⚠️ 深度 " << depth << "：DSD -m fallback 失败，继续 else_dec / build_small_tree\n";
            // 不 return，继续走后面的 else_dec / build_small_tree

        }


        std::cout << "⚠️ 深度 " << depth << "：无法双分解 → 直接建树\n";
        return build_small_tree(f);
    }



    // 打印信息
    // std::cout << "\n" << string(depth*2, ' ') << "📌 深度 " << depth << " 双分解成功：\n";
    // std::cout << string(depth*2, ' ') << "   k1=" << result.k1
    //           << "  k2=" << result.k2 << "  k3=" << result.k3 << "\n";

    // std::cout << string(depth*2, ' ') << "   Γ = { ";
    // for (int v : result.Gamma) std::cout << v << " ";
    // std::cout << "}\n";

    // std::cout << string(depth*2, ' ') << "   Θ = { ";
    // for (int v : result.Theta) std::cout << v << " ";
    // std::cout << "}\n";

    // std::cout << string(depth*2, ' ') << "   Λ = { ";
    // for (int v : result.Lambda) std::cout << v << " ";
    // std::cout << "}\n";

    // std::cout << string(depth*2, ' ') << "   F(u,v) = " << result.F01 << "\n";

        if (BD_MINIMAL_OUTPUT)
    {
        std::cout << "\n" << string(depth*2, ' ') << "深度 " << depth
                  << " 可分解真值表：" << f.f01 << "\n";
    }
    else
    {
        std::cout << "\n" << string(depth*2, ' ') << "📌 深度 " << depth << " 双分解成功：\n";
        std::cout << string(depth*2, ' ') << "   k1=" << result.k1
                  << "  k2=" << result.k2 << "  k3=" << result.k3 << "\n";
        std::cout << string(depth*2, ' ') << "   Γ = { ";
        for (int v : result.Gamma) std::cout << v << " ";
        std::cout << "}\n";
        std::cout << string(depth*2, ' ') << "   Θ = { ";
        for (int v : result.Theta) std::cout << v << " ";
        std::cout << "}\n";
        std::cout << string(depth*2, ' ') << "   Λ = { ";
        for (int v : result.Lambda) std::cout << v << " ";
        std::cout << "}\n";
        std::cout << string(depth*2, ' ') << "   F(u,v) = " << result.F01 << "\n";
    }

    // 记录变量到 FINAL_VAR_ORDER（全是原始编号）
    for (int v : result.Gamma)
        if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), v) == FINAL_VAR_ORDER.end())
            FINAL_VAR_ORDER.push_back(v);
    for (int v : result.Theta)
        if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), v) == FINAL_VAR_ORDER.end())
            FINAL_VAR_ORDER.push_back(v);
    for (int v : result.Lambda)
        if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), v) == FINAL_VAR_ORDER.end())
            FINAL_VAR_ORDER.push_back(v);

    // 准备 φ 和 ψ 的递归
    TT phi_tt = result.phi_tt;
    TT psi_tt = result.psi_tt;

    int n_phi = phi_tt.order.size();
    int n_psi = psi_tt.order.size();

    std::cout << string(depth*2, ' ') << "📌 递归分解 φ：原始变量 { ";
    for (int v : phi_tt.order) std::cout << v << " ";
    std::cout << "} → 局部编号 { ";
    for (int i = 1; i <= n_phi; i++) std::cout << i << " ";
    std::cout << "}\n";
    std::cout << string(depth*2, ' ') << "   映射关系：";
    for (int i = 0; i < n_phi; i++)
        std::cout << "位置" << (i+1) << "→变量" << phi_tt.order[i] << " ";
    std::cout << "\n";

    std::cout << string(depth*2, ' ') << "📌 递归分解 ψ：原始变量 { ";
    for (int v : psi_tt.order) std::cout << v << " ";
    std::cout << "} → 局部编号 { ";
    for (int i = 1; i <= n_psi; i++) std::cout << i << " ";
    std::cout << "}\n";
    std::cout << string(depth*2, ' ') << "   映射关系：";
    for (int i = 0; i < n_psi; i++)
        std::cout << "位置" << (i+1) << "→变量" << psi_tt.order[i] << " ";
    std::cout << "\n\n";

    // 递归分解 φ 和 ψ
    int L = bi_decomp_recursive(phi_tt, depth + 1);
    int R = bi_decomp_recursive(psi_tt, depth + 1);

    if (L < 0 || R < 0)
    return -1;

    // 创建当前节点（用 F 作为函数）
    return new_node(result.F01, {L, R});
}

// =====================================================
// 顶层调用入口（参照 run_dsd_recursive）
// =====================================================
inline bool run_bi_decomp_recursive(const std::string& binary01)
{
   const bool prev_enable_else_dec = ENABLE_ELSE_DEC;
    const bool prev_minimal_output = BD_MINIMAL_OUTPUT;
    const bool prev_only_k2_eq_0 = BD_ONLY_K2_EQ_0;
    const bool prev_dsd_mix_fallback = BD_ENABLE_DSD_MIX_FALLBACK;
    RESET_NODE_GLOBAL();
    ENABLE_ELSE_DEC = true;
    BD_MINIMAL_OUTPUT = true;

    if (!is_power_of_two(binary01.size())) {
        std::cout << "输入长度必须为 2^n\n";
         BD_MINIMAL_OUTPUT = prev_minimal_output;
        BD_ONLY_K2_EQ_0 = prev_only_k2_eq_0;
        return false;
    }

    int n = static_cast<int>(std::log2(binary01.size()));
    ORIGINAL_VAR_COUNT = n;

    TT root;
    root.f01 = binary01;
    root.order.resize(n);

    for (int i = 0; i < n; ++i)
       // root.order[i] = i + 1;  // 位置 (i+1) 对应变量 (i+1)（原始编号）
       root.order[i] = n - i;  // 位置 (i+1) 对应变量 (n - i)（高位编号大、低位编号小）

    std::cout << "======= 双分解递归开始 =======\n";
    std::cout << "输入 = " << binary01 << " (n=" << n << ")\n";
    std::cout << "初始映射：";
    for (int i = 0; i < n; i++)
        std::cout << "位置" << (i+1) << "→变量" << root.order[i] << " ";
    std::cout << "\n\n";

    NODE_LIST.clear();
    NODE_ID = 1;
    STEP_ID = 1;
    FINAL_VAR_ORDER.clear();
    
        // 预先构建所有输入节点，固定编号与变量一一对应
    for (int v = 1; v <= n; ++v)
        new_in_node(v);

    // 可选：先缩减到 support（这里用和 DSD 相同的 shrink_to_support）
    TT root_shrunk = shrink_to_support(root);
    int root_id = bi_decomp_recursive(root_shrunk, 0);

        if (root_id < 0)
    {
        ENABLE_ELSE_DEC = prev_enable_else_dec;
        BD_MINIMAL_OUTPUT = prev_minimal_output;
        BD_ONLY_K2_EQ_0 = prev_only_k2_eq_0;
        BD_ENABLE_DSD_MIX_FALLBACK = prev_dsd_mix_fallback;
        return false;
    }

    // 打印最终节点列表
    std::cout << "\n===== 最终双分解节点列表 =====\n";
    for (auto& nd : NODE_LIST)
    {
        std::cout << nd.id << " = " << nd.func;

        if (nd.func == "in")
        {
            // 输入节点：显示原始变量编号
            std::cout << "(var=" << nd.var_id << ")";
        }
        else if (!nd.child.empty())
        {
            std::cout << "(";
            for (size_t i = 0; i < nd.child.size(); ++i)
            {
                std::cout << nd.child[i];
                if (i + 1 < nd.child.size())
                    std::cout << ",";
            }
            std::cout << ")";
        }

        std::cout << "\n";
    }

    std::cout << "Root = " << root_id << "\n";

    std::cout << "FINAL_VAR_ORDER = { ";
    for (int v : FINAL_VAR_ORDER) std::cout << v << " ";
    std::cout << "}\n";

    BD_MINIMAL_OUTPUT = prev_minimal_output;
    BD_ONLY_K2_EQ_0 = prev_only_k2_eq_0;
    BD_ENABLE_DSD_MIX_FALLBACK = prev_dsd_mix_fallback;
    return true;
}

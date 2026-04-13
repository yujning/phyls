#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include <numeric>
#include <unordered_map>
#include "node_global.hpp"

// =====================================================
// 工具
// =====================================================
inline uint64_t pow2(int k) { return 1ull << k; }

// =====================================================
// 变量重排（编码映射）
// new_order / old_order: MSB -> LSB（1-based 原变量编号）
// 说明：输入/输出 mf 的存储顺序均为：下标0对应(11..1)，下标N-1对应(00..0)
// =====================================================
inline std::string reorder_tt_by_var_order(
    const std::string& mf,
    int n,
    const std::vector<int>& new_order_msb2lsb,
    const std::vector<int>& old_order_msb2lsb)
{
    const size_t N = mf.size();
    std::string out(N, '0');

    // var -> position in old_order (MSB position = 0)
    std::unordered_map<int, int> pos_in_old;
    pos_in_old.reserve(old_order_msb2lsb.size());
    for (int i = 0; i < n; ++i)
        pos_in_old[old_order_msb2lsb[i]] = i;


    for (size_t bottom_idx = 0; bottom_idx < N; ++bottom_idx)
    {
        // bottom_idx: storage order 11..1 -> 00..0
        // convert to standard index (00..0 -> 11..1)
        uint64_t new_idx_std = (uint64_t)N - 1 - bottom_idx;

        uint64_t old_idx_std = 0;

        // decode bits according to new_order (MSB->LSB)
        for (int i = 0; i < n; ++i)
        {
            int bit = (new_idx_std >> (n - 1 - i)) & 1;
            int var = new_order_msb2lsb[i]; // 1-based

            auto it = pos_in_old.find(var);
            if (it == pos_in_old.end()) continue;
            int pos = it->second;

            old_idx_std |= (uint64_t(bit) << (n - 1 - pos));
        }

        // convert back to storage order index
        uint64_t old_idx_bottom = (uint64_t)N - 1 - old_idx_std;
        out[bottom_idx] = mf[(size_t)old_idx_bottom];
    }

    return out;
}

// =====================================================
// 枚举 (x,y,z)，按 y 从小到大
// 约束：x+y ≤ K, y+z+1 ≤ K, x+y+z = n (K = LUT_PAIR_FANIN_LIMIT)
// =====================================================
inline std::vector<std::tuple<int,int,int>>
enumerate_xyz(int n)
{
    std::vector<std::tuple<int,int,int>> out;
    const int k = LUT_PAIR_FANIN_LIMIT;
    for (int y = 1; y <= std::max(1, k - 2); ++y)
    for (int x = 0; x <= k; ++x)
    {
        int z = n - x - y;
        if (z < 0) continue;
        if (x + y > k) continue;
        if (y + z + 1 > k) continue;
        out.emplace_back(x,y,z);
    }
    return out;
}

// =====================================================
// 枚举变量划分 (A,B,C)，保持 MSB->LSB
// =====================================================
inline std::vector<
    std::tuple<std::vector<int>,std::vector<int>,std::vector<int>>>
enumerate_variable_partitions(const std::vector<int>& vars_msb2lsb,
                              int x, int y, int z)
{
    std::vector<
        std::tuple<std::vector<int>,std::vector<int>,std::vector<int>>> parts;

     std::vector<int> vars = vars_msb2lsb;
    int n = (int)vars.size();

    std::vector<bool> sel_b(n,false);
    std::fill(sel_b.begin(), sel_b.begin()+y, true);

    do {
        std::vector<int> B, rest;
        for (int i=0;i<n;i++)
            (sel_b[i]?B:rest).push_back(vars[i]);

        if (x == 0) {
            parts.emplace_back(std::vector<int>{}, B, rest);
            continue;
        }

        std::vector<bool> sel_a(rest.size(), false);
        std::fill(sel_a.begin(), sel_a.begin()+x, true);

        do {
            std::vector<int> A,C;
            for (size_t i=0;i<rest.size();i++)
                (sel_a[i]?A:C).push_back(rest[i]);
            parts.emplace_back(A,B,C);
        } while (std::prev_permutation(sel_a.begin(), sel_a.end()));

    } while (std::prev_permutation(sel_b.begin(), sel_b.end()));

    return parts;
}

// =====================================================
// A|B|C → MF index (MSB→LSB)
// =====================================================
inline uint64_t mf_index(uint64_t a, uint64_t b, uint64_t c, int y, int z)
{
    return (a << (y + z)) | (b << z) | c;
}
// =====================================================
// 打印重排后的真值表（存储顺序：11...1 → 00...0）
// =====================================================
inline void print_reordered_tt(
    const std::string& mf,
    int n,
    const std::vector<int>& order_msb2lsb)
{
    std::cout << "🧮 Reordered TT (order: ";
    for (size_t i = 0; i < order_msb2lsb.size(); ++i)
    {
        if (i) std::cout << ",";
        std::cout << order_msb2lsb[i];
    }
    //std::cout << " — 11...1→00...0)\n";

    std::string seq;
    seq.reserve(mf.size());

    // for (size_t bottom_idx = 0; bottom_idx < mf.size(); ++bottom_idx)
    // {
    //     uint64_t top_idx = (uint64_t)mf.size() - 1 - bottom_idx;

    //     for (int i = 0; i < n; ++i)
    //     {
    //         int bit = (top_idx >> (n - 1 - i)) & 1;
    //         int var = order_msb2lsb[i];
    //         std::cout << "  v" << var << "=" << bit;
    //     }

    //     std::cout << " -> " << mf[bottom_idx] << "\n";
    //     seq.push_back(mf[bottom_idx]);
    // }

    std::cout << "🧮 TT sequence: " << seq << "\n";
}

// =====================================================
// 打印 Mr / MZ
// =====================================================
inline void print_MZ_delta(int x, int y)
{
    uint64_t Ix = pow2(x);
    uint64_t t  = pow2(y);

    std::cout << "🟩 MZ = I_{2^" << x << "} ⊗ Mr_{2^" << y << "}\n";
    std::cout << "🟩 Mr = δ_" << (t*t) << "[ ";
    for (uint64_t i=0;i<t;i++)
        std::cout << (i*t+i+1) << " ";
    std::cout << "]\n";

    std::cout << "🟩 MZ = δ_" << (Ix*t*t) << "[ ";
    for (uint64_t a=0;a<Ix;a++)
        for (uint64_t i=0;i<t;i++)
            std::cout << (a*t*t + i*t + i + 1) << " ";
    std::cout << "]\n";
}

// =====================================================
// Step 1：MF′ + MZ → MXY
// =====================================================
inline std::string compute_MXYX_from_MF(
    const std::string& MF,
    int x,int y,int z)
{
    uint64_t HN = pow2(x);
    uint64_t MN = pow2(y);
    uint64_t LN = pow2(z);

    uint64_t MXY_cols = pow2(x+2*y+z);
    std::string MXY(MXY_cols,'x');

    for (uint64_t a=0;a<HN;a++)
    for (uint64_t b=0;b<MN;b++)
    {
        uint64_t r = b*MN + b;
        uint64_t block = a*MN*MN + r;

        for (uint64_t c=0;c<LN;c++)
        {
            uint64_t mf_col  = mf_index(a,b,c,y,z);
            uint64_t mxy_col = block*LN + c;
            MXY[mxy_col] = MF[mf_col];
        }
    }
    return MXY;
}

// =====================================================
// Step 2：求 MX / MY
// =====================================================
// =====================================================
// Step 2：求 MX / MY
// 兼容原功能：MX 仍然输出为 0/1（x->0），用于 DAG
// 额外输出：MX_with_x 保留 x，用于打印调试
// =====================================================
inline bool solve_MX_MY_from_MXY(
    const std::string& MXY,
    int x,int y,int z,
    std::string& MX,          // 原输出：用于 DAG（x->0）
    std::string& MX_with_x,   // 新增：用于打印（保留 x）
    std::string& MY)
{
    uint64_t MYN = pow2(x+y);
    uint64_t block_size = pow2(y+z);

    // concrete blocks：会在 merge 时把 x 吃掉（尽可能具体化）
    std::vector<std::string> blocks;
    // symbolic blocks：保留首次出现的 block（含 x），merge 不修改
    std::vector<std::string> blocks_symbolic;

    MY.resize(MYN);

    for (uint64_t i=0;i<MYN;i++)
    {
        std::string blk = MXY.substr(i*block_size, block_size);

        bool hit=false;
        for (size_t k=0;k<blocks.size();k++)
        {
            bool ok=true;
            for (uint64_t j=0;j<block_size;j++)
            {
                if (blk[j]!='x' && blocks[k][j]!='x'
                    && blk[j]!=blocks[k][j]) { ok=false; break; }
            }

            if (ok)
            {
                // 只更新 concrete blocks：把已有的 x 用 blk 填起来
                for (uint64_t j=0;j<block_size;j++)
                    if (blocks[k][j]=='x') blocks[k][j]=blk[j];

                // symbolic 不动，保留 x 的结构
                MY[i]=(k==0?'1':'0');
                hit=true; break;
            }
        }

        if (!hit)
        {
            if (blocks.size()==2) return false;
            blocks.push_back(blk);
            blocks_symbolic.push_back(blk);
            MY[i]=(blocks.size()==1?'1':'0');
        }
    }

    if (blocks.size()!=2) return false;

    // ===== 输出 MX_with_x（保留 x）=====
    MX_with_x.clear();
    for (auto& b: blocks_symbolic)
        for (char c: b)
            MX_with_x.push_back(c);

    // ===== 输出 MX（兼容原逻辑：x->0，用于 DAG）=====
    MX.clear();
    for (auto& b: blocks)
        for (char c: b)
            MX.push_back(c=='x'?'0':c);

    return true;
}

// =====================================================
// children 构造（含占位符）
// =====================================================
inline std::vector<int> make_children_with_placeholder(
    const std::vector<int>& order_msb2lsb,
    int placeholder_var_id,
    int placeholder_node_id)
{
    std::vector<int> ch;
    ch.reserve(order_msb2lsb.size());

    for (int var_id : order_msb2lsb)
    {
        if (var_id == placeholder_var_id) continue;
        new_in_node(var_id);
        if (std::find(FINAL_VAR_ORDER.begin(),
                      FINAL_VAR_ORDER.end(),
                      var_id) == FINAL_VAR_ORDER.end())
            FINAL_VAR_ORDER.push_back(var_id);
    }

    for (auto it = order_msb2lsb.rbegin();
         it != order_msb2lsb.rend(); ++it)
    {
        int var_id = *it;
        if (var_id == placeholder_var_id)
            ch.push_back(placeholder_node_id);
        else
            ch.push_back(new_in_node(var_id));
    }

    return ch;
}

// =====================================================
// ★ 入口函数（最终）
// =====================================================
inline bool run_strong_bi_dec_and_build_dag(const TT& root_tt)
{
    int n = (int)root_tt.order.size();

    int max_var_id = 0;
    for (int v : root_tt.order)
    max_var_id = std::max(max_var_id, v);

    if ((size_t(1) << n) != root_tt.f01.size()) return false;

    // 原始 TT 变量顺序（MSB->LSB）
    std::vector<int> original_order = root_tt.order;

    RESET_NODE_GLOBAL();

    for (auto [x,y,z] : enumerate_xyz(n))
    {
        auto parts = enumerate_variable_partitions(original_order, x, y, z);
        for (const auto& [A,B,C] : parts)
        {
            std::vector<int> new_order;
            new_order.insert(new_order.end(), A.begin(), A.end());
            new_order.insert(new_order.end(), B.begin(), B.end());
            new_order.insert(new_order.end(), C.begin(), C.end());

            std::string MFp =
             reorder_tt_by_var_order(root_tt.f01, n, new_order, original_order);

            print_reordered_tt(MFp, n, new_order);

            std::cout<<"📐 Try "<<x<<"+"<<y<<"+"<<z
                     <<"  A={ "; for(int v:A) std::cout<<v<<" ";
            std::cout<<"} B={ "; for(int v:B) std::cout<<v<<" ";
            std::cout<<"} C={ "; for(int v:C) std::cout<<v<<" ";
            std::cout<<"}\n";

            print_MZ_delta(x,y);

            std::string MXY =
                compute_MXYX_from_MF(MFp,x,y,z);
            std::cout<<"🟨 MXY = "<<MXY<<"\n";


            std::string MX, MX_with_x, MY;
            if (!solve_MX_MY_from_MXY(MXY, x, y, z, MX, MX_with_x, MY))
            {
                std::cout<<"❌ MX/MY unsat\n";
                continue;
            }

            std::cout<<"✅ Strong Bi-Decomposition found\n";
            std::cout<<"🟦 MY = "<<MY<<"\n";
            std::cout<<"🟥 MX = "<<MX<<"\n";
            std::cout<<"🟥 MX(x) = "<<MX_with_x<<"\n";


            // ===== 构造 DAG =====
            ORIGINAL_VAR_COUNT = max_var_id;


            TT tt_my;
            tt_my.f01 = MY;
            tt_my.order.insert(tt_my.order.end(), A.begin(), A.end());
            tt_my.order.insert(tt_my.order.end(), B.begin(), B.end());

            auto children_my = make_children_from_order(tt_my);
            std::reverse(children_my.begin(), children_my.end());
            int my_node = new_node(MY, children_my);

            int placeholder_var_id = max_var_id + 1;

            std::vector<int> order_mx;
            order_mx.push_back(placeholder_var_id);
            order_mx.insert(order_mx.end(), B.begin(), B.end());
            order_mx.insert(order_mx.end(), C.begin(), C.end());

            auto children_mx = make_children_with_placeholder(
                order_mx, placeholder_var_id, my_node);
            std::reverse(children_mx.begin(), children_mx.end());

            new_node(MX, children_mx);
            return true;
        }
    }

    std::cout<<"❌ No valid strong bi-decomposition\n";
    return false;
}

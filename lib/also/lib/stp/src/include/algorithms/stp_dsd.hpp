#pragma once
#include <bits/stdc++.h>
using namespace std;

#include <set>
#include <algorithm>
#include "node_global.hpp"


#include "reorder.hpp"
#include "dsd_else_dec.hpp"
// ================================================
// kitty truth table
// ================================================
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

// ================================================
// DSD Node
// ================================================
struct DSDNode {
    int id;
    string func;      // "in", "0", "1", 或 STP 核心函数串（例如 "0111","1101"...）
    vector<int> child;
    int var_id = -1;  // 对于 "in" 节点：原始变量编号（1-based）
};
static inline bool is_prime_node(int id)
{
    if (id <= 0 || id >= (int)NODE_LIST.size()) return false;
    const auto& nd = NODE_LIST[id];
    if (nd.func == "in" || nd.func == "0" || nd.func == "1") return false;
    return nd.child.size() > 2;
}
static void replace_node_everywhere(int old_id, int new_id)
{
    if (old_id == new_id) return;

    for (auto& nd : NODE_LIST)
        for (auto& c : nd.child)
            if (c == old_id) c = new_id;

    if (ROOT_NODE_ID == old_id)
        ROOT_NODE_ID = new_id;
}
static int refine_prime_node(int node_id)
{
    const auto& nd = NODE_LIST[node_id];
    int n = log2(nd.func.size());

    TT t;
    t.f01 = nd.func;

    // child → order（只用于 Shannon / Exact）
    t.order.clear();
    for (int c : nd.child)
    {
        if (NODE_LIST[c].func == "in")
            t.order.push_back(NODE_LIST[c].var_id);
        else
            t.order.push_back(-1);
    }

    return dsd_else_decompose(t, /*depth=*/0);
}
static void refine_all_prime_nodes()
{
    vector<int> primes;

    for (auto& nd : NODE_LIST)
        if (is_prime_node(nd.id))
            primes.push_back(nd.id);

    if (primes.empty())
    {
        cout << "✅ no prime nodes to refine\n";
        return;
    }

    cout << "🔧 refining " << primes.size() << " prime nodes\n";

    for (int pid : primes)
    {
        if (!is_prime_node(pid)) continue;
        int new_root = refine_prime_node(pid);
        replace_node_everywhere(pid, new_root);
    }
}

inline bool is_terminal_tt(const TT& t)
{
    // 常量
    if (t.order.empty())
        return true;

    // 1 个或 2 个变量 → 都当叶子
    if (t.order.size() <= 2)
        return true;

    return false;
}

inline bool is_binary_constant(const std::string& f01)
{
    if (f01.empty()) return false;
    return std::all_of(f01.begin(), f01.end(),
                       [&](char c){ return c == f01[0]; });
}

// ================================================
// TT = truth table + variable order
// order[i] = 原始变量编号（1-based），对应局部位置 i+1
// ================================================



// ================================================
// make_tt_from01
// ================================================
static kitty::dynamic_truth_table make_tt_from01(const string& f01)
{
    size_t len = f01.size();
    unsigned n = log2(len);
    kitty::dynamic_truth_table tt(n);

    for (uint64_t i = 0; i < len; i++)
        if (f01[i] == '1')
            kitty::set_bit(tt, i);

    return tt;
}



// ================================================
// shrink_to_support（按 support 缩减 TT）
// ================================================
static TT shrink_to_support(const TT& in)
{
    auto tt = make_tt_from01(in.f01);

    vector<int> supp;
    for (int i = 0; i < tt.num_vars(); i++)
        if (kitty::has_var(tt, i))
            supp.push_back(i);

    int n = tt.num_vars();
    
    if (supp.size() == n) {
        cout << "✓ 所有变量都有影响，无需缩减\n\n";
        return in;
    }
    
    cout << "🔍 检测到的有效变量（Kitty位置 → 你的位置 → 原始变量）：\n";
    
    // 按Kitty顺序收集
    vector<pair<int, int>> kitty_order_retained;  // {kitty_pos, 原始变量编号}
    for (int kitty_pos : supp) {
        int your_pos = n - kitty_pos;
        int original_var = in.order[your_pos - 1];
        kitty_order_retained.push_back({kitty_pos, original_var});
        
        cout << "   Kitty位置" << kitty_pos << " → 你的位置" << your_pos 
             << " → 原始变量" << original_var << " ✓\n";
    }
    
    // 按Kitty位置升序排序
    sort(kitty_order_retained.begin(), kitty_order_retained.end(), 
         [](auto& a, auto& b) { return a.first < b.first; });
    
    cout << "\n   按Kitty顺序保留的变量：{ ";
    for (auto& p : kitty_order_retained)
        cout << p.second << " ";
    cout << "}（对应Kitty位置 ";
    for (auto& p : kitty_order_retained)
        cout << p.first << " ";
    cout << "）\n";

    // 构建缩减后的真值表（直接提取对应位）
    unsigned nv = kitty_order_retained.size();
    kitty::dynamic_truth_table new_tt(nv);

    for (uint64_t x = 0; x < (1ull << nv); x++)
    {
        uint64_t old = 0;
        for (int b = 0; b < (int)nv; b++)
        {
            uint64_t bit = (x >> b) & 1;
            int old_kitty_pos = kitty_order_retained[b].first;
            old |= (bit << old_kitty_pos);
        }
        
        if (kitty::get_bit(tt, old))
            kitty::set_bit(new_tt, x);
    }

    TT out;
    out.f01.resize(1ull << nv);
    for (uint64_t i = 0; i < out.f01.size(); i++)
        out.f01[i] = kitty::get_bit(new_tt, i) ? '1' : '0';

    cout << "   缩减后的真值表（Kitty顺序）= " << out.f01 << "\n";

    // 🔥 转换为STP变量顺序（就是把Kitty顺序反过来）
    out.order.clear();
    for (int i = nv - 1; i >= 0; i--) {
        out.order.push_back(kitty_order_retained[i].second);
    }

    cout << "   对应STP变量顺序：{ ";
    for (int v : out.order) cout << v << " ";
    cout << "}（你的位置 ";
    for (int i = 0; i < (int)nv; i++) 
        cout << (n - kitty_order_retained[nv-1-i].first) << " ";
    cout << "）\n";
    cout << "   STP真值表 = " << out.f01 << " （编码不变）\n\n";

    return out;
}

// ================================================
// mul_ui（旧 STP 模板用，保留）
// ui ∈ {"10","01","11","00"}
// ================================================
static inline string mul_ui(const string& ui, const string& w)
{
    if (ui == "10") return w;
    if (ui == "01") {
        string r; r.reserve(w.size());
        for (size_t i = 0; i + 1 < w.size(); i += 2) {
            r.push_back(w[i + 1]);
            r.push_back(w[i]);
        }
        return r;
    }
    if (ui == "11") return string(w.size(), '1');
    if (ui == "00") return string(w.size(), '0');
    return w;
}

// ================================================
// TemplateResult / run_case_once （保留给 s>1 其它 case 用）
// ================================================
struct TemplateResult {
    string MF;
    string Mphi;
    string Mpsi;
};

static TemplateResult run_case_once(
    const vector<string>& blocks,
    int s,
    const string& S0,
    const string& S1)
{
    int m = blocks.size();
    vector<string> W = blocks;

    string MF = S0 + S1;
    string Mpsi;

    auto try_u = [&](const string& u)->bool {
        if (u=="10" || u=="01") {
            vector<string> cand;
            for (int i = 0; i < m; i++) {
                if (is_constant_block(W[i])) continue;
                string c = mul_ui(u, W[i]);
                if (mul_ui(u, c) == W[i])
                    cand.push_back(c);
            }
            if (cand.empty()) return false;
            for (int k=1;k<(int)cand.size();k++)
                if (cand[k] != cand[0]) return false;
            Mpsi = cand[0];
            return true;
        }
        return false;
    };

    if (!try_u(S0) && !try_u(S1))
    {
        int p = -1;
        for (int i = 0; i < m; i++)
            if (!is_constant_block(W[i])) { p=i; break; }
        if (p < 0) p = 0;
        Mpsi = W[p];
    }

    string exp0 = mul_ui(S0, Mpsi);
    string exp1 = mul_ui(S1, Mpsi);

    string Mphi;
    Mphi.reserve(m);
    for (int i = 0; i < m; i++)
        Mphi.push_back(W[i] == exp0 ? '1' : '0');

    return { MF, Mphi, Mpsi };
}

// =====================================================
// 块辅助：判断一个 block 是否全 0 / 全 1 / 常量
// =====================================================
static bool is_all_zero(const string& b)
{
    return !b.empty() && std::all_of(b.begin(), b.end(), [](char c){return c=='0';});
}
static bool is_all_one(const string& b)
{
    return !b.empty() && std::all_of(b.begin(), b.end(), [](char c){return c=='1';});
}
static bool is_constant_block_full(const string& b)
{
    return is_all_zero(b) || is_all_one(b);
}

// =====================================================
// 统一分块语义：
//   - 对 s==1：
//       * uniq.size()==1: B0
//           MF  = B0 + B0  （比如 10 -> 1010, 11 -> 1111）
//           MΦ  = 全 '1'
//           MΨ  = B0
//       * uniq.size()==2: B0, B1
//           MF  = B0 + B1
//           MΦ[i] = '1' 若 blocks[i]==B0 否则 '0'
//           MΨ = "10"
//       * uniq.size()>2: 返回 false
//   - 对 s>1：
//       * uniq.size()==1: B0
//           MF  = B0 + B0
//           MΦ  = 全 '1'
//           MΨ  = B0
//       * uniq.size()==2 且恰好一块常数、一块非常数:
//           设常数块为 C，非常数块为 W：
//           若 C 是全 1 → MF = "1101"
//           若 C 是全 0 → MF = "1110"
//           MΦ[i] = '1' 若 blocks[i]==C, 否则 '0'
//           MΨ = W
//       * 其它情况：返回 false（交给 STP 模板 run_case_once 处理）
// =====================================================



static bool derive_block_semantics_general(
    const vector<string>& blocks,
    int s,
    string &MF,
    string &Mphi,
    string &Mpsi)
{
    vector<string> uniq;
    for (auto &b : blocks)
    {
        if (find(uniq.begin(), uniq.end(), b) == uniq.end())
            uniq.push_back(b);
    }

    if (uniq.empty())
        return false;

    // ---------- s == 1 ----------
    if (s == 1)
    {
        if (uniq.size() > 2)
            return false;

        // if (uniq.size() == 1)
        // {
        //     string B0 = uniq[0];
        //     MF = B0 + B0;
        //     Mphi.assign(blocks.size(), '1');
        //     Mpsi = B0;
        //     return true;
        // }
        if (uniq.size() == 1)
        {
            // 🔥 所有 blocks 相同 ⇒ f = Ψ
            Mpsi = uniq[0];

            MF.clear();     // 用 empty MF 作为“collapse 标记”
            Mphi.clear();   // Φ 无意义

            return true;
        }


        // uniq.size()==2
        {
            string B0 = uniq[0];
            string B1 = uniq[1];
            MF = B0 + B1;

            Mphi.resize(blocks.size());
            for (int i=0;i<(int)blocks.size();i++)
                Mphi[i] = (blocks[i] == B0 ? '1' : '0');

            Mpsi = "10";
            return true;
        }
    }

    // ---------- s > 1 ----------
    // if ((int)uniq.size() == 1)
    // {
    //     string B0 = uniq[0];
    //     MF   = B0 + B0;
    //     Mphi.assign(blocks.size(), '1');
    //     Mpsi = B0;
    //     return true;
    // }

    // ---------- s > 1 ----------
    if ((int)uniq.size() == 1)
    {
        // 🔥 所有 blocks 相同 ⇒ f = Ψ
        Mpsi = uniq[0];

        MF.clear();     // collapse 标记
        Mphi.clear();   // Φ 无意义

        return true;
    }


    if ((int)uniq.size() == 2)
    {
        string U0 = uniq[0];
        string U1 = uniq[1];

        bool U0_const = is_constant_block_full(U0);
        bool U1_const = is_constant_block_full(U1);

        // 恰好一块常数、一块非常数
        if (U0_const ^ U1_const)
        {
            string C = U0_const ? U0 : U1;  // 常数块
            string W = U0_const ? U1 : U0;  // 非常数块

            // 🔥 根据常数块的值决定 MF
            //   常数块全 0 → MF = "0010"  (00=常数, 10=恒等)
            //   常数块全 1 → MF = "1110"  (11=常数, 10=恒等)
            if (is_all_zero(C))
                MF = "0010";
            else if (is_all_one(C))
                MF = "1110";
            else
                return false;

            // MΦ：'1' 表示该块是常数块，'0' 表示非常数块
            Mphi.resize(blocks.size());
            for (int i=0;i<(int)blocks.size();i++)
                Mphi[i] = (blocks[i] == C ? '1' : '0');

            // 🔥 MΨ：直接用非常数块本身（u=10，恒等）
            Mpsi = W;

            return true;
        }

        // 其它 2 块模式（都非常数或者两个不同常数），交给 STP 模板
        return false;
    }

    // uniq.size() > 2
    return false;
}
// =====================================================
// new_node / new_in_node
// =====================================================
// static int new_node(const string& func, const vector<int>& child)
// {
//     int id = NODE_ID++;
//     NODE_LIST.push_back({ id, func, child, -1 });
//     return id;
// }

// 哈希表：func + children → node_id
inline int new_node(const std::string& func, const std::vector<int>& child)
{
    // 结构哈希 key（不 reverse）
    auto key = std::make_tuple(func, child);

    // 如果已有完全相同结构 → 直接复用
    if (NODE_HASH.count(key))
        return NODE_HASH[key];

    // 创建新节点
    int id = NODE_ID++;
    NODE_LIST.push_back({ id, func, child, -1 });

    // 加入哈希表
    NODE_HASH[key] = id;

    return id;
}


// static int new_in_node(int var_id)  // var_id = 1..n
// {
//     int id = NODE_ID++;
//     NODE_LIST.push_back({ id, "in", {}, var_id });
//     return id;
// }

// var_id → node_id


inline int new_in_node(int var_id)
{
    if (INPUT_NODE_CACHE.count(var_id))
        return INPUT_NODE_CACHE[var_id];

    int id = NODE_ID++;
    NODE_LIST.push_back({ id, "in", {}, var_id });

    INPUT_NODE_CACHE[var_id] = id;
    return id;
}

inline std::vector<int> make_children_from_order(const TT& t)
{
    std::vector<int> ch;
    ch.reserve(t.order.size());

        // 先按变量编号升序初始化输入节点，保证 in(var=i) 的编号与变量一致
    std::vector<int> sorted_vars = t.order;
    std::sort(sorted_vars.begin(), sorted_vars.end());
    sorted_vars.erase(std::unique(sorted_vars.begin(), sorted_vars.end()), sorted_vars.end());
    for (int var_id : sorted_vars)
        new_in_node(var_id);

    // FINAL_VAR_ORDER 保持原有次序

        for (int var_id : t.order)
          if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), var_id) == FINAL_VAR_ORDER.end())
            FINAL_VAR_ORDER.push_back(var_id);

    // klut 的 PI 顺序对应 kitty 变量顺序（最低位在前），所以反向绑定
    for (auto it = t.order.rbegin(); it != t.order.rend(); ++it)
    {
        int var_id = *it;
        ch.push_back(new_in_node(var_id));
    }
    
  
    return ch;
}


// =====================================================
// build_small_tree - 记录变量到 FINAL_VAR_ORDER
// =====================================================
static int build_small_tree(const TT& t)
{
    int nv = t.order.size();

    if (nv == 1)
    {
        int var_id = t.order[0];  // 原始变量编号
        int a = new_in_node(var_id);
        
        if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), var_id) 
            == FINAL_VAR_ORDER.end())
        {
            FINAL_VAR_ORDER.push_back(var_id);
        }
        
        if (t.f01 == "10") return a;                  // identity
        if (t.f01 == "01") return new_node("01",{a}); // NOT
        if (t.f01 == "00") return new_node("0",{});   // const 0
        if (t.f01 == "11") return new_node("1",{});   // const 1
        return a;
    }

    if (nv == 2)
    {
        for (int var_id : t.order)
        {
            if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), var_id) 
                == FINAL_VAR_ORDER.end())
            {
                FINAL_VAR_ORDER.push_back(var_id);
            }
        }
        
        int a = new_in_node(t.order[0]);
        int b = new_in_node(t.order[1]);
        return new_node(t.f01, { a, b });
    }
    // ⭐ 修复编号：如果不是 1bit/2bit，就要递归建子树，而不能直接返回空 children
    std::vector<int> child_ids;
    for (int i = 0; i < nv; i++)
    {
        int var_id = t.order[i];
        int leaf = new_in_node(var_id);
        child_ids.push_back(leaf);
        if (!count(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), var_id))
            FINAL_VAR_ORDER.push_back(var_id);
    }

    // 这个节点拥有 nv 个输入，所以 child 列表必须列出所有变量
    return new_node(t.f01, child_ids);

}

// =====================================================
// factor_once_with_reorder_01
//   - TT.order[i] = 变量编号（1-based）
//   - Lambda_j = STP 论文中 j = n-i
//   - MF/MΦ/MΨ：优先用 derive_block_semantics_general，
//     若其返回 false，再用 STP 模板 run_case_once
// =====================================================
static bool factor_once_with_reorder_01(
    const TT& in,
    int depth,
    std::string& MF12,
    TT& phi_tt,
    TT& psi_tt)
{
    const string& bin = in.f01;
    int len = bin.size();
    if (!is_power_of_two(len) || len <= 4)
        return false;

    int n = log2(len);
    int r = n / 2;

    // =====================================================
    // 枚举 s（优先大 s）
    // =====================================================
    for (int s = r; s >= 1; --s)
    {
        vector<bool> v(n);
        fill(v.begin(), v.begin() + s, true);

        do {
            // ---------------- Lambda_j ----------------
            vector<int> Lambda_j;
            for (int i = 0; i < n; i++)
                if (v[i]) Lambda_j.push_back(n - i);

            sort(Lambda_j.begin(), Lambda_j.end());

            // ---------------- 重排 ----------------
            vector<bool> inLam(n+1, false);
            for (int j : Lambda_j) inLam[j] = true;

            vector<int> new_order;
            for (int j = 1; j <= n; j++)
                if (!inLam[j]) new_order.push_back(j);
            for (int j : Lambda_j)
                new_order.push_back(j);

            string reordered =
                reorder_by_index_mapping(bin, n, new_order);

            int cid = theorem33_case_id(reordered, s);
            if (cid == 0) continue;

            // ---------------- 分块 ----------------
            int bl = 1 << s;
            int nb = len / bl;
            vector<string> blocks(nb);
            for (int i = 0; i < nb; i++)
                blocks[i] = reordered.substr(i * bl, bl);

            // ---------------- 生成 MF / Φ / Ψ ----------------
            string MF_use, Mphi_use, Mpsi_use;
            bool ok_block =
                derive_block_semantics_general(blocks, s,
                                               MF_use, Mphi_use, Mpsi_use);

            if (!ok_block)
            {
                // fallback 到 STP 模板
                bool has1 = false, has0 = false;
                for (auto& b : blocks)
                    if (is_constant_block(b))
                        (b[0] == '1' ? has1 : has0) = true;

                vector<pair<string,string>> S_list;
                switch (cid) {
                    case 1: S_list = {{"11","00"},{"00","11"}}; break;
                    case 2:
                        if (has1)
                            S_list = {{"11","10"},{"11","01"},{"10","11"},{"01","11"}};
                        else
                            S_list = {{"00","10"},{"00","01"},{"10","00"},{"01","00"}};
                        break;
                    case 3: S_list = {{"10","10"},{"01","01"}}; break;
                    case 4: S_list = {{"10","01"},{"01","10"}}; break;
                    default: continue;
                }

                auto R = run_case_once(blocks, s,
                                       S_list[0].first,
                                       S_list[0].second);
                MF_use   = R.MF;
                Mphi_use = R.Mphi;
                Mpsi_use = R.Mpsi;
            }

            // =================================================
            // 变量顺序恢复（⚠️ 先算 order，再做 collapse）
            // =================================================
            int n_phi = n - s;

            vector<bool> inLam_j(n + 1, false);
            for (int j : Lambda_j) inLam_j[j] = true;

            vector<int> Omega_j;
            for (int j = 1; j <= n; j++)
                if (!inLam_j[j]) Omega_j.push_back(j);

            vector<int> newPos_j = Omega_j;
            newPos_j.insert(newPos_j.end(),
                            Lambda_j.begin(),
                            Lambda_j.end());

            vector<int> newOrder_original;
            for (int j : newPos_j)
                newOrder_original.push_back(in.order[j - 1]);

            vector<int> phi_order_original(
                newOrder_original.begin(),
                newOrder_original.begin() + n_phi
            );

            vector<int> psi_order_original(
                newOrder_original.begin() + n_phi,
                newOrder_original.end()
            );

            // ====== 原有打印（一行不删） ======
            cout << STEP_ID++ << ". MF = [" << MF_use << "]\n";
            cout << "   MΦ = [" << Mphi_use << "]\n";
            cout << "   MΨ = [" << Mpsi_use << "]\n";
            // ... 后面所有你已有的 cout
            cout << "   Φ 原始变量 = { ";
            for (int v : phi_order_original) cout << v << " ";
            cout << "}  Ψ 原始变量 = { ";
            for (int v : psi_order_original) cout << v << " ";
            cout << "}\n\n";

            // =================================================
            // 🔥 在【打印之后】加 collapse 判断
            // =================================================
            if (ok_block && MF_use.empty())
            {
                cout << "   ⚠️ collapse: all blocks identical → f = Ψ\n";
                cout << "   ⚠️ collapse before shrink: Ψ = " << Mpsi_use << "\n";

                psi_tt.f01   = Mpsi_use;
                psi_tt.order = psi_order_original;

                psi_tt = shrink_to_support(psi_tt);

                cout << "   ⚠️ collapse after shrink: Ψ = "
                    << psi_tt.f01 << " vars = { ";
                for (int v : psi_tt.order) cout << v << " ";
                cout << "}\n\n";

                MF12.clear();
                phi_tt = TT{};

                return true;
            }


            // =================================================
            // 正常二叉 DSD
            // =================================================
            MF12          = MF_use;
            phi_tt.f01    = Mphi_use;
            psi_tt.f01    = Mpsi_use;
            phi_tt.order  = phi_order_original;
            psi_tt.order  = psi_order_original;

            return true;

        } while (prev_permutation(v.begin(), v.end()));
    }

    return false;
}

// dsd_factor - 递归 DSD 分解
// =====================================================
static int dsd_factor(const TT& f, int depth = 0)
{
    // =========================
    // 0) 常量
    // =========================
    if (is_binary_constant(f.f01))
        return build_small_tree(f);

    // =========================
    // 1) 1~2 输入：终点
    // =========================
    if (f.order.size() <= 2)
        return build_small_tree(f);

    // =========================
    // 2) 尝试 DSD
    // =========================
    std::string MF12;
    TT phi_tt, psi_tt;

    if (factor_once_with_reorder_01(f, depth, MF12, phi_tt, psi_tt))
    {
        // collapse: f = Ψ（继续 DSD 主线）
        if (MF12.empty())
            return dsd_factor(psi_tt, depth + 1);

        int L, R;

        // Φ
        if (is_terminal_tt(phi_tt))
            L = build_small_tree(phi_tt);
        else
            L = dsd_factor(phi_tt, depth + 1);

        // Ψ
        if (is_terminal_tt(psi_tt))
            R = build_small_tree(psi_tt);
        else
            R = dsd_factor(psi_tt, depth + 1);

        return new_node(MF12, {L, R});
    }

    // =========================
    // 3) DSD 失败：prime 节点处理
    //    - n > 4 : Shannon 1 层（dsd_else_decompose 内部已做）
    //    - n <=4 : exact 2-LUT（dsd_else_decompose 内部已做）
    // =========================
    if (ENABLE_ELSE_DEC)
        return dsd_else_decompose(f, depth);

    // 没开 -e：就把 prime 当叶子（保持原行为）
    return build_small_tree(f);
}

// =====================================================
// run_dsd_recursive
// =====================================================
inline int run_dsd_recursive(const std::string& binary01, bool enable_else_dec)

{
    RESET_NODE_GLOBAL();
     ENABLE_ELSE_DEC = enable_else_dec;
    if (!is_power_of_two(binary01.size())) {
        std::cout << "输入长度必须为 2^n\n";
        return false;
    }

    int n = static_cast<int>(std::log2(binary01.size()));
    ORIGINAL_VAR_COUNT = n;
    
    TT root;
    root.f01 = binary01;
    root.order.resize(n);

    for (int i = 0; i < n; ++i)
        root.order[i] = n - i;   // 6,5,4,3,2,1


    std::cout << "输入 = " << binary01 << " (n=" << n << ")\n";
    std::cout << "初始映射：";
    for (int i = 0; i < n; i++)
        std::cout << "位置" << (i+1) << "→变量" << root.order[i] << " ";
    std::cout << "\n\n";

    NODE_LIST.clear();
    NODE_ID = 1;
    STEP_ID = 1;
    FINAL_VAR_ORDER.clear();

        // 先创建所有输入节点，确保编号与变量一致
    for (int v = 1; v <= n; ++v)
        new_in_node(v);
        
    // 🔥 只在最开始缩减一次
    TT root_shrunk = shrink_to_support(root);
    int root_id = dsd_factor(root_shrunk);
        if (ENABLE_ELSE_DEC)
            refine_all_prime_nodes(); // 递归中不再缩减
    // int root_id = dsd_factor(root);

    // ================= 修改后的这块 =================
    std::cout << "===== 最终 DSD 节点列表 =====\n";
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
            // 任意个子节点：全部打印出来
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
    // ================= 修改结束 =================
    ROOT_NODE_ID = root_id;
    std::cout << "Root = " << root_id << "\n";

    std::cout << "FINAL_VAR_ORDER = { ";
    for (int v : FINAL_VAR_ORDER) std::cout << v << " ";
    std::cout << "}\n";

    return root_id;
}
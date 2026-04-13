#pragma once
#include <bits/stdc++.h>
using std::string;
using std::vector;
using std::set;
using std::cout;
using std::endl;

int run_dsd_recursive(const string& binary01, bool enable_else_dec = false);

//-----------------------------------------
// 判断是否为 2 的幂
//-----------------------------------------
static inline bool is_power_of_two(size_t x){
    return x && ((x & (x - 1)) == 0);
}

//-----------------------------------------
// 块相关工具
//-----------------------------------------
static bool is_constant_block(const std::string& b){
    return std::all_of(b.begin(), b.end(),
                       [&](char c){ return c == b.front(); });
}

static std::string complement_block(const std::string& b){
    std::string r = b;
    for(char& c: r) c = (c=='0' ? '1' : '0');
    return r;
}

//-----------------------------------------
// 定理 3.3 CASE 判断（完全原样）
//-----------------------------------------
static int theorem33_case_id(const std::string& binary, int s){
    const int n = std::log2(binary.size());
    if(!is_power_of_two(binary.size()) || (1<<n)!=(int)binary.size()) return 0;
    if(s<1 || s>n/2) return 0;

    const int block_len = 1<<s;
    const int num_blocks = binary.size() / block_len;

    std::vector<std::string> blocks;
    for(int i=0;i<num_blocks;++i)
        blocks.emplace_back(binary.substr(i*block_len, block_len));

    std::set<std::string> uniq_const, uniq_nonconst;
    for(auto &b: blocks)
        (is_constant_block(b) ? uniq_const : uniq_nonconst).insert(b);

    if(uniq_nonconst.empty() && uniq_const.size()==2) return 1;
    if(uniq_const.size()==1 && uniq_nonconst.size()==1) return 2;
    if(uniq_const.empty() && uniq_nonconst.size()==1) return 3;
    if(uniq_const.empty() && uniq_nonconst.size()==2){
        auto it = uniq_nonconst.begin();
        string a=*it++; string b=*it;
        if(complement_block(a)==b) return 4;
    }
    if(uniq_nonconst.empty() && uniq_const.size()==1) return 5;
    return 0;
}

// =====================================================
// ★ 核心：变量索引 → 编码映射重排
// new_order: 1-based, MSB -> LSB
// =====================================================
static inline string reorder_by_index_mapping(
    const string& binary,
    int n,
    const vector<int>& new_order)
{
    string out(binary.size(), '0');

    for (size_t new_idx = 0; new_idx < binary.size(); ++new_idx)
    {
        uint64_t old_idx = 0;

        for (int i = 0; i < n; ++i)
        {
            int bit = (new_idx >> (n - 1 - i)) & 1;
            int var = new_order[i];          // 1-based
            int pos = var - 1;               // 原始 MSB 位置
            old_idx |= (uint64_t(bit) << (n - 1 - pos));
        }

        out[new_idx] = binary[old_idx];
    }

    return out;
}

// =====================================================
// 重排主函数（等价原 STP 版本）
// =====================================================
inline void all_reorders(const string &binary)
{
    int len = binary.size();
    if(!is_power_of_two(len)){
        cout<<"输入长度必须是 2 的整数次幂\n";
        return;
    }

    int n = log2(len);
    int r = n / 2;

    for(int s=1; s<=r; ++s)
    {
        vector<bool> v(n);
        fill(v.begin(), v.begin()+s, true);

        do{
            vector<int> Lambda;
            for(int i=0;i<n;i++)
                if(v[i]) Lambda.push_back(i+1);

            // ===== 索引映射重排（Omega ⧺ Lambda）=====
            vector<int> new_order;
            vector<bool> inLam(n+1,false);
            for(int j:Lambda) inLam[j]=true;

            for(int j=1;j<=n;j++) if(!inLam[j]) new_order.push_back(j);
            for(int j:Lambda) new_order.push_back(j);

            string reordered =
                reorder_by_index_mapping(binary, n, new_order);

            int cid = theorem33_case_id(reordered, s);
            if(cid!=0){
                cout << "\n===== 重排命中：s="<<s<<" 情形("<<cid<<") =====\n";
                cout << "Λ = { ";
                for(int j : Lambda) cout<<j<<" ";
                cout << "}  => reordered: " << reordered << "\n";

                run_dsd_recursive(reordered, ENABLE_ELSE_DEC);
                return;
            }

        }while(prev_permutation(v.begin(),v.end()));
    }

    cout << "❌ 所有重排均未命中任何分解模式\n";
}

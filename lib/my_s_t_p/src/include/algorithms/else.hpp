// =====================================================
// 枚举所有 (k1,k2,k3) + 所有 Γ,Θ,Λ 重排（找到全部解）
// =====================================================
// static vector<BiDecompResult>
// enumerate_bi_decomposition_all_permutations(const TT& in)
// {
//     vector<BiDecompResult> results;

//     const string &f01 = in.f01;
//     if (f01.empty()) return results;

//     int n = (int)std::log2((double)f01.size());
//     if ((int)in.order.size() != n) return results;

//     // 枚举 k2 和 k3 的大小
//     for (int k2 = 0; k2 <= n - 2; ++k2)
//     {
//         int max_k3 = (n - k2) / 2;

//         for (int k3 = max_k3; k3 >= 1; --k3)
//         {
//             int k1 = n - k2 - k3;
//             if (k1 <= 0) continue;

//             std::cout << "\n========== 枚举 k1=" << k1 << ", k2=" << k2 << ", k3=" << k3 << " ==========\n";

//             // 枚举 Θ 的所有 C(n, k2) 种组合（按“位置”选）
//             vector<bool> theta_mask(n, false);
//             std::fill(theta_mask.begin(), theta_mask.begin() + k2, true);

//             do {
//                 vector<int> Theta_indices; // 1-based positions
//                 for (int i = 0; i < n; ++i)
//                     if (theta_mask[i])
//                         Theta_indices.push_back(i + 1);

//                 // 剩余位置用于 Γ 和 Λ
//                 vector<int> remaining;
//                 for (int i = 0; i < n; ++i)
//                     if (!theta_mask[i])
//                         remaining.push_back(i + 1);

//                 // 枚举 Λ 的所有 C(n-k2, k3) 种组合（位置）
//                 vector<bool> lambda_mask(remaining.size(), false);
//                 std::fill(lambda_mask.begin(), lambda_mask.begin() + k3, true);

//                 do {
//                     vector<int> Lambda_indices; // positions
//                     vector<int> Gamma_indices;  // positions

//                     for (size_t i = 0; i < remaining.size(); ++i)
//                     {
//                         if (lambda_mask[i])
//                             Lambda_indices.push_back(remaining[i]);
//                         else
//                             Gamma_indices.push_back(remaining[i]);
//                     }

//                     // 避免对称重复：当 k1 == k3 时，要求 Γ[0] < Λ[0]
//                     if (k1 == k3 && Gamma_indices[0] > Lambda_indices[0])
//                         continue;

//                     // 打印当前尝试的“位置分组”
//                     std::cout << "  尝试 位置 Γ={";
//                     for (int g : Gamma_indices) std::cout << g << " ";
//                     std::cout << "}, Θ={";
//                     for (int t : Theta_indices) std::cout << t << " ";
//                     std::cout << "}, Λ={";
//                     for (int l : Lambda_indices) std::cout << l << " ";
//                     std::cout << "}\n";

//                     // 根据“位置”对 f01 作重排
//                     // string reordered_f01 = apply_variable_reordering(
//                     //     f01, n, Gamma_indices, Theta_indices, Lambda_indices, k1, k2, k3);

//                     string reordered_f01 = apply_variable_reordering_swap(
//                         f01, n,
//                         Gamma_indices, Theta_indices, Lambda_indices,
//                         k1, k2, k3
//                     );


//                     std::cout << "📌 重排后的 f01（二进制） = " << reordered_f01 << "\n";

//                     // 构造重排后的 TT，order 里必须是“原始编号”，顺序为 [Γ,Θ,Λ]
//                     TT reordered_tt;
//                     reordered_tt.f01 = reordered_f01;
//                     reordered_tt.order.clear();

//                     // ★★ 这里是关键修正点：用位置去 in.order 里取“原始编号”
//                     for (int pos : Gamma_indices)
//                         reordered_tt.order.push_back(in.order[pos - 1]);
//                     for (int pos : Theta_indices)
//                         reordered_tt.order.push_back(in.order[pos - 1]);
//                     for (int pos : Lambda_indices)
//                         reordered_tt.order.push_back(in.order[pos - 1]);

//                     // 在重排后的真值表上尝试分解
//                     auto sub = enumerate_one_case(reordered_tt, k1, k2, k3);

//                     if (!sub.empty())
//                         std::cout << "    ✓ 找到解！\n";

//                     results.insert(results.end(), sub.begin(), sub.end());

//                 } while (std::prev_permutation(lambda_mask.begin(), lambda_mask.end()));

//             } while (std::prev_permutation(theta_mask.begin(), theta_mask.end()));
//         }
//     }

//     return results;
// }

// =====================================================
// 不重排变量版本：枚举所有 (k1,k2,k3) 在当前 in 上分解
// =====================================================
// static vector<BiDecompResult>
// enumerate_bi_decomposition_no_reorder(const TT& in)
// {
//     vector<BiDecompResult> results;

//     const string &f01 = in.f01;
//     if (f01.empty()) return results;

//     int n = (int)std::log2((double)f01.size());
//     if ((int)in.order.size() != n) return results;

//     for (int k2 = 1; k2 <= n - 2; ++k2)
//     {
//         int r = (n - k2) / 2;

//         for (int k3 = 1; k3 <= r; ++k3)
//         {
//             int k1 = n - k2 - k3;
//             if (k1 <= 0) continue;

//             auto sub = enumerate_one_case(in, k1, k2, k3);
//             results.insert(results.end(), sub.begin(), sub.end());
//         }
//     }

//     return results;
// }

// =====================================================
// 根据公式(34)进行真值表重排
//   Gamma_indices / Theta_indices / Lambda_indices 是位置（1-based）
//   真值表按“位置”排列，变量编号始终由 TT.order 保存
// =====================================================
// static string apply_variable_reordering(
//     const string& f01,
//     int n,
//     const vector<int>& Gamma_indices,  // 1-based positions
//     const vector<int>& Theta_indices,  // 1-based positions
//     const vector<int>& Lambda_indices, // 1-based positions
//     int k1, int k2, int k3)
// {
//     // 将 f01 转换为向量形式
//     vector<stp_data> Mf = binary_to_vec(f01);

//     // 构造交换矩阵链（从右到左按公式34的顺序）
//     vector<vector<stp_data>> swap_chain;

//     // 第一部分：W[2^k1, 2^k2]
//     swap_chain.push_back(generate_swap_vec(std::pow(2, k1), std::pow(2, k2)));

//     // 第二部分：⊗_{i=k2}^1 W[2, 2^{j_i+(k2-1)-i}]
//     for (int i = k2; i >= 1; --i)
//     {
//         int j_i = Theta_indices[i - 1];  // 位置
//         int exp = j_i + (k2 - 1) - i;
//         swap_chain.push_back(generate_swap_vec(2, std::pow(2, exp)));
//     }

//     // 第三部分：W[2^{k1+k2}, 2^k3]
//     swap_chain.push_back(generate_swap_vec(std::pow(2, k1 + k2), std::pow(2, k3)));

//     // 第四部分：⊗_{i=k3}^1 W[2, 2^{j_i+(k3-1)-i}]
//     for (int i = k3; i >= 1; --i)
//     {
//         int j_i = Lambda_indices[i - 1];  // 位置
//         int exp = j_i + (k3 - 1) - i;
//         swap_chain.push_back(generate_swap_vec(2, std::pow(2, exp)));
//     }

//     // 矩阵链乘法
//     vector<stp_data> Mperm = Vec_chain_multiply(swap_chain, false);
//     vector<stp_data> result = Vec_semi_tensor_product(Mf, Mperm);

//     // 转换回字符串（跳过第一个维度元素）
//     string reordered;
//     for (size_t i = 1; i < result.size(); ++i)
//         reordered.push_back(result[i] ? '1' : '0');

//     return reordered;
// }


// static bool
// find_first_bi_decomposition(const TT& in, BiDecompResult& out)
// {
//     const string &f01 = in.f01;
//     if (f01.empty()) return false;

//     int n = (int)std::log2((double)f01.size());
//     if ((int)in.order.size() != n) return false;

//     // =====================================================
//     // 新的搜索策略：
//     // 外层循环：偏移量 offset 从 0 开始递增
//     // 内层循环：k2 从 0 到 n-2
//     //   对于每个 k2，计算 k3 = floor((n-k2)/2) - offset
//     //   要求 k3 >= 1
//     // =====================================================
    
//     // 最大可能的 offset
//     int max_offset = n/2;  // 最坏情况 k3 减小到 0
    
//     for (int offset = 0; offset <= max_offset; ++offset)
//     {
//         std::cout << "\n================ 第 " << offset+1 << " 轮搜索 ================\n";
        
//         // k2 从小到大
//         for (int k2 = 0; k2 <= n - 2; ++k2)
//         {
//             // 计算当前偏移下的 k3
//             int k3 = (n - k2) / 2 - offset;
//             if (k3 < 1) continue;  // k3 必须至少为 1
            
//             int k1 = n - k2 - k3;
//             if (k1 <= 0) continue;

//             std::cout << "\n  尝试 k1=" << k1 
//                       << ", k2=" << k2 
//                       << ", k3=" << k3 << " (offset=" << offset << ")\n";

//             // 先试试不重排的情况（变量已经是 [Γ,Θ,Λ] 顺序）
//             auto sub = enumerate_one_case(in, k1, k2, k3);
//             if (!sub.empty())
//             {
//                 out = sub[0];
//                 std::cout << "  ✓ 不需重排即可分解！\n";
//                 return true;
//             }

//             // 枚举 Θ 的所有 C(n, k2) 种组合（按位置，1-based）
//             vector<bool> theta_mask(n, false);
//             if (k2 > 0)
//                 std::fill(theta_mask.begin(), theta_mask.begin() + k2, true);
            
//             // 如果 k2=0，只有一种情况（Θ为空）
//             bool has_theta_combination = (k2 == 0);
            
//             do {
//                 vector<int> Theta_pos;  // Θ的位置
//                 for (int i = 0; i < n; ++i)
//                     if (theta_mask[i])
//                         Theta_pos.push_back(i + 1);
                
//                 // 处理 k2=0 的特殊情况
//                 if (k2 == 0 && !has_theta_combination)
//                 {
//                     Theta_pos.clear();
//                     has_theta_combination = true;
//                 }

//                 // 剩余位置用于分配 Γ 和 Λ
//                 vector<int> remaining_pos;
//                 for (int i = 0; i < n; ++i)
//                     if (!theta_mask[i])
//                         remaining_pos.push_back(i + 1);

//                 // 枚举 Λ 的所有 C(n-k2, k3) 种组合
//                 vector<bool> lambda_mask(remaining_pos.size(), false);
//                 std::fill(lambda_mask.begin(), lambda_mask.begin() + k3, true);
                
//                 bool has_lambda_combination = false;

//                 do {
//                     vector<int> Lambda_pos;  // Λ的位置
//                     vector<int> Gamma_pos;   // Γ的位置

//                     for (size_t i = 0; i < remaining_pos.size(); ++i)
//                     {
//                         if (lambda_mask[i])
//                             Lambda_pos.push_back(remaining_pos[i]);
//                         else
//                             Gamma_pos.push_back(remaining_pos[i]);
//                     }

//                     // ⭐ 对称性剪枝：当 k1 == k3 时，要求 Γ 的首位置 < Λ 的首位置
//                     if (k1 == k3 && Gamma_pos[0] > Lambda_pos[0])
//                         continue;

//                     // 打印当前尝试
//                     std::cout << "    尝试位置：Γ={";
//                     for (int p : Gamma_pos) std::cout << p << " ";
//                     std::cout << "}, Θ={";
//                     for (int p : Theta_pos) std::cout << p << " ";
//                     std::cout << "}, Λ={";
//                     for (int p : Lambda_pos) std::cout << p << " ";
//                     std::cout << "}\n";

//                     // ⭐ 重排真值表：按 [Γ, Θ, Λ] 的位置顺序
//                     string reordered_f01 = apply_variable_reordering_swap(
//                         f01, n,
//                         Gamma_pos, Theta_pos, Lambda_pos,
//                         k1, k2, k3
//                     );

//                     std::cout << "    📌 重排后的 f01 = " << reordered_f01 << "\n";

//                     // 构造重排后的 TT，order 保存原始变量编号
//                     TT reordered_tt;
//                     reordered_tt.f01 = reordered_f01;
//                     reordered_tt.order.clear();

//                     // 按 [Γ, Θ, Λ] 顺序记录原始变量编号
//                     for (int pos : Gamma_pos)
//                         reordered_tt.order.push_back(in.order[pos - 1]);
//                     for (int pos : Theta_pos)
//                         reordered_tt.order.push_back(in.order[pos - 1]);
//                     for (int pos : Lambda_pos)
//                         reordered_tt.order.push_back(in.order[pos - 1]);

//                     // 在重排后的真值表上尝试分解
//                     sub = enumerate_one_case(reordered_tt, k1, k2, k3);

//                     if (!sub.empty())
//                     {
//                         out = sub[0];
//                         std::cout << "      ✓ 找到分解！\n";
//                         return true;
//                     }
                    
//                     has_lambda_combination = true;

//                 } while (std::prev_permutation(lambda_mask.begin(), lambda_mask.end()));

//             } while (k2 > 0 && std::prev_permutation(theta_mask.begin(), theta_mask.end()));
//         }
//     }

//     std::cout << "❌ 遍历所有 (k1,k2,k3) 和变量分组，未找到有效分解\n";
//     return false;
// }
#pragma once
#include <algorithm>
#include <map>

#include <string>
#include <tuple>
#include <vector>
#include <unordered_map>

// 前向声明，避免循环依赖
struct DSDNode;
struct TT {
    std::string f01;
    std::vector<int> order;
};
int new_node(const std::string&, const std::vector<int>&);
int new_in_node(int var_id);
// C++17 允许 inline 全局变量，只需定义一次即可，全工程自动共享
inline int NODE_ID = 1;
inline int STEP_ID = 1;
inline bool STRONG_DSD_DEBUG_PRINT = false;
inline bool ENABLE_ELSE_DEC = false;

inline int ORIGINAL_VAR_COUNT = 0;
inline int NEXT_PLACEHOLDER_ID = 0;
inline int ROOT_NODE_ID = 0;
inline std::unordered_map<int, int> PLACEHOLDER_BINDINGS;
inline std::vector<int> FINAL_VAR_ORDER;

inline std::vector<DSDNode> NODE_LIST;

inline std::map<int, int> INPUT_NODE_CACHE;

inline std::map<std::tuple<std::string, std::vector<int>>, int> NODE_HASH;

// ======================================================
// 🔥 重置所有全局节点状态（每次运行 bd/dsd 前必须调用）
// ======================================================
inline void RESET_NODE_GLOBAL()
{
    NODE_ID = 1;
    STEP_ID = 1;
    ENABLE_ELSE_DEC = false;
    ORIGINAL_VAR_COUNT = 0;
    NEXT_PLACEHOLDER_ID = 0;
    ROOT_NODE_ID = 0;
    PLACEHOLDER_BINDINGS.clear();

    NODE_LIST.clear();
    FINAL_VAR_ORDER.clear();

    INPUT_NODE_CACHE.clear();
    NODE_HASH.clear();
}


inline int allocate_placeholder_var_id(const std::unordered_map<int, int>* existing)
{
    int max_var_id = ORIGINAL_VAR_COUNT;
    if (existing)
    {
        for (const auto& [var_id, _] : *existing)
            max_var_id = std::max(max_var_id, var_id);
    }

    for (int var_id : FINAL_VAR_ORDER)
        max_var_id = std::max(max_var_id, var_id);

    for (const auto& [var_id, _] : PLACEHOLDER_BINDINGS)
        max_var_id = std::max(max_var_id, var_id);

    if (NEXT_PLACEHOLDER_ID <= max_var_id)
        NEXT_PLACEHOLDER_ID = max_var_id + 1;

    return NEXT_PLACEHOLDER_ID++;
}

inline void register_placeholder_binding(int var_id, int node_id)
{
    PLACEHOLDER_BINDINGS[var_id] = node_id;
}

inline void register_placeholder_bindings(const std::unordered_map<int, int>& bindings)
{
    for (const auto& [var_id, node_id] : bindings)
        PLACEHOLDER_BINDINGS[var_id] = node_id;
}
inline std::vector<int> make_children_from_order(const TT& t);

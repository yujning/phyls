#ifndef LUT_FUNC_CACHE_HPP
#define LUT_FUNC_CACHE_HPP

#include <string>
#include <vector>
#include <map>
#include <cstdint>

namespace alice
{

/*============================================================*
 * LUT 功能 key（用于判重）
 *============================================================*/
struct LutFuncKey
{
    uint32_t nvars;
    std::string truth01;
    uint32_t mode = 0;

    bool operator<(const LutFuncKey& rhs) const
    {
        if (nvars != rhs.nvars) return nvars < rhs.nvars;
        if (truth01 != rhs.truth01) return truth01 < rhs.truth01;
        return mode < rhs.mode;
    }
};

/*============================================================*
 * 与算法无关的 LUT DAG 节点描述
 *============================================================*/
struct CachedLutNode
{
    int id = 0;

    // "in" 或 truth table（二进制字符串）
    std::string func;

    // var_id 仅对 input 节点有效（1-based）
    int var_id = 0;

    // 子节点 id
    std::vector<int> child;
};

/*============================================================*
 * 缓存条目（完整功能 DAG）
 *============================================================*/
struct CachedResyn
{
    std::vector<CachedLutNode> nodes;
    int root_id = 0;
};

/*============================================================*
 * LUT 功能 cache（header-only）
 *============================================================*/
class LutFuncCache
{
public:
    static bool has(const LutFuncKey& key)
    {
        return cache().find(key) != cache().end();
    }

    static CachedResyn& get(const LutFuncKey& key)
    {
        return cache().at(key);
    }

    static void insert(const LutFuncKey& key, CachedResyn&& val)
    {
        cache()[key] = std::move(val);
    }

    static void clear()
    {
        cache().clear();
    }

private:
    static std::map<LutFuncKey, CachedResyn>& cache()
    {
        static std::map<LutFuncKey, CachedResyn> inst;
        return inst;
    }
};

} // namespace alice

#endif

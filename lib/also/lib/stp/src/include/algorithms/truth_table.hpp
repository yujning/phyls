// #pragma once

// #include <string>
// #include <sstream>
// #include <stdexcept>
// #include <kitty/dynamic_truth_table.hpp>
// #include <kitty/constructors.hpp>
// #include <kitty/print.hpp>

// namespace alice
// {

// inline std::string bin_to_hex(const std::string& bin)
// {
//     std::string b = bin;

//     while (b.size() % 4 != 0)
//         b = "0" + b;

//     static const char* hex_map = "0123456789abcdef";
//     std::string hex;
//     hex.reserve(b.size() / 4);

//     for (size_t i = 0; i < b.size(); i += 4)
//     {
//         int v = (b[i] - '0') * 8 +
//                 (b[i+1] - '0') * 4 +
//                 (b[i+2] - '0') * 2 +
//                 (b[i+3] - '0');
//         hex.push_back(hex_map[v]);
//     }

//     while (hex.size() > 1 && hex[0] == '0')
//         hex.erase(hex.begin());

//     return hex;
// }

// inline std::string hex_to_binary(const std::string& hex)
// {
//     std::string clean = hex;
//     if (clean.rfind("0x", 0) == 0 || clean.rfind("0X", 0) == 0)
//         clean = clean.substr(2);

//     if (clean.empty())
//         throw std::runtime_error("empty hex string");

//     unsigned bit_count = clean.size() * 4;
//     unsigned num_vars = 0;
//     while ((1u << num_vars) < bit_count)
//         num_vars++;

//     if ((1u << num_vars) != bit_count)
//         throw std::runtime_error("hex length is not 2^n bits");

//     kitty::dynamic_truth_table tt(num_vars);
//     kitty::create_from_hex_string(tt, clean);

//     std::ostringstream oss;
//     kitty::print_binary(tt, oss);
//     return oss.str();
// }

// } // namespace alice


#pragma once

#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>

namespace alice
{

// =====================================================
// FIXED: binary truth table -> hex string
// - 保留前导 0（LUT / 真值表语义必须）
// - 不再“数值化”处理
// =====================================================
inline std::string bin_to_hex(const std::string& bin)
{
    if (bin.empty())
        throw std::runtime_error("bin_to_hex: empty binary string");

    // 确保是 0/1
    for (char c : bin)
        if (c != '0' && c != '1')
            throw std::runtime_error("bin_to_hex: invalid binary character");

    std::string b = bin;

    // pad LEFT to multiple of 4 bits
    while (b.size() % 4 != 0)
        b.insert(b.begin(), '0');

    static const char* hex_map = "0123456789abcdef";
    std::string hex;
    hex.reserve(b.size() / 4);

    for (size_t i = 0; i < b.size(); i += 4)
    {
        int v =
            (b[i]   - '0') * 8 +
            (b[i+1] - '0') * 4 +
            (b[i+2] - '0') * 2 +
            (b[i+3] - '0');

        hex.push_back(hex_map[v]);
    }

    // ★ 关键修复：不再删除前导 0
    return hex;
}

// =====================================================
// FIXED: hex string -> binary truth table
// - 支持 0x / 0X
// - 支持大小写
// - 严格保证 bit 数是 2^n
// =====================================================
inline std::string hex_to_binary(const std::string& hex)
{
    std::string clean = hex;

    // 去掉 0x / 0X
    if (clean.size() >= 2 && clean[0] == '0' &&
        (clean[1] == 'x' || clean[1] == 'X'))
    {
        clean = clean.substr(2);
    }

    if (clean.empty())
        throw std::runtime_error("hex_to_binary: empty hex string");

    // 统一成小写，kitty 更稳
    std::transform(clean.begin(), clean.end(), clean.begin(), ::tolower);

    // 检查合法 hex
    for (char c : clean)
    {
        if (!((c >= '0' && c <= '9') ||
              (c >= 'a' && c <= 'f')))
        {
            throw std::runtime_error("hex_to_binary: invalid hex character");
        }
    }

    const unsigned bit_count = static_cast<unsigned>(clean.size() * 4);

    // bit_count 必须是 2^n
    unsigned num_vars = 0;
    while ((1u << num_vars) < bit_count)
        ++num_vars;

    if ((1u << num_vars) != bit_count)
        throw std::runtime_error("hex_to_binary: bit length is not 2^n");

    kitty::dynamic_truth_table tt(num_vars);
    kitty::create_from_hex_string(tt, clean);

    std::ostringstream oss;
    kitty::print_binary(tt, oss);

    return oss.str();
}

} // namespace alice

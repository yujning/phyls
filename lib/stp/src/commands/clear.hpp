#ifndef CLEAR_HPP
#define CLEAR_HPP

#include <iostream>
#include <cstdlib>
#include <alice/alice.hpp>

namespace alice
{
    // 定义一个命令类：clear_command
    class clear_command : public command
    {
    public:
        explicit clear_command(const environment::ptr &env)
            : command(env, "Clear the terminal screen")
        {
            // 可添加选项，比如确认开关；这里保持最简
        }

    protected:
        void execute() override
        {
#ifdef _WIN32
            std::system("cls");   // Windows 下清屏
#else
            std::system("clear"); // Linux / macOS 下清屏
#endif
        }
    };

    // 注册命令，命令名就是 stp 命令行里输入的 "clear"
    ALICE_ADD_COMMAND(clear, "STP")
}

#endif // CLEAR_HPP

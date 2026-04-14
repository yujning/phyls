#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include "../../src/include/algorithms/stp_dsd.hpp"

TEST_CASE("run_dsd_recursive function tests", "[dsd]") {
    SECTION("Valid input with length 2^1") {
        std::string input = "01";
        REQUIRE(run_dsd_recursive(input) == true);
    }

    SECTION("Valid input with length 2^2") {
        std::string input = "0110";
        REQUIRE(run_dsd_recursive(input) == true);
    }

    SECTION("Valid input with length 2^3") {
        std::string input = "01101001";
        REQUIRE(run_dsd_recursive(input) == true);
    }

    SECTION("Invalid input - not power of two length") {
        std::string input = "011"; // length 3
        REQUIRE(run_dsd_recursive(input) == false);
    }

    SECTION("Invalid input - empty string") {
        std::string input = "";
        REQUIRE(run_dsd_recursive(input) == false);
    }

    SECTION("Invalid input - contains non-01 characters") {
        std::string input = "0120";
        REQUIRE(run_dsd_recursive(input) == false);
    }

    SECTION("Edge case - minimum valid input") {
        std::string input = "0";
        REQUIRE(run_dsd_recursive(input) == true);
    }

    SECTION("Edge case - maximum practical input") {
        std::string input(1024, '0'); // 2^10 length
        REQUIRE(run_dsd_recursive(input) == true);
    }
}
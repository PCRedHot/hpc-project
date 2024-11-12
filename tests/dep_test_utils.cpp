#define CATCH_CONFIG_MAIN

#include <cstdio>
#include <string>
#include <cmath>
#include <functional>

#include "catch2/catch_all.hpp"
#include "exprtk.hpp"
#include "utils.hpp"


TEST_CASE("Utils Exprtk Test 1", "[utils_exprtk_1]") {
    auto func = fin_diff::expr_func_2d_from_string<double, double>("3 / 2 * sin(2 * x) * sinh(y)");

    double x, y;

    for (x = double(0); x <= double(1); x += double(0.001)) {
        for (y = double(0); y <= double(1); y += double(0.001)) {
            const double z = func(x, y);

            REQUIRE(z == 3.0 / 2.0 * std::sin(2 * x) * std::sinh(y));
        }
    }

}
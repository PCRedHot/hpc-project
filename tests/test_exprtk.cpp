#define CATCH_CONFIG_MAIN

#include <cstdio>
#include <string>
#include <cmath>

#include "catch2/catch_all.hpp"
#include "exprtk.hpp"

TEST_CASE("Exprtk Test 1", "[exprtk_1]") {
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    const std::string expression_string =
        "clamp(-1.0, sin(2 * pi * x) + cos(x / 2 * pi), +1.0)";

    double x;

    symbol_table_t symbol_table;
    symbol_table.add_variable("x", x);
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    parser.compile(expression_string, expression);

    for (x = double(-5); x <= double(+5); x += double(0.001)) {
        const double y = expression.value();

        CHECK(y == std::clamp(std::sin(2 * M_PI * x) + std::cos(x / 2 * M_PI), -1.0, 1.0));
    }
}

TEST_CASE("Exprtk Test 2", "[exprtk_2]") {
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    const std::string expression_string =
        "8 * pi^2 * (sin(2 * pi * x) * sin(2 * pi * y))";

    double x, y;

    symbol_table_t symbol_table;
    symbol_table.add_variable("x", x);
    symbol_table.add_variable("y", y);
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    parser.compile(expression_string, expression);

    for (x = double(0); x <= double(1); x += double(0.001)) {
        for (y = double(0); y <= double(1); y += double(0.001)) {
            const double z = expression.value();

            CHECK(z == 8 * M_PI * M_PI * (std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y)));
        }
    }
}

TEST_CASE("Exprtk Test 3", "[exprtk_3]") {
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    const std::string expression_string =
        "3 / 2 * sin(2 * x) * sinh(y)";

    double x, y;

    symbol_table_t symbol_table;
    symbol_table.add_variable("x", x);
    symbol_table.add_variable("y", y);
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    parser.compile(expression_string, expression);

    for (x = double(0); x <= double(1); x += double(0.001)) {
        for (y = double(0); y <= double(1); y += double(0.001)) {
            const double z = expression.value();

            CHECK(z == 3.0 / 2.0 * std::sin(2 * x) * std::sinh(y));
        }
    }
}
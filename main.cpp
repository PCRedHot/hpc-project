#include <iostream>

int main(int argc, char* argv[]) {
#ifdef __DEBUG__
    std::cout << "Debug mode is enabled" << std::endl;
#endif
    return 0;
}
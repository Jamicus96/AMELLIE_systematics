#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>

int main() {
    std::string init_str = "Long, ass, test, thinggy.";
    std::cout << init_str << std::endl;

    std::vector<std::string> str_lst;
    std::size_t pos;
    bool splitting = true;
    while (splitting) {
        pos = init_str.find(", ", 0);
        if ((unsigned int)pos < init_str.size()) {
            str_lst.push_back(init_str.substr(0, pos));
            init_str = init_str.substr(pos + 2);
        } else {
            splitting = false;
        }
    }
    str_lst.push_back(init_str);

    for (unsigned int i = 0; i < str_lst.size(); ++i) {
        std::cout << str_lst.at(i) << ", ";
    }
    std::cout << std::endl;

    return 0;
}
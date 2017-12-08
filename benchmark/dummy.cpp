#include <include/methodB.hpp>
#include <include/methodD.hpp>
#include <include/methodH.hpp>
#include <include/methodR.hpp>
#include <include/methodSH.hpp>
#include <include/methodP.hpp>

#include <iostream>
#include <random>

int main() {
    std::vector<sampling::ULONG> sample;
    sample.reserve(100000);
    sampling::HashSampling<> hs(std::random_device{}(), 1024);
    sampling::SeqDivideSampling<> s(hs, 1024, std::random_device{}());
    s.sample((1ULL << 31), 1000000, [&sample](const auto& x) { sample.push_back(x); });
    std::cout << sample.size() << std::endl;
    return 0;
}

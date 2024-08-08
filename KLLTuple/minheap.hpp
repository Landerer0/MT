#ifndef MINHEAP_HPP
#define MINHEAP_HPP

#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

class MinHeap {
public:
    MinHeap(int n_elements, int n_levels);
    MinHeap();
    void insert(pair<int64_t, int64_t> element);
    uint64_t sizeInBytes();
    void printMinHeap();
    vector<pair<int64_t, int64_t>> getHeap();
    int getLevels();
    int getCapacity();
    vector<pair<int64_t, int64_t>> getSortHeap();

private:
    vector<pair<int64_t, int64_t>> heap;
    int levels;

    bool incrementIfExists(pair<int64_t, int64_t> element); // retorna true cuando se incremento
    void heapify(int i);
};

#endif  // MINHEAP_HPP

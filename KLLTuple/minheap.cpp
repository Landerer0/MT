#include "minheap.hpp"

using namespace std;

MinHeap::MinHeap(int n_elements, int n_levels){
    if (n_elements > 0) {
        while ((1 << levels) < n_elements)
            levels++;
    } else if (n_levels > 0) {
        levels = n_levels;
    }
    heap.reserve((1 << levels) - 1);
}

void MinHeap::insert(pair<int64_t, int64_t> element) {
    //cout << "Insertar elemento (" << element.first << "," << element.second << ") " << heap.size() << " " << (1 << levels) - 1 << endl;

    if (incrementIfExists(element)) {
        // cout << "Elemento con el mismo pair.second encontrado y aumentado." << endl;
        // printMinHeap();
        return;
    }

    if (heap.size() == (1 << levels) - 1) {
        if (element.first > heap[0].first) {
            heap[0] = element;
            heapify(0);
        }

    } else {
        heap.push_back(element);
        int i = heap.size() - 1;
        while (i > 0) {
            int parent = (i - 1) / 2;
            if (heap[parent].first > heap[i].first) {
                swap(heap[i], heap[parent]);
                i = parent;
            } else {
                break;
            }
        }
    }
    //printMinHeap();
}

bool MinHeap::incrementIfExists(pair<int64_t, int64_t> element) {
    auto it = find_if(heap.begin(), heap.end(),
                           [&element](const pair<int64_t, int64_t>& pair) {
                               return pair.second == element.second;
                           });

    if (it != heap.end()) {
        it->first += element.first;
        heapify(distance(heap.begin(), it));
        return true;
    }

    return false;
}

void MinHeap::printMinHeap() {
    if (heap.empty()) {
        cout << "El Min Heap está vacío." << endl;
        return;
    }
    cout << "minHeap tiene a lo más " << heap.size() << " elementos" << endl;

    int current_level = 0;
    int level_size = 1;

    for (const auto& element : heap) {
        if (current_level == 0) {
            cout << "Nivel " << current_level << ": ";
        }

        cout << "First: " << element.first << ", Second: " << element.second << "   ";

        if (--level_size == 0) {
            cout << "\n";
            current_level++;
            level_size = 1 << current_level;
            if (level_size > heap.size() - 1) {
                break;
            }
            cout << "Nivel " << current_level << ": ";
        }
    }
}

uint64_t MinHeap::sizeInBytes(){
    uint64_t size = heap.size() * sizeof(heap[0]);
    size+= sizeof(levels);
    return size;
}

vector<pair<int64_t,int64_t>> MinHeap::getSortHeap() {
    if (heap.empty()) {
        cout << "El Min Heap está vacío." << endl;
    }

    vector<pair<int64_t,int64_t>> toReturn = heap;

    std::sort(toReturn.begin(), toReturn.end(), [](const std::pair<int64_t, int64_t>& a, const std::pair<int64_t, int64_t>& b) {
        return a.first > b.first;
    });

    return toReturn;
}


int MinHeap::getLevels(){
    return levels;
}


int MinHeap::getCapacity(){
    return heap.size();
}

vector<pair<int64_t, int64_t>> MinHeap::getHeap()  {
    return heap;
}

void MinHeap::heapify(int i) {
    int smallest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < heap.size() && heap[left].first < heap[smallest].first) {
        smallest = left;
    }

    if (right < heap.size() && heap[right].first < heap[smallest].first) {
        smallest = right;
    }

    if (smallest != i) {
        std::swap(heap[i], heap[smallest]);
        heapify(smallest);
    }
}

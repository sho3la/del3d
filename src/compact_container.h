// del3d - the slot allocator that fixes element order.
//
// Cells and vertices are handed out by this container, and the order it hands
// them out in is part of del3d's output contract: callers that break ties by
// "whichever incident cell comes first" see that order directly. During an
// insertion the star-hole step deletes the conflict cells and creates new ones,
// and freed slots are reused, so the sequence cells come out in is a function
// of the container's reuse policy - not of the geometry.
//
// The policy is therefore specified rather than left to an allocator:
//
//   * blocks grow by a fixed schedule - the first holds 14 elements and each
//     subsequent block 16 more than the last (14, 30, 46, 62, ...);
//   * a new block's slots are pushed onto the free list in REVERSE index
//     order, so that popping them yields ascending memory order;
//   * destroy() pushes the freed slot onto the head of the free list, making
//     reuse last-in-first-out: the slot freed most recently is the next one
//     handed out;
//   * iteration is memory order - blocks in allocation order, slots in index
//     order - skipping free slots.
//
// Here a "pointer" is a stable index into one flat array. Appending a block
// appends its slots to that array, so index order *is* memory order and
// iteration is a linear scan with a liveness test. Indices stay valid across
// block growth because nothing is ever moved or shrunk.
#ifndef DEL3D_COMPACT_CONTAINER_H
#define DEL3D_COMPACT_CONTAINER_H

#include <cstddef>
#include <vector>

namespace del3d {
namespace detail {

template <class T>
class CompactContainer {
public:
    using Index = int;
    static constexpr Index kNull = -1;

    CompactContainer() = default;

    /// Allocates one element: pop the free-list head, growing by a fresh block
    /// first if the list is empty. The element is value-initialised.
    Index create() {
        if (free_list_ == kNull) allocate_new_block();
        const Index ret = free_list_;
        free_list_ = next_free_[ret];
        used_[ret] = true;
        items_[ret] = T();
        ++size_;
        return ret;
    }

    /// Releases one element: mark the slot free and push it onto the free-list
    /// head (LIFO). The slot keeps its index; nothing is moved.
    void destroy(Index i) {
        used_[i] = false;
        next_free_[i] = free_list_;
        free_list_ = i;
        --size_;
    }

    bool is_used(Index i) const { return used_[i]; }
    /// Live elements.
    std::size_t size() const { return size_; }
    /// Allocated slots, live and free - i.e. the valid index range.
    std::size_t capacity() const { return items_.size(); }

    T&       operator[](Index i)       { return items_[i]; }
    const T& operator[](Index i) const { return items_[i]; }

    /// Calls f(index) for every live element in memory order.
    template <class F>
    void for_each(F&& f) const {
        for (Index i = 0; i < Index(items_.size()); ++i)
            if (used_[i]) f(i);
    }

    void clear() {
        items_.clear(); used_.clear(); next_free_.clear();
        free_list_ = kNull; size_ = 0; block_size_ = kFirstBlockSize;
    }

private:
    static constexpr std::size_t kFirstBlockSize = 14;
    static constexpr std::size_t kBlockIncrement = 16;

    /// Appends a block of block_size_ slots and threads them onto the free
    /// list back to front, so that the subsequent pops walk the block front to
    /// back and allocation order matches iteration order.
    void allocate_new_block() {
        const Index base = Index(items_.size());
        items_.resize(items_.size() + block_size_);
        used_.resize(items_.size(), false);
        next_free_.resize(items_.size(), kNull);
        for (std::size_t k = block_size_; k >= 1; --k) {
            const Index slot = base + Index(k - 1);
            next_free_[slot] = free_list_;
            free_list_ = slot;
        }
        block_size_ += kBlockIncrement;
    }

    std::vector<T>     items_;
    std::vector<char>  used_;      ///< liveness flag per slot
    std::vector<Index> next_free_; ///< free-list link, valid for free slots
    Index              free_list_ = kNull;
    std::size_t        size_ = 0;
    std::size_t        block_size_ = kFirstBlockSize;
};

} // namespace detail
} // namespace del3d

#endif // DEL3D_COMPACT_CONTAINER_H

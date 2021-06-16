// MIT License
//
// Copyright (c) 2021 Oliver Schick
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <mpo19/cluster.hpp>
#include <cassert>
#include <cstring>

cluster::item_t& cluster::item_t::get_referenced_item()
{
    using boost::get;
    item_t** referenced_item = get<item_t*>(std::addressof(value_));
    //If value_ does not contain an item_t*, then it refers to self
    if(nullptr == referenced_item) return *this;

    for(;;) {
        item_t** check =
            get<item_t*>(std::addressof((**referenced_item).value_));
        if(nullptr == check) break;
        referenced_item = check;
    }
    return **referenced_item;
}

cluster::item_t& cluster::union_t::get_left()
{
    assert(nullptr != left);
    return left->get_referenced_item();
}

cluster::item_t& cluster::union_t::get_right()
{
    assert(nullptr != right);
    return right->get_referenced_item();
}

cluster& merge(cluster& left, cluster& right, size_t rank)
{
    assert(std::addressof(left) != std::addressof(right)
           && "Cannot merge the same cluster with itself");
           
    //Right is invalid, so there is nothing to merge
    //If left is invalid the returned cluster will be invalid too
    if(right.item_ == nullptr)
        return left;
    
    //Only right  is valid
    if(left.item_ == nullptr){
        assert(right.item_ != nullptr);
        left = std::move(right);
        return left;
    }

    cluster::item_t& left_item = left.item_->get_referenced_item();
    cluster::item_t& right_item = right.item_->get_referenced_item();

    assert(std::addressof(left_item) != std::addressof(right_item)
           && "left and right should not point to same cluster through indirection");

    assert(nullptr == boost::get<cluster::item_t*>(std::addressof(left_item.value_))
           && "left_item should not be a pointer to another item");

    assert(nullptr == boost::get<cluster::item_t*>(std::addressof(right_item.value_))
           && "right_item should not be a pointer to another item");

    left_item.value_ = cluster::union_t{
        std::make_unique<cluster::item_t>(std::move(left_item.value_))
        , std::make_unique<cluster::item_t>(std::move(right_item.value_))
        , rank
    };

    //Update indirections to point to left_item
    right_item.value_ = std::addressof(left_item);
    if(std::addressof(right_item) != right.item_.get())
        right.item_->value_ = std::addressof(left_item);

    return left;
}

std::vector<index_t> get_indices(cluster const& c)
{
    std::vector<cluster::leaf_t> res;
    struct visitor : boost::static_visitor<void> {

        visitor(std::vector<cluster::leaf_t>& v) : indices{v} {}

        result_type operator()(cluster::union_t& u)
        {
            //The rank gives a good hint for the amount inserted into
            //the vector, as the vector size always end up being smaller
            // or equal to u.rank + 1
            indices.reserve(u.rank + 1);
            u.get_left().visit(*this);
            u.get_right().visit(*this);

            //Since ranks of children are strictly decreasing, u.rank + 1 > indices.capacity()
            //is only false for top-level function, if no space is reallocated in the recursive call
            //(which would be an error). If we are in the top-level function, we check if indices
            //contain indeed less or as much elements as u.rank + 1
            assert( (indices.capacity() > u.rank + 1 || indices.size() <= u.rank + 1)
                    && "indices should be smaller than u.rank + 1");
        }

        result_type operator()(cluster::leaf_t const& leaf)
        {
            indices.emplace_back(leaf);
        }

        std::vector<cluster::leaf_t>& indices;

    } v{res};

    c.visit(v);
    return res;
}

std::ostream& operator<<(std::ostream& os, cluster const& c)
{
    struct visitor : boost::static_visitor<void> {
        visitor(std::ostream& os) : os{os} {}
        result_type operator()(cluster::union_t& u)
        {
            os << '{';
            u.get_left().visit(*this);
            os << " U_" << u.rank << " ";
            u.get_right().visit(*this);
            os << '}';
        }

        result_type operator()(cluster::leaf_t const& l)
        {
            os << "{" << l << "}";
        }

        std::ostream& os;
    } v{os};
    c.visit(v);
    return os;
}

std::istream& operator>>(std::istream& is, cluster& c)
{
    char ch = is.get();
    //Use to prevent compiler warning
    (void) ch;
    assert(ch == '{');
    
    if(is.peek() == '{'){
        constexpr size_t union_sep_size = 3;
        cluster left_cluster, right_cluster;
        //reserve one additional character for null termination
        char buf[union_sep_size + 1];
        size_t rank;
        
        is >> left_cluster;
        is.read(buf, union_sep_size);
        buf[union_sep_size] = '\0';
        assert(strcmp(buf, " U_") == 0);
        
        is >> rank;
        
        is >> std::ws >> right_cluster;
        
        c = std::move(merge(left_cluster, right_cluster, rank));
        
        is >> ch;
        assert(ch == '}');
    }
    else{
        cluster::leaf_t leaf;
        is >> leaf;
        c = cluster{leaf};
        
        ch = is.get();
        (void) ch;
        assert(ch == '}');
    }
    return is;
}

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

#ifndef MPO19_CLUSTER_13072020
#define MPO19_CLUSTER_13072020

#include <mpo19/typedefs.hpp>

#include <vector>
#include <type_traits>
#include <iosfwd>

#include <boost/variant.hpp>

struct cluster {

    struct item_t;

    //Base case, when no clusters have been merged
    using leaf_t = index_t;

    struct union_t {
        std::unique_ptr<item_t> left;
        std::unique_ptr<item_t> right;
        size_t rank;

        item_t& get_left();

        item_t& get_right();
    };

    struct item_t {

        using value_type = boost::variant<leaf_t, union_t, item_t*>;

        item_t(leaf_t leaf)
            : value_{leaf} {}

        item_t(union_t&& u)
            : value_{std::move(u)} {}

        item_t(item_t* i)
            : value_{i} {}

        item_t(value_type&& v)
            : value_{std::move(v)} {}

        item_t& get_referenced_item();
        value_type& get_value()
        {
            return get_referenced_item().value_;
        }

        //Visit the item, with indirections being automatically resolved.
        template<typename Visitor>
        auto visit(Visitor&& v)
        {
            visitor_<Visitor> visitor{std::forward<Visitor>(v)};
            return boost::apply_visitor(visitor, get_value());
        }

        friend cluster& merge(cluster& left, cluster& right, size_t rank);

        friend struct cluster;

    private:

        template<typename Visitor>
        struct visitor_ {
            using result_type = typename std::decay_t<Visitor>::result_type;
            Visitor&& v;

            template<typename T>
            decltype(auto) operator()(T&& t)
            {
                return v(std::forward<T>(t));
            }

            result_type operator()(item_t* i){
                return boost::apply_visitor(*this, i->get_value());
            }
        };

        value_type value_;
    };

    cluster() = default;

    cluster(cluster&& rhs) = default;

    cluster(cluster const&) = delete;

    cluster(leaf_t const& leaf)
        : item_{std::make_unique<item_t>(leaf)} {}

    cluster& operator=(cluster&&) = default;

    cluster& operator=(cluster const&) = delete;

    ~cluster() = default;

    bool is_cluster() const
    {
        return nullptr != item_&& nullptr == boost::get<item_t*>(std::addressof(item_->value_));
    }

    friend cluster& merge(cluster& left, cluster& right, size_t rank);

    template<typename Visitor>
    auto visit(Visitor&& v) const
    {
        assert(nullptr != item_);
        return item_->visit(std::forward<Visitor>(v));
    }

private:

    //We use a unique_ptr, so that the address of the underlying
    //item_t remains the same after moving the cluster.
    //This is important, so that the indirections pointing to an item
    //remain valid.
    std::unique_ptr<item_t> item_ = nullptr;
};

std::vector<cluster::leaf_t> get_indices(cluster const& c);

std::ostream& operator<<(std::ostream& os, cluster const& c);
std::istream& operator>>(std::istream& is, cluster& c);

#endif //MPO19_CLUSTER_13072020

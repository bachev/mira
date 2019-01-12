/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2011 and later by Bastien Chevreux
 *
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 */


#ifndef _util_stlimprove_h
#define _util_stlimprove_h

#include "boost/concept/assert.hpp"
#include <numeric>      // std::accumulate
#include <algorithm>

#if __has_include(<parallel/algorithm>)   // for parallel sort
#  include <parallel/algorithm>
#  define HAVEGNUPARALLEL
#endif


/*************************************************************************
 *
 * MSTD namespace mainly for saving typing and making things easier
 *  to understand by creating version of the standard library functions
 *  which work on complete containers and do not need a first and last
 *  iterator.
 * Simple example:
 *   mstd::fill
 *
 *************************************************************************/

namespace mstd {
  template <class Container>
  inline auto accumulate(const Container & c) {
    return std::accumulate(begin(c),end(c),0);
  }

  template <class Container, class T>
  inline auto binary_search(Container & c, const T& val) {
    return std::binary_search(begin(c),end(c),val);
  }
  template <class Container, class T, class Compare>
  inline auto binary_search(Container & c, const T& val, Compare comp) {
    return std::binary_search(begin(c),end(c),val,comp);
  }
  template <class Container, class T>
  inline auto binary_search(const Container & c, const T& val) {
    return std::binary_search(begin(c),end(c),val);
  }
  template <class Container, class T, class Compare>
  inline auto binary_search(const Container & c, const T& val, Compare comp) {
    return std::binary_search(begin(c),end(c),val,comp);
  }

  template <class Container, class T>
  inline bool contains(Container & c, const T& val) {
    return (std::find(begin(c),end(c),val) != end(c));
  }

  template <class Container, class OutputIterator>
  inline auto copy(Container & c, OutputIterator result) {
    return std::copy(begin(c),end(c),result);
  }

  template <class Container, class OutputIterator>
  inline auto copy(const Container & c, OutputIterator result) {
    return std::copy(begin(c),end(c),result);
  }

  template <class Container, class T>
  void fill(Container & c, const T& val) {
    std::fill(begin(c),end(c),val);
  }

  template <class Container, class T>
  inline auto lower_bound(Container & c, const T& val) {
    return std::lower_bound(begin(c),end(c),val);
  }
  template <class Container, class T, class Compare>
  inline auto lower_bound(Container & c, const T& val, Compare comp) {
    return std::lower_bound(begin(c),end(c),val,comp);
  }
  template <class Container, class T>
  inline auto lower_bound(const Container & c, const T& val) {
    return std::lower_bound(begin(c),end(c),val);
  }
  template <class Container, class T, class Compare>
  inline auto lower_bound(const Container & c, const T& val, Compare comp) {
    return std::lower_bound(begin(c),end(c),val,comp);
  }

  template <class Container>
  inline auto & max_element(Container & c){ return std::max_element(begin(c),end(c));}
  template <class Container>
  inline auto & max_element(const Container & c){ return std::max_element(begin(c),end(c));}
  template <class Container>
  inline auto & min_element(Container & c){ return std::min_element(begin(c),end(c));}
  template <class Container>
  inline auto & min_element(const Container & c){ return std::min_element(begin(c),end(c));}

  template <class Container>
  void reverse(Container & c) {
    // this does not work for HDeque. Why?
    // std::reverse(begin(c),end(c));
    std::reverse(std::begin(c),std::end(c));
  }

  // standard sort (non-parallel)
  template <class Container>
  void ssort(Container & c) {
    std::sort(begin(c),end(c));
  }
  template <class Container, class Compare>
  void ssort(Container & c, Compare comp) {
    std::sort(begin(c),end(c),comp);
  }

#ifdef HAVEGNUPARALLEL
  template <class Container>
  void psort(Container & c) {
    __gnu_parallel::sort(begin(c),end(c));
  }
  template <class Container, class Compare>
  void psort(Container & c, Compare comp) {
    __gnu_parallel::sort(begin(c),end(c),comp);
  }
#else
  // no GNU parallel, use standard sort
  template <class Container>
  void psort(Container & c) {
    std::sort(begin(c),end(c));
  }
  template <class Container, class Compare>
  void psort(Container & c, Compare comp) {
    std::sort(begin(c),end(c),comp);
  }
#endif

  template <class Container>
  inline auto unique(Container & c) {
    return std::unique(begin(c),end(c));
  }

  template <class Container, class T>
  inline auto upper_bound(Container & c, const T& val) {
    return std::upper_bound(begin(c),end(c),val);
  }
  template <class Container, class T, class Compare>
  inline auto upper_bound(Container & c, const T& val, Compare comp) {
    return std::upper_bound(begin(c),end(c),val,comp);
  }
  template <class Container, class T>
  inline auto upper_bound(const Container & c, const T& val) {
    return std::upper_bound(begin(c),end(c),val);
  }
  template <class Container, class T, class Compare>
  inline auto upper_bound(const Container & c, const T& val, Compare comp) {
    return std::upper_bound(begin(c),end(c),val,comp);
  }
}

#endif

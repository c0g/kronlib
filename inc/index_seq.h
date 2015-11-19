#ifndef REDI_INDEX_SEQ_H
#define REDI_INDEX_SEQ_H

// Copyright Jonathan Wakely 2012-2015
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#if __cplusplus < 201402L
# include <cstddef>
#else
# include <utility>
#endif

namespace redi
{
#if __cplusplus < 201402L

  /// A type that represents a parameter pack of zero or more integers.
  template<std::size_t... Indices>
    struct index_sequence
    {
      /// Generate an index_sequence with an additional element.
      template<std::size_t N>
        using append = index_sequence<Indices..., N>;
    };

  /// Unary metafunction that generates an index_sequence containing [0, Size)
  template<std::size_t Size>
    struct make_index_sequence
    {
      using type
        = typename make_index_sequence<Size-1>::type::template append<Size-1>;
    };

  // Terminal case of the recursive metafunction.
  template<>
    struct make_index_sequence<0u>
    {
      typedef index_sequence<> type;
    };

  template<typename... Types>
    using index_sequence_for
      = typename make_index_sequence<sizeof...(Types)>::type;

#else
  using std::index_sequence;
  using std::make_index_sequence;
  using std::index_sequence_for;
#endif

}  // namespace redi

#endif  // REDI_INDEX_SEQ_H

// vi: set ft=cpp sw=2:

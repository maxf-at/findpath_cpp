


#include <common.hpp>

#include <postprocessing.hpp>

// #include <single_findpath.hpp>
// #include <single_findpath_test.hpp>
// #include <single_findpath_pool.hpp>

// #include <single_findpath_s.hpp>

// #include <single_findpath_s_old.hpp>
#include <single_findpath_h.hpp>
// #include <old/single_findpath_h.hpp>

// #include <single_findpath_hash.hpp>

// #include <single_findpath_s4hash.hpp>
#include <single_findpath_i.hpp>


#include <multi_findpath.hpp>


#include <merge_findpath.hpp>
// #include <merge_findpath_hash.hpp>
// #include <merge_findpath_benchmark.hpp>


#include <mfe_preprocessing.hpp>

// // emplace back without copy: construct it inplace, how emplace is meant
// // to be used
// auto& s = v.emplace_back("Hello");
// auto &s = v.back(); // not needed
// s += ", ";
// s += "world";

// // push back without copy
// string s = "Hello";
// s +"World";
// s + get_lots_of_stuff();
// v.push_back(std::move(s));

// // Create a deque containing integers
// std::deque<int> d = {7, 5, 16, 8};
//  // Add an integer to the beginning and end of the deque
// d.push_front(13);
// d.push_back(25);

void testfunc(const char* seq, const char* s1, const char* s2, float search_width_multiplier);
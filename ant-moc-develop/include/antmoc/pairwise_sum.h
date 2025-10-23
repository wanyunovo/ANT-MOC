/**
 * @file pairwise_sum.h
 * @brief Utility function for the accurate pairwise sum of a list of floating
 *        point numbers.
 * @author William Boyd (wboyd@mit.edu)
 * @date June 13, 2013
 */

#include <cstddef>

namespace antmoc {

/**
 * @brief Performs a pairwise sum of an array of numbers.
 * @details This type of summation uses a divide-and-conquer algorithm which
 *          is necessary to bound the error for summations of large sequences
 *          of numbers.
 * @param vector an array of numbers
 * @return the sum of all numbers in the array
 */
template <typename T>
extern double pairwise_sum(T* vector, size_t length);

} /* namespace antmoc */

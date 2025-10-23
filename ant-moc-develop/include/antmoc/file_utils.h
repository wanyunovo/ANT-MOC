/// \file include/file_utils.h
/// \brief Utilities for file manipulation (POSIX)
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef FILE_UTILES_H_
#define FILE_UTILES_H_

#include <string>

namespace antmoc {

namespace fileutils {


/// \brief Get the global output directory
std::string getOutputDirectory();

/// \brief Set the global output directory
void setOutputDirectory(const std::string &dir);

/// \brief Check the existence of a regular file
bool existsFile(const std::string &path);

/// \brief Check the existence of a directory
bool existsDirectory(const std::string &path);

/// \brief Create a directory recursively
/// \param path Path to the directory
void createDirectory(const std::string &path);

/// \brief Get file extension from file name.
/// \param file File name.
/// \return File extension if succeed, "" otherwise.
std::string getExtension(const std::string &file);


} // namespace fileutils

} // namespace antmoc

#endif  // FILE_UTILES_H_

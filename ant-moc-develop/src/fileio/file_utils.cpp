#include "antmoc/file_utils.h"
#include "antmoc/log.h"

#include <sys/stat.h>
#include <errno.h>
#include <cstring>

#ifdef MAX_PATH_LEN
#undef MAX_PATH_LEN
#endif
#define MAX_PATH_LEN 256

namespace antmoc {

namespace fileutils {


/// Global output directory
static std::string _output_directory = "./";


std::string getOutputDirectory() {
  return _output_directory;
}


void setOutputDirectory(const std::string &dir) {
  _output_directory = dir;
}


bool existsFile(const std::string &path) {
  struct stat path_stat;
  return !stat(path.c_str(), &path_stat) && S_ISREG(path_stat.st_mode);
}


bool existsDirectory(const std::string &path) {
  struct stat path_stat;
  return !stat(path.c_str(), &path_stat) && S_ISDIR(path_stat.st_mode);
}


/// \details Thanks to https://gist.github.com/JonathonReinhart/8c0d90191c38af2dcadb102c4e202950
void createDirectory(const std::string &path) {

  // Check if the directory exists
  if (existsDirectory(path))
    return;

  // Check the string length
  if (path.size() >= MAX_PATH_LEN) {
    log::error("Failed to create directory '{}': exceeded maximum length {}", path, MAX_PATH_LEN);
  }

  errno = 0;

  char _path[MAX_PATH_LEN];
  char *p;

  // Copy string to make it mutable
  std::strcpy(_path, path.c_str());

  // Iterate the string
  for (p = _path + 1; *p; p++) {
	if (*p == '/') {
	  // Temporarily truncate the string to get a directory path
	  *p = '\0';

	  if (mkdir(_path, S_IRWXU) != 0) {
		if (errno != EEXIST)
          log::error("Failed to create directory '{}': {}", path, errno);
	  }

	  *p = '/';
	}
  }

  if (mkdir(_path, S_IRWXU) != 0) {
	if (errno != EEXIST)
      log::error("Failed to create directory '{}': {}", path, errno);
  }

  return;
}


std::string getExtension(const std::string &file) {
    // Find the last occurrence of '.'
    std::size_t pos = file.rfind('.');

    if (pos != std::string::npos) {
        return file.substr(pos);
    }

    return "";
}


} // namespace fileutils

} // namespace antmoc

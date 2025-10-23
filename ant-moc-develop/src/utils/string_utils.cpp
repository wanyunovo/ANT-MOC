#include "antmoc/string_utils.h"
#include "antmoc/log.h"

#include <algorithm>
#include <cctype>
#include <sstream>

namespace antmoc
{


namespace stringutils {


  //--------------------------------------------------------------------
  // Case conversion
  //--------------------------------------------------------------------
  /// \details Use transform algorithm to process characters one by one
  std::string &toUpper(std::string &s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::toupper(c); }
                  );
    return s;
  }

  /// \details Use transform algorithm to process characters one by one
  std::string toUpper(const std::string &s) {
    std::string newstr(s);
    std::transform(newstr.begin(), newstr.end(), newstr.begin(),
                   [](unsigned char c) { return std::toupper(c); }
                  );
    return newstr;
  }

  /// \details Use transform algorithm to process characters one by one
  std::string &toLower(std::string &s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); }
                  );
    return s;
  }

  /// \details Use transform algorithm to process characters one by one
  std::string toLower(const std::string &s) {
    std::string newstr(s);
    std::transform(newstr.begin(), newstr.end(), newstr.begin(),
                   [](unsigned char c) { return std::tolower(c); }
                  );
    return newstr;
  }


  //--------------------------------------------------------------------
  // String spliting and concatenating
  //--------------------------------------------------------------------
  /// \details Split a string into words. The delimiter can be any string.
  /// \        During spliting, all occurences of the delimiter are skipped.
  StringVec splitString(std::string input, const std::string &delimiter)
  {
    size_t last_pos = 0;
    size_t next_pos = 0;
    StringVec result;

    // Find the delimiter in the string
    while ( (next_pos = input.find(delimiter, last_pos))
            != std::string::npos ) {
      if (next_pos != last_pos) {
        // Skip the leading and adjacent delimiter
        result.push_back(input.substr(last_pos, next_pos - last_pos));
      }
      // Skip the delimiter we found
      last_pos = next_pos + delimiter.size();
    }
    // Push the rest back to the vector
    if (last_pos < input.size()) {
      result.push_back(input.substr(last_pos));
    }
    return result;
  }


  /// \details An overloaded function for split strings into words.
  ///          This is only for space-delimited strings.
  StringVec splitString(std::string input)
  {
    StringVec words;
    std::stringstream ss(input);
    std::string s;

    while (ss >> s) {
      words.push_back(s);
    }

    return words;
  }


  //--------------------------------------------------------------------
  // Trimming methods, thanks to
  // https://stackoverflow.com/questions/216823/whats-the-best-way-to-
  // trim-stdstring
  //--------------------------------------------------------------------

  /// \details Erase spaces from the start of the string
  std::string& ltrim(std::string &s) {
    s.erase(s.begin(),
            std::find_if(s.begin(), s.end(),
                         [](int ch) { return !std::isspace(ch); }
                        )
           );
    return s;
  }

  /// \details Erase spaces from the start of the string
  std::string ltrim(const std::string &s) {
    std::string str(s);
    return ltrim(str);
  }

  /// \details Erase spaces from the end of the string
  std::string& rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         [](int ch) { return !std::isspace(ch); }
                        ).base(),
            s.end()
           );
    return s;
  }

  /// \details Erase spaces from the end of the string
  std::string rtrim(const std::string &s) {
    std::string str(s);
    return rtrim(str);
  }

  std::string& trim(std::string &s) {
    return rtrim(ltrim(s));
  }

  std::string trim(const std::string &s) {
    std::string str(s);
    return trim(str);
  }


  //--------------------------------------------------------------------
  // String checking
  //--------------------------------------------------------------------
  /// \details Iterate the string and check each character
  bool isSpaces(const std::string &s) {
    return std::all_of(s.cbegin(), s.cend(),
                       [](const int ch)
                          { return std::isspace(ch); }
                      );
  }


  //--------------------------------------------------------------------
  // Substitution
  //--------------------------------------------------------------------
  std::string& spaceToUnderscore(std::string &s) {
    trim(s);
    std::replace(s.begin(), s.end(), ' ', '_');
    return s;
  }


  std::string spaceToUnderscore(const std::string &s) {
    std::string str(s);
    return spaceToUnderscore(str);
  }


  std::string& underscoreToSpace(std::string &s) {
    trim(s);
    std::replace(s.begin(), s.end(), '_', ' ');
    return s;
  }


  std::string underscoreToSpace(const std::string &s) {
    std::string str(s);
    return underscoreToSpace(str);
  }


  /// \details This algorithm first splits the string with the given delimiter,
  ///          and then trys to convert substrings to integers which represent
  ///          start, stop and an optional step.
  ///          If start == stop, only 1 item will be returned.
  ///          If start != stop, the step must be a valid number so that we can
  ///          establish an arithmetic progression.
  //           The delimiter defaults to ':', which is just like Python slice.
  std::set<int> toIntegerSet(std::string s, std::string delimiter) {
    if (s.empty())
      return std::set<int>(); // empty set

    auto strvec = splitString(s, delimiter);
    if (strvec.empty() || strvec.size() > 3) {
      log::ferror("Ill-formed integer range: '%s', correct range is "
                        "represented as 'start[%sstop[%sstep]]'", s.c_str(),
                        delimiter.c_str(), delimiter.c_str());
    }

    std::set<int> int_set;

    try {

      int start = std::stoi(strvec[0]);
      int stop = start;
      int step = 1; // defaults to 1

      // start:stop
      if (strvec.size() > 1)
        stop = std::stoi(strvec[1]);

      // start:stop:step
      if (strvec.size() == 3)
        step = std::stoi(strvec[2]);

      int_set.insert(start);

      if (start != stop) {
        if ((stop - start) * step <= 0) {
          throw std::logic_error("Ill-formed integer range");
        } else {
          int num_items = (stop - start) / step;
          for (int i = 1; i <= num_items; ++i)
            int_set.insert(start + i * step);
        }
      }

    } catch (...) {
      log::ferror("Ill-formed integer range: '%s',delimited by '%s'",
                        s.c_str(), delimiter.c_str());
    }

    return int_set;
  }

} // namespace stringutils

} // namespace antmoc

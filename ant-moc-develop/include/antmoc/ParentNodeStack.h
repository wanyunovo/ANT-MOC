/// \file include/ParentNodeStack.cpp
/// \brief A class for reading XML and TOML options.
/// \date Oct 18, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef PARENTNODESTACK_H_
#define PARENTNODESTACK_H_

#include "antmoc/string_utils.h"
#include "antmoc/xml_utils.h"

#include <toml.hpp>

#include <set>
#include <stack>
#include <string>

namespace antmoc {

///---------------------------------------------------------------------
/// \class A simple stack hold and manipulate parent nodes for diffenrent
///        file formats.
///---------------------------------------------------------------------
struct ParentNodeStack {

  ParentNodeStack() = default;

  virtual ~ParentNodeStack() = default;

  /// \brief Detect the file format (XML or TOML) and initialize the nodes.
  virtual bool initialize(const std::string &file) = 0;

  /// \brief Query the stack top for a child node.
  /// \param node_name Name of the child node.
  /// \param ret The value of the child node as a string.
  /// \return True if the child was found, false otherwise.
  virtual bool queryTop(const std::string &node_name, std::string &ret) = 0;

  /// \brief Push a parent node into the stack
  /// \details If the stack top is nullptr or the node value is empty,
  ///          nullptr will be pushed onto the stack.
  /// \param node_name Name of the child node.
  virtual void pushNextParent(const std::string &node_name) = 0;

  /// \brief Return the stack size.
  virtual size_t size() = 0;

  /// \brief Remove the stack top.
  virtual void pop() = 0;

  /// \brief Retrieve options and save it into the buffer.
  /// \details The top element of the stack will be queried.
  /// \param buffer A buffer to save the results.
  /// \param boolean_options A list of boolean option names.
  /// \param name_in_file Name of the child node.
  /// \param name_in_cli Name of the associated command line option name.
  void saveOption(StringVec &buffer, const std::set<std::string> &boolean_options,
                  const std::string &name_in_file, const std::string &name_in_cli);

  /// \brief Swap the top element of the stack.
  /// \param node_name Name of the child node.
  void swapTop(const std::string &node_name);
};


///---------------------------------------------------------------------
/// \class A simple stack hold and manipulate parent nodes for XML.
///---------------------------------------------------------------------
struct ParentNodeStackXML : public ParentNodeStack {

  XMLElement *PLACE_HOLDER = nullptr;

  ///< A document.
  tinyxml2::XMLDocument *doc;

  ///< A stack for recording parent nodes.
  std::stack<XMLElement *> parents;

  ParentNodeStackXML();

  ~ParentNodeStackXML();

  bool initialize(const std::string &file) override;
  void pushNextParent(const std::string &node_name) override;
  bool queryTop(const std::string &node_name, std::string &ret) override;

  size_t size() override  { return parents.size(); }
  void pop() override     { parents.pop(); }

};


///---------------------------------------------------------------------
/// \class A simple stack hold and manipulate parent nodes for TOML.
///---------------------------------------------------------------------
struct ParentNodeStackTOML : public ParentNodeStack {

  toml::value PLACE_HOLDER = toml::value("");

  ///< A stack for recording parent nodes.
  std::stack<toml::value> parents;

  ParentNodeStackTOML() = default;

  bool initialize(const std::string &file) override;
  void pushNextParent(const std::string &node_name) override;
  bool queryTop(const std::string &node_name, std::string &ret) override;

  size_t size() override { return parents.size(); }
  void pop() override { parents.pop(); }

};


} // namespace antmoc

#endif  // PARENTNODESTACK_H_

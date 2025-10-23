#include "antmoc/ParentNodeStack.h"
#include "antmoc/log.h"
#include "antmoc/toml_utils.h"

#include <sstream>

namespace antmoc {


//----------------------------------------------------------------------
// ParentNodeStack
//----------------------------------------------------------------------
void ParentNodeStack::saveOption(StringVec &buffer,
                                 const std::set<std::string> &boolean_options,
                                 const std::string &name_in_file,
                                 const std::string &name_in_cli) {
  std::string node_value;
  bool success = queryTop(name_in_file, node_value);

  if (success) {

    // Add the option to the buffer
    std::stringstream arg;
    arg << "--" << name_in_cli;
    if (boolean_options.count(name_in_file)) {
      arg << "=" << node_value;
      buffer.push_back(arg.str());
    }
    else {
      buffer.push_back(arg.str());
      buffer.push_back(node_value);
    }

    // Debugging
    log::debug("Parent node stack: saw option in file: '{}' '{}'", name_in_file, node_value);
  }
}


/// \brief Swap the top element of the stack.
void ParentNodeStack::swapTop(const std::string &node_name) {
  if (size() == 0)
    log::error("Internal failure: swapTop failed");

  pop();
  pushNextParent(node_name);
}


//----------------------------------------------------------------------
// ParentNodeStackXML
//----------------------------------------------------------------------
ParentNodeStackXML::ParentNodeStackXML() {
  doc = new tinyxml2::XMLDocument();
}

ParentNodeStackXML::~ParentNodeStackXML() {
  delete doc;
}

/// \brief Initialize the node stack.
bool ParentNodeStackXML::initialize(const std::string &file) {

  loadFileXML(*this->doc, file, true);

  // Get the root element, return false if there is no root
  tinyxml2::XMLElement *root = xmlutils::queryFirstChild(this->doc, "settings");

  if (root) {
    this->parents.push(root);
    return true;
  }

  return false;
}


bool ParentNodeStackXML::queryTop(const std::string &node_name, std::string &ret) {

  bool success = false;

  if ( auto parent = this->parents.top() ) {
    if ( auto node_value = xmlutils::queryNodeString(parent, node_name) ) {
      ret = std::string(node_value);
      success = true;

      log::debug("Parent node stack({}): queried top for '{}'", size(), node_name);
    }
  }
  return success;
}

void ParentNodeStackXML::pushNextParent(const std::string &node_name) {
  // Skip empty node names, which indicate parent tag.
  if (node_name.empty()) {
    log::debug("Parent node stack({}): skipping parent '{}'", size(), node_name);
  }
  else {

    auto parent = this->parents.top();
    if (parent == PLACE_HOLDER) {
      log::debug("Parent node stack({}): pushing place-holder parent", size());
      this->parents.push(PLACE_HOLDER);
    }
    else {
      log::debug("Parent node stack({}): pushing parent '{}'", size(), node_name);
      this->parents.push(xmlutils::queryFirstChild(parent, node_name));
    }
  }
}


//----------------------------------------------------------------------
// ParentNodeStackTOML
//----------------------------------------------------------------------
/// \brief Initialize the node stack.
bool ParentNodeStackTOML::initialize(const std::string &file) {

  auto root = toml::parse(file);
  this->parents.push(root);

  return true;
}


bool ParentNodeStackTOML::queryTop(const std::string &node_name, std::string &ret) {
  bool success = false;

  const auto &parent = this->parents.top();
  if (parent.is_table()) {
    if (parent.contains(node_name)) {
      const auto &node = toml::find(parent, node_name);
      ret = tomlutils::toString(node);
      success = true;

      log::debug("Parent node stack({}): queried top for '{}'", size(), node_name);
    }
  }
  else {
    log::error("Internal failure: wrongly formatted TOML node at '{}'", node_name);
  }
  return success;
}


void ParentNodeStackTOML::pushNextParent(const std::string &node_name) {
  // Skip empty node names, which indicate parent tag.
  if (node_name.empty()) {
    log::debug("Parent node stack({}): skipping parent '{}'", size(), node_name);
  }
  else {

    const auto &parent = this->parents.top();

    if (parent == PLACE_HOLDER) {
      log::debug("Parent node stack({}): pushing place-holder parent", size());
      this->parents.push(PLACE_HOLDER);
    }
    else if (parent.contains(node_name)) {
      log::debug("Parent node stack: pushing parent '{}'", node_name);
      const auto &node = parent.at(node_name);
      this->parents.push(node);
    }
  }
}


} // namespace antmoc

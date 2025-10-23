#include "antmoc/xml_utils.h"
#include "antmoc/log.h"
#include "antmoc/file_utils.h"
#include "antmoc/string_utils.h"

#include <limits>

namespace antmoc {


inline namespace xmlutils {


  void loadFileXML(XMLDocument &doc, const std::string &file,
                   bool enforced) {
    if (!fileutils::existsFile(file)) {
      if (enforced) {
        log::error("File I/O: Cannot find XML file: '{}'", file);
      }
    }
    else {
      auto status = doc.LoadFile(file.c_str());
      if(enforced && status != XMLError::XML_SUCCESS) {
        log::error("File I/O: Failed to load XML file '{}':\n{}", file, doc.ErrorStr());
      }
    }
  }


  void saveFileXML(XMLDocument &doc, const std::string &file) {
    if (fileutils::existsFile(file)) {
      log::warn_once("File I/O: Overwriting file '{}'", file);
    }

    if (XMLError::XML_SUCCESS != doc.SaveFile(file.c_str())) {
      log::warn_once("File I/O: Failed to save XML file '{}:\n{}", file, doc.ErrorStr());
    }
  }


  /// \brief Extract the first child from a given element
  /// \param parent the parent element
  /// \param child name of the child to be extracted
  /// \param enforced determine if the entry is mandatory
  /// \param info extra info attached to error messages
  /// \return a pointer to child element
  template <typename T>
  XMLElement* queryFirstChild(T *parent,
                              const std::string &child,
                              bool enforced,
                              const std::string &info) {

    if (!parent)
      log::error("Could not query an empty XML object");

    XMLElement *e = parent->FirstChildElement(child.c_str());
    if (enforced && !e) {
      log::ferror("When parsing '%s' of '%s', sub-element may not exist. '%s'",
                  child, parent->Value(), info);
    }
    return e;
  }

  // Explicit instantiation
  template XMLElement* queryFirstChild(XMLDocument*, const std::string&,
                                       bool, const std::string&);
  template XMLElement* queryFirstChild(XMLElement*, const std::string&,
                                       bool, const std::string&);


  /// \brief Extract a string attribute from a given element
  /// \param elem_xml the parent element
  /// \param attr_name name of the attribute to be extracted
  /// \param elem_name name of the parent element (for debugging)
  /// \param enforced determine if the entry is mandatory
  /// \return queried value if succeed, or nullptr if failed
  const char* queryStringAttr(XMLElement *elem_xml,
                              const std::string &attr_name,
                              const std::string &elem_name,
                              bool enforced) {
    auto attr = elem_xml->Attribute(attr_name.c_str());
    if (!attr && enforced) {
      log::ferror("An error occured when parsing '%s' of '%s': "
                  "attribute may not exist", attr_name, elem_name);
    }
    return attr;
  }


  /// \brief Extract a floating point attribute from a given element
  /// \param elem_xml the parent element
  /// \param attr_name name of the attribute to be extracted
  /// \param elem_name name of the parent element (for debugging)
  /// \param enforced determine if the entry is mandatory
  /// \return queried value if succeed, or infinity if failed
  double queryDoubleAttr(XMLElement *elem_xml,
                         const std::string &attr_name,
                         const std::string &elem_name,
                         bool enforced) {
    double def = std::numeric_limits<double>::infinity();
    double attr = elem_xml->DoubleAttribute(attr_name.c_str(), def);
    if ((attr == def) && enforced) {
      log::ferror("An error occured when parsing '%s' of '%s': "
                  "attribute may not exist", attr_name, elem_name);
    }
    return attr;
  }


  /// \brief Extract a floating point attribute from a given element
  /// \param elem_xml the parent element
  /// \param attr_name name of the attribute to be extracted
  /// \param elem_name name of the parent element (for debugging)
  /// \param enforced determine if the entry is mandatory
  /// \return queried value if succeed, or min() of int if failed
  int queryIntAttr(XMLElement *elem_xml,
                   const std::string &attr_name,
                   const std::string &elem_name,
                   bool enforced) {
    int def = std::numeric_limits<int>::min();
    int attr = elem_xml->IntAttribute(attr_name.c_str(), def);
    if ((attr == def) && enforced) {
      log::ferror("An error occured when parsing '%s' of '%s': "
                  "attribute may not exist", attr_name, elem_name);
    }
    return attr;
  }


  /// \details A node is either an attribute or an sub-element. If two nodes
  ///          with the same name exist, an error occurred.
  bool existNode(const XMLElement *parent, const std::string &child) {
    auto attr = parent->Attribute(child.c_str());
    auto elem = parent->FirstChildElement(child.c_str());

    if ( elem &&
         ( attr ||
           elem->NextSiblingElement(child.c_str())
         )
       ) {
      log::ferror("When parsing '%s', detected duplicated attribute "
                  "or sub-element '%s'", parent->Name(), child);
    }

    return (attr || elem);
  }


  /// \details If a node is enforced, it has to exist.
  const char* queryNodeString(XMLElement *parent,
                              const std::string &child,
                              bool enforced,
                              const std::string &info) {

    if (!existNode(parent, child)) {
      if (enforced) {
        log::ferror("When parsing '%s' of '%s', attribute/sub-element "
                    "may not exist. '%s'", child, parent->Name(), info);
      }
      else {
        return defaultQueriedString();
      }
    }

    auto attr = parent->Attribute(child.c_str());
    auto elem = parent->FirstChildElement(child.c_str());

    return (attr) ? attr : elem->GetText();
  }


  /// \details If a node is enforced, it has to exist.
  double queryNodeDouble(XMLElement *parent,
                         const std::string &child,
                         bool enforced,
                         const std::string &info) {

    if (!existNode(parent, child)) {
      if (enforced) {
        log::ferror("When parsing '%s' of '%s', attribute/sub-element "
                    "may not exist. '%s'", child, parent->Name(), info);
      }
      else {
        return defaultQueriedDouble();
      }
    }

    auto attr = parent->Attribute(child.c_str());
    auto elem = parent->FirstChildElement(child.c_str());
    double result = defaultQueriedDouble();

    try {
      result = std::stod( (attr ? attr : elem->GetText()) );
    } catch (...) {
      log::ferror("When parsing '%s' of '%s', type double is required. '%s'",
                   child, parent->Name(), info);
    }

    return result;
  }


  /// \details If a node is enforced, it has to exist.
  int queryNodeInt(XMLElement *parent,
                   const std::string &child,
                   bool enforced,
                   const std::string &info) {

    if (!existNode(parent, child)) {
      if (enforced) {
        log::ferror("When parsing '%s' of '%s', attribute/sub-element "
                    "may not exist. '%s'", child, parent->Name(), info);
      }
      else {
        return defaultQueriedInt();
      }
    }

    auto attr = parent->Attribute(child.c_str());
    auto elem = parent->FirstChildElement(child.c_str());
    double result = defaultQueriedInt();

    try {
      result = std::stoi( (attr ? attr : elem->GetText()) );
    } catch (...) {
      log::ferror("When parsing '%s' of '%s', type int is required. '%s'",
                   child, parent->Name(), info);
    }

    return result;
  }


  const char* defaultQueriedString() {
    return nullptr;
  }


  double defaultQueriedDouble() {
    return std::numeric_limits<double>::infinity();
  }


  int defaultQueriedInt() {
    return std::numeric_limits<int>::min();
  }


} // inline namespace xmlutils


} // namespace antmoc

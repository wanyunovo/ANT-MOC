/// \file include/xml_utils.h
/// \date Sep 2, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef XML_UTILS_H_
#define XML_UTILS_H_

#include <string>

#include <tinyxml2.h>

namespace antmoc
{

// Alias for tinyxml2 types
using tinyxml2::XMLElement;
using tinyxml2::XMLDocument;
using tinyxml2::XMLError;
using tinyxml2::XMLPrinter;


///---------------------------------------------------------------------
/// \namespace xmlutils
/// \details Useful XML interfaces for antmoc.
///---------------------------------------------------------------------
inline namespace xmlutils {

  /// \brief Load an XML file.
  /// \param doc      XML document to load the file.
  /// \param file     Path to the XML file.
  /// \param enforced Treat any error as exception (true) or not (false).
  void loadFileXML(XMLDocument &doc, const std::string &file, bool enforced = false);

  /// \brief Save an XML document to a file.
  /// \param doc  XML document.
  /// \param file Path to the XML file.
  void saveFileXML(XMLDocument &doc, const std::string &file);

  template <typename T>
  XMLElement* queryFirstChild(T *,
                              const std::string &,
                              bool enforced = false,
                              const std::string &info = "");
  // Explicit instantiation
  extern template XMLElement* queryFirstChild(XMLDocument*, const std::string&,
                                              bool, const std::string&);
  extern template XMLElement* queryFirstChild(XMLElement*, const std::string&,
                                              bool, const std::string&);

  const char* queryStringAttr(XMLElement *,
                              const std::string &,
                              const std::string &,
                              bool enforced = false);

  double queryDoubleAttr(XMLElement *,
                         const std::string &,
                         const std::string &,
                         bool enforced = false);

  int queryIntAttr(XMLElement *,
                   const std::string &,
                   const std::string &,
                   bool enforced = false);


  /// \brief Check if a node exists. A node is either an attribute or
  ///        a sub-element.
  /// \param parent the parent element
  /// \param child name of the node
  bool existNode(const XMLElement *parent, const std::string &child);


  /// \brief Extract a string attribute/sub-element from a given element
  /// \param parent the parent element
  /// \param child name of the attribute or element to be read
  /// \param enforced determine if the entry is mandatory
  /// \param info extra info attached to error messages
  /// \return queried value if succeed, or nullptr if failed
  const char* queryNodeString(XMLElement *parent,
                              const std::string &child,
                              bool enforced = false,
                              const std::string &info = "");


  /// \brief Extract a floating point attribute/sub-element from a given element
  /// \param parent the parent element
  /// \param child name of the attribute or element to be read
  /// \param enforced determine if the entry is mandatory
  /// \param info extra info attached to error messages
  /// \return queried value if succeed, or infinity if failed
  double queryNodeDouble(XMLElement *parent,
                         const std::string &child,
                         bool enforced = false,
                         const std::string &info = "");


  /// \brief Extract a floating point attribute/sub-element from a given element
  /// \param parent the parent element
  /// \param child name of the attribute or element to be read
  /// \param enforced determine if the entry is mandatory
  /// \param info extra info attached to error messages
  /// \return queried value if succeed, or min() of int if failed
  int queryNodeInt(XMLElement *parent,
                   const std::string &child,
                   bool enforced = false,
                   const std::string &info = "");


  /// \brief The default value when a String query fails
  const char* defaultQueriedString();

  /// \brief The default value when a floating-point number query fails
  double defaultQueriedDouble();

  /// \brief The default value when an integer number query fails
  int defaultQueriedInt();


} // inline namespace xmlutils


} // namespace antmoc

#endif  // XML_UTILS_H_

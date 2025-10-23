/**
 * @file test_MaterialIdConversion
 * @brief Test material id conversion
 * @date April 24, 2019
 * @author Gen Wang, USTB
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#include "testing/test_utils.h"
#include "antmoc/Material.h"

using namespace antmoc;
namespace{

// Test fixture
class test_MaterialIdConversion : public testing::Test{
  // do nothing
};

/* TEST Examples in the Test fixture */
TEST_F(test_MaterialIdConversion,convertStringIdtoInt32Id){
  ROOT_ONLY();

  EXPECT_EQ(1073741824, Material::convertStringToId("0-0"));
  EXPECT_EQ(1073741825, Material::convertStringToId("0-1"));
  EXPECT_EQ(1074790400, Material::convertStringToId("1-0"));
  EXPECT_EQ(1074790401, Material::convertStringToId("1-1"));
  EXPECT_EQ(1084227584, Material::convertStringToId("10-0"));
  EXPECT_EQ(1084227585, Material::convertStringToId("10-1"));
  EXPECT_EQ(1084227594, Material::convertStringToId("10-10"));
  EXPECT_EQ(1085276170, Material::convertStringToId("11-10"));
}


TEST_F(test_MaterialIdConversion,convertInt32IdtoString){
  ROOT_ONLY();

  EXPECT_EQ("0-0", Material::convertIdToString(1073741824));
  EXPECT_EQ("0-1", Material::convertIdToString(1073741825));
  EXPECT_EQ("1-0", Material::convertIdToString(1074790400));
  EXPECT_EQ("1-1", Material::convertIdToString(1074790401));
  EXPECT_EQ("10-0", Material::convertIdToString(1084227584));
  EXPECT_EQ("10-1", Material::convertIdToString(1084227585));
  EXPECT_EQ("10-10", Material::convertIdToString(1084227594));
  EXPECT_EQ("11-10", Material::convertIdToString(1085276170));
}


//TEST_F(test_MaterialIdConversion,convertStringWithManyDelimeters){
//  ROOT_ONLY();
//
//  EXPECT_EQ(0, Material::convertStringToId("0-0-0"));
//  EXPECT_EQ(1, Material::convertStringToId("0-0-1"));
//  EXPECT_EQ(100, Material::convertStringToId("0-1-0"));
//  EXPECT_EQ(10101, Material::convertStringToId("1-1-1"));
//}

TEST_F(test_MaterialIdConversion,convertStringIdtoInt32IdWithoutDash){
  ROOT_ONLY();

  EXPECT_EQ(0, Material::convertStringToId("0"));
  EXPECT_EQ(1212, Material::convertStringToId("1212"));
}

TEST_F(test_MaterialIdConversion, testRobustness){
  ROOT_ONLY();

  EXPECT_EQ(-1, Material::convertStringToId("1*"));
  EXPECT_EQ(-1, Material::convertStringToId("A1"));
  EXPECT_EQ(-1, Material::convertStringToId("A"));
  EXPECT_EQ(-1, Material::convertStringToId("1A"));
  EXPECT_EQ(-1, Material::convertStringToId("-1A"));
  EXPECT_EQ(-1, Material::convertStringToId("1-A1"));

  const char *c_str = "11-10";
  EXPECT_EQ(1085276170, Material::convertStringToId(c_str));

}

} /*namespace*/













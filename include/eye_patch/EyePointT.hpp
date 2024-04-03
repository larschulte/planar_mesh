#pragma once
#include <pcl/point_types.h> // for using macros PCL_ADD_

struct EyePointT{
    PCL_ADD_POINT4D
    int id;
    
    EyePointT(){
        x = 0;
        y = 0;
        z = 0;
        id = 0;
    }
};

// need this for pcl/conversion to work
POINT_CLOUD_REGISTER_POINT_STRUCT (EyePointT,
                                   (float, x, x)
                                   (float, y, y)
                                   (float, z, z)
                                   (int, id, id)
)



// ================================================================================



/** NOTE
 * - include <pcl/point_types.h> after <pcl/impl/point_types.hpp> causes error
 *   thus here include <pcl/point_types.h> directly
 * 
 * - use array as oppose to Eigen::Matrix for minimum data size 
 * 
 * - inline means to replace the caller code with the code inside{}
 */


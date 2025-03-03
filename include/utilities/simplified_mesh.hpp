#include <pcl/PolygonMesh.h>
#include "MeshObject/Surface.hpp"

pcl::PolygonMesh create_simplified_mesh(const std::shared_ptr<Surface>& surface, bool with_color = false);

template <typename PointT>
pcl::PolygonMesh create_simplified_mesh_impl(const std::shared_ptr<Surface>& surface);

void merge_polygon_mesh(pcl::PolygonMesh& mesh1, pcl::PolygonMesh& mesh2);
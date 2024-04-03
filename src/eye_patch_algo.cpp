#include "eye_patch/eye_patch_algo.hpp"

using namespace eye_patch;

EyePatchAlgo::EyePatchAlgo(const EyePatchAlgoParams &p_algo):p_algo_(p_algo)
{
    polarSearch_data_ = boost::make_shared<pcl::PointCloud<EyePointT>>();
    polarSearch_.setCloudPointer(polarSearch_data_);
    polarSearch_.setResolution(p_algo_.projection_search_resolution_h, p_algo_.projection_search_resolution_v);
}

EyePatchAlgo::~EyePatchAlgo() {}

void EyePatchAlgo::add_pointcloud(const pcl::PointCloud<VilensPointT> &callback_lidar_pointcloud, const Eigen::Isometry3f &tf_lidar2world_eigen_f)
{   
    // transform to global
    pcl::PointCloud<EyePointT> surfel_pointcloud_global_frame;
    pcl::transformPointCloudWithNormals(surfel_pointcloud_local_frame, surfel_pointcloud_global_frame, tf_lidar2world_eigen_f);
    std::cout << "eye_patch - transform surfels to global frame \n";

    // add to octree
    addSurfelCloudToEyePatch(surfel_pointcloud_global_frame, tf_lidar2world_eigen_f.inverse());
    std::cout << "eye_patch - points added to global map \n";
}

// add to global map
void EyePatchAlgo::addSurfelCloudToEyePatch(pcl::PointCloud<EyePointT>& input_surfel_pointcloud, const Eigen::Isometry3f& tf_world2lidar_eigen_f)
{
    polarSearch_.prepareSearch(tf_world2lidar_eigen_f);

    // ------ [LOOP FOR L] ------
    for (auto iter = input_surfel_pointcloud.begin(); iter < input_surfel_pointcloud.end(); iter++)
    {
        // surfel l
        EyePointT& surfel_l = *iter;
        EyePointT surfel_l_backup = surfel_l; // copy surfel l -> overlapped using copied surfel l to maximum its influenced area)

        // radius search
        pcl::Indices nearby_old_points_indices = polarSearch_.extractNearbyIndices(surfel_l, tf_world2lidar_eigen_f);
        
        // check number of returns
        if (nearby_old_points_indices.size() == 0)
        {
            polarSearch_.addPoint(surfel_l, tf_world2lidar_eigen_f);
            continue;
        }

        // ------ [LOOP FOR K] ------
        for (int index : nearby_old_points_indices)
        {
            EyePointT &surfel_k = polarSearch_data_->at(index); // surfel k

             // skip invalid k
            if (surfel_k.valid == 0.f)
            {
                continue;
            }

            // skip k that is too far from l
            float normal_distance = abs((surfel_l.getVector3fMap() - surfel_k.getVector3fMap()).dot(surfel_l.getNormalVector3fMap()));
            if (normal_distance > p_algo_.projection_search_normal_distance)
            {
                continue;
            }

            // overlap surfels
            EyePointT _;
            EyePointT new_surfel_o;
            overlapSurfels(surfel_l, surfel_k, surfel_l, _, _); // update original surfel_l
            overlapSurfels(surfel_l_backup, surfel_k, _, surfel_k, new_surfel_o); // update surfel_k and surfel_overlapped 

            if (new_surfel_o.valid == 1.f) 
            {
                polarSearch_.addPoint(new_surfel_o, tf_world2lidar_eigen_f);
            }
        }

        // obtain color of l
        if (surfel_l.valid == 1.f)
        {
            polarSearch_.addPoint(surfel_l, tf_world2lidar_eigen_f);
        }
    }
}


void EyePatchAlgo::getResults(pcl::PointCloud<EyePointT>::Ptr output_surfel_pointcloud_ptr)
{
    // extract valid surfels
    pcl::PassThrough<EyePointT> pass;
    pass.setInputCloud(polarSearch_data_);
    pass.setFilterFieldName ("valid");
    pass.setFilterLimits(0.9f, 1.1f);
    pass.filter(*output_surfel_pointcloud_ptr);

    // update polarsearch data
    *polarSearch_data_ = *output_surfel_pointcloud_ptr;

    std::cout << "eye_patch - get result of " << output_surfel_pointcloud_ptr->size() << " surfels" << std::endl;
}
# read pcd file in folder and publish to ros topic

import rospy
from sensor_msgs.msg import PointCloud2
from sensor_msgs.msg import PointField
import sensor_msgs.point_cloud2 as pc2
import numpy as np
import os
from pypcd4 import PointCloud as pypcd4
from pypcd4.pointcloud2 import NPTYPE_TO_PFTYPE

import tf2_ros
import geometry_msgs.msg

def file_names_from_folder(folder_path):
    file_names = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            file_names.append(file)
    return file_names

def time_path_dict_from_vilens_file_names(folder_path):
    file_names = file_names_from_folder(folder_path)
    pcd_path_dict = {}
    for file_name in file_names:
        # read time from file name
        # cloud_1711460868_234132000.pcd
        time = file_name.split("_")[1] + " " + file_name.split("_")[2].split(".")[0]

        # Construct the file path
        file_path = os.path.join(folder_path, file_name)

        # add to dictionary
        pcd_path_dict[time] = file_path

    return pcd_path_dict

def time_path_dict_from_pcl_ros_file_names(folder_path):
    file_names = file_names_from_folder(folder_path)
    pcd_path_dict = {}
    for file_name in file_names:
        # read time from file name
        # 1711460868.234132000.pcd
        time = file_name.split(".")[0] + " " + file_name.split(".")[1]

        # Construct the file path
        file_path = os.path.join(folder_path, file_name)

        # add to dictionary
        pcd_path_dict[time] = file_path

    return pcd_path_dict

def lines_from_file(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()
    return lines

def time_pose_dict_from_slam_pose_graph(file_path):
    # get lines
    lines = lines_from_file(file_path)

    # dictionary to store pose and time
    pose_dict = {}
    for line in lines:
        data = line.split(" ")

        # VERTEX_SE3:QUAT_TIME id x y z qx qy qz qw sec nsec
        if data[0] == "VERTEX_SE3:QUAT_TIME":
            
            # pose
            pose = [float(data[2]), float(data[3]), float(data[4]), float(data[5]), float(data[6]), float(data[7]), float(data[8])]

            # time as string, without newline character
            time = data[9] + " " + data[10][:-1] # :-1 to remove new line character

            # add to dictionary
            pose_dict[time] = pose
    
    return pose_dict

def time_pose_dict_from_state_csv(file_path):
    # get lines
    lines  = lines_from_file(file_path)

    #   ADD_TO_CSV(stream, sec);
    #   ADD_TO_CSV(stream, nsec);
    #   ADD_TO_CSV(stream, state.W_r_WB[0]);
    #   ADD_TO_CSV(stream, state.W_r_WB[1]);
    #   ADD_TO_CSV(stream, state.W_r_WB[2]);
    #   ADD_TO_CSV(stream, state.q_WB.x());
    #   ADD_TO_CSV(stream, state.q_WB.y());
    #   ADD_TO_CSV(stream, state.q_WB.z());
    #   ADD_TO_CSV(stream, state.q_WB.w());
    #   ADD_TO_CSV(stream, state.B_v_WB[0]);
    #   ADD_TO_CSV(stream, state.B_v_WB[1]);
    #   ADD_TO_CSV(stream, state.B_v_WB[2]);
    #   ADD_TO_CSV(stream, state.B_w_WB[0]);
    #   ADD_TO_CSV(stream, state.B_w_WB[1]);
    #   ADD_TO_CSV(stream, state.B_w_WB[2]);
    #   ADD_TO_CSV(stream, state.biasAcc[0]);
    #   ADD_TO_CSV(stream, state.biasAcc[1]);
    #   ADD_TO_CSV(stream, state.biasAcc[2]);
    #   ADD_TO_CSV(stream, state.biasGyr[0]);
    #   ADD_TO_CSV(stream, state.biasGyr[1]);
    #   ADD_TO_CSV(stream, state.biasGyr[2]);
    #   ADD_TO_CSV(stream, state.tBiasAng[0]);
    #   ADD_TO_CSV(stream, state.tBiasAng[1]);
    #   ADD_TO_CSV(stream, state.tBiasAng[2]);
    #   ADD_TO_CSV(stream, state.tBiasLin[0]);
    #   ADD_TO_CSV(stream, state.tBiasLin[1]);
    #   ADD_TO_CSV_EOL(stream, state.tBiasLin[2]);

    # dictionary to store pose and time
    pose_dict = {}
    for line in lines:
        data = line.split(",")

        # sec,nsec,x,y,z,qx,qy,qz,qw, ...
        
        # time
        time = data[0] + " " + data[1]

        # pose
        pose = [float(data[2]), float(data[3]), float(data[4]), float(data[5]), float(data[6]), float(data[7]), float(data[8])]

        # add to dictionary
        pose_dict[time] = pose
    
    return pose_dict

def common_time_from_dicts(dict1, dict2):
    common_time = []
    
    # get keys
    keys1 = dict1.keys()
    keys2 = dict2.keys()

    # sort keys
    keys1 = sorted(keys1)
    keys2 = sorted(keys2)

    # find common time
    for time in keys1:
        if time in keys2:
            common_time.append(time)
    return common_time
    
    
def pointcloud2_msg_from_pcd(pcd_path, secs, nsecs, frame_id):
    # use pypcd4 object as intermediate
    pypcd4_object = pypcd4.from_path(pcd_path)

    # header
    header = rospy.Header()
    header.stamp.secs = secs
    header.stamp.nsecs = nsecs
    header.frame_id = frame_id

    # fields
    fields = []
    offset = 0
    for i, name in enumerate(pypcd4_object.fields):
        fields.append(PointField(name, offset, NPTYPE_TO_PFTYPE[np.dtype(pypcd4_object.types[i])], pypcd4_object.count[i]))
        offset += np.dtype(pypcd4_object.types[i]).itemsize * pypcd4_object.count[i]

    
    # data
    data = pypcd4_object.pc_data

    # return msg
    return pc2.create_cloud(header, fields, data)

def secs_nsecs_from_path(file_path):
    # /home/jiahao/Desktop/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460854_847538000.pcd
    file_name = file_path.split("/")[-1]
    secs = int(file_name.split("_")[1])
    nsecs = int(file_name.split("_")[2].split(".")[0])
    return secs, nsecs

def secs_nsecs_from_time(time):
    # 1711460854 847538000
    secs = int(time.split(" ")[0])
    nsecs = int(time.split(" ")[1])
    return secs, nsecs

def tf_from_pose(pose, secs, nsecs, frame_id, child_frame_id):
    t = geometry_msgs.msg.TransformStamped()
    t.header.frame_id = frame_id
    t.header.stamp.secs = secs
    t.header.stamp.nsecs = nsecs
    t.child_frame_id = child_frame_id
    t.transform.translation.x = pose[0]
    t.transform.translation.y = pose[1]
    t.transform.translation.z = pose[2]
    t.transform.rotation.x = pose[3]
    t.transform.rotation.y = pose[4]
    t.transform.rotation.z = pose[5]
    t.transform.rotation.w = pose[6]
    return t

def static_tf_from_pose(pose, frame_id, child_frame_id):
    t = geometry_msgs.msg.TransformStamped()
    t.header.frame_id = frame_id
    t.child_frame_id = child_frame_id
    t.transform.translation.x = pose[0]
    t.transform.translation.y = pose[1]
    t.transform.translation.z = pose[2]
    t.transform.rotation.x = pose[3]
    t.transform.rotation.y = pose[4]
    t.transform.rotation.z = pose[5]
    t.transform.rotation.w = pose[6]
    return t

# main
if __name__ == "__main__":
    
    # Initialize ROS node
    rospy.init_node("vilens_results_publisher", anonymous=True)
    pub_pc = rospy.Publisher("/vilens_pointcloud", PointCloud2, queue_size=1, latch=True)
    pub_tf = tf2_ros.TransformBroadcaster()


    # dataset 
    # (can use "rosrun pcl_ros bag_to_pcd *.bag folder/" to convert bag to pcd files as well)
    folder_path = "/home/jiahao/Desktop/osney power station/2024-03-26_13-47-27_rec004_osney_power_station_raw/pcd"
    time_pcd_path_dict = time_path_dict_from_pcl_ros_file_names(folder_path)
    state_csv_path = "/home/jiahao/Desktop/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/odometry/state.csv"
    time_pose_dict = time_pose_dict_from_state_csv(state_csv_path)
    # addtional base to lidar required, if using state.csv as pose since it records map to base not map to lidar
    # the transform is obtained from runtime_nl/config/rad002/vilens/vilens_rad002.yaml
    # B_r_BL: [-0.015152, 0.023212, 0.082800] # temp - same as XT32
    # q_BL: [ 0., 0., 0.7071068, 0.7071068 ] # from the URDF
    pub_static_tf = tf2_ros.StaticTransformBroadcaster()
    tf_BL = [-0.015152, 0.023212, 0.082800, 0., 0., 0.7071068, 0.7071068]
    msg_BL_tf = static_tf_from_pose(tf_BL, "base", "lidar")
    pub_static_tf.sendTransform(msg_BL_tf)

    # # dataset
    # folder_path = "/home/jiahao/Desktop/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds"
    # time_pcd_path_dict = time_path_dict_from_vilens_file_names(folder_path) 
    # slam_file_path = "/home/jiahao/Desktop/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_pose_graph.slam"
    # time_pose_dict = time_pose_dict_from_slam_pose_graph(slam_file_path)
    
    # common time
    common_time_list = common_time_from_dicts(time_pcd_path_dict, time_pose_dict)
    
    # publish
    rate = rospy.Rate(1)
    i = 0
    while not rospy.is_shutdown():
                
        # time
        time = common_time_list[i]
        secs, nsecs = secs_nsecs_from_time(time)

        # print 
        print(i, time)
        
        # publish pointcloud
        pcd_path = time_pcd_path_dict[time]
        msg_pc = pointcloud2_msg_from_pcd(pcd_path, secs, nsecs, "lidar")
        pub_pc.publish(msg_pc)

        # publish pose
        pose = time_pose_dict[time]
        msg_tf = tf_from_pose(pose, secs, nsecs, "map", "base") # change this along with changing dataset
        # msg_tf = tf_from_pose(pose, secs, nsecs, "map", "lidar") # change this along with changing dataset
        pub_tf.sendTransform(msg_tf)

        # control rate
        rate.sleep()
        i += 1
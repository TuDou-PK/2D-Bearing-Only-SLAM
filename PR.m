close all
clear
clc

addpath "./tools/g2o_wrapper"
addpath "./tools/visualization"
source "basicFunctions.m"
source "LS.m"
source "./tools/utilities/geometry_helpers_2d.m"

#------------------------------------------------Loading and pre-processing data------------------------------------------------------------#
# Loading data
[landmarks_gt, poses_gt, transitions_gt, observations_gt] = loadG2o("02-BearingOnlySLAM/slam2D_bearing_only_ground_truth.g2o");
[landmarks_ig, poses_ig, transitions_ig, observations_ig] = loadG2o("02-BearingOnlySLAM/slam2D_bearing_only_initial_guess.g2o");

# Pre-processing ground truth data
landmarks_gt = landmarks_gt(2: end);
poses_gt = poses_gt(2: end);
transitions_gt  = transitions_gt(2: end);
observations_gt = observations_gt(2: end);

#Pre-processing initial guess data
landmarks_ig = landmarks_ig(2: end);
poses_ig = poses_ig(2: end);
transitions_ig  = transitions_ig(2: end);
observations_ig = observations_ig(2: end);
#-------------------------------------------------------------------------------------------------------------------------------------------#



#----------------------------------------------------Prepare for initialing data------------------------------------------------------------#
#-----------------Get length of data--------------------------

length_poses_gt = length(poses_gt);	#length of ground truth(gt) poses data
length_poses_ig = length(poses_ig);	#length of initial guess(ig) poses data

length_ldmk = length(landmarks_gt);	#length of gt landmark data

length_obs_ig = 0;			#length of ig observation data
for i = 1: length(observations_ig)
	length_obs_ig += length(observations_ig(i).observation);
endfor

length_trans_ig = length(transitions_ig);	#length of ig transitions data


#---------------Initial array for storing data-----------------
#pose data
pose_array_gt = zeros(3, 3, length_poses_gt);
pose_array_ig = zeros(3, 3, length_poses_ig);

pose_id2index_gt = zeros(1, length_poses_gt);
pose_id2index_ig = zeros(1, length_poses_ig);

#landmark data
ldmk_array = zeros(2, length_ldmk);
ldmk_id2index_gt = zeros(1, length_ldmk);

#Observation data
obs_array_ig = zeros(1, length_obs_ig);
association_obs_ig = zeros(2, length_obs_ig);

trans_array_ig = zeros(3, 3, length_trans_ig);
association_trans_ig = zeros(2, length_trans_ig);
#-------------------------------------------------------------------------------------------------------------------------------------------#



#-----------------------------------------------------------Process pose data---------------------------------------------------------------#
for i = 1: length_poses_gt
	pose_array_gt(:, :, i) = v2t([poses_gt(i).x, poses_gt(i).y, poses_gt(i).theta]);
	pose_id2index_gt(i) = poses_gt(i).id;
	
	pose_array_ig(:, :, i) = v2t([poses_ig(i).x, poses_ig(i).y, poses_ig(i).theta]);
	pose_id2index_ig(i) = poses_ig(i).id;
endfor
#-------------------------------------------------------------------------------------------------------------------------------------------#



#--------------------------------------------------------process landmark data--------------------------------------------------------------#
for i = 1:length_ldmk
	ldmk_array(:, i) = [landmarks_gt(i).x_pose, landmarks_gt(i).y_pose];
	ldmk_id2index_gt(i) = landmarks_gt(i).id;
endfor
#-------------------------------------------------------------------------------------------------------------------------------------------#





#--------------------------------------------------------process observation data-----------------------------------------------------------#
k = 1;
for i = 1: length(observations_ig)     #for all observation times
	observations = observations_ig(i);   # for each observation
	pose_id = observations.pose_id;      # get id
	pose_index = find(pose_id2index_ig == pose_id);    # pose_gt's id index is same with pose id
	num_observations = length(observations.observation);    # get number of observations for each time

	for j = 1: num_observations
		obs = observations.observation(j);    
		ldmk_id = obs.id;
		ldmk_index = find(ldmk_id2index_gt == ldmk_id);
		obs_array_ig(k) = obs.bearing;
		association_obs_ig(:, k) = [pose_index; ldmk_index];
		k++;
	endfor
endfor
#-------------------------------------------------------------------------------------------------------------------------------------------#


#--------------------------------------------------------Odometry data----------------------------------------------------------------------#
m = 1;
for i = 1: length_trans_ig
	transition = transitions_ig(i);
	pose_from_id = transition.id_from;
	pose_to_id = transition.id_to;
	pose_from_index = find(pose_id2index_ig == pose_from_id);
	pose_to_index = find(pose_id2index_ig == pose_to_id);
	association_trans_ig(:, i) = [pose_from_index; pose_to_index];
	trans_array_ig(:, :, i) = v2t(transition.v);
endfor
#-------------------------------------------------------------------------------------------------------------------------------------------#


#-------------------------------------------------------------Get landmarks-----------------------------------------------------------------#
ldmk_ids = unique(association_obs_ig(2, :));
ldmk_nums = zeros(size(ldmk_ids));
num_ldmk_ids = length(ldmk_ids);

for i = 1: num_ldmk_ids
	ldmk_id = ldmk_ids(i);
	index_of_id = find(association_obs_ig(2, :) == ldmk_id);
	bearings = obs_array_ig(index_of_id);	
	num_bearings = length(bearings);
	pose_from_id = association_obs_ig(1, :)(index_of_id);	
	poses_array = zeros(3, length(bearings));
	
	for j = 1: num_bearings
		poses_array(:, j) = t2v(pose_array_ig(:, :, pose_from_id(j)));		
		bearings(j) = bearings(j) + poses_array(3, j);   
	endfor	
	
	[pose1, pose2] = findPose(num_bearings, bearings);                                            

	if num_bearings == 1
		ldmk_array_ig(:, ldmk_id) = poses_array(1:2, 1) + [cos(poses_array(3, 1)); sin(poses_array(3, 1))];
	else
		p11 = poses_array(1:2, pose1);
 		theta1 = bearings(pose1);
 		p12 = p11 + [cos(theta1); sin(theta1)];
 		
 		p21 = poses_array(1:2, pose2);
 		theta2 = bearings(pose2);
 		p22 = p21 + [cos(theta2); sin(theta2)];
 		
 		ldmk_array_ig(:, ldmk_id) = intersect_lines(p11, p12, p21, p22);
	endif
	
endfor

#-------------------------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------Least Square Algorithm------------------------------------------------------------#
iterations = 5;
damping = 1e-5;

[pose_array_pr, ldmk_array_pr] = LS(pose_array_ig, ldmk_array_ig, obs_array_ig, association_obs_ig, trans_array_ig, association_trans_ig, damping, iterations);
disp("Program run complete!");
#-------------------------------------------------------------------------------------------------------------------------------------------#


#------------------------------------------------------------Result Figures-----------------------------------------------------------------#
f1 = figure(1);
position=get(f1,"position");
set(f1,"position",[10 10 2000 900]);


subplot(1, 2, 1);
plot(ldmk_array(1,:),ldmk_array(2,:),'b*');
hold on;
plot(squeeze(pose_array_gt(1,3,:)),squeeze(pose_array_gt(2,3,:)),'b*-');
hold on;
plot(ldmk_array_ig(1,:),ldmk_array_ig(2,:),'m*');
hold on;
plot(squeeze(pose_array_ig(1,3,:)),squeeze(pose_array_ig(2,3,:)),'m*-');
hold on;
title("Grand Truth & Initial Guess", "fontsize", 16);
l1 = legend("GT landmarks", "GT poses", "IG landmarks", "IG poses");
set(l1, "fontsize", 12);


subplot(1, 2, 2);
plot(ldmk_array(1,:),ldmk_array(2,:),'b*');
hold on;
plot(squeeze(pose_array_gt(1,3,:)),squeeze(pose_array_gt(2,3,:)),'b*-');
hold on;
plot(ldmk_array_pr(1,:),ldmk_array_pr(2,:),'r*');
hold on;
plot(squeeze(pose_array_pr(1,3,:)),squeeze(pose_array_pr(2,3,:)),'r*-');
title("Grand Truth & Computed Trajectory", "fontsize", 16);
l2 = legend("GT landmarks", "GT poses", "CT landmarks", "CT poses");
set(l2, "fontsize", 12);
pause(100);



















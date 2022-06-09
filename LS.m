1;

source "basicFunctions.m"

function [pose_array_pr, ldmk_array_pr] = LS(pose_array_ig, ldmk_array_ig, obs_array_ig, association_obs_ig, trans_array_ig, association_trans_ig, damping, iterations)
    
    length_poses_ig = size(pose_array_ig)(3);	# size is (3, 3, 301), so Nr=301
    length_ldmk_ig = size(ldmk_array_ig)(2);	# size is (2, 141), so Nl = 141

    state_dim = 3 * length_poses_ig + 2 * length_ldmk_ig;

    disp("LS iterations start: ");

    for i=1:iterations
        printf("Iteration: %d\n", i);

        [H_ldmk, b_ldmk] = buildLinearSystemLandmarks(state_dim, size(trans_array_ig)(3), association_trans_ig, trans_array_ig, pose_array_ig, length_poses_ig, length_ldmk_ig);
        
	[H_pose, b_pose] = buildLinearSystemPoses(state_dim, size(obs_array_ig)(2), association_obs_ig, obs_array_ig, pose_array_ig, ldmk_array_ig, length_poses_ig, length_ldmk_ig);

        H=H_ldmk + H_pose + eye(state_dim,state_dim)*damping;
        b=b_ldmk + b_pose;

        H((length_poses_ig-1)*3+1:length_poses_ig*3,:)=[];
        H(:,(length_poses_ig-1)*3+1:length_poses_ig*3)=[];
        b((length_poses_ig-1)*3+1:length_poses_ig*3)=[];

        dx = zeros(state_dim,1);
        dx=H\(-b);
        dx = [dx(1:(length_poses_ig-1)*3); zeros(3,1); dx((length_poses_ig-1)*3+1:end)];
        [pose_array_pr, ldmk_array_pr] = box_plus(pose_array_ig, ldmk_array_ig, dx, length_poses_ig, length_ldmk_ig);
    endfor   
endfunction

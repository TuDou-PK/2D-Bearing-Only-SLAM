1;

global pose_dim=3;
global landmark_dim=2;

function A = v2t(v)
	c = cos(v(3));
	s = sin(v(3));

	A=[c, -s, v(1);
	   s,  c, v(2);
	   0,  0,   1];
endfunction;


function v = t2v(A)
	v(1:2,1) = A(1:2,3);
	v(3,1) = atan2(A(2,1), A(1,1));
endfunction;



function delta_theta = diff_(theta1, theta2)
	difference = theta2-theta1;
    if (difference>pi)
        delta_theta = abs(difference-2*pi);
    elseif (difference<-pi)
        delta_theta = abs(difference+2*pi);
    else
        delta_theta = abs(difference);
    endif
endfunction;


function [XR, XL] = box_plus(XR, XL, dX, Nr, Nl)
	pose_dim=3;
	landmark_dim=2;
	
        for i=1:Nr
                pose_i = poseMatrixIndex(i, Nr, Nl);
                dXr = dX(pose_i:pose_i+pose_dim-1,:);
                XR(:,:,i)=v2t(dXr)*XR(:,:,i);
        end
        for i=1:Nl
                land_i = landmarkMatrixIndex(i, Nr, Nl);
                dXl = dX(land_i:land_i+landmark_dim-1,:);
                XL(:,i)=dXl + XL(:,i);
        end
endfunction


function point = intersect_lines(p11,p12,p21,p22)
        k1 = (p12(2)-p11(2))/(p12(1)-p11(1));
        k2 = (p22(2)-p21(2))/(p22(1)-p21(1));
        w21 = p11(2)-k1*p11(1);
        w22 = p21(2)-k2*p21(1);
        x=(w22-w21)/(k1-k2);
        y=w21+k1*x;
        point = [x;y];
endfunction;


function [err, Ji, Jj]=poseErrorAndJacobian(Xi,Xj,Z)
        Ri=Xi(1:2,1:2);
        Rj=Xj(1:2,1:2);

        ti=Xi(1:2,3);
        tj=Xj(1:2,3);

        dR_0=[0 -1;1 0];

        Z_hat=eye(3);
        Z_hat(1:2,1:2)=Ri'*Rj;
        Z_hat(1:2,3)=Ri'*(tj-ti);

        err = [(Z_hat-Z)(1:2,1); (Z_hat-Z)(1:2,2); (Z_hat-Z)(1:2,3)];

        Ji=zeros(6,3);
        Jj=zeros(6,3);

        Jj(1:4,3)=reshape(Ri'*dR_0*Rj, 4, 1);
        
        Jj(5:6,1:2)=Ri';

        Jj(5:6,3)=Ri'*dR_0*tj;
        Ji=Ji-Jj;
endfunction


function [err, Jr_i, Jl_j] = bearingErrorAndJacobian(Xr_i, Xl_j, z);
       R=Xr_i(1:2,1:2);
       dR_0=[0, -1;1, 0];
       t=Xr_i(1:2,3);

       l_hat = R'*(Xl_j-t);
       z_hat = atan2(l_hat(2),l_hat(1));

       err=z_hat-z;
       err=atan2(sin(err),cos(err));

       Jr_i = (1./(l_hat(1:2)'*l_hat(1:2))*[-l_hat(2) l_hat(1)]) * [-R', R'*dR_0'*Xl_j];
       Jl_j = (1./(l_hat(1:2)'*l_hat(1:2))*[-l_hat(2) l_hat(1)]) * R';
endfunction

function [pose1, pose2] = findPose(num_bearings, bearings)
	diff_angle = zeros(num_bearings, num_bearings);
	for i = 1: num_bearings - 1
		for j = i + 1: num_bearings
			diff_angle(i, j) = diff_(bearings(i), bearings(j));
			if diff_angle(i, j) > pi/2
				diff_angle(i, j) = pi - diff_angle(i, j);
			endif	
		endfor
	endfor
	[pose1, pose2] = find(diff_angle == max(max(diff_angle)));
endfunction



function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)
	global pose_dim;
	global landmark_dim;

	if (pose_index>num_poses)
		v_idx=-1;
	return;
	endif;
	v_idx=1+(pose_index-1)*pose_dim;
endfunction



function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
	global pose_dim;
	global landmark_dim;
	if (landmark_index>num_landmarks)
		v_idx=-1;
		return;
	endif;
	v_idx=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
endfunction



function [H, b] = buildLinearSystemLandmarks(state_dim, length_trans_ig, association_trans_ig, trans_array_ig, pose_array_ig, length_poses_ig, length_ldmk_ig)
        H = zeros(state_dim, state_dim);
        b = zeros(state_dim,1);
        
        %handle poses
        for m=1:length_trans_ig
            pose_i = association_trans_ig(1,m);
            pose_j = association_trans_ig(2,m);
            z = trans_array_ig(:,:,m);
            X_i = pose_array_ig(:,:,pose_i);
            X_j = pose_array_ig(:,:,pose_j);
            [err, Ji, Jj] = poseErrorAndJacobian(X_i, X_j, z);

            pose_i_H=poseMatrixIndex(pose_i, length_poses_ig, length_ldmk_ig);
            pose_j_H=poseMatrixIndex(pose_j, length_poses_ig, length_ldmk_ig);

            H(pose_i_H:pose_i_H+2, pose_i_H:pose_i_H+2)+=Ji'*Ji;
            H(pose_i_H:pose_i_H+2, pose_j_H:pose_j_H+2)+=Ji'*Jj;
            H(pose_j_H:pose_j_H+2, pose_i_H:pose_i_H+2)+=Jj'*Ji;
            H(pose_j_H:pose_j_H+2, pose_j_H:pose_j_H+2)+=Jj'*Jj;

            b(pose_i_H:pose_i_H+2)+=Ji'*err;
            b(pose_j_H:pose_j_H+2)+=Jj'*err;
        endfor
endfunction;



function [H, b] = buildLinearSystemPoses(state_dim, length_obs_ig, association_obs_ig, obs_array_ig, pose_array_ig, ldmk_array_ig, length_poses_ig, length_ldmk_ig)
        H = zeros(state_dim, state_dim);
        b = zeros(state_dim,1);
        #handle landmarks
        for m=1:length_obs_ig
            pose_i = association_obs_ig(1,m);
            land_j = association_obs_ig(2,m);
            z = obs_array_ig(1,m);
            Xr_i = pose_array_ig(:,:,pose_i);
            Xl_j = ldmk_array_ig(:,land_j);
            [err, Jr_i, Jl_j] = bearingErrorAndJacobian(Xr_i, Xl_j, z);

            pose_i_H=poseMatrixIndex(pose_i, length_poses_ig, length_ldmk_ig);
            land_j_H=landmarkMatrixIndex(land_j, length_poses_ig, length_ldmk_ig);

            H(pose_i_H:pose_i_H+2, pose_i_H:pose_i_H+2)+=Jr_i'*Jr_i;
            H(pose_i_H:pose_i_H+2, land_j_H:land_j_H+1)+=Jr_i'*Jl_j;
            H(land_j_H:land_j_H+1, land_j_H:land_j_H+1)+=Jl_j'*Jl_j;
            H(land_j_H:land_j_H+1, pose_i_H:pose_i_H+2)+=Jl_j'*Jr_i;

            b(pose_i_H:pose_i_H+2)+=Jr_i'*err;
            b(land_j_H:land_j_H+1)+=Jl_j'*err;
        endfor
endfunction

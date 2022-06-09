source "../utilities/geometry_helpers_3d.m"

%projects a world point onto an image plane
function [p_img, is_outside] = projectPoint(K,p) 
  p_img=K*p;
  p_img*=1./p_img(3);
  is_outside=(P(3)<0 || P_img(1)<0 ||P_img(1)>c || P_img(2)<0
             || P_img(2)>r)
endfunction;


%projects a set of world points onto the image plane
function [P_img, is_outside]=projectPoints(K,P)
  P=K*P;
  P_img=P(1:2,:)./P(3,:);
  r=2*K(1,3);
  c=2*K(2,3);
  is_outside= ((P(3,:)<0)
  | (P_img(1,:)<0)
  | (P_img(1,:)>c)
  | (P_img(2,:)<0)
  | (P_img(2,:)>r));
endfunction;

function [P1, P2] = pruneProjectedPointPairs(P1, is_outside1,  P2, is_outside2)
  if (size(P1,2)!=size(P2,2))
    disp("points do not match");
  endif;
  n=size(P1,2);
  k=1;
  for (i=1:n)
    if (is_outside1(i) || is_outside2(i))
      continue;
    endif;
    P1(:,k)=P1(:,i);
    P2(:,k)=P2(:,i);
    ++k;
  endfor;
  P1=P1(:,1:k);
  P2=P2(:,1:k);
endfunction;

% triangulates a point, passing through two lines
% one passing through the origin, and having
% direction vector d1
% one passing through a point p2, and having
% direction d2
function [success, p, e]=triangulatePoint(p2, d1, d2)
  p=zeros(3,1);
  success=false;
  e=-1;
                      
  D=[-d1, d2];         #assemble system matrix to find ascissa 
  s=-(D'*D)\(D'*p2);          #s: ascissa of closest point on p1 and p2
  if (s(1)<0 || s(2)<0)
    return;
  endif;
  success=true;
  p1_triangulated=d1*s(1);   # point on 1st line
  p2_triangulated=d2*s(2)+p2; # point on 2nd line
  e=norm(p1_triangulated-p2_triangulated); #difference between the points
  p=0.5*(p1_triangulated+p2_triangulated);               #midpoint
endfunction;

# triangulates a batch of points in image coords,
# X: is the pose of the world w.r.t the 2nd camera
function [n_success, P, errors] = triangulatePoints(K,X,P1_img,P2_img)
  #initialize vars
  n_points=size(P1_img,2);
  P=zeros(3,n_points);
  errors=zeros(1,n_points);
  n_success= 0;
  
  #inverse transform
  iX=inv(X);
  #inverse camera matrix
  iK=inv(K);
  #inverse rotation * inverse camera matix
  iRiK=iX(1:3,1:3)*iK;

  #express the points in camera coordinates
  #rotate the direction vector of P2 in world frame
  D1_cam=iK*[P1_img; ones(1,n_points)];
  D2_cam=iRiK*[P2_img; ones(1,n_points)];
  p2=iX(1:3,4);
  for (i=1:n_points)
    p1_cam=D1_cam(:,i);
    p2_cam=D2_cam(:,i);
    [success, p, e]=triangulatePoint(p2, p1_cam, p2_cam);
    if (success==true)
      ++n_success;
      P(:,n_success)=p;
      errors(n_success)=e;
    endif;
  endfor;
endfunction

function  E = transform2essential(X)
  E=X(1:3,1:3)'*skew(X(1:3,4));
endfunction;


function [X1, X2]=essential2transform(E)
  W=[0, -1,  0;
     1,  0,  0;
     0,  0,  1];

  [U,S,V]=svd(E);
  R1=V*W*U';
  if (det(R1)<0) #right handed condition
    [U,S,V]=svd(-E);
    R1=V*W*U';
  endif;
  #1st solution for the rotation
  X1=eye(4);
  X1(1:3,1:3)=R1;
  t_cross=R1*E;
  X1(1:3,4)= [t_cross(3,2)-t_cross(2,3);
              t_cross(1,3)-t_cross(3,1);
              t_cross(2,1)-t_cross(1,2)];

  #2nd solution for the rotation
  R2=V*W'*U';
  X2=eye(4);
  X2(1:3,1:3)=R2;
  t_cross=R2*E;
  X2(1:3,4)= [t_cross(3,2)-t_cross(2,3);
              t_cross(1,3)-t_cross(3,1);
              t_cross(2,1)-t_cross(1,2)];
endfunction;

function [X1,X2] = fundamental2transform(K,F)
  E=K'*F*K;
  [X1,X2] = essential2transform(E);
endfunction;

function F=transform2fundamental(K,X)
  iK = inv(K);
  F=iK'*transform2essential(X)*iK;
endfunction


#estimate fundamental matrix
function F = estimateFundamentalSimple(P1_img, P2_img)
  H=zeros(9,9);
  n_points=size(P1_img,2);
  for (i=1:n_points)
    p1_img=[P1_img(:,i); 1];
    p2_img=[P2_img(:,i); 1];
    A=reshape(p1_img*p2_img',1,9);
    H+=A'*A;
  endfor;
  [V,lambda]=eig(H);
  F=reshape(V(:,1),3,3);
endfunction

#computes a preconditioning matrix, that scales the points around the
#origin
#A(1:2,1:2): inverse sigma of the points
#A(1:2,3)  : -inverse sigma*mean
function A=computePreconditioningMatrix(P_img)
  n_points=size(P_img,2);
  P_img=[P_img; ones(1,n_points)];
  s=sum(P_img,2);
  s2=P_img*P_img';
  mu=s/n_points;
  sigma=s2/n_points-mu*mu';
  A=eye(3);
  A(1:2,1:2)=inv(chol(sigma(1:2,1:2)));
  A(1:2,3)=-A(1:2,1:2)*mu(1:2);
endfunction

#estimate fundamental matrix
function F = estimateFundamental(P1_img, P2_img, use_preconditioning=false)
  n_points=size(P1_img,2);
  A1=A2=eye(3);
  if (use_preconditioning)
    A1=computePreconditioningMatrix(P1_img);
    A2=computePreconditioningMatrix(P2_img);
  endif;
  AP1=A1(1:2,1:2)*P1_img+repmat(A1(1:2,3),1,n_points);
  AP2=A2(1:2,1:2)*P2_img+repmat(A2(1:2,3),1,n_points);
  Fa=estimateFundamentalSimple(AP1,AP2);
  F=A1'*Fa*A2;
endfunction

#estimates a transform from a set of image projections
#with camera matrix K
function X=estimateTransform(K,P1_img, P2_img, use_preconditioning=false)
  # compute fundamental
  F=estimateFundamental(P1_img,P2_img, use_preconditioning);

  # extract essential
  E=K'*F*K;

  #extract transforms from essential
  [X1,X2]=essential2transform(E);

  X=X1;
  n_in_front=0;
  #for each transform pick the best
  X_test=X1;
  [n_test, P]=triangulatePoints(K, X_test,P1_img, P2_img);
  n_test
  if (n_test>n_in_front)
    X=X_test;
    n_in_front=n_test;
  endif;
  X_test(1:3,4)=-X_test(1:3,4);
  [n_test, P]=triangulatePoints(K, X_test,P1_img, P2_img);
  if (n_test>n_in_front)
    X=X_test;
    n_in_front=n_test;
  endif;
  X_test=X2;
  [n_test, P]=triangulatePoints(K, X_test,P1_img, P2_img);
  if (n_test>n_in_front)
    X=X_test;
    n_in_front=n_test;
  endif;
  X_test(1:3,4)=-X_test(1:3,4);
  [n_test, P]=triangulatePoints(K, X_test,P1_img, P2_img);
  if (n_test>n_in_front)
    X=X_test;
    n_in_front=n_test;
  endif;
endfunction;

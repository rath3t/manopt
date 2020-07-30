function [detjacobian,invjacobian]=JacobianOfMesh(nnel,dshapedxi,dshapedeta,xcoord,ycoord)
%------------------------------------------------------------------------
%  Purpose:
%     determine the Jacobian for two-dimensional mapping
%
%  Synopsis:
%     [detjacobian,invjacobian]=Jacobian(nnel,dshapedxi,dshapedeta,xcoord,ycoord) 
%
%  Variable Description:
%     jacobian - Jacobian 
%     nnel - number of nodes per element   
%     dshapedxi - derivative of shape functions w.r.t. natural coordinate xi
%     dshapedeta - derivative of shape functions w.r.t. natural coordinate eta
%     xcoord - x axis coordinate values of nodes
%     ycoord - y axis coordinate values of nodes
%------------------------------------------------------------------------
 jacobian=zeros(2,2);
%  for i=1:nnel
 jacobian(1,1) = dshapedxi.'*xcoord;
 jacobian(1,2) = dshapedxi.'*ycoord;
 jacobian(2,1) = dshapedeta.'*xcoord;
 jacobian(2,2) = dshapedeta.'*ycoord;
%  end
 detjacobian = det(jacobian) ;  % Determinant of Jacobian matrix
 invjacobian = inv(jacobian) ;  % Inverse of Jacobian matrix
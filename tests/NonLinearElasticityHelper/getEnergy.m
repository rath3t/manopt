function [energy,store] = getEnergy(X,coordinates, nodes,nel,nnel,d,dirichnodes,mu,store)

energy=0.0;
% index=zeros(edof,1);                     % index vector


% Order of Gauss Quadrature
%--------------------------------------------------------------------------
[pointb,weightb]=  GaussQuadrature('second');     % sampling points & weights for bending


nglb=size(pointb,1);                     % 2x2 Gauss-Legendre quadrature for bending

node =zeros(nnel,1);
xx =zeros(nnel,1);
yy =zeros(nnel,1);
zz =zeros(nnel,1);

gradDisplacements= zeros(d,2);

nodeDisplacement = zeros(nnel,d);
nnodes=length(coordinates);
Dusefull=zeros(nnodes,d);
% if ~isfield(store, 'ready')
%     for iel=1:nel                        % loop for the total number of elements
nodecounter=1;
for nn=1:nnodes
    % loop for total number of nodes for single element
    
    %             node(i)=nodes(iel,i);               % extract connected node for (iel)^th element
    
    fixedorfree = find(dirichnodes==nn,1);
    if fixedorfree ~= 0
        Dusefull(nodecounter,:) =  X.displacementFixed(:,fixedorfree);
        %                 fixCounter=fixCounter+1;
        nodecounter=nodecounter+1;
    else
        A=dirichnodes>nn;
        idx=A==0;
        out=sum(idx(:));
        Dusefull(nodecounter,:) =  X.displacement(:,nn-out);
       nodecounter=nodecounter+1;
    end
    
    
end

store.Xusefull=Dusefull;



for iel=1:nel                        % loop for the total number of elements
    
    for i=1:nnel                        % loop for total number of nodes for single element
        
        node(i)=nodes(iel,i);               % extract connected node for (iel)^th element
        
        %         fixedorfree = find(dirichnodes==node(i),1);
        %         if fixedorfree ~= 0
        nodenodeDisplacements(i,:) = Dusefull(node(i),:);
        %                 fixCounter=fixCounter+1;
        %         else
        %             A=dirichnodes>node(i);
        %             idx=A==0;
        %             out=sum(idx(:));
        %             nodeDirectors(i,:) = X.Free(node(i)-out,:);
        %             %                 freeCounter=freeCounter+1;
        %         end
        
        elementstresses(1:d,1:d,i)= X.Stresses(1:d,1:d,node(i));
%         elementrightCauchy(1:d,1:d,i)= X.RightCauchyGreen{node(i)}(1:d,1:d);
        
        
        xx(i)=coordinates(node(i),1);       % extract x value of the node
        yy(i)=coordinates(node(i),2);       % extract y value of the node
        
    end
    elementrightCauchy= X.RightCauchyGreen{iel};
    
    %--------------------------------------------------------------------------
    %  numerical integration
    %--------------------------------------------------------------------------
    
    for int=1:nglb                        % nglb is sampling points =4
        xi=pointb(int,1);
        wt=weightb(int,1);
        eta=pointb(int,2);
        
        
        [shape,dshapedxi,dshapedeta]=Shapefunctions(xi,eta);
        % compute shape functions and derivatives at sampling point
        
        [detjacobian,invjacobian]=JacobianOfMesh(nnel,dshapedxi,dshapedeta,xx,yy);  % compute Jacobian
        
        [dshapehdx,dshapehdy]=ShapefunctionDerivatives(nnel,dshapedxi,dshapedeta,invjacobian);
        % derivatives w.r.t. physical coordinate
        
        gradx= nodenodeDisplacements.'*dshapehdx;
        grady= nodenodeDisplacements.'*dshapehdy;
        
        S=interpolate(shape,elementstresses);
%         C=interpolate(shape,elementrightCauchy);
            C=elementrightCauchy;
%         C=C/det(C)^(1/d);
        F=eye(d,d)+[gradx  grady];
        
        Cu=F.'*F;
        
        psi= mu/2*(trace(C)-2);
        langrange= 0.5*trace(S*(Cu-C));
        
%         RS=norm(Cu-C)*0.5
%         RC=0.5*S*F*
        
        
        

        energy=energy+  (psi+ langrange)*wt*detjacobian;
        
    end                      % end of numerical integration loop
    
    
    
end


end
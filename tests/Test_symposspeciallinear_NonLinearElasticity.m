function Test_symposspeciallinear_NonLinearElasticity()


clear all;
close all;
clc;
if ~exist('sizedom1', 'var') || isempty(sizedom1)
    % Dimension of the embedding space: R^d
    sizedom1 = 20;
end

if ~exist('sizedom2', 'var') || isempty(sizedom2)
    % Dimension of the embedding space: R^d
    sizedom2 = 20;
end
eleT1 = 2 ;
eleT2 = 2 ;

[coordinates,nodesMat]=createMeshRectanglular(sizedom1,sizedom2,eleT1,eleT2);

nel = length(nodesMat) ;                  % number of elements
nnel=4;                                % number of nodes per element
ndof=5;                                % number of dofs per node
nnode = length(coordinates) ;          % total number of nodes in system
sdof=nnode*ndof;                       % total system dofs
edof=nnel*ndof;                        % degrees of freedom per element

% nnode=nnodefree;
nu= 0.5;
E=1000;
mu =E/(2*(1+nu)); % AMat is a material parameter exhange energy coefficient



%fix all node at xl<x<xr and yb<y<yt
xl1=-inf;
xr1=inf;
yb1=-0.001;
yt1=0.001;

xl2=-0.001;
xr2=0.001;
yb2=-inf;
yt2=inf;

xl3=-0.001+sizedom1;
xr3=0.001+sizedom1;
yb3=-inf;
yt3=inf;


xl4=-inf;
xr4=inf;
yb4=-0.001+sizedom2;
yt4=0.001+sizedom2;


dirchnod=zeros(0);
dirchnodU=zeros(0);
dirchnodL=zeros(0);
dirchnodR=zeros(0);
dirchnodT=zeros(0);

activeOrInactive=[1 0 0 0];


freenod=zeros(0);
counterDirichline1=0;
counterDirichline2=0;
counterDirichline3=0;
counterDirichline4=0;
for curnnode=1:nnode
    xcur= coordinates(curnnode,1);
    ycur= coordinates(curnnode,2);
    if xl1<xcur && xr1>xcur && yb1<ycur && yt1 >ycur && activeOrInactive(1) ==1
        dirchnod= [dirchnod; curnnode];
        dirchnodU= [dirchnodU; curnnode];
        counterDirichline1=counterDirichline1+1;
    elseif  xl2<xcur && xr2>xcur && yb2<ycur && yt2 >ycur && activeOrInactive(2) ==1
        dirchnod= [dirchnod; curnnode];
        dirchnodL= [dirchnodL; curnnode];
        counterDirichline2=counterDirichline2+1;
    elseif  xl3<xcur && xr3>xcur && yb3<ycur && yt3 >ycur && activeOrInactive(3) ==1
        dirchnod= [dirchnod; curnnode];
        dirchnodR= [dirchnodR; curnnode];
        counterDirichline3=counterDirichline3+1;
    elseif  xl4<xcur && xr4>xcur && yb4<ycur && yt4 >ycur && activeOrInactive(4) ==1
        dirchnod= [dirchnod; curnnode];
        dirchnodT= [dirchnodT; curnnode];
        counterDirichline4=counterDirichline4+1;
    else
        freenod= [freenod; curnnode];
    end
end

ndirichnodes=size(dirchnod,1);
if ~exist('d', 'var') || isempty(d)
    % Dimension of the embedding space: R^d
    d = 2;
end




%     dirichletnode=
nnodefree=nnode-ndirichnodes;


A = zeros(ndirichnodes,d);
curnodecounter=1;

if d==2
    mU=[1,0];
    mL=[0,-1];
    mR = [0,1];
    mT = [-1,0];
else
    mU=[0,1,0];
    mL=[1,0,0];
    mR = [0,0,-1];
    mT = [0,0,1];
end



%     for curnod=1:size(dirchnod,1)
%         if any(dirchnod(curnod)==dirchnod1(:))
%             A(curnodecounter,:) = [-1,0,0];
%             curnodecounter=curnodecounter+1;
%         else
%             A(curnodecounter,:) = [0,0,1];
%             curnodecounter=curnodecounter+1;
%         end
%     end
% end

manifoldRightCauchyGreen =powermanifold(symposspeciallinear(d),nnode);

%Ein C pro Element
manifoldRightCauchyGreen =powermanifold(symposspeciallinear(d),nel);

% productmanifold(manifoldRightCauchyGreen);

manifoldStresses = symmetricfactory(d,nnode);
manifoldDisplacementFree= euclideanfactory(d, nnodefree);
manifoldFixed= constantfactory(zeros(d,ndirichnodes));

manifolds = struct();
manifolds.Stresses = manifoldStresses;
manifolds.RightCauchyGreen = manifoldRightCauchyGreen;
manifolds.displacement = manifoldDisplacementFree;
manifolds.displacementFixed = manifoldFixed;

manifold= productmanifold(manifolds);
% manifold = elliptopefactory(n, d);

% Generate a random initial guess if none was given.
if ~exist('X0', 'var') || isempty(X0)
    X0 = manifold.rand();
end

for i=1:nnode
    X0.Stresses(1:d,1:d,i)=zeros(d,d);
end

for i=1:nel
    X0.RightCauchyGreen{i}=eye(d,d);
end
for i=1:nnodefree
    X0.displacement(1:d,i)=zeros(d,1);
end


% Problem

nodecounter=1;

for nn=1:nnode
    % loop for total number of nodes for single element
    
    %             node(i)=nodes(iel,i);               % extract connected node for (iel)^th element
    
    fixedorfree = sum(find(dirchnod==nn,1));
    if fixedorfree ~= 0
        displacementsUseFull(nodecounter,:) =  X0.displacementFixed(:,fixedorfree);
        %                 fixCounter=fixCounter+1;
        nodecounter=nodecounter+1;
    else
        A=dirchnod>nn;
        idx=A==0;
        out=sum(idx(:));
        displacementsUseFull(nodecounter,:) =  X0.displacement(:,nn-out);
        nodecounter=nodecounter+1;
    end
    
    
end

% for i=1:nnode
%     X0RightCauchyGreenUseFull(1:d,i)=X0.RightCauchyGreen{i}(1:d,1);
%     X0RightCauchyGreenUseFull(d+1:2*d,i)=X0.RightCauchyGreen{i}(1:d,2);
% end

% figureReturn= PlotMesh(coordinates,nodesMat);
% set(gcf,'Position',[100 100 1000 1000])
% title('Mesh with RightCauchy Green Columns') ;
% set(get(gca,'title'),'Position',[10 20.5 1.00011])
% %     h = legend('circle', 'Plus', 'Location', 'NorthEast');
% % set(h, 'FontSize', 14)
% if d==2
%     %     axis([0 sizedom1+1 -1 sizedom2+1])
% else
%     axis([-3 sizedom1+3 -3 sizedom2+3 -1 1])
% end
% % axlim = get(gca, 'XLim');                                           % Get ‘XLim’ Vector
% % aylim = get(gca, 'YLim');                                           % Get ‘YLim’ Vector
% % legpos = get(h, 'Position');                                        % Get ‘legend’ 'Position' Vector
% % x_txt = 0.3;      % Set ‘x’-Coordinate Of Text Object
% % y_txt = (sizedom2)*1.1;      % Set ‘y’-Coordinate Of Text Object
% % text(x_txt, y_txt, resultstring(jj));
% 
% hold on;
% 
% if d==3
%     view(60,30);
% else
%     axis manual
% end
% if d==2
%     quiver(coordinates(:,1),coordinates(:,2),X0RightCauchyGreenUseFull(1,:)',X0RightCauchyGreenUseFull(2,:)',0,'LineWidth',2,'color',[0 0 1],'MaxHeadSize',1);
%     quiver(coordinates(:,1),coordinates(:,2),X0RightCauchyGreenUseFull(3,:)',X0RightCauchyGreenUseFull(4,:)',0,'LineWidth',2,'color',[0 0 1],'MaxHeadSize',1);
% else
%     %     quiver3(coordinates(:,1),coordinates(:,2),zeros(nnode,1),Xusefull(:,1),Xusefull(:,2),Xusefull(:,3),'LineWidth',2,'color',[0 0 1],'MaxHeadSize',2);
% end
% zoom(0.9)


% Create the manifold structure
problem.M =manifold;



% cost description
problem.cost = @cost;


    function [f, store] = cost(X, store)
        if ~isfield(store, 'ready')
            %             XXt = X*X';
            %             store.Xusefull=struct;
            
            [f,store]=getEnergy(X,coordinates, nodesMat,nel,nnel,d,dirchnod,mu,store);
            store.f = f;
            %             store.Xusefull = Xusefull;
            store.ready = true;
            
            % Shift the exponentials by the maximum value to reduce
            % numerical trouble due to possible overflows.
            %             s = max(max(triu(XXt, 1)));
            %             expXXt = exp((XXt-s)/epsilon);
            %             % Zero out the diagonal
            %             expXXt(1:(n+1):end) = 0;
            %             u = sum(sum(triu(expXXt, 1)));
            %             store.XXt = XXt;
            %             store.s = s;
            %             store.expXXt = expXXt;
            
            
            %         u = store.u;
            %         s = store.s;
            f = store.f;
            
        end
    end


% % gradient description
% problem.grad = @(X) problem.M.egrad2rgrad(X, egrad(X));
%     function g = egrad(X)
%         g = (C*C*X);
%     end


% Hessian description
%     problem.hess = @(X, U) problem.M.ehess2rhess(X, egrad(X), ehess(X, U), U);
%     function Hess = ehess(X, eta)
%         Hess = C*C*eta;
%     end

% Initialization

% Check numerically whether gradient and Ressian are correct
%             X0 = [1 0; 0 2]
%            d0 =  [  -0.128326596995917  -0.242841769548604;  -0.242841769548604   0.256653193991834];
%             trace(X0\d0)
figure
checkgradient(problem);
drawnow;
%     pause;
% checkhessian(problem);
%     drawnow;
%     pause;

% %         X0 = [1 1; 1 2]
% X0=problem.M.rand()


% Options (not mandatory)
options.maxiter = 10;
options.maxinner = 30;
options.maxtime = 120;
options.tolgradnorm = 1e-12;

% Pick an algorithm to solve the problem
% [Xopt, costopt, info] = trustregions(problem, X0, options);
        [Xopt costopt info] = conjugategradient(problem, X0, options);
%         [Xopt costopt info] = steepestdescent(problem, X0, options);
% figureReturn= PlotMesh(coordinates,nodesMat);
% set(gcf,'Position',[100 100 1000 1000])
% title('Mesh with RightCauchy Green Columns') ;
% set(get(gca,'title'),'Position',[10 20.5 1.00011])
% %     h = legend('circle', 'Plus', 'Location', 'NorthEast');
% % set(h, 'FontSize', 14)
% if d==2
%     %     axis([0 sizedom1+1 -1 sizedom2+1])
% else
%     axis([-3 sizedom1+3 -3 sizedom2+3 -1 1])
% end
% % axlim = get(gca, 'XLim');                                           % Get ‘XLim’ Vector
% % aylim = get(gca, 'YLim');                                           % Get ‘YLim’ Vector
% % legpos = get(h, 'Position');                                        % Get ‘legend’ 'Position' Vector
% % x_txt = 0.3;      % Set ‘x’-Coordinate Of Text Object
% % y_txt = (sizedom2)*1.1;      % Set ‘y’-Coordinate Of Text Object
% % text(x_txt, y_txt, resultstring(jj));
% 
% hold on;
% 
% if d==3
%     view(60,30);
% else
%     axis manual
% end
% if d==2
%     quiver(coordinates(:,1),coordinates(:,2),X0RightCauchyGreenUseFull(1,:)',X0RightCauchyGreenUseFull(2,:)',0,'LineWidth',2,'color',[0 0 1],'MaxHeadSize',1);
%     quiver(coordinates(:,1),coordinates(:,2),X0RightCauchyGreenUseFull(3,:)',X0RightCauchyGreenUseFull(4,:)',0,'LineWidth',2,'color',[0 0 1],'MaxHeadSize',1);
% else
%     %     quiver3(coordinates(:,1),coordinates(:,2),zeros(nnode,1),Xusefull(:,1),Xusefull(:,2),Xusefull(:,3),'LineWidth',2,'color',[0 0 1],'MaxHeadSize',2);
% end
% zoom(0.9)
end


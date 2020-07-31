function Test_symposspeciallinear()
% function Test_positivedefinite()
% Test for sympositivedefinite geometry (matrix completion)
%

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, August 29, 2013.
% Contributors:
% Change log:

    % Problem
    n = 10;
    B = randn(n, n)/10;
%     B = eye(n,n);
    B(1,1)=B(1,1)+1;
    B(2,1)=B(2,1)+1;
    C = B'*B
    
    
    % Create the manifold structure
    problem.M = symposspeciallinear(n);
    problem.M.transp = problem.M.paralleltransp;
    problem.M.retr = problem.M.exp; %man kann nicht expm benutzen für
%     compelx step!
    % creating a symmetric matrix of ones and zeros
%     f = 0.1; % % fraction of ones
%     length_of_vector = n*(n-1)/2;
%     test =(rand(length_of_vector,1) > f); % a vector with fraction f of zeros
%     train = logical(ones(length_of_vector,1) - test); % a vector with fraction f of ones
%     mask = sparse(squareform(train));
    
        
    % cost description
    problem.cost = @cost;
    function f = cost(X)
        f = .5*norm(C*X, 'fro')^2;
    end
    
    
    % gradient description
    problem.grad = @(X) problem.M.egrad2rgrad(X, egrad(X));
    function g = egrad(X)
        g = (C*C*X);
    end
    
    
%     % Hessian description
    problem.hess = @(X, U) problem.M.ehess2rhess(X, egrad(X), ehess(X, U), U);
    function Hess = ehess(X, eta)
        Hess = C*C*eta;
    end
    
      % Initialization

    % Check numerically whether gradient and Ressian are correct
%             X0 = [1 0; 0 2]
%            d0 =  [  -0.128326596995917  -0.242841769548604;  -0.242841769548604   0.256653193991834];
%             trace(X0\d0) 
%     checkgradient(problem);
%     drawnow;
%     pause;
    checkhessian(problem);
%     drawnow;
%     pause;
    
% %         X0 = [1 1; 1 2]
X0=eye(n);
det(X0)
eta0=problem.M.randvec(X0);
% trace(X0\eta0)

    det(X0)
  
    % Options (not mandatory)
    options.maxiter = 20;
    options.maxinner = 1000;
    options.maxtime = 120;
    options.tolgradnorm = 1e-10;
    
    % Pick an algorithm to solve the problem
    [Xopt, costopt, info] = trustregions(problem, X0, options);
%         [Xopt costopt info] = conjugategradient(problem, X0, options);
%         [Xopt costopt info] = steepestdescent(problem, X0, options);
    Xopt
    detC =det(Xopt)
    [eigCvec,eigCval] = eig(Xopt)
    C
end


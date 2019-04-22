function [L, S, LL, SS] = RobustPCA(X, lambda, mu, tol, max_iter)
    % - X is a data matrix (of the size N x M) to be decomposed
    %   X can also contain NaN's for unobserved values
    % - lambda - regularization parameter, default = 1/sqrt(max(N,M))
    % - mu - the augmented lagrangian parameter, default = 10*lambda
    % - tol - reconstruction error tolerance, default = 1e-6
    % - max_iter - maximum number of iterations, default = 1000

    [M, N] = size(X);
    unobserved = isnan(X);
    X(unobserved) = 0;
    normX = norm(X, 'fro');

    % default arguments
    if nargin < 2
        lambda = 1 / sqrt(max(M,N));
    end
    if nargin < 3
        mu = 10*lambda;
    end
    if nargin < 4
        tol = 1e-10;
    end
    if nargin < 5
        max_iter = 1000;
    end
    
    % initial solution
    L = zeros(M, N);
    S = zeros(M, N);
    Y = zeros(M, N);
    
    myiter = 1;
    for iter = (1:max_iter)
        % ADMM step: update L and S
        L = Do(1/mu, X - S + (1/mu)*Y);
        S = So(lambda/mu, X - L + (1/mu)*Y);
        % and augmented lagrangian multiplier
        Z = X - L - S;
        Z(unobserved) = 0; % skip missing values
        Y = Y + mu*Z;
        
%  				for im = 1:size(S,1)
%                     for j = 1:size(S,2)/3
%                    SA = sqrt((S(im,3*(j-1)+1))^2+(S(im,3*(j-1)+2))^2+(S(im,3*(j-1)+3))^2)+0.0001;
% %                    LA = sqrt((L(im,1))^2+(L(im,2))^2+(L(im,3))^2); 
% % 					if  SA < LA
% % 						S(im,1) = 10;
% % 						S(im,2) = 10;
% % 						S(im,3) = 10;
% % 					else
%  						S(im,3*(j-1)+1) = (1/SA)*S(im,3*(j-1)+1);
%  						S(im,3*(j-1)+2) = (1/SA)*S(im,3*(j-1)+2);
%  						S(im,3*(j-1)+3) = (1/SA)*S(im,3*(j-1)+3);
% % 					end
%                     end
%  				end
				
        err = norm(Z, 'fro') / normX;
         
        if (iter == 1) || (mod(iter, 10) == 0) || (err < tol)
            fprintf(1, 'iter: %04d\terr: %f\trank(L): %d\tcard(S): %d\n', ...
                    iter, err, rank(L), nnz(S(~unobserved)));
         LL{myiter} = L;
         SS{myiter} = S;
         myiter = myiter + 1;
        end
        if (err < tol) break; end
    end
end

function r = So(tau, X)
    % shrinkage operator
    r = sign(X) .* max(abs(X) - tau, 0);
end

function r = Do(tau, X)
    % shrinkage operator for singular values
    [U, S, V] = svd(X, 'econ');
    r = U*So(tau, S)*V';
end

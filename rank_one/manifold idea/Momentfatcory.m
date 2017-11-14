function M = Momentfatcory(n,A)
% Manifold of n-by-n symmetric positive semidefinite matrices X of rank k
% satisfying trace(A_i*X) = b_i for all i.
%
% function M = SDPfactory(n, k, A, b)
%
% A point X on the manifold is parameterized as YY^T where Y is a matrix of
% size nxk. As such, X is symmetric, positive semidefinite. We restrict to
% full-rank Y's, such that X has rank exactly k. The point X is numerically
% represented by Y (this is more efficient than working with X, which may
% be big). Tangent vectors are represented as matrices of the same size as
% Y, call them Ydot, so that Xdot = Y Ydot' + Ydot Y and trace(Xdot) == 0.
% The metric is the canonical Euclidean metric on Y.
%
%
% Note that this geometry formally breaks down at rank-deficient Y's.
% As an alternative, you may use the sphere manifold (it has larger
% dimension (by 1), but does not break down at rank drop.)
%
% The geometry is taken from the 2010 paper:
% M. Journee, P.-A. Absil, F. Bach and R. Sepulchre,
% "Low-Rank Optimization on the Cone of Positive Semidefinite Matrices".
% Paper link: http://www.di.ens.fr/~fbach/journee2010_sdp.pdf
% 
% 
% Please cite the Manopt paper as well as the research paper:
%     @Article{journee2010low,
%       Title   = {Low-rank optimization on the cone of positive semidefinite matrices},
%       Author  = {Journ{\'e}e, M. and Bach, F. and Absil, P.-A. and Sepulchre, R.},
%       Journal = {SIAM Journal on Optimization},
%       Year    = {2010},
%       Number  = {5},
%       Pages   = {2327--2351},
%       Volume  = {20},
%       Doi     = {10.1137/080731359}
%     }
% 
%
% See also: spherefactory elliptopefactory  spectahedronfactory 
%           symfixedrankYYfactory

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, July 11, 2013.
% Contributors: Nicolas Boumal
% Change log:
%
%   April 30, 2017 (NB):
%       Adapted from spectahedronfactory by JL
    
    % Number of constraints
    nCon = length(A);
    
    M.name = @() sprintf('Moment manifold of %dx1 vector with moment constraints',n);
    
    M.dim = @() n + 1; 
    
    % Euclidean metric on the total space
    M.inner = @(Y, eta, zeta) eta(:)'*zeta(:);
    
    M.norm = @(Y, eta) sqrt(M.inner(Y, eta, eta));
    
    M.dist = @(Y, Z) error('spectrahedronfactory.dist not implemented yet.');
    
    M.typicaldist = @() 1;
    
    M.proj = @projection;
    
    % construct AA
    Atemp = [];
    for ii = 1:nCon,
        Atemp = [Atemp;A{ii}];
    end
    AA = Atemp * Atemp';
    
    % Define projection function
    function etaproj = projection(Y, eta)
        % Projection onto the tangent space, i.e., on the tangent space of
        % trace(A_i*Y) = b_i forall i
        
        gamma = eta;
        gamma(1) = 0;
        
        % construct Ydiag
        Ydiag = [];
        for i = 1:nCon,
            Ydiag = blkdiag(Ydiag,Y);
        end
        
        BB = Ydiag'*Atemp*gamma;
        temp = Ydiag'*AA*Ydiag;
%         [rank(temp),size(temp)]   % check if temp is full rank
        alpha = temp\BB;
        etaproj = gamma - Atemp' * Ydiag * alpha;
    end
    
    M.tangent = M.proj;
    M.tangent2ambient = @(Y, eta) eta;
    
%     %%%%%%%%%%%% compute symbolic retraction map %%%%%%%
%     yy = [1;sym('y',[n-1,1],'real')];
%     zz = [0;sym('z',[n-1,1],'real')];
%     alpha1 = sym('a',[nCon,1],'real');
%     ww = yy + zz;
%     eqn = [];
%     for m = 1:nCon,
%         temp = ww'*A{m}*ww;
%         for i = 1:nCon,
%             temp = temp - alpha1(i)*yy'*A{i}*A{m}*ww - alpha1(i)*ww'*A{m}*A{i}'*yy;
%             for j = 1:nCon,
%                 temp = temp + alpha1(i)*alpha1(j)*yy'*A{i}*A{m}*A{j}'*yy;
%             end
%         end
%         eqn = [eqn;temp];
%     end
%     sol = solve(eqn,alpha1);
%     ALPHA = [sol.a1';sol.a2';sol.a3'];  % will change this hard-coding later
%     alpha1 = matlabFunction(ALPHA,'Vars',[yy;zz]);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    M.retr = @retractionPROJ;
    function Ynew = retractionPROJ(Y, eta, t)
        if nargin < 3
            t = 1.0;
        end
        Ybar = Y + t*eta;
        w = Ybar;
        w(1) = 1;
        alpha = sym('a',[nCon,1],'real');
        eqn = [];
        for m = 1:nCon,
            temp = w'*A{m}*w;
            for i = 1:nCon,
                temp = temp - alpha(i)*Y'*A{i}*A{m}*w - alpha(i)*w'*A{m}*A{i}'*Y;
                for j = 1:nCon,
                    temp = temp + alpha(i)*alpha(j)*Y'*A{i}*A{m}*A{j}'*Y;
                end
            end
            eqn = [eqn;temp];
        end
        
        sol = vpasolve(eqn,alpha);
        alpha = double([sol.a1';sol.a2';sol.a3']);
        if size(alpha,2)==2,
            res1 = w;
            res2 = w;
            for i = 1:nCon,
                res1 = res1 - A{i}'*Y*alpha(i,1);
                res2 = res2 - A{i}'*Y*alpha(i,2);
            end

            temp = [norm(w-res1),norm(w-res2)];
            if find(temp==min(temp))==1,
                Ynew = res1;
            else
                Ynew = res2;
            end
        else
            Ynew = w;
            for i = 1:nCon,
                Ynew = Ynew - A{i}'*Y*alpha(i,1);
            end
        end
        [Ynew';Y']
    end
    
    M.egrad2rgrad = @projection;
    
    M.ehess2rhess = @ehess2rhess;  % DON'T USE HESSIAN: NOT IMPLEMENTED YET!!!!!!!!!!!!!!!!!!!
    function Hess = ehess2rhess(Y, egrad, ehess, eta)
        YtEta = Y.'*eta;
        
        A_ly_Om = Y.'*Y;
        Q_ly_Om = -(YtEta - YtEta.');
        Omega = lyap(A_ly_Om, Q_ly_Om);
        
        A_ly_DOm = A_ly_Om;
        Q_ly_DOm = -(comm(eta.',egrad) + comm(Y.',ehess) + comm(Omega.',YtEta + YtEta.'));
        DOmega = lyap(A_ly_DOm, Q_ly_DOm);
        
        Hess = ehess  - eta*Omega - Y*DOmega;
        
        for i = 1:nCon
            AY = A{i}*Y;
            Aeta = A{i}*eta;
            trYA2Y = (AY(:).'*AY(:));
            
            alpha = eta(:).'*AY(:)/trYA2Y;
            Dalpha = 1/trYA2Y*(ehess(:).'*AY(:)+egrad(:).'*Aeta(:)) - ...
                     eta(:).'*AY(:)/trYA2Y^2 * (2*Aeta(:).'*AY(:));
                 
            Hess = Hess - alpha*Aeta - Dalpha*AY;
        end
        
        % Project on the horizontal space
        Hess = M.proj(Y, Hess);
    end
    
    M.exp = @exponential;
    function Ynew = exponential(Y, eta, t)
        if nargin < 3
            t = 1.0;
        end
        
        Ynew = retractionPROJ(Y, eta, t);
        warning('manopt:spectrahedronfactory:exp', ...
            ['Exponential for Moment ' ...
            'manifold not implenented yet. Used retraction instead.']);
    end
    
    % Notice that the hash of two equivalent points will be different...
    M.hash = @(Y) ['z' hashmd5(Y(:))];
    
    M.rand = @random;
    
    function Y = random()
        y = randn(1);
        Y = [1;y];
        for i = 2:n-1
            Y = [Y;y^i];
        end
    end
    
    M.randvec = @randomvec;
    function eta = randomvec(Y)
        eta = randn(n, 1);
        eta = projection(Y, eta);
        nrm = M.norm(Y, eta);
        eta = eta / nrm;
    end
    
    M.lincomb = @matrixlincomb;
    
    M.zerovec = @(Y) zeros(n, 1);
    
    M.transp = @(Y1, Y2, d) projection(Y2, d);
    
    M.vec = @(Y, u_mat) u_mat(:);
    M.mat = @(Y, u_vec) reshape(u_vec, [n, 1]);
    M.vecmatareisometries = @() true;
end

function res = comm(x,y)
    xy = x*y;
    res = xy - xy.';
end
function res = anticomm(x,y)
    xy = x*y;
    res = xy + xy.';
end
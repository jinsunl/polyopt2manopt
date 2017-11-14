function hess = getHessian_aug(problem, x, d, storedb, key)
% Computes the Hessian of the cost function at x along d.
%
% function hess = getHessian(problem, x, d)
% function hess = getHessian(problem, x, d, storedb)
% function hess = getHessian(problem, x, d, storedb, key)
%
% Returns the Hessian at x along d of the cost function described in the
% problem structure.
%
% storedb is a StoreDB object, key is the StoreDB key to point x.
%
% If an exact Hessian is not provided, an approximate Hessian is returned
% if possible, without warning. If not possible, an exception will be
% thrown. To check whether an exact Hessian is available or not (typically
% to issue a warning if not), use canGetHessian.
%
% See also: getPrecon getApproxHessian canGetHessian

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
%
%   April 3, 2015 (NB):
%       Works with the new StoreDB class system.

    % Allow omission of the key, and even of storedb.
    if ~exist('key', 'var')
        if ~exist('storedb', 'var')
            storedb = StoreDB();
        end
        key = storedb.getNewKey();
    end
    
    
    hess = problem.rhess_aug(x,d);
    
    
end

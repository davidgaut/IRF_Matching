function y = irf_comp(dr, e1, long, drop, replic, iorder,M_,options_)

% function y = irf(dr, e1, long, drop, replic, iorder)
% Computes impulse response functions
%
% INPUTS
%    dr:     structure of decisions rules for stochastic simulations
%    e1:     exogenous variables value in time 1 after one shock
%    long:   number of periods of simulation
%    drop:   truncation (in order 2)
%    replic: number of replications (in order 2)
%    iorder: first or second order approximation
%
% OUTPUTS
%    y:      impulse response matrix
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.



if M_.maximum_lag >= 1
    temps = repmat(dr.ys,1,M_.maximum_lag);
else
    temps = zeros(M_.endo_nbr, 1); % Dummy values for purely forward models
end
y       = 0;

local_order = iorder;
if M_.hessian_eq_zero && local_order~=1
    local_order = 1;
end

if local_order == 1
    y1 = repmat(dr.ys,1,long);
    ex2 = zeros(long,M_.exo_nbr);
    ex2(1,:) = e1';
    y2 = simult_comp(temps,dr,ex2,local_order,M_,options_);
    y = y2(:,M_.maximum_lag+1:end)-y1;
else
    % eliminate shocks with 0 variance
    i_exo_var = setdiff([1:M_.exo_nbr],find(diag(M_.Sigma_e) == 0 ));
    nxs = length(i_exo_var);
    ex1 = zeros(long+drop,M_.exo_nbr);
    ex2 = ex1;
    chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));
    for j = 1: replic
        ex1(:,i_exo_var) = randn(long+drop,nxs)*chol_S;
        ex2 = ex1;
        ex2(drop+1,:) = ex2(drop+1,:)+e1';
        y1 = simult_(temps,dr,ex1,local_order);
        y2 = simult_(temps,dr,ex2,local_order);
        y = y+(y2(:,M_.maximum_lag+drop+1:end)-y1(:,M_.maximum_lag+drop+1:end));
    end
    y=y/replic;
end


function [fval,dlik] = irf_matching(Model,oo_,options_,lnprior)

Var_IRF        = Model.Var_IRF;
invVarNorm     = Model.invVarNorm;
logdetVarNorm  = Model.logdetVarNorm;
horizon_est    = Model.horizon_est;
mod_var_list   = Model.mod_var_list;
mod_shock_list = Model.mod_shock_list;
SR             = Model.SR;

% Get impulse responses
warning off;
var_list     = mod_var_list;
[i_var,~]    = varlist_indices(var_list,Model.endo_names);

iter_ = max(options_.periods,1);
if Model.exo_nbr > 0
    oo_.exo_simul= ones(iter_ + Model.maximum_lag + Model.maximum_lead,1) * oo_.exo_steady_state';
end

SS(Model.exo_names_orig_ord,Model.exo_names_orig_ord) = Model.Sigma_e+1e-14*eye(Model.exo_nbr);
cs = transpose(chol(SS));

[i_var_exo,nvar_exo] = varlist_indices(mod_shock_list,Model.exo_names);

% Compute Impulse Response
y2 = zeros(Model.endo_nbr,horizon_est,nvar_exo);
for jj=1:nvar_exo        
    y2(:,:,jj) = irf_comp(oo_.dr,cs(Model.exo_names_orig_ord,i_var_exo(jj)), horizon_est, options_.drop, ...
                    options_.replic, options_.order,Model,options_) * 1;
end

% Stack model impulse responses
DSGE_IRF = permute(y2(i_var,:,:),[2,1,3]);
DSGE_IRF = DSGE_IRF(1:horizon_est,:,:) * 1;
psitheta = DSGE_IRF(:);

% evaluate criterion and form likelihood
criterion  = ((Var_IRF-psitheta))'*invVarNorm*((Var_IRF-psitheta));
likelihood = (size(Var_IRF,1)/2*log(1/2/3.141592653589793)-1/2*logdetVarNorm-1/2*criterion); %in logs
dlik       = (Var_IRF-psitheta).^2.*diag(invVarNorm);

if isnan(criterion); 
            fval = Inf;
return
end

% ------------------------------------------------------------------------------
% Form posterior for impulse response matching
% ------------------------------------------------------------------------------
fval = - (likelihood + lnprior);




clc, close all

// Create fake data series to trick Dynare
X1 = [1:10];
save fsdat X1

// Add codes
addpath(genpath('Codes'))

// Model
%--------------------------------------------------------------------------
var X1 X2;
varexo eps1 eps2;
parameters f11 f12 f21 f22 b12 b21;

// Calibrated params
f11 = 0.4;
f12 = 0.4;
f21 = 0.4;
f22 = 0.4;
b12 = 1.;
b21 = 0.6;

// Model is a VAR(1) in [X1 X2]
model;
    X1 = f11*X1(-1) + f12*X2(-1) + eps1         + (1/b12)*eps2;
    X2 = f21*X1(-1) + f22*X2(-1) + (1/b21)*eps1 +         eps2;
end;

steady_state_model;
    X2 = 0;
    X1 = 0;
end;

shocks;
    var eps1; stderr 1;
    var eps2; stderr 1;
end;

steady;
stoch_simul(order=1);

// Save calibrated params for later
calib = M_.params;

// Create fake VAR IRFs
%--------------------------------------------------------------------------
clear X1_resp_eps1 X1_resp_eps2 X2_resp_eps1 X2_resp_eps2 Var_IRF

// True VAR parameters
F11 = 0.3;
F12 = 0.2;
F21 = 0.5;
F22 = 0.1;
B12 = 1;
B21 = 0.5;
sd = 0.5;

// Construct IRFs
X1_resp_eps1(:,1) = 1+sd*randn(1e3,1);
X2_resp_eps1(:,1) = (1/B21)*X1_resp_eps1(:,1);
X2_resp_eps2(:,1) = 1+sd*randn(1e3,1);
X1_resp_eps2(:,1) = (1/B12)*X2_resp_eps2(:,1);
for horiz = 2 : 10
    X1_resp_eps1(:,horiz) = F11*X1_resp_eps1(:,horiz-1) + F12*X2_resp_eps1(:,horiz-1) ;
    X2_resp_eps1(:,horiz) = F21*X1_resp_eps1(:,horiz-1) + F22*X2_resp_eps1(:,horiz-1);

    X1_resp_eps2(:,horiz) = F11*X1_resp_eps2(:,horiz-1) + F12*X2_resp_eps2(:,horiz-1) ;
    X2_resp_eps2(:,horiz) = F21*X1_resp_eps2(:,horiz-1) + F22*X2_resp_eps2(:,horiz-1);
end

// IRF matching settings
%--------------------------------------------------------------------------

// Priors (center priors on true VAR parameters)
estimated_params;
    f11, normal_pdf, F11, 1;
    f12, normal_pdf, F12, 1;
    f21, normal_pdf, F21, 1;
    f22, normal_pdf, F22, 1;
    b12, normal_pdf, B12, 1;
    b21, normal_pdf, B21, 1;
end;
             
// Declare observed variable (only to trick Dynare)
varobs X1;

// Settings
mod_var_list = {'X2' 'X1'};                                                % List variables to match
mod_shock_list = {'eps1' 'eps2'};                                          % Shocks conditional on which to match
options_.irfs_match_estimation = 1;                                        % IRF Matching
horizon_est  = 10;                                                         % Horizon to match


// Plot VAR IRFs and compare them with calibrated DSGE
%--------------------------------------------------------------------------
FigSize
subplot(2,2,1)
H = PlotSwathe(mean(X1_resp_eps1),squeeze(prctile(X1_resp_eps1,[5 95]))); hold on; h = plot([oo_.irfs.X1_eps1(1:horizon_est)'],'--','LineWidth',1.5,'Color',cmap(2)); 
title('X1 to eps1'); legend([H.patch h],{'VAR';'Calibrated DSGE'})
subplot(2,2,2)
H = PlotSwathe(mean(X1_resp_eps2),squeeze(prctile(X1_resp_eps2,[5 95]))); hold on; h = plot([oo_.irfs.X1_eps2(1:horizon_est)'],'--','LineWidth',1.5,'Color',cmap(2)); 
title('X1 to eps2'); legend([H.patch h],{'VAR';'Calibrated DSGE'})
subplot(2,2,3)
H = PlotSwathe(mean(X2_resp_eps1),squeeze(prctile(X2_resp_eps1,[5 95]))); hold on; h = plot([oo_.irfs.X2_eps1(1:horizon_est)'],'--','LineWidth',1.5,'Color',cmap(2)); 
title('X2 to eps1'); legend([H.patch h],{'VAR';'Calibrated DSGE'})
subplot(2,2,4)
H = PlotSwathe(mean(X2_resp_eps2),squeeze(prctile(X2_resp_eps2,[5 95]))); hold on; h = plot([oo_.irfs.X2_eps2(1:horizon_est)'],'--','LineWidth',1.5,'Color',cmap(2)); 
title('X2 to eps2'); legend([H.patch h],{'VAR';'Calibrated DSGE'})


// Reshape IRFs for IRF matching 
%--------------------------------------------------------------------------
// Reshape
X1_resp_eps1 = permute(X1_resp_eps1,[2 , 3 , 4 ,1]);
X1_resp_eps2 = permute(X1_resp_eps2,[2 , 3 , 4 ,1]);
X2_resp_eps1 = permute(X2_resp_eps1,[2 , 3 , 4 ,1]);
X2_resp_eps2 = permute(X2_resp_eps2,[2 , 3 , 4 ,1]);

// Stack by (horizon * var * shock * obs)
Var_IRF(:,:,1,:) = [X2_resp_eps1 X1_resp_eps1];
Var_IRF(:,:,2,:) = [X2_resp_eps2 X1_resp_eps2];

// Settings
StdNorm          = std(Var_IRF,0,4);
Var_IRF          = mean(Var_IRF,4);
Var_IRF          = Var_IRF(:);
StdNorm          = StdNorm(:);
VarNorm          = StdNorm.^2;
VarNorm          = diag(VarNorm);
invVarNorm       = inv(VarNorm);
logdetVarNorm    = logdet(VarNorm);
M_.Var_IRF         = Var_IRF;
M_.invVarNorm      = invVarNorm;
M_.logdetVarNorm   = logdetVarNorm;
M_.horizon_est     = horizon_est;
M_.invVarNorm      = invVarNorm;
M_.mod_var_list    = mod_var_list;
M_.mod_shock_list  = mod_shock_list;
M_.SR              = [];                 

// Estimation
%--------------------------------------------------------------------------
estimation(datafile='fsdat.mat',mh_replic=0,mh_nblocks=2,mh_jscale=0.8,mode_compute = 7,plot_priors = 1);
xparam = get_posterior_parameters('mode',M_,estim_params_,oo_,options_);
M_ = set_all_parameters(xparam,estim_params_,M_);
stoch_simul(order=1);

// Plot VAR IRFs and compare them with estimated DSGE
%--------------------------------------------------------------------------
FigSize
subplot(2,2,1)
H = PlotSwathe(mean(X1_resp_eps1,4),squeeze(prctile(X1_resp_eps1,[5 95],4))); hold on; h = plot([oo_.irfs.X1_eps1(1:horizon_est)'],'--','LineWidth',1.5,'Color',cmap(2)); 
title('X1 to eps1'); legend([H.patch h],{'VAR';'Estimated DSGE'})
subplot(2,2,2)
H = PlotSwathe(mean(X1_resp_eps2,4),squeeze(prctile(X1_resp_eps2,[5 95],4))); hold on; h = plot([oo_.irfs.X1_eps2(1:horizon_est)'],'--','LineWidth',1.5,'Color',cmap(2)); 
title('X1 to eps2'); legend([H.patch h],{'VAR';'Estimated DSGE'})
subplot(2,2,3)
H = PlotSwathe(mean(X2_resp_eps1,4),squeeze(prctile(X2_resp_eps1,[5 95],4))); hold on; h = plot([oo_.irfs.X2_eps1(1:horizon_est)'],'--','LineWidth',1.5,'Color',cmap(2)); 
title('X2 to eps1'); legend([H.patch h],{'VAR';'Estimated DSGE'})
subplot(2,2,4)
H = PlotSwathe(mean(X2_resp_eps2,4),squeeze(prctile(X2_resp_eps2,[5 95],4))); hold on; h = plot([oo_.irfs.X2_eps2(1:horizon_est)'],'--','LineWidth',1.5,'Color',cmap(2)); 
title('X2 to eps2'); legend([H.patch h],{'VAR';'Estimated DSGE'})


// Print estimates and compare then with calibrated/priors
%--------------------------------------------------------------------------
aux = [F11; F12; F21; F22; B12; B21;];
aux = [calib aux xparam oo_.prior.mean ]
in.cnames = strvcat('calibrated', 'true','posterior','prior');
in.rnames = ['   '; M_.param_names];
disp('-----------')
disp('Estimation results:')
mprint(aux,in)



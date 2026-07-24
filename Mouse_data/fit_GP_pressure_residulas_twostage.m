Mou = [1,3,4,5];
addpath('Results\')
save_on = false;

for mm = 1:4
    mouse = num2str(Mou(mm));    
    data{mm} = importdata(['../Mouse_data/pH',mouse,'_1024.dat']);
end


for mm = 1:4
    mouse = num2str(Mou(mm));
%     for ii = 1:R
%         ResM{ii,mm} = load(['NumOpt_mousedata',mouse,'_5p_f1f2vaso2x_stnd_results',num2str(ii),'_relaxedbounds_1vaso.mat']);    
%     end

%     for ii = 1:20
%         ResMrel{ii,mm} = load(['NumOpt_mousedata',mouse,'_5p_f1f2vaso2x_stnd_results',num2str(ii),'_relaxedbounds.mat']);    
%     end

    for ii = 1:20
        ResMold{ii,mm} = load(['NumOpt_mousedata',mouse,'_5p_f1f2vaso2x_stnd_results',num2str(ii),'.mat']);    
    end
  
end

% Res = [ResMold;ResMrel];
Res = [ResMrel];

% fvalOld = cellfun(@(xx) xx.fval, ResMold);
fvalRel = cellfun(@(xx) xx.fval, ResMrel);
% fval = [fvalOld;fvalRel];
fval = [fvalRel];
[~,indmin] = min(fval);

pmin = cellfun(@(xx) xx.p(:,1),Res,'UniformOutput',false);
indmin_ind = indmin + size(Res,1) * ((1:size(Res,2)) - 1);

pmin_ind = pmin(indmin_ind);
pmin_ind = cell2mat(pmin_ind);
data = cell2mat(data);
data = data(1:2:end,:);

pressure_fit = pmin_ind;

save('Pressure_fit_Data.mat','pressure_fit',"data")
writematrix([pressure_fit,data],'Pressure_fit_Data.csv')


residulas = pmin_ind - data;

figure(100)
set(0,'defaulttextInterpreter','latex') %latex axis labels
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(gcf,'units','normalized','outerposition',[0.2 0.2 0.6 0.7]);

subplot(2,2,1)
plot(residulas,LineWidth=2)
xlim([1,512])
% legend([repelem('M',4)',num2str(Mou')])
legend({'Mouse 1','Mouse 3', 'Mouse 4', 'Mouse 5'})
title('(a) Pressure residuals (model fit - data)')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Step 1
% Get the mean f(t) and variance var(t) as functions of time t from
% the four pressure discrepancy profiles (measurements minus model fit)
% you have. Since all four profiles are smooth functions of t, 
% and the average of smooth functions is also smooth,
% these are smooth functions. 

res_mean = mean(residulas,2);
res_var = var(residulas,[],2);
res_std = sqrt(res_var);

figure(44)
plot(residulas,LineWidth=2)
xlim([1,512])
% legend([repelem('M',4)',num2str(Mou')])
title('Residuals (fit - data)')
hold on
plot(res_mean,'k',LineWidth=2)
plot(res_mean+sqrt(res_var),'k--',LineWidth=2)
plot(res_mean-sqrt(res_var),'k--',LineWidth=2)

figure(45)
plot(res_std,LineWidth=2)
xlim([1,512])

%% Step 2
% Fit a GP with squared exponential kernel to your pressure discrepancy
% profiles and obtain the lengthscale L.
X = (1:512)';
for mm = 1:4
    Y = residulas(:,mm);
    GPfit{1,mm} =  fitrgp(X,Y,'KernelFunction','ardsquaredexponential','standardize',1);
end

GPfit{1,1}.KernelInformation.KernelParameters
GPfit{1,1}.KernelInformation.KernelParameterNames


LenSc = cellfun(@(xx) xx.KernelInformation.KernelParameters(1,1), GPfit)';
MagnSig = cellfun(@(xx) xx.KernelInformation.KernelParameters(2,1), GPfit)';


LenScR = 0.2895; % From R package hetGP


%% Step 3
% Take the kernel defined in eq.(2) of the paper by Heinonen et al.
% where x in that equation corresponds to time t, l(t)=l(t')=L, 
% and sigma(t)= sqrt(var(t)), as obtained in Steps 1 and 2.

ker_heinonen = @(xx,xx2,ll,ll2) sqrt(2.*ll.*ll2./(ll.^2+ll2.^2)) .* exp(-(pdist2(xx,xx2).^2)./(ll.^2+ll2.^2));

if false
figure(56)
subplot(1,2,1)
hold on
plot(ker_heinonen(X,mean(X),1e4*LenScR,1e4*LenScR),LineWidth=2)
plot(ker_heinonen(X,X(1),1e4*LenScR,1e4*LenScR),LineWidth=2)
plot(ker_heinonen(X,X(end),1e4*LenScR,1e4*LenScR),LineWidth=2)
xlim([1,512])
legend({'Input corr with middle intput','Input corr with first input','Input corr with last input'},'FontSize',12)
title('Heinonen kernel modulo sigma')

subplot(1,2,2)
hold on
plot(ker_heinonen(res_std.*X,res_std(round((1+512)/2))*X(round((1+512)/2)),1e4*LenScR,1e4*LenScR),LineWidth=2)
plot(ker_heinonen(res_std.*X,res_std(1)*X(1),1e4*LenScR,1e4*LenScR),LineWidth=2)
plot(ker_heinonen(res_std.*X,res_std(end)*X(end),1e4*LenScR,1e4*LenScR),LineWidth=2)
xlim([1,512])
title('Heinonen kernel scaled by heterogenous sigma')
legend({'Input corr with middle intput','Input corr with first input','Input corr with last input'},'FontSize',12)
end


ker_heinonen_scaled = @(xx,xx2,ll,ll2) res_std(xx).*res_std(xx2).*sqrt(2.*ll.*ll2./(ll.^2+ll2.^2)) .* exp(-(pdist2(xx,xx2).^2)./(ll.^2+ll2.^2));
ker_heinonen_scaled(X,X(1),1e4*LenScR,1e4*LenScR)

%% Step 4
% Sample noise realisations by sampling signals from 
% a zero-mean GP with the covariance matrix from Step 3, 
% and then adding these signals to f(t).

% Kstar = ker_heinonen_scaled(X,X,1e4*LenScR,1e4*LenScR);
Kstar = ker_heinonen(X,X,1e4*LenScR,1e4*LenScR);
issymmetric(Kstar)
det(Kstar)

ScaleSigma = res_std*res_std';
issymmetric(ScaleSigma)
det(ScaleSigma)

Kstar = Kstar.*ScaleSigma;
issymmetric(Kstar)
det(Kstar)

Nsamples = 100;
CholKstar = chol(Kstar+ diag(1e-9*ones(512,1)));
draw = randn(512,Nsamples);
draw = res_mean + CholKstar'*draw;
    


figure(57)
plot(residulas,'Color',[0.7,0.7,0.7],LineWidth=2)
hold on
plot(draw(:,1:10),'Color',[0,0.7,0.9],LineWidth=1)
xlim([1,512])
title('Grey: residuals, Blue: samples')


figure(58)
hold on
plot(draw,LineWidth=1)
plot(residulas,'k--',LineWidth=2)
xlim([1,512])
title('Black dashed: residuals, Color: samples')

if save_on
    save('Results/Pressure_draws_GPprior.mat',"draw","res_mean","res_var",...
        "CholKstar","ker_heinonen","LenSc","LenScR","MagnSig")
end


figure(100)
subplot(2,2,2)
hold on
plot(draw(:,1:50),LineWidth=1)
plot(residulas,'k--',LineWidth=2)
xlim([1,512])
title('(b) Pressure residuals (black) and GP samples (colour)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Use the fitted GP model for the pressure to obtain a GP model for the flows

q2min = cellfun(@(xx) xx.q(:,2),Res,'UniformOutput',false);
q3min = cellfun(@(xx) xx.q(:,3),Res,'UniformOutput',false);
% indmin_ind = indmin + size(Res,1) * ((1:size(Res,2)) - 1);

q2min_ind = q2min(indmin_ind);
q2min_ind = cell2mat(q2min_ind);
q3min_ind = q3min(indmin_ind);
q3min_ind = cell2mat(q3min_ind);

flow2_fit = q2min_ind;
flow3_fit = q3min_ind;

figure(88)
set(gcf,'units','normalized','outerposition',[0.30 0.36 0.57 0.53]);
subplot(1,3,1)
plot(pressure_fit,LineWidth=2)
xlim([1,512])
title('Pressure fitted')
subplot(1,3,2)
plot(flow2_fit,LineWidth=2)
xlim([1,512])
title('Flow v2 fitted')
subplot(1,3,3)
plot(flow3_fit,LineWidth=2)
xlim([1,512])
title('Flow v3 fitted')
if save_on
    name = ['Results/pressure_flow_fit.png'];
    set(gcf,'PaperPositionMode','auto');
    print(gcf,name,'-dpng','-r0')    
end

var(pressure_fit)
var(flow2_fit)
var(flow3_fit)
 
scale_flow2 = std(flow2_fit,[],2)./std(pressure_fit,[],2);
scale_flow3 = std(flow3_fit,[],2)./std(pressure_fit,[],2);


flow2_mean = mean(flow2_fit,2);
flow3_mean = mean(flow3_fit,2);


figure(61)
set(gcf,'units','normalized','outerposition',[0.30 0.36 0.57 0.53]);
subplot(1,3,1)
plot(std(pressure_fit,[],2),LineWidth=2)
xlim([1,512])
title('std(Pressure fitted)')
subplot(1,3,2)
plot(std(flow2_fit,[],2),LineWidth=2)
xlim([1,512])
title('std(Flow v2 fitted)')
subplot(1,3,3)
plot(std(flow3_fit,[],2),LineWidth=2)
xlim([1,512])
title('std(Flow v3 fitted)')
if save_on
    name = ['Results/std_pressure_flow_fit.png'];
    set(gcf,'PaperPositionMode','auto');
    print(gcf,name,'-dpng','-r0')    
end

figure(62)
plot(scale_flow2,LineWidth=2)
hold on 
plot(scale_flow3,LineWidth=2)
legend('flow 2','flow 3')
title('Scaling for GP sigma for flow fits')
xlim([1,512])
if save_on
    name = ['Results/gp_time_var_scaling.png'];
    set(gcf,'PaperPositionMode','auto');
    print(gcf,name,'-dpng','-r0')  
end



Kstar0 = ker_heinonen(X,X,1e4*LenScR,1e4*LenScR);

ScaleSigma_flow2 = (scale_flow2.*res_std)*(scale_flow2'.*res_std');
ScaleSigma_flow3 = (scale_flow3.*res_std)*(scale_flow3'.*res_std');
issymmetric(ScaleSigma_flow2)
det(ScaleSigma_flow2)

Kstar_flow2 = Kstar0.*ScaleSigma_flow2;
Kstar_flow3 = Kstar0.*ScaleSigma_flow3;
issymmetric(Kstar_flow2)
issymmetric(Kstar_flow3)
det(Kstar)


Nsamples = 100;
CholKstar_flow2 = chol(Kstar_flow2 + diag(1e-8*ones(512,1)));
CholKstar_flow3 = chol(Kstar_flow3 + diag(1e-9*ones(512,1)));
draw_flow2 = randn(512,Nsamples);
draw_flow2 = CholKstar_flow2'*draw_flow2;
    
draw_flow3 = randn(512,Nsamples);
draw_flow3 = CholKstar_flow3'*draw_flow3;
    

figure(63)
set(gcf,'units','normalized','outerposition',[0.30 0.36 0.57 0.53]);
subplot(1,2,1)
hold on
plot(draw_flow2,LineWidth=1)
% plot(residulas,'k--',LineWidth=2)
title('Zero-mean noise for flow v2')
xlim([1,512])
subplot(1,2,2)
hold on
plot(draw_flow3,LineWidth=1)
% plot(residulas,'k--',LineWidth=2)
title('Zero-mean noise for flow v3')
xlim([1,512])
% title('Black dashed: residuals, Color: samples')
if save_on
    name = ['Results/gp_noise_samples_flow.png'];
    set(gcf,'PaperPositionMode','auto');
    print(gcf,name,'-dpng','-r0')  
end


figure(64)
set(gcf,'units','normalized','outerposition',[0.30 0.36 0.57 0.53]);
subplot(1,2,1)
hold on
plot(flow2_mean + draw_flow2,LineWidth=1)
plot(flow2_fit,'k--',LineWidth=2)
title('Flow v2')
xlim([1,512])
subplot(1,2,2)
hold on
plot(flow3_mean + draw_flow3,LineWidth=1)
plot(flow3_fit,'k--',LineWidth=2)
title('Flow v3')
xlim([1,512])
sgtitle('Black dashed: residuals, Color: samples')
    name = ['Results/flow_gp_samples.png'];
    set(gcf,'PaperPositionMode','auto');
    print(gcf,name,'-dpng','-r0')  


figure(100)
subplot(2,2,3)
hold on
% plot(mean(flow2_fit,2)+ draw_flow2_scalar,LineWidth=1)
plot(flow2_mean + draw_flow2(:,1:50),LineWidth=1)
plot(flow2_fit,'k--',LineWidth=2)
title('(c) Fitted LPA flow (black) and GP samples (colour)')
xlim([1,512])
subplot(2,2,4)
hold on
plot(flow3_mean + draw_flow3(:,1:50),LineWidth=1)
plot(flow3_fit,'k--',LineWidth=2)
title('(d) Fitted RPA flow (black) and GP samples (colour)')
xlim([1,512])
if save_on
    name = ['Results/pressure_residulas_GPsamples_pressure_flow.pdf'];
    set(gcf,'PaperPositionMode','auto');
%     print(gcf,name,'-dpng','-r0')
    exportgraphics(gcf,name, 'ContentType', 'vector');
end


%% SAVE FOR SAMPLING
if save_on
    data_pressure = data;
    save('Results/Pressure_flow_sampling.mat',...
        'data_pressure','CholKstar','CholKstar_flow2','CholKstar_flow3','res_mean','flow2_mean','flow3_mean');
end
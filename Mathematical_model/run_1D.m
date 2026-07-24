% Run the 1D ST code in the 21 vessel mouse geometry
% clear; close all;
% !chmod +x sor06
!make clean
!make


% Exponential stiffness
f1   = 7e+4;%5e+6;
f2   = -10; %10
f3   = 1e4;%2.5e4;%1e+4;%8e+4;
fs1  = f1;%5e+6;%f1;
fs2  = f2;%-20;%f2;
fs3  = f3;%*10;%1e+6;%f3;
Z0  = 0;%1e2;


alpha = 0.88; %Alpha
beta  = 0.68; %Beta
rm    = 0.005;%minimum radius
lrr   = 17;   %length-to-radius ratio
%%

vaso_L = 1.0;
vaso_R = 0.6;

pars = [f1 f2 f3 fs1 fs2 fs3 alpha beta lrr rm Z0 vaso_L vaso_R];
pars_str = mat2str(pars);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
out = unix(sprintf('sor06.exe  %s',pars_str(2:end-1)));
toc
%%

if out == 0
    fname = strcat('pu_ALL.2d');
    data = load(fname);
    [t,x,p,q,a,c] = gnuplot(data);
end

figure; plot(p);
figure; plot(q);

% look at flow split

figure(99);
subplot(1,2,1); hold on; plot(q(:,2),'LineWidth',3);
subplot(1,2,2); hold on; plot(q(:,3),'LineWidth',3);



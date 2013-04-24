function gmres_diffBC
%let the residual tol very small, while change MG sweeps large, say 2--10. then compare the convergence behavior
%check the log(res) Vs itations
%KSP 
% fgmres

% 3D case
% the below is 3D case: 64*64*64, true norm is used
%Here we use the real residual norm, not KSP residual
% All Dirichlet
%resid_1=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_1VcylePerIt_finegrid.dat');

figure
% all Dirichlet, Vel MG using different bottom smooth steps
%resid_1=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_1VcylePerIt_finegrid.dat');  %  
resid_1=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_gmres_finegridGS_16botsmooth.dat');  %  
%
% % all Dirichlet, Pres MG
%resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_slip_1VcylePerIt_finegrid.dat'); %
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_slip_gmres_finegridGS_16botsmooth.dat'); %
%resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_periodic_1VcylePerIt_finegrid.dat'); %
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_periodic_gmres_finegridGS_16botsmooth.dat'); %
 
%plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b*'); 
%hold on;
plot([0;resid_1(1:end,1)], [0;log10(resid_1(1:end,2))],  'b*'); 
hold on;
plot([0;resid_2(1:end,1)], [0;log10(resid_2(1:end,2))],  'r*'); 
hold on;
plot([0;resid_3(1:end,1)], [0;log10(resid_3(1:end,2))],  'k*');  
hold on;
%plot([0;resid_4(:,1)], [0;log10(resid_4(:,2))],  'ro'); 
%hold on; 
%plot([0;resid_5(:,1)], [0;log10(resid_5(:,2))],  'ks'); 
%hold on;
%plot([0;resid_6(:,1)], [0;log10(resid_6(:,2))],  'c^'); 
%hold on;
%plot([0;resid_7(1:end,1)], [0;log10(resid_7(1:end,2))],  'r^'); 
%hold on; 
%plot([0;resid_8(1:end,1)], [0;log10(resid_8(1:end,2))],  'k^'); 
%hold on; 


%hold on;
% legend('Diri','Diri Vel MG','Diri Pre MG',  'x-Peri', 'xy-Peri', 'Peri', 'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(10, -14.93, 'total no. of MG Sweeps');
% axis square 

% legend('Diri Vel MG','Diri Pre MG',  'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(4, -14.93, 'total no. of MG Sweeps');
legend('Diri','Slip', 'Per' ,'Location','NorthEast')
 text(10, -15.13, 'total no. of MG Sweeps');
axis square 
%legend('Pure Dirichlet 1 MG/it', 'x- Periodic 1 MG/it', 'x,y -Periodic 1 MG/it', 'Pure Periodic 1 MG/it', 'Pure Periodic 2 MG/it', 'Pure Periodic 4 MG/it', 'Pure Periodic 10 MG/it', 'Location','NorthEastOutside')
%text(10, -12.43, 'total no. of MG Sweeps');

 hold on;
 text(-5.5, -7.15, 'log(r_n/r_0)');
 
 
figure  
%plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b*'); 
%hold on;
resid_1=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_finegrid25Vcyle_gmres.dat'); %
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_slip_finegrid25Vcyle_gmres.dat'); %
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_periodic_finegrid25Vcyle_gmres.dat'); %

plot([0;resid_1(1:end,1)/25], [0;log10(resid_1(1:end,2))],  'b*'); 
hold on;
plot([0;resid_2(1:end,1)/25], [0;log10(resid_2(1:end,2))],  'r*'); 
hold on;
plot([0;resid_3(1:end,1)/25], [0;log10(resid_3(1:end,2))],  'k*');  
hold on;

%hold on;
% legend('Diri','Diri Vel MG','Diri Pre MG',  'x-Peri', 'xy-Peri', 'Peri', 'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(10, -14.93, 'total no. of MG Sweeps');
% axis square 

% legend('Diri Vel MG','Diri Pre MG',  'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(4, -14.93, 'total no. of MG Sweeps');
legend('Diri','Slip', 'Per' ,'Location','NorthEast')
text(10, -16.13, 'total no. GMRES It');
axis square 
%legend('Pure Dirichlet 1 MG/it', 'x- Periodic 1 MG/it', 'x,y -Periodic 1 MG/it', 'Pure Periodic 1 MG/it', 'Pure Periodic 2 MG/it', 'Pure Periodic 4 MG/it', 'Pure Periodic 10 MG/it', 'Location','NorthEastOutside')
%text(10, -12.43, 'total no. of MG Sweeps');

hold on;
text(-5.5, -7.15, 'log(r_n/r_0)');
 
  
 

figure  % Dirichlet using different V-Cycle numbers

 % all Dirichlet, Vel MG
resid_1=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_1VcylePerIt_finegrid.dat');  %  
%
% % all Dirichlet, Pres MG
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_2VcylePerIt_finegrid.dat'); %
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_4VcylePerIt_finegrid.dat'); %

resid_4=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_dirichelt_finegrid25Vcyle_gmres.dat'); %

%resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/3d_periodic_1VcylePerIt_finegrid.dat'); %


%plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b*'); 
%hold on;
plot([0;resid_1(1:end,1)/1], [0;log10(resid_1(1:end,2))],  'b*'); 
hold on;
plot([0;resid_2(1:end,1)/2], [0;log10(resid_2(1:end,2))],  'r*'); 
hold on;
plot([0;resid_3(1:end,1)/4], [0;log10(resid_3(1:end,2))],  'k*');  
hold on;
plot([0;resid_4(1:end,1)/25], [0;log10(resid_4(1:end,2))],  'c*');  
hold on;


%hold on;
% legend('Diri','Diri Vel MG','Diri Pre MG',  'x-Peri', 'xy-Peri', 'Peri', 'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(10, -14.93, 'total no. of MG Sweeps');
% axis square 

% legend('Diri Vel MG','Diri Pre MG',  'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(4, -14.93, 'total no. of MG Sweeps');
legend('Diri 1-Vcycle/it','Diri 2-Vcycle/it', 'Diri 4-Vcycle/it', 'Diri 25-Vcycle/it', 'Location','NorthEast')
 text(10, -15.23, 'total no. GMRES It');
axis square 
%legend('Pure Dirichlet 1 MG/it', 'x- Periodic 1 MG/it', 'x,y -Periodic 1 MG/it', 'Pure Periodic 1 MG/it', 'Pure Periodic 2 MG/it', 'Pure Periodic 4 MG/it', 'Pure Periodic 10 MG/it', 'Location','NorthEastOutside')
%text(10, -12.43, 'total no. of MG Sweeps');

 hold on;
 text(-5.5, -7.15, 'log(r_n/r_0)');
 
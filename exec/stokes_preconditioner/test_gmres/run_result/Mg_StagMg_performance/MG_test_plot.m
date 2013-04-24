function MG_test_plot
%let the residual tol very small, while change MG sweeps large, say 2--10. then compare the convergence behavior
%check the log(res) Vs itations
%KSP 
% fgmres

% 3D case
% the below is 3D case: 64*64*64, true norm is used
%Here we use the real residual norm, not KSP residual

disp('MG tests using different bottom smooth steps');

figure
% all Dirichlet, Vel MG
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_2botsmooth_vel.dat');
% 
% % all Dirichlet, Pres MG
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_2botsmooth_pres.dat');
     
%plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b*'); 
%hold on;
plot([0;resid_2(1:end,1)], [0;log10(resid_2(1:end,2))],  'r*'); 
hold on;
plot([0;resid_3(1:end,1)], [0;log10(resid_3(1:end,2))],  'k*');  
hold on;

% legend('Diri Vel MG','Diri Pre MG',  'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(4, -14.93, 'total no. of MG Sweeps');
legend('Diri Vel MG','Diri Pre MG',  'Location','NorthEast')
text(4, -14.93, 'total no. of MG Sweeps');
axis square 

hold on;
text(-2.1, -7, 'log(r_n/r_0)');


figure
% all Dirichlet, Vel MG
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_4botsmooth_vel.dat');
% 
% % all Dirichlet, Pres MG
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_4botsmooth_pres.dat');
     
%plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b*'); 
%hold on;
plot([0;resid_2(1:end,1)], [0;log10(resid_2(1:end,2))],  'r*'); 
hold on;
plot([0;resid_3(1:end,1)], [0;log10(resid_3(1:end,2))],  'k*');  
hold on;
% legend('Diri Vel MG','Diri Pre MG',  'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(4, -14.93, 'total no. of MG Sweeps');
legend('Diri Vel MG','Diri Pre MG',  'Location','NorthEast')
text(4, -14.93, 'total no. of MG Sweeps');
axis square 
hold on;
text(-2.1, -7, 'log(r_n/r_0)');



figure
% all Dirichlet, Vel MG
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_8botsmooth_vel.dat');
% 
% % all Dirichlet, Pres MG
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_8botsmooth_pres.dat');
     
%plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b*'); 
%hold on;
plot([0;resid_2(1:end,1)], [0;log10(resid_2(1:end,2))],  'r*'); 
hold on;
plot([0;resid_3(1:end,1)], [0;log10(resid_3(1:end,2))],  'k*');  
hold on;
% legend('Diri Vel MG','Diri Pre MG',  'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(4, -14.93, 'total no. of MG Sweeps');
legend('Diri Vel MG','Diri Pre MG',  'Location','NorthEast')
text(4, -14.93, 'total no. of MG Sweeps');
axis square 
hold on;
text(-2.1, -7, 'log(r_n/r_0)');


figure % the following comments maybe not correct
resid_1=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_2*2vel.dat');
% 
% % all Dirichlet, Pres MG
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_4*4vel.dat');
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_8*8vel.dat');

plot([0;resid_1(1:end,1)], [0;log10(resid_1(1:end,2))],  'b*'); 
hold on;
plot([0;resid_2(1:end,1)], [0;log10(resid_2(1:end,2))],  'r*'); 
hold on;
plot([0;resid_3(1:end,1)], [0;log10(resid_3(1:end,2))],  'k*');  
hold on;
 
legend('stag min width=2','stag min width=4','stag min width=8',  'Location','NorthEast')
text(14, -14.93, 'total no. of MG Sweeps');
axis square 

%legend('Pure Dirichlet 1 MG/it', 'x- Periodic 1 MG/it', 'x,y -Periodic 1 MG/it', 'Pure Periodic 1 MG/it', 'Pure Periodic 2 MG/it', 'Pure Periodic 4 MG/it', 'Pure Periodic 10 MG/it', 'Location','NorthEastOutside')
%text(10, -12.43, 'total no. of MG Sweeps');

 hold on;
 text(-4.1, -7, 'log(r_n/r_0)');
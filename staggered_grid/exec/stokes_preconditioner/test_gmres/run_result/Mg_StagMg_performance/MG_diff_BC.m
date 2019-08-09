function MG_diff_BC

% constant density and viscosity case

figure
% all Dirichlet, Vel MG
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_8botsmooth_vel.dat');
% 
% % all Dirichlet, Pres MG
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_8botsmooth_pres.dat');
     
% all slip, Vel MG
resid_4=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_slip_precon_finegridGS_8botsmooth_vel.dat');
% 
% % all slip, Pres MG
resid_5=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_slip_precon_finegridGS_8botsmooth_pres.dat');

% all periodic, Vel MG
resid_6=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_periodic_precon_finegridGS_8botsmooth_vel.dat');
% 
% % all periodic, Pres MG
resid_7=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_periodic_precon_finegridGS_8botsmooth_pres.dat');


%plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b*'); 
%hold on;
plot([0;resid_2(1:end,1)], [0;log10(resid_2(1:end,2))],  'r*'); 
hold on;
plot([0;resid_3(1:end,1)], [0;log10(resid_3(1:end,2))],  'k*');  
hold on;
plot([0;resid_4(1:end,1)], [0;log10(resid_4(1:end,2))],  'rs'); 
hold on;
plot([0;resid_5(1:end,1)], [0;log10(resid_5(1:end,2))],  'ks');  
hold on;

plot([0;resid_6(1:end,1)], [0;log10(resid_6(1:end,2))],  'ro'); 
hold on;
plot([0;resid_7(1:end,1)], [0;log10(resid_7(1:end,2))],  'ko');  
hold on;
% legend('Diri Vel MG','Diri Pre MG',  'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(4, -14.93, 'total no. of MG Sweeps');
legend('Diri Vel MG','Diri Pre MG', 'Slip Vel MG','Slip Pre MG','Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
text(4, -14.93, 'total no. of MG Sweeps');
axis square 
hold on;
text(-2.1, -7, 'log(r_n/r_0)');


% variable density and variable viscosity case
figure
% all Dirichlet, Vel MG
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_var_finegridGS_8botsmooth_vel.dat');
% 
% % all Dirichlet, Pres MG
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_var_finegridGS_8botsmooth_pres.dat');
     
% all slip, Vel MG
resid_4=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_slip_precon_var_finegridGS_8botsmooth_vel.dat');
% 
% % all slip, Pres MG
resid_5=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_slip_precon_var_finegridGS_8botsmooth_pres.dat');

% all periodic, Vel MG
resid_6=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_periodic_precon_var_finegridGS_8botsmooth_vel.dat');
% 
% % all periodic, Pres MG
resid_7=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_periodic_precon_var_finegridGS_8botsmooth_pres.dat');


%plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b*'); 
%hold on;
plot([0;resid_2(1:end,1)], [0;log10(resid_2(1:end,2))],  'r*'); 
hold on;
plot([0;resid_3(1:end,1)], [0;log10(resid_3(1:end,2))],  'k*');  
hold on;
plot([0;resid_4(1:end,1)], [0;log10(resid_4(1:end,2))],  'rs'); 
hold on;
plot([0;resid_5(1:end,1)], [0;log10(resid_5(1:end,2))],  'ks');  
hold on;

plot([0;resid_6(1:end,1)], [0;log10(resid_6(1:end,2))],  'ro'); 
hold on;
plot([0;resid_7(1:end,1)], [0;log10(resid_7(1:end,2))],  'ko');  
hold on;
% legend('Diri Vel MG','Diri Pre MG',  'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(4, -14.93, 'total no. of MG Sweeps');
legend('Diri Vel MG','Diri Pre MG', 'Slip Vel MG','Slip Pre MG','Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
text(8, -15.93, 'total no. of MG Sweeps');
axis square 
hold on;
text(-4.1, -7, 'log(r_n/r_0)');

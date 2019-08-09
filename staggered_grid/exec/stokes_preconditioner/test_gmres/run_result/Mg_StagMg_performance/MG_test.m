function MG_test
%let the residual tol very small, while change MG sweeps large, say 2--10. then compare the convergence behavior
%check the log(res) Vs itations
%KSP 
% fgmres

% 3D case
% the below is 3D case: 64*64*64, true norm is used
%Here we use the real residual norm, not KSP residual

% All Dirichlet
%resid_1=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_1VcylePerIt_finegrid.dat');

% all Dirichlet, Vel MG
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_vel.dat');  %  
% 
% % all Dirichlet, Pres MG
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_pres.dat'); %
     
% % x periodic, yz Dirichlet num_MG_sweep=1
% %resid_4=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_x-per-yzdirichelt_1VcylePerIt_finegrid.dat');
% 
% % x, y periodic, z Dirichlet num_MG_sweep=1     
% %resid_5=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_xy-per-zdirichelt_1VcylePerIt_finegrid.dat');
% 
% % x, y, z periodic, num_MG_sweep=1
% %resid_6=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_periodic_1VcylePerIt_finegrid.dat');
% 
% % x, y, z periodic, num_MG_sweep=1, velocity
% resid_7=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_periodic_precon_finegridGS_vel.dat');
% 
% % x, y, z periodic, num_MG_sweep=1, pressure
% resid_8=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_periodic_precon_finegridGS_pre.dat');


%resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_4*4vel.dat');
%resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/slip_mg_2*2vel.dat');

% all Dirichlet, Pres MG
%resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_dirichelt_precon_finegridGS_4*4pres.dat');
%resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/slip_mg_2*2pres.dat');


% x periodic, yz Dirichlet num_MG_sweep=1
%resid_4=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_x-per-yzdirichelt_1VcylePerIt_finegrid.dat');

% x, y periodic, z Dirichlet num_MG_sweep=1     
%resid_5=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_xy-per-zdirichelt_1VcylePerIt_finegrid.dat');

% x, y, z periodic, num_MG_sweep=1
%resid_6=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_periodic_1VcylePerIt_finegrid.dat');

% x, y, z periodic, num_MG_sweep=1, velocity
%resid_7=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_periodic_precon_finegridGS_4*4vel.dat');

% x, y, z periodic, num_MG_sweep=1, pressure
%resid_8=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/Mg_StagMg_performance/3d_periodic_precon_finegridGS_4*4pre.dat');   

figure  
%plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b*'); 
%hold on;
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
legend('Diri Vel MG','Diri Pre MG',  'Location','NorthEast')
 text(4, -14.93, 'total no. of MG Sweeps');
axis square 
%legend('Pure Dirichlet 1 MG/it', 'x- Periodic 1 MG/it', 'x,y -Periodic 1 MG/it', 'Pure Periodic 1 MG/it', 'Pure Periodic 2 MG/it', 'Pure Periodic 4 MG/it', 'Pure Periodic 10 MG/it', 'Location','NorthEastOutside')
%text(10, -12.43, 'total no. of MG Sweeps');

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
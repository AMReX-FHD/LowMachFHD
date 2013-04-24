function compare_P1P2_var_vis
figure  % Dirichlet using different V-Cycle numbers

 % all Dirichlet, use 1Vcyle/it 
resid_1=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/random/3d_dirichlet_GS_P1_steady_vis1.dat');  %  
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/random/3d_dirichlet_GS_P2_steady_vis1.dat'); %
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/random/3d_dirichlet_Jacobi_P1_steady_vis1.dat'); %
resid_4=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/random/3d_dirichlet_Jacobi_P2_steady_vis1.dat'); %


%plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b*'); 
%hold on;
plot([0;resid_1(1:end,1)/1], [0;log10(resid_1(1:end,2))],  'b*'); 
hold on;
plot([0;resid_2(1:end,1)/1], [0;log10(resid_2(1:end,2))],  'r*'); 
hold on;
plot([0;resid_3(1:end,1)/1], [0;log10(resid_3(1:end,2))],  'k*');  
hold on;
plot([0;resid_4(1:end,1)/1], [0;log10(resid_4(1:end,2))],  'c*');  
hold on;


%hold on;
% legend('Diri','Diri Vel MG','Diri Pre MG',  'x-Peri', 'xy-Peri', 'Peri', 'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(10, -14.93, 'total no. of MG Sweeps');
% axis square 

% legend('Diri Vel MG','Diri Pre MG',  'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(4, -14.93, 'total no. of MG Sweeps');
legend('P1, GS smooth','P2, GS smooth', 'P1, Jacobi smooth', 'P2, Jacobi smooth','Location','NorthEast')
 text(25, -15.23, 'total no. GMRES It');
axis square 
%legend('Pure Dirichlet 1 MG/it', 'x- Periodic 1 MG/it', 'x,y -Periodic 1 MG/it', 'Pure Periodic 1 MG/it', 'Pure Periodic 2 MG/it', 'Pure Periodic 4 MG/it', 'Pure Periodic 10 MG/it', 'Location','NorthEastOutside')
%text(10, -12.43, 'total no. of MG Sweeps');

 hold on;
 text(-16.5, -7.15, 'log(r_n/r_0)');
 
 
 figure  % Dirichlet using different V-Cycle numbers

 % all Dirichlet, use 1Vcyle/it 
resid_1=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/random/3d_dirichlet_GS_P1_steady_vis2.dat');  %  
resid_2=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/random/3d_dirichlet_GS_P2_steady_vis2.dat'); %
resid_3=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/random/3d_dirichlet_Jacobi_P1_steady_vis2.dat'); %
resid_4=load('/home/mccai/src/FluctHydro/stokes_preconditioner/test_gmres/run_result/fig2/random/3d_dirichlet_Jacobi_P2_steady_vis2.dat'); %


%plot([0;resid_1(:,1)], [0;log10(resid_1(:,2))],  'b*'); 
%hold on;
plot([0;resid_1(1:end,1)/1], [0;log10(resid_1(1:end,2))],  'b*'); 
hold on;
plot([0;resid_2(1:end,1)/1], [0;log10(resid_2(1:end,2))],  'r*'); 
hold on;
plot([0;resid_3(1:end,1)/1], [0;log10(resid_3(1:end,2))],  'k*');  
hold on;
plot([0;resid_4(1:end,1)/1], [0;log10(resid_4(1:end,2))],  'c*');  
hold on;


%hold on;
% legend('Diri','Diri Vel MG','Diri Pre MG',  'x-Peri', 'xy-Peri', 'Peri', 'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(10, -14.93, 'total no. of MG Sweeps');
% axis square 

% legend('Diri Vel MG','Diri Pre MG',  'Peri Vel MG','Peri Pre MG', 'Location','NorthEast')
%  text(4, -14.93, 'total no. of MG Sweeps');
legend('P1, GS smooth','P2, GS smooth', 'P1, Jacobi smooth', 'P2, Jacobi smooth','Location','NorthEast')
 text(25, -15.23, 'total no. GMRES It');
axis square 
%legend('Pure Dirichlet 1 MG/it', 'x- Periodic 1 MG/it', 'x,y -Periodic 1 MG/it', 'Pure Periodic 1 MG/it', 'Pure Periodic 2 MG/it', 'Pure Periodic 4 MG/it', 'Pure Periodic 10 MG/it', 'Location','NorthEastOutside')
%text(10, -12.43, 'total no. of MG Sweeps');

 hold on;
 text(-16.5, -7.15, 'log(r_n/r_0)');
 
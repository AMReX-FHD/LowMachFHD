% Evaluate parameters from Lemarchand and Nowakowski
% paper "Do the internal fluctuations blur..."
clear all
k_1 = 4; k_2 = 1.37; k_3 = 1; k_m3 = 10;
D_A = 5; D_B = 50;
%D_A = 25; D_B = 250;
A_0 = 0;
B_0 = k_m3/k_3;
Delta = k_m3^2 - 4*k_1^2*k_3/k_2;
A_p = (k_m3 + sqrt(Delta))/(2*k_1);
A_m = (k_m3 - sqrt(Delta))/(2*k_1);
B_p = (k_m3 - k_1*A_p)/k_3;
B_m = (k_m3 - k_1*A_m)/k_3;
fprintf('Steady states:\n');
fprintf('A_0 = %g, B_0 = %g\n',A_0,B_0);
fprintf('A_+ = %g, B_+ = %g\n',A_p,B_p);
fprintf('A_- = %g, B_- = %g\n',A_m,B_m);
a(1,1) = - k_1 + 2*k_2*A_p*B_p;
a(1,2) = k_2*A_p^2;
a(2,1) = -2*k_2*A_p*B_p;
a(2,2) = -k_3 - k_2*A_p^2;
fprintf('Matrix a(i,j) = df_i/dx_j :\n');
disp(a)
ll_1 = sqrt(D_A/a(1,1));
ll_2 = sqrt(D_B/-a(2,2));
fprintf('Penetration lengths per species:\n');
fprintf(' ell_1 = %g, ell_2 = %g \n',ll_1, ll_2);

rates=eig(a);
ll_g=sqrt(-D_A/real(rates(1)));
fprintf('Global penetration length:\n');
fprintf(' ell_1 = %g \n', ll_g);

q_m = sqrt( 0.5*(1/ll_1^2 - 1/ll_2^2) );
ll_m = 2*pi/q_m;
fprintf('Turing wavelength = %g\n',ll_m);

ratio = ll_m / ll_g

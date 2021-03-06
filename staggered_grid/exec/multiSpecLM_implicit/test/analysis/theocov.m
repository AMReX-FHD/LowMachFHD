% Code to match the theoretical covariance with the low-mach code calculation

clear all; format long;
rhotot = 3; c_init(1,1:3) = rhotot*[0.2 0.35 0.45];
molmass(1:3) = [1.d0 2.d0 3.d0]; 
dx = 1; dy = 1; dz = 1; dV = dx*dy*dz;
variance_coef_mass = 1e-6;

covW_theo(1,1) = variance_coef_mass*(c_init(1,1)*(molmass(2)*c_init(1,1)*c_init(1,2) + molmass(3)*c_init(1,1)* ...
                 c_init(1,3) + molmass(1)*(c_init(1,2) + c_init(1,3))^2))/(dV* ...
                 (c_init(1,1) + c_init(1,2) + c_init(1,3))^4);
covW_theo(1,2) = -variance_coef_mass*((c_init(1,1)*c_init(1,2)*(-(molmass(3)*c_init(1,3)) + molmass(2)*(c_init(1,1) + ...
                 c_init(1,3)) + molmass(1)*(c_init(1,2) + c_init(1,3))))/(dV*( ...
                 c_init(1,1) + c_init(1,2) + c_init(1,3))^4));
covW_theo(1,3) = -variance_coef_mass*((c_init(1,1)*c_init(1,3)*(-(molmass(2)*c_init(1,2)) + molmass(3)*(c_init(1,1) + ...
                 c_init(1,2)) + molmass(1)*(c_init(1,2) + c_init(1,3))))/(dV*( ...
                 c_init(1,1) + c_init(1,2) + c_init(1,3))^4));
covW_theo(2,1) = covW_theo(1,2);
covW_theo(2,2) = variance_coef_mass*(c_init(1,2)*(molmass(2)*(c_init(1,1) + c_init(1,3))^2 + c_init(1,2)*(molmass(1)* ...
                 c_init(1,1) + molmass(3)*c_init(1,3))))/(dV*(c_init(1,1) + ...
                         c_init(1,2) + c_init(1,3))^4); 
covW_theo(2,3) = -variance_coef_mass*((c_init(1,2)*c_init(1,3)*(-(molmass(1)*c_init(1,1)) + molmass(3)*(c_init(1,1) + ...
                 c_init(1,2)) + molmass(2)*(c_init(1,1) + c_init(1,3))))/(dV*( ...
                 c_init(1,1) + c_init(1,2) + c_init(1,3))^4));
covW_theo(3,1) = covW_theo(1,3); 
covW_theo(3,2) = covW_theo(2,3); 
covW_theo(3,3) = variance_coef_mass*(c_init(1,3)*(molmass(3)*(c_init(1,1) + c_init(1,2))^2 + (molmass(1)*c_init(1,1) + ...
                 molmass(2)*c_init(1,2))*c_init(1,3)))/(dV*(c_init(1,1) + ...
                 c_init(1,2) + c_init(1,3))^4);

Theo_cov=[covW_theo(1,1) covW_theo(1,2) covW_theo(1,3); covW_theo(2,1) covW_theo(2,2) covW_theo(2,3); covW_theo(3,1) covW_theo(3,2) covW_theo(3,3)]



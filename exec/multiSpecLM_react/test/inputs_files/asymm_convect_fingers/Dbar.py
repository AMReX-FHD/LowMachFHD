# cgs units
# 1 = HCl, NaOH, NaCl, H2O

# tracer diffusion coefficients (from de Wit's paper)
D1 = 3.336e-5
D2 = 2.129e-5
D3 = 1.611e-5

# self-diffusion coefficient of water (from Phys. Fluids paper)
D4 = 2.3e-5

# maxwell-stefan
D14 = D1
D24 = D2
D34 = D3
D12 = D1*D2/D4
D13 = D1*D3/D4
D23 = D2*D3/D4

print "(D1,D2,D3,D4) = (%e,%e,%e,%e)" % (D1,D2,D3,D4)
print "(D12,D13,D23,D14,D24,D34) = (%e,%e,%e,%e,%e,%e)" % (D12,D13,D23,D14,D24,D34)

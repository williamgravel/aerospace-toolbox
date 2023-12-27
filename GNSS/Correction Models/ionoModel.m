function [rho_IF,I_a,I_b,TEC] = ionoModel(rho_a,f_a,rho_b,f_b)

TEC = f_a^2*f_b^2/(40.3*(f_a^2 - f_b^2))*(rho_b - rho_a);

rho_IF = (f_a^2*rho_a - f_b^2*rho_b)/(f_a^2 - f_b^2);

I_a = 40.3*TEC/f_a^2;
I_b = 40.3*TEC/f_b^2;

end

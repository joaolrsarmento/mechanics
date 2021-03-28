% Instituto Tecnológico de Aeronáutica
% Engenharia Mecânica-Aeronáutica
% MPD-11 - Dinâmica de máquinas
% Aluno: João Sarmento
% 
% Funções: Descrição do Mecanismo (mecanismo de quatro barras)

function [f] = f_quatro_barras(x, c)

theta_3 = x(1); theta_4 = x(2); r_1 = c(1); r_2 = c(2); r_3 = c(3); r_4 = c(4); theta_1 = c(5); theta_2 = c(6);

f_1 = r_1*cos(theta_1) + r_4*cos(theta_4) - r_2*cos(theta_2) - r_3*cos(theta_3);

f_2 = r_1*sin(theta_1) + r_4*sin(theta_4) - r_2*sin(theta_2) - r_3*sin(theta_3);

f = [ f_1 ; f_2 ];
end


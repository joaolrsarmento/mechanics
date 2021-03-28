% Instituto Tecnológico de Aeronáutica
% Engenharia Mecânica-Aeronáutica
% MPD-11 - Dinâmica de máquinas
% Aluno: João Sarmento
% 
% Funções: Jacobiana (quatro barras)

function [J] = j_quatro_barras(x, c)

theta_3 = x(1); theta_4 = x(2); r_1 = c(1); r_2 = c(2); r_3 = c(3); r_4 = c(4); theta_1 = c(5); theta_2 = c(6);

df_1_dx_1 =  r_3*sin(theta_3);
df_1_dx_2 = -r_4*sin(theta_4);
df_2_dx_1 = -r_3*cos(theta_3);
df_2_dx_2 =  r_4*cos(theta_4);

J = [ df_1_dx_1 , df_1_dx_2 ; ...
      df_2_dx_1 , df_2_dx_2 ];
end

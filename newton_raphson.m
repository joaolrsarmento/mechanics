% Instituto Tecnológico de Aeronáutica
% Engenharia Mecânica-Aeronáutica
% MPD-11 - Dinâmica de máquinas
% Aluno: João Sarmento
% 
% Função: Newton Raphson

function [x] = newton_raphson(F, J, x0, c)

epsilon  = 1e-8;
epsilon_f  = epsilon;
epsilon_dx = epsilon;
iter_max = 50;

iter = 0;
x = x0;
f = F(x, c);

while norm(f, 2) > epsilon_f
    dx = -J(x, c)\f;
    if norm(dx, 2) < epsilon_dx
        error("Solução não convergiu!");
    end
    x = x + dx;
    f = F(x, c);
    iter = iter + 1;
    if iter > iter_max
        error(['Solução não convergiu! Tente aumentar o ' ...
               'número máximo de iterações permitidas!']);
    end
end
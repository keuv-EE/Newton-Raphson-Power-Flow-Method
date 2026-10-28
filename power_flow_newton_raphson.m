%% ============================================================

clear; clc; close all;

%% ============================================================
%
% Cabecalho
%
%% ============================================================

% Quantidade de barras
n_barras = input("Quantidade de Barras 'n': ");

% Premissas: coletar informacoes das barras de i = 1:n_barras
% Contudo, atentar-se a dimensionalidade funcao da quantidade de barras PV

% Tensoes na barra i
V = zeros(1,n_barras);

for i = 1:n_barras
    V(i) = input(sprintf('Tensão na barra %d: ', i));
end

% Tipo da barra i
tipo_barra = string(zeros(1, n_barras));
for i = 1:n_barras
    tipo_barra(i) = input(sprintf("Qual é o tipo da barra %d? ", i), "s");
end

%% Ybus
y = zeros(n_barras, n_barras+1);     % admitancias entre as barras
for i = 1:n_barras
    for j = 1:n_barras+1
        y(i,j) = input(sprintf('y_{%d, %d}: ', i, j-1));
    end
end

% Matriz Ybus
Y = zeros(n_barras, n_barras);

for i = 1:n_barras
    for j = 1:n_barras+1
        % j-1 = 0:n_barras
        if i == j-1
            for p = 1:n_barras+1
                if i ~= p-1
                    Y(i,i) = Y(i,i) + y(i,p);
                end
            end

        elseif j-1 > 0
            Y(i,j-1) = -y(i,j);
        end
        % elseif j == 1 (Yi0 -> faca nada)
    
    end
end

%% Potencias

% sch
S = zeros(1, n_barras);     % Potencias que nao devem ser atualizadas

for i = 1:n_barras
    S(i) = input(sprintf('Potência (S_sch) na barra %d: ', i));
end

erro_ang_abs_tensao = 1e5;
sc = 1e-5;

while any(erro_ang_abs_tensao > sc)
    Pk = zeros(1, n_barras);
    Qk = zeros(1, n_barras);
    
    for i = 1:n_barras
        for j = 1:n_barras
            if tipo_barra(i) == "PQ"
                Pk(i) = Pk(i) + abs(V(i))*abs(V(j))*abs(Y(i,j))*cos(angle(Y(i,j)) - angle(V(i)) + angle(V(j)));
    
                Qk(i) = Qk(i) - abs(V(i))*abs(V(j))*abs(Y(i,j))*sin(angle(Y(i,j)) - angle(V(i)) + angle(V(j)));

            elseif tipo_barra(i) == "PV"
                Pk(i) = Pk(i) + abs(V(i))*abs(V(j))*abs(Y(i,j))*cos(angle(Y(i,j)) - angle(V(i)) + angle(V(j)));
            end
        end
    end
    
    qntdd_PV = 0;
    
    for i = 1:n_barras
        if tipo_barra(i) == "PV"
            qntdd_PV = qntdd_PV + 1;
        end
    end
    
    delt_Pk = zeros(1, n_barras);
    delt_Qk = zeros(1, n_barras);
    
    for i = 1:n_barras
        if tipo_barra(i) == "PQ"
            delt_Pk(i) = real(S(i)) - Pk(i);
            delt_Qk(i) = imag(S(i)) - Qk(i);
        elseif tipo_barra(i) == "PV"
            delt_Pk(i) = real(S(i)) - Pk(i);
        end
    end
    
    delt_Pk(:, 1, :) = [];
    
    inx_aux = 0;
    for i = 1:n_barras
        if tipo_barra(i) == "PV" || tipo_barra(i) == "Slack"
            delt_Qk(:, i - inx_aux, :) = [];
            inx_aux = inx_aux + 1;
        end
    end
    
    erro_Pk_Qk = [delt_Pk.'; delt_Qk.'];
    
    %% Matriz Jacobiana
    
        J1 = zeros(n_barras - 1, n_barras - 1);
        J2 = zeros(n_barras - 1, n_barras - 1);
        J3 = zeros(n_barras - 1, n_barras - 1);
        J4 = zeros(n_barras - 1, n_barras - 1);
    
        for i = 2:n_barras
            for j = 2:n_barras
    
                if i == j
                    for p = 1:n_barras
                        if i ~= p
                            J1(i-1,i-1) = J1(i-1,i-1) + abs(V(i))*abs(V(p))*abs(Y(i,p))*sin(angle(Y(i,p)) -  angle(V(i)) + angle(V(p)));
                            J2(i-1,i-1) = J2(i-1,i-1) + abs(V(p))*abs(Y(i,p))*cos(angle(Y(i,p)) -  angle(V(i)) + angle(V(p)));
                            J3(i-1,i-1) = J3(i-1,i-1) + abs(V(i))*abs(V(p))*abs(Y(i,p))*cos(angle(Y(i,p)) -  angle(V(i)) + angle(V(p)));
                            J4(i-1,i-1) = J4(i-1,i-1) - abs(V(p))*abs(Y(i,p))*sin(angle(Y(i,p)) -  angle(V(i)) + angle(V(p)));
                        end
                    end
                    J2(i-1,i-1) = J2(i-1,i-1) + 2*abs(V(i))*abs(Y(i,i))*cos(angle(Y(i,i)));
                    J4(i-1,i-1) = J4(i-1,i-1) - 2*abs(V(i))*abs(Y(i,i))*sin(angle(Y(i,i)));
                elseif i ~= j
                    J1(i-1,j-1) = - abs(V(i))*abs(V(j))*abs(Y(i,j))*sin(angle(Y(i,j)) - angle(V(i)) + angle(V(j)));
                    J2(i-1,j-1) =   abs(V(i))*abs(Y(i,j))*cos(angle(Y(i,j)) - angle(V(i)) + angle(V(j)));
                    J3(i-1,j-1) = - abs(V(i))*abs(V(j))*abs(Y(i,j))*cos(angle(Y(i,j)) - angle(V(i)) + angle(V(j)));
                    J4(i-1,j-1) = - abs(V(i))*abs(Y(i,j))*sin(angle(Y(i,j)) - angle(V(i)) + angle(V(j)));
                end
    
            end
        end
    
        J = [J1, J2; J3, J4];
    
    for i = n_barras:-1:1
        if tipo_barra(i) == "PQ"
            continue
        elseif tipo_barra(i) == "PV"
            J(:, n_barras-1 + i-1) = [];
            J(n_barras-1 + i-1, :) = [];
        end
    end
    
    erro_ang_abs_tensao = J \ erro_Pk_Qk;
    
    erro_ang = [0; erro_ang_abs_tensao(1:n_barras-1)];
    
    erro_abs_tensao = erro_ang_abs_tensao(n_barras:end);
    
    for i = 1:n_barras
        if tipo_barra(i) == "PV"
            erro_abs_tensao = [erro_abs_tensao(1:i-1); 0; erro_abs_tensao(i:end)];
        elseif tipo_barra(i) == "Slack"
            erro_abs_tensao = [0; erro_abs_tensao(i:end)];
        end
    end
    
    tens = zeros(1, n_barras);
    
    for i = 1:n_barras
        tens(i) = abs(V(i));
    end
    
    inx_aux = 0;
    for i = 1:n_barras
        if tipo_barra(i) == "PV" || tipo_barra(i) == "Slack"
            tens(i - inx_aux) = [];
            inx_aux = inx_aux + 1;
        end
    end
    
    ang = zeros(1, n_barras);
    
    for i = 1:n_barras
        ang(i) = angle(V(i));
    end
    
    ang(1) = []; % Barra de Slack
    
    ang_abs_tensao = [ang.'; tens.'] + erro_ang_abs_tensao;
    
    for i = 1:n_barras
        if tipo_barra(i) == "PQ"
            V(i) = abs(abs(V(i)) +  erro_abs_tensao(i) )*exp(1i * (angle(V(i)) + erro_ang(i) ));
        elseif tipo_barra(i) == "PV"
            V(i) = abs(V(i))*exp(1i * (angle(V(i)) + erro_ang(i) ));
        end
    end
    
    if all(erro_ang_abs_tensao < sc)
        % fprintf('Iter: %d\n', iter);
        for i = 1:n_barras
            fprintf('Modulo Tensao V%d: %.5f\n', i, abs(V(i)));
            fprintf('Fase Tensao V%d: %.5f\n', i, rad2deg(angle(V(i))));
        end
    end

end
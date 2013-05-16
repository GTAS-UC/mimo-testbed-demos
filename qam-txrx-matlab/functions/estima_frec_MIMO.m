function [X_corregidos,H_g_corregido] = estima_frec_MIMO(X,Pilotos,niter)
% ESTIMA_FREC: Dados unos observables (Pilotos + Info) y una estima de la
% matriz del canal corrige el offset de frecuencia en los observables y en
% la matriz del canal
%
% estima_frec:
%       [X_corregidos,H_g_corregido] = estima_frec(X,H_g,Longitud_pilotos,niter)
%           X = Observables (Pilotos + Info)
%           Pilotos = Pilotos enviados
%           H_g = Estima de la matriz del canal
%           niter = numero maximo de iteraciones
%
%           X_corregidos = Observables corregidos (Sin Pilotos)
%           H_g_corregido = Estima de la matriz del canal con el Offset
%           corregido

% David Ramirez Garcia
% G.T.A.S.
% 2008

[nr,N] = size(X);

n = 0:(N-1);

% Longitud_pilotos = [64*ones(1,5) 128*ones(1,5) 256*ones(1,5) 512*ones(1,niter-5)];
Longitud_pilotos = [64*ones(1,5) 128*ones(1,5) 256*ones(1,niter-5)];

H_g = X(:,1:Longitud_pilotos(1))*Pilotos(:,1:Longitud_pilotos(1)).'/Longitud_pilotos(1);

for k = 1:niter
    
    f1 = estimafrec(X(:,1:Longitud_pilotos(k)), Pilotos(:,1:Longitud_pilotos(k)), H_g, round(Longitud_pilotos(k)/2));
    
    % Corregimos el error
    
    for kk = 1:nr
        
        X(kk,:) = X(kk,:).*exp(-1i*f1*n);
        
    end
    
    % Reestimamos canal
    H_g = X(:,1:Longitud_pilotos(k))*Pilotos(:,1:Longitud_pilotos(k))'/Longitud_pilotos(k);
    
end

X_corregidos = X;

H_g_corregido = H_g;


function delta_f = estimafrec(x, a, H, N)
% ESTIMAFREC: Funcion que se encarga de calcular el offset de frecuencia en
% una transmision MIMO
%
% Sintaxis:
%    delta_f = estimafrec(x, a, H, N)
%        x = Pilotos Recibidos
%        a = Pilotos Enviados
%        H = Matriz de Canal
%        N = Numero de puntos para evaluar la autocorrelacion (length(x)/2)

% Gabriel García Ocejo
% G.T.A.S.
% 2007

L = length(x);
M = N;
%A sera una matriz de tamaño 2xLongitudPilotos
%En cada instante k se transmitiran A(:, k)
%H es la matriz del canal de tamaño 2x2
y = zeros(1,L);
for i=1:L
    y(i) = a(:,i)'*H'*x(:,i);
end

%Calculamos la autocorrelacion para un m;
R = zeros(1, M);
for m=1:M
    R(m)=sum(y(m+1:L).*conj(y(1:L-m)))/(L-m);
end
m = 1:M;
%plot(angle(R))
%pause
delta_f = sum(m.*(abs(R).^2).*angle(R)) / sum((m.^2).*(abs(R).^2));
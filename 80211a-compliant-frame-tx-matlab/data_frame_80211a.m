
% PAR�METROS
%   - tasa      ->	velocidad en Mbps (6,9,12,18,24,48,54)
%   - n_datos   ->	n�mero de octetos de datos a enviar

function trama = data_frame_80211a(tasa,data_hex)

if (tasa ~= 6) && (tasa ~= 9) && (tasa ~= 12) && (tasa ~= 18) && (tasa ~= 24) && (tasa ~= 36) && (tasa ~= 48) && (tasa ~= 54)
    error('La tasa no es v�lida'); 
end

% PAR�METROS SEG�N LA TASA

switch tasa
    
    case 6
        rate = '1101';          % BPSK
        Ncbps = 48;             % N� de bits por s�mbolo OFDM
        Nbpsc = 1;              % N� de bits codificados por subportadora  
        Ndbps = 24;             % N� de bits de datos por s�mbolo OFDM
        tipo_mod = 1;           % Tipo de modulaci�n
        punct = 1;              % Patr�n de perforado
    case 9
        rate = '1111';      % BPSK
        Ncbps = 48;
        Nbpsc = 1;
        Ndbps = 36;
        punct = [1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1];      % 3/4   
        tipo_mod = 1;
    case 12
        rate = '0101';      % QPSK
        Ncbps = 96;
        Nbpsc = 2;
        Ndbps = 48;
        tipo_mod = 2;
        punct = 1;
    case 18
        rate = '0111';      % QPSK
        Ncbps = 96;
        Nbpsc = 2;
        Ndbps = 72;
        punct = [1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1];      % 3/4
        tipo_mod = 2;
    case 24
        rate = '1001';      % 16QAM
        Ncbps = 192;
        Nbpsc = 4;
        Ndbps = 96;
        tipo_mod = 3;
        punct = 1;
    case 36
        rate = '1011';      % 16QAM
        Ncbps = 192;
        Nbpsc = 4;
        Ndbps = 144;
        punct = [1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1];      % 3/4
        tipo_mod = 3;
    case 48
        rate = '0001';      % 64QAM
        Ncbps = 288;
        Nbpsc = 6;
        Ndbps = 192;
        punct = [1 1 1 0 1 1 1 0 1 1 1 0];      % 2/3
        tipo_mod = 4;
    case 54
        rate = '0011';      % 64QAM
        Ncbps = 288;
        Nbpsc = 6;
        Ndbps = 216;
        punct = [1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1];      % 3/4
        tipo_mod = 4;
end

% SECUENCIAS DE ENTRENAMIENTO

% SECUENCIA CORTA
short_training_seq = ...
    sqrt(13/6)*[0 0 0 0 0 0 0 0 1+j 0 0 0 -1-j 0 0 0 1+j 0 0 0 -1-j 0 0 0 -1-j 0 0 0 1+j 0 0 0 ...
    0 0 0 0 -1-j 0 0 0 -1-j 0 0 0 1+j 0 0 0 1+j 0 0 0 1+j 0 0 0 1+j 0 0 0 0 0 0 0];
% Muestras temporales 
muestras_short = ifft(ifftshift(short_training_seq));
% Repetimos 2.5 veces
muestras_short_tot = [muestras_short muestras_short muestras_short(1:33)];
% Enventanamos
muestras_short_tot = [0.5*muestras_short_tot(1) muestras_short_tot(2:160) 0.5*muestras_short_tot(161)];

% SECUENCIA LARGA
long_training_seq = ...
    [0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 0 ...
     1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0];
% Muestras temporales
muestras_long = ifft(ifftshift(long_training_seq)); 
% A�adimos intervalo de guarda
cyclic_prefix = muestras_long(33:64);   % 2 * GI
muestras_long_tot = [cyclic_prefix muestras_long muestras_long muestras_long(1)];
muestras_long_tot = [0.5*muestras_long_tot(1) muestras_long_tot(2:160) 0.5*muestras_long_tot(161)]; 

% SECUENCIA DE ENTRENAMIENTO
training_sequence = [muestras_short_tot(1:160) muestras_short_tot(161)+ muestras_long_tot(1) muestras_long_tot(2:161)];

% CAMPO SIGNAL (modulado BPSK con R = 1/2) (sin service)

% rate      -> (4 bits)
% res       -> (1 bit) de valor 0
% long      -> (12 bits) Codifica el n�mero de bytes en la trama MAC.
%               LSB-MSB   -> Se mete en octetos
% parity    -> (1 bit) par
% tail      -> (6 bits) de valor 0
% service   -> (16 bits) se tx en el campo de datos al rate de la trama MAC,
%              los 1os 8 bits son 0

long=dec2bin(length(data_hex)/2,12);
long = double(long(end:-1:1))-48;           % long(LSB:MSB)
rate = double(rate)-48;
parity = mod(sum([rate long]),2);

% Formamos el campo SIGNAL
sig_field = [rate 0 long parity 0 0 0 0 0 0];

% CODIFICADOR CONVOLUCIONAL
trellis = poly2trellis(7,[133 171]);
sig_field_c = convenc(sig_field,trellis);

% INTERLEAVER

Ncbps_bpsk = 48;
Nbpsc_bpsk = 1;

k = 0:(Ncbps_bpsk-1);
s = max(Nbpsc_bpsk/2,1);

% 1� permutaci�n
ind_1 = (Ncbps_bpsk/16)*mod(k,16) + floor(k/16);
sig_field_c(ind_1+1) = sig_field_c;
% 2� permutaci�n
ind_1 = 0:(Ncbps_bpsk-1);
ind_2 = s*floor(ind_1/s) + mod(ind_1 + Ncbps_bpsk - floor(16*ind_1/Ncbps_bpsk),s);
sig_field_c_i(ind_2+1) = sig_field_c;

% MODULAMOS BPSK
sig_field_mod = 2*sig_field_c_i - 1;

% Completamos con 0's y pilotos

pilotos = [1 1 1 -1];

sig_field_mod = [0 0 0 0 0 0 sig_field_mod(1:5) pilotos(1) sig_field_mod(6:18) pilotos(2) ...
    sig_field_mod(19:24) 0 sig_field_mod(25:30) pilotos(3) sig_field_mod(31:43) ....
    pilotos(4) sig_field_mod(44:48) 0 0 0 0 0];

% Muestras temporales
sig_field_temp = ifft(ifftshift(sig_field_mod));
% A�adimos el prefijo c�clico
cyclic_prefix = sig_field_temp(49:64);
sig_field_temp = [0.5*cyclic_prefix(1) cyclic_prefix(2:16) sig_field_temp(1:64) 0.5*sig_field_temp(1)];

% Annex G frame
%data_hex = '0802002e006008cd37a60020d6013cf1006008ad3baf00004a6f792c2062726967687420737061726b206f6620646976696e6974792c0a4461756768746572206f6620456c797369756d2c0a466972652d696e73697265642077652074726561da5799ed';

data_hex=flipdim(reshape(data_hex,2,[]),1);
data_hex=data_hex(:);

data = zeros(1,length(data_hex)*4);

for k = 1:length(data_hex)

    aux = dec2bin(hex2dec(data_hex(k)),4);
    aux = double(aux)>48;
    aux = aux(end:-1:1);
    data((k-1)*4+1:k*4) = aux;
    
end

% Service + Data + Tail + Pad
datos = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 data 0 0 0 0 0 0];
aux = length(datos)/Ndbps;
% Calculamos el relleno;
pad = ceil(aux)*Ndbps - length(datos);
datos = [datos zeros(1,pad)];

% SCRAMBLER (secuencia de 127 bits)

secuencia = [0 1 1 0 1 1 0 0 0 0 0 1 1 0 0 1 1 0 1 0 1 0 0 1 1 ...
    1 0 0 1 1 1 1 0 1 1 0 1 0 0 0 0 1 0 1 0 1 0 1 1 1 1 1 0 1 0 0 1 0 1 ...
    0 0 0 1 1 0 1 1 1 0 0 0 1 1 1 1 1 1 1 0 0 0 0 1 1 1 0 1 1 1 1 0 0 1 ...
    0 1 1 0 0 1 0 0 1 0 0 0 0 0 0 1 0 0 0 1 0 0 1 1 0 0 0 1 0 1 1 1 0 1];

secuencia = repmat(secuencia,1,ceil(length(datos)/127));
secuencia = secuencia(1:length(datos));
datos_scr = xor(datos,secuencia);
datos_scr(16+length(data)+1:16+length(data)+6) = 0;            

% CODIFICADOR CONVOLUCIONAL

% Hacemos que el patr�n de perforado tenga longitud igual al n� de datos
% para el caso de 1/2.

if punct == 1
    punct = ones(1,length(datos_scr));
end

datos_scr_c = convenc(datos_scr,trellis,punct);

% Inicializamos el vector de datos
data_tot_temp = zeros(1,81 + 80*(length(datos_scr_c)/Ncbps-1));

% h)

for n_trama = 1:length(datos_scr_c)/Ncbps
    
    data_trama = datos_scr_c((n_trama-1)*Ncbps+1:n_trama*Ncbps);
    data_aux = zeros(1,length(data_trama));
    data_int = zeros(1,length(data_trama));
    
    % Definimos los pilotos iniciales para, posteriormente, hacer el
    % control de polaridad
    pilotos = [1 1 1 -1];
    
    % INTERLEAVER

    k = 0:(Ncbps - 1);
    s = max(Nbpsc/2,1);
    % 1� permutaci�n
    ind_1 = (Ncbps/16)*mod(k,16) + floor(k/16);
    data_aux(ind_1+1) = data_trama;
    % 2� permutaci�n
    ind_1 = 0:(Ncbps-1);
    ind_2 = s*floor(ind_1/s) + mod(ind_1 + Ncbps - floor(16*ind_1/Ncbps),s);
    data_int(ind_2+1) = data_aux;

    %% MODULAMOS

    switch tipo_mod

        case 1
            % BPSK
            data_mod = 2*data_int-1;            
            
        case 2
            % QPSK
            data_mod = (-1)*(data_int(1:2:length(data_int))==0) + 1*(data_int(1:2:length(data_int))==1) + (i)*(data_int(2:2:length(data_int))==1) + (-i)*(data_int(2:2:length(data_int))==0);
            data_mod = (1/sqrt(2))*data_mod;
            
        case 3
            % 16-QAM
            m_1 = modem.qammod('M', 16, 'SymbolOrder', 'Gray', 'InputType', 'Bit');
            data_mod = modulate(m_1,data_int');
            data_mod = (1/sqrt(10))*data_mod';
            
        case 4
            % 64-QAM
            m_2 = modem.qammod('M', 64, 'SymbolOrder', 'Gray', 'InputType', 'Bit');
            data_mod = modulate(m_2,data_int');
            data_mod = (1/sqrt(42))*data_mod';
    end

    control_polaridad = [1 1 1 -1 -1 -1 1 -1 -1 -1 -1 1 1 -1 1 -1 -1 1 1 -1 1 1 -1 1 1 1 1 1 1 -1 1 ... 
        1 1 -1 1 1 -1 -1 1 1 1 -1 1 -1 -1 -1 1 -1 1 -1 1 1 -1 -1 1 1 1 1 1 -1 -1 1 1 ... 
        -1 -1 1 -1 1 -1 1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 -1 -1 1 -1 1 1 1 1 -1 1 -1 1 -1 1 ...
        -1 -1 -1 -1 -1 1 -1 1 1 -1 1 -1 1 1 1 -1 -1 1 -1 -1 -1 1 1 1 -1 -1 -1 -1 -1 -1 -1];
    
    pilotos = pilotos * control_polaridad(n_trama);
    
    % A�adimos pilotos 
    data_mod = [0 0 0 0 0 0 data_mod(1:5) pilotos(1) data_mod(6:18) pilotos(2) ...
        data_mod(19:24) 0 data_mod(25:30) pilotos(3) data_mod(31:43) pilotos(4) ...
        data_mod(44:48) 0 0 0 0 0];

    % Sacamos muestras temporales
    data_temp = ifft(fftshift(data_mod));
    % A�adimos el prefijo c�clico
    cyclic_prefix = data_temp(49:64);
    data_temp = [0.5*cyclic_prefix(1) cyclic_prefix(2:16) data_temp(1:64) 0.5*data_temp(1)];    
    
    if n_trama==1
        data_tot_temp(1:81) = data_temp;
    else
        data_tot_temp((n_trama-1)*81 - (n_trama - 2)) = data_tot_temp((n_trama-1)*81 - (n_trama - 2)) + data_temp(1);
        data_tot_temp((n_trama-1)*81 - (n_trama - 2) + 1:(n_trama-1)*81 - (n_trama - 2) + 80) = data_temp(2:end);        
    end 

end

% TRAMA COMPLETA

trama = [training_sequence(1:320) training_sequence(321)+sig_field_temp(1) ...
    sig_field_temp(2:end-1) sig_field_temp(end)+data_tot_temp(1) ...
    data_tot_temp(2:end)];










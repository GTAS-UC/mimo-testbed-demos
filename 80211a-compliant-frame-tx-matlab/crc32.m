function crc=crc32(packet)
%generate CRC
m = uint32(hex2dec('ffffffff'));
crc_m = uint32(hex2dec('edb88320'));

if exist('crc32_table') == 0
    %disp('Creating crc32_table...');
    crc32_table = zeros(1,256,'uint32');
    for byte = 0:255
        crc = byte;
        for j = 1:8
            if bitand(crc,1) == 1
                mask = m;
            else
                mask = 0;
            end
            crc = bitxor(bitshift(crc,-1),bitand(crc_m, mask));
        end
        crc32_table(byte+1) = crc;
        dec2hex(crc32_table(byte+1));
    end
end

len = length(packet);
i = 1;
ff = uint32(hex2dec('ff'));
crc = m;

while i < len
    byte = uint32(hex2dec(packet(i:i+1)));
    t_index = bitand(bitxor(crc, byte), ff) + 1;
    crc = bitxor(bitshift(crc,-8), crc32_table(t_index));
    i = i+2;
end

crc = bitxor(crc,m);
crc_bin = dec2bin(crc,32);
crc = fliplr(crc_bin);
crc=logical((sscanf(crc,'%1d')).');

%Convert to hex, flip nibbles and invert bits within a nibble
aux=reshape(char(crc+'0'),8,[]);
crc=dec2hex(bin2dec(reshape(aux(end:-1:1,:),4,[]).')).';

% optional: force the WRONG CRC32 for the reference frame in order to
% compare with the reference in ANNEX G of 802.11-2007
%crc = logical([0 1 0 1 1 0 1 1 1 1 1 0 1 0 1 0 1 0 0 1 1 0 0 1 1 0 1 1 0 1 1 1]);
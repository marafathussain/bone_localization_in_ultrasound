%% ---- Mohammad Arafat Hussain (arafat@ece.ubc.ca) ----
function [RF,header] = RPread_arafat(filename)


fid = fopen(filename,'r');
header = fread(fid,19,'int32');
RF = zeros(header(4),header(3),floor(header(2)/2));
t = 1;
for i = 1:header(2)/2
    RF(:,:,t) = fread(fid, [header(4), header(3)],'int16');
    t = t+1;
end



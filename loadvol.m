%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loads the ultrasound volume data saved from the Propello application
%%
%% Inputs:  
%%     datapath - The path of the data to open
%%     volN   - The volume to retrieve
%%
%% Return:
%%     volume -  The volume data returned into a 3D array
%%     datatype -   Data type (0 = prescan b, 1  = postscan b, 2 = rf)
%%     numVols -    The number of volumes in the file
%%     fpV -        Frames per volume
%%     h -          Image height of the image (samples per line)
%%     w -          Image width (scan lines per frame)
%%     ss -         Sample size
%%
%% Ultrasonix Medical Corporation Nov 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [volume, datatype, numVols, fpV, h, w, ss]  = loadvol(datapath, volN)

if nargin == 1
    volN = 1;
end

fid= fopen(datapath, 'r');

if( fid == -1)
    error('Cannot open file');
end

% read the header info
volinfo = fread(fid, 7, 'int');

datatype = volinfo(1);  
numVols = volinfo(2);
fpV = volinfo(3);
w = volinfo(4);
h = volinfo(5);
ss = volinfo(6);
buffer = [];
volume = [];
    
if(datatype == 0 || datatype == 1)
  % read volume data
  volData = fread(fid, inf, 'uchar=>uchar');

   fclose(fid);

    if (numVols < volN)
        error('Error invalid volume number');
        volume = 0;
        return
    end

    buffer = volData((volN-1)*(h*w*fpV)+1 : volN*(h*w*fpV));

    % prescan b 
    if(datatype == 0)
        for i = 1 : fpV
            volume(:,:,i) = reshape(buffer((i-1)*h*w+1 : i*h*w), h, w);  
        end
    end
    
    % postscan b
    if(datatype == 1)
        for i = 1 : fpV
            temp = reshape(buffer((i-1)*h*w+1 : i*h*w), w, h);  
            volume(:,:,i) = imrotate(temp, -90);
        end
    end

% rf data
else
    other = fread(fid, (volN-1)*(h*w*fpV),'int16');
    
    for i = 1:fpV
        %the actual data...
        [temp,count] = fread(fid,h*w,'int16'); 
        volume(:,:,i) = int16(reshape(temp,h,w));
    end

    fclose(fid);
    
end
    



function Env=envelope(rf,n)
% ENVELOPE computes the envelope of an RF line/frame
%   envelope(rf) computes the envelope of the RF echo: 
%   rf, using the Hilbert transform algorithm. 
%   envelope(rf,n) uses n-point FFT to reduce computation 
%   time, if the echo length is not 2^m, m being an integer.
%
%   An advantage of the Hilbert algorithm is that, 
%   being an asynchronous detector, it does not need 
%   the center frequency.

% Author:  S. Kaisar Alam
% Email:   kalam@rrinyc.org
% Written: 06-99
% Revised: 03-22-05 (SKA)
% Version: 4.0
%
% New in this version: Now using MATLAB's hilbert function.
%
% Copyright © 1999 S. Kaisar Alam. All rights reserved.
% Questions & suggestions to S. Kaisar Alam.
%___________________________________________________________

if nargin < 1, error('Needs rf frame!'); end
if nargin > 2, error('Too many inputs! Maximum 2 necessary.'); end

nd = ndims(rf);

if nd > 3
   error('Number of dimension cannot be over 3') % CANNOT HANDLE 3-D DATA YET
elseif nd == 3 % 3-D DATA
   [row,colx,colz]=size(rf);
   if nargin == 1, n = row; end
   
   for k = 1:colz
      RF = rf(:,:,k);	% ONE SCAN PLANE
      env = abs(hilbert(RF,n)); % COMPUTE ENVELOPE 1 SCAN PLANE AT A TIME
      Env(:,:,k) = env(1:row,:);
   end   
else % 1-D & 2-D DATA
   [row,col]=size(rf);
   if nargin == 1, n = row; end
   Env = abs(hilbert(rf,n)); % COMPUTE ENVELOPE
end
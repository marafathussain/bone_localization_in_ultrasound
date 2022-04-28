function y = LogComp(x,option,fac)
%LOGCOMP Log compression for B-mode display
%  LOGCOMP(X,OPTION,FAC) performs log-compression according
%  to OPTION (output Y = 0, when input X = 0).
%  For OPTION = 'mu' (default), it performs a mu-law compression,
%  well-known in communications, using FAC as the mu value.
%  For OPTION = 'A', it performs an A-law compression, used in
%  communications in Europe, using FAC as the A value.
%  OPTION = 'E' returns max(X(:))*log10(1 + FAC/X/max(X(:)))/log10(1 + FAC).
%  OPTION = 'NONE' returns max(X(:))*log10(FAC + X)/log10(FAC + max(X(:)).
%  Values lower than zero are set to zero. FAC determines the dynamic
%  range of the output signal. Suggested values are 0.001, 0.01, 0.1,
%  and 1 (default: 0.05).
%  OPTION 'E' and 'NONE' are equivalent, with 'E' producing the same 
%  output for FAC 1000 times greater than for 'NONE'
%  X cannot be negative. X and Y have the same dynamic range.
%
%  See also ALAW,LOG10,MULAW.

% Author:  S. Kaisar Alam, Ph.D.
% Email:   kaisar.alam@ieee.org
% Written: 03-18-05
% Revised: 07-16-07 (SKA)
% Version: 1.2
%
% New in this version: changed default option to 'A' and removed error on
%          line 17
%
% Copyright © 2005 S. Kaisar Alam. All rights reserved.
% Questions & suggestions to S. Kaisar Alam.
%___________________________________________________________

if nargin == 0, error('No input signal!');
elseif nargin == 1, fac = 5; option = 'A';
elseif nargin == 2
   if ~ischar(option)
      fac = option; option = 'A';
   else
      if strcmpi(option,'mu')
         fac = 255;
      elseif strcmpi(option,'A')
         fac = 87.6;
      elseif strcmpi(option,'none')
         fac = 0.05;
      else
         error('Unknown log-compression option!');
      end
   end
end
if ~exist('option','var'), option = 'none'; end

if nargin < 3 % MU and A-LAW NEED THOSE PARAMETERS. NO DEFAULT.
   if any(x < 0), error('Input singal cannot be negative!'); end
end

if strcmpi(option,'mu')
   y = mulaw(x,fac); % FAC is MU
elseif strcmpi(option,'A')
   y = alaw(x,fac); % FAC is A
elseif strcmpi(option,'e') % Ernie's formula for straight log
   y = max(x(:))*log10(1 + fac*x/max(x(:)))/log10(1 + fac);
   % ERNIE USED x - x_MIN instead of x.
else % My formula for straight log
   y = max(x(:))*log10(1 + fac*x)/log10(1 + fac*max(x(:)));
   mask = y >= 0; y = y.*mask; % SETTING NEGATIVE VALUES TO ZERO, IF ANY.
end
% DYNAMIC RANGE OF X and Y ARE THE SAME
%% RUNME
% Computes and plots the frequency response of a plate equipped with tuned
% vibration absorbers (TVA). For more information, see
% http://modelreduction.org/index.php/Plate_with_tuned_vibration_absorbers.

% Copyright 2023 Quirin Aumann
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

clear; close all

% speficy evaluation type:
%   'rms' - evaluate root mean square z-displacement of all nodes on the
%       plate surface
%   'siso' - evaluate the z-displacement at the load location
type = 'rms';

if strcmp(type,'rms')
    load('plateTVA_n201900m1q28278')
    n_nodes = full(sum(sum(C)));
elseif strcmp(type,'siso')
    load('plateTVA_n201900m1q1')
else
    error('type must be "rms" or "siso"')
end

recompute = false;
if recompute == true
    res = zeros(1,length(s));
    for ii=1:length(s)
        fprintf('Frequency step %d, f=%.2f Hz ... ',ii,imag(s(ii))/2/pi)
        tic
        x = C * ((s(ii)^2 * M + s(ii) * E + K) \ B);
        if strcmp(type,'rms')
            res(ii) = sqrt((x'*x)/n_nodes);
        elseif strcmp(type,'siso')
            res(ii) = x;
        end
        fprintf('finished in %.2f s\n',toc)
    end
end

f = imag(s)/2/pi;
mag = 10*log10(abs(res)/1e-9);

figure('name','Transfer function')
plot(f,mag)
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

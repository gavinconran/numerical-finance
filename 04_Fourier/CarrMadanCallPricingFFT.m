% This is material illustrating the methods from the book
% Financial Modelling  - Theory, Implementation and Practice with Matlab
% source
% Wiley Finance Series
% ISBN 978-0-470-74489-5
%
% Date: 02.05.2012
%
% Authors:  Joerg Kienitz
%           Daniel Wetterau
%           Manuel Wittke
%
% Please send comments, suggestions, bugs, code etc. to
% kienitzwetterau_FinModelling@gmx.de
%
% (C) Joerg Kienitz, Daniel Wetterau, Manuel Wittke
% 
% Since this piece of code is distributed via the mathworks file-exchange
% it is covered by the BSD license 
%
% This code is being provided solely for information and general 
% illustrative purposes. The authors will not be responsible for the 
% consequences of reliance upon using the code or for numbers produced 
% from using the code. 



function [call_price_fft, delta_fft,gamma_fft] = CarrMadanCallPricingFFT(N,eta,model,S,K,T,r,d)

lnS = log(S);
lnK = log(K);

%optAlpha = optimalAlpha(model,lnS,lnK,T,r,d,varargin{:});

% sigma = varargin{:};
% d1 = (log(S/K)+(r-d+0.5*sigma^2)*T)/sigma/sqrt(T);
% optAlpha = -d1/sigma/sqrt(T);

optAlpha = 1.5;

DiscountFactor = exp(-r*T);


%-------------------------
%--- FFT Option Pricing --
%-------------------------
% from: Option Valuation Using the Fast Fourier Transform, 
%       Peter Carr, March 1999, pp 10-11
%-------------------------

% predefined parameters
FFT_N = N;% must be a power of two (2^14)
FFT_eta = eta; % spacing of psi integrand

% effective upper limit for integration (18)
% uplim = FFT_N * FFT_eta;

FFT_lambda = (2 * pi) / (FFT_N * FFT_eta); %spacing for log strike output (23)
FFT_b = (FFT_N * FFT_lambda) / 2; % (20)

uvec = 1:FFT_N;
%log strike levels ranging from lnS-b to lnS+b
ku = - FFT_b + FFT_lambda * (uvec - 1); %(19)

jvec = 1:FFT_N;
vj = (jvec-1) * FFT_eta;

%applying FFT
tmp = DiscountFactor * psi(model,vj,optAlpha,lnS,r,d,T) .* exp(1i * vj * (FFT_b)) * FFT_eta;
tmp = (tmp / 3) .* (3 + (-1).^jvec - ((jvec - 1) == 0) ); %applying simpson's rule
cpvec = real(exp(-optAlpha .* ku) .* fft(tmp) / pi); %call price vector resulting in equation 24
deltavec = real(exp(-optAlpha .* ku) .* fft((1i*vj+(optAlpha+1)).*tmp) / S / pi);
gammavec = real(exp(-optAlpha .* ku) .* fft((1i*vj+(optAlpha+1)).*((1i*vj+(optAlpha+1))-1).*tmp) / S /S / pi);

call_price_fft = real(interp1(ku,cpvec,lnK));
delta_fft = real(interp1(ku,deltavec,lnK));
gamma_fft = real(interp1(ku,gammavec,lnK));

end




%analytical formula for zhi in equation ( 6 ) of Madan's paper
function ret = psi(model,v,alpha,lnS,r,d,T)
  ret = exp((1i*v+(alpha+1))*lnS+feval(@CharacteristicFunctionLib, model, v - (alpha + 1) * 1i,T,r,d,model.params)) ./ (alpha.^2 + alpha - v.^2 + 1i * (2 * alpha + 1) .* v);
end


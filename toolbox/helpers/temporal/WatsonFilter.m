function [timp,filterParams] = WatsonFilter(filterParams,tSecs)
% Return the temporal impulse function suggested by Watson.
%
% Synopsis:
%    timp = WatsonFilter(filterParams,tSecs) 
%
% Description:
%    Return the temporal impulse response function as described in Watson's
%    1986 review chapter.
%
% Inputs:
%    filterParams           - Structure containing filter parameters.  Here
%                                     is an example that indicates the
%                                     fields, as well as the default values
%                                     if filterParams is passed as [].
%                                         filterParams.type = 'watson';            % filter type      
%                                         filterParams.tau = 6.25;                 % time constant in msec
%                                         filterParams.k = 1.33;                   % scaling factor for the time constant of the second filter (which equals k*tau)
%                                         filterParams.n1 = 9;                     % number of stages for the first filter
%                                         filterParams.n2 = 10;                    % number of stages for the second filter
%                                         filterParams.zeta = 1;                   % transience factor (0 = no adaptation/sustained; 1 = full adaptation/transient)
%                                         filterParams.xi = 22;                    % sensitivity factor (sensitivity factor or gain that scales the impulse response and amplitude response up or down in amplitude)
%   tSecs                   - Timebase in seconds for the impulse response, in seconds
%
% Outputs:
%   timp                        - The temporal impulse response
%   filterParams           - The filter parameter structure.  Useful for getting the default parameters
%
% See example in source code

% History:
%   03/11/25   dhb  Made this a callable function from code sent to me by Mengxin Wang

%{ 
% Examples:
    clear;
    filterParams.type = 'watson';            % filter type
    filterParams.tau = 6.25;                 % time constant
    filterParams.k = 1.33;                   % scaling factor for the time constant of the second filter (which equals k*tau)
    filterParams.n1 = 9;                     % number of stages for the first filter
    filterParams.n2 = 10;                    % number of stages for the second filter
    filterParams.zeta = 1;                   % transience factor (0 = no adaptation/sustained; 1 = full adaptation/transient)
    filterParams.xi = 22;                    % sensitivity factor (sensitivity factor or gain that scales the impulse response and amplitude response up or down in amplitude)

    tSecs = linspace(0,0.3,1000);
    timp = WatsonFilter(filterParams,tSecs);
    figure; plot(tSecs,timp);

    [~,defaultParams] = WatsonFilter([],[]);
%}

if (isempty(filterParams))
    filterParams.type = 'watson';            % filter type
    filterParams.tau = 6.25;                 % time constant
    filterParams.k = 1.33;                   % scaling factor for the time constant of the second filter (which equals k*tau)
    filterParams.n1 = 9;                     % number of stages for the first filter
    filterParams.n2 = 10;                    % number of stages for the second filter
    filterParams.zeta = 1;                   % transience factor (0 = no adaptation/sustained; 1 = full adaptation/transient)
    filterParams.xi = 22;                    % sensitivity factor (sensitivity factor or gain that scales the impulse response and amplitude response up or down in amplitude)
end
 
tmSecs = 1000*tSecs;
u = @(t)(0*(t<0)+1*(t>=0));
h1 = @(t)u(t)./real(filterParams.tau.*factorial(filterParams.n1-1)).*((t./filterParams.tau).^(filterParams.n1-1)).*exp(-t./filterParams.tau);
h2 = @(t)u(t)./real(filterParams.k.*filterParams.tau.*factorial(filterParams.n2-1)).*((t./(filterParams.k.*filterParams.tau)).^(filterParams.n2-1)).*exp(-t./(filterParams.k.*filterParams.tau));
h = @(t)filterParams.xi.*(h1(t)-filterParams.zeta.*h2(t));
timp = h(tmSecs);

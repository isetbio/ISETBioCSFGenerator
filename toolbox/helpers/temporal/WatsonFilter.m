

filterParams.integrationTime = 10; % integration time
filterParams.tau = 6.25; % time constant
filterParams.k = 1.33; % scaling factor for the time constant of the second filter (which equals k*tau)
filterParams.n1 = 9; % number of stages for the first filter
filterParams.n2 = 10; % number of stages for the second filter
filterParams.zeta = 1; % transience factor (0 = no adaptation/sustained; 1 = full adaptation/transient)
filterParams.xi = 22; % sensitivity factor (sensitivity factor or gain that scales the impulse response and amplitude response up or down in amplitude)
 
u = @(t)(0*(t<0)+1*(t>=0));
h1 = @(t)u(t)./real(filterParams.tau.*factorial(filterParams.n1-1)).*((t./filterParams.tau).^(filterParams.n1-1)).*exp(-t./filterParams.tau);
h2 = @(t)u(t)./real(filterParams.k.*filterParams.tau.*factorial(filterParams.n2-1)).*((t./(filterParams.k.*filterParams.tau)).^(filterParams.n2-1)).*exp(-t./(filterParams.k.*filterParams.tau));
h = @(t)filterParams.xi.*(h1(t)-filterParams.zeta.*h2(t));

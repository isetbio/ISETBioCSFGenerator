function plotLMthresholds()
    maxThreshold = 0.0401;
    theScriptName = 't_isoResponseLMplaneEllipses';
    plotCones_vs_mRGCs(maxThreshold, theScriptName);
end

function plotCones_vs_mRGCs(maxVisualizedThreshold, theScriptName)
    mRGCOutputSignalType = 'cones';             % {'cones', 'mRGCs'}
    opticsType = 'loadComputeReadyRGCMosaic';   % {'loadComputeReadyRGCMosaic', 'adaptiveOptics6MM'} 
    stimulusSpatialFrequency = 0;
    mosaicEccDegs = [0 0];
    mosaicSizeDegs = [1 1];
    orientationDegs = 90;
    presentationMode = 'counterphasemodulated';

    useMetaContrast = true;
    useConeContrast = true;
    useFixationalEMs = false;
    whichNoiseFreeNre = 'mRGCMosaic';
    whichNoisyInstanceNre = 'Gaussian';
    whichClassifierEngine = 'rceTemplateDistance';

    [figureFileBaseDir, resultsFileBaseDir] = setupFigureDirectory(theScriptName, ...
        useMetaContrast,useConeContrast,useFixationalEMs,...
        whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType);

    thresholdsDataFileName = ...
        sprintf('%sCSF_SF_%2.2fCPD_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusSpatialFrequency, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);

    load(fullfile(resultsFileBaseDir,thresholdsDataFileName), 'examinedDirectionsOnLMplane', 'thresholdContrasts');
    coneThresholds = thresholdContrasts;

    mRGCOutputSignalType = 'mRGCs';
    [~, resultsFileBaseDir] = setupFigureDirectory(theScriptName, ...
        useMetaContrast,useConeContrast,useFixationalEMs,...
        whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType);
    thresholdsDataFileName = ...
        sprintf('%sCSF_SF_%2.2fCPD_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusSpatialFrequency, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
    fullfile(resultsFileBaseDir,thresholdsDataFileName)
    load(fullfile(resultsFileBaseDir,thresholdsDataFileName), 'examinedDirectionsOnLMplane', 'thresholdContrasts');
    mRGCThresholds = thresholdContrasts;


    visualizeStimuliAndThresholdsOnLMPlane(...
            examinedDirectionsOnLMplane, ...
            coneThresholds, mRGCThresholds, ...
            'cones (L+M+S)', 'mRGCs (L+M+S)', ...
            maxVisualizedThreshold,  ...
            fullfile(figureFileBaseDir,'mRGCs_vs_cones.pdf'));

 end


function  visualizeStimuliAndThresholdsOnLMPlane(...
            examinedDirectionsOnLMplane, ...
            thresholds1, thresholds2, ...
            legend1, legend2, ...
            maxVisualizedThreshold,  ...
            pdfFileName)

    hFig = figure(2); clf;
    set(hFig, 'Color', [1 1 1]);
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard tall figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    % Plot the axes
    plot(theAxes{1,1},  maxVisualizedThreshold*[-1 1], [0 0], 'k-', 'LineWidth', 1.0);
    hold(theAxes{1,1}, 'on');
    plot(theAxes{1,1},  [0 0], maxVisualizedThreshold*[-1 1], 'k-', 'LineWidth', 1.0);

 
    x1 = cosd(examinedDirectionsOnLMplane) .* thresholds1;
    y1 = sind(examinedDirectionsOnLMplane) .* thresholds1;

    x2 = cosd(examinedDirectionsOnLMplane) .* thresholds2;
    y2 = sind(examinedDirectionsOnLMplane) .* thresholds2;
    
    if (numel(examinedDirectionsOnLMplane)>6)
        [z, a, b, rotationRadians] = fitEllipseToXYpoints(...
            [x1(:) y1(:)]', ...
            'nonLinear', false);
        %[z, a, b, rotationRadians] = fitellipse(thresholdConeContrasts, 'linear');
    
        % Plot the fitted ellipse
        % form the parameter vector
        npts = 100;
        t = linspace(0, 2*pi, npts);
    
        % Rotation matrix
        Q = [cos(rotationRadians), -sin(rotationRadians); sin(rotationRadians) cos(rotationRadians)];
        % Ellipse points
        X = Q * [a * cos(t); b * sin(t)] + repmat(z, 1, npts);

        % PLot the ellipse points
        plot(theAxes{1,1}, X(1,:), X(2,:), 'k-', 'LineWidth', 6.0);
        plot(theAxes{1,1}, X(1,:), X(2,:), 'c-', 'Color', [1 0.5 0.7], 'LineWidth',2.0);


        [z, a, b, rotationRadians] = fitEllipseToXYpoints(...
            [x2(:) y2(:)]', ...
            'nonLinear', false);
        %[z, a, b, rotationRadians] = fitellipse(thresholdConeContrasts, 'linear');
    
        % Plot the fitted ellipse
        % form the parameter vector
        npts = 100;
        t = linspace(0, 2*pi, npts);
    
        % Rotation matrix
        Q = [cos(rotationRadians), -sin(rotationRadians); sin(rotationRadians) cos(rotationRadians)];
        % Ellipse points
        X = Q * [a * cos(t); b * sin(t)] + repmat(z, 1, npts);

        % PLot the ellipse points
        plot(theAxes{1,1}, X(1,:), X(2,:), 'k-', 'LineWidth', 6.0);
        plot(theAxes{1,1}, X(1,:), X(2,:), 'c-', 'Color', [0.5 1 0.85], 'LineWidth',2.0);

    end

     % LM plane thresholds-1
    
    p1 = scatter(theAxes{1,1}, x1,y1, 300, ...
        'MarkerFaceColor', [[1 0.5 0.7]], 'MarkerEdgeColor', [1 0.5 0.7]*0.5, ...
        'MarkerFaceAlpha', 0.6, 'LineWidth', 1.5);
    hold(theAxes{1,1}, 'on');


    % LM plane thresholds-2
    
    p2 = scatter(theAxes{1,1}, x2,y2, 300, ...
        'MarkerFaceColor', [0.5 1 0.85], 'MarkerEdgeColor', [0.5 1 0.85]*0.5, ...
        'MarkerFaceAlpha', 0.6, 'LineWidth', 1.5);

   

    xLims = maxVisualizedThreshold*[-1 1];
    yLims = maxVisualizedThreshold*[-1 1];
    hold(theAxes{1,1}, 'off');
    set(theAxes{1,1}, 'XLim', maxVisualizedThreshold*[-1 1], 'YLim', maxVisualizedThreshold*[-1 1]);
    axis(theAxes{1,1}, 'square');
    box(theAxes{1,1}, 'off');
    set(theAxes{1,1}, 'Color', 'none', 'Box', 'off', 'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1]);
    set(theAxes{1,1}, 'XTick', -.1:0.02:0.1, 'XTickLabel', {'-.1', '-.08', '-.06', '-.04', '-.02', '0', '+.02', '+.04', '+.06', '+.08', '+.1'})
    set(theAxes{1,1}, 'YTick', -.1:0.02:0.1, 'YTickLabel', {'-.1', '-.08', '-.06', '-.04', '-.02', '0', '+.02', '+.04', '+.06', '+.08', '+.1'})
    
    %PublicationReadyPlotLib.offsetAxes(theAxes{1,1},ff, xLims, yLims);
    PublicationReadyPlotLib.labelAxes(theAxes{1,1},ff, 'threshold contrast (L-cone)', 'threshold contrast (M-cone)');
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

    %set(theThresholdAxes, 'FontSize', 24);

end

function [figureFileBaseDir, resultsFileBaseDir] = setupFigureDirectory(theScriptName, ...
    useMetaContrast,useConeContrast,useFixationalEMs,...
    whichNoiseFreeNre,whichNoisyInstanceNre,...
    whichClassifierEngine,mRGCOutputSignalType)

    % Make sure local/figures directory exists so we can write out our figures in peace
    projectBaseDir = ISETBioCSFGeneratorRootPath;
    if (~exist(fullfile(projectBaseDir,'local', theScriptName, 'figures'),'dir'))
        mkdir(fullfile(projectBaseDir,'local', theScriptName, 'figures'));
        fprintf('Generated figure directory at %s\n', fullfile(projectBaseDir,'local', theScriptName, 'figures'))
    end

    if (~exist(fullfile(projectBaseDir,'local', theScriptName, 'results'),'dir'))
        mkdir(fullfile(projectBaseDir,'local', theScriptName, 'results'));
        fprintf('Generated results directory at %s\n', fullfile(projectBaseDir,'local', theScriptName, 'figures'))
    end


    figureFileBaseDir = fullfile(projectBaseDir,'local',theScriptName,'figures', ...
        sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', theScriptName, ...
        useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType));

    if (~exist(figureFileBaseDir, 'dir'))
        mkdir(figureFileBaseDir);
        fprintf('Generated figure sub-directory at %s\n', figureFileBaseDir);
    end

    resultsFileBaseDir = fullfile(projectBaseDir,'local',theScriptName,'results', ...
        sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', theScriptName, ...
        useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType));
    
    if (~exist(resultsFileBaseDir, 'dir'))
        mkdir(resultsFileBaseDir);
        fprintf('Generated figure sub-directory at %s\n', resultsFileBaseDir);
    end

end


function [z, a, b, alpha] = fitEllipseToXYpoints(xyPoints, varargin)

	% Parse input
    p = inputParser;
    p.addParameter('maxIterations', 200, @isnumeric);
    p.addParameter('tolerance', 1e-5, @isnumeric);
    p.addParameter('constraint', 'bookstein', @(x)(ismember(x, {'bookstein', 'trace'})));
    p.addParameter('nonLinear', true, @islogical);
    p.parse(varargin{:});

	rows = size(xyPoints,1);
	cols = size(xyPoints,2);
	assert(rows == 2, 'xyPoints must be a 2 x N matrix');
	assert(cols >=6, 'xyPoints must have at least 6 points');

	fitParams = struct(...
		'nonLinear', p.Results.nonLinear, ...
		'constraint', p.Results.constraint, ...
		'maxIterations', p.Results.maxIterations, ...
		'tolerance', p.Results.tolerance);

	% Remove centroid
	centroid = mean(xyPoints, 2);
	xyPoints = bsxfun(@minus, xyPoints, centroid);

	% Obtain a linear estimate
	switch fitParams.constraint
    	case 'bookstein'
        	[z, a, b, alpha] = fitbookstein(xyPoints);
    	case 'trace'
       		[z, a, b, alpha] = fitggk(xyPoints);
	end % switch


	if (fitParams.nonLinear)
		% Initial conditions
	    z0     = z;
	    a0     = a;
	    b0     = b;
	    alpha0 = alpha;

	    % Fit
    	[z, a, b, alpha, converged, isCircle] = fitNonLinear(xyPoints, z0, a0, b0, alpha0, params);

    	% Return linear estimate if GN doesn't converge or if the data points fall on a circle
	    if (~converged) || (isCircle)
	        fprintf('*** FailureToConverge: Gauss-Newton did not converge, returning linear estimate.\n');
	        z = z0;
	        a = a0;
	        b = b0;
	        alpha = alpha0;
	    end
	end % if (fitParams.nonLinear)

	% Add the centroid back on
	z = z + centroid;
end


function [z, a, b, alpha] = fitbookstein(x)
	%FITBOOKSTEIN   Linear ellipse fit using bookstein constraint
	%   lambda_1^2 + lambda_2^2 = 1, where lambda_i are the eigenvalues of A

	% Convenience variables
	m  = size(x, 2);
	x1 = x(1, :)';
	x2 = x(2, :)';

	% Define the coefficient matrix B, such that we solve the system
	% B *[v; w] = 0, with the constraint norm(w) == 1
	B = [x1, x2, ones(m, 1), x1.^2, sqrt(2) * x1 .* x2, x2.^2];

	% To enforce the constraint, we need to take the QR decomposition
	[Q, R] = qr(B);

	% Decompose R into blocks
	R11 = R(1:3, 1:3);
	R12 = R(1:3, 4:6);
	R22 = R(4:6, 4:6);

	% Solve R22 * w = 0 subject to norm(w) == 1
	[U, S, V] = svd(R22);
	w = V(:, 3);

	% Solve for the remaining variables
	v = -R11 \ R12 * w;

	% Fill in the quadratic form
	A        = zeros(2);
	A(1)     = w(1);
	A([2 3]) = 1 / sqrt(2) * w(2);
	A(4)     = w(3);
	bv       = v(1:2);
	c        = v(3);

	% Find the parameters
	[z, a, b, alpha] = conic2parametric(A, bv, c);
end

function [z, a, b, alpha] = fitggk(x)
	% Linear least squares with the Euclidean-invariant constraint Trace(A) = 1
	% Convenience variables
	m  = size(x, 2);
	x1 = x(1, :)';
	x2 = x(2, :)';

	% Coefficient matrix
	B = [2 * x1 .* x2, x2.^2 - x1.^2, x1, x2, ones(m, 1)];
	v = B \ -x1.^2;

	% For clarity, fill in the quadratic form variables
	A        = zeros(2);
	A(1,1)   = 1 - v(2);
	A([2 3]) = v(1);
	A(2,2)   = v(2);
	bv       = v(3:4);
	c        = v(5);

	% find parameters
	[z, a, b, alpha] = conic2parametric(A, bv, c);
end


function [z, a, b, alpha, converged, isCircle] = fitNonLinear(x, z0, a0, b0, alpha0, params)
	% Gauss-Newton least squares ellipse fit minimising geometric distance 

	% Get initial rotation matrix
	Q0 = [cos(alpha0), -sin(alpha0); sin(alpha0) cos(alpha0)];
	m = size(x, 2);

	% Get initial phase estimates
	phi0 = angle( [1 i] * Q0' * (x - repmat(z0, 1, m)) )';
	u = [phi0; alpha0; a0; b0; z0];

	% Iterate using Gauss Newton
	converged = false;

	for nIts = 1:params.maxIterations
	    % Find the function and Jacobian
	    [f, J, isCircle] = computeJacobian(u);
    
    	if (isCircle)
    		fprintf('Ellipse is near-circular - nonlinear fit may not succeed\n.')
    	end

	    % Solve for the step and update u
	    h = -J \ f;
	    u = u + h;
    
	    % Check for convergence
	    delta = norm(h, inf) / norm(u, inf);
	    if delta < params.tolerance
	        converged = true;
	        break
	    end
	end

	alpha = u(end-4);
	a = u(end-3);
	b = u(end-2);
	z = u(end-1:end);

	% ---- Nested function ---
	function [f, J, isCircle] = computeJacobian(u)
        % Define the system of nonlinear equations and Jacobian. 

        % Tolerance for whether it is a circle
        circTol = 1e-5;
        
        % Unpack parameters from u
        phi   = u(1:end-5);
        alpha = u(end-4);
        a     = u(end-3);
        b     = u(end-2);
        z     = u(end-1:end);
        
        % If it is a circle, the Jacobian will be singular, and the
        % Gauss-Newton step won't work. 
        %TODO: This can be fixed by switching to a Levenberg-Marquardt
        %solver
        if (abs(a - b) / (a + b) < circTol)
            isCircle = true;
        else
        	isCircle = false;
        end

        % Convenience trig variables
        c = cos(phi);
        s = sin(phi);
        ca = cos(alpha);
        sa = sin(alpha);
        
        % Rotation matrices
        Q    = [ca, -sa; sa, ca];
        Qdot = [-sa, -ca; ca, -sa];

        % Preallocate function and Jacobian variables
        f = zeros(2 * m, 1);
        J = zeros(2 * m, m + 5);
        for i = 1:m
            rows = (2*i-1):(2*i);
            % Equation system - vector difference between point on ellipse
            % and data point
            f((2*i-1):(2*i)) = x(:, i) - z - Q * [a * cos(phi(i)); b * sin(phi(i))];
            
            % Jacobian
            J(rows, i) = -Q * [-a * s(i); b * c(i)];
            J(rows, (end-4:end)) = ...
                [-Qdot*[a*c(i); b*s(i)], -Q*[c(i); 0], -Q*[0; s(i)], [-1 0; 0 -1]];
        end
    end % ---- Nested function ---
end

function [z, a, b, alpha] = conic2parametric(A, bv, c)
	% Diagonalise A - find Q, D such at A = Q' * D * Q
	[Q, D] = eig(A);
	Q = Q';

	% If the determinant < 0, it's not an ellipse
	if prod(diag(D)) <= 0 
	    error('NotEllipse', 'Linear fit did not produce an ellipse');
	end

	% We have b_h' = 2 * t' * A + b'
	t = -0.5 * (A \ bv);

	c_h = t' * A * t + bv' * t + c;

	z = t;
	a = sqrt(-c_h / D(1,1));
	b = sqrt(-c_h / D(2,2));
	alpha = atan2(Q(1,2), Q(1,1));
end % conic2parametric

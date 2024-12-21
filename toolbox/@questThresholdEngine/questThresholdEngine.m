classdef questThresholdEngine < contrastThresholdEngine
    %questThresholdEngine  Adaptive contrast threshold estimation procedure
    % based on the QUEST+ routine. The class maintain an array of questData
    % ojbect and loop through them sequentially. A stop criterion can be
    % triggered if the standard error among the set of questData object
    % drops below a pre-defined a threshold and we have the minimal number of
    % trials, or maximum number of trials is reached.
    %
    % Usage:
    %   See t_spatialCSF.m
    %   Also see base class contrastThresholdEngine
    %
    %   Note. The interpretation of the stimulus units is handled by the
    %   psychometric function associated with the object. The default is
    %   qpPFWeibull, which itself works in dB units of contrast. If you
    %   want log10 units, specify the PF as @qpPFWeibullLog using the
    %   'qpPF' key/value pair. Similarly for linear contrast units, use
    %   @qpPFStandardWeibull.
    %
    %   Other PFs could be added, but they need to be in the form expected
    %   by mQUESTPlus, and need to have an inverse PF associated with them.
    %   Add to the list in the constructor to add additional PFs; follow
    %   the template there to see how to check for the PF and add its
    %   associated inverse.
    %
    % questThresholdEngine Properties:
    %   estimators     - The array of questData object.
    %   numEstimator   - The number questData object to maintain
    %   stopCriterion  - Stop criterion for deciding if we have enough
    %                    trials. See key/value pair description below.
    %   slopeRange     - Range of slopes for psychometric curve
    %   guessRate      - Range of guess rate for psychometric curve
    %   lapseRate      - Range of lapseRate rate for psychometric curve
    %   estIdx         - Current questData object being used
    %
    %   Also see base class contrastThresholdEngine for additional info,
    %   including info about key/value pair 'estDomain'.
    %
    % questThresholdEngine Methods:
    %   thresholdEstimate    - Current running estimate of threshold and
    %                          its standard error across questData object
    %   combineData          - Return all stimulus - response data we have
    %                          recorded so far
    %
    %   thresholdMLE         - Run a MLE on all data we have recorded so
    %                          far, plot data and the psychometric curve
    %
    %   Also see base class QuestThresholdEstimator
    %
    % Inputs:
    %   None.
    %
    % Outputs:
    %   questThresholdEngine Object.
    %
    % Optional key/value pairs.  Also inherits those in contrastThesholdEngine.
    %   'minTrial'        - Int. The minimal number of trials to run
    %
    %   'maxTrial'        - Int. The maximum number of trials to run
    %
    %   'estDomain'       - Array. The domain (range) for which contrast is
    %                       being estimated. For examples, linear contrast
    %                       [0:0.01:1] or log contrast [-10:0.1:0]. User
    %                       should handle the conversion externally
    %
    %   'stopCriterion'   - Either empty, a number, or a function handle.
    %                       If empty, then no stop criterion. If a number,
    %                       stops when SE across estimators is less than
    %                       the number.  If a function handle, it should
    %                       take two arguments, the current threshold
    %                       estimate, and the S.E. of the estimate among N
    %                       quest objects, and return a boolean variable if
    %                       more trials are required. If numEstimator is 1,
    %                       then the S.E. is always zero so the routine
    %                       stops when the minTrials value is reached. Same
    %                       when this is empty.  Indeed, setting to empty
    %                       means that S.E. across estimators is not used
    %                       as a stopping criterion. See
    %                       t_spatialCSF for elaboration on stopping.
    %
    %   'slopeRange'       - Array. An array of all possible slope for the
    %                        psychometric curve
    %
    %   'guessRate'        - Array. An array of all possible guess rate for
    %                        the psychometric curve. Default 0.5;
    %
    %   'lapseRate'        - Array. An array of all possible lapse rate for
    %                        the psychometric curve. Default 0.
    %
    %   'validation'       - Boolean. If set to true run the entire
    %                        psychometric curve
    %
    %   'blocked'          - Run trials blocked by contrast.  Set this to true
    %                        if a block of trials at same contrast will be run.
    %                        Often true in simulations, rarely true in experiments.
    %                        Default 0.
    %
    %   'nRepeat'          - Double. Number of trials per contrast level
    %                        when running in validation mode.  This can be
    %                        confusing, because it feels a lot the same as
    %                        the nTest parameter that we use for the
    %                        classifier engine.
    %
    %   'nOutcomes'    - Double. Number of stimulus alternatives per
    %                        trial. Default 2.
    %   'qpPF'         - Psychometric function as expected by QuestPlus.
    %                    Default: @qpPFWeibull.  Other options are
    %                    qpPFWeibullLog and qpPFStandardWeibull.
    %
    %    See also t_spatialCSF, contrastThresholdEngine, mQUESTPlus

    % History:
    %   04/22/23  dhb  Change sign of what happens from true to false for
    %                  empty stop criterion, to match what the documentation above says
    %                  will happen. 
    %             dhb  Change comparison of se to numeric to specified stop
    %                  criterion from <= to <, so that we don't throw a stop
    %                  when the criterion is zero se and there is only one
    %                  estimator.
    
   
    % Class properties
    properties %(Access = private)
        
        estimators;
        numEstimator;
        stopCriterion;
        validation;
        blocked;
        
        slopeRange;
        guessRate;
        lapseRate;
        
        estIdx;
        nRepeat;

        nOutcomes;

        qpPF;
        qpPFInv;
        
        validationTrialContrasts;
    end
    
    methods
        
        % Constructor method
        function this = questThresholdEngine(varargin)
            
            this@contrastThresholdEngine(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addParameter('numEstimator', 1, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('stopCriterion', 0.05, @(x) (isa(x, 'function_handle') || isnumeric(x)));            
            p.addParameter('slopeRange', 0.1 : 0.5 : 50);
            p.addParameter('guessRate', 0.5);
            p.addParameter('nOutcomes', 2);
            p.addParameter('lapseRate', 0.0);
            p.addParameter('validation', false, @(x)(islogical(x) && numel(x) == 1));
            p.addParameter('blocked', false, @(x)(islogical(x) && numel(x) == 1));
            p.addParameter('nRepeat', 64, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('qpPF',@qpPFWeibull,@(x) isa(x,'function_handle'));
                        
            parse(p, varargin{:});
            this.numEstimator  = p.Results.numEstimator;
            this.slopeRange = p.Results.slopeRange;
            this.guessRate = p.Results.guessRate;
            this.lapseRate = p.Results.lapseRate;
            this.validation = p.Results.validation;
            this.blocked = p.Results.blocked;
            this.nRepeat = p.Results.nRepeat;
            this.nOutcomes = p.Results.nOutcomes;
            this.qpPF = p.Results.qpPF;
            
            stopCriterion = p.Results.stopCriterion;
            if isempty(stopCriterion)
                this.stopCriterion = @(threshold, se) false;  
            elseif isnumeric(stopCriterion)
                this.stopCriterion = @(threshold, se) se < stopCriterion;
            elseif isa(stopCriterion, 'function_handle')
                this.stopCriterion = stopCriterion;
            else
                error('Input argument stopCriterion is an invalid type')
            end

            % Make sure we know about PF and set up its inverse.  Can add
            % more PFs here as desired.
            if (isequal(this.qpPF,@qpPFWeibull))
                this.qpPFInv = @qpPFWeibullInv;
            elseif (isequal(this.qpPF,@qpPFWeibullLog))
                this.qpPFInv = @qpPFWeibullLogInv;
            elseif (isequal(this.qpPF,@qpPFStandardWeibull))
                this.qpPFInv = @qpPFStandardWeibullInv;
            else
                error('Unknown qpPF specified');
            end
            
            % Initialize QUEST+ objects specified by 'numEstimator'
            this.estimators = cell(this.numEstimator, 1);
            for idx = 1 : this.numEstimator
                this.estimators{idx} = ...
                    qpInitialize('stimParamsDomainList', {this.estDomain}, ...
                    'psiParamsDomainList',  {this.estDomain, this.slopeRange, this.guessRate, this.lapseRate}, ...
                    'nOutcomes',this.nOutcomes,'qpPF',this.qpPF);
            end
            
            % Set the current estimator to #1
            this.estIdx = 1;            
            this.nTrial = 0;
            this.nextFlag = true;
            this.testCrst = qpQuery(this.estimators{this.estIdx});
            
            if this.validation
                if (this.numEstimator ~= 1)
                    error('Only run one estimator in validation mode');
                end

                % Pre-randomize the contrasts in blocked fashion in this
                % mode.
                % 
                % If trials are blocked, just set to blank as we use
                % estDomain directly in order, in this case.  Also,
                % override testCrst obtained above to start the party off
                % right.
                if (this.blocked)
                    this.validationTrialContrasts = [];
                    this.testCrst = this.estDomain(1);
                else
                    this.validationTrialContrasts = [];
                    for tt = 1:this.nRepeat
                        this.validationTrialContrasts = [this.validationTrialContrasts...
                            this.estDomain(randperm(length(this.estDomain)))];
                    end
                    this.testCrst = this.validationTrialContrasts(1);
                end
            end
        end
        
        % Running estimate of parameters of the psychometric curve
        estimates = parameterEstimate(this)
        
        % Running estimate of threshold and its standard error
        [threshold, stderr] = thresholdEstimate(this)
        
        % Combine data recorded so far from all QUEST object
        % Return a list of stimulus - response pair, and a struct of both
        [stimVec, responseVec, structVec] = combineData(this)
        
        % Record a set of trials of the experiment
        % Return next next query contrast, an indicator for new trial
        [nextCrst, nextFlag] = multiTrial(this, stimVec, responseVec)

        % Record a set of trials of the experiment
        % Return next next query contrast, an indicator for new trial
        % This version for Quest when nRepeat > 1.
        [nextCrst, nextFlag] = multiTrialQuestBlocked(this, stimVec, responseVec)
        
        % Record one trial of the experiment
        % Return next next query contrast, an indicator for new trial
        [nextCrst, nextFlag] = singleTrial(this, stim, response)

        % Record one trial of the experiment
        % Return next next query contrast, an indicator for new trial
        % This version is for cases where a block of trials (nRepeat > 1) at a single
        % contrast is run, in validation mode.
        [nextCrst, nextFlag] = singleTrialValidationBlocked(this, stim, response)
        
        % Run MLE estimate of psychometric curve parameter on combined data
        [threshold, para, psychometricDataOut] = thresholdMLE(this, varargin);
        
        % Plot data and MLE fit
        %
        % Inputs:
        %    this                   - The questThresholdEngine object.
        %    baseMarkerSize         - Base marker size.
        %    
        % Optional key/value pairs.
        %    'para'                 - Psychometric function parameters.  If passed,
        %                             these are used and the fit is skipped. Must
        %                             be matched to what the object's PF expects.
        %                             Default empty, so that fitting is done.
        %    'showPlot'             - Show a plot of data and fit. Default false.
        %    'newFigure'            - Create a new figure window for the plot. 
        %                             Otherwise the plot goes into the current
        %                             figure. Default false.
        function plotMLE(this, baseMarkerSize, varargin)
            p = inputParser;
            p.addParameter('newFigure', false);
            p.addParameter('para', []);
            parse(p, varargin{:});
            this.thresholdMLE('showPlot', true, 'newFigure', p.Results.newFigure, ...
                'para', p.Results.para, 'pointSize', baseMarkerSize);
        end
    end   
end



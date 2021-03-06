classdef questThresholdEngine < contrastThresholdEngine
    %questThresholdEngine  Adaptive contrast threshold estimation procedure
    % based on the QUEST+ routine. The class maintain an array of questData
    % ojbect and loop through them sequentially. A stop criterion can be
    % triggered if the standard error among the set of questData object
    % drops below a pre-defined a threshold and we have the minimal number of
    % trials, or maximum number of trials is reached.
    %
    % Usage:
    %   See t_thresholdEngine.m
    %   Also see base class QuestThresholdEstimator
    %
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
    %   Also see base class QuestThresholdEstimator
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
    %                       as a a stopping criterion. See
    %                       t_thresholdEngine for elaboration on stopping.
    %
    %   'slopeRange'       - Array. An array of all possible slope for the
    %                        psychometric curve
    %
    %   'guessRate'        - Array. An array of all possible guess rate for
    %                        the psychometric curve
    %
    %   'lapseRate'        - Array. An array of all possible lapse rate for
    %                        the psychometric curve
    %
    %   'validation'       - Boolean. If set to True run the entire
    %                        psychometric curve
    %
    %   'nRepeat'          - Double. Number of trials per contrast level
    %                        when running the validation
    %
    %    See also t_thresholdEngine, contrastThresholdEngine
    
   
    % Class properties
    properties %(Access = private)
        
        estimators;
        numEstimator;
        stopCriterion;
        validation;
        
        slopeRange;
        guessRate;
        lapseRate;
        
        estIdx;
        nRepeat;
        
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
            p.addParameter('lapseRate', 0.0);
            p.addParameter('validation', false, @(x)(islogical(x) && numel(x) == 1));
            p.addParameter('nRepeat', 64, @(x)(isnumeric(x) && numel(x) == 1));
                        
            parse(p, varargin{:});
            this.numEstimator  = p.Results.numEstimator;
            this.slopeRange = p.Results.slopeRange;
            this.guessRate = p.Results.guessRate;
            this.lapseRate = p.Results.lapseRate;
            this.validation = p.Results.validation;
            this.nRepeat = p.Results.nRepeat;
            
            stopCriterion = p.Results.stopCriterion;
            if isempty(stopCriterion)
                this.stopCriterion = @(threshold, se) true;  
            elseif isnumeric(stopCriterion)
                this.stopCriterion = @(threshold, se) se <= stopCriterion;
            elseif isa(stopCriterion, 'function_handle')
                this.stopCriterion = stopCriterion;
            else
                error('Input argument stopCriterion is an invalid type')
            end
            
            % Initialize QUEST+ objects specified by 'numEstimator'
            this.estimators = cell(this.numEstimator, 1);
            for idx = 1 : this.numEstimator
                this.estimators{idx} = ...
                    qpInitialize('stimParamsDomainList', {this.estDomain}, ...
                    'psiParamsDomainList',  {this.estDomain, this.slopeRange, this.guessRate, this.lapseRate});
            end
            
            % Set the current estimator to #1
            this.estIdx = 1;            
            this.nTrial = 0;
            this.nextFlag = true;
            this.testCrst = qpQuery(this.estimators{this.estIdx});
            
            if this.validation
                this.testCrst = this.estDomain(1);
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
        
        % Record one trial of the experiment
        % Return next next query contrast, an indicator for new trial
        [nextCrst, nextFlag] = singleTrial(this, stim, response)
        
        % Run MLE estimate of psychometric curve parameter on combined data
        [threshold, para, psychometricDataOut] = thresholdMLE(this, varargin);
        
        % Plot data and MLE fit
        function plotMLE(this, markerSize)
            this.thresholdMLE('showPlot', true, 'newFigure', false, 'pointSize', markerSize);
        end
    end
    
end



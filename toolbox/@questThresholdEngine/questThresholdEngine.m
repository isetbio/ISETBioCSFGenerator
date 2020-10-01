classdef questThresholdEngine < contrastThresholdEngine
    %QuestThresholdEstimator  Adaptive contrast threshold estimation procedure
    % based on the QUEST+ routine. The class maintain an array of questData
    % ojbect and loop through them sequentially. A stop criterion can be
    % triggered if the standard error among the set of questData object
    % drops below a pre-defined a threshold and we have the minimal number of
    % trials, or maximum number of trials is reached.
    %
    % Usage:
    %   See t_questAdaptiveEstimator.m
    %   Also see base class QuestThresholdEstimator
    %
    %
    % ContrastThresholdEstimator Properties:
    %   estimators     - The array of questData object.
    %   numEstimator   - The number questData object to maintain
    %   stopCriterion  - Stop criterion for deciding if we have enough trials
    %   slopeRange     - Range of slopes for psychometric curve
    %   guessRate      - Range of guess rate for psychometric curve
    %   lapseRate      - Range of lapseRate rate for psychometric curve
    %   estIdx         - Current questData object being used
    %
    %   Also see base class QuestThresholdEstimator
    %
    % QuestThresholdEstimator Methods:
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
    %   QuestThresholdEstimator Object.
    %
    % Optional key/value pairs:
    %   'numEstimator'    - Int. Number of questData object to run
    %
    %   'stopCriterion'   - Double. Should be between [0, 1]. Stop the
    %                       procedure if the standard error drops below
    %                       stopCriterion * thresholdEstimate
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
    %    Also see base class QuestThresholdEstimator
    
    
    % Class properties
    properties %(Access = private)
        
        estimators;
        numEstimator;
        stopCriterion;
        
        slopeRange;
        guessRate;
        lapseRate;
        
        estIdx;
        
    end
    
    methods
        
        % Constructor method
        function this = questThresholdEngine(varargin)
            
            this@contrastThresholdEngine(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addParameter('numEstimator', 1, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('stopCriterion', 0.05, @(x)(isnumeric(x) && numel(x) == 1));
            p.addParameter('slopeRange', 0.1 : 0.5 : 50);
            p.addParameter('guessRate', 0.5);
            p.addParameter('lapseRate', 0.0);
            
            
            parse(p, varargin{:});
            this.numEstimator  = p.Results.numEstimator;
            this.stopCriterion = p.Results.stopCriterion;
            this.slopeRange = p.Results.slopeRange;
            this.guessRate = p.Results.guessRate;
            this.lapseRate = p.Results.lapseRate;
            
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
            
        end
        
        % Running estimate of parameters of the psychometric curve 
        parameterEstimate(this)
        
        % Running estimate of threshold and its standard error
        thresholdEstimate(this)
        
        % Combine data recorded so far from all QUEST object
        combineData(this)
        
        % Record a set of trials of the experiment
        multiTrial(this, stimVec, responseVec)        
        
        % Record one trial of the experiment
        singleTrial(this, stim, response)
        
        % Run MLE estimate of psychometric curve parameter on combined data        
        thresholdMLE(this, varargin);
        
    end
    
end



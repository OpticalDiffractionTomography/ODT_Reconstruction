function cfg = pipelineConfig(pipelineVersion)
% pipelineConfig  Processing options for a named pipeline version.
%   cfg = pipelineConfig('default')  cluster pipeline
%   cfg = pipelineConfig('alice')    reproduces Alice's PC script
%                                    (backup/field_retrieval_Tomogram_reconstruction_Alice.mat)
%
% Selected via --version/-v in main.sh and tomo_process, which inject
% pipelineVersion into the MATLAB workspace before each stage runs.
    if nargin < 1 || isempty(pipelineVersion)
        pipelineVersion = 'default';
    end

    switch lower(pipelineVersion)
        case 'default'
            cfg.carrierFrame   = [];     % [] = mean peak over all bg frames
            cfg.tiltCorrection = true;   % re-centre FFT peak after bg division
            cfg.resOverride    = [];     % [] = use res stored in the .mat files
            cfg.allBackgrounds = true;   % pair every sample with every bg file
            cfg.madFactor      = 4;      % outlier thresholds: median + madFactor*MAD
            cfg.tiffFlipLR     = true;   % mirror TIFF slices left-right
        case 'alice'
            cfg.carrierFrame   = 49;     % peak of bg frame 49
            cfg.tiltCorrection = false;
            cfg.resOverride    = 0.133481157776655;  % hardcoded after loading each sample
            cfg.allBackgrounds = false;  % first bg*_Tomog.mat only
            cfg.madFactor      = 4;      % kept from default (her script used fixed 1.5 / 0.05)
            cfg.tiffFlipLR     = false;
        otherwise
            error('pipelineConfig:unknownVersion', ...
                'Unknown pipeline version "%s" (expected "default" or "alice").', pipelineVersion);
    end
    cfg.name = lower(pipelineVersion);
end

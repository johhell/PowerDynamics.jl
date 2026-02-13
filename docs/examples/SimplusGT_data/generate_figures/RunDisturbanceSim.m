%% RunDisturbanceSim.m
% Disturbance simulation script for Simplus Grid Tool.
% Run this script AFTER running UserMain.m (model must already be open).
%
% What it does:
%   1. Sets R of Branch1-1 to R_fault between t_fault_on and t_fault_off
%   2. Logs v_dq signals from DS-Bus1 through DS-Bus4
%   3. Runs the simulation for t_end seconds total (3 independent segments)
%   4. Saves results to a CSV file and plots voltage magnitude
%
% Each phase runs independently from t=0. The final Simulink state of each
% phase is passed as the initial state of the next. The time vectors are
% then stacked in post-processing to form a continuous timeline.
clc;

%% Configuration (edit these as needed)
modelName   = 'mymodel_v1';
branchPath  = [modelName '/Branch1-1'];
t_fault_on  = 0.1;     % Time to apply fault [s]
t_fault_off = 0.2;     % Time to clear fault [s]
t_end       = 30;      % Total simulation time [s]
R_fault     = '1e-4';  % Resistance during fault (string, same units as block)
busIndices  = 1:4;     % DS-Bus indices to log v_dq from

% Create figures directory with absolute path
scriptDir  = fileparts(mfilename('fullpath'));
figureDir  = fullfile(scriptDir, 'figures');
if ~exist(figureDir, 'dir')
    [status, msg] = mkdir(figureDir);
    if ~status
        error('Failed to create directory %s: %s', figureDir, msg);
    end
end
outputFile = fullfile(figureDir, 'disturbance_sim_results.csv');
fprintf('Output file will be: %s\n', outputFile);

%% Save original model settings
R_orig        = get_param(branchPath, 'Resistance');
origStopTime  = get_param(modelName, 'StopTime');
origStartTime = get_param(modelName, 'StartTime');

fprintf('Branch: %s\n', branchPath);
fprintf('Original Resistance: %s\n', R_orig);
fprintf('Disturbance: R = %s from %.3gs to %.3gs\n', R_fault, t_fault_on, t_fault_off);

%% Add To Workspace blocks for v_dq from each DS-Bus
set_param(modelName, 'ReturnWorkspaceOutputs', 'on');

activeBuses = [];
addedBlocks = {};

for k = busIndices
    dsName  = ['DS-Bus' num2str(k)];
    dsPath  = [modelName '/' dsName];
    twsName = ['ToWS_v_dq_Bus' num2str(k)];
    twsPath = [modelName '/' twsName];
    varName = ['v_dq_Bus' num2str(k)];

    if ~exist_block(modelName, dsName)
        warning('DS-Bus%d not found in model, skipping.', k);
        continue;
    end

    if ~exist_block(modelName, twsName)
        add_block('simulink/Sinks/To Workspace', twsPath);
        dsPos = get_param(dsPath, 'Position');
        twsPos = [dsPos(1), dsPos(4)+40, dsPos(3), dsPos(4)+60];
        set_param(twsPath, 'Position', twsPos);
    end
    set_param(twsPath, 'VariableName', varName);
    set_param(twsPath, 'SaveFormat', 'Timeseries');

    lh = get_param(twsPath, 'LineHandles');
    if lh.Inport(1) == -1
        add_line(modelName, [dsName '/1'], [twsName '/1']);
    end

    addedBlocks{end+1} = twsPath;
    activeBuses(end+1) = k;
end

if isempty(activeBuses)
    error('No DS-Bus blocks found. Check that the model is open.');
end
fprintf('Logging v_dq from DS-Bus: %s\n', num2str(activeBuses));

%% Configure state saving (no SaveCompleteFinalSimState to avoid timestamp issues)
set_param(modelName, 'SaveFinalState', 'on');
set_param(modelName, 'SaveCompleteFinalSimState', 'off');

%% ---- Phase 1: Normal operation (duration = t_fault_on) ----
fprintf('\n--- Phase 1: 0 to %gs (normal) ---\n', t_fault_on);
set_param(modelName, 'StartTime', '0');
set_param(modelName, 'StopTime', num2str(t_fault_on));
set_param(modelName, 'LoadInitialState', 'off');
set_param(modelName, 'FinalStateName', 'xFinal1');

simOut1 = sim(modelName);
xFinal1 = simOut1.get('xFinal1');
fprintf('Phase 1 complete.\n');

%% ---- Phase 2: Disturbance (duration = t_fault_off - t_fault_on) ----
dur2 = t_fault_off - t_fault_on;
fprintf('\n--- Phase 2: %gs to %gs (R = %s, sim duration %gs) ---\n', ...
    t_fault_on, t_fault_off, R_fault, dur2);
set_param(branchPath, 'Resistance', R_fault);
set_param(modelName, 'StartTime', '0');
set_param(modelName, 'StopTime', num2str(dur2));
set_param(modelName, 'LoadInitialState', 'on');
set_param(modelName, 'InitialState', 'xFinal1');
set_param(modelName, 'FinalStateName', 'xFinal2');

simOut2 = sim(modelName);
xFinal2 = simOut2.get('xFinal2');
fprintf('Phase 2 complete.\n');

%% ---- Phase 3: Recovery (duration = t_end - t_fault_off) ----
dur3 = t_end - t_fault_off;
fprintf('\n--- Phase 3: %gs to %gs (R restored, sim duration %gs) ---\n', ...
    t_fault_off, t_end, dur3);
set_param(branchPath, 'Resistance', R_orig);
set_param(modelName, 'StartTime', '0');
set_param(modelName, 'StopTime', num2str(dur3));
set_param(modelName, 'InitialState', 'xFinal2');
set_param(modelName, 'FinalStateName', 'xFinal3');

simOut3 = sim(modelName);
fprintf('Phase 3 complete.\n');

%% Remove temporary To Workspace blocks and restore model
for i = 1:length(addedBlocks)
    lh = get_param(addedBlocks{i}, 'LineHandles');
    if lh.Inport(1) ~= -1
        delete_line(lh.Inport(1));
    end
    delete_block(addedBlocks{i});
end

set_param(branchPath, 'Resistance', R_orig);
set_param(modelName, 'StartTime', origStartTime);
set_param(modelName, 'StopTime', origStopTime);
set_param(modelName, 'LoadInitialState', 'off');
set_param(modelName, 'SaveFinalState', 'off');
fprintf('\nModel settings restored, temporary blocks removed.\n');

%% Extract and stack signals
% Each phase ran from t=0. We offset Phase 2 by t_fault_on and Phase 3
% by t_fault_off to reconstruct the continuous timeline.
fprintf('Extracting and stacking logged signals...\n');

allTime = [];
allData = [];
colNames = {'Time'};

for k = activeBuses
    varName = ['v_dq_Bus' num2str(k)];

    ts1 = simOut1.get(varName);
    ts2 = simOut2.get(varName);
    ts3 = simOut3.get(varName);

    d1 = squeeze(ts1.Data);
    d2 = squeeze(ts2.Data);
    d3 = squeeze(ts3.Data);

    % Ensure column orientation
    if size(d1,1) == 1 && numel(d1) > 1; d1 = d1'; end
    if size(d2,1) == 1 && numel(d2) > 1; d2 = d2'; end
    if size(d3,1) == 1 && numel(d3) > 1; d3 = d3'; end

    % Stack time: offset each phase, skip first sample of Phase 2 & 3
    % to avoid duplicate points at the boundaries
    if isempty(allTime)
        t1 = ts1.Time;
        t2 = ts2.Time + t_fault_on;
        t3 = ts3.Time + t_fault_off;
        allTime = [t1; t2(2:end); t3(2:end)];
    end

    allData = [allData, [d1; d2(2:end,:); d3(2:end,:)]];

    nSig = size(d1, 2);
    if nSig == 2
        colNames = [colNames, {['v_d_Bus' num2str(k)], ['v_q_Bus' num2str(k)]}];
    else
        for j = 1:nSig
            colNames = [colNames, {sprintf('v_dq%d_Bus%d', j, k)}];
        end
    end
end

%% Write CSV
outMatrix = [allTime, allData];

% Write CSV file
fid = fopen(outputFile, 'w');
if fid == -1
    error('Cannot open file for writing: %s\nCurrent directory: %s', outputFile, pwd);
end
fprintf(fid, '%s\n', strjoin(colNames, ','));
fclose(fid);
dlmwrite(outputFile, outMatrix, '-append', 'precision', '%.10g');

fprintf('\nResults saved to: %s\n', fullfile(pwd, outputFile));
fprintf('  %d time steps, %d columns\n', size(outMatrix, 1), size(outMatrix, 2));
fprintf('  Columns: %s\n', strjoin(colNames, ', '));
fprintf('\nRun ExportFigures.m to generate and save all plots.\n');

%% Helper function
function tf = exist_block(sys, name)
    try
        get_param([sys '/' name], 'Handle');
        tf = true;
    catch
        tf = false;
    end
end

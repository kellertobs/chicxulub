clear; close all; clc;

par_Kardla;

% SET MODEL PARAMETERS
runID   = 'Kardla'; % run identifier tag
nout    = 50;       % print output every 'nop' steps
svout   = 1;        % save figures and data to file (1)

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

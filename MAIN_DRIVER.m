%% This is the main driver file for calculating differential TT and attenuation
% 
%  This file demonstrates the workflow for obtaining seismic data,
%  processing it, and computing differential travel time and differential
%  attenuation. This script is not necessarily intended to be executed all
%  at once, but could be in principle.
% 
%  Each of the constituent scripts in this workflow has various parameters
%  defined at the top which are important to check and customise according
%  to your needs.
% 
%  Z. Eilon 2016


%% ========================================================================
%% ================================  PREP  ================================  
%% ========================================================================
global wd datadir dbnam dbdir resdir
wd = 'path_to_top_directory'; % MASTER.m sits in this directory
datadir = 'path_to_data_directory'; % where the data will sit (e.g. on a separate hard drive)
resdir = 'path_to_results_directory'; % where results will be stored
dbnam = 'name_of_antelope_database';
dbdir = 'path_to_antelope_database';


%% add paths 
% path to all the necessary sub functions
addpath([wd,'/matguts'])
% not essential, but useful: path to colour maps from seizmo - GET THIS at http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/
addpath('~/Documents/MATLAB/seizmo-master/cmap/'); 

%% ========================================================================
%% =======================  GET + PROCESS DATA  ===========================
%% ========================================================================

% Download the data from IRIS
evdata_1_DOWNLOAD_IRIS;

% Rotate the components into NEZ format - only need if OBS
evdata_2_ROTATE;

% Remove instrument response using the SAC poles and zeros file that you 
evdata_3_RMRESP_SACPZ;

% Make arrival structures - parse data into handy event-based structures
evdata_4_ARRIVAL_STRUCT;

%% ========================================================================
%% ===============  CALCULATE TRAVEL TIME + ATTEN THINGS  =================  
%% ========================================================================

% Calculate differential travel times through cross correlation with manual
% windowing
calc_1_DIFFERENTIAL_TT;

% Compute power spectra and calc. spectral ratios
calc_2_SPECTRAL_RATIOS;

% Run through combs of filters and calculate attenuation from joint
% inversion of pairwise amplitude and phase spectra
calc_3_COMBSPECTRA;

% Extract the results into large [nevts*nstas] structures 
calc_4_RESULTS_EXTRACT

%Optional - instead of using the dtstar results calculated within
%calc_3_COMBSPECTRA which solves for each pair of stations and then for
%each event, we can solve for the station-wise dtstar values that fit all
%amplitude and phase measurements at once
calc_4b_ALL_IN_ONE_SOLVE_spectra

%% ========================================================================
%% ==========================  MAKE FIGURES  ==============================
%% ========================================================================

% Make figure of individual dtstar results 
mkfig_DTSTAR_BAZ_MAP

% Make figure of station-averaged dtstar results
mkfig_DTSTAR_STAV_MAP

% Make figure of attenuation across the array for a sinle event
mkfig_EVENT_ATTEN

% Plot the comb of filters used
mkfig_filter_comb






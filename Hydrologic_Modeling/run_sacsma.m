
% SACSMA Run for CALFEWS Watersheds

addpath('/home/fs02/pmr82_0001/rg727/WGEN/streamflow/matlab_module')

tic;

%create parallel pool for batch 
pc=parcluster('local')
pc.JobStorageLocation=strcat('/home/fs02/pmr82_0001/rg727/check_parallel_matlab/parallel_job_info')

parpool(pc,str2num(getenv('SLURM_CPUS_ON_NODE')))

%% HRU Information File for the whole CALFEWS area
hruinfo_calfews = load('/home/fs02/pmr82_0001/rg727/WGEN/streamflow/hruinfo/HRUinfo_ALL_CALFEWS.txt');
hru_lat_calfews = hruinfo_calfews(:,1); % HRU Lat
hru_lon_calfews = hruinfo_calfews(:,2); % HRU Long
hru_elev_calfews = hruinfo_calfews(:,4); % HRU elevation (m)


%% Calibration result file to load optimal parameters
fid = fopen('/home/fs02/pmr82_0001/rg727/WGEN/streamflow/sacramento_ga_cdec_pool_calfews_optpar.txt');
hru_par_calfews = cell(length(hru_lat_calfews),33);
while ~feof(fid)
    str = fgets(fid);
    if contains(str,'SACRAMENTO OPTIMAL PARAMETERS')
        fgets(fid);
        fgets(fid);
        n=1;
        for jj = 1:length(hru_lat_calfews)
            str = fgets(fid);
            hru_par_calfews(n,:) = textscan(str,['%s %s '...
                '%f %f %f %f %f %f %f %f %f %f '...
                '%f %f %f %f %f %f %f %f %f %f '...
                '%f %f %f %f %f %f %f %f %f %f %f ']);
            n=n+1;
        end
        break;
    end
end
fclose(fid);


%% HRU Information File for each watershed
%cdec_id = {'SHA','ORO','YRS','FOL','NML','TLG','MRC','MIL','PNF','TLG','SCC','ISB'};
cdec_id = {'TLG'};

for i=1:length(cdec_id)
    hruinfo = load(['/home/fs02/pmr82_0001/rg727/WGEN/streamflow/hruinfo/HRUinfo_',cdec_id{i},'.txt']);
    % HRU Lat
    eval(['hru_lat_',cdec_id{i},' = hruinfo(:,1);']) 
    % HRU Long
    eval(['hru_lon_',cdec_id{i},' = hruinfo(:,2);']) 
    % HRU Area fraction
    eval(['hru_area_',cdec_id{i},' = hruinfo(:,3);']) 
    % HRU flow length to the basin outlet (m)
    eval(['hru_flen_',cdec_id{i},' = hruinfo(:,5);']) 
end


%% Run SACSMA for the entire CALFEWS watersheds
% Simulation period
%sim_startdate = [1985 10 1];
%sim_enddate = [2013 9 30];
%sim_datemat = datevec(datenum(sim_startdate):datenum(sim_enddate));

sim_startdate = [1399 11 1];
sim_enddate = [2017 4 30];
sim_datemat = datevec(datenum(sim_startdate):datenum(sim_enddate));

% Module initial state
inistate = [0 0 5 5 5 0]; % for SACSMA [uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc]
snow_inistate = [0 0 0 0]; % for Snow17 [W_ice,W_liq,ATI,Deficit]


% Sacrameto HRU output initialization
hru_surf = zeros(size(sim_datemat,1),length(hru_lat_TLG));
hru_base = zeros(size(sim_datemat,1),length(hru_lat_TLG));
% tic;

%Adjust the parameter file
y=string(hru_par_calfews(:,1));
y_new=str2double(y);

[tf,idx] = ismember(hru_lat_TLG,y_new);
hru_par_calfews=hru_par_calfews(idx,:);
hru_elev_calfews=hru_elev_calfews(idx,:);

ensemble_TLG_sim=double.empty(225536,0);

for j=1:50
parfor i = 1:length(hru_lat_TLG)
    
    dates=parquetread(['/home/fs02/pmr82_0001/rg727/paleo_dates.parquetgzip']);
    year=double(table2array(dates(:,1)));
    month=double(table2array(dates(:,2)));
    day=double(table2array(dates(:,3)));
     
    precip=parquetread(['/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/perturbations/',num2str(T),'T_',num2str(C),'CC/precip/',num2str(j),'/meteo_',num2str(hru_lat_TLG(i),'%1.6f'),'_',num2str(hru_lon_TLG(i),'%1.6f'),'.parquetgzip']);
    precip_array=double(table2array(precip(:,1)));
    temp=parquetread(['/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/perturbations/',num2str(T),'T_',num2str(C),'CC/temp/',num2str(j),'/meteo_',num2str(hru_lat_TLG(i),'%1.6f'),'_',num2str(hru_lon_TLG(i),'%1.6f'),'.parquetgzip']); % yr mon day pr(mm) tas(C)
    temp_array=double(table2array(temp(:,1)));
    hru_meteo=zeros(225536,5);
    hru_meteo(:,1)=year;
    hru_meteo(:,2)=month;
    hru_meteo(:,3)=day;
    hru_meteo(:,4)=precip_array;
    hru_meteo(:,5)=temp_array;
    
    
    sind = find(hru_meteo(:,1)==sim_startdate(1)&hru_meteo(:,2)==sim_startdate(2)&hru_meteo(:,3)==sim_startdate(3)); %start simulation indices
    eind = find(hru_meteo(:,1)==sim_enddate(1)&hru_meteo(:,2)==sim_enddate(2)&hru_meteo(:,3)==sim_enddate(3)); % end simulation indices 
    pr = hru_meteo(sind:eind,4); % HRU pr for the simulation period
    tas = hru_meteo(sind:eind,5); % HRU tas for the simulation period
    
    selind = i;
    % Module parameters
    Coeff = hru_par_calfews{selind,3};
    SnowPar = [hru_par_calfews{selind,20},hru_par_calfews{selind,21},hru_par_calfews{selind,22},hru_par_calfews{selind,23},hru_par_calfews{selind,24},hru_par_calfews{selind,25},hru_par_calfews{selind,26},hru_par_calfews{selind,27},hru_par_calfews{selind,28},hru_par_calfews{selind,29}];
    SMA_Par = [hru_par_calfews{selind,4},hru_par_calfews{selind,5},hru_par_calfews{selind,6},hru_par_calfews{selind,7},hru_par_calfews{selind,8},hru_par_calfews{selind,9},hru_par_calfews{selind,10}...
        ,hru_par_calfews{selind,11},hru_par_calfews{selind,12},hru_par_calfews{selind,13},hru_par_calfews{selind,14},hru_par_calfews{selind,15},hru_par_calfews{selind,16},hru_par_calfews{selind,17},hru_par_calfews{selind,18},hru_par_calfews{selind,19}];
    
    % Run modules: Hamon -> Snow17 -> SACSMA
    pet = pet_hamon(sim_startdate, sim_enddate, tas, hru_lat_calfews(selind), Coeff);
    [snow_outflow, snow_melt, snow_swe, snow_inistate_new] = snow_snow17(sim_startdate, sim_enddate, pr, tas, hru_elev_calfews(selind), SnowPar, snow_inistate);
    [surf, base] = sma_sacramento(pet,snow_outflow,SMA_Par, inistate);
    
    hru_surf(:,i) = surf;
    hru_base(:,i) = base;
 
end

% toc % This step takes about 130 sec;

%% Run Lohmann Routing for each watershed
% tic;
for k = 1:length(cdec_id)

    eval(['hru_flen_sta = hru_flen_',cdec_id{k},';'])
    eval(['hru_lat_sta = hru_lat_',cdec_id{k},';'])
    eval(['hru_lon_sta = hru_lon_',cdec_id{k},';'])
    eval(['hru_area_sta = hru_area_',cdec_id{k},';'])
    
    hru_flow_sta = zeros(size(sim_datemat,1),length(hru_lat_sta));

    [~,minind] = min(hru_flen_sta);
    isOutlet = zeros(length(hru_flen_sta),1);
    isOutlet(minind) = 1;

     for i = 1:length(hru_lat_sta)
        selind = i;
        route_par = [hru_par_calfews{selind,30},hru_par_calfews{selind,31},hru_par_calfews{selind,32},hru_par_calfews{selind,33}];
        [totflow, baseflow] = rout_lohmann(hru_surf(:,selind), hru_base(:,selind), hru_flen_sta(i), route_par, isOutlet(i));
        hru_flow_sta(:,i) = totflow * hru_area_sta(i)/sum(hru_area_sta);
     end

    eval(['simflow_',cdec_id{k},' = sum(hru_flow_sta,2);'])
    
end
ensemble_TLG_sim(:,j)=simflow_TLG;
end
% toc; % This takes about 320 sec
%filename = sprintf('E:/sacsma_calfews_cornell/sacsma_calfews_cornell/Tuolumne_Ensemble_paleo.txt');
%writetable(table(sim_datemat,ensemble_TLG_sim),filename,'Delimiter','\t')

%filename = sprintf('E:/sacsma_calfews_cornell/sacsma_calfews_cornell/ANOVA/2T.txt');
%writetable(table(sim_datemat,ensemble_TLG_sim),filename,'Delimiter','\t')

%filename = sprintf('/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/streamflow/%dT_%.1fCC/streamflow_%dT_%.1fCC.txt',T,C,T,C);
%writetable(table(sim_datemat,ensemble_TLG_sim),filename,'Delimiter','\t')

filename = sprintf('/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/streamflow/%dT_%dCC/streamflow_%dT_%dCC_full.txt',T,C,T,C);
writetable(table(sim_datemat,ensemble_TLG_sim),filename,'Delimiter','\t')


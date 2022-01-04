%% Filter ONC Endeavour catalog for events to use as training data
% 10 minute windows

clear all;


%% Get list of all earthquakes good enough for usability

dirOutput = '/Users/zkrauss/Documents/GitHub/Axial-Endeavour-AutoLocate/EndeavourOutput/multi_it5';

startDate = datenum(2016,08,01);
endDate = datenum(2021,03,01);

locOut = load_location(dirOutput,startDate,endDate);
locOut = locOut([locOut.status]==0);

% List of P-wave SNR
p_snr = []; s_snr = [];
for i=1:length(locOut)
    p_sub = [];
    s_sub = [];
    for j = 1:length(locOut(i).pick.Ppeak)
        p_sub = [p_sub locOut(i).pick.Ppeak(j).dB];
    end
    p_snr = [p_snr p_sub];
    for j = 1:length(locOut(i).pick.Speak)
        s_sub = [s_sub locOut(i).pick.Speak(j).dB];
    end
    s_snr = [s_snr s_sub];
end

min_psnr = mean(p_snr);
min_ssnr = mean(s_snr);

figure(2);clf;
subplot(2,1,1);
histogram(p_snr);
subplot(2,1,2);
histogram(s_snr);

% Clear up some space
clear p_snr; clear s_snr;

% We want earthquakes that have both a P- and S-pick on the same station
% with both picks larger than the mean SNR (in dB) with weights >= 0.5
count = 0; template_data = [];
for i=1:length(locOut)
    % Get list of stations that have both P and S wave for this earthquake
    ppick = locOut(i).pick.Ppeak;
    spick = locOut(i).pick.Speak;
    if isempty(ppick) | isempty(spick)
        continue
    end
    [sta,iP,iS] = intersect([locOut(i).pick.Ppeak.station],[locOut(i).pick.Speak.station]);
    for j = 1:length(sta)
        % Make sure both picks have weight >= 0.25
        if isnan(ppick(iP(j)).weight) | isnan(spick(iS(j)).weight)
            continue
        elseif ppick(iP(j)).weight < 0.5 | spick(iS(j)).weight<0.5
            continue
        % Make sure both picks have SNR > mean
        elseif ppick(iP(j)).dB < min_psnr | spick(iS(j)).dB <min_ssnr
            continue
        end
        % If the station/earthquake passes the test, save its pick info
        % Station, earthquake ot, earthquake lat, earthquake lon, earthquake mag, nwr, P and S pick time, P and S pick weight, P
        % and S SNR in dB
        count = count + 1;
        
        template_data(count).station = sta(j);
        template_data(count).label = 'earthquake';
        template_data(count).eq_ot = locOut(i).nlloc_hypo.ot;
        template_data(count).eq_lat = locOut(i).nlloc_hypo.lat;
        template_data(count).eq_lon = locOut(i).nlloc_hypo.lon;
        template_data(count).ID = locOut(i).ID;
        template_data(count).nwr = length(locOut(i).nlloc_hypo.picks.stations);
        template_data(count).ptime = locOut(i).on+ppick(iP(j)).time1Kurt/86400;
        template_data(count).stime = locOut(i).on+spick(iS(j)).time1Kurt/86400;
        template_data(count).psnr_db = ppick(iP(j)).dB;
        template_data(count).ssnr_db = spick(iS(j)).dB;
        template_data(count).pweight = ppick(iP(j)).weight;
        template_data(count).sweight = spick(iS(j)).weight;
    end
end
            
            
%% Whittle down list to only earthquakes that are 2 minutes from another
% After all, PhaseNet only takes in 30 s


% Notes- 20 min on either side gives 5887 keepers, 10 min gives 7557, 5 min
% gives 8939, 2 min gives 10586

keep=[];
for i=1:length(template_data)
    tm = template_data(i).eq_ot;
    others = nonzeros([locOut.on] > datenum(tm-minutes(2)) & [locOut.on] < datenum(tm+minutes(2)));
    if length(others)<2
        keep = [keep i];
    end
end

template_data = template_data(keep);
%% Priority 2 - whale times

% Start by finding timestamps and stations for when whales are active at increments of
% 20 minutes ( > 10 calls per 20 minutes)

% Loop over days so as to not pull in too much data at once
day_bands = (datetime(2016,08,01):caldays(1):datetime(2021,03,01));
whale_counts = [];whale_times = [];count = 0;
for i=1:length(day_bands)
    % Load one day of events
    event = load_event(dirOutput,datenum(day_bands(i)),datenum(day_bands(i)+caldays(1)));
    if isempty(event)
        continue
    end
    % Consider only those events which are classified as whales
    event = event(strcmp('whale',{event.type}));
    time_bands = day_bands(i):minutes(20):(day_bands(i)+caldays(1));
    for j=1:length(time_bands)
        whales = event(nonzeros([event.on] > datenum(time_bands(j)) & [event.on] < datenum(time_bands(j)+minutes(20))));
        % If there were more than 10 whale events in that 20-minute band,
        % save the code of the stations that recorded them
        if length(whales)>10
            sta = unique([whales.station]);
            for k=1:length(sta)
                count = count + 1;
                whale_times(count).station = sta(k);
                whale_times(count).timeband = time_bands(j);

            end
        end
    end
    clear event; clear time_bands;
end

%% Loop through those same bands and make sure there are no earthquakes
% within the whale time band - if so, discard whale band

wkeep = [];
for i=1:length(whale_times)
    within = nonzeros([locOut.on] > datenum(whale_times(i).timeband) & [locOut.on] < datenum(whale_times(i).timeband + minutes(20)));
    if length(within)<1
        wkeep = [wkeep i];
    end
end
        
whale_times = whale_times(wkeep);
%% Save whale bands

% timestarts = rand([1 length(whale_times)]);
count = length(template_data);
for i=1:length(whale_times)
    count = count + 1;
    template_data(count).noise_start = whale_times(i).timeband;
    template_data(count).label = 'whale';
end

    


%% Priority 3 - noise at other times

% Find timestamps for periods that do not have whale calls and do not have
% earthquakes

% Random number seed
rng(1);
% Define how many timebands you want:
num_noise = 10;


search_bands = (datetime(2016,08,01):minutes(20):datetime(2021,03,01));
% Shuffle the search bands 
search_bands = search_bands(randperm(length(search_bands)));

% Now check each timeband to see if it doesn't contain an earthquake or
% whale band until you have the specified number of noise bands saved:
count = length(template_data);noise_count = 0;
for i=1:length(search_bands)
    band = search_bands(i);
    band_day = dateshift(band,'start','day');
    if any([template_data.noise_start]==band)
        continue
    elseif any([locOut.on] > datenum(band) & [locOut.on] < datenum(band + minutes(20)))
        continue
    % If time band passes, get a list of stations operating at that time
    else
        pickOut = load_pick(dirOutput,datenum(band_day),datenum(band_day+caldays(1)));
        if ~isempty(pickOut)
            sta = unique([pickOut(1).Speak.station]);
            for j=1:length(sta)
                count = count + 1;
                template_data(count).noise_start = band;
                template_data(count).label = 'noise';
                template_data(count).station = sta(j);
                noise_count = noise_count + 1;
              
            end
        end
    end
    if noise_count==num_noise
        break
    end
end
    
%% Next step- vary where the P-pick starts in the window (consult the PhaseNet training data here)
% This should probably be done in python. All we want from this is a list
% of times and stations.
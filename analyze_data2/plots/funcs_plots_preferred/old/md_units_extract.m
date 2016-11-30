
function myspikesmm = md_units_extract(md)
%     This function extracts the unit time series from the data in md. This is probably
%     not useful, since it doesn't take into account the different ctg schemes. I might
%     should instead just use the run_analysis script which outputs spike data to mode5.1.
%     This does include organization by category type.

    units_per_file = cellfun(@length, {md(:).unit_names});
    unit_nums = arrayfunu(@(N) 1:N,units_per_file);
    myspikes = cellfunu(@(md,unit_nums) arrayfunu(@(curr_unit) get_spike_mean_ts(curr_unit,5,md), unit_nums )  ,num2cell(md),unit_nums);
    myspikesm = cellfunu(@(x) cell2mat(x),myspikes);
    myspikesm = cellfunu(@(x) x / get_dt, myspikesm);
    myspikesmm = cell2mat(myspikesm);
    
end



function [currspike] = get_spike_mean_ts(varargin)

%     [currspike] = get_spike_mean_ts(curr_unit,curr_stage,md)
%     INPUT
%     unit_ind - Units indices
%     curr_unit - Current unit to look at. The code will automagically pair this up with
%         the correct electrode based on metadata
%     curr_stage - Current stage of recording to look at (sample, delay, both, etc)
%     md - Metadata structure
% 
%     OUTPUT
%     currspike - Time series for spikes

    apply_filter=1;

    [currspike] = get_spike_ts(varargin{:});
    
    currspike = mean(currspike,2);
    
    if apply_filter
        dt = get_dt;
        currspike = sgolayfilt(currspike,3,round(0.151/dt));   % Apply 3rd-order filter, 151 ms 
    end
    
end




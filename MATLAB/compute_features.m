clear all
close all
clc

path = 'C:\git\acosize_ts';

freq    = '70khz';
ramping = {'slow','fast'};
orientation = {'dorsal'};

data_path   = fullfile(path,'data');
results_path    = fullfile(path,'results','measurements');

windowing = 'variableWindow'; % fixedWindow variableWindow

load(fullfile(data_path,'mat',strcat('acoSize_TS_',freq,'_',windowing,'.mat')))

for idx = 1:length(outStruct)
    char(outStruct(idx).fish_id)
    
    for idxTS = 1:size(outStruct(idx).TS,2)
        if ~isempty(outStruct(idx).TS(idxTS).TSMat)
            outStruct(idx).TS(idxTS).f_vec_inter = ceil(min(outStruct(idx).TS(idxTS).f_vec)*1e-3)*1e3:1e3:floor(max(outStruct(idx).TS(idxTS).f_vec)*1e-3)*1e3;
            outStruct(idx).TS(idxTS).TSMat_inter = interp1(outStruct(idx).TS(idxTS).f_vec,outStruct(idx).TS(idxTS).TSMat',outStruct(idx).TS(idxTS).f_vec_inter)';
            outStruct(idx).TS(idxTS).TSkde_inter = interp1(outStruct(idx).TS(idxTS).f_vec,outStruct(idx).TS(idxTS).TSkde',outStruct(idx).TS(idxTS).f_vec_inter)';
            outStruct(idx).TS(idxTS).TSStats_inter = interp1(outStruct(idx).TS(idxTS).f_vec,outStruct(idx).TS(idxTS).TSStats',outStruct(idx).TS(idxTS).f_vec_inter)';
        end
    end
end

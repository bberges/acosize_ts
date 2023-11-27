clear all
close all
clc

path = 'G:\git\acosize_ts';

freq    = '200khz';
ramping = {'slow','fast'};
orientation = {'ventral'};

data_path   = fullfile(path,'data');
results_path    = fullfile(path,'results','measurements');

load(fullfile(data_path,'mat',strcat('acoSize_TS_',freq,'_cut.mat')))

measurement_summary = struct(   'fish_id',{outStruct.fish_id}, ...
                                'type',{outStruct.type}, ...
                                'mode',{outStruct.mode}, ...
                                'ramping',{outStruct.ramping}, ...
                                'resolution',{outStruct.resolution}, ...
                                'orientation',{outStruct.orientation}, ...
                                'fish_length',{outStruct.fish_length}, ...
                                'length_4',{outStruct.length_4}, ...
                                'width',{outStruct.width}, ...
                                'length_1',{outStruct.length_1}, ...
                                'length_2',{outStruct.length_2}, ...
                                'length_3',{outStruct.length_3}, ...
                                'length_PVA',{outStruct.length_PVA}, ...
                                'SNR',{outStruct.SNR});
summaryTab = struct2table(measurement_summary);
writetable(summaryTab,fullfile(results_path,'summaryTab.csv'),'Delimiter',',')

write_path    = fullfile(path,'results','measurements',freq);

% export TS into CSV files
for idx = 1:length(outStruct)
    idx/length(outStruct)
    currentName = char(strcat( outStruct(idx).fish_id,'_', ...
                            outStruct(idx).mode,'_', ...
                            outStruct(idx).ramping,'_', ...
                            outStruct(idx).orientation,'_', ...
                            outStruct(idx).resolution,'_', ...
                            freq));
                        
	for idxTS = 1:size(outStruct(idx).TS,2)
        if ~isempty(outStruct(idx).TS(idxTS).TSMat)
            outStruct(idx).TS(idxTS).f_vec_inter = ceil(min(outStruct(idx).TS(idxTS).f_vec)*1e-3)*1e3:1e3:floor(max(outStruct(idx).TS(idxTS).f_vec)*1e-3)*1e3;
            outStruct(idx).TS(idxTS).TSMat_inter = interp1(outStruct(idx).TS(idxTS).f_vec,outStruct(idx).TS(idxTS).TSMat',outStruct(idx).TS(idxTS).f_vec_inter)';
            outStruct(idx).TS(idxTS).TSkde_inter = interp1(outStruct(idx).TS(idxTS).f_vec,outStruct(idx).TS(idxTS).TSkde',outStruct(idx).TS(idxTS).f_vec_inter)';
            outStruct(idx).TS(idxTS).TSStats_inter = interp1(outStruct(idx).TS(idxTS).f_vec,outStruct(idx).TS(idxTS).TSStats',outStruct(idx).TS(idxTS).f_vec_inter)';
        end
    end

%             plot(outStruct(idx).TS(idxTS).f_vec,outStruct(idx).TS(idxTS).TSStats(2,:), ...
%                 outStruct(idx).TS(idxTS).f_vec_inter,outStruct(idx).TS(idxTS).TSStats_inter(2,:))
    
    % write TS correlation
    TSCorr = array2table([outStruct(idx).headingPercentilesVec', outStruct(idx).TSCorr'], 'VariableNames', {'angle', 'TSCorr'});
    writetable(TSCorr,fullfile(write_path,strcat(currentName,'_TSCorr.csv')),'Delimiter',',')
    
    idxTemp=1;
    while isempty(outStruct(idx).TS(idxTemp).f_vec) && idxTemp < length(outStruct(idx).TS)
        idxTemp = idxTemp+1;
    end
    
    colNamesTSStats = [cellstr('angle'); cellstr('N_spectra');cellstr('percentiles');cellstr(num2str(outStruct(idx).TS(idxTemp).f_vec_inter'))];
    TSStatsTab = array2table(   zeros(0,size(outStruct(idx).TS(idxTemp).TSStats_inter,2)+3), ...
                                'VariableNames', colNamesTSStats);
                            
	colNamesTSkde = [cellstr('angle'); cellstr('N_spectra');cellstr('TSmesh');cellstr(num2str(outStruct(idx).TS(idxTemp).f_vec_inter'))];
    TSkdeTab = array2table(   zeros(0,size(outStruct(idx).TS(idxTemp).TSkde_inter,2)+3), ...
                                'VariableNames', colNamesTSkde);

    % make TS densities
    flagFirst = true;
    for idxAngle = 1:length(outStruct(idx).headingPercentilesVec)
        if ~isempty(outStruct(idx).TS(idxAngle).TSStats)
            TSStats = [repmat(outStruct(idx).headingPercentilesVec(idxAngle),3,1), ...
                        repmat(size(outStruct(idx).TS(idxAngle).TSMat_inter,1),3,1), ...
                        [25 50 75]', ...
                        outStruct(idx).TS(idxAngle).TSStats_inter];
            TSStatsCurrent = array2table(TSStats, 'VariableNames', colNamesTSStats);
            TSStatsTab = [TSStatsTab;TSStatsCurrent];
            
            TSkde = [repmat(outStruct(idx).headingPercentilesVec(idxAngle),length(outStruct(idx).TS(idxAngle).xmesh),1), ...
                        repmat(size(outStruct(idx).TS(idxAngle).TSMat_inter,1),length(outStruct(idx).TS(idxAngle).xmesh),1), ...
                        outStruct(idx).TS(idxAngle).xmesh', ...
                        outStruct(idx).TS(idxAngle).TSkde_inter];
                    
            TSkdeCurrent = array2table(TSkde,'VariableNames', colNamesTSkde);
            TSkdeTab = [TSkdeTab;TSkdeCurrent];
        end
    end
    
    writetable(TSStatsTab,fullfile(write_path,strcat(currentName,'_TSStats.csv')),'Delimiter',',')
    writetable(TSkdeTab,fullfile(write_path,strcat(currentName,'_TSkde.csv')),'Delimiter',',')
end

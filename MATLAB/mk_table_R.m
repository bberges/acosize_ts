clear all
close all
clc

path = 'V:\acosize_ts';

freq    = '200kHz';
ramping = {'slow','fast'};
orientation = {'dorsal'};

data_path   = fullfile(path,'data');
results_path    = fullfile(path,'results','measurements');

windowing = 'variableWindow'; % fixedWindow variableWindow

load(fullfile(data_path,'mat',strcat('acoSize_TS_',freq,'_',windowing,'.mat')))

 

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
%     idx/length(outStruct)
    char(outStruct(idx).fish_id)
    currentName = char(strcat( outStruct(idx).fish_id,'_', ...
                            outStruct(idx).mode,'_', ...
                            outStruct(idx).ramping,'_', ...
                            outStruct(idx).orientation,'_', ...
                            outStruct(idx).resolution,'_', ...
                            windowing, '_', ...
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

    colNamesTSMat = [cellstr('fish_id'); ...
                     cellstr('SNR'); ...
                     cellstr('type'); ...
                     cellstr('mode'); ...
                     cellstr('ramping'); ...
                     cellstr('resolution'); ...
                     cellstr('orientation'); ...
                     cellstr('fish_length'); ...
                     cellstr('width'); ...
                     cellstr('angle'); ...
                     cellstr('spectra_id'); ...
                     cellstr('idxPing'); ...
                     cellstr('nullsN'); ...
                     cellstr('nullsProm'); ...
                     cellstr('nullsWidthMean'); ...
                     cellstr('nullsWidthStd'); ...
                     cellstr('nullsMean'); ...
                     cellstr(num2str(outStruct(idx).TS(idxTemp).f_vec_inter'))];
    colNamesTSMatTemp = [cellstr('fish_id'); ...
                         cellstr('SNR'); ...
                         cellstr('type'); ...
                         cellstr('mode'); ...
                         cellstr('ramping'); ...
                         cellstr('resolution'); ...
                         cellstr('orientation'); ...
                         cellstr('fish_length'); ...
                         cellstr('width'); ...
                         cellstr('angle'); ...
                         cellstr('spectra_id'); ...
                         cellstr('nullsN'); ...
                         cellstr('nullsProm'); ...
                         cellstr('nullsWidthMean'); ...
                         cellstr('nullsWidthStd'); ...
                         cellstr('nullsMean')];
    TSMatTab = array2table(   zeros(0,size(outStruct(idx).TS(idxTemp).TSStats_inter,2)+17), ...
                                'VariableNames', colNamesTSMat);
    
    colNamesTSStats = [cellstr('angle'); cellstr('N_spectra');cellstr('percentiles');cellstr(num2str(outStruct(idx).TS(idxTemp).f_vec_inter'))];
    TSStatsTab = array2table(   zeros(0,size(outStruct(idx).TS(idxTemp).TSStats_inter,2)+3), ...
                                'VariableNames', colNamesTSStats);
                            
	colNamesTSkde = [cellstr('angle'); cellstr('N_spectra');cellstr('TSmesh');cellstr(num2str(outStruct(idx).TS(idxTemp).f_vec_inter'))];
    TSkdeTab = array2table(   zeros(0,size(outStruct(idx).TS(idxTemp).TSkde_inter,2)+3), ...
                                'VariableNames', colNamesTSkde);
                            
    colNamesTimeMatTemp = [cellstr('fish_id'); ...
                         cellstr('SNR'); ...
                         cellstr('type'); ...
                         cellstr('mode'); ...
                         cellstr('ramping'); ...
                         cellstr('resolution'); ...
                         cellstr('orientation'); ...
                         cellstr('fish_length'); ...
                         cellstr('width'); ...
                         cellstr('angle'); ...
                         cellstr('idxPing'); ...
                         cellstr('meanVal'); ...
                         cellstr('SNR_sig'); ...
                         cellstr('mu'); ...
                         cellstr('sigma'); ...
                         cellstr('kurt'); ...
                         cellstr('skew'); ...
                         cellstr('Npeaks'); ...
                         cellstr('widthPeaks'); ...
                         cellstr('promPeaks'); ...
                         cellstr('stdPromPeaks'); ...
                         cellstr('stdWidthPeaks'); ...
                         cellstr('meanPeaks'); ...
                         cellstr('stdPeaks')];
                     
    TimeMatTab = array2table(   zeros(0,length(colNamesTimeMatTemp)), ...
                                'VariableNames', colNamesTimeMatTemp);

    % make TS densities
    for idxAngle = 1:length(outStruct(idx).headingPercentilesVec)
        if ~isempty(outStruct(idx).TS(idxAngle).TSStats)
            TSMatCurrent  = cell2table(cell(size(outStruct(idx).TS(idxAngle).TSMat_inter,1),length(colNamesTSMatTemp)),...
                                        'VariableNames', colNamesTSMatTemp);
            TSMatCurrent.fish_id(:) = outStruct(idx).fish_id;
            TSMatCurrent.SNR(:) = {outStruct(idx).SNR};
            TSMatCurrent.type(:) = outStruct(idx).type;
            TSMatCurrent.mode(:) = outStruct(idx).mode;
            TSMatCurrent.ramping(:) = outStruct(idx).ramping;
            TSMatCurrent.resolution(:) = outStruct(idx).resolution;
            TSMatCurrent.orientation(:) = outStruct(idx).orientation;
            TSMatCurrent.fish_length(:) = {outStruct(idx).fish_length};
            TSMatCurrent.width(:) = {outStruct(idx).width};
            TSMatCurrent.angle(:) = {outStruct(idx).headingPercentilesVec(idxAngle)};
            TSMatCurrent.spectra_id(:) = num2cell(1:size(outStruct(idx).TS(idxAngle).TSMat_inter,1));
            TSMatCurrent.idxPing(:) = num2cell(outStruct(idx).TS(idxAngle).idxTS);
            
            for idxSpectra = 1:size(outStruct(idx).TS(idxAngle).TSMat,1)
                sig = (outStruct(idx).TS(idxAngle).TSMat(idxSpectra,:)-mean(outStruct(idx).TS(idxAngle).TSMat(idxSpectra,:)))/ ...
                        std(outStruct(idx).TS(idxAngle).TSMat(idxSpectra,:));
                [~,locs,w,p]  = findpeaks( 	-sig, ...
                                            outStruct(idx).TS(idxAngle).f_vec, ...
                                            'WidthReference','halfheight');
                                        
                idxFilt = p>4 & outStruct(idx).TS(idxAngle).TSMat(idxSpectra,ismember(outStruct(idx).TS(idxAngle).f_vec,locs))<-70;
                w = w(idxFilt);
                locs = locs(idxFilt);
                p = p(idxFilt);
                valNull = outStruct(idx).TS(idxAngle).TSMat(idxSpectra,ismember(outStruct(idx).TS(idxAngle).f_vec,locs));
                Npeaks = length(locs);
                
                if Npeaks >= 1
                    TSMatCurrent.nullsN(idxSpectra) = {Npeaks};
                    TSMatCurrent.nullsProm(idxSpectra) = {mean(p)};
                    TSMatCurrent.nullsWidthMean(idxSpectra) = {mean(w)};
                    TSMatCurrent.nullsWidthStd(idxSpectra) = {std(w)};
                    TSMatCurrent.nullsMean(idxSpectra) = {mean(valNull)};
                else
                    TSMatCurrent.nullsN(idxSpectra) = {0};
                    TSMatCurrent.nullsProm(idxSpectra) = {0};
                    TSMatCurrent.nullsWidthMean(idxSpectra) = {0};
                    TSMatCurrent.nullsWidthStd(idxSpectra) = {0};
                    TSMatCurrent.nullsMean(idxSpectra) = {0};
                end
            end
            
            temp = array2table(outStruct(idx).TS(idxAngle).TSMat_inter, ...
                                'VariableNames', cellstr(num2str(outStruct(idx).TS(idxTemp).f_vec_inter')));
            TSMatCurrent = [TSMatCurrent temp];
            TSMatTab = [TSMatTab;TSMatCurrent];
            
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
            
            % time features
            TimeMatCurrent  = cell2table(cell(length(outStruct(idx).TS(idxAngle).idxTS),length(colNamesTimeMatTemp)),...
                                        'VariableNames', colNamesTimeMatTemp);
                                    
            TimeMatCurrent.fish_id(:) = outStruct(idx).fish_id;
            TimeMatCurrent.SNR(:) = {outStruct(idx).SNR};
            TimeMatCurrent.type(:) = outStruct(idx).type;
            TimeMatCurrent.mode(:) = outStruct(idx).mode;
            TimeMatCurrent.ramping(:) = outStruct(idx).ramping;
            TimeMatCurrent.resolution(:) = outStruct(idx).resolution;
            TimeMatCurrent.orientation(:) = outStruct(idx).orientation;
            TimeMatCurrent.fish_length(:) = {outStruct(idx).fish_length};
            TimeMatCurrent.width(:) = {outStruct(idx).width};
            TimeMatCurrent.angle(:) = {outStruct(idx).headingPercentilesVec(idxAngle)};
            TimeMatCurrent.idxPing(:) = num2cell(outStruct(idx).TS(idxAngle).idxTS);
            
            for idxPingTS = 1:length(outStruct(idx).TS(idxAngle).idxTS)
                idxCurrent = outStruct(idx).TS(idxAngle).idxTS(idxPingTS);
                minRange = min(outStruct(idx).TS(idxAngle).rangeTS(idxPingTS,1));
                maxRange   = max(outStruct(idx).TS(idxAngle).rangeTS(idxPingTS,2));
                [~,idxStartRange] = min(abs(outStruct(idx).range - minRange));
                [~,idxEndRange] = min(abs(outStruct(idx).range - maxRange));
                
                rangeSel = outStruct(idx).range(idxStartRange:idxEndRange);
                sig = outStruct(idx).SvMat(idxCurrent,idxStartRange:idxEndRange);
                
                sig(isinf(sig)) = NaN;
                %pd = fitdist(10*log10(abs(timeSig(idxTrans).data(:,idxCell))),'Normal');
                pd = fitdist(sig','Normal');
                TimeMatCurrent.mu(idxPingTS) = {pd.mu};
                TimeMatCurrent.sigma(idxPingTS) = {pd.sigma};
                TimeMatCurrent.kurt(idxPingTS) = {kurtosis(sig)};
                TimeMatCurrent.skew(idxPingTS) = {skewness(sig)};
                
                [~,locs,w,p]  = findpeaks( 	(sig-mean(sig))/std(sig), ...
                                            rangeSel, ...
                                            'WidthReference','halfheight');
                                        
                idxFilt = p>1 & sig(ismember(rangeSel,locs))>-50;
                w = w(idxFilt);
                locs = locs(idxFilt);
                pnorm = p(idxFilt);
                p = sig(ismember(rangeSel,locs));
                Npeaks = length(locs);
                
                idxSNRsig = [];
                for idxPeakSel = 1:Npeaks
                    [~,idxPeakStart] = min(abs(rangeSel-(locs(idxPeakSel)-w(idxPeakSel)/2)));
                    [~,idxPeakend] = min(abs(rangeSel-(locs(idxPeakSel)+w(idxPeakSel)/2)));
                    
                    idxSNRsig = [idxSNRsig,idxPeakStart:idxPeakend];
                end
                
                idxSNRnoise = find(~ismember(1:length(rangeSel),idxSNRsig));
                
                sigSNR = 10.^(sig(idxSNRsig)/20);
                Ps = sum(abs(sigSNR).^2)/length(sigSNR);
                noise = 10.^(sig(idxSNRnoise)/20);
                Pn = sum(abs(noise).^2)/length(noise);
                
                SNR_sig = 10*log10(Ps/Pn);
                
                if Npeaks > 1
                    widthPeaks     = mean(w);
                    promPeaks      = mean(pnorm);
                    meanPeaks      = mean(p);
                    stdWidthPeaks  = std(w);
                    stdPromPeaks   = std(pnorm);
                    stdPeaks       = std(p);
                elseif Npeaks == 1
                    widthPeaks     = mean(w);
                    promPeaks      = mean(pnorm);
                    meanPeaks      = mean(p);
                    stdWidthPeaks  = std(w);
                    stdPromPeaks   = std(pnorm);
                    stdPeaks       = std(p);
                else
                    widthPeaks     = 0;
                    promPeaks      = 0;
                    meanPeaks      = mean(p);
                    stdWidthPeaks  = 0;
                    stdPromPeaks   = 0;
                    stdPeaks       = std(p);
                end
                
                TimeMatCurrent.meanVal(idxPingTS) = {mean(sig)};
                TimeMatCurrent.Npeaks(idxPingTS) = {Npeaks}; 
                TimeMatCurrent.widthPeaks(idxPingTS) = {widthPeaks};
                TimeMatCurrent.promPeaks(idxPingTS) = {promPeaks};
                TimeMatCurrent.stdPromPeaks(idxPingTS) = {stdPromPeaks};
                TimeMatCurrent.stdWidthPeaks(idxPingTS) = {stdWidthPeaks};
                TimeMatCurrent.meanPeaks(idxPingTS) = {meanPeaks};
                TimeMatCurrent.stdPeaks(idxPingTS) = {stdPeaks};
                TimeMatCurrent.SNR_sig(idxPingTS) = {SNR_sig};
            end
            TimeMatTab = [TimeMatTab;TimeMatCurrent];
        end
    end
    
    writetable(TimeMatTab,fullfile(write_path,strcat(currentName,'_TimeMat.csv')),'Delimiter',',')
    writetable(TSMatTab,fullfile(write_path,strcat(currentName,'_TSMat.csv')),'Delimiter',',')
    writetable(TSStatsTab,fullfile(write_path,strcat(currentName,'_TSStats.csv')),'Delimiter',',')
    writetable(TSkdeTab,fullfile(write_path,strcat(currentName,'_TSkde.csv')),'Delimiter',',')
end

function Features = computeTimeFeaturesWMR_v4(timeSig, range, mask,maxCellDist,nAveragingCells)

nWindows = length(range)/size(timeSig(1).data,1);
nPings = size(timeSig(1).data,2)/nWindows;

maskMat = reshape(mask, nWindows,nPings);

R = reshape(range, size(timeSig(1).data,1),nWindows);
R = repmat(R,1,nPings);

%% compute features
SvMat       = zeros(length(timeSig),size(timeSig(1).data,2));
SNRMat      = zeros(length(timeSig),size(timeSig(1).data,2));
muMat       = zeros(length(timeSig),size(timeSig(1).data,2));
sigmaMat    = zeros(length(timeSig),size(timeSig(1).data,2));
kurtMat     = zeros(length(timeSig),size(timeSig(1).data,2));
skewMat     = zeros(length(timeSig),size(timeSig(1).data,2));
NpeaksMat         = zeros(length(timeSig),size(timeSig(1).data,2));
widthPeaksMat     = zeros(length(timeSig),size(timeSig(1).data,2));
promPeaksMat      = zeros(length(timeSig),size(timeSig(1).data,2));
stdWidthPeaksMat  = zeros(length(timeSig),size(timeSig(1).data,2));
stdPromPeaksMat   = zeros(length(timeSig),size(timeSig(1).data,2));

for idxTrans = 1:length(timeSig)
    idxTransOther = find([1:length(timeSig)]~=idxTrans);
    % % % % % % % computing features
    % combining pixels in windows
    % SvVal = 10*log10(sum(abs(timeSig(idxTrans).data).^2.*R.^2./sum(R.^2,1)));
    SvVal = 10*log10(sum(abs(timeSig(idxTrans).data).^2.*R.^2./repmat(sum(R.^2,1), size(timeSig(idxTrans).data, 1), 1))); % R2015b compatible	
    SvValMat = reshape(SvVal, nWindows,nPings);
    
    [ImatIdxNoise,JmatIdxNoise] = ind2sub([nWindows nPings],find(isnan(mask)));

    
    % SNR
    % What if < 3 windows, need to keep fail-safe in case of temporary 
    % weird (e.g. empty sonar data erroneous bottom detection etc.    
    % SNR = SvValMat - mean(SvValMat(3,:)); 
%     SNR = SvValMat - mean(SvValMat(min(3,end),:)); 
%     SNR = SNR(:)';
    
    mu = zeros(size(timeSig(idxTrans).data,2),1);
    sig = zeros(size(timeSig(idxTrans).data,2),1);
    kurt = zeros(size(timeSig(idxTrans).data,2),1);
    skew = zeros(size(timeSig(idxTrans).data,2),1);
    Npeaks          = zeros(size(timeSig(idxTrans).data,2),1);
	widthPeaks      = zeros(size(timeSig(idxTrans).data,2),1);
	promPeaks       = zeros(size(timeSig(idxTrans).data,2),1);
	stdWidthPeaks   = zeros(size(timeSig(idxTrans).data,2),1);
	stdPromPeaks    = zeros(size(timeSig(idxTrans).data,2),1);
    SNR             = zeros(size(timeSig(idxTrans).data,2),1);
    for idxCell = 1:size(timeSig(idxTrans).data,2)
        if mask(idxCell) == 1
            
            [I,J] = ind2sub([nWindows nPings],idxCell);
            
            SNR(idxCell) = SvVal(idxCell) - mean(SvVal((I > ImatIdxNoise) & (J == JmatIdxNoise)));
            % normal distribution fit
            % JS 2018/1018: bugfix on inf values
            values = 10*log10(abs(timeSig(idxTrans).data(:,idxCell)));
            values(isinf(values)) = NaN;
            %pd = fitdist(10*log10(abs(timeSig(idxTrans).data(:,idxCell))),'Normal');
            pd = fitdist(values,'Normal');
            mu(idxCell) = pd.mu;
            sig(idxCell) = pd.sigma;
            kurt(idxCell) = kurtosis(values);
            skew(idxCell) = skewness(values);
            
            % peak detection
            [~,locs,w,p]  = findpeaks( 	double(abs(timeSig(idxTrans).data(:,idxCell))/max(abs(timeSig(idxTrans).data(:,idxCell)))), ...
                                        R(:,idxCell), ...
                                        'WidthReference','halfheight');

            w = w(p>0.25);
            locs = locs(p>0.25);
            p = p(p>0.25);
            Npeaks(idxCell) = length(locs);

            if Npeaks(idxCell) > 1
                widthPeaks(idxCell)     = mean(w);
                promPeaks(idxCell)      = mean(p);
                stdWidthPeaks(idxCell)  = std(w);
                stdPromPeaks(idxCell)   = std(p);
            elseif Npeaks(idxCell) == 1
                widthPeaks(idxCell)     = mean(w);
                promPeaks(idxCell)      = mean(p);
                stdWidthPeaks(idxCell)  = std(w);
                stdPromPeaks(idxCell)   = std(p);
            else
                widthPeaks(idxCell)     = 0;
                promPeaks(idxCell)      = 0;
                stdWidthPeaks(idxCell)  = 0;
                stdPromPeaks(idxCell)   = 0;
            end
        end
    end
%     close(h)
    
    % % averaging features
    
    DataBuffers(1).data = SvValMat;
    DataBuffers(2).data = reshape(SNR, nWindows,nPings);
    DataBuffers(3).data = reshape(mu, nWindows,nPings);
    DataBuffers(4).data = reshape(sig, nWindows,nPings);
    DataBuffers(5).data = reshape(kurt, nWindows,nPings);
    DataBuffers(6).data = reshape(skew, nWindows,nPings);
    DataBuffers(7).data = reshape(Npeaks, nWindows,nPings);
    DataBuffers(8).data = reshape(widthPeaks, nWindows,nPings);
    DataBuffers(9).data = reshape(promPeaks, nWindows,nPings);
    DataBuffers(10).data = reshape(stdWidthPeaks, nWindows,nPings);
    DataBuffers(11).data = reshape(stdPromPeaks, nWindows,nPings);
    DataBuffersAvg      = performLocalAveraging(DataBuffers, maskMat, nAveragingCells, maxCellDist); % Default median smoothing
    
    
    % % % % % % % storing features
    SvMat(idxTrans,:)   = DataBuffersAvg(1).data(:);
    SNRMat(idxTrans,:)  = DataBuffersAvg(2).data(:);
    muMat(idxTrans,:)   = DataBuffersAvg(3).data(:);
    sigmaMat(idxTrans,:)= DataBuffersAvg(4).data(:);
    kurtMat(idxTrans,:) = DataBuffersAvg(5).data(:);
    skewMat(idxTrans,:) = DataBuffersAvg(6).data(:);
    NpeaksMat(idxTrans,:)           = DataBuffersAvg(7).data(:);
    widthPeaksMat(idxTrans,:)       = DataBuffersAvg(8).data(:);
    promPeaksMat(idxTrans,:)        = DataBuffersAvg(9).data(:);
    stdWidthPeaksMat(idxTrans,:)    = DataBuffersAvg(10).data(:);
    stdPromPeaksMat(idxTrans,:)     = DataBuffersAvg(11).data(:);
    
%     out.Buffers(idxTrans).SvVal             = SvVal;
%     out.Buffers(idxTrans).SNR               = SNR;
%     out.Buffers(idxTrans).mu                = mu;
%     out.Buffers(idxTrans).sigma             = sig;
%     out.Buffers(idxTrans).kurt              = kurt;
%     out.Buffers(idxTrans).skew              = skew;
end

%% saving

Features.Output.data = [    SvMat; ...              % 1:3
                            SNRMat; ...             % 4:6
                            muMat; ...              % 7:9
                            sigmaMat; ...           % 10:12
                            kurtMat; ...            % 13:15
                            skewMat; ...            % 16:18
                            NpeaksMat; ...          % 19:21
                            widthPeaksMat; ...      % 22:24
                            promPeaksMat; ...       % 25:27
                            stdPromPeaksMat; ...    % 28:31
                            stdWidthPeaksMat; ...   % 32:33
                            mean(SvMat); ...        % 34
                            mean(SNRMat); ...       % 35
                            mean(muMat); ...        % 36
                            mean(sigmaMat); ...     % 37
                            mean(kurtMat); ...      % 38
                            mean(skewMat); ...      % 39
                            mean(NpeaksMat); ...    % 40
                            mean(widthPeaksMat); ...% 41
                            mean(promPeaksMat); ... % 42
                            mean(stdPromPeaksMat); ...  %43
                            mean(stdWidthPeaksMat)];    % 44
                        
Features.Output.names = {   'SvMat1', ...
                            'SvMat2', ...
                            'SvMat3', ...
                            'SNRMat1', ...
                            'SNRMat2', ...
                            'SNRMat3', ...
                            'muMat1', ...
                            'muMat2', ...
                            'muMat3', ...
                            'sigmaMat1', ...
                            'sigmaMat2', ...
                            'sigmaMat3', ...
                            'kurtMat1', ...
                            'kurtMat2', ...
                            'kurtMat3', ...
                            'skewMat1', ...
                            'skewMat2', ...
                            'skewMat3', ...
                            'NpeaksMat1', ...
                            'NpeaksMat2', ...
                            'NpeaksMat3', ...
                            'widthPeaksMat1', ...
                            'widthPeaksMat2', ...
                            'widthPeaksMat3', ...
                            'promPeaksMat1', ...
                            'promPeaksMat2', ...
                            'promPeaksMat3', ...
                            'stdPromPeaksMat1', ...
                            'stdPromPeaksMat2', ...
                            'stdPromPeaksMat3', ...
                            'stdWidthPeaksMat1', ...
                            'stdWidthPeaksMat2', ...
                            'stdWidthPeaksMat3', ...
                            'SvMatMean', ...
                            'SNRMatMean', ...
                            'muMatMean', ...
                            'sigmaMatMean', ...
                            'kurtMatMean', ...
                            'skewMatMean', ...
                            'NpeaksMatMean', ...
                            'widthPeaksMatMean', ...
                            'promPeaksMatMean', ...
                            'stdPromPeaksMatMean', ...
                            'stdWidthPeaksMatMean'};

% timeFeaturesWMR.Buffers
% % % % % % % storing features
% out.Buffers(idxTrans).SvVal             = SvMat;
% out.Buffers(idxTrans).SNR               = SNRMat;
% out.Buffers(idxTrans).mu                = muMat;
% out.Buffers(idxTrans).sigma             = sigmaMat;
% out.Buffers(idxTrans).kurt              = kurtMat;
% out.Buffers(idxTrans).skew              = kurtMat;

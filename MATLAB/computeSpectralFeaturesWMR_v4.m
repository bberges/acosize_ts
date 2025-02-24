function Features = computeSpectralFeaturesWMR_v4(  ClassificationModel, ...
                                                    data, ...
                                                    timeSig, ...
                                                    range, ...
                                                    mask, ...
                                                    maxCellDist, ...
                                                    nAveragingCells)

nWindows = length(range)/size(timeSig(1).data,1);
nPings = size(timeSig(1).data,2)/nWindows;

maskMat = reshape(mask, nWindows,nPings);

newU            = ClassificationModel.newU;
df              = ClassificationModel.df;
df_low          = ClassificationModel.df_low;
interp_inter    = ClassificationModel.interp_inter;

S_f = data.data;
f_vec = data.freqs;

S_interp    = [];
f_interp    = [];
S_interp_low    = [];
f_interp_low    = [];
for idxInter = 1:length(interp_inter)
	% add 2kHz to make sure one includes the boundaries for the
	% interpolation
	idxCrop = find( f_vec >= interp_inter(idxInter,1)-df*2 & ...
                  	f_vec <= interp_inter(idxInter,2)+df*2);
                
	f_interp_temp = interp_inter(idxInter,1):df:interp_inter(idxInter,2);
	S_interp_temp = interp1(f_vec(idxCrop),S_f(idxCrop,:),f_interp_temp);
    
	f_interp_temp_low = interp_inter(idxInter,1):df_low:interp_inter(idxInter,2);
	S_interp_temp_low = interp1(f_vec(idxCrop),S_f(idxCrop,:),f_interp_temp_low);
                
	f_interp    = [f_interp f_interp_temp];
	S_interp    = [S_interp; S_interp_temp];
    
	f_interp_low    = [f_interp_low f_interp_temp_low];
	S_interp_low    = [S_interp_low; S_interp_temp_low];
end


%% compute features
polParam1   = zeros(2,size(S_f,2));
deltaMat    = zeros(1,size(S_f,2));
std_S_f     = zeros(1,size(S_f,2));
S_f_PCA_red	= zeros(size(newU,2),size(S_f,2));
S_low       = zeros(length(f_interp_low),size(S_f,2));

for idxCell = 1:size(S_f,2)
	if mask(idxCell) == 1
        
        S_interp_smooth(:,idxCell) = smooth(abs(S_interp(:,idxCell)),8);
        S_interp_smooth(:,idxCell) = (S_interp_smooth(:,idxCell) - mean(S_interp_smooth(:,idxCell)))/std(S_interp_smooth(:,idxCell));
            
        [polParam6,polyS,muPoly] = polyfit(   f_vec(:), ...
                                              10*log10(abs(S_f(:,idxCell))), ...
                                              6);
        [f_y,delta] = polyval(polParam6,f_vec,polyS,muPoly);
        deltaMat(idxCell) = mean(delta);

        [polParam1(:,idxCell),~,~] = polyfit( 	f_vec(:), ...
                                               	10*log10(abs(S_f(:,idxCell))), ...
                                              	1);
    %     [f_y,delta] = polyval(polParam2,f_vec,polyS2,muPoly2);
        % spectral feature 2:
        % spread of detrended signasl.
        % 1 value, std of detrend signal
        S_f_detrend = 10*log10(abs(S_f(:,idxCell))) - f_y(:); % Detrended data
        std_S_f(idxCell) = std(S_f_detrend);

        [polParam1(:,idxCell),~,~] = polyfit(    f_vec(:), ...
                                                            10*log10(abs(S_f(:,idxCell))), ...
                                                            1);

        S_f_detrend = 10*log10(abs(S_f(:,idxCell))) - f_y(:); % Detrended data
        std_S_f(idxCell) = std(S_f_detrend);

        % spectral feature 3:
        % spectrum transformed to PC axis
        S_f_PCA_red(:,idxCell) = newU'*S_interp_smooth(:,idxCell);
        
        idxRef = f_interp_low == ClassificationModel.ref_low;
        S_low(:,idxCell) = abs(S_interp_low(:,idxCell))/abs(S_interp_low(idxRef,idxCell));
    end
end

DataBuffers(1).data = reshape(std_S_f,nWindows,nPings);

count = 1;
for idx = 2:size(S_f_PCA_red,1)+1
    DataBuffers(idx).data = reshape(S_f_PCA_red(count,:),nWindows,nPings);
    count = count +1;
end

count = 1;
for idx = size(S_f_PCA_red,1)+2:size(S_f_PCA_red,1)+size(S_low,1)+1
    DataBuffers(idx).data = reshape(S_low(count,:),nWindows,nPings);
    count = count +1;
end

DataBuffersAvg      = performLocalAveraging(DataBuffers, maskMat, nAveragingCells, maxCellDist); % Default median smoothing

for idx = 1:size(S_f_PCA_red,1)
    S_f_PCA_red_av(idx,:) = DataBuffersAvg(idx+1).data(:);
end

for idx = 1:size(S_low,1)
    S_low_av(idx,:) = DataBuffersAvg(idx+size(S_f_PCA_red,1)+1).data(:);
end

%% saving

Features.Output.data = [    polParam1; ...  % 34:35
                            DataBuffersAvg(1).data(:)'; ...    % 36
                            S_f_PCA_red_av; ...  % 37:45
                            S_low_av];  % 46:62
                        
                        
Features.Output.names = {   'polParam1_1', ...
                            'polParam1_2', ...
                            'std_S_f', ...
                            'S_f_PCA_red1', ...
                            'S_f_PCA_red2', ...
                            'S_f_PCA_red3', ...
                            'S_f_PCA_red4', ...
                            'S_f_PCA_red5', ...
                            'S_f_PCA_red6', ...
                            'S_f_PCA_red7', ...
                            'S_f_PCA_red8', ...
                            'S_f_PCA_red9', ...
                            'S_f_PCA_red10', ...
                            strcat('S_low_',num2str(f_interp_low(1)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(2)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(3)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(4)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(5)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(6)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(7)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(8)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(9)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(10)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(11)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(12)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(13)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(14)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(15)*1e-3)), ...
                            strcat('S_low_',num2str(f_interp_low(16)*1e-3))};

% % % % % % % storing features
% spectralFeaturesWMR.polParam1      	= polParam1;
% spectralFeaturesWMR.delta           = delta;
% spectralFeaturesWMR.std_S_f       	= std_S_f;
% spectralFeaturesWMR.S_f_PCA_red   	= S_f_PCA_red;


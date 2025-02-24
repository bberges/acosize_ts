function Features = computeSpectralFeaturesWMR_v1(ClassificationModel, data, mask)


newU        = ClassificationModel.newU;
freqsPCA    = ClassificationModel.fPCA;

S_f = data.data;
f_vec = data.freqs;

polParam1   = zeros(2,size(S_f,2));
deltaMat       = zeros(1,size(S_f,2));
std_S_f     = zeros(1,size(S_f,2));
S_f_PCA_red	= zeros(size(newU,2),size(S_f,2));

for idxCell = 1:size(S_f,2)
	if mask(idxCell) == 1
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
        S_f_PCA = smooth(interp1(f_vec,abs(S_f(:,idxCell)),freqsPCA),10);
        S_f_PCA_red(:,idxCell) = newU'*S_f_PCA;
    end
end

Features.Output.data = [    polParam1; ...  % 34:35
                            std_S_f; ...    % 36
                            S_f_PCA_red];   % 37:45
                        
                        
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
                            'S_f_PCA_red9'};

% % % % % % % storing features
% spectralFeaturesWMR.polParam1      	= polParam1;
% spectralFeaturesWMR.delta           = delta;
% spectralFeaturesWMR.std_S_f       	= std_S_f;
% spectralFeaturesWMR.S_f_PCA_red   	= S_f_PCA_red;


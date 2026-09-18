%% Chapter 4 Results And Discussion Script
clear;close all;clc;
set(0,'DefaultAxesFontSize',13);
set(0,'DefaultLineLineWidth',1.35);
set(0,'DefaultFigureColor','w');
set(0,'DefaultFigurePosition',[100 100 1100 760]);

%% 0. GLOBAL CONFIGURATION
CFG.ChapterNumber=4;
CFG.OutDir=fullfile(pwd,'Chapter4_Figures');
if ~exist(CFG.OutDir,'dir');mkdir(CFG.OutDir);end
CFG.NumofAntenna=10;
CFG.Kd=pi;
CFG.ModuRadians=pi/4;
CFG.theta_tx=deg2rad(60);
CFG.theta_Int1=deg2rad(30);
CFG.theta_Int2=deg2rad(120);
CFG.SIR1_dB=3;
CFG.SIR2_dB=3;
CFG.lambda=0.75;
CFG.delta=1e-2;
CFG.EbNo_dB=-12:2:15;
CFG.SNR_dB=CFG.EbNo_dB+10*log10(2);
CFG.SNR_demo_dB=12;
CFG.BitsPerBlock=20000;
CFG.MinErrors=20;
CFG.MaxBlocks=15;
CFG.figIdx=0;
CFG.Style.FontSize=11;
CFG.Style.TitleSize=14;
CFG.Style.LineWidth=1.5;
logfile=fullfile(CFG.OutDir,'Figure_List_Chapter4.txt');
if exist(logfile,'file');delete(logfile);end
fprintf('=== Chapter 4 automated results generation ===\n');
fprintf('Output folder: %s\n\n',CFG.OutDir);

%% 1. FIGURE 4.1 - SYSTEM MODEL / ARRAY GEOMETRY
fig=figure('Name','System Model');
hold on;box on;grid on;
Nshow=CFG.NumofAntenna;
elementsX=(0:Nshow-1)*0.5;
origin=mean(elementsX);
arrowLen=3.0;
draw_doa_ray(origin,CFG.theta_tx,arrowLen,[0 0.55 0]);
draw_doa_ray(origin,CFG.theta_Int1,arrowLen,[0.75 0 0]);
draw_doa_ray(origin,CFG.theta_Int2,arrowLen,[0.75 0 0]);
plot(elementsX,zeros(1,Nshow),'ks','MarkerFaceColor','k','MarkerSize',9);
for k=1:Nshow
text(elementsX(k),-0.35,sprintf('%d',k-1),'HorizontalAlignment','center','FontSize',9);
end
labelR=arrowLen*1.18;
text(origin+labelR*cos(CFG.theta_tx),labelR*sin(CFG.theta_tx),...
sprintf('Desired User (%d^{\\circ})',round(rad2deg(CFG.theta_tx))),...
'Color',[0 0.5 0],'FontWeight','bold','FontSize',11,'HorizontalAlignment','center','VerticalAlignment','bottom');
text(origin+labelR*cos(CFG.theta_Int1),labelR*sin(CFG.theta_Int1),...
sprintf('Interferer 1 (%d^{\\circ})',round(rad2deg(CFG.theta_Int1))),...
'Color','r','FontWeight','bold','FontSize',11,'HorizontalAlignment','left','VerticalAlignment','bottom');
text(origin+labelR*cos(CFG.theta_Int2),labelR*sin(CFG.theta_Int2),...
sprintf('Interferer 2 (%d^{\\circ})',round(rad2deg(CFG.theta_Int2))),...
'Color','r','FontWeight','bold','FontSize',11,'HorizontalAlignment','right','VerticalAlignment','bottom');
text(origin,-0.85,'Antenna Element Index','HorizontalAlignment','center','FontSize',9,'Color',[0.3 0.3 0.3]);
axis equal;
xlim([-1.0 elementsX(end)+2.3]);ylim([-1.1 4.15]);
xlabel('Array Aperture (\lambda)');
title('Uniform Linear Array Geometry and Signal Directions of Arrival');
set(gca,'YTick',[]);
set(gca,'FontSize',CFG.Style.FontSize);
CFG=save_thesis_figure(fig,CFG,...
'Uniform linear array (N=10, d=\lambda/2) with the desired user and the two co-channel interferers at their respective angles of arrival.');

%% 2. RUN THE BASELINE (DEMO) RLS SIMULATION
rng(1);
bits_demo=max(CFG.BitsPerBlock*4,20000);
[Rx_Sig,QPSK_TxSig,uncoded_bits_Tx]=generate_scenario(...
CFG.NumofAntenna,CFG.Kd,CFG.theta_tx,CFG.theta_Int1,CFG.theta_Int2,...
CFG.SIR1_dB,CFG.SIR2_dB,CFG.SNR_demo_dB,bits_demo,CFG.ModuRadians);
[y_bf,e_bf,w_time,w_final]=rls_beamformer(Rx_Sig,QPSK_TxSig,...
CFG.NumofAntenna,CFG.lambda,CFG.delta);
y_ref=Rx_Sig(1,:);

%% 3. FIGURE 4.2 - CONSTELLATION: BEFORE vs AFTER RLS BEAMFORMING
fig=figure('Name','Constellation Before/After Beamforming');
border=2.5;
plotIdx=1:min(12000,length(QPSK_TxSig));
BER_ref_demo=mean(abs(uncoded_bits_Tx-qpsk_demod(y_ref)));
BER_bf_demo=mean(abs(uncoded_bits_Tx-qpsk_demod(y_bf)));
EVM_ref_pct=100*sqrt(mean(abs(y_ref-QPSK_TxSig).^2)/mean(abs(QPSK_TxSig).^2));
EVM_bf_pct=100*sqrt(mean(abs(y_bf-QPSK_TxSig).^2)/mean(abs(QPSK_TxSig).^2));
subplot(1,2,1);
plot(real(y_ref(plotIdx)),imag(y_ref(plotIdx)),'.','Color',[0.85 0.33 0.1],'MarkerSize',3);
axis(border*[-1 1-1 1]);axis square;grid on;box on;
line(border*[-1 1],[0 0],'Color','k');line([0 0],border*[-1 1],'Color','k');
title('Single Element: No Beamforming','FontSize',CFG.Style.TitleSize);
xlabel('In-Phase');ylabel('Quadrature');
text(-2.35,2.15,sprintf('BER = %.2e\nEVM = %.1f%%',BER_ref_demo,EVM_ref_pct),...
'FontSize',10,'BackgroundColor','w','EdgeColor',[0.75 0.75 0.75]);
subplot(1,2,2);
plot(real(y_bf(plotIdx)),imag(y_bf(plotIdx)),'.','Color',[0 0.45 0.74],'MarkerSize',3);
axis(border*[-1 1-1 1]);axis square;grid on;box on;
line(border*[-1 1],[0 0],'Color','k');line([0 0],border*[-1 1],'Color','k');
title('Array Output: RLS Beamforming','FontSize',CFG.Style.TitleSize);
xlabel('In-Phase');ylabel('Quadrature');
text(-2.35,2.15,sprintf('BER = %.2e\nEVM = %.1f%%',BER_bf_demo,EVM_bf_pct),...
'FontSize',10,'BackgroundColor','w','EdgeColor',[0.75 0.75 0.75]);
sgtitle_compat('Received Symbol Constellation: Impact of RLS Adaptive Beamforming');
CFG=save_thesis_figure(fig,CFG,...
'QPSK constellation observed at a single array element versus the RLS beamformer output at the same SNR, including BER and error-vector-magnitude indicators.');

%% 4. FIGURE 4.3 - RLS LEARNING CURVE BASED ON MEAN-SQUARE ERROR
NumOfSamples=min(1500,length(e_bf));
MSE_inst_dB=10*log10(abs(e_bf(1:NumOfSamples)).^2+eps);
MSE_smooth_dB=10*log10(moving_average(abs(e_bf(1:NumOfSamples)).^2,35)+eps);
tailStart=max(1,NumOfSamples-300);
steadyMSE_demo=mean(MSE_smooth_dB(tailStart:end));
convMSE_demo=find(MSE_smooth_dB<=steadyMSE_demo+1,1,'first');
if isempty(convMSE_demo);convMSE_demo=NumOfSamples;end
fig=figure('Name','RLS Learning Curve');
plot(MSE_inst_dB,'Color',[0.75 0.82 0.92]);hold on;grid on;box on;
plot(MSE_smooth_dB,'Color',[0 0.35 0.7],'LineWidth',1.8);
line(convMSE_demo*[1 1],ylim,'LineStyle','--','Color',[0.2 0.2 0.2]);
xlim([1 NumOfSamples]);
title('RLS Learning Curve from the A-Priori Error Signal');
xlabel('Sample Index (n)');ylabel('Mean-Square Error (dB)');
legend('Instantaneous |e(n)|^2','35-sample moving average',...
sprintf('Convergence: n = %d',convMSE_demo),'Location','best');
CFG=save_thesis_figure(fig,CFG,...
'A-priori RLS error shown as instantaneous and smoothed MSE, replacing the weak raw error trace with a standard learning-curve metric.');

%% 5. FIGURE 4.4 - ARRAY OUTPUT SINR IMPROVEMENT vs ITERATION NUMBER
a_tx=steering_vector(CFG.NumofAntenna,CFG.Kd,CFG.theta_tx);
a_i1=steering_vector(CFG.NumofAntenna,CFG.Kd,CFG.theta_Int1);
a_i2=steering_vector(CFG.NumofAntenna,CFG.Kd,CFG.theta_Int2);
Pint1=1/10^(CFG.SIR1_dB/10);
Pint2=1/10^(CFG.SIR2_dB/10);
N0=1/10^(CFG.SNR_demo_dB/10);
L_full=size(w_time,2);
idxSamples=unique(round(linspace(1,L_full,min(900,L_full))));
SINR_dB=zeros(1,length(idxSamples));
for ii=1:length(idxSamples)
wk=w_time(:,idxSamples(ii));
SINR_dB(ii)=output_sinr_db(wk,a_tx,a_i1,a_i2,Pint1,Pint2,N0);
end
SINR_smooth_dB=moving_average(SINR_dB,25);
w_ref=zeros(CFG.NumofAntenna,1);w_ref(1)=1;
SINR_ref_dB=output_sinr_db(w_ref,a_tx,a_i1,a_i2,Pint1,Pint2,N0);
fig=figure('Name','SINR Improvement');
plot(idxSamples,SINR_dB,'Color',[0.70 0.86 0.72]);hold on;grid on;box on;
plot(idxSamples,SINR_smooth_dB,'Color',[0.1 0.5 0.2],'LineWidth',1.8);
line([idxSamples(1) idxSamples(end)],SINR_ref_dB*[1 1],'LineStyle','--','Color',[0.55 0.1 0.1]);
xlim([1 min(5000,L_full)]);
title('Array Output SINR Convergence with RLS Beamforming');
xlabel('Sample Index (n)');ylabel('Output SINR (dB)');
legend('Instantaneous RLS weights','Smoothed RLS trend',...
sprintf('Single-element baseline (%.1f dB)',SINR_ref_dB),'Location','best');
CFG=save_thesis_figure(fig,CFG,...
'Array output SINR computed from the evolving RLS weight vector and compared with the single-element receiver baseline.');

%% 6. FIGURE 4.5 - WEIGHT-VECTOR NORM AND UPDATE MAGNITUDE
NumOfSamples2=min(1500,size(w_time,2));
w_norm=sqrt(sum(abs(w_time(:,1:NumOfSamples2)).^2,1));
w_update=[0,sqrt(sum(abs(diff(w_time(:,1:NumOfSamples2),1,2)).^2,1))];
fig=figure('Name','Weight Vector Stability');
subplot(2,1,1);
plot(moving_average(w_norm,25),'Color',[0 0.45 0.74]);grid on;box on;
xlim([1 NumOfSamples2]);
title('RLS Weight-Vector Norm');
ylabel('||w(n)||_2');
subplot(2,1,2);
semilogy(moving_average(w_update+eps,25),'Color',[0.64 0.08 0.18]);grid on;box on;
xlim([1 NumOfSamples2]);
title('RLS Weight-Update Magnitude');
xlabel('Sample Index (n)');ylabel('||w(n)-w(n-1)||_2');
CFG=save_thesis_figure(fig,CFG,...
'Weight-vector norm and update magnitude, showing adaptation stability without relying on visually noisy individual coefficient phase traces.');

%% 7. FIGURE 4.6 - ERROR-VECTOR MAGNITUDE DISTRIBUTION
err_ref=abs(y_ref-QPSK_TxSig);
err_bf=abs(y_bf-QPSK_TxSig);
[cdf_ref_x,cdf_ref_y]=empirical_cdf(err_ref);
[cdf_bf_x,cdf_bf_y]=empirical_cdf(err_bf);
fig=figure('Name','Error Vector Distribution');
plot(cdf_ref_x,cdf_ref_y,'Color',[0.85 0.33 0.1]);hold on;grid on;box on;
plot(cdf_bf_x,cdf_bf_y,'Color',[0 0.45 0.74],'LineWidth',1.8);
xlim([0 prctile_compat([err_ref err_bf],99.5)]);
ylim([0 1]);
title('Empirical CDF of Symbol Error-Vector Magnitude');
xlabel('|y(n)-s(n)|');ylabel('Cumulative Probability');
legend(sprintf('No beamforming, EVM = %.1f%%',EVM_ref_pct),...
sprintf('RLS beamforming, EVM = %.1f%%',EVM_bf_pct),'Location','southeast');
CFG=save_thesis_figure(fig,CFG,...
'Empirical distribution of the symbol error-vector magnitude before and after RLS beamforming, quantifying constellation tightening.');

%% 8. FIGURES 4.7 - 4.9 - Rx ARRAY RADIATION PATTERN AND DIRECTIONAL GAINS
theta_scan=linspace(0,180,721);
AF=exp(1j*CFG.Kd*(0:CFG.NumofAntenna-1)'*cos(deg2rad(theta_scan)));
RadPattern=abs(w_final.'*AF);
RadPattern_dB=20*log10(RadPattern/max(RadPattern));
[~,idxTx]=min(abs(theta_scan-rad2deg(CFG.theta_tx)));
[~,idxI1]=min(abs(theta_scan-rad2deg(CFG.theta_Int1)));
[~,idxI2]=min(abs(theta_scan-rad2deg(CFG.theta_Int2)));
dirGain_dB=[RadPattern_dB(idxTx),RadPattern_dB(idxI1),RadPattern_dB(idxI2)];
fig=figure('Name','Radiation Pattern - Cartesian');
plot(theta_scan,RadPattern_dB);grid on;hold on;box on;
xlim([0 180]);ylim([-40 5]);
yl=[-40 5];
line(rad2deg(CFG.theta_tx)*[1 1],yl,'LineStyle','--','Color',[0 0.6 0]);
line(rad2deg(CFG.theta_Int1)*[1 1],yl,'LineStyle','--','Color','r');
line(rad2deg(CFG.theta_Int2)*[1 1],yl,'LineStyle','--','Color','r');
text(rad2deg(CFG.theta_Int1)+2,dirGain_dB(2)-2,sprintf('%.1f dB',dirGain_dB(2)),...
'Color','r','FontSize',10,'BackgroundColor','w');
text(rad2deg(CFG.theta_Int2)+2,dirGain_dB(3)-2,sprintf('%.1f dB',dirGain_dB(3)),...
'Color','r','FontSize',10,'BackgroundColor','w');
title('Receive Array Gain Pattern - RLS Steady-State Weights');
xlabel('Angle (degrees)');ylabel('Normalized Gain (dB)');
legend('Array Pattern',...
sprintf('Desired User @ %d^{\\circ}',round(rad2deg(CFG.theta_tx))),...
sprintf('Interferer 1 @ %d^{\\circ}',round(rad2deg(CFG.theta_Int1))),...
sprintf('Interferer 2 @ %d^{\\circ}',round(rad2deg(CFG.theta_Int2))),...
'Location','best');
CFG=save_thesis_figure(fig,CFG,...
'Receive array gain pattern formed by the converged RLS weight vector, showing the main beam steered toward the desired user and nulls placed toward the two interferers.');
fig=figure('Name','Radiation Pattern - Polar');
pax=polar_pattern_plot(deg2rad(theta_scan),RadPattern/max(RadPattern),...
CFG.theta_tx,CFG.theta_Int1,CFG.theta_Int2);
title('Receive Array Polar Radiation Pattern - RLS Steady-State Weights');
CFG=save_thesis_figure(fig,CFG,...
'Polar representation of the receive array radiation pattern after RLS adaptation.');

%% 9. FIGURE 4.9 - DIRECTIONAL GAIN SUMMARY
fig=figure('Name','Directional Gain Summary');
b=bar(dirGain_dB,'FaceColor','flat');grid on;box on;
b.CData=[0 0.55 0;0.8 0.1 0.1;0.8 0.1 0.1];
set(gca,'XTickLabel',{'Desired 60 deg','Interferer 30 deg','Interferer 120 deg'});
ylabel('Normalized Array Gain (dB)');
title('Directional Response of the Converged RLS Beamformer');
ylim([min(dirGain_dB)-8 3]);
for kk=1:numel(dirGain_dB)
text(kk,dirGain_dB(kk)+0.7,sprintf('%.1f dB',dirGain_dB(kk)),...
'HorizontalAlignment','center','FontSize',10);
end
CFG=save_thesis_figure(fig,CFG,...
'Normalized gain at the desired-user direction and at the two interferer directions, summarizing main-beam preservation and adaptive null formation.');

%% 10. BER SWEEP OVER SNR (WITH & WITHOUT BEAMFORMING)
fprintf('\n--- Running BER sweep over SNR (Fig. 4.10 / 4.11) ---\n');
nSNR=length(CFG.SNR_dB);
BER_bf=nan(1,nSNR);
BER_nobf=nan(1,nSNR);
rng(2);
for si=1:nSNR
[BER_bf(si),BER_nobf(si),~]=ber_monte_carlo(...
CFG.NumofAntenna,CFG.Kd,CFG.theta_tx,CFG.theta_Int1,CFG.theta_Int2,...
CFG.SIR1_dB,CFG.SIR2_dB,CFG.SNR_dB(si),CFG.lambda,CFG.delta,...
CFG.ModuRadians,CFG.BitsPerBlock,CFG.MinErrors,CFG.MaxBlocks);
fprintf('  SNR = %6.2f dB | BER(RLS) = %.3e | BER(no BF) = %.3e\n',...
CFG.SNR_dB(si),BER_bf(si),BER_nobf(si));
end

%% 11. FIGURE 4.10 - BER vs SNR: RLS BEAMFORMING vs THEORETICAL BOUNDS
EbNo_th=-12:1:20;
SNR_th=EbNo_th+10*log10(2);
EbN0Lin=10.^(EbNo_th/10);
BER_Rayleigh=0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));
BER_AWGN=0.5*erfc(sqrt(EbN0Lin));
fig=figure('Name','BER vs SNR - Theoretical Bounds');
semilogy(CFG.SNR_dB,BER_bf,'o-','Color',[0 0.45 0.74]);hold on;grid on;box on;
semilogy(SNR_th,BER_Rayleigh,'--s','Color',[0.85 0.33 0.1],'MarkerSize',4);
semilogy(SNR_th,BER_AWGN,'-','Color',[0.2 0.6 0.2]);
xlim([CFG.SNR_dB(1) CFG.SNR_dB(end)]);ylim([1e-6 1]);
legend('RLS Adaptive Beamforming (Simulated)','Theoretical Rayleigh Fading',...
'Theoretical AWGN','Location','southwest');
title('BER Performance of QPSK with RLS Adaptive Beamforming');
xlabel('SNR (dB)');ylabel('Bit Error Rate');
CFG=save_thesis_figure(fig,CFG,...
'Simulated BER of the RLS-beamformed QPSK link versus SNR, compared against the theoretical AWGN and Rayleigh-fading bounds for an equivalent single-branch receiver.');

%% 12. FIGURE 4.11 - BER: WITH vs WITHOUT BEAMFORMING
fig=figure('Name','BER With vs Without Beamforming');
semilogy(CFG.SNR_dB,BER_bf,'o-','Color',[0 0.45 0.74]);hold on;grid on;box on;
semilogy(CFG.SNR_dB,BER_nobf,'s--','Color',[0.64 0.08 0.18]);
xlim([CFG.SNR_dB(1) CFG.SNR_dB(end)]);ylim([1e-6 1]);
legend('With RLS Beamforming (N=10)','Without Beamforming (Single Element)',...
'Location','southwest');
title('Effect of RLS Adaptive Beamforming on BER Performance');
xlabel('SNR (dB)');ylabel('Bit Error Rate');
CFG=save_thesis_figure(fig,CFG,...
'BER comparison between the RLS-beamformed array output and a single reference element with no interference suppression, quantifying the benefit of adaptive beamforming.');

%% 13. EFFECT OF NUMBER OF ANTENNAS (used by Fig. 4.12 / 4.14 / Table 4.3)
fprintf('\n--- Running antenna-count sweep (Fig. 4.12 / 4.14 / Table 4.3) ---\n');
N_list=[4 6 8 10 12];
BER_vs_N=nan(1,length(N_list));
OutputSINR_vs_N_dB=nan(1,length(N_list));
NullDepth1_dB=nan(1,length(N_list));
NullDepth2_dB=nan(1,length(N_list));
HPBW_vs_N_deg=nan(1,length(N_list));
RP_all=cell(1,length(N_list));
rng(3);
for kk=1:length(N_list)
Nk=N_list(kk);
[BER_vs_N(kk),~,~]=ber_monte_carlo(Nk,CFG.Kd,CFG.theta_tx,CFG.theta_Int1,CFG.theta_Int2,...
CFG.SIR1_dB,CFG.SIR2_dB,CFG.SNR_demo_dB,CFG.lambda,CFG.delta,...
CFG.ModuRadians,CFG.BitsPerBlock,CFG.MinErrors,CFG.MaxBlocks);
bits_k=max(CFG.BitsPerBlock*2,8000);
[RxK,TxK]=generate_scenario(Nk,CFG.Kd,CFG.theta_tx,CFG.theta_Int1,CFG.theta_Int2,...
CFG.SIR1_dB,CFG.SIR2_dB,CFG.SNR_demo_dB,bits_k,CFG.ModuRadians);
[~,~,~,wk_final]=rls_beamformer(RxK,TxK,Nk,CFG.lambda,CFG.delta);
AFk=exp(1j*CFG.Kd*(0:Nk-1)'*cos(deg2rad(theta_scan)));
RPk=abs(wk_final.'*AFk);
RPk_dB=20*log10(RPk/max(RPk));
a_tx_k=steering_vector(Nk,CFG.Kd,CFG.theta_tx);
a_i1_k=steering_vector(Nk,CFG.Kd,CFG.theta_Int1);
a_i2_k=steering_vector(Nk,CFG.Kd,CFG.theta_Int2);
OutputSINR_vs_N_dB(kk)=output_sinr_db(wk_final,a_tx_k,a_i1_k,a_i2_k,Pint1,Pint2,N0);
[~,iTx]=min(abs(theta_scan-rad2deg(CFG.theta_tx)));
[~,iInt1]=min(abs(theta_scan-rad2deg(CFG.theta_Int1)));
[~,iInt2]=min(abs(theta_scan-rad2deg(CFG.theta_Int2)));
NullDepth1_dB(kk)=RPk_dB(iInt1);
NullDepth2_dB(kk)=RPk_dB(iInt2);
HPBW_vs_N_deg(kk)=half_power_beamwidth(theta_scan,RPk_dB,iTx);
if Nk==CFG.NumofAntenna
RPk_dB_ref=RPk_dB;
end
RP_all{kk}=RPk_dB;
end
fig=figure('Name','Array Size Performance Metrics');
subplot(3,1,1);
plot(N_list,OutputSINR_vs_N_dB,'o-','Color',[0 0.45 0.74]);grid on;box on;
xlim([N_list(1)-0.5 N_list(end)+0.5]);
ylabel('SINR (dB)');
title(sprintf('Effect of Array Size on RLS Beamforming Metrics (SNR = %.1f dB)',CFG.SNR_demo_dB));
subplot(3,1,2);
plot(N_list,HPBW_vs_N_deg,'s-','Color',[0.49 0.18 0.56]);grid on;box on;
xlim([N_list(1)-0.5 N_list(end)+0.5]);
ylabel('HPBW (deg)');
subplot(3,1,3);
plot(N_list,NullDepth1_dB,'^-','Color',[0.85 0.33 0.1]);hold on;
plot(N_list,NullDepth2_dB,'v--','Color',[0.64 0.08 0.18]);grid on;box on;
xlim([N_list(1)-0.5 N_list(end)+0.5]);
xlabel('Number of Antenna Elements (N)');ylabel('Null Depth (dB)');
legend('Interferer 1 at 30 deg','Interferer 2 at 120 deg','Location','best');
CFG=save_thesis_figure(fig,CFG,...
'Effect of array size on output SINR, half-power beamwidth, and interferer null depth, replacing an unstable single-SNR BER-only plot.');

%% 14. EFFECT OF THE RLS FORGETTING FACTOR (Fig. 4.13 / Table 4.4)
fprintf('\n--- Running forgetting-factor sweep (Fig. 4.13 / Table 4.4) ---\n');
lambda_list=[0.90 0.95 0.98 0.995];
bits_lambda=max(CFG.BitsPerBlock*4,20000);
MSE_curves=cell(1,length(lambda_list));
SteadyMSE_dB=nan(1,length(lambda_list));
ConvSamples=nan(1,length(lambda_list));
rng(4);
[RxL,TxL]=generate_scenario(CFG.NumofAntenna,CFG.Kd,CFG.theta_tx,CFG.theta_Int1,CFG.theta_Int2,...
CFG.SIR1_dB,CFG.SIR2_dB,CFG.SNR_demo_dB,bits_lambda,CFG.ModuRadians);
for kk=1:length(lambda_list)
[~,e_l,~,~]=rls_beamformer(RxL,TxL,CFG.NumofAntenna,lambda_list(kk),CFG.delta);
inst_MSE_dB=10*log10(moving_average(abs(e_l).^2,25)+eps);
MSE_curves{kk}=inst_MSE_dB;
SteadyMSE_dB(kk)=mean(inst_MSE_dB(end-min(500,floor(length(inst_MSE_dB)/4))+1:end));
thresh=SteadyMSE_dB(kk)+1;
convIdx=find(inst_MSE_dB<=thresh,1,'first');
if isempty(convIdx);convIdx=length(inst_MSE_dB);end
ConvSamples(kk)=convIdx;
end
fig=figure('Name','Effect of Forgetting Factor');
subplot(2,1,1);hold on;grid on;box on;
colorsL=lines(length(lambda_list));
legL=cell(1,length(lambda_list));
plotLen=min(800,length(MSE_curves{1}));
for kk=1:length(lambda_list)
plot(MSE_curves{kk}(1:plotLen),'Color',colorsL(kk,:));
legL{kk}=sprintf('\\lambda = %.3f',lambda_list(kk));
end
title('Effect of the RLS Forgetting Factor on Convergence Speed');
xlabel('Sample Index (n)');ylabel('Smoothed MSE (dB)');
legend(legL,'Location','best');
subplot(2,1,2);
yyaxis left;
bar(lambda_list,SteadyMSE_dB,0.45,'FaceColor',[0 0.45 0.74]);
ylabel('Steady-State MSE (dB)');
yyaxis right;
plot(lambda_list,ConvSamples,'s-','Color',[0.64 0.08 0.18],'LineWidth',1.6);
ylabel('Convergence Samples');
grid on;box on;
xlabel('Forgetting Factor \lambda');
title('Steady-State Error and Convergence-Speed Summary');
CFG=save_thesis_figure(fig,CFG,...
'Smoothed learning curves of the RLS algorithm for several values of the forgetting factor \lambda, illustrating the convergence-speed/steady-state-error trade-off.');

%% 15. FIGURE 4.14 - RADIATION PATTERN COMPARISON FOR DIFFERENT N
fig=figure('Name','Radiation Pattern vs N');
hold on;grid on;box on;
colorsN=lines(length(N_list));
legN=cell(1,length(N_list));
for kk=1:length(N_list)
plot(theta_scan,RP_all{kk},'Color',colorsN(kk,:));
legN{kk}=sprintf('N = %d',N_list(kk));
end
xlim([0 180]);ylim([-40 5]);
yl=[-40 5];
line(rad2deg(CFG.theta_tx)*[1 1],yl,'LineStyle',':','Color','k');
line(rad2deg(CFG.theta_Int1)*[1 1],yl,'LineStyle',':','Color',[0.6 0 0]);
line(rad2deg(CFG.theta_Int2)*[1 1],yl,'LineStyle',':','Color',[0.6 0 0]);
title('Radiation Pattern Comparison for Different Array Sizes');
xlabel('Angle (degrees)');ylabel('Normalized Gain (dB)');
legend(legN,'Location','best');
CFG=save_thesis_figure(fig,CFG,...
'RLS steady-state radiation patterns for different numbers of array elements, showing the improvement in beamwidth and null depth with increasing array size.');

%% 16. FIGURE 4.15 - EFFECT OF SIR ON BER PERFORMANCE
fprintf('\n--- Running SIR-dependent BER curves (Fig. 4.15) ---\n');
SIR_curve_list=[0 3 10];
sirSNR_idx=1:2:length(CFG.SNR_dB);
SNR_sir_dB=CFG.SNR_dB(sirSNR_idx);
BER_vs_SIR_curve=nan(length(SIR_curve_list),length(SNR_sir_dB));
rng(5);
for kk=1:length(SIR_curve_list)
for si=1:length(SNR_sir_dB)
[BER_vs_SIR_curve(kk,si),~,~]=ber_monte_carlo(CFG.NumofAntenna,CFG.Kd,CFG.theta_tx,...
CFG.theta_Int1,CFG.theta_Int2,SIR_curve_list(kk),SIR_curve_list(kk),SNR_sir_dB(si),...
CFG.lambda,CFG.delta,CFG.ModuRadians,CFG.BitsPerBlock,CFG.MinErrors,CFG.MaxBlocks);
end
end
fig=figure('Name','BER vs SNR for SIR Levels');
hold on;grid on;box on;
colorsSIR=lines(length(SIR_curve_list));
legSIR=cell(1,length(SIR_curve_list));
for kk=1:length(SIR_curve_list)
semilogy(SNR_sir_dB,BER_vs_SIR_curve(kk,:),'o-','Color',colorsSIR(kk,:));
legSIR{kk}=sprintf('SIR_1 = SIR_2 = %d dB',SIR_curve_list(kk));
end
xlim([SNR_sir_dB(1) SNR_sir_dB(end)]);ylim([1e-6 1]);
title('BER Sensitivity to Co-Channel Interference Level');
xlabel('SNR (dB)');ylabel('Bit Error Rate');
legend(legSIR,'Location','southwest');
CFG=save_thesis_figure(fig,CFG,...
'BER curves across SNR for multiple signal-to-interference ratios, replacing the weak single-SNR SIR plot with a full sensitivity comparison.');

%% 16. FIGURE 4.16 - LEARNING CURVE COMPARISON: LMS vs NLMS vs RLS
fprintf('\n--- Running algorithm comparison (LMS vs NLMS vs RLS) ---\n');
mu_LMS=0.002;
mu_NLMS=0.1;
[~,e_lms,~,~]=lms_beamformer(Rx_Sig,QPSK_TxSig,CFG.NumofAntenna,mu_LMS);
[~,e_nlms,~,~]=nlms_beamformer(Rx_Sig,QPSK_TxSig,CFG.NumofAntenna,mu_NLMS);
MSE_LMS_dB=10*log10(moving_average(abs(e_lms).^2,35)+eps);
MSE_NLMS_dB=10*log10(moving_average(abs(e_nlms).^2,35)+eps);
MSE_RLS_dB=MSE_smooth_dB;
fig=figure('Name','Learning Curves: LMS vs NLMS vs RLS');
plot(MSE_LMS_dB(1:NumOfSamples),'Color',[0.85 0.33 0.1],'LineWidth',1.5);hold on;grid on;box on;
plot(MSE_NLMS_dB(1:NumOfSamples),'Color',[0.49 0.18 0.56],'LineWidth',1.5);
plot(MSE_RLS_dB(1:NumOfSamples),'Color',[0 0.45 0.74],'LineWidth',1.5);
xlim([1 NumOfSamples]);
title('Convergence Speed Comparison: LMS vs NLMS vs RLS');
xlabel('Sample Index (n)');ylabel('Smoothed MSE (dB)');
legend('LMS (\mu = 0.002)','NLMS (\mu = 0.1)',sprintf('RLS (\\lambda = %.2f)',CFG.lambda),'Location','northeast');
CFG=save_thesis_figure(fig,CFG,...
'Smoothed learning curves comparing the convergence speed and steady-state error of LMS, NLMS, and RLS algorithms.');
steadyMSE_lms=mean(MSE_LMS_dB(tailStart:end));
steadyMSE_nlms=mean(MSE_NLMS_dB(tailStart:end));
conv_lms=find(MSE_LMS_dB<=steadyMSE_lms+1,1,'first');
conv_nlms=find(MSE_NLMS_dB<=steadyMSE_nlms+1,1,'first');
if isempty(conv_lms);conv_lms=NumOfSamples;end
if isempty(conv_nlms);conv_nlms=NumOfSamples;end

%% 17. FIGURE 4.17 - BER COMPARISON: LMS vs NLMS vs RLS
BER_lms=nan(1,nSNR);
BER_nlms=nan(1,nSNR);
rng(6);
for si=1:nSNR
T_err_lms=0;T_err_nlms=0;blocks=0;
while (T_err_lms<CFG.MinErrors||T_err_nlms<CFG.MinErrors)&&(blocks<CFG.MaxBlocks)
blocks=blocks+1;
[Rx_S,QPSK_Tx,uncoded_Tx]=generate_scenario(...
CFG.NumofAntenna,CFG.Kd,CFG.theta_tx,CFG.theta_Int1,CFG.theta_Int2,...
CFG.SIR1_dB,CFG.SIR2_dB,CFG.SNR_dB(si),CFG.BitsPerBlock,CFG.ModuRadians);
[y_lms,~,~,~]=lms_beamformer(Rx_S,QPSK_Tx,CFG.NumofAntenna,mu_LMS);
[y_nlms,~,~,~]=nlms_beamformer(Rx_S,QPSK_Tx,CFG.NumofAntenna,mu_NLMS);
bits_rx_lms=qpsk_demod(y_lms);
bits_rx_nlms=qpsk_demod(y_nlms);
T_err_lms=T_err_lms+sum(abs(uncoded_Tx-bits_rx_lms));
T_err_nlms=T_err_nlms+sum(abs(uncoded_Tx-bits_rx_nlms));
end
BER_lms(si)=T_err_lms/(blocks*length(uncoded_Tx));
BER_nlms(si)=T_err_nlms/(blocks*length(uncoded_Tx));
end
fig=figure('Name','BER: LMS vs NLMS vs RLS');
semilogy(CFG.SNR_dB,BER_lms,'s-','Color',[0.85 0.33 0.1]);hold on;grid on;box on;
semilogy(CFG.SNR_dB,BER_nlms,'^-','Color',[0.49 0.18 0.56]);
semilogy(CFG.SNR_dB,BER_bf,'o-','Color',[0 0.45 0.74]);
xlim([CFG.SNR_dB(1) CFG.SNR_dB(end)]);ylim([1e-6 1]);
legend('LMS Beamforming','NLMS Beamforming','RLS Beamforming','Location','southwest');
title('BER Performance Comparison of Adaptive Algorithms');
xlabel('SNR (dB)');ylabel('Bit Error Rate');
CFG=save_thesis_figure(fig,CFG,...
'BER comparison between LMS, NLMS, and RLS beamforming, evaluating overall performance.');
fprintf('\n=== All 17 figures generated (41.png ... 417.png) ===\n');

%% Tables Exported As Figures And Data Files

%% TABLE 4.1 - SIMULATION PARAMETERS SUMMARY
T1_headers={'Parameter','Value'};
T1_rows={
'Adaptive algorithm','Recursive Least Squares (RLS)';
'Array configuration','Uniform Linear Array (ULA)';
'Number of antenna elements, N',num2str(CFG.NumofAntenna);
'Element spacing','\lambda / 2';
'Modulation','QPSK (Gray-coded)';
'Desired user DoA',sprintf('%d deg',round(rad2deg(CFG.theta_tx)));
'Interferer 1 DoA / SIR_1',sprintf('%d deg / %.1f dB',round(rad2deg(CFG.theta_Int1)),CFG.SIR1_dB);
'Interferer 2 DoA / SIR_2',sprintf('%d deg / %.1f dB',round(rad2deg(CFG.theta_Int2)),CFG.SIR2_dB);
'RLS forgetting factor, \lambda',num2str(CFG.lambda);
'RLS regularization, \delta',num2str(CFG.delta);
'SNR sweep range',sprintf('%.1f to %.1f dB',CFG.SNR_dB(1),CFG.SNR_dB(end));
};
T1=cell2table(T1_rows,'VariableNames',{'Parameter','Value'});
CFG=save_thesis_table(T1,T1_headers,T1_rows,CFG,...
'Table 4.1','Simulation Parameters',...
'Summary of the system-model and RLS algorithm parameters used throughout the simulations of this chapter.');

%% TABLE 4.2 - BER RESULTS SUMMARY (WITH vs WITHOUT BEAMFORMING)
T2_headers={'SNR (dB)','BER - With RLS BF','BER - Without BF','Improvement Factor'};
T2_rows=cell(nSNR,4);
for si=1:nSNR
T2_rows{si,1}=sprintf('%.2f',CFG.SNR_dB(si));
T2_rows{si,2}=sprintf('%.3e',BER_bf(si));
T2_rows{si,3}=sprintf('%.3e',BER_nobf(si));
T2_rows{si,4}=sprintf('%.1f x',(BER_nobf(si)+eps)/(BER_bf(si)+eps));
end
T2=cell2table(T2_rows,'VariableNames',{'SNR_dB','BER_With_BF','BER_Without_BF','Improvement_Factor'});
CFG=save_thesis_table(T2,T2_headers,T2_rows,CFG,...
'Table 4.2','BER Performance Summary',...
'Bit error rate versus SNR with and without RLS beamforming, and the resulting improvement factor.');

%% TABLE 4.3 - EFFECT OF NUMBER OF ANTENNAS
T3_headers={'N','Output SINR (dB)','HPBW (deg)','Null @ 30 deg (dB)','Null @ 120 deg (dB)'};
T3_rows=cell(length(N_list),5);
for kk=1:length(N_list)
T3_rows{kk,1}=num2str(N_list(kk));
T3_rows{kk,2}=sprintf('%.2f',OutputSINR_vs_N_dB(kk));
T3_rows{kk,3}=sprintf('%.2f',HPBW_vs_N_deg(kk));
T3_rows{kk,4}=sprintf('%.2f',NullDepth1_dB(kk));
T3_rows{kk,5}=sprintf('%.2f',NullDepth2_dB(kk));
end
T3=cell2table(T3_rows,'VariableNames',{'N_elements','OutputSINR_dB','HPBW_deg','NullDepth30_dB','NullDepth120_dB'});
CFG=save_thesis_table(T3,T3_headers,T3_rows,CFG,...
'Table 4.3','Effect of Array Size',...
'Output SINR, half-power beamwidth, and interferer null depth as a function of the number of antenna elements.');

%% TABLE 4.4 - EFFECT OF RLS FORGETTING FACTOR
T4_headers={'Forgetting Factor \lambda','Steady-State MSE (dB)','Conv. Samples'};
T4_rows=cell(length(lambda_list),3);
for kk=1:length(lambda_list)
T4_rows{kk,1}=sprintf('%.3f',lambda_list(kk));
T4_rows{kk,2}=sprintf('%.2f',SteadyMSE_dB(kk));
T4_rows{kk,3}=num2str(ConvSamples(kk));
end
T4=cell2table(T4_rows,'VariableNames',{'Lambda','SteadyStateMSE_dB','ConvergenceSamples'});
CFG=save_thesis_table(T4,T4_headers,T4_rows,CFG,...
'Table 4.4','Effect of the RLS Forgetting Factor',...
'Steady-state MSE and approximate convergence speed for several values of the RLS forgetting factor.');

%% TABLE 4.5 - ALGORITHM COMPARISON SUMMARY
T5_headers={'Algorithm','Steady-State MSE (dB)','Convergence Samples'};
T5_rows={
'LMS',sprintf('%.2f',steadyMSE_lms),num2str(conv_lms);
'NLMS',sprintf('%.2f',steadyMSE_nlms),num2str(conv_nlms);
'RLS',sprintf('%.2f',steadyMSE_demo),num2str(convMSE_demo);
};
T5=cell2table(T5_rows,'VariableNames',{'Algorithm','SteadyStateMSE_dB','ConvergenceSamples'});
CFG=save_thesis_table(T5,T5_headers,T5_rows,CFG,...
'Table 4.5','Adaptive Algorithms Comparison',...
'Comparison of convergence speed and steady-state MSE for LMS, NLMS, and RLS.');
fprintf('\n=== All 5 tables generated (as PNG + CSV + XLSX) ===\n');
fprintf('=== Chapter 4 generation complete. See "%s" ===\n',CFG.OutDir);

%% Local Functions
function s=qpsk_gray_modulate(b1,b2,moduRad)
s=((b1==0).*(b2==0)*exp(1i*moduRad)+...
(b1==0).*(b2==1)*exp(3i*moduRad)+...
(b1==1).*(b2==1)*exp(5i*moduRad)+...
(b1==1).*(b2==0)*exp(7i*moduRad));
end

%% QPSK Demodulation Helper
function bits_rx=qpsk_demod(y)
B4=(real(y)<0);
B3=(imag(y)<0);
bits_rx=zeros(1,2*length(y));
bits_rx(1:2:end)=B3;
bits_rx(2:2:end)=B4;
end

%% ULA Steering Vector Helper
function v=steering_vector(N,kd,theta)
v=exp(1i*kd*(0:N-1)'*cos(theta));
end

%% Output SINR Calculation Helper
function sinr_dB=output_sinr_db(w,a_tx,a_i1,a_i2,Pint1,Pint2,N0)
Sdes=abs(w.'*a_tx)^2;
Iint=Pint1*abs(w.'*a_i1)^2 + Pint2*abs(w.'*a_i2)^2;
Nout=N0*(w'*w);
sinr_dB=10*log10(Sdes/(Iint+Nout+eps));
end

%% Scenario Generation Helper
function [Rx_Sig,QPSK_TxSig,uncoded_bits_Tx]=generate_scenario(...
N,kd,theta_tx,theta_i1,theta_i2,SIR1_dB,SIR2_dB,SNR_dB,bit_count,moduRad)
uncoded_bits_Tx=round(rand(1,bit_count));
B1=uncoded_bits_Tx(1:2:end);B2=uncoded_bits_Tx(2:2:end);
QPSK_TxSig=qpsk_gray_modulate(B1,B2,moduRad);
bits1=round(rand(1,bit_count));
B1=bits1(1:2:end);B2=bits1(2:2:end);
Inter1=1/10^(SIR1_dB/10);
QPSK_Int1=sqrt(Inter1/2)*qpsk_gray_modulate(B1,B2,moduRad);
bits2=round(rand(1,bit_count));
B1=bits2(1:2:end);B2=bits2(2:2:end);
Inter2=1/10^(SIR2_dB/10);
QPSK_Int2=sqrt(Inter2/2)*qpsk_gray_modulate(B1,B2,moduRad);
a_tx=steering_vector(N,kd,theta_tx);
a_i1=steering_vector(N,kd,theta_i1);
a_i2=steering_vector(N,kd,theta_i2);
Tx_Sig=a_tx*QPSK_TxSig;
Int1=a_i1*QPSK_Int1;
Int2=a_i2*QPSK_Int2;
N0=1/10^(SNR_dB/10);
Noise=sqrt(N0/2)*(randn(N,length(QPSK_TxSig))+1i*randn(N,length(QPSK_TxSig)));
Rx_Sig=Tx_Sig+Int1+Int2+Noise;
end

%% RLS Beamformer Helper
function [y,e,w_time,w_final]=rls_beamformer(Rx_Sig,desired,N,lambda,delta)
w=zeros(N,1);
P=(1/delta)*eye(N);
L=size(Rx_Sig,2);
y=zeros(1,L);
e=zeros(1,L);
w_time=zeros(N,L);
for n=1:L
x=Rx_Sig(:,n);
y(n)=w.' * x;
e(n)=desired(n)-y(n);
alpha=desired(n)-x.'*w;
g1=P*conj(x);
g2=lambda+x.'*P*conj(x);
g=g1/g2;
P=(1/lambda)*P-g*x.'*(1/lambda)*P;
w=w+alpha*g;
w_time(:,n)=w;
end
w_final=w_time(:,end);
end

%% LMS Beamformer Helper
function [y,e,w_time,w_final]=lms_beamformer(Rx_Sig,desired,N,mu)
w=zeros(N,1);
L=size(Rx_Sig,2);
y=zeros(1,L);
e=zeros(1,L);
w_time=zeros(N,L);
for n=1:L
x=Rx_Sig(:,n);
y(n)=w.' * x;
e(n)=desired(n)-y(n);
w=w+mu*conj(x)*e(n);
w_time(:,n)=w;
end
w_final=w_time(:,end);
end

%% NLMS Beamformer Helper
function [y,e,w_time,w_final]=nlms_beamformer(Rx_Sig,desired,N,mu)
w=zeros(N,1);
L=size(Rx_Sig,2);
y=zeros(1,L);
e=zeros(1,L);
w_time=zeros(N,L);
for n=1:L
x=Rx_Sig(:,n);
y(n)=w.' * x;
e(n)=desired(n)-y(n);
w=w+(mu/(x'*x + 1e-6)) * conj(x) * e(n);
w_time(:,n)=w;
end
w_final=w_time(:,end);
end

%% Moving Average Helper
function y=moving_average(x,span)
x=x(:).';
n=length(x);
y=zeros(1,n);
half=floor(span/2);
cs=[0,cumsum(x)];
for i=1:n
lo=max(1,i-half);
hi=min(n,i+half);
y(i)=(cs(hi+1)-cs(lo))/(hi-lo+1);
end
end

%% Empirical CDF Helper
function [xs,ps]=empirical_cdf(x)
xs=sort(x(:).');
n=numel(xs);
if n==0
ps=[];
else
ps=(1:n)/n;
end
end

%% Percentile Helper
function p=prctile_compat(x,pct)
x=sort(x(:));
if isempty(x)
p=NaN;
return;
end
q=max(0,min(100,pct))/100;
idx=1+(numel(x)-1)*q;
lo=floor(idx);hi=ceil(idx);
if lo==hi
p=x(lo);
else
p=x(lo)+(idx-lo)*(x(hi)-x(lo));
end
end

%% Half Power Beamwidth Helper
function bw=half_power_beamwidth(theta_deg,pattern_dB,peakIdx)
threshold=pattern_dB(peakIdx)-3;
left=peakIdx;
while left>1&&pattern_dB(left)>=threshold
left=left-1;
end
right=peakIdx;
while right<numel(pattern_dB)&&pattern_dB(right)>=threshold
right=right+1;
end
bw=theta_deg(right)-theta_deg(left);
end

%% Monte Carlo BER Helper
function [BER_bf,BER_nobf,total_bits]=ber_monte_carlo(...
N,kd,theta_tx,theta_i1,theta_i2,SIR1_dB,SIR2_dB,SNR_dB,...
lambda,delta,moduRad,bitsPerBlock,minErrors,maxBlocks)
T_err_bf=0;T_err_nobf=0;T_bits=0;blocks=0;
while (T_err_bf<minErrors)&&(blocks<maxBlocks)
blocks=blocks+1;
[Rx_Sig,QPSK_TxSig,uncoded_bits_Tx]=generate_scenario(...
N,kd,theta_tx,theta_i1,theta_i2,SIR1_dB,SIR2_dB,SNR_dB,bitsPerBlock,moduRad);
[y,~,~,~]=rls_beamformer(Rx_Sig,QPSK_TxSig,N,lambda,delta);
bits_rx_bf=qpsk_demod(y);
T_err_bf=T_err_bf+sum(abs(uncoded_bits_Tx-bits_rx_bf));
y_ref=Rx_Sig(1,:);
bits_rx_nobf=qpsk_demod(y_ref);
T_err_nobf=T_err_nobf+sum(abs(uncoded_bits_Tx-bits_rx_nobf));
T_bits=T_bits+length(uncoded_bits_Tx);
end
if T_err_bf==0
BER_bf=0.5/T_bits;
else
BER_bf=T_err_bf/T_bits;
end
if T_err_nobf==0
BER_nobf=0.5/T_bits;
else
BER_nobf=T_err_nobf/T_bits;
end
total_bits=T_bits;
end

%% Polar Pattern Plot Helper
function pax=polar_pattern_plot(theta_rad,rho,theta_tx,theta_i1,theta_i2)
hold on;axis equal;axis off;
rMax=1.05;
for rr=[0.25 0.5 0.75 1.0]
ang=linspace(0,pi,181);
plot(rr*cos(ang),rr*sin(ang),'Color',[0.85 0.85 0.85]);
end
for a_deg=0:30:180
a=deg2rad(a_deg);
plot([0 rMax*cos(a)],[0 rMax*sin(a)],'Color',[0.85 0.85 0.85]);
text(1.08*cos(a),1.08*sin(a),sprintf('%d^{\\circ}',a_deg),...
'HorizontalAlignment','center','FontSize',9);
end
hPattern=plot(rho.*cos(theta_rad),rho.*sin(theta_rad),'Color',[0 0.45 0.74],'LineWidth',1.6);
hTx=plot(1.0*cos(theta_tx),1.0*sin(theta_tx),'o','MarkerFaceColor',[0 0.6 0],'MarkerEdgeColor','k','MarkerSize',8);
hI1=plot(1.0*cos(theta_i1),1.0*sin(theta_i1),'*','Color','r','MarkerSize',10,'LineWidth',1.5);
hI2=plot(1.0*cos(theta_i2),1.0*sin(theta_i2),'*','Color','r','MarkerSize',10,'LineWidth',1.5);
legend([hPattern hTx hI1 hI2],{'Array Pattern','Desired User','Interferer 1','Interferer 2'},...
'Location','southoutside','Orientation','horizontal');
xlim([-1.15 1.15]);ylim([-0.15 1.15]);
pax=gca;
end

%% Direction Ray Drawing Helper
function draw_doa_ray(originX,theta,len,rgbColor)
xTip=originX+len*cos(theta);
yTip=len*sin(theta);
plot([originX xTip],[0 yTip],'-','Color',rgbColor,'LineWidth',1.8);
headLen=0.18*len;
headAng=deg2rad(18);
a1=theta+pi-headAng;
a2=theta+pi+headAng;
p1=[xTip+headLen*cos(a1),yTip+headLen*sin(a1)];
p2=[xTip+headLen*cos(a2),yTip+headLen*sin(a2)];
fill([xTip p1(1) p2(1)],[yTip p1(2) p2(2)],rgbColor,'EdgeColor',rgbColor);
end

%% Figure Title Compatibility Helper
function sgtitle_compat(str)
try
sgtitle(str);
catch
annotation('textbox',[0 0.93 1 0.06],'String',str,...
'HorizontalAlignment','center','EdgeColor','none',...
'FontWeight','bold','FontSize',12);
end
end

%% Automatic Figure Saving Helper
function CFG=save_thesis_figure(fig,CFG,captionStr)
CFG.figIdx=CFG.figIdx+1;
fname=sprintf('%d%d',CFG.ChapterNumber,CFG.figIdx);
outfile=fullfile(CFG.OutDir,[fname '.png']);
set(fig,'Color','w');
try
exportgraphics(fig,outfile,'Resolution',300);
catch
set(fig,'PaperPositionMode','auto');
print(fig,outfile,'-dpng','-r300');
end
logfile=fullfile(CFG.OutDir,'Figure_List_Chapter4.txt');
fid=fopen(logfile,'a');
fprintf(fid,'Figure %d.%d  (%s.png)\n    %s\n\n',CFG.ChapterNumber,CFG.figIdx,fname,captionStr);
fclose(fid);
fprintf('  -> Saved %s.png\n',fname);
end

%% Automatic Table Saving Helper
function CFG=save_thesis_table(T,headers,dataRows,CFG,tableLabel,shortName,captionStr)
fig=draw_table_figure([tableLabel ' - ' shortName],headers,dataRows);
CFG=save_thesis_figure(fig,CFG,[tableLabel ' (' shortName '): ' captionStr]);
baseName=sprintf('Table_%d_%s',CFG.ChapterNumber,matlab.lang.makeValidName(shortName));
try
writetable(T,fullfile(CFG.OutDir,[baseName '.csv']));
writetable(T,fullfile(CFG.OutDir,[baseName '.xlsx']));
catch ME
warning('Could not export %s to CSV/XLSX: %s',baseName,ME.message);
end
end

%% Table Figure Drawing Helper
function fig=draw_table_figure(titleStr,headers,dataRows)
nCols=numel(headers);
nRows=size(dataRows,1);
colW=190;
rowH=32;
figW=60+colW*nCols;
figH=110+rowH*(nRows+1);
fig=figure('Color','w','Units','pixels','Position',[100 100 figW figH]);
ax=axes('Parent',fig,'Position',[0 0 1 1]);hold(ax,'on');
axis(ax,[0 nCols 0 (nRows+2)]);
set(ax,'YDir','reverse');
axis(ax,'off');
text(nCols/2,0.55,titleStr,'FontWeight','bold','FontSize',13,...
'HorizontalAlignment','center','Parent',ax);
rectangle('Position',[0 1 nCols 1],'FaceColor',[0.20 0.29 0.42],...
'EdgeColor','none','Parent',ax);
for c=1:nCols
text(c-0.5,1.5,headers{c},'FontWeight','bold','FontSize',10.5,...
'Color','w','HorizontalAlignment','center','Parent',ax);
end
for r=1:nRows
y0=1+r;
if mod(r,2)==0
rectangle('Position',[0 y0 nCols 1],'FaceColor',[0.93 0.95 0.98],...
'EdgeColor','none','Parent',ax);
end
for c=1:nCols
text(c-0.5,y0+0.5,dataRows{r,c},'FontSize',9.5,...
'HorizontalAlignment','center','Parent',ax);
end
end
for r=0:(nRows+1)
line([0 nCols],[1+r 1+r],'Color',[0.75 0.75 0.75],'Parent',ax);
end
for c=0:nCols
line([c c],[1 nRows+2],'Color',[0.75 0.75 0.75],'Parent',ax);
end
end

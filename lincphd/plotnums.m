%--- plot estimated number of targets vs time

%estimated number of targets
subplot(2,1,1); hold;
stairs([1:K],N_true,'k');
plot([1:K],hat_N,'k.');

set(gca, 'XLim',[1 K]);
set(gca, 'YLim',[0 max(N_true)+1]);
legend(gca,'True','Estimated');
xlabel('Time','fontsize',10);
ylabel('Number of Targets','fontsize',10); grid on;

% %cardinality dn plot: mean +/- 1 std dev
% subplot(2,1,2); hold;
% stairs([1:K],N_true,'b');
% plot([1:K],cdn_update_mean,'k');
% plot([1:K],cdn_update_mean+sqrt(cdn_update_var),'k:',...
%      [1:K],cdn_update_mean-sqrt(cdn_update_var),'k:');
% 
% set(gca, 'XLim',[1 K]);
% set(gca, 'YLim',[0 max(N_true)+1]);
% legend(gca,'True','Mean','Std Dev');
% xlabel('time step','fontsize',10); 
% ylabel('cardinality distribution','fontsize',10);

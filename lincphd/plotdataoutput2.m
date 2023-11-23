%--- coordinate plot with respect to time
% you need to define model.C_posn in order to work properly

[X_track,k_birth,k_death]= extract_tracks(X,track_list,total_tracks);

subplot(211); box on;

for k=1:K
    hlined= line(k*ones(size(Z{k},2),1),Z{k}(1,:)','LineStyle','none','Marker','x',...
        'Markersize',5,'Color',0.7*ones(1,3));
    if ~isempty(hat_X{k}),
        P= model.C_posn*hat_X{k};

    
        hline2= line(k*ones(size(hat_X{k},2),1),P(1,:)','LineStyle','none','Marker','.',...
            'Markersize',8,'Color',0*ones(1,3));
    end;
end;

for i=1:total_tracks
    P= model.C_posn*X_track(:,k_birth(i):1:k_death(i),i);
    hline1= line(k_birth(i):1:k_death(i),P(1,:),'LineStyle','-','Marker','none',...
        'LineWidth',1,'Color',0*ones(1,3));
end;

subplot(212); box on;
    
for k=1:K
    yhlined= line(k*ones(size(Z{k},2),1),Z{k}(2,:)','LineStyle','none','Marker','x',...
        'Markersize',5,'Color',0.7*ones(1,3));
    if ~isempty(hat_X{k}),
        P= model.C_posn*hat_X{k};

        yhline2= line(k*ones(size(hat_X{k},2),1),P(2,:)','LineStyle','none','Marker','.',...
            'Markersize',8,'Color',0*ones(1,3));
    end;
end;

for i=1:total_tracks
    P= model.C_posn*X_track(:,k_birth(i):1:k_death(i),i);
    yhline1= line(k_birth(i):1:k_death(i),P(2,:),'LineStyle','-','Marker','none',...
        'LineWidth',1,'Color',0*ones(1,3));
end;

subplot(211); xlabel('Time'); ylabel('x-coordinate (m)');
set(gca, 'XLim',[1 K]); set(gca, 'YLim',range_c(1,:));


subplot(212); xlabel('Time'); ylabel('y-coordinate (m)');
set(gca, 'XLim',[1 K]); set(gca, 'YLim',range_c(2,:));
legend([yhline2 yhline1 yhlined],'Estimates          ','True tracks','Measurements');

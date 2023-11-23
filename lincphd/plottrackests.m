% Position plot for the signal

%--- range of the plots
limit= 1000*[ -1 1; -1 1 ]; limit2= 1000*[ -1 1 -1 1 ];

%---
[X_track,k_birth,k_death]= extract_tracks(X,track_list,total_tracks);

%plot tracks
figure(1); clf;
for i=1:total_tracks
    P= model.C_posn*X_track(:,k_birth(i):1:k_death(i),i);
    plot(P(1,:),P(2,:),'k-'); hold on;
end;

for k=1:K
    plot(hat_X{k}(1,:),hat_X{k}(3,:),'LineStyle','none','Marker','o','Markersize',4,'Color',0*ones(1,3));
end
xlabel('x coordinate (m)','fontsize',15); ylabel('y coordinate (m)','fontsize',15);
axis(limit2); axis square;


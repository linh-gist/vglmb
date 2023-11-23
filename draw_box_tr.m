% Draw Box
function draw_box(t,frame, state, id)
% state : sn By N
% sn: state dimension
% N: the number of sample
% state(:,1) = [center X, center Y, width, height]'

% Draw Frame
imshow(uint8(frame),'Border','tight');
hold on;
text(5, 65, 'Tracking', 'Color','r', 'FontWeight','bold', 'FontSize',25);
text(5, 18, num2str(t), 'Color','r', 'FontWeight','bold', 'FontSize',25);


[sn N] = size(state);

for i=1:N
    W = sqrt(state(3,i)*state(4,i)); 
    H = state(3,i)./W;
%     X = state(1,i);
%     Y = state(2,i);
    X = state(1,i) - W/2;
    Y = state(2,i) - H/2;

    plot(X,Y,'o','linewidth',2); % Center
%     X = state(1,i); Y= state(2,i);
%     W = state(3,i); H = state(4,i);
    R1 = [X + floor(W/2);Y + floor(H/2)];
    R2 = [X - floor(W/2);Y + floor(H/2)];
    R3 = [X - floor(W/2);Y - floor(H/2)];
    R4 = [X + floor(W/2);Y - floor(H/2)];
    R5 = [X + floor(W/2);Y + floor(H/2)];
%     R1 = [X ; Y];
%     R2 = [X + floor(W/2);Y];
%     R3 = [X ;Y + floor(H/2)];
%     R4 = [X - floor(W/2);Y + floor(H/2)];
%     R5 = [X - floor(W/2);Y - floor(H/2)];
    BOX = [R1 R2 R3 R4 R5];
    
    idx = getColorFromID(id(i));
    line(BOX(1,:),BOX(2,:),'linewidth',2, 'Color',idx);
%     line(BOX(1,:),BOX(2,:),'linewidth',3);
end
% drawnow
% hold off;
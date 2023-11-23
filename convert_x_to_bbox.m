
%convert_x_to_bbox

function bbox = convert_x_to_bbox(x)
w = sqrt(x(3,:)*x(4,:));
h = x(3,:)./w;
bbox = [x(1,:)-w./2,x(2,:)-h./2,x(1,:)+w./2,x(2,:)+h./2]';
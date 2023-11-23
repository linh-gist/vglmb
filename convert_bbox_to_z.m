
%convert_bbox_to_z

function z = convert_bbox_to_z(bbox)

w = bbox(3,:)-bbox(1,:);
h = bbox(4,:)-bbox(2,:);
x = bbox(1,:)+w./2;
y = bbox(2,:)+h./2;
s = w.*h;
r = w./h;
sc = bbox(5,:);
z = [x;y;s;r;sc];
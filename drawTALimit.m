function drawTALimit(trackingArea,camPar)
% 
% (C) Anton Andriyenko, 2012
%
% The code may be used free of charge for non-commercial and
% educational purposes, the only requirement is that this text is
% preserved within the derivative work. For any other purpose you
% must contact the authors for permission. This code may not be
% redistributed without written permission from the authors.

% global sceneInfo

c1=trackingArea([1 3]);
c2=trackingArea([2 3]);
c3=trackingArea([2 4]);
c4=trackingArea([1 4]);


% camPar=camPar;
[mR mT]=getRoTran(camPar);

[slx sly]=worldToImg(c1(1),c1(2),0,mR,mT,camPar.mInt,camPar.mGeo);
x(1)=slx; y(1)=sly;

[slx sly]=worldToImg(c2(1),c2(2),0,mR,mT,camPar.mInt,camPar.mGeo);
x(2)=slx; y(2)=sly;

[slx sly]=worldToImg(c3(1),c3(2),0,mR,mT,camPar.mInt,camPar.mGeo);
x(3)=slx; y(3)=sly;

[slx sly]=worldToImg(c4(1),c4(2),0,mR,mT,camPar.mInt,camPar.mGeo);
x(4)=slx; y(4)=sly;

line([x x(1)],[y y(1)],'linewidth',2,'color','w','linestyle','--');

end
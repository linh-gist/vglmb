function detections3D=projToGP(detections,camPar)
    detections3D = zeros(2,size(detections,2));
    F=length(detections);
%     if opt.track3d
        heightweight=zeros(F,0);
        % height prior:
        muh=1.7; sigmah=.7; factorh=1/sigmah/sqrt(2*pi);

        [mR mT]=getRoTran(camPar);

%         for t=1:length(detections)
            ndet=size(detections,2);
%             detections(t).xw=zeros(1,ndet);
%             detections(t).yw=zeros(1,ndet);

            for det=1:ndet
                [xw yw zw]=imageToWorld(detections(1,det), detections(2,det), camPar);
                detections3D(1,det)=xw;
                detections3D(2,det)=yw;

                % one meter
                xi=detections(1,det)+detections(3,det)/2; yi=detections(2,det)+detections(4,det)/2;
                [xiu yiu]=worldToImage(xw,yw,1000,mR,mT,camPar.mInt,camPar.mGeo);
                onemeteronimage=norm([xi yi]-[xiu yiu]);
                worldheight=detections(4,det)/onemeteronimage; % in meters
                weight=normpdf(worldheight,muh,sigmah)/factorh;

%                 detections3D(3,det)=detections(5,det)*weight;
            end

            %         detections(t).xp=detections(t).xw;
            %         detections(t).yp=detections(t).yw;
%         end
%     end
end
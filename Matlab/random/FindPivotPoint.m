function [pivotpoint,minv] = FindPivotPoint(HelicalAxes_point,HelicalAxes_vect,r0)

% r0: 1x3 vector die startpositie aangeeft
% HelicalAxes_point is Mx3 matrix met punten op de rechten
% HelicalAxes_vect is Mx3 matrix met vectoren volgens rechten

% normaliseer helicalaxes_vect
for i = 1:size(HelicalAxes_vect,1)
    HelicalAxes_vect(i,:) = HelicalAxes_vect(i,:)./norm(HelicalAxes_vect(i,:));
end
minv = DistanceFromVertexToLine2(r0,HelicalAxes_point,HelicalAxes_vect);

% Checken of er naburige waarden zijn die een kleinere ssd hebben, indien
% ja: naburige waarde wordt nieuwe startwaarde
for step = [50 10 5 1 .1 .01 .001 .0001]
    next = 0;
    while(~next)
        minxneg = DistanceFromVertexToLine2(r0+[-step 0 0],HelicalAxes_point,HelicalAxes_vect);
        minxpos = DistanceFromVertexToLine2(r0+[step 0 0],HelicalAxes_point,HelicalAxes_vect);
        minyneg = DistanceFromVertexToLine2(r0+[0 -step 0],HelicalAxes_point,HelicalAxes_vect);
        minypos = DistanceFromVertexToLine2(r0+[0 step 0],HelicalAxes_point,HelicalAxes_vect);
        minzneg = DistanceFromVertexToLine2(r0+[0 0 -step],HelicalAxes_point,HelicalAxes_vect);
        minzpos = DistanceFromVertexToLine2(r0+[0 0 step],HelicalAxes_point,HelicalAxes_vect);
        minvect = [minxneg; minxpos; minyneg; minypos; minzneg; minzpos];
        
        [min_tmp, indx] = min(minvect);
        if(min_tmp < minv)
            minv = min_tmp;
            switch indx
                case 1, r0 = r0 + [-step 0 0];
                case 2, r0 = r0 + [step 0 0];
                case 3, r0 = r0 + [0 -step 0];
                case 4, r0 = r0 + [0 step 0];
                case 5, r0 = r0 + [0 0 -step];
                case 6, r0 = r0 + [0 0 step];
            end
        else
            next = 1;
        end
    end
    
end
pivotpoint = r0;
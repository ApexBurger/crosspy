


nodeCoords = [  193, 94; 
                193, 190;
                193, 286;
                289, 94;
                289, 190;
                289, 286;
                385, 94; 
                385, 190;
                385, 286; ];
            Position_X = reshape(nodeCoords(:,2),[3 3])'
            Position_Y = reshape(nodeCoords(:,1),[3 3])'

displacementsP = [  3.6, -2;
                   3.3, -1.9;
                   3.1, -1.8;
                   3.6, -1.8;
                   3.4, -1.7;
                   3.1, -1.6;
                   3.8, -1.7;
                   3.5, -1.5;
                   3.2, -1.4; ];
 Shift_X  = reshape(displacementsP(:,2),[3 3])';
 Shift_Y  = reshape(displacementsP(:,1),[3 3])';
%                
 displacementsP = [ 0 0; 0 0; 0 2; 0 2; 0 0; 0 0;0 2; 0 0; 0 0];
  Shift_X = [0 0 0; 0 0 0; 0 0 0] 
  Shift_Y = [0 0 0; 0 0 0; 2 2 2]
 
 

[strain_11, strain_22, strain_12] = fDIC_StrainCalc(Shift_X, Shift_Y, Position_X, Position_Y)
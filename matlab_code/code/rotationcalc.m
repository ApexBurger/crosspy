function [ params1 ] = rotationcalc( x_shift,y_shift,x_pos,y_pos,params0,UB,LB )
options = optimset('lsqnonlin');
options.Display = 'on';
options.UseParallel='never';

params1 = lsqnonlin(@(params0) rotationobj(params0,x_shift,y_shift,x_pos,y_pos),params0,LB,UB,options);

end


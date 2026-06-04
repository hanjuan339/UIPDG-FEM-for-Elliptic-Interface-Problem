function nx=upfemitfnx(lvsfunx,lvsfuny,x,y)
%Find the x component of the unit outward normal to the interface in the unit square.
%   lvsfun: its zero level set is the interface
%   lvsfunx: the partial derivative of the lvsfun w.r.t. x
%   lvsfuny: the partial derivative of the reflvsfun w.r.t. y
%   x,y: may be matrices of the same size
%   nx: the x-component of the unit outward normal to the interface.
%   
zx=lvsfunx(x,y);
zy=lvsfuny(x,y);
nx=zx./sqrt(zx.^2+zy.^2);

function nx=upfemrefitfnx(reflvsfunx,reflvsfuny,x,y)
%Find the x component of the unit outward normal to the interface in the unit square.
%   reflvsfun: its zero level set is the interface
%   reflvsfunx: the partial derivative of the reflvsfun w.r.t. x
%   reflvsfuny: the partial derivative of the reflvsfun w.r.t. y
%   x,y: may be matrices of the same size
%   nx: the x-component of the unit outward normal to the interface.
%   
zx=reflvsfunx(x,y);
zy=reflvsfuny(x,y);
nx=zx./sqrt(zx.^2+zy.^2);

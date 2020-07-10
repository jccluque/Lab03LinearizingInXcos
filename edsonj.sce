function xdot=edsonj(u1,u2,u3,u4,u5,u6,u7,u8)
// This is nonlinear EDSON-J model

// Load the parameters
exec('edsonjParameters.sce', -1);

// state variables
u=u1;		
v=u2;
r=u3;
x=u4;
y=u5;
apsi=u6;

// control variables
nc=u7;	// modo común (n1+n2)/2
nd=u8;	// modo diferencial (n1-n2)/2	

//

vare=[u v r]';


M=[m-Xudot     0          0;
    0       m-Yvdot   m*xg-Yrdot;
    0     m*xg-Yrdot   Iz-Nrdot];
//
C=[     0                  -m*r               -m*xg*r+Yvdot*v+((Yrdot+Nvdot)/2)*r ;
       m*r                 0                -Xudot*u;
  m*xg*r-Yvdot*v-((Yrdot-Nvdot)/2)*r   Xudot*u       0];
//
D=[ -Xu-Xuu*abs(u)      0                 0;
    0            -Yv-Yvv*abs(v)           0;
    0                  0          -Nr-Nrr*abs(r)];

// propeller force starboard and portboard, respectively 

tau_c=[Xunc*u*nc+Xncnc*nc*nc+Xrnd*r*nd+Xndnd*nd;
       Yvnc*v*nc+Yrnc*r*nc;
       Nvnc*v*nc+Nrnc*r*nc+Nund*u*nd+Nncnd*nc*nd];
    //tau_c=[Fs+Fp;
    //        %0;
    //   %(Fs-Fp)*0.9];
// The dynamics are contained in the equation below
// forces = 
forces=-C*vare-D*vare+tau_c;
varedot=inv(M)*forces;
// review with FOSSEN
etadot=[u-sin(apsi)*v;
        sin(apsi)*u+cos(apsi)*v;
         r];
// For little apsi angles sin(apsi)~= apsi and cos(apsi)~= 1

xdot=[varedot;
      etadot];
endfunction

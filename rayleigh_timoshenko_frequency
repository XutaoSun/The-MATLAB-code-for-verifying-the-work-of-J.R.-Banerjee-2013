clear all;
clc;
format long;
E=216e9; %Young's modulus
G=7.894736e10; %shear modulus
k=5/6; %section shape factor (rectangular cross-section)
I=0.025*0.0078^3/12; %Second moment of area
Ip=1.11449E-8; %polar second moment of area
A=0.025*0.0078; %Cross-sectional area
Rho=7850; %density
b=0.025; %Cross-sectional width
h=0.0078; %Cross-sectional height
nodes=4;
sup1=[1 2 3]; %left boundary condition
sup2=[]; %right boundary condition
L=[0.08 0.12];
Xi=0.2; %Crack depth to sectional height ratio
Mu=0.28; %Poisson's ratio
F(1,1)=exp(1)^(1/(1-Xi))*(-0.326584e-5*Xi+1.455190*Xi^2-0.984690*Xi^3+4.895396*Xi^4-6.501832*Xi^5+12.792091*Xi^6-26.723556*Xi^7+35.073593*Xi^8-34.954632*Xi^9+9.054062*Xi^10);
F(2,2)=exp(1)^(1/(1-Xi))*(-0.326018e-6*Xi+1.454954*Xi^2-1.455784*Xi^3-0.421981*Xi^4-0.279522*Xi^5+0.455399*Xi^6-2.432830*Xi^7+5.427219*Xi^8-6.643057*Xi^9+4.466758*Xi^10);
F(3,3)=exp(1)^(1/(1-Xi))*(-0.219628e-4*Xi+52.37903*Xi^2-130.2483*Xi^3+308.442769*Xi^4-602.445544*Xi^5+939.044538*Xi^6-1310.95029*Xi^7+1406.52368*Xi^8-1067.4998*Xi^9+391.536356*Xi^10);
c11=(1-Mu^2)/(E*b)*F(1,1);
c22=(1-Mu^2)/(E*b)*F(2,2);
c33=(1-Mu^2)/(E*b*h^2)*F(3,3);
reqmodes=[1 2 3 4 5 6]; %Required natural frequencies, must be arranged in ascending order of modes
fr=length(reqmodes); %Total number of required natural frequencies
w=0.000001; %Trial value for 'w'
J0=WWalgorithm(w,nodes,sup1,sup2,E,G,k,I,Ip,A,Rho,L,Mu,c11,c22,c33);
for r=reqmodes %Calculating the required natural frequencies, the orders of which are defined by 'reqmodes'
    while r>J0 %Establishing the lower and upper limits around the required natural frequency
	    wl=w; 
		w=w*2; 
		wu=w;
		J0=WWalgorithm(w,nodes,sup1,sup2,E,G,k,I,Ip,A,Rho,L,Mu,c11,c22,c33);
	end
	while (wu-wl)>0.0000001 %Required accuracy
	    w=(wu+wl)/2; %Bisection
		J0=WWalgorithm(w,nodes,sup1,sup2,E,G,k,I,Ip,A,Rho,L,Mu,c11,c22,c33);
		if r>J0
		    wl=w;
		else
		    wu=w;
		end
	end
	wr(r)=w;
end
fHz(1:fr,1)=wr(reqmodes)./(2*pi); %Natural frequencies in Hz
ApproxiHz(1:fr,1)=fix(fHz*1000)/1000; %Up to three decimal places
function[WW]=WWalgorithm(w,nodes,sup1,sup2,E,G,k,I,Ip,A,Rho,L,Mu,c11,c22,c33)
%Code to be stored in a file named 'WWalgorithm.m'.
%'w' is the trial frequency in rad/s.
%'L' is a matrix of the length of the intact beam segment
Jm=0; ng=0; %Initial values for the terms of the Wittrick-Williams algorithm, where 'Jm' is the number of fixed end natural frequencies below the trial frequency 'w', 'ng' is the number of negative leading diagonal elements of the upper triangular matrix formed from the global stiffness matrix after performing Gaussian elimination without row interchange.
DFGstiff=zeros(3*nodes,3*nodes); %Global dynamic stiffness matrix of the cracked beam.
crackstiff=[1/c11 0 0 -1/c11 0 0; 0 1/c22 0 0 -1/c22 0; 0 0 1/c33 0 0 -1/c33; -1/c11 0 0 1/c11 0 0; 0 -1/c22 0 0 1/c22 0; 0 0 -1/c33 0 0 1/c33]; %Dynamic stiffness matrix of the crack element
for i=1:2 %Two intact beam segments
	Alpha2(i)=Rho*A*w^2*L(i)^2/(E*A);
	Beta2(i)=Rho*Ip*Mu^2*w^2/(E*A);
	Gamma(i)=(Alpha2(i)/(1-Beta2(i)))^0.5;
	bb(i)=Rho*A*w^2*L(i)^4/(E*I);
	rr(i)=I/(A*L(i)^2);
	ss(i)=E*I/(k*A*G*L(i)^2);
	Phi(i)=(bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5)^0.5;
	Lambda(i)=(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5)^0.5;
	Z(i)=Phi(i)-bb(i)*ss(i)/Phi(i);
	Eta(i)=Z(i)/(Lambda(i)+bb(i)*ss(i)/Lambda(i));
	Tau(i)=(Lambda(i)+Eta(i)*Phi(i))/(2*Eta(i)*(1-cos(Phi(i))*cosh(Lambda(i)))+(1-Eta(i)^2)*sin(Phi(i))*sinh(Lambda(i)));
	a1(i)=E*A/L(i)*Gamma(i)*(1-Beta2(i))*cot(Gamma(i));
	a2(i)=-E*A/L(i)*Gamma(i)*(1-Beta2(i))*csc(Gamma(i));
	d1(i)=E*I/L(i)^3*bb(i)*Tau(i)*(cos(Phi(i))*sinh(Lambda(i))+Eta(i)*sin(Phi(i))*cosh(Lambda(i)))/(Lambda(i)*Phi(i));
	d2(i)=E*I/L(i)^2*Z(i)*Tau(i)*((Phi(i)+Eta(i)*Lambda(i))*sin(Phi(i))*sinh(Lambda(i))-(Lambda(i)-Eta(i)*Phi(i))*(1-cos(Phi(i))*cosh(Lambda(i))))/(Lambda(i)+Eta(i)*Phi(i));
	d3(i)=E*I/L(i)*Tau(i)*(sin(Phi(i))*cosh(Lambda(i))-Eta(i)*cos(Phi(i))*sinh(Lambda(i)));
	d4(i)=-E*I/L(i)^3*bb(i)*Tau(i)*(sinh(Lambda(i))+Eta(i)*sin(Phi(i)))/(Lambda(i)*Phi(i));
	d5(i)=E*I/L(i)^2*Z(i)*Tau(i)*(cosh(Lambda(i))-cos(Phi(i)));
	d6(i)=E*I/L(i)*Tau(i)*(Eta(i)*sinh(Lambda(i))-sin(Phi(i)));
	Jm=Jm+int16(fix(Gamma(i)/pi));
	DFGstiff(((6*i)-5):(6*i),((6*i)-5):(6*i))=[a1(i) 0 0 a2(i) 0 0; 0 d1(i) d2(i) 0 d4(i) d5(i); 0 d2(i) d3(i) 0 -d5(i) d6(i); a2(i) 0 0 a1(i) 0 0; 0 d4(i) -d5(i) 0 d1(i) -d2(i); 0 d5(i) d6(i) 0 -d2(i) d3(i)];
end
DFGstiff(4:9,4:9)=DFGstiff(4:9,4:9)+crackstiff;
%Assembling the global dynamic stiffness matrix of the cracked beam.
for i=1:2
	Jm=Jm+int16(fix(Phi(i)/pi))-0.5*(2-sign(d3(i))-sign(d3(i)-d6(i)^2/d3(i)));
end
if isempty(sup1)&&isempty(sup2) %In the case of no supports or restraints.
    DFGstiff=genre(DFGstiff); %Perfom Gaussian elimination using the function file 'genre.m'.
    diagonal=diag(DFGstiff);
    for i=1:length(diagonal);
        if diagonal(i)<0
	        ng=ng+1; %Sign counts of the diagonal terms.
	    end
    end
    WW=ng+Jm; %Summation of terms of the Wittrick-Williams algorithm.
else
    DFGstiff([sup2],:)=[];
    DFGstiff(:,[sup2])=[];
	DFGstiff([sup1],:)=[];
    DFGstiff(:,[sup1])=[];%Delete the rows and columns in the stiffness matrix, corresponding to the suppressed degrees of freedom.
    DFGstiff=genre(DFGstiff);
    diagonal=diag(DFGstiff);
    for i=1:length(diagonal);
        if diagonal(i)<0
	        ng=ng+1;
	    end
    end
    WW=ng+Jm;
end
end
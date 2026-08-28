% six bar linkage
% Static Equilibrium

clc;
clear;

% define the joints
A = [7 4 0];
B = [5 16 0];
C = [25 25 0];
D = [23 10 0];
E = [18 35 0];
F = [43 32 0];
G = [45 17 0];

% Define the length of links
    AB = norm(B - A);
    BC = norm(C - B);
    CD = norm(D - C);
    BE = norm(E - B);
    EF = norm(F - E);
    FG = norm(G - F);

%Define Weight of each link

WAB = [0 -1 0];
WBEC = [ 0 -1 0];
WCD = [0 -1 0];
WEF = [ 0 -1 0];
WFG = [ 0 -1 0];


%center of mass of each link
S1 = (A+B)/2;
S2 = (B+C+E)/3;
S3 = (C+D)/2; 
S4 = (E+F)/2;
S5 = (F+G)/2;




syms FAx FAy FBy FBx FCy FCx FDy FDx FFx FFy FGx FGy FEx FEy Tin

ForceA = [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
ForceE = [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG = [FGx FGy 0];
InputTorque = [ 0 0 Tin];

AppliedFroce = [50 0 0];

%Sum of Forces = 0
%Fa+Fb + WeightofAB = 0

eqn1 = ForceA + ForceB 

%Moment Equation


eqn2 = cross(A-S1,ForceA)+cross(B-S1,ForceB) + InputTorque == 0;


%Equations for Link BEC
%Sum of Forces = 0
% -Fb + Fc + Fe +WBEC = 0

eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;

%Sum of Moments = 0
eqn4 = cross(B-S2, -ForceB) + cross(C-S2, ForceC) + cross(E-S2, ForceE) == 0;



%Equations for Link CD
%Sum of Forces = 0 for Link CD
% -Fc + Fd + WCD = 0
eqn5 = -ForceC + ForceD + WCD == 0;

%Sum of Moments = 0 for Link CD
eqn6 = cross(C-S3, -ForceC) + cross(D-S3, ForceD) == 0;


%Equations for Link EF
%Sum of Forces = 0 for Link EF
% -Fe + Ff + WEF = 0
eqn7 = -ForceE + ForceF + WEF == 0;

%Sum of Moments = 0 for Link EF
eqn8 = cross(E-S4, -ForceE) + cross(F-S4, ForceF) == 0;

%Sum of Forces = 0 for Link FG
% -Ff + Fg + WFG = 0
eqn9 = -ForceF + ForceG + WFG + AppliedFroce == 0;

%Sum of Moments = 0 
% for Link FG
eqn10 = cross(F-S5, -ForceF) + cross(G-S5, ForceG) == 0;


%Solve the 10 equations

eqns = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10];
StaticSolution = solve(eqns, [FAx, FAy, FBx, FBy, FCx, FCy, FDx, FDy, FEx, FEy, FFx, FFy, FGx, FGy, Tin]);

Force_Ax = double(StaticSolution.FAx);
Force_Ay = double(StaticSolution.FAy);
Force_Bx = double(StaticSolution.FBx);
Force_By = double(StaticSolution.FBy);
Force_Cy = double(StaticSolution.FCy);
Force_Cx = double(StaticSolution.FCx);
Force_Dy = double(StaticSolution.FDy);
Force_Dx = double(StaticSolution.FDx);
Force_Ex = double(StaticSolution.FEx);
Force_Ey = double(StaticSolution.FEy);
Force_Gx = double(StaticSolution.FGx);
Force_Gy = double(StaticSolution.FGy);
Input_Torque = double(StaticSolution.Tin);

disp('ForceA');
disp([Force_Ax, Force_Ay]);
disp('ForceB');
disp([Force_Bx, Force_By]);
% Display remaining forces and input torque
disp('ForceC');
disp([Force_Cx, Force_Cy]);
disp('ForceD');
disp([Force_Dx, Force_Dy]);
disp('ForceE');
disp([Force_Ex, Force_Ey]);
disp('ForceF');
disp([FFx, FFy]);
disp('ForceG');
disp([Force_Gx, Force_Gy]);
disp('Input Torque');
disp(Input_Torque);

%Angular velocity Calculations
%Loop ABCDA

syms wBEC wCD 
omega_AB = [0 0 1];
omega_BEC = [0 0 wBEC];
omega_CD = [0 0 wCD];

eqn11 = cross(omega_AB,B-A) + cross(omega_BEC,C-B) + cross(omega_CD,D-C) ==0;

loop1Solution = solve(eqn11,[wBEC wCD]);

% Calculate angular velocities
wBEC_value = double(loop1Solution.wBEC)
wCD_value = double(loop1Solution.wCD)

omegaBEC = [0 0 wBEC_value]
omegaCD = [0 0 wCD_value]

%Loop DCEFGD
syms  wEF wFG

omega_EF = [0 0 wEF];
omega_FG = [0 0 wFG];

% Loop DCEFGD equations
eqn12 = cross(omegaCD,C-D) + cross(omegaBEC,E-C) + cross(omega_EF,F-E) + cross(omega_FG,G-F) == 0;

loop2solution = solve(eqn12,[wEF wFG]);

% Calculate angular velocities
wEF_value = double(loop2solution.wEF)
wFG_value = double(loop2solution.wFG)



% Angular Acceleration
% loop AbCDA
% Angular acceleration calculations
syms aBEC aCD 
alpha_AB= [0 0 0];
alpha_BEC=[0 0 aBEC];
alpha_CD=[0 0 aCD];

a_B_A = cross(alpha_AB,B-A)+cross(omega_AB,cross(omega_AB,B-A));
a_C_B = cross(alpha_BEC,C-B)+cross(omegaBEC,cross(omegaBEC,C-B));
a_D_C = cross(alpha_CD,D-C)+cross(omegaCD,cross(omegaCD,D-C));

eqn13 = a_B_A + a_C_B + a_D_C == 0;

loop1AccSolution = solve(eqn13,[aBEC aCD]);
alphaBEC = double(loop1AccSolution.aBEC)
alphaCD = double(loop1AccSolution.aCD)

alphaBEC_vector = [ 0 0 alphaBEC];
alphaCD_vector = [0 0 alphaCD];



%Loop 2 angular acceleration DCEFGD

syms aEF aFG 

alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];


a_C_D = cross(alphaCD_vector,C-D) + cross(omegaCD,cross(omegaCD,C-D));
a_E_C = cross(alphaBEC_vector,E-C) + cross(omegaBEC,cross(omegaBEC,E-C));

angVel_EF = [0 0 wEF_value];
angVel_FG = [0 0 wFG_value];

a_F_E = cross(alpha_EF,F-E)+cross(angVel_EF,cross(angVel_EF,F-E));

a_G_F = cross(alpha_FG,G-F)+cross(angVel_FG,cross(angVel_FG,G-F));


eqn14 = a_C_D + a_E_C + a_F_E + a_G_F == 0;

loop2AccSolution = solve(eqn14,[aEF aFG]);


% Extract angular accelerations from the solution
aEF_value = double(loop2AccSolution.aEF)
aFG_value = double(loop2AccSolution.aFG)


loop2AccSolution = solve(eqn14, [aEF aFG]);


%Velocity at Joints


vB_A = cross(omega_AB,B-A);



%vE_A = V_E_B + V_B_A;

v_E_B = cross(omegaBEC,E-B);

vE_A = v_E_B + vB_A;

%V_S4/G = V_S4 + V_F_G

V_S4_F = cross(angVel_EF,S4-F);

V_F_G = cross(angVel_FG,F-G);

vS4_G = V_S4_F + V_F_G

%% Velocity at Joints (remainder)

% Velocity of Joint C (through link CD, from fixed D)
vC_D = cross(omegaCD, C - D);

% Velocity of Joint F (through link FG, from fixed G)
vF_G = cross(angVel_FG, F - G);

disp('vC ='); disp(vC_D);
disp('vF ='); disp(vF_G);


%% Velocity at Center of Mass of Each Link

vS1 = cross(omega_AB, S1 - A);                       
vS2 = vB_A + cross(omegaBEC, S2 - B);                 
vS3 = cross(omegaCD, S3 - D);                         
% vS4_G already computed above (Link EF, from F)
vS5 = cross(angVel_FG, S5 - G);                       

disp('vS1 ='); disp(vS1);
disp('vS2 ='); disp(vS2);
disp('vS3 ='); disp(vS3);
disp('vS4 ='); disp(vS4_G);
disp('vS5 ='); disp(vS5);

%% Acceleration at Joints (remainder)

% Acceleration of Joint B 
aB = a_B_A;    

% Acceleration of Joint C 
aC = a_C_D;    

% Acceleration of Joint E 
aE = aB + cross(alphaBEC_vector, E - B) + cross(omegaBEC, cross(omegaBEC, E - B));

% Acceleration of Joint F 
aF_G = cross([0 0 aFG_value], F - G) + cross(angVel_FG, cross(angVel_FG, F - G));

disp('aA = aD = aG (fixed pivots) = [0 0 0]');
disp('aB ='); disp(aB);
disp('aC ='); disp(aC);
disp('aE ='); disp(aE);
disp('aF ='); disp(aF_G);


%% Acceleration at Center of Mass of Each Link
% (need aB, aC, aF as intermediate joint accelerations first)

aB = a_B_A;   
aC = a_C_D;    
aF_G = cross([0 0 aFG_value], F - G) + cross(angVel_FG, cross(angVel_FG, F - G));

aS1 = cross(alpha_AB, S1 - A) + cross(omega_AB, cross(omega_AB, S1 - A));
aS2 = aB + cross(alphaBEC_vector, S2 - B) + cross(omegaBEC, cross(omegaBEC, S2 - B));
aS3 = cross(alphaCD_vector, S3 - D) + cross(omegaCD, cross(omegaCD, S3 - D));
aS4 = aF_G + cross([0 0 aEF_value], S4 - F) + cross(angVel_EF, cross(angVel_EF, S4 - F));
aS5 = cross([0 0 aFG_value], S5 - G) + cross(angVel_FG, cross(angVel_FG, S5 - G));

disp('aS1 ='); disp(aS1);
disp('aS2 ='); disp(aS2);
disp('aS3 ='); disp(aS3);
disp('aS4 ='); disp(aS4);
disp('aS5 ='); disp(aS5);



%Newton's Second Law Implementation

MassAB = 1;
MassBEC = 1;
MassCD = 1;
MassEF = 1;
MassFG = 1;


%Mass Moment of Inertia

J_AB = 1;
J_BEC = 1;
J_CD = 1;
J_EF = 1;
J_FG = 1;

syms NFAx NFAy NBx NBy NCx NCy NDx NDy NEx NEy NFx NFy NGx NGy 

%define forces

NForceA = [NFAx NFAy 0];
NForceB = [NBx NBy 0];
NForceC = [NCx NCy 0];
NForceD = [NDx NDy 0];
NForceE = [NEx NEy 0];
NForceF = [NFx NFy 0];
NForceG = [NGx NGy 0];
NInputTorque = [0 0 Tin];

alphaEF_vector = [0 0 aEF_value];
alphaFG_vector = [0 0 aFG_value];

%equations for Link AB
eqn15 = NForceA +NForceB + WAB == MassAB * aS1

%sum of moments = 0
eqn16 = cross(A-S1,NForceA) + cross(B-S1,NForceB) + NInputTorque == J_AB * alpha_AB;

%equations for Link BEC

eqn17 = NForceC + NForceE -NForceB +WBEC == MassBEC * aS2

eqn18 = cross(C-S2,NForceC) + cross(E-S2,NForceE) + cross(B-S2,-NForceB) == J_BEC *  alphaBEC_vector;

%equations for link CD

eqn19 = -NForceC + NForceD +WCD == MassCD *aS3

eqn20 = cross(C-S3,-NForceC) + cross(D-S3,NForceD) == J_CD * alphaCD_vector;

%equations for link EF

eqn21 = -NForceE +NForceF +WEF == MassEF * aS4

% Sum of moments = 0 for link EF
eqn22 = cross(E-S4, -NForceE) + cross(F-S4, NForceF) == J_EF * alphaEF_vector;

% Equations for link FG
eqn23 = -NForceF + NForceG + WFG == MassFG * aS5;

% Sum of moments = 0 for link FG
eqn24 = cross(F-S5, -NForceF) + cross(G-S5, NForceG) == J_FG * alphaFG_vector;

%Solving Equations
eqns = [eqn15, eqn16, eqn17, eqn18, eqn19, eqn20, eqn21, eqn22, eqn23, eqn24];
NDynamicSolution = solve(eqns, [NFAx, NFAy, NBx, NBy, NCx, NCy, NDx, NDy, NEx, NEy, NFx, NFy, NGx, NGy, Tin]);

vecs = {aS1, aS2, aS3, aS4, aS5, alpha_AB, alphaBEC_vector, alphaCD_vector, alphaEF_vector, alphaFG_vector};
names = {'aS1','aS2','aS3','aS4','aS5','alpha_AB','alphaBEC_vector','alphaCD_vector','alphaEF_vector','alphaFG_vector'};
for i = 1:numel(vecs)
    v = vecs{i};
    fprintf('%s: %s | finite: %d | real: %d\n', names{i}, mat2str(v), all(isfinite(v)), isreal(v));
end



%Extracting values
% Extracting values from the dynamic solution
NForce_Ax = double(NDynamicSolution.NFAx)
NForce_Ay = double(NDynamicSolution.NFAy)
NForce_Bx = double(NDynamicSolution.NBx)
NForce_By = double(NDynamicSolution.NBy)
NForce_Cx = double(NDynamicSolution.NCx)
NForce_Cy = double(NDynamicSolution.NCy)
NForce_Dx = double(NDynamicSolution.NDx)
NForce_Dy = double(NDynamicSolution.NDy)
NForce_Ex = double(NDynamicSolution.NEx)
NForce_Ey = double(NDynamicSolution.NEy)
NForce_Fx = double(NDynamicSolution.NFx)
NForce_Fy = double(NDynamicSolution.NFy)
NForce_Gx = double(NDynamicSolution.NGx)
NForce_Gy = double(NDynamicSolution.NGy)
Input_Torque_N = double(NDynamicSolution.Tin)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---

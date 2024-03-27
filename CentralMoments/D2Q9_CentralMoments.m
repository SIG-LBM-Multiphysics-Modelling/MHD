clear all
clc
addpath('/Applications/Matlab_R2019b.app/toolbox/symbolic/symbolic')
%% Initialize some symbolic variables
syms U V R omega f0 f1 f2 f3 f4 f5 f6 f7 f8 BX BY real
%% Define lattice directions, weights and other useful quantities of the D2Q9 model
cx = [0, -1, -1, -1,  0,  1, 1, 1, 0];
cy = [0,  1,  0, -1, -1, -1, 0, 1, 1];
w = [4./9., 1./36., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36., 1./9.];
cs = 1/sqrt(3);
cs2 = cs^2;
cs3 = cs^3;
cs4 = cs^4;
cs6 = cs^6;
cs8 = cs^8;

f = [f0 f1 f2 f3 f4 f5 f6 f7 f8]'; % hydrodynamic populations
feq = sym(zeros(9,1));

Force = sym(zeros(9,1)); % generic force vector
T = sym(zeros(9,9)); %transformation matrix for central moments
M = zeros(9,9);  %transformation matrix for raw moments
Lambda = diag([1, 1, 1, 1, omega, omega, 1, 1, 1]);
Id = eye(9,9); % identity matrix
for i=1:9
    % build the complete equilibria
    first_order = 1/cs2*(U*cx(i)+V*cy(i));
    second_order = 1/(2*cs4)*((cx(i)*cx(i)-1/3)*U^2+...
                              (cy(i)*cy(i)-1/3)*V^2+...
                              2*cx(i)*cy(i)*U*V);
    third_order = 1/(2*cs6)*((cx(i)^2-1/3)*cy(i)*U*U*V+(cy(i)^2-1/3)*cx(i)*U*V*V);
    fourth_order = 1/(4*cs8)*((cx(i)^2-1/3)*(cy(i)^2-1/3)*U*U*V*V);
    termMHD = w(i)/(2*cs4)*(0.5*(BX*BX+BY*BY)*(cx(i)*cx(i)+cy(i)*cy(i))-(cx(i)*BX+cy(i)*BY)^2);
    feq(i) = w(i)*R*(1+first_order+second_order+third_order+fourth_order)+termMHD;
    
    % build the transformation matrix T 
    CX = cx(i)-U;
    CY = cy(i)-V;
    T(1,i) = 1;
    T(2,i) = CX;
    T(3,i) = CY;
    T(4,i) = CX*CX+CY*CY;
    T(5,i) = CX*CX-CY*CY;
    T(6,i) = CX*CY;
    T(7,i) = CX*CX*CY;
    T(8,i) = CX*CY*CY;
    T(9,i) = CX*CX*CY*CY;
    
    % build the tranformation matrix M
    CX = cx(i);
    CY = cy(i);
    M(1,i) = 1;
    M(2,i) = CX;
    M(3,i) = CY;
    M(4,i) = CX*CX+CY*CY;
    M(5,i) = CX*CX-CY*CY;
    M(6,i) = CX*CY;
    M(7,i) = CX*CX*CY;
    M(8,i) = CX*CY*CY;
    M(9,i) = CX*CX*CY*CY;
end
T = simplify(T);
%T = M;
N = T*M^(-1); %shift matrix
%N = Id;
%% HYDRO
syms k0_pre k1_pre k2_pre k3_pre k4_pre k5_pre k6_pre k7_pre k8_pre real
syms k0_star k1_star k2_star k3_star k4_star k5_star k6_star k7_star k8_star real
K_pre = simplify(T*f); %pre-collision central moments
K_eq = simplify(T*feq); % equilibrium central moments
K_pre(5) = k4_pre;
K_pre(6) = k5_pre;
K_star = simplify((Id-Lambda)*K_pre + Lambda*K_eq) %post-collision central moments

%post-collision populations
K_sym = [R k1_star k2_star k3_star k4_star k5_star k6_star k7_star k8_star];
for i=1:9
    if(K_star(i)~=sym(0))
        K_star(i) = K_sym(i);
    end
end
f_post_collision_onestep = collect(simplify(T \ K_star), K_star);

% two-steps approach
raw_moments = simplify(N^(-1)*K_star)
syms r0 r1 r2 r3 r4 r5 r6 r7 r8 real
r = [r0 r1 r2 r3 r4 r5 r6 r7 r8]'; %symbolic raw moments
f_post_collision_twosteps = collect(simplify(M\r),K_star)




% %% MHD
% syms BX BY g0 g1 g2 g3 g4 g5 g6 g7 g8 real
% g = [g0 g1 g2 g3 g4]'; % MHD populations
% w = [1/3, 1/6, 1/6, 1/6, 1/6];
% geqx = sym(zeros(5,1));
% geqy = sym(zeros(5,1));
% T = sym(zeros(5,5)); %transformation matrix for central moments
% M = zeros(5,5);  %transformation matrix for raw moments
% Id = eye(5,5);
% for i=1:5
%     geqx(i) = w(i)*(BX+3*cy(i)*(BX*V-U*BY));
%     geqy(i) = w(i)*(BY+3*cx(i)*(U*BY-BX*V));
%     CX = cx(i)-U;
%     CY = cy(i)-V;
%     T(1,i) = 1;
%     T(2,i) = CX;
%     T(3,i) = CY;
%     T(4,i) = CX*CX;
%     T(5,i) = CY*CY;
%    
%     % build the tranformation matrix M
%     CX = cx(i);
%     CY = cy(i);
%     M(1,i) = 1;
%     M(2,i) = CX;
%     M(3,i) = CY;
%     M(4,i) = CX*CX;
%     M(5,i) = CY*CY;
% end
% T = simplify(T);
% N = simplify(T*M^(-1));
% Lambda = diag([1, omega, omega, 1, 1]);
% K_pre = simplify(T*g); %pre-collision central moments
% K_eqx = simplify(T*geqx); % equilibrium central moments
% K_pre(2) = k1_pre;
% K_pre(3) = k2_pre;
% K_star = simplify((Id-Lambda)*K_pre + Lambda*K_eqx) %post-collision central moments
% K_sym = [BX k1_star k2_star k3_star k4_star];
% for i=1:5
%     if(K_star(i)~=sym(0))
%         K_star(i) = K_sym(i);
%     end
% end
% gx_post_collision_onestep = collect(simplify(T \ K_star), K_star);
% % two-steps approach
% raw_moments = simplify(N^(-1)*K_star)
% r = [r0 r1 r2 r3 r4]'; %symbolic raw moments
% gx_post_collision_twosteps = collect(simplify(M\r),K_star)
% 
% K_eqy = simplify(T*geqy); % equilibrium central moments
% K_pre(2) = k1_pre;
% K_pre(3) = k2_pre;
% K_star = simplify((Id-Lambda)*K_pre + Lambda*K_eqy) %post-collision central moments
% K_sym = [BY k1_star k2_star k3_star k4_star];
% for i=1:5
%     if(K_star(i)~=sym(0))
%         K_star(i) = K_sym(i);
%     end
% end
% gy_post_collision_onestep = collect(simplify(T \ K_star), K_star);
% % two-steps approach
% raw_moments = simplify(N^(-1)*K_star)
% r = [r0 r1 r2 r3 r4]'; %symbolic raw moments
% gy_post_collision_twosteps = collect(simplify(M\r),K_star)
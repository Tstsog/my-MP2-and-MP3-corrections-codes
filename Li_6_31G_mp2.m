% This Matlab code computes the ground state energy for lithium (Li) atom by solving the Pople-Nesbet equation using the
% 6-31G basis set ([3s,2p]) (unrestricted Hartree-Fock (uhf) calculation).
% The code also computes the correction to the ground state energy from the
% Second-order Møller–Plesset Perturbation Theory (MP2).
%
% The core Hamiltonian matrix (H_core), overlap matrix (S_ov) and two-electron integrals (tei) (Li_6_31G_tei.txt) are computed 
% by my own developing code. An obtained total energy is compared with
% those of the Psi4 computational chemistry package and of the https://cccbdb.nist.gov/energy2x.asp 
% 
% Ref: A. Szabo and N. S. Ostlund "Modern Quantum Chemistry" book 
% Ref: https://psicode.org/psi4manual/master/dfmp2.html
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
%
% May 28, 2024 & University of North Dakota 
%
%
function [] = Li_6_31G_mp2

[f_a, f_b, En_0, tei_a, tei_b, tei_ab] = Li_UHF_6_31G;
%
[f_a]; % f_a are orbital energies with alpha-spin
[f_b]; % f_b are orbital energies with beta-spin
[En_0]; % En_0 is the ground state energy from the UHF calculation with 3-21G basis set 
%
f_a = [
  -2.477475902375693 % hole state
  -0.195763084033764 % hole 
   0.037332387065282 % virtual (particle) state
   0.037332387065378 % virtual
   0.037332387065378 % virtual
   0.196702036056173 % virtual
   0.205836031095390 % virtual
   0.205836031095390 % virtual
   0.205836031095391  % virtual
   ];
%
hole_a = [1, 2];  % hole states with alpha-spin 
virt_a = [3, 4, 5, 6, 7, 8, 9];     % virtial states with alpha-spin
%
%
f_b = [
  -2.461082149276856 % hole state
   0.032189118907499 % virtual (particle) state
   0.071534765401699 % virtual
   0.071534765401699 % virtual
   0.071534765401988 % virtual
   0.233224746080843 % virtual
   0.233224746083242 % virtual
   0.233224746083242 % virtual
   0.237145234578376 % virtual
   ];
%
hole_b = [1];           % hole states with beta-spin
virt_b = [2, 3, 4, 5, 6, 7, 8, 9];  % virtual states with beta-spin 
%
En_0 = -7.431235806090976; % Ground state energy from the UHF calculation, in which the Pople-Nesbet equation is solved. 

%%%%%%%%%% MP2 calculation
%
%%% opposite-spin contribution & d_E2_os_ab
%
%              (i-alpha a-alpha | j-beta b-beta)*(i-alpha a-alpha | j-beta b-beta)
%d_E2_os_ab = --------------------------------------------------------------------
%               e^alpha_i + e^beta_j - e^alpha_a - e^beta_b            
%
d_E2_os_ab = 0.; 
for i = hole_a
    for j = hole_b
        for a = virt_a
            for b = virt_b
                denom_ab = f_a(i) + f_b(j) - f_a(a) - f_b(b);  
                %
                val_os_ab = tei_ab(i,a,j,b) * tei_ab(i,a,j,b);
                d_E2_os_ab = d_E2_os_ab + val_os_ab./denom_ab;
            end
        end
    end
end
d_E2_os_ab; % -2.664493087241966e-04 vs -0.0002652110436367 = Psi4


%%% same-spin contribution
% d_E2_ss_a = 0.;               % alpha-alpha spin contribution
%
%            1 ( (i-alpha a-alpha | j-alpha b-alpha) - (i-alpha b-alpha | j-alpha a-alpha)) * (i-alpha a-alpha | j-alpha b-alpha)
%d_E2_ss_a = - ----------------------------------------------------------------------------------------------------------------
%            2   e^alpha_i + e^alpha_j - e^alpha_a - e^alpha_b 
%
%
d_E2_ss_a = 0.;               % alpha-alpha spin contribution
for i = hole_a
    for j = hole_a
        for a = virt_a
            for b = virt_a
                denom_a = f_a(i) + f_a(j) - f_a(a) - f_a(b);  
                %
                val_ss_a = (tei_a(i,a,j,b) - tei_a(i,b,j,a))*tei_a(i,a,j,b); 
                %
                d_E2_ss_a = d_E2_ss_a + 0.5*val_ss_a./denom_a;                
            end
        end
    end
end
%
%%%
%d_E2_ss_b = 0.;                % beta-beta spin contribution
%
%            1 ( (i-beta a-beta | j-beta b-beta) - (i-beta b-beta | j-beta a-beta)) * (i-beta a-beta | j-beta b-beta)
%d_E2_ss_b = - -----------------------------------------------------------------------------------------------------
%            2   e^beta_i + e^beta_j - e^beta_a - e^beta_b 
%
%
d_E2_ss_b = 0.;                % beta-beta spin contribution
for i = hole_b
    for j = hole_b
        for a = virt_b
            for b = virt_b
                denom_b = f_b(i) + f_b(j) - f_b(a) - f_b(b);  
                %
                val_ss_b = (tei_b(i,a,j,b) - tei_b(i,b,j,a))*tei_b(i,a,j,b);
                %
                d_E2_ss_b = d_E2_ss_b + 0.5*val_ss_b./denom_b; 
            end
        end
    end
end
d_E2_ss_ab = d_E2_ss_a + d_E2_ss_b;       % -1.897684420158066e-05 vs -0.0000189743988670 = Psi4
%
E_mp2 = En_0 + d_E2_os_ab + d_E2_ss_ab     % -7.431521232243902 vs -7.4315209967904092 = Psi4 
%                                                               vs -7.431520 from https://cccbdb.nist.gov/energy2x.asp 



%%%
return
end


%%%
function [f_a, f_b, En_0, tei_a, tei_b, tei_ab] = Li_UHF_6_31G

clc; format long
itermax = 100;
alf = 1.000;
tol = 1e-8;

tei_n = 6561;
read_tei_data = fopen('Li_6_31G_tei.txt', 'r');
tei_data_n6 = textscan(read_tei_data, '%d %d %d %d %f');
%
%
p = zeros(tei_n,1); q = zeros(tei_n,1); r = zeros(tei_n,1); s = zeros(tei_n,1); vals = zeros(tei_n,1);
p(1:tei_n) = tei_data_n6{1};
q(1:tei_n) = tei_data_n6{2};
r(1:tei_n) = tei_data_n6{3};
s(1:tei_n) = tei_data_n6{4};
vals(1:tei_n) = tei_data_n6{5};
for i = 1:tei_n
    tei(p(i),q(i),r(i),s(i)) = vals(i);
%    tei(q(i),p(i),r(i),s(i)) = vals(i);    
%    tei(p(i),q(i),s(i),r(i)) = vals(i);    
%    tei(q(i),p(i),s(i),r(i)) = vals(i);   
    %
%    tei(r(i),s(i),p(i),q(i)) = vals(i);    
%    tei(s(i),r(i),p(i),q(i)) = vals(i);        
%    tei(r(i),s(i),q(i),p(i)) = vals(i);        
%    tei(s(i),r(i),q(i),p(i)) = vals(i);            
end
%
%
H0 =     [-4.44526922 -0.5288766   0.          0.          0.         -0.73461839  0.          0.          0.        ;
          -0.5288766  -0.9669833   0.          0.          0.         -0.86304442  0.          0.          0.        ;
           0.          0.         -0.78005727  0.          0.          0.         -0.53536894  0.          0.        ;
           0.          0.          0.         -0.78005727  0.          0.          0.         -0.53536894  0.        ;
           0.          0.          0.          0.         -0.78005727  0.          0.          0.         -0.53536894;
          -0.73461839 -0.86304442  0.          0.          0.         -0.85390487  0.          0.          0.        ;
           0.          0.         -0.53536894  0.          0.          0.         -0.51532696  0.          0.        ;
           0.          0.          0.         -0.53536894  0.          0.          0.         -0.51532696  0.        ;
           0.          0.          0.          0.         -0.53536894  0.          0.          0.         -0.51532696];
%
S_ov =  [1.         0.14525828 0.         0.         0.         0.18574238 0.         0.         0.        ;
         0.14525828 1.         0.         0.         0.         0.90551548 0.         0.         0.        ;
         0.         0.         1.         0.         0.         0.         0.80206344 0.         0.        ;
         0.         0.         0.         1.         0.         0.         0.         0.80206344 0.        ;
         0.         0.         0.         0.         1.         0.         0.         0.         0.80206344;
         0.18574238 0.90551548 0.         0.         0.         1.         0.         0.         0.        ;
         0.         0.         0.80206344 0.         0.         0.         1.         0.         0.        ;
         0.         0.         0.         0.80206344 0.         0.         0.         1.         0.        ;
         0.         0.         0.         0.         0.80206344 0.         0.         0.         1.        ];


%
%
Q_pqrs = tei;
H_core = H0;
%
dim = 9;
Nel = 3.;
%
N_a = 2;
N_b = Nel - N_a;
%
P_old_a = 0.5.*ones(dim,dim);
P_old_b = 0.5.*ones(dim,dim);
P_T = 0.5.*ones(dim,dim);

%%% Fock hamiltonian
%
for iter = 1:itermax
    iter;
    %
    P_a = P_old_a;
    P_b = P_old_b; 
    %
    F_a = H_core;
    F_b = H_core;    
    for p = 1:dim
        for q = 1:dim
            for r = 1:dim
                for s = 1:dim
                    F_a(p,q) = F_a(p,q) + (P_T(r,s) * Q_pqrs(p,q,r,s) -  P_a(r,s) * Q_pqrs(p,r,q,s));
                    F_b(p,q) = F_b(p,q) + (P_T(r,s) * Q_pqrs(p,q,r,s) -  P_b(r,s) * Q_pqrs(p,r,q,s));                    
                end
    
            end
    
        end
    end
    Fock_a = F_a ;
    S_mat_fock = S_ov;
    [Vec_a,En_a] = eig(Fock_a,S_mat_fock);                                     % Eigenvalue problem
    En_a = diag(En_a);
    [foo, ij] = sort(En_a);
    En_a = En_a(ij);
    [En_a(1:dim)'];

    Vec_a = Vec_a(:,ij);                       % The unnormalized eigenfunctions
%
   for i = 1:dim
        norm = 0.;
        for p = 1:dim
            for q = 1:dim
                norm = norm + Vec_a(p,i) * Vec_a(q,i) * S_ov(p,q);
            end
        end
        Vec_a(:,i) = Vec_a(:,i)/sqrt(norm);
    end
    %
    P_new_a = zeros(dim,dim);
    for i = 1:N_a
        for pp = 1:dim
            for qq = 1:dim
                P_new_a(pp,qq) = P_new_a(pp,qq) + Vec_a(pp,i)*Vec_a(qq,i);
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Fock_b = F_b ;
    S_mat_fock = S_ov;
    [Vec_b,En_b] = eig(Fock_b,S_mat_fock);                                     % Eigenvalue problem
    En_b = diag(En_b);
    [foo, ij] = sort(En_b);
    En_b = En_b(ij);
    [En_b(1:dim)'];

    Vec_b = Vec_b(:,ij);                       % The unnormalized eigenfunctions
%
   for i = 1:dim
        norm = 0.;
        for p = 1:dim
            for q = 1:dim
                norm = norm + Vec_b(p,i) * Vec_b(q,i) * S_ov(p,q);
            end
        end
        Vec_b(:,i) = Vec_b(:,i)/sqrt(norm);
    end
    %
    P_new_b = zeros(dim,dim);
    for i = 1:N_b
        for pp = 1:dim
            for qq = 1:dim
                P_new_b(pp,qq) = P_new_b(pp,qq) + Vec_b(pp,i)*Vec_b(qq,i);
            end
        end
    end
    %
    P_T = P_new_a + P_new_b;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P_n_a = alf.*P_new_a + (1.-alf).*P_a;
    P_n_b = alf.*P_new_b + (1.-alf).*P_b;

     if ( abs(sum(sum(P_n_a - P_old_a))) && abs(sum(sum(P_n_b - P_old_b))) < tol)
        break 
     end
%        
    P_old_a = P_n_a;
    P_old_b = P_n_b;    

end
%%%
En_0 = sum(0.5*diag( P_T(:,:)*(H_core(:,:)) +  P_a(:,:)*F_a(:,:) + P_b(:,:)*F_b(:,:) ));  % -7.431235806090976 vs  -7.4312368113479055 = Psi4
%
f_a = En_a;
f_b = En_b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
P_a = Vec_a'; % charge density matrix
P_b = Vec_b'; % charge density matrix
%
tei_a = zeros(dim,dim,dim,dim);  % two-electron integral (tei) in molecular orbital basis sets: (alpha alpha| alpha alpha)
tei_b = zeros(dim,dim,dim,dim);  % two-electron integral (tei) in molecular orbital basis sets: (beta beta| beta beta)
tei_ab = zeros(dim,dim,dim,dim);  % two-electron integral (tei) in molecular orbital basis sets: (alpha alpha| beta beta)
for ii = 1:dim
    for jj = 1:dim
        for kk = 1:dim
            for ll = 1:dim
                for mm = 1:dim
                    for nn = 1:dim
                        for oo = 1:dim
                            for pp = 1:dim
                                tei_a(ii,jj,kk,ll) = tei_a(ii,jj,kk,ll) + P_a(ii,mm)*P_a(jj,nn)*P_a(kk,oo)*P_a(ll,pp)*Q_pqrs(mm,nn,oo,pp);
                                tei_b(ii,jj,kk,ll) = tei_b(ii,jj,kk,ll) + P_b(ii,mm)*P_b(jj,nn)*P_b(kk,oo)*P_b(ll,pp)*Q_pqrs(mm,nn,oo,pp);     
                                tei_ab(ii,jj,kk,ll) = tei_ab(ii,jj,kk,ll) + P_a(ii,mm)*P_a(jj,nn)*P_b(kk,oo)*P_b(ll,pp)*Q_pqrs(mm,nn,oo,pp); 
                            end
                        end
                    end
                end
            end
        end
    end
end
%

%%%
return
end

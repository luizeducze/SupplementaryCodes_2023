%%% Pseudopotencial Interface Plana

clear all
close all

% Geometria do dominio
Nx = 120;
Ny = 200;

% Geometria da interface plana
radius = 25; % largura da fase liquido

% Tempo de relaxacao
tau_v = 1;
tau_rho = 1;
tau_j = 1;
tau_e = 1;
tau_s = 1;
tau_q = 1/1.99;
lambda = [ tau_rho^(-1) tau_e^(-1) tau_s^(-1) tau_j^(-1) tau_q^(-1) tau_j^(-1) tau_q^(-1) tau_v^(-1) tau_v^(-1) ];

% Matriz de relaxacao de Gram-Scmidt
M = [ 1, 1, 1, 1, 1, 1, 1, 1, 1; ...
     -4,-1,-1,-1,-1, 2, 2, 2, 2; ...
      4,-2,-2,-2,-2, 1, 1, 1, 1; ...
      0, 1, 0,-1, 0, 1,-1,-1, 1; ...
      0,-2, 0, 2, 0, 1,-1,-1, 1; ...
      0, 0, 1, 0,-1, 1, 1,-1,-1; ...
      0, 0,-2, 0, 2, 1, 1,-1,-1; ...
      0, 1,-1, 1,-1, 0, 0, 0, 0; ...
      0, 0, 0, 0, 0, 1,-1, 1,-1];
  
%Parametros Shan-Chen
G = -1;
%kappa_1 = 0.9821830475483113;
%kappa_2 = -1.6316993517870626;
%epsilon = 1.917596060230674;
kappa_1 = 0.8772583900683957;
kappa_2 = -0.9483889821543865;
epsilon = 1.947427650061337;
sigma = - epsilon*kappa_1/G/16;


% Parametros do modelo multifasico
a = 0.25;
b = 4;
R = 1;

Tc = a/b/R/2.6502;
%T = 0.80*Tc;
T = 0.60*Tc;

%rho_1 = 0.02172865606779925; % densidade da fase 1
%rho_2 = 0.30717839640781797; % densidade da fase 2
rho_1 = 0.0030824221257919567;
rho_2 = 0.40619262842616505;

% Parametros do lattice
dx = 1; % unidades de lattice
dt = 1; % unidades de lattice
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
w_SC = [4/3, 1/3, 1/3, 1/3, 1/3, 1/12, 1/12, 1/12, 1/12];
cx = [0, 1, 0,-1, 0, 1,-1,-1, 1];
cy = [0, 0, 1, 0,-1, 1, 1,-1,-1];
c = 1;
cs2 = 1/3;

% Campos de densidade e velocidade
rho = zeros(Nx,Ny);
Ux = zeros(Nx,Ny);
Vy = zeros(Nx,Ny);

% Funcao de distribuicao
f = zeros(Nx,Ny,9);
m = zeros(Nx,Ny,9);
meq = zeros(Nx,Ny,9);
mf = zeros(Nx,Ny,9);
mC = zeros(Nx,Ny,9); % Termo fonte

% Estoque Resultados
Est_1 = zeros(Nx,Ny,10);
Est_2 = zeros(Nx,Ny,41);

%*************************************************************************
% Inicializacao da Densidade
%*************************************************************************
for i = 1:Nx
    for j = 1:Ny
        
        x_pos_1 = i - 0.5 - Nx/2;
        y_pos_1 = j - 0.5 - 0.25*Ny;
        pos_1 = sqrt( x_pos_1^2 + y_pos_1^2 );
        
        x_pos_2 = i - 0.5 - Nx/2;
        y_pos_2 = j - 0.5 - 0.75*Ny;
        pos_2 = sqrt( x_pos_2^2 + y_pos_2^2 );
        
        pos = min( pos_1, pos_2 );
        
        rho(i,j) = (rho_1+rho_2)/2 + (rho_2-rho_1)/2*tanh( 2*(radius-pos)/5 );
    end
end
%*************************************************************************


%*************************************************************************
%INICIALIZACAO DO CAMPO DE FORCA
%*************************************************************************
% Calculo do pseudopotencial
Peos = Pressure_EOS(rho,T,a,b,R);
psi = sqrt( 2*( Peos - rho*cs2 )/G/c^2 );
% Calculo da Forca de Shan-Chen
F_SC_x = 0;
F_SC_y = 0;
for k = 1:9
    F_SC_x = F_SC_x - G./dt.*psi.*w_SC(k).*cx(k).*circshift(psi,[-cx(k),-cy(k),0]);
    F_SC_y = F_SC_y - G./dt.*psi.*w_SC(k).*cy(k).*circshift(psi,[-cx(k),-cy(k),0]);
end
Fx = F_SC_x;
Fy = F_SC_y;
% Termo fonte para tensao superficial
term_1 = 0;
term_2_xx = 0;
term_2_xy = 0;
term_2_yy = 0;
for k = 1:9
    term_1 = term_1 + 2*psi.*w_SC(k).*( circshift(psi,[-cx(k),-cy(k),0]) - psi );
    term_2_xx = term_2_xx + 3*psi.*w_SC(k).*( cx(k).*cx(k) - cs2 ).*circshift(psi,[-cx(k),-cy(k),0]);
    term_2_xy = term_2_xy + 3*psi.*w_SC(k).*( cx(k).*cy(k) ).*circshift(psi,[-cx(k),-cy(k),0]);
    term_2_yy = term_2_yy + 3*psi.*w_SC(k).*( cy(k).*cy(k) - cs2 ).*circshift(psi,[-cx(k),-cy(k),0]);
end
A = 5/12*( kappa_1 - 1 )*c^4;
B = 1/6*( 1 - kappa_1 )*c^4;
Qxx = G*( A*term_1 + B*term_2_xx );
Qxy = G*( B*term_2_xy );
Qyy = G*( A*term_1 + B*term_2_yy );
%
%Qxx = 0;
%Qxy = 0;
%Qyy = 0;
for k = 1:9
    Qxx = Qxx + kappa_2*G./2.*psi.*w_SC(k).*cx(k).*cx(k).*( circshift(psi,[-cx(k),-cy(k),0]) - psi );
    Qxy = Qxy + kappa_2*G./2.*psi.*w_SC(k).*cx(k).*cy(k).*( circshift(psi,[-cx(k),-cy(k),0]) - psi );
    Qyy = Qyy + kappa_2*G./2.*psi.*w_SC(k).*cy(k).*cy(k).*( circshift(psi,[-cx(k),-cy(k),0]) - psi );
end
%
%*************************************************************************


%*************************************************************************
% Inicializacao da Funcao de Distribuicao
%*************************************************************************
t2 = Ux.*Ux + Vy.*Vy;
for k = 1:9
    t1 = cx(k)*Ux + cy(k)*Vy;
    f(:,:,k) = w(k).*rho.*( 1 + 3.*t1 + 4.5.*t1.*t1 - 1.5.*t2 );
end
%*************************************************************************


%*************************************************************************
% Loop Principal
%*************************************************************************
contador = 1;
contador_2 = 1;
contador_3 = 1;
tempo = 10000;
t_print = 6000;

n_step = 100;

err = 1;
tol = 1e-6;
rho_old = rho;

while ( contador <= tempo )

    %---------------------------------------------------------------------
    % Calculo do Esquema de Forcas
    %---------------------------------------------------------------------
    mf(:,:,1) = 0;
    mf(:,:,2) = 6*( Ux.*Fx + Vy.*Fy ) + 12*sigma*(Fx.^2+Fy.^2)./psi.^2./(tau_e-0.5);
    mf(:,:,3) = - 6*( Ux.*Fx + Vy.*Fy ) - 12*sigma*(Fx.^2+Fy.^2)./psi.^2/(tau_s-0.5);
    mf(:,:,4) = Fx;
    mf(:,:,5) = - Fx;
    mf(:,:,6) = Fy;
    mf(:,:,7) = - Fy;
    mf(:,:,8) = 2*( Ux.*Fx - Vy.*Fy );
    mf(:,:,9) = Ux.*Fy + Vy.*Fx;
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Calculo do Termo Fonte
    %---------------------------------------------------------------------
    % Esquema do termo fonte
    mC(:,:,1) = 0;
    mC(:,:,2) = 1.5*lambda(2)*( Qxx + Qyy );
    mC(:,:,3) = - 1.5*lambda(3)*( Qxx + Qyy );
    mC(:,:,4) = 0;
    mC(:,:,5) = 0;
    mC(:,:,6) = 0;
    mC(:,:,7) = 0;
    mC(:,:,8) = - lambda(8)*( Qxx - Qyy );
    mC(:,:,9) = - lambda(9)*Qxy;
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Calculo da Funcao de Equilibrio
    %---------------------------------------------------------------------
    meq(:,:,1) = rho;
    meq(:,:,2) = rho.*( - 2 + 3*( Ux.*Ux + Vy.*Vy ) );
    meq(:,:,3) = rho.*( 1 - 3*( Ux.*Ux + Vy.*Vy ) );
    meq(:,:,4) = rho.*Ux;
    meq(:,:,5) = - rho.*Ux;
    meq(:,:,6) = rho.*Vy;
    meq(:,:,7) = - rho.*Vy;
    meq(:,:,8) = rho.*( Ux.*Ux - Vy.*Vy );
    meq(:,:,9) = rho.*Ux.*Vy;
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Conversao das funcoes para momentos
    %---------------------------------------------------------------------
    m(:,:,1) = f(:,:,1) + f(:,:,2) + f(:,:,3) + f(:,:,4) + f(:,:,5) + f(:,:,6) + f(:,:,7) + f(:,:,8) + f(:,:,9);
    m(:,:,2) = -4*f(:,:,1) - f(:,:,2) - f(:,:,3) - f(:,:,4) - f(:,:,5) + 2*f(:,:,6) + 2*f(:,:,7) + 2*f(:,:,8) + 2*f(:,:,9);
    m(:,:,3) = 4*f(:,:,1) - 2*f(:,:,2) - 2*f(:,:,3) - 2*f(:,:,4) - 2*f(:,:,5) + f(:,:,6) + f(:,:,7) + f(:,:,8) + f(:,:,9);
    m(:,:,4) = f(:,:,2) - f(:,:,4) + f(:,:,6) - f(:,:,7) - f(:,:,8) + f(:,:,9);
    m(:,:,5) = - 2*f(:,:,2) + 2*f(:,:,4) + f(:,:,6) - f(:,:,7) - f(:,:,8) + f(:,:,9);
    m(:,:,6) = f(:,:,3) - f(:,:,5) + f(:,:,6) + f(:,:,7) - f(:,:,8) - f(:,:,9);
    m(:,:,7) = -2*f(:,:,3) + 2*f(:,:,5) + f(:,:,6) + f(:,:,7) - f(:,:,8) - f(:,:,9);
    m(:,:,8) = f(:,:,2) - f(:,:,3) + f(:,:,4) - f(:,:,5);
    m(:,:,9) = f(:,:,6) - f(:,:,7) + f(:,:,8) - f(:,:,9);
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Colisao
    %---------------------------------------------------------------------
    m(:,:,1) = m(:,:,1) - lambda(1).*( m(:,:,1) - meq(:,:,1) ) + mf(:,:,1) - lambda(1)./2.*mf(:,:,1) + mC(:,:,1);
    m(:,:,2) = m(:,:,2) - lambda(2).*( m(:,:,2) - meq(:,:,2) ) + mf(:,:,2) - lambda(2)./2.*mf(:,:,2) + mC(:,:,2);
    m(:,:,3) = m(:,:,3) - lambda(3).*( m(:,:,3) - meq(:,:,3) ) + mf(:,:,3) - lambda(3)./2.*mf(:,:,3) + mC(:,:,3);
    m(:,:,4) = m(:,:,4) - lambda(4).*( m(:,:,4) - meq(:,:,4) ) + mf(:,:,4) - lambda(4)./2.*mf(:,:,4) + mC(:,:,4);
    m(:,:,5) = m(:,:,5) - lambda(5).*( m(:,:,5) - meq(:,:,5) ) + mf(:,:,5) - lambda(5)./2.*mf(:,:,5) + mC(:,:,5);
    m(:,:,6) = m(:,:,6) - lambda(6).*( m(:,:,6) - meq(:,:,6) ) + mf(:,:,6) - lambda(6)./2.*mf(:,:,6) + mC(:,:,6);
    m(:,:,7) = m(:,:,7) - lambda(7).*( m(:,:,7) - meq(:,:,7) ) + mf(:,:,7) - lambda(7)./2.*mf(:,:,7) + mC(:,:,7);
    m(:,:,8) = m(:,:,8) - lambda(8).*( m(:,:,8) - meq(:,:,8) ) + mf(:,:,8) - lambda(8)./2.*mf(:,:,8) + mC(:,:,8);
    m(:,:,9) = m(:,:,9) - lambda(9).*( m(:,:,9) - meq(:,:,9) ) + mf(:,:,9) - lambda(9)./2.*mf(:,:,9) + mC(:,:,9);
    %---------------------------------------------------------------------

    
    %---------------------------------------------------------------------
    % Reconversao dos momentos em funcoes
    %---------------------------------------------------------------------
    m4 = 1/4;
    m6 = 1/6;
    m9 = 1/9;
    m12 = 1/12;
    m18 = 1/18;
    m36 = 1/36;
    f(:,:,1) = m9*m(:,:,1) - m9*m(:,:,2) + m9*m(:,:,3);
    f(:,:,2) = m9*m(:,:,1) - m36*m(:,:,2) - m18*m(:,:,3) + m6*m(:,:,4) - m6*m(:,:,5) + m4*m(:,:,8);
    f(:,:,3) = m9*m(:,:,1) - m36*m(:,:,2) - m18*m(:,:,3) + m6*m(:,:,6) - m6*m(:,:,7) - m4*m(:,:,8);
    f(:,:,4) = m9*m(:,:,1) - m36*m(:,:,2) - m18*m(:,:,3) - m6*m(:,:,4) + m6*m(:,:,5) + m4*m(:,:,8);
    f(:,:,5) = m9*m(:,:,1) - m36*m(:,:,2) - m18*m(:,:,3) - m6*m(:,:,6) + m6*m(:,:,7) - m4*m(:,:,8);
    f(:,:,6) = m9*m(:,:,1) + m18*m(:,:,2) + m36*m(:,:,3) + m6*m(:,:,4) + m12*m(:,:,5) + m6*m(:,:,6) + m12*m(:,:,7) + m4*m(:,:,9);
    f(:,:,7) = m9*m(:,:,1) + m18*m(:,:,2) + m36*m(:,:,3) - m6*m(:,:,4) - m12*m(:,:,5) + m6*m(:,:,6) + m12*m(:,:,7) - m4*m(:,:,9);
    f(:,:,8) = m9*m(:,:,1) + m18*m(:,:,2) + m36*m(:,:,3) - m6*m(:,:,4) - m12*m(:,:,5) - m6*m(:,:,6) - m12*m(:,:,7) + m4*m(:,:,9);
    f(:,:,9) = m9*m(:,:,1) + m18*m(:,:,2) + m36*m(:,:,3) + m6*m(:,:,4) + m12*m(:,:,5) - m6*m(:,:,6) - m12*m(:,:,7) - m4*m(:,:,9);
    %---------------------------------------------------------------------
    

    %---------------------------------------------------------------------
    % Propagacao
    %---------------------------------------------------------------------
    for k = 1:9
        f(:,:,k) = circshift(f(:,:,k),[cx(k),cy(k),0]);
    end
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Calculo da Massa
    %---------------------------------------------------------------------
    rho = sum(f,3);
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Calculo da Forca
    %---------------------------------------------------------------------
    % Calculo do pseudopotencial
    Peos = Pressure_EOS(rho,T,a,b,R);
    psi = sqrt( 2*( Peos - rho*cs2 )/G/c^2 );
    % Calculo da Forca de Shan-Chen
    F_SC_x = 0;
    F_SC_y = 0;
    for k = 1:9
        F_SC_x = F_SC_x - G./dt.*psi.*w_SC(k).*cx(k).*circshift(psi,[-cx(k),-cy(k),0]);
        F_SC_y = F_SC_y - G./dt.*psi.*w_SC(k).*cy(k).*circshift(psi,[-cx(k),-cy(k),0]);
    end
    Fx = F_SC_x;
    Fy = F_SC_y;
    % Termo fonte para tensao superficial
    term_1 = 0;
    term_2_xx = 0;
    term_2_xy = 0;
    term_2_yy = 0;
    for k = 1:9
        term_1 = term_1 + 2*psi.*w_SC(k).*( circshift(psi,[-cx(k),-cy(k),0]) - psi );
        term_2_xx = term_2_xx + 3*psi.*w_SC(k).*( cx(k).*cx(k) - cs2 ).*circshift(psi,[-cx(k),-cy(k),0]);
        term_2_xy = term_2_xy + 3*psi.*w_SC(k).*( cx(k).*cy(k) ).*circshift(psi,[-cx(k),-cy(k),0]);
        term_2_yy = term_2_yy + 3*psi.*w_SC(k).*( cy(k).*cy(k) - cs2 ).*circshift(psi,[-cx(k),-cy(k),0]);
    end
    A = 5/12*( kappa_1 - 1 )*c^4;
    B = 1/6*( 1 - kappa_1 )*c^4;
    Qxx = G*( A*term_1 + B*term_2_xx );
    Qxy = G*( B*term_2_xy );
    Qyy = G*( A*term_1 + B*term_2_yy );
    %
    %Qxx = 0;
    %Qxy = 0;
    %Qyy = 0;
    for k = 1:9
        Qxx = Qxx + kappa_2*G./2.*psi.*w_SC(k).*cx(k).*cx(k).*( circshift(psi,[-cx(k),-cy(k),0]) - psi );
        Qxy = Qxy + kappa_2*G./2.*psi.*w_SC(k).*cx(k).*cy(k).*( circshift(psi,[-cx(k),-cy(k),0]) - psi );
        Qyy = Qyy + kappa_2*G./2.*psi.*w_SC(k).*cy(k).*cy(k).*( circshift(psi,[-cx(k),-cy(k),0]) - psi );
    end
    %
    
    % Forca de atracao
    if contador < 2000
        g = 0;
    else 
        g = 0.5e-5;
    end
    rhom = sum(sum(rho))/Nx/Ny;
    Fg = [ ( rho(:,1:Ny/2) - rhom )*g, -( rho(:,Ny/2+1:Ny) - rhom )*g ];

    Fy = Fy + Fg;
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Calculo das Variaveis
    %---------------------------------------------------------------------
    Ux = ( sum(f(:,:,[2 6 9]),3) - sum(f(:,:,[4 7 8]),3) + 0.5*Fx )./rho;
    Vy = ( sum(f(:,:,[3 6 7]),3) - sum(f(:,:,[5 8 9]),3) + 0.5*Fy )./rho;
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Plotar Figura
    %---------------------------------------------------------------------
    if mod(contador,n_step) == 0
        % erro
        err = sum( sum( abs(rho-rho_old)./rho_old ) )/Nx/Ny;
        rho_old = rho;
        disp('erro =');
        disp(err);
    end
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Plotar Figura
    %---------------------------------------------------------------------
    if mod(contador,400) == 0
        figure(1)
        contourf(rho)
        pbaspect([Ny Nx 1])
        drawnow
        %
        disp(contador);
        rho_leve = rho(1,1);
        rho_pesado = rho(Nx/2,Ny/2);
    end
    %---------------------------------------------------------------------
    
    %---------------------------------------------------------------------
    % Guardar Resultado
    %---------------------------------------------------------------------
    if mod(contador,1000) == 0
        Est_1(:,:,contador_2) = rho;
        contador_2 = contador_2 + 1;
    end
    
    if mod(contador,t_print) == 0
        Est_2(:,:,contador_3) = rho;
        t_print = t_print + 20;
        contador_3 = contador_3 + 1;
        if contador_3 > 41
            contador_3 = 41;
        end
    end
    %---------------------------------------------------------------------

contador = contador + 1;

%*************************************************************************
%*************************************************************************

end

% Salvando estoque de resultados
%save('Coalescencia_Pseudopotencial','Est_1','Est_2','a','b','R','epsilon','kappa_1','kappa_2')

%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
% Equacao de Estado
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
function [Peos] = Pressure_EOS(rho,T,a,b,R)

n = b*rho/4;

T1 = rho.*R.*T.*( 1 + n + n.^2 - n.^3 )./( 1 - n ).^3;

Peos = T1 - a*rho.^2;

end
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


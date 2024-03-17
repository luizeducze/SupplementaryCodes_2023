% Pseudopotencial abordagem hibrida (Forca + Pressao)

clear all
close all

% Geometria do dominio
Nx = 1; % largura do dominio
Ny = 100; % altura do dominio
% Tempo de relaxacao
tau_v = 1;
tau_b = 1;
tau_0 = 1;
tau_2 = 1;
tau_4 = 1;
lambda = [tau_0^(-1) tau_b^(-1) tau_2^(-1) tau_0^(-1) tau_4^(-1) tau_0^(-1) tau_4^(-1) tau_v^(-1) tau_v^(-1) ];
%A = 0.3;
A = 0;
% Matriz de relaxacao de Gram-Schmidt
M = [ 1, 1, 1, 1, 1, 1, 1, 1, 1; ...
     -4,-1,-1,-1,-1, 2, 2, 2, 2; ...
      4,-2,-2,-2,-2, 1, 1, 1, 1; ...
      0, 1, 0,-1, 0, 1,-1,-1, 1; ...
      0,-2, 0, 2, 0, 1,-1,-1, 1; ...
      0, 0, 1, 0,-1, 1, 1,-1,-1; ...
      0, 0,-2, 0, 2, 1, 1,-1,-1; ...
      0, 1,-1, 1,-1, 0, 0, 0, 0; ...
      0, 0, 0, 0, 0, 1,-1, 1,-1];
Minv = inv(M);
% Geometria da interface plana
radius = 25; % largura da fase liquido
thickness = 5;
% Parametros da equacao de estado
a = 0.5;
b = 4;
R = 1;
T = 0.65 * ( a / b / R / 2.6502 );
% Parametros do modelo multifasico
rho_1 = 0.0055881;
rho_2 = 0.3823228;
kappa = 0.5;
% Parametros do lattice
dx = 1; % unidades de lattice
dt = 1; % unidades de lattice
w = [4/9 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
cx = [0, 1, 0,-1, 0, 1,-1,-1, 1];
cy = [0, 0, 1, 0,-1, 1, 1,-1,-1];
cs2 = 1/3;
% Campos de densidade e velocidade
rho = zeros(Nx,Ny);
Ux = zeros(Nx,Ny);
Vy = zeros(Nx,Ny);
% Funcao de distribuicao
f = zeros(Nx,Ny,9);
feq = zeros(Nx,Ny,9);
Fk = zeros(Nx,Ny,9);
m = zeros(Nx,Ny,9);
meq = zeros(Nx,Ny,9);
mf = zeros(Nx,Ny,9);

%*************************************************************************
% Inicializacao da Densidade
%*************************************************************************
for i = 1:Nx
    for j = 1:Ny
        x_pos = 0;
        %x_pos = i - 0.5 - Nx/2;
        y_pos = j - 0.5 - Ny/2;
        pos = sqrt( x_pos^2 + y_pos^2 );
        
        rho(i,j) = (rho_1+rho_2)/2 + (rho_2-rho_1)/2*tanh( 2*(radius-pos)/5 );
        %rho(i,j) = (rho_1+rho_2)/2 + (rho_2-rho_1)/2*tanh( (radius-pos)/sqrt(2)/(epsilon) );
    end
end
%*************************************************************************

%*************************************************************************
% Inicializacao da Forca
%*************************************************************************
% Laplaciano da densidade
d2_rho = 0;
for k = 1:9
    d2_rho = d2_rho + 2./cs2./(dx.^2).*w(k).*( circshift(rho,[-cx(k),-cy(k),0]) - rho );
end

mu_0 = 0;
%v_rho = (rho-rho_c)/rho_c;
%mu = 4*p_c/rho_c*v_rho.*(v_rho.^2-beta_x_tauw) + mu_0 - kappa*d2_rho;
mu_b = chemical_potential( rho, a, b, R, T );
mu = mu_b + mu_0 - kappa*d2_rho;

% Derivadas do potencial quimico
dx_mu = 0;
dy_mu = 0;
dx_rho = 0;
dy_rho = 0;
for k = 1:9
    dx_rho = dx_rho + 1./cs2./(dx.^2./dt).*w(k).*cx(k).*circshift(rho,[-cx(k),-cy(k),0]);
    dy_rho = dy_rho + 1./cs2./(dx.^2./dt).*w(k).*cy(k).*circshift(rho,[-cx(k),-cy(k),0]);
    dx_mu = dx_mu + 1./cs2./(dx.^2./dt).*w(k).*cx(k).*circshift(mu,[-cx(k),-cy(k),0]);
    dy_mu = dy_mu + 1./cs2./(dx.^2./dt).*w(k).*cy(k).*circshift(mu,[-cx(k),-cy(k),0]);
end
Fx = - rho.*dx_mu;
Fy = - rho.*dy_mu;
%*************************************************************************

%*************************************************************************
% Inicializacao da Funcao de Distribuicao
%*************************************************************************
t2 = Ux.*Ux + Vy.*Vy;
for k = 1:9
    t1 = cx(k)*Ux + cy(k)*Vy;
    feq(:,:,k) = w(k).*rho.*( 1 + 3.*t1 + 4.5.*t1.*t1 - 1.5.*t2 );
    f(:,:,k) = feq(:,:,k);
end
%*************************************************************************

%load('solution_50K','rho','Ux','Vy','f','feq','Fx','Fy','Fk')
%save('solution_50K','rho','Ux','Vy','f','feq','Fx','Fy','Fk')

%*************************************************************************
% Loop Principal
%*************************************************************************
contador = 1;
tempo = 50000;
err = 1;
tol = 1e-6;
rho_old = rho;

while ( contador <= tempo )
    

    %---------------------------------------------------------------------
    % Calculo do Esquema de Forcas
    %---------------------------------------------------------------------
    for k = 1:9
        c_k2 = ( cx(k)^2 + cy(k)^2 );
        
        Fk(:,:,k) = w(k)*( cx(k)/cs2*Fx + cy(k)/cs2*Fy ...
            + ( cx(k)*cx(k) - cs2 )/cs2^2*( Ux.*( Fx + cs2*dx_rho ) ) ...
            + ( cx(k)*cy(k) )/cs2^2*( Ux.*( Fy + cs2*dy_rho ) ...
            + Vy.*( Fx + cs2*dx_rho ) ) ...
            + ( cy(k)*cy(k) - cs2 )/cs2^2*( Vy.*( Fy + cs2*dy_rho ) ) ...
            + 0.5*( c_k2/cs2 - 2 )*( Ux.*dx_rho + Vy.*dy_rho ) );
    end
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Calculo da Funcao de Equilibrio
    %---------------------------------------------------------------------
    % Derivadas da velocidade
    dx_Ux = 0;
    dy_Ux = 0;
    dx_Vy = 0;
    dy_Vy = 0;
    for k = 1:9
        dx_Ux = dx_Ux + 1./cs2./(dx.^2./dt).*w(k).*cx(k).*circshift(Ux,[-cx(k),-cy(k),0]);
        dy_Ux = dy_Ux + 1./cs2./(dx.^2./dt).*w(k).*cy(k).*circshift(Ux,[-cx(k),-cy(k),0]);
        dx_Vy = dx_Vy + 1./cs2./(dx.^2./dt).*w(k).*cx(k).*circshift(Vy,[-cx(k),-cy(k),0]);
        dy_Vy = dy_Vy + 1./cs2./(dx.^2./dt).*w(k).*cy(k).*circshift(Vy,[-cx(k),-cy(k),0]);
    end
    
    Sxx = dx_Ux + dx_Ux;
    Sxy = dy_Ux + dx_Vy;
    %Syx = Sxy;
    Syy = dy_Vy + dy_Vy;
    
    for k = 2:9
        feq(:,:,k) = w(k)*rho.*( cx(k)/cs2*Ux + cy(k)/cs2*Vy ...
            + ( cx(k)*cx(k) - cs2 )/2/cs2^2*( Ux.*Ux + cs2*A*dt*Sxx ) ...
            + ( cx(k)*cy(k) )/cs2^2*( Ux.*Vy + cs2*A*dt*Sxy ) ...
            + ( cy(k)*cy(k) - cs2 )/2/cs2^2*( Vy.*Vy + cs2*A*dt*Syy ) );
    end
    feq(:,:,1) = rho - w(1)*rho.*( Ux.*Ux + Vy.*Vy )/2/cs2 - w(1)*rho.*A.*dt.*( dx_Ux + dy_Vy ); 
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Conversao da funcoes para momentos
    %---------------------------------------------------------------------
    m = zeros(Nx,Ny,9);
    meq = zeros(Nx,Ny,9);
    mf = zeros(Nx,Ny,9);
    %
    m(:,:,1) = f(:,:,1) + f(:,:,2) + f(:,:,3) + f(:,:,4) + f(:,:,5) + f(:,:,6) + f(:,:,7) + f(:,:,8) + f(:,:,9);
    m(:,:,2) = -4*f(:,:,1) - f(:,:,2) - f(:,:,3) - f(:,:,4) - f(:,:,5) + 2*f(:,:,6) + 2*f(:,:,7) + 2*f(:,:,8) + 2*f(:,:,9);
    m(:,:,3) = 4*f(:,:,1) - 2*f(:,:,2) - 2*f(:,:,3) - 2*f(:,:,4) - 2*f(:,:,5) + f(:,:,6) + f(:,:,7) + f(:,:,8) + f(:,:,9);
    m(:,:,4) = f(:,:,2) - f(:,:,4) + f(:,:,6) - f(:,:,7) - f(:,:,8) + f(:,:,9);
    m(:,:,5) = - 2*f(:,:,2) + 2*f(:,:,4) + f(:,:,6) - f(:,:,7) - f(:,:,8) + f(:,:,9);
    m(:,:,6) = f(:,:,3) - f(:,:,5) + f(:,:,6) + f(:,:,7) - f(:,:,8) - f(:,:,9);
    m(:,:,7) = -2*f(:,:,3) + 2*f(:,:,5) + f(:,:,6) + f(:,:,7) - f(:,:,8) - f(:,:,9);
    m(:,:,8) = f(:,:,2) - f(:,:,3) + f(:,:,4) - f(:,:,5);
    m(:,:,9) = f(:,:,6) - f(:,:,7) + f(:,:,8) - f(:,:,9);
    %
    meq(:,:,1) = feq(:,:,1) + feq(:,:,2) + feq(:,:,3) + feq(:,:,4) + feq(:,:,5) + feq(:,:,6) + feq(:,:,7) + feq(:,:,8) + feq(:,:,9);
    meq(:,:,2) = -4*feq(:,:,1) - feq(:,:,2) - feq(:,:,3) - feq(:,:,4) - feq(:,:,5) + 2*feq(:,:,6) + 2*feq(:,:,7) + 2*feq(:,:,8) + 2*feq(:,:,9);
    meq(:,:,3) = 4*feq(:,:,1) - 2*feq(:,:,2) - 2*feq(:,:,3) - 2*feq(:,:,4) - 2*feq(:,:,5) + feq(:,:,6) + feq(:,:,7) + feq(:,:,8) + feq(:,:,9);
    meq(:,:,4) = feq(:,:,2) - feq(:,:,4) + feq(:,:,6) - feq(:,:,7) - feq(:,:,8) + feq(:,:,9);
    meq(:,:,5) = - 2*feq(:,:,2) + 2*feq(:,:,4) + feq(:,:,6) - feq(:,:,7) - feq(:,:,8) + feq(:,:,9);
    meq(:,:,6) = feq(:,:,3) - feq(:,:,5) + feq(:,:,6) + feq(:,:,7) - feq(:,:,8) - feq(:,:,9);
    meq(:,:,7) = -2*feq(:,:,3) + 2*feq(:,:,5) + feq(:,:,6) + feq(:,:,7) - feq(:,:,8) - feq(:,:,9);
    meq(:,:,8) = feq(:,:,2) - feq(:,:,3) + feq(:,:,4) - feq(:,:,5);
    meq(:,:,9) = feq(:,:,6) - feq(:,:,7) + feq(:,:,8) - feq(:,:,9);
    %
    mf(:,:,1) = Fk(:,:,1) + Fk(:,:,2) + Fk(:,:,3) + Fk(:,:,4) + Fk(:,:,5) + Fk(:,:,6) + Fk(:,:,7) + Fk(:,:,8) + Fk(:,:,9);
    mf(:,:,2) = -4*Fk(:,:,1) - Fk(:,:,2) - Fk(:,:,3) - Fk(:,:,4) - Fk(:,:,5) + 2*Fk(:,:,6) + 2*Fk(:,:,7) + 2*Fk(:,:,8) + 2*Fk(:,:,9);
    mf(:,:,3) = 4*Fk(:,:,1) - 2*Fk(:,:,2) - 2*Fk(:,:,3) - 2*Fk(:,:,4) - 2*Fk(:,:,5) + Fk(:,:,6) + Fk(:,:,7) + Fk(:,:,8) + Fk(:,:,9);
    mf(:,:,4) = Fk(:,:,2) - Fk(:,:,4) + Fk(:,:,6) - Fk(:,:,7) - Fk(:,:,8) + Fk(:,:,9);
    mf(:,:,5) = - 2*Fk(:,:,2) + 2*Fk(:,:,4) + Fk(:,:,6) - Fk(:,:,7) - Fk(:,:,8) + Fk(:,:,9);
    mf(:,:,6) = Fk(:,:,3) - Fk(:,:,5) + Fk(:,:,6) + Fk(:,:,7) - Fk(:,:,8) - Fk(:,:,9);
    mf(:,:,7) = -2*Fk(:,:,3) + 2*Fk(:,:,5) + Fk(:,:,6) + Fk(:,:,7) - Fk(:,:,8) - Fk(:,:,9);
    mf(:,:,8) = Fk(:,:,2) - Fk(:,:,3) + Fk(:,:,4) - Fk(:,:,5);
    mf(:,:,9) = Fk(:,:,6) - Fk(:,:,7) + Fk(:,:,8) - Fk(:,:,9);
    %
    %{
    for i = 1:9
        for j = 1:9
            m(:,:,i) = m(:,:,i) + M(i,j)*f(:,:,j);
            meq(:,:,i) = meq(:,:,i) + M(i,j)*feq(:,:,j);
            mf(:,:,i) = mf(:,:,i) + M(i,j)*Fk(:,:,j);
        end
    end
    %}
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Colisao
    %---------------------------------------------------------------------
    m(:,:,1) = m(:,:,1) - lambda(1).*(m(:,:,1)-meq(:,:,1)) + mf(:,:,1)-lambda(1)./2.*mf(:,:,1);
    m(:,:,2) = m(:,:,2) - lambda(2).*(m(:,:,2)-meq(:,:,2)) + mf(:,:,2)-lambda(2)./2.*mf(:,:,2);
    m(:,:,3) = m(:,:,3) - lambda(3).*(m(:,:,3)-meq(:,:,3)) + mf(:,:,3)-lambda(3)./2.*mf(:,:,3);
    m(:,:,4) = m(:,:,4) - lambda(4).*(m(:,:,4)-meq(:,:,4)) + mf(:,:,4)-lambda(4)./2.*mf(:,:,4);
    m(:,:,5) = m(:,:,5) - lambda(5).*(m(:,:,5)-meq(:,:,5)) + mf(:,:,5)-lambda(5)./2.*mf(:,:,5);
    m(:,:,6) = m(:,:,6) - lambda(6).*(m(:,:,6)-meq(:,:,6)) + mf(:,:,6)-lambda(6)./2.*mf(:,:,6);
    m(:,:,7) = m(:,:,7) - lambda(7).*(m(:,:,7)-meq(:,:,7)) + mf(:,:,7)-lambda(7)./2.*mf(:,:,7);
    m(:,:,8) = m(:,:,8) - lambda(8).*(m(:,:,8)-meq(:,:,8)) + mf(:,:,8)-lambda(8)./2.*mf(:,:,8);
    m(:,:,9) = m(:,:,9) - lambda(9).*(m(:,:,9)-meq(:,:,9)) + mf(:,:,9)-lambda(9)./2.*mf(:,:,9);
    %---------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------
    % Reconversao dos momentos em funcoes
    %---------------------------------------------------------------------
    f = zeros(Nx,Ny,9);
    %
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
    %
    %{
    for i = 1:9
        for j = 1:9
            f(:,:,i) = f(:,:,i) + Minv(i,j)*m(:,:,j);
        end
    end
    %}
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
    % Laplaciano da densidade
    d2_rho = 0;
    for k = 1:9
        d2_rho = d2_rho + 2./cs2./(dx.^2).*w(k).*( circshift(rho,[-cx(k),-cy(k),0]) - rho );
    end
    
    mu_0 = 0;
    %v_rho = (rho-rho_c)/rho_c;
    %mu = 4*p_c/rho_c*v_rho.*(v_rho.^2-beta_x_tauw) + mu_0 - kappa*d2_rho;
    mu_b = chemical_potential( rho, a, b, R, T );
    mu = mu_b + mu_0 - kappa*d2_rho;
    
    % Derivadas do potencial quimico
    dx_mu = 0;
    dy_mu = 0;
    dx_rho = 0;
    dy_rho = 0;
    for k = 1:9
        dx_rho = dx_rho + 1./cs2./(dx.^2./dt).*w(k).*cx(k).*circshift(rho,[-cx(k),-cy(k),0]);
        dy_rho = dy_rho + 1./cs2./(dx.^2./dt).*w(k).*cy(k).*circshift(rho,[-cx(k),-cy(k),0]);
        dx_mu = dx_mu + 1./cs2./(dx.^2./dt).*w(k).*cx(k).*circshift(mu,[-cx(k),-cy(k),0]);
        dy_mu = dy_mu + 1./cs2./(dx.^2./dt).*w(k).*cy(k).*circshift(mu,[-cx(k),-cy(k),0]);
    end
    Fx = - rho.*dx_mu;
    Fy = - rho.*dy_mu;
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
    if mod(contador,100) == 0
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
    if mod(contador,100) == 0
        figure(1)
        y = linspace(0.5,Ny-0.5,Ny);
        %rho_teorico = rho_c + rho_c*sqrt( beta_x_tauw )*tanh( ( radius - sqrt( (y-Ny/2).^2 ) )/sqrt(2)/epsilon );
        %plot(y,rho(1,:),'ko',y,rho_teorico,'r');
        plot(y,rho(1,:),'ko')
        %contourf(rho)
        %pbaspect([Nx Ny 1])
        drawnow
        %
        disp(contador);
        rho_leve = rho(1,1)
        rho_pesado = rho(1,Ny/2)
    end
    %---------------------------------------------------------------------
    
    contador = contador + 1;
    
end
%*************************************************************************

erro = 100*abs( rho_leve - rho_1 )/rho_1

%{
rho_in = rho(Nx/2,Ny/2)
rho_out = rho(1,1)

Pin = Pressure_EOS(rho_in,rho_c,p_c,beta_x_tauw);
Pout = Pressure_EOS(rho_out,rho_c,p_c,beta_x_tauw);
% Raio
rhom = (rho_1+rho_2)/2;
A = contour(rho,[rhom rhom]);
pbaspect([Ny Nx 1])
diam = max(A(2,2:end)) - min(A(2,2:end));
dP = Pin - Pout
raio_num = diam/2

Surf_Tension = (Pin-Pout)*raio_num
%}

%save('Perfil_Numerico_7','y','rho','rho_teorico','Ny','epsilon');
%

%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
function [mu_b] = chemical_potential( rho, a, b, R, T )
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
n = rho*b/4;

aux_1 = R*T*( 1 + n + n.^2 - n.^3 )./( 1 - n ).^3;
aux_2 = - 2*a*rho;
aux_3 = R*T*log(rho);
aux_4 = 2*R*T./( 1 - n );
aux_5 = R*T./( 1 - n ).^2;

mu_b = aux_1 + aux_2 + aux_3 + aux_4 + aux_5;

end
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
% Equacao de Estado
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
function [Peos] = Pressure_EOS(rho,rho_c,p_c,beta_x_tauw)

v_rho = (rho-rho_c)/rho_c;

Peos = p_c*( v_rho + 1 ).^2.*( 3*v_rho.^2 - 2*v_rho + 1 - 2*beta_x_tauw );

end
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo





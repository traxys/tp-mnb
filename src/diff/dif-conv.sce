Nx=200;
Nt=200;
kappa=0.001;
exec("my_cholesky.sce")

function y=phi_0(x)
	if x < 0.25
		y = 0
	else 
		if x < 0.375
			y = 2*(x-0.25)
		else 
			if x < 0.5
				y = 2*(0.5-x)
			else
				y = 0
			end
		end
	end
endfunction

function y=conv(x)
	y = 0.4*(x-0.25)
endfunction

phi_i = zeros(Nx, 1)
for i = 0:(Nx-1)
	phi_i(i+1) = phi_0(i/Nx)
end
phi = phi_i

//Coeficients hors de la diagonale de N
N_nd = -kappa*Nx^2/Nt
//Coeficients diagonaux de N
N_d = 1+2*kappa*Nx^2/Nt

//Coefficents plus que diagonaux de M
function y=M_pd(x)
	y = conv(x)*Nx/2/Nt*(conv(x)*Nx/Nt-1)
endfunction
//Coefficents moins que diagonaux de M
function y=M_ld(x)
	y = conv(x)*Nx/2/Nt*(conv(x)*Nx/Nt+1)
endfunction
//Coefficients diagonaux de M
function y=M_d(x)
	y = 1 - conv(x)^2*Nx^2/Nt^2
endfunction

N = zeros(Nx, Nx)
N(1,1) = N_d
N(1,2) = N_nd
N(1, Nx) = N_nd

N(Nx,Nx) = N_d
N(Nx, 1) = N_nd
N(Nx, Nx-1) = N_nd

row = [N_nd, N_d, N_nd]
for i = 2:(Nx-1)
	N(i, (i-1):(i+1)) = row
end

M = zeros(Nx, Nx)
M(1,1) = M_d(0)
M(1,2) = M_pd(0)
M(1,Nx) = M_ld(0)

last_x = (Nx-1)/Nx
M(Nx, Nx) = M_d(last_x)
M(Nx, Nx-1) = M_ld(last_x)
M(Nx, 1) = M_pd(last_x)

for i = 2:(Nx-1)
	x = i/Nx
	M(i,i-1) = M_ld(x)
	M(i,i) = M_d(x)
	M(i,i+1) = M_pd(x)
end

fin=Nt;
for i=1:fin
    phi=my_cholesky(N,M*phi);
end

maillage = linspace(0,1,Nx)
scf;
plot(maillage,[phi_i phi]);

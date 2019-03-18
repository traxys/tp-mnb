
function dt=calcul_dt(U,dx)
  dt=dx/max(abs(U));
endfunction

function vort=solveur_1D(vort, Ux, Nx, kappa, dt, dx)
	//Coeficients hors de la diagonale de N
	N_nd = -kappa*Nx^2*dt
	//Coeficients diagonaux de N
	N_d = 1+2*kappa*Nx^2*dt

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
	
	//Coefficents plus que diagonaux de M
	M_pd = (Ux*Nx/2*dt).*(Ux*Nx*dt-1)
	//Coefficents moins que diagonaux de M
	M_ld = (Ux*Nx/2*dt).*(Ux*Nx*dt+1)
	//Coefficients diagonaux de M
	M_d = 1 - Ux^2*Nx^2*dt^2

	M = zeros(Nx, Nx)
	M(1,1) = M_d(1)
	M(1,2) = M_pd(1)
	M(1,Nx) = M_ld(1)

	M(Nx, Nx) = M_d(Nx-1)
	M(Nx, Nx-1) = M_ld(Nx-1)
	M(Nx, 1) = M_pd(Nx-1)

	for i = 2:(Nx-1)
		M(i,i-1) = M_ld(i)
		M(i,i) = M_d(i)
		M(i,i+1) = M_pd(i)
	end
	vort = umfpack(sparse(N),"\",M*vort)
endfunction

function vort=solveur_2D(vort, Ux, Uy, Nx, Ny, kappa, dt, dx, dy)
	for y = 1:Ny
		trans_vort = vort(y,:)'
		vort(y,:) = solveur_1D(trans_vort, Ux(y,:), Nx, kappa, dt, dx)'
	end
	for x = 1:Nx
		vort(:,x) = solveur_1D(vort(:,x), Uy(:,y), Ny, kappa, dt, dy)
	end
endfunction

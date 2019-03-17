
dx=Lx/Nx;
dy=Ly/Ny;

exec("dif-conv-f.sce")

maillage_x=linspace(0,(Nx-1)/Nx*Lx,Nx)';
maillage_y=linspace(0,(Ny-1)/Ny*Ly,Ny)';

function [cx, cy, phi_i] = get_data(Nx, Ny)

	cx=zeros(Ny,Nx); //composante x de la vitesse de convection
	cy=zeros(Ny,Nx); //composante y de la vitesse de convection

	phi=zeros(Ny,Nx);   //fonction à calculer
	phi_i=zeros(Ny,Nx); //condtion initiale

	for x = 1:Nx
		for y = 1:Ny
			c = conv(y*dy,x*dx)
			cx(y,x) = c(2)
			cy(y,x) = c(1)
		
			phi_i(y,x) = phi_0(y*dy,x*dx)
		end
	end
endfunction

function phi = advance_time(cx, cy, Nx, Ny, nu, phi_i, Tf)
	phi = phi_i
	dt=min(calcul_dt(cx,dx),calcul_dt(cy,dy));
	Nt=floor(Tf/dt);
	for k=1:Nt
    	phi = solveur_2D(phi, cx, cy, Nx, Ny, nu, dt, dx, dy)
		disp(k)
	end
endfunction


// Retourne la fréquence d'échantillonage de la transformée de Fourier discrète
function [freqs]=fftfreq(N, L)
    // TODO: Calculer les fréquences d'échantillonage en fonction de L et de la parité de N
	if modulo(N, 2) == 0 then
		k = cat(2,linspace(0, N/2-1, N/2), linspace(-N/2, -1, N/2))
	else
		k = cat(2,linspace(0, (N-1)/2, (N+1)/2), linspace(-(N-1)/2, -1, (N-1)/2))
	end
	freqs = 2*%i*%pi/L * k
endfunction


// Résolution de l'équation de Poisson en dimension 2 en utilisant la FFT
//    laplacien(psi) = f
// Entrée: f de taille (Ny,Nx) sur une domaine de taille (Ly,Lx)
// Sortie: psi, solution de l'équation
function [psi]=poisson_2d(f, Nx, Ny, Lx, Ly)
    // calcul des nombres d'ondes
    kx = fftfreq(Nx, Lx)
    ky = fftfreq(Ny, Ly)
    // transformation de fourier de f
    f_hat = fft(f, "nonsymmetric")
    // calcul de psi_hat
    psi_hat = zeros(Ny, Nx)
    for p=1:Ny
        for q=1:Nx
            if (p == 1 & q == 1) then
                psi_hat(p,q) = 0
            else
                psi_hat(p,q) = f_hat(p,q)/(kx(q)^2 + ky(p)^2)
            end
        end
    end
    // psi_hat(1,1) a été laissé à 0
    // transformation de fourier inverse de psi_hat
    psi = real(ifft(psi_hat, "nonsymmetric"))
endfunction

// Résolution de l'équation de Poisson avec rot en dimension 2 en utilisant la FFT
//    laplacien(Ux) = -dW/dy
//    laplacien(Uy) = +dW/dx
// Entrée: champs de vorticité W de taille (Ny,Nx) sur un domaine de taille (Ly,Lx)
// Sortie: Ux et Uy, vitesses solution des équations
function [Ux,Uy]=poisson_curl_2d(W, Nx, Ny, Lx, Ly)
    // calcul des nombres d'ondes
    kx = fftfreq(Nx, Lx)
    ky = fftfreq(Ny, Ly)
    // transformation de fourier de la vorticité
    w_hat = fft(W, "nonsymmetric")
    //calcul de Nx_hat et Ny_hat
    Ux_hat = zeros(Ny, Nx)
    Uy_hat = zeros(Ny, Nx)
    for p=1:Ny
        for q=1:Nx
            if (p==1 & q==1) then
                Ux_hat(p,q) = 0
                Uy_hat(p,q) = 0
            else
                kxq = kx(q)
                kyp = ky(p)
                squares = kxq*kxq + kyp*kyp
                Ux_hat(p,q) = w_hat(p,q)*kyp/squares
                Uy_hat(p,q) = w_hat(p,q)*kxq/squares
            end
    end
    // transformations de fourier inverses
    Ux = real(ifft(Ux_hat))
    Uy = real(ifft(Uy_hat))
endfunction


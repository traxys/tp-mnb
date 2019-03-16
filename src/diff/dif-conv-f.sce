
function dt=calcul_dt(U,dx)
  dt=dx/max(abs(U));
endfunction

function vort=solveur_1D(vort, Ux, Nx, kappa, dt, dx)
    //----------------------------------------------
    //TODO implémenter une itération temporelle 
    // de l'algo 1D codé dans le script dif-conv.sce
    //   vort = champ transporté (correspond à phi dans le sujet)
    //   Ux = vitesse sur la composante x
    //   Nx = nombre de points de discrétisation spatiale en x et y.
    //   kappa = constante de diffusion dynamique
    //   dt = pas de temps
    //   dx = pas d'espace
    //-----------------------------------------------
endfunction

function vort=solveur_2D(vort, Ux, Uy, Nx, Ny, kappa, dt, dx, dy)
    //----------------------------------------------
    //TODO implémenter une itération temporelle 
    // de l'algo de splitting 2D utilisant l'algo 1D
    //   vort = champ 2D transporté (correspond à phi dans le sujet)
    //   Ux = vitesse sur la composante x
    //   Uy = vitesse sur la composante y
    //   (Nx, Ny) nombre de points de discrétisation spatiale en x et y.
    //   kappa = constante de diffusion dynamique
    //   dt    = pas de temps
    //   (dx, dy) pas d'espaces dans chaque direction
    //-----------------------------------------------
endfunction

// Initialise w(x,y)
function [w]=init_field(y,x)
    w = 8*%pi*%pi*cos(2*%pi*x)*cos(2*%pi*y)
endfunction


// Solution de référence du problème laplacien(u(x,y)) = rot(w(x,y))
function [refx]=solution_fieldx(y,x)
    refx = -2*%pi*cos(2*%pi*x)*sin(2*%pi*y)
endfunction

function [refy]=solution_fieldy(y,x)
    refy = 2*%pi*sin(2*%pi*x)*cos(2*%pi*y)
endfunction

// Affichage de la fonction f, de la solution de référence et 
// de la solution obtenue, ainsi que de l'erreur commise
// W   = w(x,y)         -- fonction testée
// Refx= ux_alpha(x,y) -- solution analytique de ux
// Refy= uy_alpha(x,y) -- solution analytique de uy
// Ux = ux_star(x,y)  -- solution du solveur pour ux
// Uy = uy_star(x,y)  -- solution du solveur pour uy
function plot_error(W, Refx, Refy, Ux,Uy)  
    fig = gcf()
    subplot(421)
    colorbar(min(W),max(W))
    plot3d(Y,X,W)
    subplot(423)
    colorbar(min(Refx),max(Refx))
    plot3d(Y,X,Refx)
    subplot(425)
    colorbar(min(Ux),max(Ux))
    plot3d(Y,X,Ux)
    subplot(424)
    colorbar(min(Refy),max(Refy))
    plot3d(Y,X,Refy)
    subplot(426)
    colorbar(min(Uy),max(Uy))
    plot3d(Y,X,Uy)
    subplot(427)
    colorbar(min(abs(Refx-Ux)),max(abs(Refx-Ux)))
    plot3d(Y,X,abs(Refx-Ux))
    subplot(428)
    colorbar(min(abs(Refy-Uy)),max(abs(Refy-Uy)))
    plot3d(Y,X,abs(Refy-Uy))
    xs2png(fig, "poisson_curl_error.png")
endfunction

// Fonction de test pour le solveur de Poisson_curl
function test_poisson_curl(Lx, Ly, Nx, Ny)
    printf("::Testing poisson operator::")
    printf("\n  Domain size:    [%0.2f, %0.2f]", Lx, Ly)
    printf("\n  Discretization: [%i, %i]", Nx, Ny)
    
    // X[i] = i*dx avec dx = Lx/Nx et i=0..Nx-1
    // Y[i] = j*dy avec dy = Ly/Ny et j=0..Ny-1
    X = linspace(0.0, Lx*(Nx-1)/Nx, Nx)
    Y = linspace(0.0, Ly*(Ny-1)/Ny, Ny)
    
    printf("\n\n  Initializing field W(x,y).")
    W   = feval(Y, X, init_field)
    
    printf("\n  Initializing reference solution Refx(x,y).")
    Refx = feval(Y, X, solution_fieldx)
    printf("\n  Initializing reference solution Refy(x,y).")
    Refy = feval(Y, X, solution_fieldy)

    dir  = get_absolute_file_path("test_poisson_curl.sce")
    file = dir+"poisson.sce" 
    printf("\n\n  Loading poisson_2d function from file %s%s%s.", char(39), file, char(39))
    exec(file, -1)

    printf("\n\n  Computing Poisson curl solution Ux(x,y) et Uy(x,y).")
    [Ux,Uy] = poisson_curl_2d(W, Nx, Ny, Lx, Ly)

 
    

    printf("\n  Computing error |Ux-Refx|(x,y).")
    Errx = abs(Ux-Refx)
    printf("\n  Computing error |Uy-Refy|(x,y).")
    Erry = abs(Uy-Refy)
    
    file = pwd()+"/poisson_curl_error.png"
    printf("\n\n  Plotting everything to %s%s%s.", char(39), file, char(39))
    plot_error(W, Refx, Refy, Ux,Uy)
    
    printf("\n\n")
    mErrx = max(Errx)
    mErry = max(Erry)
      
    printf("  Maximal error on Ux is %.10ef, Maximal error on Uy is %.10ef.\n", mErrx, mErry)
    exit(0)
endfunction


// Taille du domaine
Lx = 1.0
Ly = 1.0

// Discretisation du domaine
Nx = 64
Ny = 32


test_poisson_curl(Lx, Ly, Nx, Ny)

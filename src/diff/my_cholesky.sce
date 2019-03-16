
function [A]=cholesky_fact(A)
	n = size(A)(1,1)
	
	//Calcule la premiere colonne
	A(1,1) = sqrt(A(1,1))
	A(2:n,1) = A(2:n,1)/A(1,1)
	
	for p = 2:n
		A(1:(p-1),p) = 0
		//calcul le terme diagonal
		A(p,p) = sqrt(A(p,p) - A(p,1:(p-1))*A(p,1:(p-1))')
		//calcul le reste de la colonne
		A((p+1):n,p) = (A((p+1):n,p) - A((p+1):n, 1:(p-1))*A(p, 1:(p-1))')/A(p,p)
	end
endfunction

function [y]=up_sweep_cholesky(A,x)
  	[m,n]=size(A);
  	if (m~=n) then
   		print(%io(2), "error, not a square matrix");
  	else
  		y = zeros(n, 1)
		y(n) = x(n) / A(n,n)
		for i = -(-(n-1):-1)
			y(i) = (x(i) - A(i,(i+1):n)*y((i+1):n))/A(i,i)
		end
	end
endfunction

function [y]=down_sweep_cholesky(A,x)
  	[m,n]=size(A);
	if (m~=n) then
    	print(%io(2), "error, not a square matrix");
  	else
		y = zeros(n, 1)
		y(1) = x(1) / A(1,1)
		for i = 2:n
			y(i) = (x(i) - A(i,1:(i-1))*y(1:(i-1)))/A(i,i) 
		end
 	end
endfunction

function [U]=my_cholesky(N,S)
	T = cholesky_fact(N)
	y = down_sweep_cholesky(T, S)
	U = up_sweep_cholesky(T', y)
endfunction

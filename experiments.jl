using Metaheuristics
include("tools.jl")
include("testTrajectories.jl")

function funcError(p, ctr, punpres, distr)
    θ1 = 0
    D  = 6

    #Esta es la variable donde se almacenan las violaciones de las restricciones
    contx = 0

    if ctr[1] > 0
        contx += ctr[1] 
    end
    if ctr[2] > 0
        contx += ctr[2]
    end
    if ctr[3] > 0
        contx += ctr[3]
    end
    if ctr[4] > 0
        contx += ctr[4]
    end



    # Se ve que el individuo sea factible. Si es factible se calcula la funcion objetivo, si no lo es se coloca un valor
    # muy grande en el valor de dicha funcion. F.O=1000. Existe un caso especial
    # en el que el individuo es factible pero produce un mecanismo de doble
    # balancin, en ese caso se penaliza la funcion objetivo.

    N_points = size(punpres,1)

    flag = false
    if contx == 0
        my_error = 0

        X = zeros(N_points, 2)
        i = 1

        Q = floor(Int, p[11])

        θ2 = 0

        if Q != 0
        	distr2 = p[10]*(distr)
        else
        	distr2 = linspace(0, p[10], N_points) 
        end

        distr2 *= sign(p[end])
        
        for  θ2 = distr2

            A1= 2p[3] * (p[2] * cos(θ2)-p[1] * cos(θ1))
            B1= 2p[3] * (p[2] * sin(θ2)-p[1] * sin(θ1))
            C1= p[1]^2 + p[2]^2 + p[3]^2-p[4]^2 - 2p[1] * p[2] * cos(θ2-θ1)
            
            Diss = B1^2+A1^2-C1^2

            if Diss < 0
            	flag = true
				my_error = 1000
                break
    			# continue
            end

            θ3=2.0atan((-B1+sqrt(Diss))/(C1-A1))
            
            Cxr = p[2] * cos(θ2) + p[5] * cos(θ3)-p[6] * sin(θ3)
            Cyr = p[2] * sin(θ2) + p[5] * sin(θ3)+p[6] * cos(θ3)
            Cx = Cxr*cos(p[7])-Cyr*sin(p[7])+p[8]
            Cy = Cxr*sin(p[7])+Cyr*cos(p[7])+p[9]
            
            X[i, 1] = Cx
            X[i, 2] = Cy
            my_error+= norm(X[i,:] - punpres[i,:]) #(punpres[i, 1] - Cx)^2 + (punpres[i, 2]-Cy) ^ 2

            i += 1
        end

        # my_error += hausdorffDistance(X, punpres)


        

        
        sumr = contx
    else
        my_error=1000
        sumr = contx
    end

    return my_error, sumr

end


function myerror(x, error_func, precision_points, distr)
    f, g = error_func(x, grashof(x), precision_points, distr)

    return f + 1000g
end

function getDist(Pts,q=2)
	n = size(Pts, 1)

	ds = zeros(n)
	ds[1] = norm( Pts[1,:]- Pts[2,:],q)
	
	for i = 2:n-1
		ds[i] = ds[i-1]+ norm( Pts[i,:]- Pts[i+1,:],q)
	end

	ds /= maximum(ds)

	return ds
end

function distanceMatrix(X, Y)
    Nx = size(X,1)
    Ny = size(Y,1)

    D = zeros(Nx, Ny)

    for i = 1:Nx
        for j = 1:Ny
            D[i, j] = sum((X[i,:]- Y[j,:]) .^ 2)
        end
    end

    return D
end

function hausdorffDistance(X, Y)
    D = distanceMatrix(X, Y)

    h(D) = maximum( minimum(D, 1) )

    return max(h(D), h(D'))


end



function experimentsInfo(ncase, npts)
    bounds =[ 0 60;
          0 60;
          0 60;
          0 60;
        -60 60;
        -60 60;
          0 2π;
        -60 60;
        -60 60;
         0 2π;
         -10 10;
         -1 1;
        ]'
	D = size(bounds,2)

    precision_points = getTrayectory(ncase, npts)

     if ncase == 12
		bounds =[ 0 5;
		          0 5;
		          0 5;
		          0 5;
		        -5 5;
		        -5 5;
		          0 2π;
		        -5 5;
		        -5 5;
		         0 2π;
		         0 20;
		        ]'
     end

    
    error_func = funcError

    return D, bounds, precision_points, error_func, getDist(precision_points)
end

function findStructure(ncase, npts=50)
    # case information
    D, bounds, precision_points, error_func, distr = experimentsInfo(ncase, npts)

    # fitness
    fitness(x) = myerror(x, error_func, precision_points, distr)
    
    # optimize
    x, f = eca(fitness, D; saveLast = "last.csv",
                           saveConvergence = "conv.csv",
                           max_evals=20000D,
                           # N = 200,
                           # η_max=4,
                           limits=bounds)

    return x, f

end

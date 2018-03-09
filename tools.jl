# part of this code was taken from
# MATLAB code of Eric Valentin 
# CIDETEC 2018
# EMAIL: e.santiag.valentin[at]gmail.com

function grashof(p)
    gra = zeros(4)
    gra[1] = p[1]+p[2]-p[3]-p[4]
    gra[2] = p[2]-p[3]
    gra[3] = p[3]-p[4]
    gra[4] = p[4]-p[1]

    return gra

end

function secuencia(KC, n=6)
    # Esta función evalua la secuencia de los angulos de la manivela
    # del mecanismo de cuatro barras. Verifica que dicha secuencia sea ascendente.
    # El valor que regresa la funcion es cero si dicha secuencia es valida y
    # diferente de cero si no lo es. En este ultimo caso el valor representa la
    # suma de iolacion de restricciones de la secuencia. Un valor mayor indicca
    # que la secuencia es menos ordenada en ascendencia.

    # La entrada de la función es el conjunto de angulos de la manivela
    # y tiene el siguiente orden: [teta21, teta22, teta23, teta24, teta25, teta26]


    # Esta es la variable donde se almacenan las violaciones de las restricciones
    contx = 0;

    # Se busca el valor para el corrimiento de la secuencia y generar una copia
    # de dcho vector pero recorrida.
    offset = KC[1];

    #Se comienza a generar el nuevo vector.
    nvector = zeros(n)
    nvector[1]=0;

    # Se determina donde comienza la secuencia y se recorren los demas
    # elementos.
    if offset>=π  && offset<=2π 
        #offset=2π -offset
        for i=2:6
            if KC[i]>=π  && KC[i] <= 2π 
               nvector[i] = KC[i]-offset; 
            else
               nvector[i] = KC[i]+offset;
            end
        end
    else
        for i=2:6
            nvector[i] = KC[i]-offset;
        end
    end
              
    # Se verifica que se cumpla la ascendencia de la nueva secuencia.
    contx = 0
    for i=1 : n-1
        if nvector[i] - nvector[i+1]>=0
            contx += (nvector[i] - nvector[i+1]);
        end
    end

    # Se verifica que el ultimo valor no sobrepase al primer valor de la
    # secuencia original.
    if nvector[n] >= 2π 
        contx += 2π ;
    end
        
    return contx
end

# caso 2
function error_caso2(p, ctr, punpres)
    """
    Esta función evalua la funcion objetivo del problema de optimización
    del mecanismo de cuatro barras. Dicha funcion es el error cuadratico
    acumulado y es el valor que regresa esta rutina

    La entrada de la función es el conjunto de las 15 variables de diseño, las
    5 restricciones y los puntos de precision.
     y tiene el siguiente orden: p=[r1, r2, r3, r4, rcx, rcy θ_0, x0, y0, θ_21, θ_22, θ_23, θ_24, θ_25, θ_26]
     ctr=[4-grashof, 1-secuencia], punpres=[(c1x, c1y), (c2x, c2y),(c3x, c3y),(c4x, c4y),(c5x, c5y),(c6x, c6y)]
    """
    θ1 = 0
    D = 6
    pp = zeros(11)
    pp[1:6] = p
    pp[7:11] = [π/6, π/4, π/3, 10π/24, π/2]
    p = pp

    #Esta es la variable donde se almacenan las violaciones de las restricciones
    contx=0

    if ctr[1]>0
        contx += ctr[1] 
    end
    if ctr[2]>0
        contx += ctr[2]
    end
    if ctr[3]>0
        contx += ctr[3]
    end
    if ctr[4]>0
        contx += ctr[4]
    end


    # Se ve que el individuo sea factible. Si es factible se calcula la funcion objetivo, si no lo es se coloca un valor
    # muy grande en el valor de dicha funcion. F.O=1000. Existe un caso especial
    # en el que el individuo es factible pero produce un mecanismo de doble
    # balancin, en ese caso se penaliza la funcion objetivo.

    if contx == 0
        my_error = 0
        for i=1:5
            A1= 2p[3] * (p[2] * cos(p[D+i])-p[1] * cos(θ1))
            B1= 2p[3] * (p[2] * sin(p[D+i])-p[1] * sin(θ1))
            C1= p[1]^2 + p[2]^2 + p[3]^2-p[4]^2 - 2p[1] * p[2] * cos(p[D+i]-θ1)
            
            θ3=2 * atan((-B1+sqrt(B1^2+A1^2-C1^2))/(C1-A1))
            
            Cxr = p[2] * cos(p[D+i]) + p[5] * cos(θ3)-p[6] * sin(θ3)
            Cyr = p[2] * sin(p[D+i]) + p[5] * sin(θ3)+p[6] * cos(θ3)
            Cx  = Cxr
            Cy  = Cyr
            
            my_error+= (punpres[i, 1] - Cx)^2 + (punpres[i, 2]-Cy) ^ 2
        end
        
        # verifica si produce imaginarios la funcion objetivo.
        # Caso en el que se produce un doble balancin.
        if !isreal(my_error)  
            my_error=1000
        end
        
        sumr = contx
    else
        my_error=1000
        sumr = contx
    end

    return my_error, sumr

end

function error_caso1(p,ctr,punpres)
    """
    Esta función evalua la funcion objetivo del problema de optimización
    del mecanismo de cuatro barras. Dicha funcion es el error cuadratico
    acumulado y es el valor que regresa esta rutina

    La entrada de la función es el conjunto de las 15 variables de diseño, las
    5 restricciones y los puntos de precision.
     y tiene el siguiente orden: p=[r1, r2, r3, r4, rcx, rcy teta0, x0, y0, teta21, teta22, teta23, teta24, teta25, teta26]
     ctr=[4-grashof, 1-secuencia], punpres=[(c1x, c1y), (c2x, c2y),(c3x, c3y),(c4x, c4y),(c5x, c5y),(c6x, c6y)]
    """
    θ1=0;
    D=9;

    npoints = size(punpres, 1)

    tmp = zeros(5)
    tmp[1:4] = ctr
    tmp[5]   = secuencia(p[10:10+npoints-1])
    ctr = tmp

    #Esta es la variable donde se almacenan las violaciones de las restricciones
    # contx = sum( ctr[ ctr .> 0 ] )
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
    if ctr[5] > 0
        contx += ctr[5]
    end


    # Se ve que el individuo sea factible. Si es factible se calcula la funcion 
    # objetivo, si no lo es se coloca un valor
    # muy grande en el valor de dicha funcion. F.O=1000. Existe un caso especial
    # en el que el individuo es factible pero produce un mecanismo de doble
    # balancin, en ese caso se penaliza la funcion objetivo.

    if contx == 0
        my_error = 0
        for i=1:npoints
            A1 = 2p[3] * (p[2] * cos(p[D+i]) - p[1]*cos(θ1))
            B1 = 2p[3] * (p[2] * sin(p[D+i]) - p[1]*sin(θ1))
            C1 = p[1]^2+p[2]^2+p[3]^2-p[4]^2- 2p[1]*p[2]*cos(p[D+i]-θ1)
            
            θ3 = 2atan((-B1+sqrt(B1^2+A1^2-C1^2))/(C1-A1))
            
            Cxr = p[2]*cos(p[D+i])+p[5]*cos(θ3)-p[6]*sin(θ3)
            Cyr = p[2]*sin(p[D+i])+p[5]*sin(θ3)+p[6]*cos(θ3)
            
            Cx = Cxr*cos(p[7])-Cyr*sin(p[7])+p[8];
            Cy = Cxr*sin(p[7])+Cyr*cos(p[7])+p[9];
                      
            my_error+= (punpres[i, 1] - Cx)^2 + (punpres[i, 2]-Cy) ^ 2
        end
        
        # verifica si produce imaginarios la funcion objetivo.
        # Caso en el que se produce un doble balancin.
        if !isreal(my_error)
            my_error=1000
        end
        
        sumr = contx
    else
        my_error=1000
        sumr = contx
    end

    return my_error, sumr

end

function error_caso3(p,ctr,punpres1, punpres2)
    my_error1, sum1 = error_caso1(p, ctr, punpres1)
    my_error2, sum2 = error_caso1(p, ctr, punpres2)

    return my_error1 + my_error2, sum1 + sum2
end

function case_info(ncase)
    # bounds and precision points
    if ncase == 1
        BOUNDS =
        [ 0 60;
          0 60;
          0 60;
          0 60;
        -60 60;
        -60 60;
          0 2π;
        -60 60;
        -60 60;
          0 2π;
          0 2π;
          0 2π;
          0 2π;
          0 2π;
          0 2π
        ]'

        PRECISION_POINTS = 
        [
            20 20;
            20 25;
            20 30;
            20 35;
            20 40;
            20 45
        ]
        error_func = error_caso1
    elseif ncase == 2
        BOUNDS = 
        [   0 50;
            0 50;
            0 50;
            0 50;
          -50 50;
          -50 50
        ]'

        PRECISION_POINTS = 
        [       3 3;
            2.759 3.363;
            2.372 3.663;
            1.890 3.862;
            1.355 3.943
        ]

        error_func = error_caso2
    else
        BOUNDS = [
            0 60;
            0 60;
            0 60;
            0 60;
          -60 60;
          -60 60;
            0 2*pi;
          -60 60;
          -60 60;
            0 2*pi;
            0 2*pi;
            0 2*pi;
            0 2*pi;
            0 2*pi;
            0 2*pi;
            0 2*pi;
            0 2*pi;
            0 2*pi;
            0 2*pi
        ]'

        PRECISION_POINTS1 = [
            1.768 2.3311;
            1.947 2.6271;
            1.595 2.7951;
            1.019 2.7241;
            0.479 2.4281;
            0.126 2.0521;
           -0.001 1.720;
            0.103 1.514;
            0.442 1.549;
            1.055 1.905
        ]

        PRECISION_POINTS2 = [
           1.9592 2.44973;
            2.168 2.675;
            1.821 2.804;
            1.244 2.720;
            0.705 2.437;
            0.346 2.104;
            0.195 1.833;
            0.356 1.680;
            0.558 1.742;
            1.186 2.088
        ]

        f(p, ctr, punpres) = error_caso3(p,ctr, PRECISION_POINTS1, PRECISION_POINTS2)
        error_func = f
        PRECISION_POINTS = [PRECISION_POINTS1; PRECISION_POINTS2]

    end

    D = size(BOUNDS,2)
    return D, BOUNDS, PRECISION_POINTS, error_func
end

function generateTrayectory(p, nframes=25)
    # Save mechanism dynamics 
    r1,r2,r3,r4,rcx,rcy,θ0,x0,y0 = p
    θ1 = 0

    X0 = zeros(nframes, 2)
    X1 = zeros(nframes, 2)
    X2 = zeros(nframes, 2)
    X3 = zeros(nframes, 2)
    
    C = zeros(nframes, 2)
    i = 1
    for θ2 = linspace(0, 2π, nframes)
        
        

        A1 = 2p[3] * (p[2] * cos(θ2) - p[1]*cos(θ1))
        B1 = 2p[3] * (p[2] * sin(θ2) - p[1]*sin(θ1))
        C1 = p[1]^2+p[2]^2+p[3]^2-p[4]^2- 2p[1]*p[2]*cos(θ2 - θ1)

        θ3 = 2atan((-B1+sqrt(B1^2+A1^2-C1^2))/(C1-A1))


        Cxr = r2*cos(θ2)+p[5]*cos(θ3) - rcy*sin(θ3)
        Cyr = r2*sin(θ2)+p[5]*sin(θ3) + rcy*cos(θ3)

        C[i, 1] = Cxr*cos(θ0)-Cyr*sin(θ0) + x0;
        C[i, 2] = Cxr*sin(θ0)+Cyr*cos(θ0) + y0;
        
        X0[i, 1] = x0
        X0[i, 2] = y0

        X1[i, 1] = x0 + r1*cos(θ0)
        X1[i, 2] = y0 + r1*sin(θ0)

        X2[i, 1] = x0 + r2*cos(θ2+θ0)
        X2[i, 2] = y0 + r2*sin(θ2+θ0)

        X3[i, 1] = X2[i, 1] + r3*cos(θ3+θ0)
        X3[i, 2] = X2[i, 2] + r3*sin(θ3+θ0)
        
        i += 1
    end
    

    return C, X0, X1, X2, X3

end

function bestEver(ncase)
    if ncase == 1
        p = [39.35432733966176,8.738597335076454,27.91020567579063,38.41653206619049,36.576759454887835,
            18.058222622003274,3.922197401424182,-9.015069651004032,59.8594568458345,1.6275724377055465,
            2.42129538238649,2.93021895068062,3.4158306289051543,3.9401612344766854,5.428565504554411]
        my_error = 1e-20
        return p, my_error
    end
end
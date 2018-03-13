using DifferentialEquations


# Solve using Runge-Kutta (4th order)
function myode45(f::Function, tspan::Tuple, x0::Vector{Float64}, KC::Vector{Float64})
    
    # ODE function
    ddx(dx::Vector{Float64}, x::Vector{Float64}, p, t::Float64 ) = begin
        
        a = f(x,t, KC)
        
        for i = 1:length(x)
            dx[i] = a[i]
        end
        
        return dx
    end 

    # problem definition
    prob = ODEProblem(ddx, x0,tspan)

    # find solution
    sol = solve(prob, DP5())

    # returns time and solution matrix
    return sol.t, hcat(sol.u...)'
end


function mod4b(x::Vector{Float64},t::Float64, K::Vector{Float64})
    # Modulo para resolver el sistema de ecuaciones del modelo
    # dinamico del mecanismo de cuatro barras.
    # 12 de junio de 2012

    Kp, Ki, Kd = K

    # Parametros de la ecuacion diferencial del mecanismo de cuatro barras
    L1 = 0.5593
    L2 = 0.102
    L3 = 0.610
    L4 = 0.406
    r2 = 0.0
    r3 = 0.305
    r4 = 0.203
    m2 = 1.362
    m3 = 1.362
    m4 = 0.2041
    J2 = 0.00071
    J3 = 0.0173
    J4 = 0.00509
    ϕ2 = 0.0
    ϕ3 = 0.0
    ϕ4 = 0.0
    θ1 = 0.0
    
    # Constantes del resorte y amortiguador
    C  = 0.0
    k  = 0.0
    θ40= 0.0
    
    # Parametros de la ecuacion diferencial del motor de C.D.
    R = 0.4
    L = 0.05
    Kf= 0.678
    Kb= 0.678
    J = 0.056
    TL= 0.0
    B = 0.226
    
    # Parametros de la caja de engranes del motor de C.D.
    N1 = 1.0
    N2 = 1.0
    n = N2/N1

    θ2 = x[1]
    θ2p= x[2]

    # Ecuaciones del modelo cinematico del mecanismo de cuatro barras
    AA = 2.0L3 * (L2 * cos(θ2)-L1 * cos(θ1))
    BB = 2.0L3 * (L2 * sin(θ2)-L1 * sin(θ1))
    CC = L1^2+L2^2+L3^2-L4^2-2 * L1 * L2 * cos(θ1-θ2)
    θ3 = 2.0atan((-BB+(AA^2+BB^2-CC^2)^0.5)/(CC-AA))

    DD = 2.0L4 * (L1 * cos(θ1)-L2 * cos(θ2))
    EE = 2.0L4 * (L1 * sin(θ1)-L2 * sin(θ2))
    FF = L1^2+L2^2+L4^2-L3^2-2 * L1 * L2 * cos(θ1-θ2)
    θ4 = 2.0atan((-EE-(DD^2+EE^2-FF^2)^0.5)/(FF-DD))

    γ3 = (L2 * sin(θ4-θ2))/(L3 * sin(θ3-θ4))
    γ4 = (L2 * sin(θ3-θ2))/(L4 * sin(θ3-θ4))
    α2 = -r2 * sin(θ2-ϕ2)
    α3 = -L2 * sin(θ2)-r3 * γ3 * sin(θ3-ϕ3)
    α4 = -r4 * γ4 * sin(θ4-ϕ4)
    β2 =  r2 * cos(θ2+ϕ2)
    β3 =  L2 * cos(θ2)+r3 * γ3 * cos(θ3+ϕ3)
    β4 =  r4 * γ4 * cos(θ4+ϕ4)

    # Ecuaciones auxiliares del modelo dinámico
    C0  = J2+m2 * r2^2+m3 * L2^2
    C1  = J3+m3 * r3^2
    C2  = J4+m4 * r4^2
    C3  = 2.0m3 * L2 * r3
    Ax1 = C0+C1 * γ3^2+C2 * γ4^2+C3 * γ3 * cos(θ2-θ3-ϕ3)

    D1  = (γ4-1) * sin(θ3-θ4) * cos(θ4-θ2)
    D2  = (γ4-γ3) * sin(θ4-θ2) * cos(θ3-θ4)
    D3  = (γ3-1) * sin(θ3-θ4) * cos(θ3-θ2)
    D4  = (γ4-γ3) * sin(θ3-θ2) * cos(θ3-θ4)
    dg3 = (L2 * (D1+D2))/(L3 * (sin(θ3-θ4))^2)
    dg4 = (L2 * (D3+D4))/(L4 * (sin(θ3-θ4))^2)
    dAx1= 2.0C1 * γ3 * dg3 + 2.0C2 * γ4 * dg4+C3 * (dg3 * cos(θ2-θ3-ϕ3)-γ3 * (1-γ3) * sin(θ2-θ3-ϕ3))


    A0 = 1.0/(Ax1+n^2 * J)
    A1 = -(0.5) * dAx1
    A2 = -(C * γ4^2+n^2 * B)
    A3 = -(n * TL)-k * γ4 * (θ4-θ40)

    # Ley de control PID
    xref = 30.0
    myerror = xref-θ2p
    x2_p = A0 * (A1 * θ2p^2+A2 * θ2p+n * Kf * x[3]+A3)

    #u = 30
    u = Kp * myerror * x[5] + Ki * x[4] + Kd * (-x2_p)
    x3_p = (1/L) * (u-n * Kb * θ2p-R * x[3])
    x6_p = ((1/n * Kf) * ((0.5) * dAx1 * x2_p^2+n^2 * B * x2_p^2+n * TL))^2

    #Ecuaciones de estado
    #Las primeras 3 son las ecuaciones de estado, la siguiente  1 es la
    #integral del myerror,la siguiente 1 es la integral de la velocidad deseada,
    #la siguiente 1 es la integral de la corriente funcion objetivo.
    return [θ2p, x2_p, x3_p, myerror, xref, x6_p]
end

function merito1(KC::Vector{Float64}, onlydelta::Bool=false)
    #Esta función evalua la funcion de variacion de la velocidad de entrada
    #al mecanismo de cuatro barras.
    #La entrada de la función es el individuo a evaluar, donde el vector
    #individuo tiene el siguiente orden: [Kp,Ki,Kd]
    ########### ESTA FUNCION SE MINIMIZA ############


    #Esta es la variable donde se almacenan la violaciones de las restricciones
    contx = 0.0
    #evaluacion de la aptitud de los individuos
    #tiempo de simulacion y paro del algoritmo.
    t0 = 0.0
    tf = 4π/30
    TSPAN=(t0, tf)
    
    # Condiciones iniciales de las variables de estado del sistema    
    θ2  = 0.0
    θ2p = 0.0
    x30 = 0.0 # DC motor
    x40 = 0.0 # ∫error
    x50 = 0.0 # ∫v (desired)
    x60 = 0.0 # ∫I^2

    # initial conditions
    x0 = [θ2, θ2p, x30, x40, x50, x60]
      
    # Solucion del sistema dinamico
    t, x = myode45(mod4b,TSPAN,x0,KC) 

    #Código para evaluar la restriccion del rise time. 
    #El rise time se asigna en trise.

    trise = 0.1

    i = 1
    ren, col = size(x)

    #evaluacion de la primera restriccion
    #el rise time maximo es de 0.1
    while x[i, 2] < 30.0 && i < ren
      i += 1
    end
    if i == ren
      contx = 0.3
    else
      if t[i] > trise
         contx = t[i]-trise
      end
    end

    #evaluacion de la segunda restriccion, el overshoot menor del 1.7#

    if i==ren
        ov = 100.0
        contx += ov
    else
        ov = ((maximum(x[:,2]) - 30.0)/30.0)*100.0
        if ov >= 1.7
            contx += ov
        end
    end


    # Codigo para evaluar la función de mérito.
    # Si es un individuo factible, evalua la función,
    # Si no es factible le asigna el valor 1000 como valor de merito

    if (contx == 0.0)
       Δx = abs(maximum(x[:,2])-minimum(x[i:ren,2]))

       fit = Δx
       viorest = contx

    else
       fit = 1000.0
       viorest = contx
    end

    if onlydelta
        return abs(maximum(x[:,2])-minimum(x[i:ren,2]))
    end

    return fit, viorest
    #Eso es todo...

end
function soln = realSpaceSolver()
    NSteps = 20000;

    veloSpread = .05;
    posSpread = 1/(2*veloSpread);

    tau = .2596;

    n = 4;
    int = 3.0871;
    fBragg = 4*n;
    wm = 8*n;

    ti = 0;
    tf = 88*tau;
    t = linspace(ti,tf,NSteps);
    t1 = 10*tau;
    t2 = 32*tau;
    t3 = 42*tau;
    t4 = 64*tau;

    xi = -12*posSpread;
    xf = 70*posSpread;
    x = linspace(xi,xf, NSteps);

    xInit = 15*posSpread;

    m = 0;

    tic
    sol = pdepe(m,@pdeParams,@initialState,@pdeBCs,x,t);
    toc
    % Extract the first solution component as u.
    u = sqrt(abs(sol(:,:,1)).^2);

    soln = abs(sol(:,:,1)).^2;
    
    % surf(x,t,u) 
    figure
    pcolor(x,t,u)
    hold on
    shading interp 
    colormap jet
    xlim([xi,xf])
    ylim([ti,tf])
    title('SCI Interferometer')
    xlabel('Distance x')
    ylabel('Time t')
    hold off

    % A solution profile can also be illuminating.
    figure
    plot(x,u(end,:).^2)
    hold on
    title('Solution at t = tf')
    xlabel('Distance x')
    ylabel('u(x,tf)')
    hold off
    % --------------------------------------------------------------
    function [c,f,s] = pdeParams(x,t,u,DuDx)
        c = 1;
        f = -1i*DuDx;
        s = (-1i*(2*pi*int)).*(exp(-(t-t1).^2./(tau.^2)).*cos(fBragg.*t - 2.*x)+exp(-(t-t2).^2./(tau.^2)).*cos(fBragg.*t - 2.*x)+2.*exp(-(t-t3).^2./(tau.^2)).*cos(fBragg.*t - 2.*x).*cos(wm.*t)+2.*exp(-(t-t4).^2./(tau.^2)).*cos(fBragg.*t - 2.*x).*cos(wm.*t)).*u;
    end
    % --------------------------------------------------------------
    function u0 = initialState(x)
        u0 = sqrt(pi./posSpread).*exp(-(1/2).*((x-xInit)./posSpread).^2);
    end
    % --------------------------------------------------------------
    function [pl,ql,pr,qr] = pdeBCs(xl,ul,xr,ur,t)
        pl = 0;
        ql = 1;
        pr = 0;
        qr = 1;
    end
    
end
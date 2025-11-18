function h = plot_stars()

[example, beta_choice, alpha, x0, y0, a,b,...
    eps_disk,m_polar,x0_polar,y0_polar,r0_polar,r1_polar,...
    n_polar, bc] = read_vars();

th = linspace(0,2*pi,512);
for i = 1:m_polar
    xc = x0_polar(i);
    yc = y0_polar(i);
        
    r0 = r0_polar(i);
    r1 = r1_polar(i);
    n = n_polar(i);
    
    x = r0*(1 + r1*cos(n*th)).*cos(th) + xc;
    y = r0*(1 + r1*cos(n*th)).*sin(th) + yc;
    
    h(i) = plot(x,y,'k','linewidth',2);
    hold on;        
end

end
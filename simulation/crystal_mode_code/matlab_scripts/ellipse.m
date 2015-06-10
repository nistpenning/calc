function [x,y] = ellipse(a,b,x0,y0,col,edge,flag)
  % (x0,y0) - center
  % (a,b) - ellipse axes
  t = linspace(0,2*pi,1000);
  x = x0 + a*cos(t);
  y = y0 + b*sin(t);
  if flag==1
    patch(x,y,col,'EdgeColor',edge)
  end
  if flag == 2
    plot(x,y,col,'Color',edge)
  end
end
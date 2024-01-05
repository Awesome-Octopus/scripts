[x, y params] = eprload;
y0thorder = y-(mean(y(1:150))); % subtract a flat baseline from the linear part of the spectrum;
x = x(280:end-350);
y0thorder = y0thorder(280:end-350);
plot (x, y0thorder);
intspec = cumtrapz(x, y0thorder);
hold on;
plot(x, intspec);
deltay = intspec(end)-intspec(1);
q = numel(intspec);
for n=1:q
    linCorrection(n) = intspec(n)-(deltay*n/q);
end
plot (x,linCorrection);
totalArea = cumtrapz(x,linCorrection);
plot(x, totalArea);
txt = sprintf("the background corrected integrated area of your spectrum is: %f\n", totalArea(end));
disp(txt);
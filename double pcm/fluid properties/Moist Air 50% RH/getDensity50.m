function density = getDensity50(temp)
density = 3.899e-12*temp.^4 - 8.854e-9*temp.^3 + 7.946e-6*temp.^2 -3.828e-3*temp + 1.238;
end
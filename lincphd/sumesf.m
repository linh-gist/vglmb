function sums = sumesf(numbers)

m = length(numbers);
if m==0
    sums = 1;
else
    result = poly(numbers);
    sums = abs(result)';
end
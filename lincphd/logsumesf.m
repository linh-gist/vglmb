function logsums = logsumesf(lognumbers)

%Calculates log transformed values for elementary symmetric functions
%Inputs are log-transformed values
%Expands with an iterated convolution performed entirely in log-domain
%with logconv which can be computationally expensive since the logconv
%function does not use an fft shortcut to perform convolution
%WARNING: ensure that input numbers are strictly positive real
%so that log(numbers) has imaginary parts that are strictly multiples of pi
%thus the calculation can ignore coef signs in intermediate calculations
%since the end result is also strictly positive real

m = length(lognumbers);
if m==0
    logsums = 0;
else
    
    result = [0; lognumbers(1)];
    for idx=2:m
        result = logconv([0; lognumbers(idx)],result);
    end
    
    logsums = result;
end
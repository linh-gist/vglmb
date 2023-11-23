function logsums = logsumesfsw(numbers)

%Calculates log transformed values for elementary symmetric functions
%Inputs are untransformed values
%When an overflow is suspected, algorithm algorithm divides the root set 
%into a number smaller sets and expands the polynomials for each first,
%then calls the logconv operation on their log-transforms to get the result
%i.e. this minimizes calls to the logconv operation which is very expensive
%WARNING: ensure that input numbers are strictly positive real
%so that log(numbers) has imaginary parts that are strictly multiples of pi
%thus the calculation can ignore coef signs in intermediate calculations
%since the end result is also strictly positive real

cutlength = 100;
m = length(numbers);
if m==0
    logsums = 0;
elseif m<=cutlength
    logsums = log(abs(poly(numbers))+eps(0))';
else
    nparts = floor(m/cutlength);
    logbits = zeros(cutlength+1,nparts);
    for cp=1:nparts
        coefs = numbers((cp-1)*cutlength+1:cp*cutlength);
        logbits(:,cp) = log(abs(poly(coefs))+eps(0))';
    end
    remcoefs = numbers(nparts*cutlength+1:end);
    logremdr = log(abs(poly(remcoefs))+eps(0))';

    result = logbits(:,1);
    for cp=2:nparts
        result = logconv(logbits(:,cp),result);
    end
    
    logsums = logconv(logremdr,result);
end
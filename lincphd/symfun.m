function combs = symfun(numbers)

m = length(numbers);
combs = cell(m,1);

%i=1
combs{1} = numbers;
highidx_old = numbers;
disp([num2str(1),' ',num2str(sum(numbers))]);

%i=2:m
for i=2:m
    old = combs{i-1};
    new = zeros(nchoosek(m,i),1);
      
    new_idx = 1;
    for old_idx=1:length(old)
            for pt = highidx_old(old_idx)+1:m
                new(new_idx) = old(old_idx)*numbers(pt);
                highidx_new(new_idx) = pt;
                new_idx = new_idx+1;
            end
    end
    highidx_old = highidx_new;
    clear highidx_new;
combs{i} = new;
disp([num2str(i),' ',num2str(sum(new))]);
end
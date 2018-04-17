function t = firstDigit(typ)
    t = typ;
    for i = 1:length(typ)
        while(typ(i) > 0)
            t(i) = typ(i);
            typ(i) = fix(typ(i)/10);
        end
    end
end
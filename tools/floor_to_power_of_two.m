function result = floor_to_power_of_two(number)
    if number <= 0
        result = 0;
    else
        result = 2^floor(log2(number));
    end
end
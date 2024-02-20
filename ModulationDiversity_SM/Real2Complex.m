function Complex_X = Real2Complex(Real_X)
    Nt = size(Real_X, 1)/2;
    totalTime = size(Real_X, 2);
    Complex_X = zeros(Nt, totalTime, size(Real_X, 3));

    for timeSlot = 1:totalTime
        for ant=1:Nt
            Complex_X(ant, timeSlot, :) = Real_X(ant, timeSlot, :)+1j*Real_X(ant+Nt, timeSlot, :);
        end
    end
end
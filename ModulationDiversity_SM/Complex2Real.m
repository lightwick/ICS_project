function Real_X = Complex2Real(Complex_X)
    Real_X = [real(Complex_X);
        imag(Complex_X)];
end
function retVal = transistorComparison( DTIRPath )

    retVal.d = readtable(DTIRPath);

    figure();
    plot( retVal.d.time, retVal.d.n7 );
    title( "S-Param input" );
    figure();
    plot( retVal.d.time, retVal.d.n16 );
    title( "S-Param output" );

    
    f_s.d = 1 / retVal.d.time( 2 );
    f_s.d = f_s.d;

    time.d = retVal.d.time(end) + retVal.d.time( 2 );
    time.d = time.d;


    %freqSequence.d = 0 : 1 / time.d : f_s.d - 1 / time.d;
    %freqSequence.v = 0 : 1 / time.v : f_s.v - 1 / time.v;
    
    freqSequence.d = - f_s.d/2 : 1 / time.d : f_s.d/2 - 1 / time.d;

    len = length(retVal.d.n7);
    N7.d = fftshift( fft( retVal.d.n7 ) ) / len;
    N16.d = fftshift(fft( retVal.d.n16 )) / len;
    
    figure();
    plot( freqSequence.d, abs( fftshift( fft( retVal.d.n9 ) ) ) / len );
    title( "FFT Transistor input" );
    set( gca, 'YScale', 'log' );
    
    figure();
    plot( freqSequence.d, abs( N7.d ) );
    title( "FFT S-Param input DTIR" );
    set( gca, 'YScale', 'log' );
    
    figure();
    plot( freqSequence.d, abs( N16.d ) );
    title( "FFT S-Param output DTIR" );
    set( gca, 'YScale', 'log' );
    
    figure();
    plot( freqSequence.d, abs( N16.d ) ./ abs( N7.d ) );
    title( "FFT S-Param Trans DTIR" );
    set( gca, 'YScale', 'log' );
    disp("done");

end
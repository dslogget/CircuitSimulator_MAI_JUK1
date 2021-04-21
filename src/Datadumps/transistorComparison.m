function retVal = transistorComparison( DTIRPath, VFPath )

    retVal.d = readtable(DTIRPath);
    retVal.v = readtable(VFPath);
    sametime = false;
    if (length( retVal.d.time ) == length( retVal.v.time )) && all(retVal.d.time == retVal.v.time)
        sametime = true;
    end

    figure();
    plot( retVal.d.time, retVal.d.n7 ); hold on;
    plot( retVal.v.time, retVal.v.n7, '--' ); hold off;
    title( "S-Param input" );
    figure();
    plot( retVal.d.time, retVal.d.n16 ); hold on;
    plot( retVal.v.time, retVal.v.n16, '--' ); hold off;
    title( "S-Param output" );
    if sametime
        figure();
        plot( retVal.d.time, retVal.d.n7 - retVal.v.n7 );
        title( "S-Param input difference" );
        figure();
        plot( retVal.d.time, retVal.d.n16 - retVal.v.n16 );
        title( "S-Param output difference" );
    end
    
    f_s.d = 1 / retVal.d.time( 2 );
    f_s.v = 1 / retVal.v.time( 2 );
    f_s.d = f_s.d;
    f_s.v = f_s.v;

    time.d = retVal.d.time(end) + retVal.d.time( 2 );
    time.v = retVal.v.time(end) + retVal.v.time( 2 );
    time.d = time.d;
    time.v = time.v;

    %freqSequence.d = 0 : 1 / time.d : f_s.d - 1 / time.d;
    %freqSequence.v = 0 : 1 / time.v : f_s.v - 1 / time.v;
    
    freqSequence.d = - f_s.d/2 : 1 / time.d : f_s.d/2 - 1 / time.d;
    freqSequence.v = - f_s.v/2 : 1 / time.v : f_s.v/2 - 1 / time.v;

    len = length(retVal.d.n7);
    N7.d = fftshift( fft( retVal.d.n7 ) ) / len;
    N16.d = fftshift(fft( retVal.d.n16 )) / len;
    N7.v = fftshift(fft( retVal.v.n7 )) / len;
    N16.v = fftshift(fft( retVal.v.n16 )) / len;
    
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
    plot( freqSequence.v, abs( N7.v ) );
    title( "FFT S-Param input VF" );
    set( gca, 'YScale', 'log' );
    
    figure();
    plot( freqSequence.v, abs( N16.v ) );
    title( "FFT S-Param output VF" );
    set( gca, 'YScale', 'log' );
    disp("done");
    
    figure();
    plot( freqSequence.d, abs( N16.d ) ./ abs( N7.d ) );
    title( "FFT S-Param Trans DTIR" );
    set( gca, 'YScale', 'log' );
    disp("done");
    
    figure();
    plot( freqSequence.v, abs( N16.v ) ./ abs( N7.v ) );
    title( "FFT S-Param Trans VF" );
    set( gca, 'YScale', 'log' );
    disp("done");
    
    

end
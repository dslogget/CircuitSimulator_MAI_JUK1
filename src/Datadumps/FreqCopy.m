t1 = readtable("../datadump.txt");
%% sect
offset = 1500100;
step = 0.001;
periodNS = 1;
period = periodNS / step;
nPeriods = 100;
fLower = 0;
fUpper = 200;

%time = t1.time(end) + t1.time(2) - offset*t1.time(2);
f_s = 1/t1.time(2);
freqSequence = - f_s/2 : 1 / (nPeriods*periodNS) : f_s/2 - 1 / (nPeriods*periodNS);

f1 = abs( fftshift( fft( t1.n7((offset):(offset + nPeriods*period - 1))) ) );
f1 = f1 / length( f1 );
p1 = ( 2 * abs( f1 ).^2 ) / (2*50);

%f2 = abs( fftshift( fft( t1.n16((offset):(offset + nPeriods*period - 1))) ) );
%p2 = ( f2 .* conj(f2) ) / length( f2 );
%f2 = f2 ./ length(f2);

fi = abs( fftshift( fft( t1.n9((offset):(offset + nPeriods*period - 1))) ) );
fi = fi ./ length(fi);
pin = ( 2 * abs( fi ).^2 ) / (2*50);

dbm1 = 10 * log10( p1 / (1E-3));
%dbm2 = 10 * log10( p2 / (1E-3));
dbmi = 10 * log10( pin / (1E-3));

[~, idxLower] = min( abs( freqSequence - fLower ) );
[~, idxUpper] = min( abs( freqSequence - fUpper ) );

figure();
hold on;
base = min([dbmi;dbm1]) - 10;
stem( freqSequence, dbmi, 'BaseValue', base );
stem( freqSequence, dbm1, 'BaseValue', base );
%stem( freqSequence, dbm2 );
%set( gca, 'YScale', 'log' );
hold off;
dlmwrite("fftoutput.txt", [freqSequence( idxLower:idxUpper ).', dbm1( idxLower:idxUpper ), dbmi( idxLower:idxUpper ) ],'delimiter',' ','precision','%20.12f');
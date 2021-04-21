function  [S, z_ref, numPorts] = CPPVectFitAdaptor( filePath, N, Niter )  
    if nargin < 3 || isempty(Niter)
        Niter = 5;
        if nargin < 2 || isempty(N)
            N = 100;
        end
    end
    s2p = sparameters( filePath );

    freqSamples = s2p.Frequencies / 1e9;
    Ns = length( freqSamples );
    z_ref = s2p.Impedance;
    numPorts = s2p.NumPorts;
    
    freqSamplesLong = [ freqSamples; (freqSamples(end) + freqSamples(2:end-1))];
    freqSamplesLong = [ freqSamplesLong; (freqSamplesLong(end) + freqSamplesLong(2:end-1))];
    freqSamplesLong = [ freqSamplesLong; (freqSamplesLong(end) + freqSamplesLong(2:end-1))];
    freqSamplesLong = [ freqSamplesLong; (freqSamplesLong(end) + freqSamplesLong(2:end-1))];
    freqSamplesLong = [ freqSamplesLong; (freqSamplesLong(end) + freqSamplesLong(2:end-1)); 100000];
    
    
    fitData = struct("poles", {}, "residues", {}, "remainder", {});
    S = repmat( fitData, s2p.NumPorts, s2p.NumPorts );
    

    %Complex starting poles :
    bet=linspace(2*pi*freqSamples(1),2*pi*freqSamples(Ns),N/2);

    
    opts.relax=0;      %Use vector fitting with relaxed non-triviality constraint
    opts.stable=1;     %Enforce stable poles
    opts.asymp=2;      %Include D in fitting    
    opts.skip_pole=0;  %Do NOT skip pole identification
    opts.skip_res=0;   %DO NOT skip identification of residues (C,D,E) 
    opts.cmplx_ss=1;   %Create real-only state space model

    opts.spy1=0;       %No plotting for first stage of vector fitting
    opts.spy2=0;       %Create magnitude plot for fitting of f(s) 
    opts.logx=0;       %Use linear abscissa axis
    opts.logy=0;       %Use logarithmic ordinate axis 
    opts.errplot=1;    %Include deviation in magnitude plot
    opts.phaseplot=0;  %Do NOT produce plot of phase angle
    opts.legend=1;     %Include legends in plots
    
    
    
    % Find a complex rational model for freq domain data
    for a = 1 : 1 : s2p.NumPorts
        for b = 1 : 1 : s2p.NumPorts
            poles=[];
            for n=1:length(bet)
              alf=-bet(n)*1e-2;
              poles = [poles (alf-1j*bet(n)) (alf+1j*bet(n)) ]; 
            end
            
            posData = squeeze(s2p.Parameters( a, b, : ));
            t3 = 2*pi*1j*freqSamples.';%[freqSamples; (freqSamples(end) + freqSamples(2:end-1) ) ];
            t4 = posData.';%[posData; flip(conj(posData(2:end-1))) ];
            
            
            for iter=1:Niter
                [SER,poles,rmserr,fit,opts] = vectfit3( t4, t3, poles, ones(1,Ns), opts );
            end
            S(a, b).poles = complex(poles);
            S(a, b).residues = complex(SER.C);
            S(a, b).remainder = complex(SER.D);
            disp(rmserr);
            
            Dk=zeros(length(freqSamplesLong),N); 
            for m=1:N
                Dk(:,m)=1./(2*pi*1j*freqSamplesLong-S(a, b).poles(m));
            end 
            fitlonger = (Dk*S(a, b).residues.').' + S(a, b).remainder;
            
            
            
            if ( max( abs(fitlonger) ) > 1 )
               close all;
               [S, z_ref, numPorts] = CPPVectFitAdaptor( filePath, N-2 );
               return;
            end
            if ( max( abs( S(a, b).residues ) ) > 1000 )
                disp( "Large residue warning" );
            end
            
            figure();
            plot(freqSamples, abs(t4)); hold on; plot(freqSamplesLong, abs(fitlonger) ); hold off;
            
            filepath = sprintf("S(%d,%d).txt",a,b);
            dlmwrite(filepath, [freqSamplesLong(1:floor(length(freqSamplesLong)/2)), abs(fitlonger(1:floor(length(freqSamplesLong)/2)).')],'delimiter',' ','precision','%20.12f') 
        end
    end
    fprintf( "%d poles", N );
    S = S.';
end
function [d, v] = comp( DTIRPath, VFPath, nodeNum )

    d = readtable(DTIRPath);
    v = readtable(VFPath);
    sametime = false;
    if d.time == v.time
        sametime = true;
    end

    figure();
    plot( d.time, d.(nodeNum) ); hold on;
    plot( v.time, v.(nodeNum), '--' ); hold off;
    if sametime
        figure();
        plot( d.time, d.(nodeNum) - v.(nodeNum) );
    end

end
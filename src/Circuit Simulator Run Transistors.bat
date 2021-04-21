@echo off
CircuitSimulator.exe "Netlists/FClass_Amplifier10G_VF.netlist" > "Datadumps/FClass_Amplifier10G_VF.out"
timeout 5
CircuitSimulator.exe "Netlists/FClass_Amplifier10G_DTIR.netlist" > "Datadumps/FClass_Amplifier10G_DTIR.out"
timeout 5
CircuitSimulator.exe "Netlists/FClass_Amplifier10G_DTIR_Trunc.netlist" > "Datadumps/FClass_Amplifier10G_DTIR_Trunc.out"
pause
@echo off
CircuitSimulator.exe "Netlists\Benchmarks\NUDTIR\1GNUDTIR10.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\NUDTIR\1GNUDTIR50.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\NUDTIR\1GNUDTIR100.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\NUDTIR\1GNUDTIR500.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\NUDTIR\1GNUDTIR1000.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\NUDTIR\1GNUDTIR5000.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\NUDTIR\1GNUDTIR10000.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\NUDTIR\1GNUDTIR50000.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\NUDTIR\1GNUDTIR100000.netlist" > nul

timeout 5

CircuitSimulator.exe "Netlists\Benchmarks\VFRC\1GVFRC10.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\VFRC\1GVFRC50.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\VFRC\1GVFRC100.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\VFRC\1GVFRC500.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\VFRC\1GVFRC1000.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\VFRC\1GVFRC5000.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\VFRC\1GVFRC10000.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\VFRC\1GVFRC50000.netlist" > nul
timeout 5
CircuitSimulator.exe "Netlists\Benchmarks\VFRC\1GVFRC100000.netlist" > nul

pause

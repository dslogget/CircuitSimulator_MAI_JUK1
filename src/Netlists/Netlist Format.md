# Netlist format
Most netlist components follow the format:
```
<Component  Prefix><ID> <n1> <n2> ... <nN> <value1> <value2> ... <valueN>
```
# Core Components
## Voltage Sources
### DC source
```
V<id> <n1> <n2> <Voltage (V)>
```
### Sine Source
```
VS<id> <n1> <n2> ( <Amplitude (V)=1> <Frequency (GHz)=1> <DC offset (V)=0> <Phase offset (degrees)=0> )
```
### Time series source
```
VT<id> <n1> <n2> <Time Multiplier> <File path>
```
## Current Sources
### DC current source
```
I<id> <n1> <n2> <Current (A)>
```
### NLCurrentSource
```
IN<id> <n1> <n2> <r1+> <r1-> <r2+> <r2->
```
represents the equation described by
```
constexpr T alpha = 1.3;
constexpr T beta0 = 0.42;
constexpr T gamma = 0.0005;
constexpr T delta = 0.3;
constexpr T xi = 0.06;
constexpr T lambda = 1.5;
constexpr T mu = 0.0;
constexpr T zeta = 0.18;
constexpr T Vto = -2.4;

using ADT = AD::DiffVar<T, 2>;
ADT V_gs(r1, 1, 0);
ADT V_ds(r2, 0, 1);

auto beta = beta0;
auto Vgst = V_gs - (1 + beta * beta) * Vto + gamma * V_ds;
auto Veff = 0.5 * (Vgst + AD::sqrt(AD::pow(Vgst, 2) + delta * delta));
auto power = lambda / (1 + mu * AD::pow(V_ds, 2) + xi * Veff);
auto area = alpha * V_ds * (1 + zeta * Veff);
auto f1 = AD::tanh(area);
auto Ids_lim = beta * AD::pow(Veff, power);
auto Idrain = Ids_lim * f1;
auto I_ds = Idrain[0] - Idrain[1] * r1 - Idrain[2] * r2;
```

## Resistor
```
R<id> <n1> <n2> <Resistance (Ohms)> (<Any character to indicate group 2>)
```
## Capacitors
### Linear Capacitor
```
C<id> <n1> <n2> <Capacitance (nF)>
```
### Non-linear Capacitor
```
CN<id> <n1> <n2> <C_p (nF)> <C_o (nF)> <P_10> <P_11>
```
Where the capacitance follows the following equation:
```
C = C_p + C_o * (1.0 + std::tanh(P_10 + P_11 * u));
```
## Inductor
```
L<id> <n1> <n2> <Inductance (nH)>
```
## S-parameter Blocks
### NUDTIR block
```
S<id> <Pruning threshold % of max val (recommend 0)> <number of ports> <p1+> <p1-> ... <pn+> <pn-> <S-param file path>
```
### VF block (PRR)
```
SV<id> <number of ports> <p1+> <p1-> ... <pn+> <pn-> <PRR file path>
```
### VF block (Matlab)
```
SVF<id> <number of ports> <p1+> <p1-> ... <pn+> <pn-> <S-param file path>
```

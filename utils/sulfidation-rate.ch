set terminal postscript enhanced color eps "Times, 18" size 8cm,8cm
system "mkdir -p output"
set output 'output/sulfidation-rate.eps'

T=500 # in K
R=8.314 # Gas constant

sulfid_gas_met(x)= 10**((0.00473*T)-5.645+((x+63562)/(R*T)))  # in nm/s
sulfid_gas_mkw(x)= 0.01372 + 0.04356*(exp(x/(R*T)))  # in nm/s
sulfid_gas_pht(x)= exp(-(11766/T)-0.6478)*1E9 # in nm/s
sulfid_gas_pyr(x)= 7.45E8 * exp(-(98400/(R*T))) # in nm/s
sulfid_liq_met(x)= 0.0666*(1 + ((x+63562)/(R*T)))  # in nm/s
sulfid_liq_mkw(x)= 0.1332*(1 + (2*(x+63562)/(R*T)))  # in nm/s
sulfid_liq_pht(x)= 2.416  # in nm/s
sulfid_liq_pyr(x)= 0.003543 # in nm/s

set xlabel '{/Symbol m}_S (kJ/mol)'
set ylabel 'Log (Growth rate/ms^{-1})'

set xrange [-75000:25000]

set xtics('-75' -75000,'-50' -50000,'-25' -25000,'0' 0, '25' 25000)
set xtics out nomirror
set ytics 5 out nomirror

set key top left

plot log10(sulfid_gas_met(x)*1E-9) w l lt 1 lw 3 lc rgb 'red' title 'Met-Gas', \
     log10(sulfid_gas_mkw(x)*1E-9) w l lt 1 lw 3 lc rgb 'magenta' title 'Mkw-Gas', \
     log10(sulfid_gas_pht(x)*1E-9) w l lt 1 lw 3 lc rgb 'blue' title 'Pht-Gas', \
     log10(sulfid_gas_pyr(x)*1E-9) w l lt 1 lw 3 lc rgb 'black' title 'Pyr-Gas', \
     log10(sulfid_liq_met(x)*1E-9) w l lt 2 lw 3 lc rgb 'red' title 'Met-Liq', \
     log10(sulfid_liq_mkw(x)*1E-9) w l lt 2 lw 3 lc rgb 'magenta' title 'Mkw-Liq', \
     log10(sulfid_liq_pht(x)*1E-9) w l lt 2 lw 3 lc rgb 'blue' title 'Pht-Liq', \
     log10(sulfid_liq_pyr(x)*1E-9) w l lt 2 lw 3 lc rgb 'black' title 'Pyr-Liq'


quit
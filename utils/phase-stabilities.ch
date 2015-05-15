set terminal postscript enhanced color eps "Times, 18" size 8cm,8cm
system "mkdir -p output"
set output 'output/phase-stabilities.eps'

set border 3
set xtics 25 out nomirror
set ytics 25 out nomirror

set xlabel 'Sulfur chemical potential, {/Symbol m}_S, (kJ/mol)'
set ylabel 'Free energy, G, (kJ/mol)'

set xrange [-75:25]

set xtics 25 out nomirror
set ytics 50 out nomirror

plot 'free-energies.txt' u ($1/1000):($2/1000) w l lt 1 lw 3 lc rgb 'blue' title 'Metal', \
     'free-energies.txt' u ($1/1000):($3/1000) w l lt 1 lw 3 lc rgb 'red' title 'Mackinawite', \
     'free-energies.txt' u ($1/1000):($4/1000) w l lt 1 lw 3 lc rgb 'black' title 'Pyrrhotite', \
     'free-energies.txt' u ($1/1000):($5/1000) w l lt 1 lw 3 lc rgb 'magenta' title 'Pyrite', \
     'free-energies.txt' u ($1/1000):($6/1000) w l lt 1 lw 3 lc rgb 'dark-green' title 'Environment'


     
quit

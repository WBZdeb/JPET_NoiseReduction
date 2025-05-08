# Ustawienia wyjścia
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'activity_plot_norm.png'

# Tytuły i etykiety
set title 'Zależność średniej liczby randomów od aktywności'
set xlabel 'Aktywność'
set ylabel 'Średnia liczba randomów na okno'

# Skala logarytmiczna osi X (opcjonalnie, jeśli duży rozrzut danych)
#set logscale x

# Dane i dopasowanie funkcji kwadratowej
f(x) = a*(x/1e6)**2 + b*(x/1e6) + c
fit f(x) 'act_output.txt' using 1:2 via a, b, c

# Wykres
plot 'act_output.txt' using 1:2 with points pointtype 7 pointsize 1.5 title 'Dane', \
     f(x) with lines lw 2 lc rgb 'red' title 'Dopasowanie: f(x) = a*x² + b*x + c'

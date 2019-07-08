set terminal postscript enhanced color
#set terminal pdf
set logscale xy
set style line 1 lw 3

set xlabel "Redshift"
set ylabel "Fractional Abundance of Species"

set xrange [*:*] reverse
set yrange [*:1.5]

set format y "10^{%T}"
set format x "10^{%T}"

set key inside right bottom

#FILE NUMBER
i = "avg"

    set output "".i.".ps"
    #set output "".i.".pdf"

    # H, H+, e, He, He+
    set title "H, H+, e, He, He+"
    plot "grid.".i.".data" u 1:2 title "H" w lines ls 1 lt rgb "red",\
    '' u 1:3 title "H+" w lines ls 1 lt rgb "blue",\
    '' u 1:4 title "e" w lines ls 1 lt rgb "green",\
    '' u 1:5 title "He" w lines ls 1 lt rgb "orange",\
    '' u 1:6 title "He+" w lines ls 1 lt rgb "cyan"\

    # H-, D, D+, D-
    set title "H-, D, D+, D-"
    plot "grid.".i.".data" u 1:7 title "H-" w lines ls 1 lt rgb "red",\
    '' u 1:8 title "D" w lines ls 1 lt rgb "blue",\
    '' u 1:9 title "D+" w lines ls 1 lt rgb "green",\
    '' u 1:10 title "D-" w lines ls 1 lt rgb "orange"

    # Li, Li+, Li-, LiH+, LiH
    set title "Li, Li+, Li-, LiH+"
    plot "grid.".i.".data" u 1:11 title "Li" w lines ls 1 lt rgb "red",\
    '' u 1:12 title "Li+" w lines ls 1 lt rgb "blue",\
    '' u 1:13 title "Li-" w lines ls 1 lt rgb "green",\
    '' u 1:23 title "LiH+" w lines ls 1 lt rgb "orange"\

    # H2+, H2, HeH+, H3+, He2+
    set title "H2+, H2, HeH+, H3+, He2+"
    plot "grid.".i.".data" u 1:14 title "H2+" w lines ls 1 lt rgb "red",\
    '' u 1:15 title "H2" w lines ls 1 lt rgb "blue",\
    '' u 1:16 title "HeH+" w lines ls 1 lt rgb "green",\
    '' u 1:17 title "H3+" w lines ls 1 lt rgb "orange",\
    '' u 1:22 title "He2+" w lines ls 1 lt rgb "cyan"\

    # HD+, HD, HeD+, H2D+
    set title "HD+, HD, HeD+, H2D+"
    plot "grid.".i.".data" u 1:18 title "HD+" w lines ls 1 lt rgb "red",\
    '' u 1:19 title "HD" w lines ls 1 lt rgb "blue",\
    '' u 1:20 title "HeD+" w lines ls 1 lt rgb "green",\
    '' u 1:21 title "H2D+" w lines ls 1 lt rgb "orange"\


    # D2, D2+, HD2+, D3+
    set title "D2, D2+, HD2+, D3+"
    plot "grid.".i.".data" u 1:25 title "D2" w lines ls 1 lt rgb "red",\
    '' u 1:26 title "D2+" w lines ls 1 lt rgb "blue",\
    '' u 1:27 title "HD2+" w lines ls 1 lt rgb "green",\
    '' u 1:28 title "D3+" w lines ls 1 lt rgb "orange",\
    '' u 1:24 title "LiH" w lines ls 1 lt rgb "cyan",

    # H*, D*
    set title "H*, D*"
    plot "grid.".i.".data" u 1:29 title "H*" w lines ls 1 lt rgb "red",\
    '' u 1:30 title "D*" w lines ls 1 lt rgb "blue"

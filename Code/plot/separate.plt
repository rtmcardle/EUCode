#set terminal postscript enhanced color
set terminal pdf
set logscale y
set style line 1 lw 3

set title font 'helvetica,35'

set xlabel "Redshift (z)" font 'helvetica,25' offset 0,-1.5
set ylabel "Fractional Abundance of Species" font 'helvetica,25' offset -2.0,0

set xrange [31:10]
set yrange [*:1.5]

set format y "10^{%T}"
#set format x "10^{%T}"

set key at 25,10**-2
set key spacing 1.5

#FILE NUMBER
i = 10
f = "avg"

    #set output "plot.".i.".ps"
    set output "plot.".i.".pdf"

    # H, H+, e, He, He+
    set title "H, H^+, e, He, He^+"
    set yrange [10**-3.5:10**0]
    plot "grid.".f.".data" u 1:2 title "H" w lines ls 1 lt rgb "red",\
    '' u 1:3 title "H^+" w lines ls 1 lt rgb "blue",\
    '' u 1:4 title "e" w lines ls 1 lt rgb "green",\
    '' u 1:5 title "He" w lines ls 1 lt rgb "orange",\
    '' u 1:6 title "He^+" w lines ls 1 lt rgb "cyan",\

    i = i+1
    #set output "plot.".i.".ps"
    set output "plot.".i.".pdf"

    # H-, D, D+, D-
    set title "H^-, D, D^+, D^-"
    set yrange[10**-25:10**-2.5]
    set key at 25,10**-18
    plot "grid.".f.".data" u 1:7 title "H^-" w lines ls 1 lt rgb "red",\
    '' u 1:8 title "D" w lines ls 1 lt rgb "blue",\
    '' u 1:9 title "D^+" w lines ls 1 lt rgb "green",\
    '' u 1:10 title "D^-" w lines ls 1 lt rgb "orange"

    i = i+1
    #set output "plot.".i.".ps"
    set output "plot.".i.".pdf"

    # Li, Li^+, Li^-, LiH^+, LiH
    set title "Li, Li^+, Li^-, LiH^+"
    set yrange[10**-25:10**-7.5]
    set key at 25,10**-20
    plot "grid.".f.".data" u 1:11 title "Li" w lines ls 1 lt rgb "red",\
    '' u 1:12 title "Li^+" w lines ls 1 lt rgb "blue",\
    '' u 1:13 title "Li^-" w lines ls 1 lt rgb "green",\
    '' u 1:23 title "LiH^+" w lines ls 1 lt rgb "orange",\

    i = i+1
    #set output "plot.".i.".ps"
    set output "plot.".i.".pdf"

    # H2+, H2, HeH+, H3+, He2+
    set title "H@_2^+, H_2, HeH^+, H@_3^+, He@_2^+"
    set yrange[10**-25:10**-5]
    set key at 25,10**-18
    plot "grid.".f.".data" u 1:14 title "H@_2^+" w lines ls 1 lt rgb "red",\
    '' u 1:15 title "H_2" w lines ls 1 lt rgb "blue",\
    '' u 1:16 title "HeH^+" w lines ls 1 lt rgb "green",\
    '' u 1:17 title "H@_3^+" w lines ls 1 lt rgb "orange",\
    '' u 1:22 title "He@_2^+" w lines ls 1 lt rgb "cyan",\

    i = i+1
    #set output "plot.".i.".ps"
    set output "plot.".i.".pdf"

    # HD+, HD, HeD+, H2D+
    set title "HD^+, HD, HeD^+, H_2D^+"
    set yrange[10**-25:10**-7.5]
    set key at 25,10**-20
    plot "grid.".f.".data" u 1:18 title "HD^+" w lines ls 1 lt rgb "red",\
    '' u 1:19 title "HD" w lines ls 1 lt rgb "blue",\
    '' u 1:20 title "HeD^+" w lines ls 1 lt rgb "green",\
    '' u 1:21 title "H_2D^+" w lines ls 1 lt rgb "orange",\


    ## D2, D2+, HD2+, D3+
    #set title "D_2, D@_2^+, HD2^+, D@_3^+, LiH"
    #plot "grid.".f.".data" u 1:25 title "D_2" w lines ls 1 lt rgb "red",\
    #'' u 1:26 title "D@_2^+" w lines ls 1 lt rgb "blue",\
    #'' u 1:27 title "HD@_2^+" w lines ls 1 lt rgb "green",\
    #'' u 1:28 title "D@_3^+" w lines ls 1 lt rgb "orange",\
    #'' u 1:24 title "LiH" w lines ls 1 lt rgb "cyan",
    #
    ## H*, D*
    #set title "H*, D*"
    #plot "fort.".f.".16" u 1:29 title "H*" w lines ls 1 lt rgb "red",\
    #'' u 1:30 title "D*" w lines ls 1 lt rgb "blue"
    #

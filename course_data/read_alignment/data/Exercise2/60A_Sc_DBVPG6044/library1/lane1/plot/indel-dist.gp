
        set terminal png size 600,400 truecolor
        set output "plot/indel-dist.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style increment user
        set ylabel "Indel count [log]"
        set xlabel "Indel length"
        set y2label "Insertions/Deletions ratio"
        set log y
        set y2tics nomirror
        set ytics nomirror
        set title "lane1.stats.txt" noenhanced
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	5649
2	1286
3	1015
4	540
5	475
6	463
7	103
8	69
9	69
10	102
11	54
12	63
13	15
14	10
15	36
16	22
17	4
18	21
19	7
20	8
21	7
22	2
23	0
24	2
25	3
26	1
27	3
28	7
29	4
30	3
31	0
32	0
34	0
36	0
40	0
48	0
end
1	7092
2	1781
3	997
4	859
5	555
6	349
7	404
8	118
9	206
10	104
11	27
12	225
13	48
14	80
15	29
16	147
17	42
18	37
19	38
20	23
21	12
22	8
23	4
24	21
25	2
26	0
27	3
28	2
29	2
30	1
31	28
32	3
34	7
36	2
40	12
48	1
end
1	0.796531
2	0.722066
3	1.018054
4	0.628638
5	0.855856
6	1.326648
7	0.254950
8	0.584746
9	0.334951
10	0.980769
11	2.000000
12	0.280000
13	0.312500
14	0.125000
15	1.241379
16	0.149660
17	0.095238
18	0.567568
19	0.184211
20	0.347826
21	0.583333
22	0.250000
23	0.000000
24	0.095238
25	1.500000
26	0.000000
27	1.000000
28	3.500000
29	2.000000
30	3.000000
31	0.000000
32	0.000000
34	0.000000
36	0.000000
40	0.000000
48	0.000000
end

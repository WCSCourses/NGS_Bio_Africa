
            set terminal png size 600,400 truecolor
            set output "plot/insert-size.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set rmargin 5
            set label sprintf("%d",273) at 273+10,1625
            set ylabel  "Number of pairs"
            set xlabel  "Insert Size"
            set title "lane1.stats.txt" noenhanced
            plot \
                '-' with lines lc rgb 'black' title 'All pairs', \
                '-' with lines title 'Inward', \
                '-' with lines title 'Outward', \
                '-' with lines title 'Other'
        0	167
1	0
2	11
3	3
4	1
5	2
6	0
7	0
8	0
9	0
10	0
11	1
12	1
13	0
14	1
15	1
16	0
17	0
18	0
19	7
20	5
21	2
22	0
23	2
24	0
25	3
26	1
27	4
28	2
29	2
30	0
31	2
32	3
33	0
34	2
35	5
36	2
37	2
38	0
39	0
40	6
41	4
42	0
43	3
44	3
45	7
46	4
47	3
48	2
49	1
50	2
51	7
52	5
53	4
54	5
55	4
56	3
57	3
58	3
59	4
60	7
61	4
62	7
63	7
64	4
65	3
66	2
67	9
68	5
69	7
70	6
71	7
72	8
73	7
74	6
75	7
76	10
77	5
78	13
79	4
80	17
81	14
82	7
83	10
84	13
85	8
86	12
87	8
88	7
89	9
90	6
91	11
92	14
93	19
94	20
95	18
96	15
97	17
98	14
99	15
100	22
101	24
102	23
103	26
104	19
105	33
106	32
107	26
108	24
109	30
110	49
111	25
112	41
113	40
114	46
115	60
116	44
117	41
118	53
119	59
120	55
121	60
122	52
123	66
124	79
125	71
126	60
127	83
128	64
129	85
130	83
131	89
132	93
133	111
134	100
135	95
136	81
137	109
138	108
139	116
140	93
141	124
142	104
143	117
144	114
145	116
146	125
147	122
148	129
149	147
150	123
151	128
152	115
153	120
154	129
155	147
156	134
157	136
158	130
159	105
160	127
161	129
162	131
163	126
164	133
165	134
166	165
167	129
168	134
169	154
170	143
171	146
172	158
173	135
174	171
175	146
176	151
177	131
178	162
179	151
180	145
181	160
182	153
183	157
184	156
185	144
186	150
187	160
188	183
189	150
190	139
191	144
192	171
193	146
194	163
195	164
196	174
197	163
198	167
199	155
200	160
201	168
202	199
203	185
204	167
205	190
206	170
207	168
208	182
209	186
210	183
211	183
212	174
213	191
214	173
215	188
216	226
217	179
218	214
219	212
220	220
221	205
222	209
223	224
224	199
225	223
226	242
227	218
228	230
229	250
230	258
231	264
232	275
233	294
234	308
235	298
236	334
237	331
238	356
239	357
240	439
241	460
242	478
243	508
244	514
245	600
246	648
247	685
248	687
249	815
250	789
251	851
252	895
253	915
254	1004
255	1103
256	1105
257	1124
258	1220
259	1229
260	1268
261	1284
262	1303
263	1354
264	1313
265	1442
266	1466
267	1505
268	1531
269	1493
270	1474
271	1496
272	1527
273	1625
274	1535
275	1546
276	1553
277	1515
278	1541
279	1557
280	1579
281	1576
282	1614
283	1560
284	1610
285	1495
286	1608
287	1563
288	1606
289	1577
290	1543
291	1537
292	1511
293	1560
294	1589
295	1577
296	1578
297	1517
298	1467
299	1535
300	1523
301	1474
302	1489
303	1474
304	1448
305	1412
306	1359
307	1390
308	1379
309	1318
310	1257
311	1304
312	1274
313	1255
314	1228
315	1157
316	1103
317	1096
318	1021
319	907
320	989
321	928
322	812
323	785
324	779
325	710
326	736
327	688
328	674
329	644
330	610
331	589
332	543
333	502
334	522
335	487
336	421
337	463
338	421
339	385
340	352
341	353
342	362
343	352
344	307
345	301
346	257
347	250
348	249
349	234
350	238
351	212
352	204
353	227
354	189
355	177
356	163
357	149
358	160
359	128
360	131
361	105
362	103
363	94
364	87
365	96
366	81
367	71
368	83
369	65
370	55
371	55
372	56
373	43
374	43
375	33
376	40
377	37
378	34
379	29
380	30
381	24
382	22
383	19
384	25
385	22
386	19
387	13
388	16
389	14
390	14
391	10
392	15
393	7
394	12
395	11
396	6
397	8
398	8
399	9
400	9
401	13
402	4
403	8
404	4
405	5
406	8
407	11
408	3
409	6
410	4
411	9
412	2
413	4
414	4
415	6
416	2
417	4
418	4
419	7
end
0	0
1	0
2	0
3	0
4	0
5	0
6	0
7	0
8	0
9	0
10	0
11	0
12	0
13	0
14	0
15	0
16	0
17	0
18	0
19	5
20	2
21	1
22	0
23	0
24	0
25	0
26	0
27	2
28	0
29	0
30	0
31	0
32	1
33	0
34	1
35	3
36	0
37	1
38	0
39	0
40	1
41	1
42	0
43	0
44	2
45	3
46	0
47	0
48	1
49	1
50	0
51	1
52	3
53	2
54	0
55	3
56	0
57	0
58	0
59	0
60	1
61	0
62	3
63	3
64	0
65	1
66	2
67	3
68	1
69	3
70	3
71	0
72	3
73	2
74	1
75	3
76	2
77	3
78	3
79	2
80	7
81	7
82	5
83	4
84	6
85	2
86	7
87	3
88	2
89	5
90	4
91	8
92	7
93	10
94	8
95	7
96	7
97	9
98	7
99	4
100	16
101	9
102	14
103	9
104	7
105	15
106	17
107	5
108	14
109	28
110	47
111	24
112	40
113	40
114	43
115	59
116	44
117	36
118	51
119	58
120	54
121	57
122	50
123	64
124	77
125	68
126	60
127	79
128	63
129	83
130	81
131	89
132	91
133	106
134	99
135	93
136	80
137	108
138	106
139	116
140	90
141	122
142	102
143	114
144	112
145	111
146	123
147	118
148	126
149	145
150	122
151	125
152	114
153	120
154	129
155	146
156	132
157	135
158	121
159	103
160	127
161	128
162	130
163	123
164	133
165	134
166	163
167	129
168	134
169	150
170	140
171	146
172	154
173	133
174	170
175	146
176	151
177	131
178	160
179	151
180	142
181	156
182	152
183	157
184	154
185	143
186	149
187	159
188	180
189	149
190	138
191	144
192	171
193	145
194	161
195	162
196	173
197	163
198	167
199	154
200	160
201	167
202	199
203	183
204	167
205	190
206	169
207	167
208	180
209	185
210	183
211	183
212	174
213	190
214	173
215	185
216	224
217	178
218	212
219	212
220	220
221	204
222	209
223	224
224	199
225	221
226	241
227	217
228	228
229	250
230	257
231	264
232	274
233	290
234	306
235	297
236	334
237	331
238	355
239	356
240	439
241	459
242	478
243	508
244	514
245	600
246	646
247	685
248	687
249	815
250	789
251	850
252	895
253	915
254	1004
255	1103
256	1105
257	1124
258	1220
259	1229
260	1267
261	1284
262	1303
263	1353
264	1313
265	1442
266	1465
267	1504
268	1531
269	1493
270	1474
271	1496
272	1526
273	1625
274	1535
275	1546
276	1553
277	1515
278	1541
279	1557
280	1578
281	1576
282	1614
283	1560
284	1610
285	1495
286	1608
287	1563
288	1606
289	1577
290	1543
291	1537
292	1511
293	1560
294	1589
295	1577
296	1578
297	1517
298	1467
299	1535
300	1523
301	1474
302	1489
303	1474
304	1448
305	1412
306	1359
307	1390
308	1378
309	1318
310	1257
311	1304
312	1274
313	1255
314	1228
315	1157
316	1103
317	1096
318	1021
319	907
320	988
321	928
322	812
323	785
324	779
325	710
326	736
327	688
328	674
329	644
330	610
331	589
332	543
333	502
334	522
335	487
336	421
337	463
338	421
339	385
340	351
341	353
342	362
343	352
344	307
345	301
346	257
347	250
348	249
349	234
350	238
351	212
352	204
353	227
354	189
355	177
356	163
357	149
358	160
359	128
360	131
361	105
362	103
363	94
364	87
365	96
366	81
367	71
368	83
369	65
370	55
371	55
372	56
373	43
374	43
375	33
376	40
377	37
378	34
379	28
380	30
381	24
382	22
383	19
384	24
385	22
386	19
387	12
388	16
389	14
390	14
391	10
392	15
393	7
394	12
395	11
396	6
397	8
398	8
399	8
400	9
401	13
402	4
403	8
404	4
405	5
406	8
407	10
408	3
409	6
410	4
411	8
412	2
413	4
414	4
415	6
416	2
417	4
418	4
419	7
end
0	0
1	0
2	1
3	0
4	0
5	0
6	0
7	0
8	0
9	0
10	0
11	1
12	0
13	0
14	1
15	1
16	0
17	0
18	0
19	1
20	2
21	0
22	0
23	1
24	0
25	3
26	1
27	1
28	1
29	0
30	0
31	1
32	1
33	0
34	0
35	0
36	1
37	1
38	0
39	0
40	3
41	2
42	0
43	2
44	0
45	3
46	2
47	1
48	1
49	0
50	2
51	3
52	2
53	2
54	3
55	0
56	1
57	1
58	2
59	4
60	4
61	4
62	4
63	4
64	3
65	1
66	0
67	3
68	2
69	3
70	3
71	2
72	4
73	4
74	2
75	2
76	6
77	2
78	7
79	2
80	7
81	4
82	2
83	5
84	5
85	6
86	2
87	3
88	3
89	3
90	1
91	2
92	5
93	8
94	9
95	11
96	6
97	7
98	4
99	8
100	6
101	14
102	8
103	14
104	11
105	17
106	15
107	17
108	6
109	0
110	0
111	0
112	0
113	0
114	1
115	0
116	0
117	0
118	0
119	0
120	0
121	1
122	0
123	0
124	0
125	1
126	0
127	0
128	0
129	0
130	0
131	0
132	0
133	0
134	0
135	0
136	0
137	0
138	0
139	0
140	0
141	0
142	0
143	0
144	0
145	1
146	0
147	0
148	0
149	0
150	0
151	0
152	0
153	0
154	0
155	0
156	0
157	0
158	0
159	0
160	0
161	0
162	0
163	0
164	0
165	0
166	0
167	0
168	0
169	0
170	0
171	0
172	0
173	0
174	0
175	0
176	0
177	0
178	0
179	0
180	0
181	0
182	0
183	0
184	0
185	0
186	0
187	0
188	0
189	0
190	0
191	0
192	0
193	0
194	0
195	0
196	0
197	0
198	0
199	1
200	0
201	0
202	0
203	1
204	0
205	0
206	0
207	0
208	0
209	0
210	0
211	0
212	0
213	0
214	0
215	0
216	1
217	0
218	0
219	0
220	0
221	1
222	0
223	0
224	0
225	0
226	1
227	0
228	1
229	0
230	0
231	0
232	0
233	0
234	2
235	0
236	0
237	0
238	0
239	0
240	0
241	0
242	0
243	0
244	0
245	0
246	0
247	0
248	0
249	0
250	0
251	0
252	0
253	0
254	0
255	0
256	0
257	0
258	0
259	0
260	0
261	0
262	0
263	1
264	0
265	0
266	0
267	0
268	0
269	0
270	0
271	0
272	0
273	0
274	0
275	0
276	0
277	0
278	0
279	0
280	0
281	0
282	0
283	0
284	0
285	0
286	0
287	0
288	0
289	0
290	0
291	0
292	0
293	0
294	0
295	0
296	0
297	0
298	0
299	0
300	0
301	0
302	0
303	0
304	0
305	0
306	0
307	0
308	1
309	0
310	0
311	0
312	0
313	0
314	0
315	0
316	0
317	0
318	0
319	0
320	0
321	0
322	0
323	0
324	0
325	0
326	0
327	0
328	0
329	0
330	0
331	0
332	0
333	0
334	0
335	0
336	0
337	0
338	0
339	0
340	1
341	0
342	0
343	0
344	0
345	0
346	0
347	0
348	0
349	0
350	0
351	0
352	0
353	0
354	0
355	0
356	0
357	0
358	0
359	0
360	0
361	0
362	0
363	0
364	0
365	0
366	0
367	0
368	0
369	0
370	0
371	0
372	0
373	0
374	0
375	0
376	0
377	0
378	0
379	0
380	0
381	0
382	0
383	0
384	0
385	0
386	0
387	0
388	0
389	0
390	0
391	0
392	0
393	0
394	0
395	0
396	0
397	0
398	0
399	1
400	0
401	0
402	0
403	0
404	0
405	0
406	0
407	0
408	0
409	0
410	0
411	0
412	0
413	0
414	0
415	0
416	0
417	0
418	0
419	0
end
0	167
1	0
2	10
3	3
4	1
5	2
6	0
7	0
8	0
9	0
10	0
11	0
12	1
13	0
14	0
15	0
16	0
17	0
18	0
19	1
20	1
21	1
22	0
23	1
24	0
25	0
26	0
27	1
28	1
29	2
30	0
31	1
32	1
33	0
34	1
35	2
36	1
37	0
38	0
39	0
40	2
41	1
42	0
43	1
44	1
45	1
46	2
47	2
48	0
49	0
50	0
51	3
52	0
53	0
54	2
55	1
56	2
57	2
58	1
59	0
60	2
61	0
62	0
63	0
64	1
65	1
66	0
67	3
68	2
69	1
70	0
71	5
72	1
73	1
74	3
75	2
76	2
77	0
78	3
79	0
80	3
81	3
82	0
83	1
84	2
85	0
86	3
87	2
88	2
89	1
90	1
91	1
92	2
93	1
94	3
95	0
96	2
97	1
98	3
99	3
100	0
101	1
102	1
103	3
104	1
105	1
106	0
107	4
108	4
109	2
110	2
111	1
112	1
113	0
114	2
115	1
116	0
117	5
118	2
119	1
120	1
121	2
122	2
123	2
124	2
125	2
126	0
127	4
128	1
129	2
130	2
131	0
132	2
133	5
134	1
135	2
136	1
137	1
138	2
139	0
140	3
141	2
142	2
143	3
144	2
145	4
146	2
147	4
148	3
149	2
150	1
151	3
152	1
153	0
154	0
155	1
156	2
157	1
158	9
159	2
160	0
161	1
162	1
163	3
164	0
165	0
166	2
167	0
168	0
169	4
170	3
171	0
172	4
173	2
174	1
175	0
176	0
177	0
178	2
179	0
180	3
181	4
182	1
183	0
184	2
185	1
186	1
187	1
188	3
189	1
190	1
191	0
192	0
193	1
194	2
195	2
196	1
197	0
198	0
199	0
200	0
201	1
202	0
203	1
204	0
205	0
206	1
207	1
208	2
209	1
210	0
211	0
212	0
213	1
214	0
215	3
216	1
217	1
218	2
219	0
220	0
221	0
222	0
223	0
224	0
225	2
226	0
227	1
228	1
229	0
230	1
231	0
232	1
233	4
234	0
235	1
236	0
237	0
238	1
239	1
240	0
241	1
242	0
243	0
244	0
245	0
246	2
247	0
248	0
249	0
250	0
251	1
252	0
253	0
254	0
255	0
256	0
257	0
258	0
259	0
260	1
261	0
262	0
263	0
264	0
265	0
266	1
267	1
268	0
269	0
270	0
271	0
272	1
273	0
274	0
275	0
276	0
277	0
278	0
279	0
280	1
281	0
282	0
283	0
284	0
285	0
286	0
287	0
288	0
289	0
290	0
291	0
292	0
293	0
294	0
295	0
296	0
297	0
298	0
299	0
300	0
301	0
302	0
303	0
304	0
305	0
306	0
307	0
308	0
309	0
310	0
311	0
312	0
313	0
314	0
315	0
316	0
317	0
318	0
319	0
320	1
321	0
322	0
323	0
324	0
325	0
326	0
327	0
328	0
329	0
330	0
331	0
332	0
333	0
334	0
335	0
336	0
337	0
338	0
339	0
340	0
341	0
342	0
343	0
344	0
345	0
346	0
347	0
348	0
349	0
350	0
351	0
352	0
353	0
354	0
355	0
356	0
357	0
358	0
359	0
360	0
361	0
362	0
363	0
364	0
365	0
366	0
367	0
368	0
369	0
370	0
371	0
372	0
373	0
374	0
375	0
376	0
377	0
378	0
379	1
380	0
381	0
382	0
383	0
384	1
385	0
386	0
387	1
388	0
389	0
390	0
391	0
392	0
393	0
394	0
395	0
396	0
397	0
398	0
399	0
400	0
401	0
402	0
403	0
404	0
405	0
406	0
407	1
408	0
409	0
410	0
411	1
412	0
413	0
414	0
415	0
416	0
417	0
418	0
419	0
end

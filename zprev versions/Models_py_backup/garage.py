#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""
#import csv
#with open('Amanullah.txt') as f:
#    reader = csv.reader(f, delimiter="\t")
#    d = list(reader)
#print(d[5][:]) # 248

#import pandas as pd
#pd.read_csv('Amanullah.txt', delim_whitespace=True)

import numpy as np
zpicks1 =[0.028488		
        ,0.050043		
        ,0.052926		
        ,0.070086		
        ,0.062668		
        ,0.087589		
        ,0.078577		
        ,0.017227		
        ,0.042233		
        ,0.045295		
        ,0.03648		
        ,0.019599		
        ,0.100915		
        ,0.027342		
        ,0.074605		
        ,0.026489		
        ,0.049922		
        ,0.030604		
        ,0.016345641	
        ,0.0154363		
        ,0.030529		
        ,0.024525		
        ,0.023953		
        ,0.026038		
        ,0.048948		
        ,0.024314		
        ,0.015166		
        ,0.03572		
        ,0.048818		
        ,0.0219800059146
        ,0.0275			
        ,0.1244			
        ,0.036			
        ,0.01673		
        ,0.016321		
        ,0.021793		
        ,0.01645		
        ,0.023208		
        ,0.036457		
        ,0.019264		
        ,0.017605		
        ,0.031528		
        ,0.023536		
        ,0.016743		
        ,0.05371		
        ,0.016991		
        ,0.027865		
        ,0.017173		
        ,0.029955		
        ,0.016559		
        ,0.015			
        ,0.0544			
        ,0.1561			
        ,0.0393			
        ,0.1241			
        ,0.1441			
        ,0.1299			
        ,0.0784			
        ,0.0583			
        ,0.0309			
        ,0.0406			
        ,0.0152			
        ,0.0224			
        ,0.016			
        ,0.0362			
        ,0.0173			
        ,0.0312			
        ,0.0221			
        ,0.016			
        ,0.0249			
        ,0.0303			
        ,0.0283			
        ,0.0152			
        ,0.0345			
        ,0.036			
        ,0.0248			
        ,0.0292			
        ,0.0163			
        ,0.0187			
        ,0.0195			
        ,0.0256			
        ,0.0337			
        ,0.0546			
        ,0.024			
        ,0.0336			
        ,0.0341			
        ,0.0261			
        ,0.0211			
        ,0.0321			
        ,0.0221			
        ,0.0298			
        ,0.0334			
        ,0.0284			
        ,0.0341			
        ,0.045			
        ,0.0421			
        ,0.0576			
        ,0.033			
        ,0.0753			
        ,0.0204			
        ,0.0205			
        ,0.0402			
        ,0.026			
        ,0.0259			
        ,0.0268			
        ,0.0239			
        ,0.069			
        ,0.0651			
        ,0.0229			
        ,0.018			
        ,0.0315			
        ,0.0215			
        ,0.0255			
        ,0.0325			
        ,0.0843			
        ,0.0308			
        ,0.0327			
        ,0.0423			
        ,0.0684			
        ,0.0153			
        ,0.0233			
        ,0.0491			
        ,0.0425			
        ,0.0192			
        ,0.0308			
        ,0.0212			
        ,0.0277			
        ,0.0335			
        ,0.0208			
        ,0.0173			
        ,0.036			
        ,0.0233			
        ,0.0589			
        ,0.0583			
        ,0.0688			
        ,0.0321			
        ,0.0522			
        ,0.0308			
        ,0.0329			
        ,0.023			
        ,0.015			
        ,0.0321			
        ,0.0643			
        ,0.032			
        ,0.0209			
        ,0.0219			
        ,0.032			
        ,0.0151			
        ,0.0192			
        ,0.0266			
        ,0.0377			
        ,0.0247			
        ,0.0242			
        ,0.0366			
        ,0.0229			
        ,0.0312			
        ,0.015			
        ,0.0341			
        ,0.0251			
        ,0.0189			
        ,0.065			
        ,0.147			
        ,0.13			
        ,0.104			
        ,0.244			
        ,0.3			
        ,0.045			
        ,0.114			
        ,0.258			
        ,0.297			
        ,0.382			
        ,0.147			
        ,0.274			
        ,0.299			
        ,0.38			
        ,0.382			
        ,0.303			
        ,0.35			
        ,0.087			
        ,0.262			
        ,0.217			
        ,0.119			
        ,0.183			
        ,0.359			
        ,0.143			
        ,0.262			
        ,0.234			
        ,0.153			
        ,0.095			
        ,0.288			
        ,0.195			
        ,0.148			
        ,0.213			
        ,0.181			
        ,0.265			
        ,0.193			
        ,0.34			
        ,0.118			
        ,0.143			
        ,0.162			
        ,0.29			
        ,0.124			
        ,0.265			
        ,0.127			
        ,0.174			
        ,0.165			
        ,0.251			
        ,0.259			
        ,0.108			
        ,0.161			
        ,0.206			
        ,0.245			
        ,0.25			
        ,0.23			
        ,0.087			
        ,0.063			
        ,0.279			
        ,0.277			
        ,0.156			
        ,0.332			
        ,0.363			
        ,0.332			
        ,0.146			
        ,0.39			
        ,0.175			
        ,0.301			
        ,0.117			
        ,0.22			
        ,0.121			
        ,0.156			
        ,0.41			
        ,0.179			
        ,0.252			
        ,0.253			
        ,0.393			
        ,0.13			
        ,0.401			
        ,0.311			
        ,0.172			
        ,0.28			
        ,0.31			
        ,0.187			
        ,0.067			
        ,0.318			
        ,0.259			
        ,0.3			
        ,0.272			
        ,0.281			
        ,0.294			
        ,0.19			
        ,0.267			
        ,0.125			
        ,0.184			
        ,0.314			
        ,0.311			
        ,0.09			
        ,0.404			
        ,0.202			
        ,0.328			
        ,0.213			
        ,0.181			
        ,0.094			
        ,0.304			
        ,0.11			
        ,0.087			
        ,0.204			
        ,0.198			
        ,0.128			
        ,0.216			
        ,0.322			
        ,0.219			
        ,0.191			
        ,0.094			
        ,0.381			
        ,0.212			
        ,0.368			
        ,0.422			
        ,0.259			
        ,0.185			
        ,0.214			
        ,0.361			
        ,0.395			
        ,0.116			
        ,0.145			
        ,0.254			
        ,0.389			
        ,0.35			
        ,0.257			
        ,0.218			
        ,0.43			
        ,0.62			
        ,0.57			
        ,0.3			
        ,0.38			
        ,0.43			
        ,0.24			
        ,0.44			
        ,0.5			
        ,0.97			
        ,0.479			
        ,0.83			
        ,0.416			
        ,0.581			
        ,0.45			
        ,0.579			
        ,0.32			
        ,0.657			
        ,0.43			
        ,0.472			
        ,0.374			
        ,0.18			
        ,0.55			
        ,0.592			
        ,0.172			
        ,0.526			
        ,0.763			
        ,0.58			
        ,0.43			
        ,0.45			
        ,0.656			
        ,0.495			
        ,0.49			
        ,0.57			
        ,0.388			
        ,0.45			
        ,0.48			
        ,0.615			
        ,0.4			
        ,0.655			
        ,0.498			
        ,0.465			
        ,0.453			
        ,0.425			
        ,0.514			
        ,0.423			
        ,0.859			
        ,0.936			
        ,0.528			
        ,0.978			
        ,0.885			
        ,0.815			
        ,0.698			
        ,0.568			
        ,0.711			
        ,0.3396			
        ,0.3965			
        ,0.812			
        ,0.799			
        ,0.882			
        ,0.833			
        ,0.874			
        ,0.772			
        ,0.178			
        ,0.26			
        ,0.186			
        ,0.269			
        ,0.215			
        ,0.543			
        ,0.75			
        ,0.64			
        ,0.43			
        ,0.64			
        ,0.497			
        ,0.44			
        ,0.355			
        ,0.78			
        ,0.54			
        ,0.86			
        ,0.468			
        ,0.84			
        ,0.96			
        ,0.8218			
        ,0.93			
        ,0.451			
        ,0.61			
        ,0.83			
        ,0.707			
        ,0.415			
        ,0.557			
        ,0.791			
        ,0.695			
        ,0.633			
        ,0.2486			
        ,0.532			
        ,0.331			
        ,0.346			
        ,0.961			
        ,0.613			
        ,0.3402			
        ,0.983			
        ,0.71			
        ,0.73			
        ,0.47			
        ,0.62			
        ,0.521			
        ,0.369			
        ,0.571			
        ,0.604			
        ,0.9271			
        ,0.285			
        ,0.2912			
        ,0.548			
        ,0.868			
        ,0.496			
        ,0.811			
        ,0.756			
        ,0.817			
        ,0.752			
        ,0.5516			
        ,0.3578			
        ,1.01			
        ,0.741			
        ,0.43			
        ,0.526			
        ,0.592			
        ,0.905			
        ,0.949			
        ,0.4607			
        ,0.3709			
        ,0.8			
        ,0.679			
        ,0.5817			
        ,0.55			
        ,0.81			
        ,0.95			
        ,0.3373			
        ,0.91			
        ,0.263			
        ,0.643			
        ,0.691			
        ,0.357			
        ,0.721			
        ,0.581			
        ,0.6268			
        ,0.818			
        ,0.4627			
        ,0.449			
        ,0.688			
        ,0.87			
        ,0.5043			
        ,0.591			
        ,0.426			
        ,0.329			
        ,0.583			
        ,0.519			
        ,0.401			
        ,0.205			
        ,0.34			
        ,0.436			
        ,0.363			
        ,0.436			
        ,0.309			
        ,0.342			
        ,0.159			
        ,0.332			
        ,0.469			
        ,0.239			
        ,0.352			
        ,0.612			
        ,0.631			
        ,0.645			
        ,0.429			
        ,0.497			
        ,0.539			
        ,0.561			
        ,0.41			
        ,0.412			
        ,0.599			
        ,0.619			
        ,0.422			
        ,0.54			
        ,0.401			
        ,0.218			
        ,0.633			
        ,0.383			
        ,0.302			
        ,0.34			
        ,0.51			
        ,0.421			
        ,0.399			
        ,0.493			
        ,0.687			
        ,0.687			
        ,0.495			
        ,0.603			
        ,0.421			
        ,0.348			
        ,0.213			
        ,0.344			
        ,0.271			
        ,0.564			
        ,0.274			
        ,0.181			
        ,0.582			
        ,0.68			
        ,0.401			
        ,0.416			
        ,0.286			
        ,0.562			
        ,0.266			
        ,0.314			
        ,0.581			
        ,0.463			
        ,0.341			
        ,0.631			
        ,0.522			
        ,0.368			
        ,0.309			
        ,0.528			
        ,0.216			
        ,0.284			
        ,0.508			
        ,0.781			
        ,0.613			
        ,0.278			
        ,0.477			
        ,0.95			
        ,1.057			
        ,0.816			
        ,0.455			
        ,1.02			
        ,1.14			
        ,0.854			
        ,1.37			
        ,0.975			
        ,0.97			
        ,0.74			
        ,1.39			
        ,0.46			
        ,1.02			
        ,1.12			
        ,1.23			
        ,1.19			
        ,0.839			
        ,1.01			
        ,0.521			
        ,0.475			
        ,0.95			
        ,1.3			
        ,1.4			
        ,1.305			
        ,0.216			
        ,0.735			
        ,1.14			
        ,1.307			
        ,1.265			
        ,0.67			
        ,0.64			
        ,1.34			
        ,0.84			
        ,0.935			
        ,0.953			
        ,1.124			
        ,0.552			
        ,0.671			
        ,0.511			
        ,1.03]

zpicks2 = [0.226143925635
        ,0.167114251513
        ,0.155866408242
        ,0.158669048234
        ,0.156270379521
        ,0.189597973592
        ,0.155790553548
        ,0.199535140164
        ,0.167648033398
        ,0.165271499987
        ,0.170400148449
        ,0.184736519598
        ,0.167818465426
        ,0.175738040875
        ,0.160110641392
        ,0.19147874559
        ,0.162480236617
        ,0.173537588946
        ,0.14402997371
        ,0.149773066605
        ,0.0905547451176
        ,0.10395529982
        ,0.108670309216
        ,0.110164462529
        ,0.175878781605
        ,0.184772695826
        ,0.21780640315
        ,0.17455423199
        ,0.163660035553
        ,0.190634594378
        ,0.179523627555
        ,0.167742220247
        ,0.171121785922
        ,0.212053128352
        ,0.207622958381
        ,0.231469837469
        ,0.250007930671
        ,0.230453431037
        ,0.216601021623
        ,0.239750047895
        ,0.270777214079
        ,0.224419994663
        ,0.233717264297
        ,0.247451826306
        ,0.220796341332
        ,0.308392297269
        ,0.221390433935
        ,0.247302669853
        ,0.223360506658
        ,0.25040542174
        ,0.160930618134
        ,0.0851859276563
        ,0.0831803243516
        ,0.0983616480332
        ,0.110959981085
        ,0.156387993518
        ,0.129040616704
        ,0.0863760987398
        ,0.199697931638
        ,0.176677068679
        ,0.164649900237
        ,0.209050001205
        ,0.232580853069
        ,0.215338744788
        ,0.164155647932
        ,0.210066600914
        ,0.173133567729
        ,0.183198718389
        ,0.201710026772
        ,0.187033320457
        ,0.169350596107
        ,0.172739361028
        ,0.236356060735
        ,0.205416560749
        ,0.171522593533
        ,0.177535782419
        ,0.167312607114
        ,0.206892325236
        ,0.190583331504
        ,0.191040914154
        ,0.179228108521
        ,0.172930780679
        ,0.16950794454
        ,0.189072007408
        ,0.181737645865
        ,0.167946113517
        ,0.185600393489
        ,0.18146401803
        ,0.166313999169
        ,0.182214269008
        ,0.209882805071
        ,0.168522663763
        ,0.168292526178
        ,0.165324028128
        ,0.155012862381
        ,0.162357046505
        ,0.153692702221
        ,0.164165671037
        ,0.152378958895
        ,0.188558589288
        ,0.184564810559
        ,0.162967150273
        ,0.176878840961
        ,0.171154394323
        ,0.17097815303
        ,0.174177026864
        ,0.169085209011
        ,0.154689426092
        ,0.178350941227
        ,0.190930640521
        ,0.163304382923
        ,0.180036324178
        ,0.189515900553
        ,0.162662308732
        ,0.193449865152
        ,0.171179303431
        ,0.161712348214
        ,0.16403631389
        ,0.162821861461
        ,0.207155778089
        ,0.176317048358
        ,0.168279067414
        ,0.188970108359
        ,0.191424646957
        ,0.165778427196
        ,0.18474579249
        ,0.176243118914
        ,0.162762916979
        ,0.195361399576
        ,0.215091501139
        ,0.169897387897
        ,0.178698403151
        ,0.157634614976
        ,0.160351225061
        ,0.191426088364
        ,0.167577472063
        ,0.182864275271
        ,0.164890394315
        ,0.162681757891
        ,0.179840694371
        ,0.211404589303
        ,0.166296974492
        ,0.158132776968
        ,0.186494110705
        ,0.187650195008
        ,0.181514604995
        ,0.165934850147
        ,0.207587590171
        ,0.190576794911
        ,0.170306255652
        ,0.290373531881
        ,0.174379321042
        ,0.172887454113
        ,0.16122683881
        ,0.175847772194
        ,0.166284241086
        ,0.207266601483
        ,0.170422924262
        ,0.171895520078
        ,0.20380344245
        ,0.120378117628
        ,0.121080634003
        ,0.135784293243
        ,0.12064387065
        ,0.150780305961
        ,0.215824546738
        ,0.127260196442
        ,0.120017551192
        ,0.157083138554
        ,0.215462354547
        ,0.220813701877
        ,0.126442356391
        ,0.164190866001
        ,0.204243298637
        ,0.197647141922
        ,0.199423960287
        ,0.246333608197
        ,0.218357189469
        ,0.117239118491
        ,0.142760435008
        ,0.140960118616
        ,0.114480614211
        ,0.120171462823
        ,0.202932960788
        ,0.12033040226
        ,0.149344667277
        ,0.152379820049
        ,0.114718538439
        ,0.118788973643
        ,0.134501367232
        ,0.127181263793
        ,0.112499009039
        ,0.16065563488
        ,0.116442235995
        ,0.128294661027
        ,0.121960815999
        ,0.17167768821
        ,0.114590986149
        ,0.114524011818
        ,0.119565350744
        ,0.150676952536
        ,0.122100457699
        ,0.140962475796
        ,0.120227176144
        ,0.114282533804
        ,0.113957220599
        ,0.133182501728
        ,0.136851565925
        ,0.118564236278
        ,0.115619678811
        ,0.124355782411
        ,0.12515209647
        ,0.126539062794
        ,0.122886919759
        ,0.113367509506
        ,0.116292465656
        ,0.170515825332
        ,0.143749208555
        ,0.116301202618
        ,0.157824605437
        ,0.167338045039
        ,0.153791748359
        ,0.113634577753
        ,0.202376651268
        ,0.129902757627
        ,0.152785577721
        ,0.113090217803
        ,0.117918087643
        ,0.111553736186
        ,0.112833879905
        ,0.208370772279
        ,0.125612605347
        ,0.121669684179
        ,0.119228104024
        ,0.171808246679
        ,0.110314772736
        ,0.20749905315
        ,0.135622441007
        ,0.117223091751
        ,0.13095448346
        ,0.177971613506
        ,0.12145112418
        ,0.118306255874
        ,0.192756486723
        ,0.132143781265
        ,0.155174252389
        ,0.146266612538
        ,0.139677729246
        ,0.142023222928
        ,0.120475217309
        ,0.125621852074
        ,0.11193654747
        ,0.114893026219
        ,0.140263606877
        ,0.154509467451
        ,0.118243252284
        ,0.212015765825
        ,0.133809098038
        ,0.144862922481
        ,0.11609784707
        ,0.117868646113
        ,0.115343190964
        ,0.135822682884
        ,0.120481707591
        ,0.118814561252
        ,0.129617980863
        ,0.13579473086
        ,0.124812931046
        ,0.117827563103
        ,0.135100374854
        ,0.12534004292
        ,0.120534772995
        ,0.119858620288
        ,0.148158100576
        ,0.121635351197
        ,0.199892747232
        ,0.203900920127
        ,0.1381090777
        ,0.124767998966
        ,0.132159093773
        ,0.176630102674
        ,0.152039612546
        ,0.114183512043
        ,0.117222451448
        ,0.145717293983
        ,0.188189991384
        ,0.178077262302
        ,0.152073477348
        ,0.154102535988
        ,0.36344840931
        ,0.395579278552
        ,0.395227786467
        ,0.320340946007
        ,0.333518201962
        ,0.466023796791
        ,0.408470067632
        ,0.32528185368
        ,0.322646999342
        ,0.818432977093
        ,0.363960881605
        ,0.476221516489
        ,0.566269702273
        ,0.511751148392
        ,0.45952604589
        ,0.653193752978
        ,0.426633125026
        ,0.668671821574
        ,0.619060682936
        ,0.522213466148
        ,0.943967733848
        ,0.452691548423
        ,1.02320057236
        ,0.727363486932
        ,0.436282215602
        ,0.521010819611
        ,0.909043307134
        ,0.526483896531
        ,0.465243749852
        ,0.598757469559
        ,0.635534544683
        ,0.449365583478
        ,0.452457762216
        ,0.476394763516
        ,0.471680427567
        ,0.51587600329
        ,0.529360126953
        ,0.564010349483
        ,0.483871492231
        ,0.503703328639
        ,0.652414722917
        ,0.610115606889
        ,0.532970328234
        ,0.494543625705
        ,0.506403802203
        ,0.2534244988
        ,0.302526442652
        ,0.717608197762
        ,0.250124518086
        ,0.298988010763
        ,0.29259857127
        ,0.746941209578
        ,0.440949943442
        ,0.295573375025
        ,0.37487817081
        ,0.232651192244
        ,0.207508245548
        ,0.386258428909
        ,0.238374333285
        ,0.597150330975
        ,0.531348208869
        ,0.395564414276
        ,0.223384059982
        ,0.243792731173
        ,0.207808738715
        ,0.197085550242
        ,0.264144092838
        ,0.201887086183
        ,0.0965091021297
        ,0.137773519339
        ,0.189960701162
        ,0.148525009826
        ,0.194442979286
        ,0.168225439785
        ,0.106490971154
        ,0.203457956158
        ,0.171350462953
        ,0.11058193483
        ,0.172897098508
        ,0.165250052616
        ,0.225125911159
        ,0.278620780131
        ,0.215627670046
        ,0.293679174465
        ,0.144017146978
        ,0.15737107879
        ,0.211176417394
        ,0.259277417037
        ,0.141995153329
        ,0.158128342293
        ,0.203347956975
        ,0.207672991039
        ,0.171303725502
        ,0.156001611413
        ,0.167913176451
        ,0.13852607255
        ,0.153627203737
        ,0.345401297016
        ,0.166247667787
        ,0.136350904009
        ,0.441320917548
        ,0.185960509148
        ,0.196706401553
        ,0.154395011149
        ,0.171049818244
        ,0.162167199818
        ,0.15383794776
        ,0.175419559321
        ,0.160870192314
        ,0.290419777723
        ,0.137402734867
        ,0.147860252643
        ,0.190445811927
        ,0.249257746544
        ,0.154020408449
        ,0.212840468683
        ,0.194569634979
        ,0.212381883478
        ,0.20904577048
        ,0.150278973145
        ,0.136956742249
        ,0.380891045794
        ,0.225350070805
        ,0.152228805806
        ,0.168434413662
        ,0.173662227482
        ,0.260816306234
        ,0.288339749303
        ,0.180730669073
        ,0.155882980598
        ,0.217280907113
        ,0.197639177097
        ,0.160984148284
        ,0.156209026367
        ,0.215181910354
        ,0.29513371916
        ,0.138055963607
        ,0.270806993269
        ,0.135415118781
        ,0.167688951886
        ,0.256534768993
        ,0.13857499344
        ,0.188595219767
        ,0.156453875009
        ,0.158306777046
        ,0.27125777036
        ,0.148626981409
        ,0.158124901298
        ,0.166956993834
        ,0.281601232201
        ,0.15061141667
        ,0.308183442634
        ,0.22595475705
        ,0.273050416654
        ,0.288536338632
        ,0.316956305777
        ,0.23053144449
        ,0.219199816494
        ,0.229672711685
        ,0.217391698785
        ,0.20891536552
        ,0.2294063655
        ,0.238102446996
        ,0.215459944632
        ,0.246358907571
        ,0.240709017167
        ,0.283065430909
        ,0.205242395733
        ,0.225170628892
        ,0.284735745242
        ,0.245972086709
        ,0.23539370649
        ,0.216241856053
        ,0.222751847877
        ,0.243223346147
        ,0.307329015793
        ,0.259435785458
        ,0.315186986894
        ,0.386382183638
        ,0.249438141428
        ,0.232044816512
        ,0.284960259957
        ,0.371779826612
        ,0.224122410473
        ,0.261565759588
        ,0.253453216931
        ,0.275320135034
        ,0.232764787176
        ,0.227776883863
        ,0.325882220619
        ,0.296231218043
        ,0.265794738732
        ,0.286826360161
        ,0.274033552724
        ,0.233155755193
        ,0.280738355979
        ,0.23465390113
        ,0.218430762841
        ,0.225411924471
        ,0.210577314192
        ,0.217182931008
        ,0.296179216596
        ,0.209762707255
        ,0.229723298662
        ,0.327843886551
        ,0.292732159101
        ,0.336024813008
        ,0.304687254477
        ,0.260077754247
        ,0.336087850286
        ,0.2525260008
        ,0.243685251902
        ,0.394367550265
        ,0.268565103689
        ,0.225338954944
        ,0.235138397674
        ,0.248621683092
        ,0.207836005267
        ,0.215775825102
        ,0.303416211397
        ,0.217189018845
        ,0.207444974165
        ,0.232554896589
        ,0.325104262678
        ,0.261568700364
        ,0.20080587099
        ,0.185467817051
        ,0.30169867661
        ,0.240606195136
        ,0.466035003635
        ,0.232259103398
        ,0.208732914262
        ,0.216820067565
        ,0.213134015411
        ,0.244915950922
        ,0.198885573363
        ,0.282255398884
        ,0.183311419886
        ,0.238383165249
        ,0.196808697491
        ,0.2262922033
        ,0.21187088553
        ,0.222899260788
        ,0.236655407187
        ,0.209549982904
        ,0.322525273833
        ,0.184037573567
        ,0.248548172881
        ,0.222117064104
        ,0.217324926249
        ,1.07792101581
        ,0.249576425145
        ,0.231050242645
        ,0.185908921392
        ,0.364467879127
        ,0.265400797402
        ,0.224124266884
        ,0.186583218843
        ,0.186744572784
        ,0.267856909216
        ,0.194732639155
        ,0.215358476354
        ,0.971680645776
        ,0.198690905622
        ,0.105933336834
        ,0.121926562564
        ,0.0899118005512
        ,0.142160782046]

mag = [
        35.3355510594
        ,36.6754415683
        ,36.8168806729
        ,37.4403208027
        ,37.4803312569
        ,38.2312880884
        ,37.4880153326
        ,34.6528548205
        ,36.3351409577
        ,36.6329198059
        ,35.9028802813
        ,34.5849105142
        ,38.4564792878
        ,35.0742634622
        ,37.5857951062
        ,35.4739511994
        ,36.5637135925
        ,35.5474875749
        ,34.0357608093
        ,33.9288622487
        ,35.5987234724
        ,35.0588114943
        ,34.9628931602
        ,35.3604256119
        ,36.7224923466
        ,35.0989067311
        ,34.0945752677
        ,35.9784201813
        ,36.3736376383
        ,34.8477997627
        ,35.6443099344
        ,39.0474440635
        ,35.8166920463
        ,34.2137334689
        ,33.9968666808
        ,34.9662747364
        ,34.1782747061
        ,35.0760231692
        ,36.1310606977
        ,34.9290818261
        ,34.3341878854
        ,35.7209576638
        ,35.1655185904
        ,34.0023106368
        ,36.47359007
        ,34.3727214476
        ,35.079691262
        ,34.2567548437
        ,35.9695297947
        ,34.3359548469
        ,34.1611242178
        ,36.9482868597
        ,39.227736887
        ,36.3296338828
        ,38.8089373233
        ,38.8286718588
        ,38.9767208561
        ,37.6784334873
        ,37.0227081887
        ,35.9151154293
        ,36.3635574234
        ,34.0128062018
        ,34.9604436109
        ,34.1665118574
        ,35.9818330401
        ,34.2481675909
        ,35.6190668111
        ,34.8880671183
        ,33.8142603068
        ,34.7844502384
        ,35.6227445028
        ,35.5095000297
        ,34.2611733597
        ,35.9597476335
        ,35.6729412003
        ,35.2550362369
        ,35.9909847314
        ,34.428922484
        ,35.0398652996
        ,34.7556470033
        ,35.6821627437
        ,35.8268020378
        ,36.6034159881
        ,35.1708581962
        ,36.0054974141
        ,35.8323574252
        ,35.3540611638
        ,34.6494717372
        ,35.8862461317
        ,34.9213754902
        ,35.506995688
        ,35.8612155366
        ,35.5704019595
        ,35.942753551
        ,36.5991688426
        ,36.4032354598
        ,37.1076369832
        ,35.9894165447
        ,37.6492880383
        ,34.7125139232
        ,34.5973526852
        ,36.3669644792
        ,35.3618288289
        ,35.4093397274
        ,35.3032674902
        ,35.0233196812
        ,37.5632893538
        ,37.3353162366
        ,35.1846393132
        ,34.3947741135
        ,35.6474839975
        ,34.9211199992
        ,35.6747779958
        ,35.7976725473
        ,38.0527926414
        ,35.6124448286
        ,36.0723318412
        ,36.3883376162
        ,37.7252349619
        ,34.6566486547
        ,34.8677153894
        ,36.7199743885
        ,35.9144831533
        ,34.7285878322
        ,35.7737272495
        ,34.8346147484
        ,35.685542453
        ,35.9646644449
        ,34.7963237588
        ,34.2290778645
        ,36.1419516923
        ,35.192149373
        ,37.1375683747
        ,37.0511345128
        ,37.4743372416
        ,35.624004108
        ,36.6643549794
        ,35.5802476433
        ,35.9289937858
        ,35.066402753
        ,34.3704032433
        ,35.8655549482
        ,37.1708016197
        ,35.8135243381
        ,34.6890069038
        ,34.8412253332
        ,35.5819126039
        ,34.5091434423
        ,34.4813665248
        ,35.3141855009
        ,35.7917282129
        ,34.9096518953
        ,35.197313794
        ,35.969795049
        ,35.130731713
        ,35.880989798
        ,34.0863548459
        ,35.7604333404
        ,34.9393048575
        ,34.3767780517
        ,37.3262250234
        ,39.5592698132
        ,38.9057647374
        ,38.5020998563
        ,40.1393375624
        ,41.0666243477
        ,36.393636275
        ,38.5632028326
        ,40.6227568239
        ,41.1521879097
        ,41.6800684404
        ,39.0624777172
        ,40.7342862612
        ,40.786167145
        ,41.5943818036
        ,41.2798844609
        ,41.4958048947
        ,41.3186791849
        ,37.9849524764
        ,40.5403068912
        ,40.3106668232
        ,38.5741014523
        ,39.6049215973
        ,41.3392339457
        ,39.2819503516
        ,40.8215974583
        ,40.2129770735
        ,39.1872750155
        ,38.1877012201
        ,41.0542374079
        ,39.9663270466
        ,39.3027549618
        ,40.5699103652
        ,39.6471841809
        ,40.7841092312
        ,40.0333916493
        ,41.3167357808
        ,38.7531726515
        ,39.1375855648
        ,39.3037024249
        ,40.8596631914
        ,38.8075567467
        ,40.5592893806
        ,38.7165439115
        ,39.5156068266
        ,39.4205055796
        ,40.7874215484
        ,40.6700406667
        ,38.6596242466
        ,39.3773509139
        ,40.0315745848
        ,40.2148561368
        ,40.2851903681
        ,40.2680210961
        ,37.9819349544
        ,37.1535798244
        ,40.8489855129
        ,40.7532140999
        ,39.3285968248
        ,41.2409920908
        ,41.396805305
        ,41.0563346847
        ,39.3150656135
        ,41.6178771066
        ,39.5466952992
        ,40.8576976883
        ,38.7723677225
        ,40.2077951229
        ,38.7645868009
        ,39.3526221585
        ,41.8584359733
        ,40.0988547901
        ,40.7555023446
        ,40.5530673076
        ,41.5147227504
        ,38.8646925888
        ,41.7503175801
        ,41.0341628563
        ,39.4855106255
        ,40.7832381081
        ,40.9195402038
        ,39.8366280383
        ,37.3759156237
        ,41.4006012087
        ,40.6017587956
        ,41.045986724
        ,40.6829741028
        ,40.5154942024
        ,40.9271885736
        ,39.7779949203
        ,40.5175395642
        ,38.7102646797
        ,39.8151415007
        ,41.0285113011
        ,41.210916386
        ,37.8397856085
        ,41.9224021245
        ,39.8812126015
        ,40.9849963429
        ,39.9909774034
        ,39.7491170234
        ,38.2966544602
        ,41.0429300945
        ,38.6660439788
        ,38.0070616987
        ,39.9944125695
        ,39.9334159336
        ,38.9203810159
        ,40.2168645316
        ,41.1205696642
        ,40.2294909454
        ,39.9496020384
        ,38.1541733477
        ,41.4455391404
        ,39.9869791346
        ,41.8848542117
        ,42.1513901027
        ,40.6523430162
        ,39.7133625937
        ,40.0983712623
        ,41.1643837956
        ,41.7971037106
        ,38.7000698518
        ,39.1892893551
        ,40.5264639914
        ,41.9150400049
        ,41.3091071435
        ,40.4975321919
        ,40.3796999961
        ,41.3250227805
        ,43.2414539718
        ,42.5033997022
        ,40.978473029
        ,42.0839004042
        ,42.4243368959
        ,40.7410932476
        ,42.0682788206
        ,42.3772235339
        ,42.8096709251
        ,42.3725071675
        ,43.5528793931
        ,42.4507674693
        ,42.0658249423
        ,41.8473340603
        ,43.2099933879
        ,41.2583245282
        ,42.984995669
        ,41.7782581789
        ,41.9753335602
        ,43.2038030699
        ,40.2161567859
        ,44.3775825486
        ,44.1658958719
        ,39.3060887655
        ,41.9606135946
        ,44.5058457473
        ,43.3165540203
        ,41.8021270066
        ,42.2798419395
        ,43.163757796
        ,42.1349159624
        ,41.7890062516
        ,42.6854322147
        ,42.2304967413
        ,42.4252888348
        ,42.1742925555
        ,42.5543829775
        ,42.3331570972
        ,42.3212910372
        ,42.9963723488
        ,41.8196468352
        ,42.8441573957
        ,41.2107609985
        ,42.8072073246
        ,41.5784932046
        ,44.1117949333
        ,43.3110517598
        ,42.4736471276
        ,43.5101906132
        ,44.207875513
        ,44.0859183588
        ,43.7942852172
        ,42.7280926842
        ,43.5871786373
        ,41.0979300094
        ,41.4943125288
        ,43.6655368654
        ,43.3737484585
        ,43.3697792983
        ,43.6994307822
        ,43.2887235894
        ,43.5326743178
        ,39.4462680572
        ,40.8362001457
        ,39.7205596334
        ,40.7878025381
        ,40.3913025086
        ,42.4924155173
        ,43.2578441337
        ,42.7777247214
        ,42.1996725697
        ,43.1807709647
        ,42.3326264116
        ,42.012529312
        ,41.3484969437
        ,43.6120619833
        ,42.4390506937
        ,43.9375562435
        ,42.5493325602
        ,43.9018471523
        ,43.6192033448
        ,43.8124852428
        ,43.560407391
        ,41.7962473711
        ,42.9023716786
        ,44.0535603584
        ,43.3156514414
        ,41.8744251286
        ,42.5769689992
        ,43.5871216709
        ,43.2116658024
        ,43.0819715543
        ,40.6239072462
        ,42.5758707909
        ,41.0697589526
        ,41.3543412802
        ,44.2645401128
        ,42.9845810521
        ,41.3549032931
        ,44.1893249256
        ,43.0248980987
        ,43.3070687497
        ,42.1421923658
        ,43.0397842919
        ,42.1985723808
        ,41.6294017233
        ,42.4289650526
        ,42.5524984575
        ,43.9666494559
        ,40.8511146401
        ,40.8220972223
        ,42.2963775812
        ,43.5171093106
        ,42.2096742771
        ,43.4274096531
        ,43.841972097
        ,43.6787360537
        ,43.3334915384
        ,42.2898806519
        ,41.4214717885
        ,44.0371993632
        ,43.7160273544
        ,41.8054832197
        ,42.4095309248
        ,42.5652698448
        ,43.6621900525
        ,43.4319714347
        ,42.0972776673
        ,41.6640301668
        ,43.7102538752
        ,43.4861035947
        ,42.6983720673
        ,42.2953525972
        ,43.3737778575
        ,43.9924930818
        ,41.3217710555
        ,44.3268198751
        ,40.6244185037
        ,43.0199766985
        ,43.0987782021
        ,41.4503101861
        ,43.1874760163
        ,42.7733707088
        ,42.787445377
        ,43.4040974169
        ,42.0541814484
        ,42.0196027616
        ,43.0704958563
        ,44.2406074622
        ,42.365523828
        ,43.2254104475
        ,41.7735214129
        ,41.3426796309
        ,42.4054033796
        ,43.0843340799
        ,41.7052599054
        ,39.9080744849
        ,41.2288219413
        ,41.8946716216
        ,41.5682803687
        ,41.9353580779
        ,41.1833648256
        ,41.3961782167
        ,39.4107379798
        ,41.2676756269
        ,42.3702690971
        ,40.2574128205
        ,41.4318861711
        ,42.8246855768
        ,42.3865265731
        ,42.8351146004
        ,41.9059782715
        ,42.0884631073
        ,42.2377984208
        ,42.885843734
        ,41.3581075779
        ,41.4232510759
        ,42.7610577503
        ,43.0710915246
        ,41.7458788277
        ,42.5245932467
        ,42.5627722403
        ,40.0720006034
        ,42.210722156
        ,41.6620032351
        ,41.2957407355
        ,41.0820238262
        ,41.8887223665
        ,42.1362463861
        ,41.4864269384
        ,42.1528623126
        ,43.0153077416
        ,42.8464372814
        ,42.259591209
        ,42.657813727
        ,42.1973742088
        ,41.5970773005
        ,40.1001985435
        ,41.1905064352
        ,40.5328864644
        ,42.3924561345
        ,40.7263236177
        ,39.6956112205
        ,43.1797675167
        ,42.9107934322
        ,41.9343745013
        ,41.5401633967
        ,41.2113198625
        ,43.069622149
        ,40.3605248738
        ,41.2482461469
        ,43.6987176723
        ,41.9498405771
        ,41.0072769026
        ,42.8981977617
        ,42.6951558939
        ,41.4880607389
        ,40.876793648
        ,42.3751743659
        ,40.4037818654
        ,40.8493617042
        ,42.2037743798
        ,43.4540256112
        ,42.6305005668
        ,40.5752670828
        ,42.0650099983
        ,43.6483165605
        ,44.1666838278
        ,43.7075622463
        ,42.3383621215
        ,44.2946061494
        ,44.2601276384
        ,43.6206227915
        ,45.0317825245
        ,44.2680083865
        ,44.481706307
        ,43.3016559562
        ,44.8354256258
        ,42.143982369
        ,44.0808776922
        ,44.438447316
        ,44.9361444895
        ,44.2836832301
        ,43.4104143734
        ,44.835860554
        ,42.387525323
        ,42.1139203511
        ,43.8364762181
        ,44.8687556563
        ,42.9259149014
        ,44.6420998553
        ,40.5426584948
        ,43.1000774659
        ,44.2011896159
        ,45.1191615825
        ,44.8535903247
        ,43.1550267844
        ,42.9104214351
        ,45.0024532701
        ,43.5102466754
        ,43.5500467052
        ,44.2912529866
        ,44.5764033273
        ,42.5302511229
        ,42.9914166016
        ,42.388163454
        ,44.2516558833]

import matplotlib.pyplot as plt
plt.figure()
plt.scatter(zpicks1,mag)
plt.title('zpicks1')

plt.figure()
plt.scatter(zpicks2,mag)
plt.title('zpicks2')


zpicks1, mag = (list(t) for t in zip(*sorted(zip(zpicks1, mag))))

plt.figure()
plt.plot(zpicks1,mag)
plt.title('zpicks1 sorted')

mag = np.asarray(mag)

output1 = mag, zpicks1

# Relative path of output folder.
save_path = './data/'+'sorted1'
    
import pickle
pickle.dump(output1, open(save_path, 'wb'))


import results
mag1, zpicks1 = results.load('./data', 'sorted1')


plt.figure()
plt.title('model: sorted1'
          +'\n Evolution of magnitude with redshift')
#data = plt.errorbar(zpicks1, mag2, yerr=0.7, fmt='.', alpha=0.3)
best_fit = plt.scatter(zpicks1, mag1, marker='.', lw='1', c='xkcd:tomato')
plt.ylabel('magnitude')
plt.xlabel('z')
plt.show(block=False)
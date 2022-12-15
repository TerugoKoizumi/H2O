# 水分子の全エネルギー収束について

まずカットオフエネルギーについて収束を求めます。

```
import numpy as np
from ase.parallel import paropen
from gpaw import GPAW, PW, FermiDirac, Davidson
from ase.build import molecule

with paropen('a.dat', 'w') as fd:
    for ecut in range(200, 1201, 40):
        a = 15
        h2o = molecule('H2O')
        h2o.set_cell((a, a, a))
        h2o.center()
        calc = GPAW(mode=PW(ecut),
                    nbands=10,
                    kpts={'gamma':True},
                    xc='PBE',
                    occupations=FermiDirac(0.0, fixmagmom=True),
                    eigensolver='rmm-diis',
                    convergence={'energy': 2.0e-11},
                    txt=f'h2o-{ecut}.txt')
        h2o.calc = calc
        e = h2o.get_potential_energy()
        calc.write(f'h2o-{ecut}.gpw')

        print(ecut, e, file=fd)
```

このスクリプトを利用し、それぞれのecutにおけるエネルギーをプロットした図が以下になります。

![グラフ](/スクリーンショット%202022-12-15%2016.44.42.jpg)

平面波の数は`Number of coefficients`で合ってますでしょうか。
ecutそれぞれの出力から得られた結果は以下になります。

```
h2o-1000.txt:  Number of coefficients: 242361 (reduced to 121181)
h2o-1040.txt:  Number of coefficients: 257031 (reduced to 128516)
h2o-1080.txt:  Number of coefficients: 271985 (reduced to 135993)
h2o-1120.txt:  Number of coefficients: 287369 (reduced to 143685)
h2o-1160.txt:  Number of coefficients: 302659 (reduced to 151330)
h2o-1200.txt:  Number of coefficients: 318721 (reduced to 159361)
h2o-200.txt:  Number of coefficients: 21823 (reduced to 10912)
h2o-240.txt:  Number of coefficients: 28425 (reduced to 14213)
h2o-280.txt:  Number of coefficients: 35825 (reduced to 17913)
h2o-320.txt:  Number of coefficients: 43819 (reduced to 21910)
h2o-360.txt:  Number of coefficients: 52299 (reduced to 26150)
h2o-400.txt:  Number of coefficients: 61445 (reduced to 30723)
h2o-440.txt:  Number of coefficients: 70751 (reduced to 35376)
h2o-480.txt:  Number of coefficients: 80725 (reduced to 40363)
h2o-520.txt:  Number of coefficients: 90879 (reduced to 45440)
h2o-560.txt:  Number of coefficients: 101505 (reduced to 50753)
h2o-600.txt:  Number of coefficients: 112451 (reduced to 56226)
h2o-640.txt:  Number of coefficients: 124097 (reduced to 62049)
h2o-680.txt:  Number of coefficients: 135883 (reduced to 67942)
h2o-720.txt:  Number of coefficients: 148381 (reduced to 74191)
h2o-760.txt:  Number of coefficients: 160467 (reduced to 80234)
h2o-800.txt:  Number of coefficients: 173541 (reduced to 86771)
h2o-840.txt:  Number of coefficients: 186623 (reduced to 93312)
h2o-880.txt:  Number of coefficients: 200237 (reduced to 100119)
h2o-920.txt:  Number of coefficients: 213583 (reduced to 106792)
h2o-960.txt:  Number of coefficients: 228007 (reduced to 114004)
```

この結果よりセルのサイズはカットオフを640eVにして、サイズごとにエネルギーを求めました。



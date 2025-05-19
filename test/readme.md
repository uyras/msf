# Файлы для тестирования программы

Запускать так:
```
msf -f cyrrhus.mfsys -o cyrrhus.dat  --qmin=-7 --qmax=7 --qstep=0.4
```
После завершения вычислений эти данные нужно визуализировать.
Примерный код на python для визуализации:

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.colors as col

def sq(arr):
    return arr.reshape((-1,int(np.sqrt(arr.shape[0]))))

fname="cyrrhus.dat"
arr = np.loadtxt(fname)
qx = sq(arr[:,1])
qy = sq(arr[:,2])
re = sq(arr[:,3])
im = sq(arr[:,4])

# draw real part of MSF
plt.figure(dpi = 300)
plt.gca().set_aspect('equal')
plt.xlabel("$q_x$",labelpad=-8, x=0.6) 
plt.ylabel("$q_y$",labelpad=-8, y=0.6)
plt.pcolormesh(qx,qy,re,cmap='jet')
cbr = plt.colorbar(pad=0)
plt.savefig("cyrrhus_re.png", bbox_inches='tight')

# draw imaginary part of MSF
plt.figure(dpi = 300)
plt.gca().set_aspect('equal')
plt.xlabel("$q_x$",labelpad=-8, x=0.6) 
plt.ylabel("$q_y$",labelpad=-8, y=0.6)
plt.pcolormesh(qx,qy,im,cmap='jet')
cbr = plt.colorbar(pad=0)
plt.savefig("cyrrhus_im.png", bbox_inches='tight')
```

## todo

[ ] добавить в тест txt-файл с системой
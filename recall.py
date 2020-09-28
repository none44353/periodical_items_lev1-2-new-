import matplotlib.pyplot as plt
import numpy as np

CM=[0.335, 0.382, 0.432, 0.479, 0.525, 0.568, 0.6, 0.649, 0.69, 0.746, 0.795, 0.82, 0.882, 0.916, 0.942, 0.962, 0.986, 0.995, 1, 1, 1]
CU=[0.548, 0.505, 0.543, 0.572, 0.641, 0.683, 0.742, 0.786, 0.833, 0.878, 0.921, 0.959, 0.985, 1, 1, 1, 1, 1, 1, 1, 1]

plt.plot(CM,color='m',linewidth=1,linestyle='-',label='CM T=3')
plt.plot(CU,color='m',linewidth=1,linestyle='--',label='CU T=3')

plt.legend(loc=2)#标签展示位置，数字代表标签具位置

plt.title('recall for top1000')

plt.xlabel('memory/10KB')
plt.xticks(np.arange(21),('5','6','7','8','9','10','11','12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25'))#更改图表X标签为制定内容

plt.ylim(0,1)#Y轴标签范围

plt.show()
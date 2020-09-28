import matplotlib.pyplot as plt
import numpy as np

CMt1=[0.561888, 0.271383, 0.169245, 0.147825, 0.100676, 0.08581, 0.090209, 0.082289, 0.047986, 0.059673, 0.040393, 0.045131, 0.030214, 0.038099, 0.032636, 0.040724, 0.0276, 0.026083, 0.040394, 0.034317]
CMt2=[0.807827, 0.23631, 0.106896, 0.064582, 0.037389, 0.028304, 0.01814, 0.017712, 0.012509, 0.010158, 0.008389, 0.008141, 0.007254, 0.007943, 0.005095, 0.004332, 0.003669, 0.002798, 0.003226, 0.002749]
CMt3=[1.151924, 0.229616, 0.093756, 0.043369, 0.028343, 0.018734, 0.012654, 0.008738, 0.007238, 0.006602, 0.003988, 0.003844, 0.003317, 0.002286, 0.002671, 0.002296, 0.001744, 0.001459, 0.001341, 0.001127]
CMt4=[1.688441, 0.267915, 0.086018, 0.042228, 0.02278, 0.015969, 0.011048, 0.008137, 0.006016, 0.004867, 0.003888, 0.003231, 0.003073, 0.002118, 0.001829, 0.001663, 0.001407, 0.001143, 0.001063, 0.001118]
CMt5=[2.39399, 0.316061, 0.093186, 0.042947, 0.025577, 0.017284, 0.011226, 0.009077, 0.007549, 0.005593, 0.004698, 0.003534, 0.003007, 0.002645, 0.002062, 0.001768, 0.001693, 0.001432, 0.001107, 0.001017]
CMt6=[3.323052,0.391474,0.108911,0.050183,0.028866,0.019552,0.014607,0.010703,0.008402,0.006577,0.005376,0.004495,0.003437,0.003383,0.002743,0.002186,0.002093,0.001607,0.001464,0.00126]

CUt1=[0.561888, 0.271383, 0.169245, 0.147825, 0.100676, 0.08581, 0.090209, 0.082289, 0.047986, 0.059673, 0.040393, 0.045131, 0.030214, 0.038099, 0.032636, 0.040724, 0.0276, 0.026083, 0.040394, 0.034317]
CUt2=[0.247958, 0.07061, 0.032976, 0.020274, 0.011737, 0.00904, 0.006317, 0.005822, 0.004303, 0.003272, 0.002893, 0.002926, 0.002423, 0.002529, 0.001733, 0.001553, 0.001328, 0.000997, 0.001105, 0.001076]
CUt3=[0.265828, 0.05672, 0.025676, 0.01016, 0.00644, 0.004056, 0.002613, 0.002008, 0.001505, 0.001343, 0.000821, 0.000604, 0.000634, 0.000378, 0.000344, 0.000309, 0.000296, 0.000319, 0.000203, 0.000163]
CUt4=[0.353144, 0.060986, 0.019283, 0.008276, 0.003419, 0.002198, 0.001363, 0.00083, 0.000648, 0.000445, 0.000332, 0.000205, 0.000189, 0.000155, 0.000069, 0.000105, 0.000088, 0.000039, 0.000035, 0.000036]
CUt5=[0.514279,0.065767,0.017742,0.005768,0.002816,0.001329,0.00086,0.000605,0.00046,0.000218,0.000219,0.000143,0.000166,0.000112,0.00007,0.000066,0.000051,0.000058,0.00003,0.00004]
CUt6=[0.781064,0.077542,0.017268,0.005433,0.002265,0.00109,0.00077,0.000443,0.000296,0.000159,0.000221,0.000118,0.000146,0.000147,0.000061,0.000058,0.000052,0.000047,0.000067,0.000023]

HGh1=[0.330755,0.156214,0.097622,0.070932,0.058547,0.044407,0.040237,0.03537,0.029261,0.024362,0.02882,0.027162,0.01829,0.019776,0.021326,0.015808,0.012787,0.010049,0.021263,0.012581]
HGh2=[0.339709,0.087785,0.045592,0.019992,0.018515,0.01102,0.008121,0.005217,0.005567,0.006721,0.00385,0.002351,0.0024,0.002744,0.00174,0.000591,0.000823,0.001704,0.000411,0.00237]
HGh3=[0.33639,0.082595,0.031681,0.017679,0.010813,0.005202,0.003786,0.002375,0.003247,0.001489,0.001284,0.00102,0.001204,0.000668,0.000593,0.000456,0.000456,0.000352,0.000178,0.000208]
HGh4=[0.343339,0.075186,0.029244,0.013114,0.008525,0.005308,0.004004,0.002263,0.001858,0.001332,0.001196,0.000887,0.000951,0.00081,0.000466,0.000402,0.000305,0.000275,0.000184,0.000212]
HGh5=[0.350815,0.075474,0.02847,0.013429,0.008533,0.005041,0.003764,0.002226,0.001842,0.001373,0.001178,0.000888,0.000907,0.000591,0.000452,0.000363,0.000309,0.00028,0.000164,0.000196]
HGh6=[0.33931,0.073807,0.028091,0.013513,0.009039,0.004759,0.003371,0.002188,0.001836,0.001439,0.001175,0.0009,0.000901,0.000547,0.000437,0.000356,0.000322,0.000254,0.00017,0.000201]

plt.plot(CMt1,color='b',linewidth=1,linestyle='-',label='CM/CU T=1')
'''
plt.plot(CMt2,color='g',linewidth=1,linestyle='-',label='CM T=2')
plt.plot(CMt3,color='m',linewidth=1,linestyle='-',label='CM T=3')
plt.plot(CMt4,color='y',linewidth=1,linestyle='-',label='CM T=4')
plt.plot(CMt5,color='r',linewidth=1,linestyle='-',label='CM T=5')
plt.plot(CMt6,color='c',linewidth=1,linestyle='-',label='CM T=6')
'''
plt.plot(CUt2,color='g',linewidth=1,linestyle='--',label='CU T=2')
plt.plot(CUt3,color='m',linewidth=1,linestyle='--',label='CU T=3')
plt.plot(CUt4,color='y',linewidth=1,linestyle='--',label='CU T=4')
plt.plot(CUt5,color='r',linewidth=1,linestyle='--',label='CU T=5')
plt.plot(CUt6,color='c',linewidth=1,linestyle='--',label='CU T=6')

plt.plot(HGh1,color='b',linewidth=1,linestyle='-.',label='HG H=1')
plt.plot(HGh2,color='g',linewidth=1,linestyle='-.',label='HG H=2')
plt.plot(HGh3,color='m',linewidth=1,linestyle='-.',label='HG H=3')
plt.plot(HGh4,color='y',linewidth=1,linestyle='-.',label='HG H=4')
plt.plot(HGh5,color='r',linewidth=1,linestyle='-.',label='HG H=5')
plt.plot(HGh6,color='c',linewidth=1,linestyle='-.',label='HG H=6')

plt.legend(loc=2)#标签展示位置，数字代表标签具位置

plt.title('ARE for topK items')

plt.xlabel('memory/100KB')
plt.xticks(np.arange(20),('1','2','3','4','5','6','7','8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20'))#更改图表X标签为制定内容

plt.ylim(0,2)#Y轴标签范围

plt.show()
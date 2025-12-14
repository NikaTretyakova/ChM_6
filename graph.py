import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

# Чтение данных из файлов
data1 = np.loadtxt('signal_data.txt')
x = data1[:, 0]
y_original = data1[:, 1]
y_filtered = data1[:, 2]

# Интерполяция сплайнами
x_smooth = np.linspace(x.min(), x.max(), 300)
spline_original = make_interp_spline(x, y_original, k=3)
spline_filtered = make_interp_spline(x, y_filtered, k=3)
y_original_smooth = spline_original(x_smooth)
y_filtered_smooth = spline_filtered(x_smooth)

# Построение графика
plt.figure(figsize=(12, 6))
plt.plot(x_smooth, y_original_smooth, 'b-', linewidth=2, alpha=0.7, label='Исходный сигнал')
plt.plot(x_smooth, y_filtered_smooth, 'r-', linewidth=2, label='После фильтрации')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Фильтрация зашумленного сигнала')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('filtering_result.png', dpi=150)
plt.show()


data2 = np.loadtxt('signal2_data.txt')
x2 = data2[:, 0]
y2_original = data2[:, 1]

x2_smooth = np.linspace(x2.min(), x2.max(), 300)
spline2 = make_interp_spline(x2, y2_original, k=3)
y2_smooth = spline2(x2_smooth)

plt.figure(figsize=(12, 6))
plt.plot(x2_smooth, y2_smooth, 'g-', linewidth=2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Сигнал с разрывами')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('signal2_plot.png', dpi=150)
plt.show()
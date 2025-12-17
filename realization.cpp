#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <direct.h>
#include <iomanip> 

using namespace std;
using Complex = complex<double>;
using Signal = vector<Complex>;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Реализация DFT и IDFT
Signal dft(const Signal& x) {
    int N = static_cast<int>(x.size());
    Signal X(N);
    for (int k = 0; k < N; ++k) {
        X[k] = Complex(0, 0);
        for (int n = 0; n < N; ++n) {
            double angle = -2.0 * M_PI * k * n / N;
            X[k] += x[n] * Complex(cos(angle), sin(angle));
        }
    }
    return X;
}

Signal idft(const Signal& X) {
    int N = static_cast<int>(X.size());
    Signal x(N);
    for (int n = 0; n < N; ++n) {
        x[n] = Complex(0, 0);
        for (int k = 0; k < N; ++k) {
            double angle = 2.0 * M_PI * k * n / N;
            x[n] += X[k] * Complex(cos(angle), sin(angle));
        }
        x[n] /= static_cast<double>(N);
    }
    return x;
}

// Вспомогательная рекурсивная функция FFT
void computeFFTRecursive(Signal& x, int start, int step, int N) {
    if (N == 1) return;

    int halfN = N / 2;

    computeFFTRecursive(x, start, step * 2, halfN);
    computeFFTRecursive(x, start + step, step * 2, halfN);

    double angleIncrement = -2.0 * M_PI / N;
    Complex rotationStep(cos(angleIncrement), sin(angleIncrement));
    Complex currentRotation(1.0, 0.0);

    for (int k = 0; k < halfN; ++k) {
        int evenIndex = start + step * (2 * k);
        int oddIndex = start + step * (2 * k + 1);

        Complex evenValue = x[evenIndex];
        Complex oddValue = x[oddIndex];

        Complex rotatedOdd = currentRotation * oddValue;

        x[evenIndex] = evenValue + rotatedOdd;
        x[oddIndex] = evenValue - rotatedOdd;

        currentRotation *= rotationStep;
    }
}

// Реализация FFT 
Signal fft(const Signal& x) {
    int N = static_cast<int>(x.size());
    if (N <= 1) return x;

    if ((N & (N - 1)) != 0) {
        cerr << "Размер сигнала должен быть степенью двойки для FFT" << endl;
        return Signal();
    }

    Signal result = x;  
    computeFFTRecursive(result, 0, 1, N);
    return result;
}

Signal ifft(const Signal& X) {
    int N = static_cast<int>(X.size());
    Signal conjX;
    for (const auto& val : X) conjX.push_back(conj(val));

    Signal x = fft(conjX);

    for (auto& val : x) {
        val = conj(val) / static_cast<double>(N);
    }
    return x;
}

// Генерация зашумлённого сигнала
Signal generateNoisySignal(int N, double A, double omega1, double phi, double B, double omega2) {
    Signal z(N);
    for (int j = 0; j < N; ++j) {
        double value = A * cos(2.0 * M_PI * omega1 * j / N + phi) +
            B * cos(2.0 * M_PI * omega2 * j / N);
        z[j] = Complex(value, 0);
    }
    return z;
}

// Генерация сигнала 
Signal generateSignal2(int N, double A, double B, double omega2) {
    Signal z(N);
    int N4 = N / 4;
    int N2 = N / 2;
    int N34 = 3 * N / 4;

    for (int j = 0; j < N; ++j) {
        double value = 0.0;
        if (j >= N4 && j <= N2) {
            value = A + B * cos(2.0 * M_PI * omega2 * j / N);
        }
        else if (j > N34) {
            value = A + B * cos(2.0 * M_PI * omega2 * j / N);
        }
        z[j] = Complex(value, 0);
    }
    return z;
}

// Вывод таблицы спектров
void printSpectrumTable(const Signal& z_original, const Signal& Z, const string& title) {
    int N = static_cast<int>(Z.size());
    cout << "\n" << title << ":\n";
    cout << setw(3) << "m" << " | "
        << setw(12) << "Re(z)" << " | "
        << setw(12) << "Re(Z)^" << " | "
        << setw(12) << "Im(Z)^" << " | "
        << setw(15) << "Амплитуда" << " | "
        << setw(12) << "Фаза" << endl;

    cout << scientific << setprecision(6);

    for (int m = 0; m < N; ++m) {
        double amp = abs(Z[m]);
        if (amp > 1e-10) {
            double phase = arg(Z[m]);

            cout << setw(3) << m << " | "
                << setw(12) << z_original[m].real() << " | "
                << setw(12) << Z[m].real() << " | "
                << setw(12) << Z[m].imag() << " | "
                << setw(15) << amp << " | "
                << setw(12) << phase << endl;
        }
    }
    cout << defaultfloat;
}

// Обнуление высокочастотных шумовых компонентов
Signal removeHighFrequencyNoise(const Signal& Z) {
    int N = static_cast<int>(Z.size());
    Signal Z_filtered = Z;

    int m_noise1 = 185;
    int m_noise2 = N - 185; // 1024 - 185 = 839

    Z_filtered[m_noise1] = Complex(0, 0);
    Z_filtered[m_noise2] = Complex(0, 0);

    int bandwidth = 2; // ±2 отсчета
    for (int delta = -bandwidth; delta <= bandwidth; ++delta) {
        int m1 = m_noise1 + delta;
        int m2 = m_noise2 + delta;

        if (m1 >= 0 && m1 < N) Z_filtered[m1] = Complex(0, 0);
        if (m2 >= 0 && m2 < N) Z_filtered[m2] = Complex(0, 0);
    }

    return Z_filtered;
}

// Сохранение данных для построения графиков
void saveToFile(const string& filename, const Signal& signal1, const Signal& signal2 = Signal()) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Ошибка открытия файла: " << filename << endl;
        return;
    }

    int N = static_cast<int>(signal1.size());
    for (int i = 0; i < N; ++i) {
        file << i << " " << signal1[i].real();
        if (!signal2.empty() && static_cast<int>(signal2.size()) == N) {
            file << " " << signal2[i].real();
        }
        file << endl;
    }
    file.close();
}

int main() {
    setlocale(0, "");

    int n = 10;
    int N = static_cast<int>(pow(2, n));

    double A = 2.15;
    double omega1 = 2.0;
    double phi = M_PI / 4;
    double B = 0.18;
    double omega2 = 185.0;

    // Путь для сохранения файлов 
    string desktopPath = "C:\\Users\\Veronika\\Desktop\\График\\";
    _mkdir(desktopPath.c_str());

    // Генерация сигнала
    Signal z = generateNoisySignal(N, A, omega1, phi, B, omega2);

    // DFT и FFT с замером времени
    cout << "   DFT... ";
    auto start = chrono::high_resolution_clock::now();
    Signal Z_dft = dft(z);
    auto end = chrono::high_resolution_clock::now();
    auto dft_time = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << dft_time << " мкс\n";

    cout << "   FFT... ";
    start = chrono::high_resolution_clock::now();
    Signal Z_fft = fft(z);
    end = chrono::high_resolution_clock::now();
    auto fft_time = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << fft_time << " мкс\n";

    cout << "\n   Сравнение производительности:\n";
    cout << "   DFT: " << dft_time << " микросекунд\n";
    cout << "   FFT: " << fft_time << " микросекунд\n";
    if (fft_time > 0) {
        double speedup = static_cast<double>(dft_time) / fft_time;
        cout << "   Ускорение FFT: " << speedup << " раз\n";
    }

    // Таблицы спектров
    printSpectrumTable(z, Z_dft, "Спектр DFT");
    printSpectrumTable(z, Z_fft, "Спектр FFT");

    // Проверка, что DFT и FFT дают одинаковые результаты
    double max_diff = 0.0;
    for (int i = 0; i < N; ++i) {
        double diff = abs(Z_dft[i] - Z_fft[i]);
        if (diff > max_diff) max_diff = diff;
    }
    cout << "Максимальная разница между DFT и FFT: " << max_diff << endl;
    if (max_diff < 1e-10) {
        cout << "DFT и FFT дают идентичные результаты" << endl;
    }

    // Обнуление конкретных шумовых компонентов
    Signal Z_filtered = removeHighFrequencyNoise(Z_fft);

    // Табличка спектра после фильтрации
    printSpectrumTable(z, Z_filtered, "Спектр после фильтрации");

    Signal z_filtered = ifft(Z_filtered);

    // Сохранение спектров для сравнения
    saveToFile(desktopPath + "spectrum_before.txt", Z_fft);
    saveToFile(desktopPath + "spectrum_after.txt", Z_filtered);

    // Сохранение сигналов
    string file1 = desktopPath + "signal_data.txt";
    saveToFile(file1, z, z_filtered);

    // Генерация сигнала с разрывами
    Signal z2 = generateSignal2(N, A, B, omega2);

    // Сохранение данных для второй задачи
    string file2 = desktopPath + "signal2_data.txt";
    saveToFile(file2, z2);

    return 0;
}

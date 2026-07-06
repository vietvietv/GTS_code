# Lời giải tìm giá trị riêng trội

> Các phép tính bên trong dùng giá trị gốc. Ma trận trình bày 4 chữ số thập phân; các số liệu khác 7 chữ số thập phân.

## Bài toán

Cho ma trận

$$
A=\begin{pmatrix}
2.0000 & 1.0000 & 0.0000\\
1.0000 & 2.0000 & 1.0000\\
0.0000 & 1.0000 & 2.0000
\end{pmatrix}
$$

Dùng phương pháp lũy thừa để tìm giá trị riêng trội và véc tơ riêng tương ứng.

Các tham số sử dụng:

- Sai số: \(\varepsilon=0.0000001\)
- Số bước lặp tối đa: \(N=30\)
- Chế độ xuất lời giải: `exam`

## Phương pháp lũy thừa

Ta chọn véc tơ ban đầu

$$
x^{(0)}=\begin{pmatrix}
1.0000000\\
1.0000000\\
1.0000000
\end{pmatrix}
$$

Tại mỗi bước, tính

$$
y^{(k)}=Ax^{(k)}.
$$

Chọn chỉ số \(s\) sao cho

$$
|y_s^{(k)}|=\max_i |y_i^{(k)}|,
$$

rồi chuẩn hóa

$$
x^{(k+1)}=\frac{y^{(k)}}{y_s^{(k)}}.
$$

| \(k\) | \(y^{(k)}=Ax^{(k)}\) | Chuẩn hóa | \(x^{(k+1)}\) | \(\lambda^{(k)}\) | Sai số |
|---:|---|---|---|---:|---:|
| 1 | \((3.0000000, 4.0000000, 3.0000000)^T\) | \(s=2,\ y_s=4.0000000\) | \((0.7500000, 1.0000000, 0.7500000)^T\) | 4.0000000 | 0.8660254 |
| 2 | \((2.5000000, 3.5000000, 2.5000000)^T\) | \(s=2,\ y_s=3.5000000\) | \((0.7142857, 1.0000000, 0.7142857)^T\) | 3.5000000 | 0.1237179 |
| 3 | \((2.4285714, 3.4285714, 2.4285714)^T\) | \(s=2,\ y_s=3.4285714\) | \((0.7083333, 1.0000000, 0.7083333)^T\) | 3.4285714 | 0.0206197 |
| 4 | \((2.4166667, 3.4166667, 2.4166667)^T\) | \(s=2,\ y_s=3.4166667\) | \((0.7073171, 1.0000000, 0.7073171)^T\) | 3.4166667 | 0.0035204 |
| 5 | \((2.4146341, 3.4146341, 2.4146341)^T\) | \(s=2,\ y_s=3.4146341\) | \((0.7071429, 1.0000000, 0.7071429)^T\) | 3.4146341 | 0.0006035 |
| 6 | \((2.4142857, 3.4142857, 2.4142857)^T\) | \(s=2,\ y_s=3.4142857\) | \((0.7071130, 1.0000000, 0.7071130)^T\) | 3.4142857 | 0.0001035 |
| 7 | \((2.4142259, 3.4142259, 2.4142259)^T\) | \(s=2,\ y_s=3.4142259\) | \((0.7071078, 1.0000000, 0.7071078)^T\) | 3.4142259 | 0.0000178 |
| 8 | \((2.4142157, 3.4142157, 2.4142157)^T\) | \(s=2,\ y_s=3.4142157\) | \((0.7071070, 1.0000000, 0.7071070)^T\) | 3.4142157 | 0.0000030 |
| 9 | \((2.4142139, 3.4142139, 2.4142139)^T\) | \(s=2,\ y_s=3.4142139\) | \((0.7071068, 1.0000000, 0.7071068)^T\) | 3.4142139 | 0.0000005 |
| 10 | \((2.4142136, 3.4142136, 2.4142136)^T\) | \(s=2,\ y_s=3.4142136\) | \((0.7071068, 1.0000000, 0.7071068)^T\) | 3.4142136 | 0.0000001 |
| 11 | \((2.4142136, 3.4142136, 2.4142136)^T\) | \(s=2,\ y_s=3.4142136\) | \((0.7071068, 1.0000000, 0.7071068)^T\) | 3.4142136 | 0.0000000 |

Dãy lặp hội tụ theo sai số đã chọn.

Suy ra giá trị riêng trội xấp xỉ

$$
\lambda_1\approx 3.4142136.
$$

Véc tơ riêng tương ứng có thể lấy là

$$
v_1\approx \begin{pmatrix}
0.7071068\\
1.0000000\\
0.7071068
\end{pmatrix}.
$$

## Kiểm tra trường hợp

Ta kiểm tra lần lượt theo thứ tự:

$$
\text{TH1} \rightarrow \text{TH2} \rightarrow \text{TH3}.
$$

### Kiểm tra TH1

TH1 xảy ra khi dãy lặp hội tụ trực tiếp, tức là \(x^{(k+1)}\approx x^{(k)}\) và \(\lambda^{(k)}\) ổn định.

Ở bước cuối, sai số kiểm tra xấp xỉ:

$$
E_1=0.0000000,\qquad \varepsilon_1=0.0000010.
$$

Vì \(E_1\le \varepsilon_1\), bài toán thuộc TH1.

## Phương pháp xuống thang

Sau khi tìm được \((\lambda_1,v_1)\), tìm véc tơ riêng trái \(w_1\) từ \(A^T w_1=\lambda_1w_1\), rồi đặt

$$
A_1=A-\frac{\lambda_1}{w_1^Tv_1}v_1w_1^T.
$$

### Xuống thang lần 1

Ta có

$$
\lambda_1\approx 3.4142136,\qquad v_1\approx \begin{pmatrix}
0.7071068\\
1.0000000\\
0.7071068
\end{pmatrix}
$$

Véc tơ riêng trái xấp xỉ:

$$
w_1\approx \begin{pmatrix}
0.7071068\\
1.0000000\\
0.7071068
\end{pmatrix}
$$

Mẫu số

$$
w_1^Tv_1=2.0000000.
$$

Ma trận sau khi xuống thang:

$$
A_1=\begin{pmatrix}
1.1464 & -0.2071 & -0.8536\\
-0.2071 & 0.2929 & -0.2071\\
-0.8536 & -0.2071 & 1.1464
\end{pmatrix}
$$

Tiếp tục áp dụng phương pháp lũy thừa cho \(A_1\), ta thu được giá trị riêng trội tiếp theo

$$
\lambda_2\approx 2.0000000.
$$

Véc tơ riêng tương ứng trong bài toán sau xuống thang là

$$
\tilde v_2\approx \begin{pmatrix}
-1.0000000\\
0.0000000\\
1.0000000
\end{pmatrix}.
$$

## Kết luận

Giá trị riêng trội tìm được là

$$
\lambda_1\approx 3.4142136.
$$

Véc tơ riêng tương ứng là

$$
v_1\approx \begin{pmatrix}
0.7071068\\
1.0000000\\
0.7071068
\end{pmatrix}.
$$

Các giá trị riêng tiếp theo tìm được bằng xuống thang:

- \(\lambda_2\approx 2.0000000\)

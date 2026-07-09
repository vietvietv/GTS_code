# Lời giải tìm giá trị riêng trội

> Các phép tính bên trong dùng giá trị gốc. Ma trận trình bày 4 chữ số thập phân; các số liệu khác 7 chữ số thập phân.

## Bài toán

Cho ma trận

$$
A=\begin{pmatrix}
86.0000 & 72.0000 & 149.0000\\
72.0000 & 194.0000 & 221.0000\\
149.0000 & 221.0000 & 330.0000
\end{pmatrix}
$$

Dùng phương pháp lũy thừa để tìm giá trị riêng trội và véc tơ riêng tương ứng.

Các tham số sử dụng:

- Sai số: \(\varepsilon=0.0000001\)
- Số bước lặp tối đa: \(N=300\)
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

Cột sai số trong bảng được tính bằng chuẩn của vector dư:

$$
r^{(k)}=\left\|Ax^{(k+1)}-\lambda^{(k)}x^{(k+1)}\right\|_2.
$$

| \(k\) | \(y^{(k)}=Ax^{(k)}\) | Chuẩn hóa | \(x^{(k+1)}\) | \(\lambda^{(k)}\) | Sai số |
|---:|---|---|---|---:|---:|
| 1 | \((307.0000000, 487.0000000, 700.0000000)^T\) | \(s=3,\ y_s=700.0000000\) | \((0.4385714, 0.6957143, 1.0000000)^T\) | 700.0000000 | 193.8783165 |
| 2 | \((236.8085714, 387.5457143, 549.1000000)^T\) | \(s=3,\ y_s=549.1000000\) | \((0.4312667, 0.7057835, 1.0000000)^T\) | 549.1000000 | 1.8274639 |
| 3 | \((236.9053516, 388.9732029, 550.2368967)^T\) | \(s=3,\ y_s=550.2368967\) | \((0.4305516, 0.7069195, 1.0000000)^T\) | 550.2368967 | 0.2231980 |
| 4 | \((236.9256392, 389.1420988, 550.3813956)^T\) | \(s=3,\ y_s=550.3813956\) | \((0.4304754, 0.7070408, 1.0000000)^T\) | 550.3813956 | 0.0238541 |
| 5 | \((236.9278196, 389.1601412, 550.3968467)^T\) | \(s=3,\ y_s=550.3968467\) | \((0.4304673, 0.7070537, 1.0000000)^T\) | 550.3968467 | 0.0025437 |
| 6 | \((236.9280522, 389.1620652, 550.3984944)^T\) | \(s=3,\ y_s=550.3984944\) | \((0.4304664, 0.7070551, 1.0000000)^T\) | 550.3984944 | 0.0002712 |
| 7 | \((236.9280770, 389.1622703, 550.3986701)^T\) | \(s=3,\ y_s=550.3986701\) | \((0.4304663, 0.7070553, 1.0000000)^T\) | 550.3986701 | 0.0000289 |
| 8 | \((236.9280796, 389.1622922, 550.3986888)^T\) | \(s=3,\ y_s=550.3986888\) | \((0.4304663, 0.7070553, 1.0000000)^T\) | 550.3986888 | 0.0000031 |
| 9 | \((236.9280799, 389.1622945, 550.3986908)^T\) | \(s=3,\ y_s=550.3986908\) | \((0.4304663, 0.7070553, 1.0000000)^T\) | 550.3986908 | 0.0000003 |
| 10 | \((236.9280799, 389.1622948, 550.3986910)^T\) | \(s=3,\ y_s=550.3986910\) | \((0.4304663, 0.7070553, 1.0000000)^T\) | 550.3986910 | 0.0000000 |
| 11 | \((236.9280799, 389.1622948, 550.3986910)^T\) | \(s=3,\ y_s=550.3986910\) | \((0.4304663, 0.7070553, 1.0000000)^T\) | 550.3986911 | 0.0000000 |

Dãy lặp hội tụ theo sai số đã chọn.

Suy ra giá trị riêng trội xấp xỉ

$$
\lambda_1\approx 550.3986911.
$$

Véc tơ riêng tương ứng có thể lấy là

$$
v_1\approx \begin{pmatrix}
0.3315964\\
0.5446581\\
0.7703190
\end{pmatrix}.
$$

## Kiểm tra trường hợp

Ta kiểm tra lần lượt theo thứ tự:

$$
\text{TH1} \rightarrow \text{TH2} \rightarrow \text{TH3}.
$$

### Kiểm tra TH1

TH1 xảy ra khi dãy lặp hội tụ trực tiếp, tức là \(x^{(k+1)}\approx x^{(k)}\) và \(\lambda^{(k)}\) ổn định.

Ta dùng đại lượng kiểm tra

$$
E_1=\max\left\{\frac{|\lambda^{(k)}-\lambda^{(k-1)}|}{\max(1,|\lambda^{(k)}|)},\frac{\left\|Ax^{(k+1)}-\lambda^{(k)}x^{(k+1)}\right\|_2}{\max(1,|\lambda^{(k)}|)}\right\}.
$$

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

Từ TH1, ta lấy trị riêng và véc tơ riêng phải:

$$
\lambda_1\approx 550.3986911,\qquad v_1\approx \begin{pmatrix}
0.3315964\\
0.5446581\\
0.7703190
\end{pmatrix}
$$

Tìm véc tơ riêng trái từ hệ \((A^T-\lambda I)w=0\):

$$
w_1\approx \begin{pmatrix}
0.3315964\\
0.5446581\\
0.7703190
\end{pmatrix}
$$

Mẫu số

$$
w_1^Tv_1=1.0000000.
$$

Ma trận sau khi xuống thang:

$$
A_1=\begin{pmatrix}
25.4803 & -27.4057 & 8.4089\\
-27.4057 & 30.7229 & -9.9256\\
8.4089 & -9.9256 & 3.3982
\end{pmatrix}
$$

Tiếp tục áp dụng phương pháp lũy thừa cho \(A_1\), ta thu được giá trị riêng trội tiếp theo

$$
\lambda_2\approx 58.6894315.
$$

Véc tơ riêng tương ứng trong bài toán sau xuống thang là

$$
\tilde v_2\approx \begin{pmatrix}
-0.6534213\\
0.7215563\\
-0.2289043
\end{pmatrix}.
$$

### Xuống thang lần 2

Từ TH1, ta lấy trị riêng và véc tơ riêng phải:

$$
\lambda_2\approx 58.6894315,\qquad v_2\approx \begin{pmatrix}
-0.6534213\\
0.7215563\\
-0.2289043
\end{pmatrix}
$$

Tìm véc tơ riêng trái từ hệ \((A^T-\lambda I)w=0\):

$$
w_2\approx \begin{pmatrix}
0.6534213\\
-0.7215563\\
0.2289043
\end{pmatrix}
$$

Mẫu số

$$
w_2^Tv_2=-1.0000000.
$$

Ma trận sau khi xuống thang:

$$
A_2=\begin{pmatrix}
0.4223 & 0.2652 & -0.3693\\
0.2652 & 0.1666 & -0.2320\\
-0.3693 & -0.2320 & 0.3230
\end{pmatrix}
$$

Tiếp tục áp dụng phương pháp lũy thừa cho \(A_2\), ta thu được giá trị riêng trội tiếp theo

$$
\lambda_3\approx 0.9118774.
$$

Véc tơ riêng tương ứng trong bài toán sau xuống thang là

$$
\tilde v_3\approx \begin{pmatrix}
0.6805031\\
0.4274390\\
-0.5951566
\end{pmatrix}.
$$

## Kết luận

Giá trị riêng trội tìm được là

$$
\lambda_1\approx 550.3986911.
$$

Véc tơ riêng tương ứng là

$$
v_1\approx \begin{pmatrix}
0.3315964\\
0.5446581\\
0.7703190
\end{pmatrix}.
$$

Các giá trị riêng tiếp theo tìm được bằng xuống thang:

- \(\lambda_2\approx 58.6894315\)
- \(\lambda_3\approx 0.9118774\)


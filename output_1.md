# Lời giải khai triển kỳ dị SVD

> Các phép tính bên trong dùng giá trị gốc. Ma trận trình bày 4 chữ số thập phân; các số liệu khác lấy 7 chữ số thập phân.

## Bài toán

Cho ma trận

$$
A=\begin{pmatrix}
3.0000 & 3.0000 & 7.0000\\
3.0000 & 9.0000 & 10.0000\\
8.0000 & 2.0000 & 10.0000\\
2.0000 & 10.0000 & 9.0000
\end{pmatrix}
$$

Ở đây \(A\in\mathbb{R}^{4\times 3}\).

## Bước 1. Quy SVD về bài toán trị riêng

Vì \(m\ge n\), ta lập

$$
B=A^TA.
$$

Khi đó các véc tơ kỳ dị phải \(v_i\) là véc tơ riêng của \(B\):

$$
Bv_i=\lambda_i v_i,\qquad \sigma_i=\sqrt{\lambda_i},\qquad u_i=\frac{Av_i}{\sigma_i}.
$$

Ta có

$$
B=\begin{pmatrix}
86.0000 & 72.0000 & 149.0000\\
72.0000 & 194.0000 & 221.0000\\
149.0000 & 221.0000 & 330.0000
\end{pmatrix}
$$

## Bước 2. Tìm trị riêng bằng phương pháp lũy thừa và xuống thang

Áp dụng phương pháp lũy thừa cho ma trận đối xứng ở trên, sau mỗi trị riêng thì xuống thang

Trước khi xuống thang, véc tơ riêng được chuẩn hóa theo chuẩn 2:

$$
q_k\leftarrow \frac{q_k}{\|q_k\|_2},\qquad \|q_k\|_2=1.
$$

Do đó công thức xuống thang dùng trong bài là

$$
M_{k+1}=M_k-\lambda_k q_kq_k^T
$$

để tìm trị riêng tiếp theo.

Sai số dư trong bảng được tính bởi

$$
r_i=\|Mq_i-\lambda_iq_i\|_2.
$$

Nếu \(r_i\) càng nhỏ thì cặp \((\lambda_i,q_i)\) càng gần với trị riêng và véc tơ riêng đúng.

| \(i\) | \(\lambda_i\) | \(\sigma_i=\sqrt{\lambda_i}\) | Số lặp | Sai số dư |
|---:|---:|---:|---:|---:|
| 1 | 550.3986911 | 23.4605774 | 6 | 0.0000199 |
| 2 | 58.6894315 | 7.6609028 | 5 | 0.0000002 |
| 3 | 0.9118774 | 0.9549227 | 2 | 0.0000000 |

## Bước 3. Giá trị kỳ dị và hạng

Với mỗi trị riêng của ma trận đối xứng vừa xét:

$$
\lambda_i>0\Rightarrow \sigma_i=\sqrt{\lambda_i}>0,\qquad \lambda_i=0\Rightarrow \sigma_i=0.
$$

Số giá trị kỳ dị khác 0 chính là hạng của ma trận.

Các giá trị kỳ dị khác 0 là

$$
\sigma_1=23.4605774,\quad \sigma_2=7.6609028,\quad \sigma_3=0.9549227.
$$

Suy ra \(\operatorname{rank}(A)=3\).

## Bước 4. Vector kỳ dị

Với \(\sigma_1=23.4605774\), ta có

$$
v_1\approx \begin{pmatrix}
0.3315963\\
0.5446581\\
0.7703190
\end{pmatrix},\qquad u_1\approx \begin{pmatrix}
0.3418925\\
0.5796917\\
0.4878515\\
0.5559388
\end{pmatrix}.
$$

Các véc tơ này đã được chuẩn hóa:

$$
\|v_1\|_2=1.0000000,\qquad \|u_1\|_2=1.0000000.
$$

Với \(\sigma_2=7.6609028\), ta có

$$
v_2\approx \begin{pmatrix}
0.6534214\\
-0.7215561\\
0.2289046
\end{pmatrix},\qquad u_2\approx \begin{pmatrix}
0.1824757\\
-0.2930065\\
0.7927663\\
-0.5023659
\end{pmatrix}.
$$

Các véc tơ này đã được chuẩn hóa:

$$
\|v_2\|_2=1.0000000,\qquad \|u_2\|_2=1.0000000.
$$

Với \(\sigma_3=0.9549227\), ta có

$$
v_3\approx \begin{pmatrix}
0.6805032\\
0.4274389\\
-0.5951566
\end{pmatrix},\qquad u_3\approx \begin{pmatrix}
0.8820294\\
0.0660855\\
-0.3637337\\
-0.2921553
\end{pmatrix}.
$$

Các véc tơ này đã được chuẩn hóa:

$$
\|v_3\|_2=1.0000000,\qquad \|u_3\|_2=1.0000000.
$$

## Bước 5. Khai triển SVD rút gọn

Đặt

$$
U_r=\begin{pmatrix}
0.3419 & 0.1825 & 0.8820\\
0.5797 & -0.2930 & 0.0661\\
0.4879 & 0.7928 & -0.3637\\
0.5559 & -0.5024 & -0.2922
\end{pmatrix}
$$

$$
\Sigma_r=\begin{pmatrix}
23.4606 & 0.0000 & 0.0000\\
0.0000 & 7.6609 & 0.0000\\
0.0000 & 0.0000 & 0.9549
\end{pmatrix}
$$

$$
V_r=\begin{pmatrix}
0.3316 & 0.6534 & 0.6805\\
0.5447 & -0.7216 & 0.4274\\
0.7703 & 0.2289 & -0.5952
\end{pmatrix}
$$

Khi đó

$$
A=U_r\Sigma_rV_r^T.
$$

Tương đương

$$
A=\sigma_1u_1v_1^T+\sigma_2u_2v_2^T+\sigma_3u_3v_3^T.
$$

## Kiểm tra tái tạo

Từ các thành phần SVD đã tìm được, ta tính lại

$$
\widehat A=U_r\Sigma_rV_r^T=\begin{pmatrix}
4.1463 & 3.7200 & 5.9974\\
3.0859 & 9.0539 & 9.9249\\
7.5273 & 1.7031 & 10.4134\\
1.6203 & 9.7615 & 9.3321
\end{pmatrix}
$$

Sai số tái tạo theo chuẩn Frobenius là

$$
E=\|A-\widehat A\|_F=1.9098455.
$$

Vì đã lấy đủ số thành phần cần thiết, về lý thuyết ta có khai triển đúng của ma trận.

## Kiểm tra trực chuẩn

Ta kiểm tra các cột của \(U_r\) và \(V_r\):

$$
U_r^TU_r=\begin{pmatrix}
1.0000 & 0.0000 & 0.0000\\
0.0000 & 1.0000 & 0.0000\\
0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

$$
V_r^TV_r=\begin{pmatrix}
1.0000 & 0.0000 & 0.0000\\
0.0000 & 1.0000 & 0.0000\\
0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

## Nghịch đảo suy rộng

Theo công thức Moore-Penrose,

$$
A^\dagger=V_r\Sigma_r^{-1}U_r^T.
$$

Suy ra

$$
A^\dagger=\begin{pmatrix}
0.6490 & 0.0303 & -0.1847 & -0.2432\\
0.3856 & 0.0706 & -0.2262 & -0.0706\\
-0.5330 & -0.0309 & 0.2664 & 0.1853
\end{pmatrix}
$$

## Số điều kiện

Nếu các giá trị kỳ dị khác 0 phủ đủ số chiều cần xét thì

$$
\operatorname{cond}(A)=\frac{\sigma_{\max}}{\sigma_{\min}}.
$$

Nếu tồn tại \(\sigma_i=0\), ma trận bị suy biến theo nghĩa SVD và số điều kiện bằng \(+\infty\).

Trong bài này ta có

$$
\operatorname{cond}(A)=\frac{\sigma_{\max}}{\sigma_{\min}}=24.5680370.
$$

## Kết luận

Khai triển kỳ dị rút gọn của ma trận là

$$
A=U_r\Sigma_rV_r^T
$$

với các ma trận \(U_r,\Sigma_r,V_r\) đã tính ở trên; sai số tái tạo là \(E=1.9098455\).


# Lời giải khai triển kỳ dị SVD

> Các phép tính bên trong dùng giá trị gốc. Ma trận trình bày 4 chữ số thập phân; các số liệu khác lấy 7 chữ số thập phân.

## Bài toán

Cho ma trận

$$
A=\begin{pmatrix}
1.0000 & 2.0000\\
2.0000 & 4.0000\\
3.0000 & 6.0000
\end{pmatrix}
$$

Ở đây \(A\in\mathbb{R}^{3\times 2}\).

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
14.0000 & 28.0000\\
28.0000 & 56.0000
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
| 1 | 70.0000000 | 8.3666003 | 2 | 0.0000000 |
| 2 | 0.0000000 | 0.0000000 | 1 | 0.0000000 |

## Bước 3. Giá trị kỳ dị và hạng

Với mỗi trị riêng của ma trận đối xứng vừa xét:

$$
\lambda_i>0\Rightarrow \sigma_i=\sqrt{\lambda_i}>0,\qquad \lambda_i=0\Rightarrow \sigma_i=0.
$$

Số giá trị kỳ dị khác 0 chính là hạng của ma trận.

Các giá trị kỳ dị khác 0 là

$$
\sigma_1=8.3666003.
$$

Suy ra \(\operatorname{rank}(A)=1\).

## Bước 4. Vector kỳ dị

Với \(\sigma_1=8.3666003\), ta có

$$
v_1\approx \begin{pmatrix}
0.4472136\\
0.8944272
\end{pmatrix},\qquad u_1\approx \begin{pmatrix}
0.2672612\\
0.5345225\\
0.8017837
\end{pmatrix}.
$$

Các véc tơ này đã được chuẩn hóa:

$$
\|v_1\|_2=1.0000000,\qquad \|u_1\|_2=1.0000000.
$$

## Bước 5. Khai triển SVD rút gọn

Đặt

$$
U_r=\begin{pmatrix}
0.2673\\
0.5345\\
0.8018
\end{pmatrix}
$$

$$
\Sigma_r=\begin{pmatrix}
8.3666
\end{pmatrix}
$$

$$
V_r=\begin{pmatrix}
0.4472\\
0.8944
\end{pmatrix}
$$

Khi đó

$$
A=U_r\Sigma_rV_r^T.
$$

Tương đương

$$
A=\sigma_1u_1v_1^T.
$$

## Kiểm tra tái tạo

Từ các thành phần SVD đã tìm được, ta tính lại

$$
\widehat A=U_r\Sigma_rV_r^T=\begin{pmatrix}
1.0000 & 2.0000\\
2.0000 & 4.0000\\
3.0000 & 6.0000
\end{pmatrix}
$$

Sai số tái tạo theo chuẩn Frobenius là

$$
E=\|A-\widehat A\|_F=0.0000000.
$$

Vì đã lấy đủ số thành phần cần thiết, về lý thuyết ta có khai triển đúng của ma trận.

## Nghịch đảo suy rộng

Theo công thức Moore-Penrose,

$$
A^\dagger=V_r\Sigma_r^{-1}U_r^T.
$$

Suy ra

$$
A^\dagger=\begin{pmatrix}
0.0143 & 0.0286 & 0.0429\\
0.0286 & 0.0571 & 0.0857
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
\operatorname{cond}(A)=\frac{\sigma_{\max}}{\sigma_{\min}}=+\infty.
$$

## Kết luận

Khai triển kỳ dị rút gọn của ma trận là

$$
A=U_r\Sigma_rV_r^T
$$

với các ma trận \(U_r,\Sigma_r,V_r\) đã tính ở trên; sai số tái tạo là \(E=0.0000000\).

# Lời giải khai triển kỳ dị SVD

> Các phép tính bên trong dùng giá trị gốc. Ma trận trình bày 4 chữ số thập phân; các số liệu khác lấy 7 chữ số thập phân.

## Bài toán

Cho ma trận

$$
A=\begin{pmatrix}
3.0000 & 4.0000 & 0.0000\\
0.0000 & 0.0000 & 2.0000
\end{pmatrix}
$$

Ở đây \(A\in\mathbb{R}^{2\times 3}\).

## Bước 1. Quy SVD về bài toán trị riêng

Vì \(m<n\), ta lập

$$
C=AA^T.
$$

Khi đó các véc tơ kỳ dị trái \(u_i\) là véc tơ riêng của \(C\):

$$
Cu_i=\lambda_i u_i,\qquad \sigma_i=\sqrt{\lambda_i},\qquad v_i=\frac{A^Tu_i}{\sigma_i}.
$$

Ta có

$$
C=\begin{pmatrix}
25.0000 & 0.0000\\
0.0000 & 4.0000
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
| 1 | 25.0000000 | 5.0000000 | 10 | 0.0000005 |
| 2 | 4.0000000 | 2.0000000 | 2 | 0.0000000 |

## Bước 3. Giá trị kỳ dị và hạng

Với mỗi trị riêng của ma trận đối xứng vừa xét:

$$
\lambda_i>0\Rightarrow \sigma_i=\sqrt{\lambda_i}>0,\qquad \lambda_i=0\Rightarrow \sigma_i=0.
$$

Số giá trị kỳ dị khác 0 chính là hạng của ma trận.

Các giá trị kỳ dị khác 0 là

$$
\sigma_1=5.0000000,\quad \sigma_2=2.0000000.
$$

Suy ra \(\operatorname{rank}(A)=2\).

## Bước 4. Vector kỳ dị

Với \(\sigma_1=5.0000000\), ta có

$$
v_1\approx \begin{pmatrix}
0.6000000\\
0.8000000\\
0.0000000
\end{pmatrix},\qquad u_1\approx \begin{pmatrix}
1.0000000\\
0.0000000
\end{pmatrix}.
$$

Các véc tơ này đã được chuẩn hóa:

$$
\|v_1\|_2=1.0000000,\qquad \|u_1\|_2=1.0000000.
$$

Với \(\sigma_2=2.0000000\), ta có

$$
v_2\approx \begin{pmatrix}
0.0000002\\
0.0000003\\
-1.0000000
\end{pmatrix},\qquad u_2\approx \begin{pmatrix}
0.0000001\\
-1.0000000
\end{pmatrix}.
$$

Các véc tơ này đã được chuẩn hóa:

$$
\|v_2\|_2=1.0000000,\qquad \|u_2\|_2=1.0000000.
$$

## Bước 5. Khai triển SVD rút gọn

Đặt

$$
U_r=\begin{pmatrix}
1.0000 & 0.0000\\
0.0000 & -1.0000
\end{pmatrix}
$$

$$
\Sigma_r=\begin{pmatrix}
5.0000 & 0.0000\\
0.0000 & 2.0000
\end{pmatrix}
$$

$$
V_r=\begin{pmatrix}
0.6000 & 0.0000\\
0.8000 & 0.0000\\
0.0000 & -1.0000
\end{pmatrix}
$$

Khi đó

$$
A=U_r\Sigma_rV_r^T.
$$

Tương đương

$$
A=\sigma_1u_1v_1^T+\sigma_2u_2v_2^T.
$$

## Kiểm tra tái tạo

Từ các thành phần SVD đã tìm được, ta tính lại

$$
\widehat A=U_r\Sigma_rV_r^T=\begin{pmatrix}
3.0000 & 4.0000 & 0.0000\\
0.0000 & 0.0000 & 2.0000
\end{pmatrix}
$$

Sai số tái tạo theo chuẩn Frobenius là

$$
E=\|A-\widehat A\|_F=0.0000006.
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
0.1200 & 0.0000\\
0.1600 & 0.0000\\
0.0000 & 0.5000
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
\operatorname{cond}(A)=\frac{\sigma_{\max}}{\sigma_{\min}}=2.5000000.
$$

## Kết luận

Khai triển kỳ dị rút gọn của ma trận là

$$
A=U_r\Sigma_rV_r^T
$$

với các ma trận \(U_r,\Sigma_r,V_r\) đã tính ở trên; sai số tái tạo là \(E=0.0000006\).

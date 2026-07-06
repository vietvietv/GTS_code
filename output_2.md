# Lời giải khai triển kỳ dị SVD

> Các phép tính bên trong dùng giá trị gốc. Ma trận trình bày 4 chữ số thập phân; các số liệu khác lấy 7 chữ số thập phân.

## Bài toán

Cho ma trận

$$
A=\begin{pmatrix}
3.0000 & 1.0000 & 0.0000 & 2.0000 & 1.0000 & 0.0000\\
0.0000 & 2.0000 & 4.0000 & 1.0000 & 0.0000 & 3.0000\\
1.0000 & 0.0000 & 2.0000 & 3.0000 & 2.0000 & 1.0000\\
2.0000 & 3.0000 & 1.0000 & 0.0000 & 4.0000 & 2.0000
\end{pmatrix}
$$

Ở đây \(A\in\mathbb{R}^{4\times 6}\).

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
15.0000 & 4.0000 & 11.0000 & 13.0000\\
4.0000 & 30.0000 & 14.0000 & 16.0000\\
11.0000 & 14.0000 & 19.0000 & 14.0000\\
13.0000 & 16.0000 & 14.0000 & 34.0000
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
| 1 | 63.4756862 | 7.9671630 | 12 | 0.0000026 |
| 2 | 19.1887544 | 4.3804971 | 34 | 0.0000017 |
| 3 | 11.8534265 | 3.4428806 | 15 | 0.0000004 |
| 4 | 3.4821329 | 1.8660474 | 2 | 0.0000000 |

## Bước 3. Giá trị kỳ dị và hạng

Với mỗi trị riêng của ma trận đối xứng vừa xét:

$$
\lambda_i>0\Rightarrow \sigma_i=\sqrt{\lambda_i}>0,\qquad \lambda_i=0\Rightarrow \sigma_i=0.
$$

Số giá trị kỳ dị khác 0 chính là hạng của ma trận.

Các giá trị kỳ dị khác 0 là

$$
\sigma_1=7.9671630,\quad \sigma_2=4.3804971,\quad \sigma_3=3.4428806,\quad \sigma_4=1.8660474.
$$

Suy ra \(\operatorname{rank}(A)=4\).

## Bước 4. Vector kỳ dị

Với \(\sigma_1=7.9671630\), ta có

$$
v_1\approx \begin{pmatrix}
0.3375277\\
0.4157898\\
0.4609393\\
0.3157807\\
0.4753503\\
0.4183940
\end{pmatrix},\qquad u_1\approx \begin{pmatrix}
0.3182168\\
0.5329751\\
0.4488226\\
0.6428328
\end{pmatrix}.
$$

Các véc tơ này đã được chuẩn hóa:

$$
\|v_1\|_2=1.0000000,\qquad \|u_1\|_2=1.0000000.
$$

Với \(\sigma_2=4.3804971\), ta có

$$
v_2\approx \begin{pmatrix}
0.5211247\\
0.0229691\\
-0.6133095\\
0.0621068\\
0.4732233\\
-0.3520729
\end{pmatrix},\qquad u_2\approx \begin{pmatrix}
0.4985235\\
-0.7764897\\
0.0171671\\
0.3850239
\end{pmatrix}.
$$

Các véc tơ này đã được chuẩn hóa:

$$
\|v_2\|_2=1.0000000,\qquad \|u_2\|_2=1.0000000.
$$

Với \(\sigma_3=3.4428806\), ta có

$$
v_3\approx \begin{pmatrix}
0.1894097\\
-0.4459856\\
0.1559358\\
0.7990175\\
-0.2439316\\
-0.2073000
\end{pmatrix},\qquad u_3\approx \begin{pmatrix}
0.4288115\\
-0.0264633\\
0.6399207\\
-0.6371200
\end{pmatrix}.
$$

Các véc tơ này đã được chuẩn hóa:

$$
\|v_3\|_2=1.0000000,\qquad \|u_3\|_2=1.0000000.
$$

Với \(\sigma_4=1.8660474\), ta có

$$
v_4\approx \begin{pmatrix}
0.5702079\\
0.4348468\\
-0.0466880\\
-0.0909117\\
-0.6893499\\
0.0111028
\end{pmatrix},\qquad u_4\approx \begin{pmatrix}
0.6828859\\
0.3351137\\
-0.6235103\\
-0.1805562
\end{pmatrix}.
$$

Các véc tơ này đã được chuẩn hóa:

$$
\|v_4\|_2=1.0000000,\qquad \|u_4\|_2=1.0000000.
$$

## Bước 5. Khai triển SVD rút gọn

Đặt

$$
U_r=\begin{pmatrix}
0.3182 & 0.4985 & 0.4288 & 0.6829\\
0.5330 & -0.7765 & -0.0265 & 0.3351\\
0.4488 & 0.0172 & 0.6399 & -0.6235\\
0.6428 & 0.3850 & -0.6371 & -0.1806
\end{pmatrix}
$$

$$
\Sigma_r=\begin{pmatrix}
7.9672 & 0.0000 & 0.0000 & 0.0000\\
0.0000 & 4.3805 & 0.0000 & 0.0000\\
0.0000 & 0.0000 & 3.4429 & 0.0000\\
0.0000 & 0.0000 & 0.0000 & 1.8660
\end{pmatrix}
$$

$$
V_r=\begin{pmatrix}
0.3375 & 0.5211 & 0.1894 & 0.5702\\
0.4158 & 0.0230 & -0.4460 & 0.4348\\
0.4609 & -0.6133 & 0.1559 & -0.0467\\
0.3158 & 0.0621 & 0.7990 & -0.0909\\
0.4754 & 0.4732 & -0.2439 & -0.6893\\
0.4184 & -0.3521 & -0.2073 & 0.0111
\end{pmatrix}
$$

Khi đó

$$
A=U_r\Sigma_rV_r^T.
$$

Tương đương

$$
A=\sigma_1u_1v_1^T+\sigma_2u_2v_2^T+\sigma_3u_3v_3^T+\sigma_4u_4v_4^T.
$$

## Kiểm tra tái tạo

Từ các thành phần SVD đã tìm được, ta tính lại

$$
\widehat A=U_r\Sigma_rV_r^T=\begin{pmatrix}
3.0000 & 1.0000 & 0.0000 & 2.0000 & 1.0000 & 0.0000\\
0.0000 & 2.0000 & 4.0000 & 1.0000 & 0.0000 & 3.0000\\
1.0000 & 0.0000 & 2.0000 & 3.0000 & 2.0000 & 1.0000\\
2.0000 & 3.0000 & 1.0000 & 0.0000 & 4.0000 & 2.0000
\end{pmatrix}
$$

Sai số tái tạo theo chuẩn Frobenius là

$$
E=\|A-\widehat A\|_F=0.0000015.
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
0.3050 & 0.0311 & -0.1343 & -0.0172\\
0.1228 & 0.1053 & -0.2047 & 0.0760\\
-0.0491 & 0.1300 & 0.0681 & -0.0411\\
0.0859 & -0.0124 & 0.1969 & -0.1081\\
-0.2098 & -0.1740 & 0.2136 & 0.1918\\
-0.0451 & 0.0940 & -0.0201 & 0.0401
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
\operatorname{cond}(A)=\frac{\sigma_{\max}}{\sigma_{\min}}=4.2695395.
$$

## Kết luận

Khai triển kỳ dị rút gọn của ma trận là

$$
A=U_r\Sigma_rV_r^T
$$

với các ma trận \(U_r,\Sigma_r,V_r\) đã tính ở trên; sai số tái tạo là \(E=0.0000015\).


# Lời giải bằng phương pháp Danielevsky

> Các phép tính bên trong sử dụng giá trị gốc, không dùng giá trị đã làm tròn để tính tiếp.
> Khi trình bày, ma trận được làm tròn 4 chữ số thập phân; các số liệu khác làm tròn 7 chữ số thập phân.

## Bài toán

Cho ma trận

$$
A=\begin{pmatrix}
4.0000 & -1.0000 & -1.0000\\
0.0000 & 1.0000 & 0.0000\\
2.0000 & 2.0000 & 1.0000
\end{pmatrix}
$$

Dùng phương pháp Danielevsky tìm đa thức đặc trưng, giá trị riêng và véc tơ riêng của ma trận \(A\).

## Quá trình biến đổi

Đặt

$$
A^{(1)}=A.
$$

### 1. Bước đưa hàng 3 về dạng Frobenius

Xét hàng 3, ta có \(a_{3,2}=2.0000000\ne0\).

Do đó thuộc trường hợp 1 của phương pháp Danielevsky.

Chọn \(M\) sao cho hàng 2 của \(M\) bằng phần tương ứng của hàng 3.

Khi đó \(A_{\mathrm{mới}}=MA_{\mathrm{cũ}}M^{-1}\).

$$
A_{\mathrm{cũ}}=\begin{pmatrix}
4.0000 & -1.0000 & -1.0000\\
0.0000 & 1.0000 & 0.0000\\
2.0000 & 2.0000 & 1.0000
\end{pmatrix}
$$

$$
M=\begin{pmatrix}
1.0000 & 0.0000 & 0.0000\\
2.0000 & 2.0000 & 1.0000\\
0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

$$
M^{-1}=\begin{pmatrix}
1.0000 & 0.0000 & 0.0000\\
-1.0000 & 0.5000 & -0.5000\\
0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

$$
A_{\mathrm{mới}}=\begin{pmatrix}
5.0000 & -0.5000 & -0.5000\\
8.0000 & 1.0000 & -2.0000\\
0.0000 & 1.0000 & 0.0000
\end{pmatrix}
$$

### 2. Bước đưa hàng 2 về dạng Frobenius

Xét hàng 2, ta có \(a_{2,1}=8.0000000\ne0\).

Do đó thuộc trường hợp 1 của phương pháp Danielevsky.

Chọn \(M\) sao cho hàng 1 của \(M\) bằng phần tương ứng của hàng 2.

Khi đó \(A_{\mathrm{mới}}=MA_{\mathrm{cũ}}M^{-1}\).

$$
A_{\mathrm{cũ}}=\begin{pmatrix}
5.0000 & -0.5000 & -0.5000\\
8.0000 & 1.0000 & -2.0000\\
0.0000 & 1.0000 & 0.0000
\end{pmatrix}
$$

$$
M=\begin{pmatrix}
8.0000 & 1.0000 & -2.0000\\
0.0000 & 1.0000 & 0.0000\\
0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

$$
M^{-1}=\begin{pmatrix}
0.1250 & -0.1250 & 0.2500\\
0.0000 & 1.0000 & 0.0000\\
0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

$$
A_{\mathrm{mới}}=\begin{pmatrix}
6.0000 & -11.0000 & 6.0000\\
1.0000 & 0.0000 & 0.0000\\
0.0000 & 1.0000 & 0.0000
\end{pmatrix}
$$

Sau các phép biến đổi đồng dạng, thu được ma trận Frobenius

$$
F=\begin{pmatrix}
6.0000 & -11.0000 & 6.0000\\
1.0000 & 0.0000 & 0.0000\\
0.0000 & 1.0000 & 0.0000
\end{pmatrix}
$$

Ma trận biến đổi tổng hợp là

$$
P=\begin{pmatrix}
10.0000 & 2.0000 & -1.0000\\
2.0000 & 2.0000 & 1.0000\\
0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

và

$$
P^{-1}=\begin{pmatrix}
0.1250 & -0.1250 & 0.2500\\
-0.1250 & 0.6250 & -0.7500\\
0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

Khi đó

$$
F=PAP^{-1}.
$$

## Đa thức đặc trưng

Với ma trận Frobenius dạng

$$
F=\begin{pmatrix}-p_1&-p_2&\cdots&-p_n\\1&0&\cdots&0\\0&1&\cdots&0\\\vdots&\vdots&\ddots&\vdots\end{pmatrix},
$$

ta có

$$
p(\lambda)=\lambda^n+p_1\lambda^{n-1}+\cdots+p_n.
$$

\(p_1=-6.0000000\).
\(p_2=11.0000000\).
\(p_3=-6.0000000\).

Do đó

$$
p(\lambda)=\lambda^{3} - 6.0000000\lambda^{2} + 11.0000000\lambda - 6.0000000.
$$

## Giá trị riêng và véc tơ riêng

Giải phương trình

$$
\lambda^{3} - 6.0000000\lambda^{2} + 11.0000000\lambda - 6.0000000=0
$$

thu được các giá trị riêng:

- \(\lambda_1=1.0000000\)
- \(\lambda_2=2.0000000\)
- \(\lambda_3=3.0000000\)

Với mỗi giá trị riêng \(\lambda\), véc tơ riêng của khối Frobenius có dạng

$$
u_\lambda=\begin{pmatrix}\lambda^{n-1}\\\lambda^{n-2}\\\vdots\\1\end{pmatrix}.
$$

Vì \(F=PAP^{-1}\), véc tơ riêng của ma trận ban đầu là

$$
x_\lambda=P^{-1}u_\lambda.
$$

Với \(\lambda_1=1.0000000\):

$$
u_1=\begin{pmatrix}
1.0000000\\
1.0000000\\
1.0000000
\end{pmatrix}
$$

$$
x_1=P^{-1}u_1=\begin{pmatrix}
0.2500000\\
-0.2500000\\
1.0000000
\end{pmatrix}
$$

Với \(\lambda_2=2.0000000\):

$$
u_2=\begin{pmatrix}
4.0000000\\
2.0000000\\
1.0000000
\end{pmatrix}
$$

$$
x_2=P^{-1}u_2=\begin{pmatrix}
0.5000000\\
0.0000000\\
1.0000000
\end{pmatrix}
$$

Với \(\lambda_3=3.0000000\):

$$
u_3=\begin{pmatrix}
9.0000000\\
3.0000000\\
1.0000000
\end{pmatrix}
$$

$$
x_3=P^{-1}u_3=\begin{pmatrix}
1.0000000\\
0.0000000\\
1.0000000
\end{pmatrix}
$$

## Kết luận

Đa thức đặc trưng của ma trận là

$$
p(\lambda)=\lambda^{3} - 6.0000000\lambda^{2} + 11.0000000\lambda - 6.0000000.
$$

Các cặp giá trị riêng và véc tơ riêng tương ứng là:

$$
\lambda_1=1.0000000,\qquad x_1=\begin{pmatrix}
0.2500000\\
-0.2500000\\
1.0000000
\end{pmatrix}.
$$

$$
\lambda_2=2.0000000,\qquad x_2=\begin{pmatrix}
0.5000000\\
0.0000000\\
1.0000000
\end{pmatrix}.
$$

$$
\lambda_3=3.0000000,\qquad x_3=\begin{pmatrix}
1.0000000\\
0.0000000\\
1.0000000
\end{pmatrix}.
$$

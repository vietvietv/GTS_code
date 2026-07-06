# Lời giải bằng phương pháp Danielevsky

> Các phép tính bên trong sử dụng giá trị gốc, không dùng giá trị đã làm tròn để tính tiếp.
> Khi trình bày, ma trận được làm tròn 4 chữ số thập phân; các số liệu khác làm tròn 7 chữ số thập phân.

## Bài toán

Cho ma trận

$$
A=\begin{pmatrix}
1.0000 & 2.0000 & 7.0000 & 8.0000\\
3.0000 & 4.0000 & 9.0000 & 10.0000\\
0.0000 & 0.0000 & 5.0000 & -6.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000
\end{pmatrix}
$$

Dùng phương pháp Danielevsky tìm đa thức đặc trưng, giá trị riêng và véc tơ riêng của ma trận \(A\).

## Quá trình biến đổi

Đặt

$$
A^{(1)}=A.
$$

### 1. Bước đưa hàng 4 về dạng Frobenius

Xét hàng 4, ta có \(a_{4,3}=1.0000000\ne0\).

Do đó thuộc trường hợp 1 của phương pháp Danielevsky.

Chọn \(M\) sao cho hàng 3 của \(M\) bằng phần tương ứng của hàng 4.

Khi đó \(A_{\mathrm{mới}}=MA_{\mathrm{cũ}}M^{-1}\).

$$
A_{\mathrm{cũ}}=\begin{pmatrix}
1.0000 & 2.0000 & 7.0000 & 8.0000\\
3.0000 & 4.0000 & 9.0000 & 10.0000\\
0.0000 & 0.0000 & 5.0000 & -6.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000
\end{pmatrix}
$$

$$
M=\begin{pmatrix}
1.0000 & 0.0000 & 0.0000 & 0.0000\\
0.0000 & 1.0000 & 0.0000 & 0.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000\\
0.0000 & 0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

$$
M^{-1}=\begin{pmatrix}
1.0000 & 0.0000 & 0.0000 & 0.0000\\
0.0000 & 1.0000 & 0.0000 & 0.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000\\
0.0000 & 0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

$$
A_{\mathrm{mới}}=\begin{pmatrix}
1.0000 & 2.0000 & 7.0000 & 8.0000\\
3.0000 & 4.0000 & 9.0000 & 10.0000\\
0.0000 & 0.0000 & 5.0000 & -6.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000
\end{pmatrix}
$$

### 2. Bước tách khối tại hàng 3

Xét hàng 3, ta có \(a_{3,1}=\cdots=a_{3,2}=0\).

Do đó thuộc trường hợp 3 của phương pháp Danielevsky.

Ma trận khi đó có dạng khối tam giác trên \(\begin{pmatrix}A_1&B\\0&F_1\end{pmatrix}\).

Suy ra đa thức đặc trưng bằng tích đa thức đặc trưng của các khối đường chéo.

Ta tách khối Frobenius cấp 2 ở góc dưới bên phải và tiếp tục xử lý ma trận con cấp 2 ở phía trên.

$$
A_{\mathrm{hiện\ tại}}=\begin{pmatrix}
1.0000 & 2.0000 & 7.0000 & 8.0000\\
3.0000 & 4.0000 & 9.0000 & 10.0000\\
0.0000 & 0.0000 & 5.0000 & -6.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000
\end{pmatrix}
$$

$$
F_1=\begin{pmatrix}
5.0000 & -6.0000\\
1.0000 & 0.0000
\end{pmatrix}
$$

### 3. Bước đưa hàng 2 về dạng Frobenius

Xét hàng 2, ta có \(a_{2,1}=3.0000000\ne0\).

Do đó thuộc trường hợp 1 của phương pháp Danielevsky.

Chọn \(M\) sao cho hàng 1 của \(M\) bằng phần tương ứng của hàng 2.

Khi đó \(A_{\mathrm{mới}}=MA_{\mathrm{cũ}}M^{-1}\).

$$
A_{\mathrm{cũ}}=\begin{pmatrix}
1.0000 & 2.0000 & 7.0000 & 8.0000\\
3.0000 & 4.0000 & 9.0000 & 10.0000\\
0.0000 & 0.0000 & 5.0000 & -6.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000
\end{pmatrix}
$$

$$
M=\begin{pmatrix}
3.0000 & 4.0000 & 0.0000 & 0.0000\\
0.0000 & 1.0000 & 0.0000 & 0.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000\\
0.0000 & 0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

$$
M^{-1}=\begin{pmatrix}
0.3333 & -1.3333 & 0.0000 & 0.0000\\
0.0000 & 1.0000 & 0.0000 & 0.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000\\
0.0000 & 0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

$$
A_{\mathrm{mới}}=\begin{pmatrix}
5.0000 & 2.0000 & 57.0000 & 64.0000\\
1.0000 & 0.0000 & 9.0000 & 10.0000\\
0.0000 & 0.0000 & 5.0000 & -6.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000
\end{pmatrix}
$$

Sau các phép biến đổi đồng dạng, thu được ma trận Frobenius dạng khối

$$
F=\begin{pmatrix}
5.0000 & 2.0000 & 57.0000 & 64.0000\\
1.0000 & 0.0000 & 9.0000 & 10.0000\\
0.0000 & 0.0000 & 5.0000 & -6.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000
\end{pmatrix}
$$

Ma trận biến đổi tổng hợp là

$$
P=\begin{pmatrix}
3.0000 & 4.0000 & 0.0000 & 0.0000\\
0.0000 & 1.0000 & 0.0000 & 0.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000\\
0.0000 & 0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

và

$$
P^{-1}=\begin{pmatrix}
0.3333 & -1.3333 & 0.0000 & 0.0000\\
0.0000 & 1.0000 & 0.0000 & 0.0000\\
0.0000 & 0.0000 & 1.0000 & 0.0000\\
0.0000 & 0.0000 & 0.0000 & 1.0000
\end{pmatrix}
$$

Khi đó

$$
F=PAP^{-1}.
$$

## Đa thức đặc trưng

Vì ma trận cuối có dạng khối tam giác trên, nên đa thức đặc trưng bằng tích đa thức đặc trưng của các khối đường chéo.

Khối \(F_1\) cấp 2:

$$
F_1=\begin{pmatrix}
5.0000 & 2.0000\\
1.0000 & 0.0000
\end{pmatrix}
$$

$$
p_1(\lambda)=\lambda^{2} - 5.0000000\lambda - 2.0000000.
$$

Khối \(F_2\) cấp 2:

$$
F_2=\begin{pmatrix}
5.0000 & -6.0000\\
1.0000 & 0.0000
\end{pmatrix}
$$

$$
p_2(\lambda)=\lambda^{2} - 5.0000000\lambda + 6.0000000.
$$

Do đó

$$
p(\lambda)=\lambda^{4} - 10.0000000\lambda^{3} + 29.0000000\lambda^{2} - 20.0000000\lambda - 12.0000000.
$$

## Giá trị riêng và véc tơ riêng

Giải phương trình

$$
\lambda^{4} - 10.0000000\lambda^{3} + 29.0000000\lambda^{2} - 20.0000000\lambda - 12.0000000=0
$$

thu được các giá trị riêng:

- \(\lambda_1=-0.3722813\)
- \(\lambda_2=2.0000000\)
- \(\lambda_3=3.0000000\)
- \(\lambda_4=5.3722813\)

Với mỗi giá trị riêng \(\lambda\), ta tìm véc tơ riêng \(u_\lambda\) của ma trận Frobenius khối bằng cách giải

$$
(F-\lambda I)u_\lambda=0.
$$

Vì \(F=PAP^{-1}\), véc tơ riêng của ma trận ban đầu là

$$
x_\lambda=P^{-1}u_\lambda.
$$

Với \(\lambda_1=-0.3722813\):

$$
u_1=\begin{pmatrix}
-0.3722813\\
1.0000000\\
0.0000000\\
0.0000000
\end{pmatrix}
$$

$$
x_1=P^{-1}u_1=\begin{pmatrix}
-1.4574271\\
1.0000000\\
0.0000000\\
0.0000000
\end{pmatrix}
$$

Với \(\lambda_2=2.0000000\):

$$
u_2=\begin{pmatrix}
-1.0000000\\
-0.2281553\\
0.0388350\\
0.0194175
\end{pmatrix}
$$

$$
x_2=P^{-1}u_2=\begin{pmatrix}
-0.0291262\\
-0.2281553\\
0.0388350\\
0.0194175
\end{pmatrix}
$$

Với \(\lambda_3=3.0000000\):

$$
u_3=\begin{pmatrix}
-1.0000000\\
-0.2066752\\
0.0308087\\
0.0102696
\end{pmatrix}
$$

$$
x_3=P^{-1}u_3=\begin{pmatrix}
-0.0577664\\
-0.2066752\\
0.0308087\\
0.0102696
\end{pmatrix}
$$

Với \(\lambda_4=5.3722813\):

$$
u_4=\begin{pmatrix}
1.0000000\\
0.1861407\\
0.0000000\\
0.0000000
\end{pmatrix}
$$

$$
x_4=P^{-1}u_4=\begin{pmatrix}
0.0851458\\
0.1861407\\
0.0000000\\
0.0000000
\end{pmatrix}
$$

## Kết luận

Đa thức đặc trưng của ma trận là

$$
p(\lambda)=\lambda^{4} - 10.0000000\lambda^{3} + 29.0000000\lambda^{2} - 20.0000000\lambda - 12.0000000.
$$

Các cặp giá trị riêng và véc tơ riêng tương ứng là:

$$
\lambda_1=-0.3722813,\qquad x_1=\begin{pmatrix}
-1.4574271\\
1.0000000\\
0.0000000\\
0.0000000
\end{pmatrix}.
$$

$$
\lambda_2=2.0000000,\qquad x_2=\begin{pmatrix}
-0.0291262\\
-0.2281553\\
0.0388350\\
0.0194175
\end{pmatrix}.
$$

$$
\lambda_3=3.0000000,\qquad x_3=\begin{pmatrix}
-0.0577664\\
-0.2066752\\
0.0308087\\
0.0102696
\end{pmatrix}.
$$

$$
\lambda_4=5.3722813,\qquad x_4=\begin{pmatrix}
0.0851458\\
0.1861407\\
0.0000000\\
0.0000000
\end{pmatrix}.
$$

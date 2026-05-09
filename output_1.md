=== CHƯƠNG TRÌNH GIẢI HỆ PHƯƠNG TRÌNH PHI TUYẾN (NEWTON CẢI TIẾN) ===


============================================================
### 1. Thông tin hệ phương trình và Thiết lập ban đầu
* **Hệ phương trình:** $F(X) = \left[\begin{matrix}x^{2} + 2 y^{2} - 8 & x y + 5 y^{2} - 4\end{matrix}\right]^T$
* **Ma trận Jacobi tổng quát:** $J(X) = \left[\begin{matrix}2 x & 4 y\\y & x + 10 y\end{matrix}\right]$
* **Giá trị khởi đầu:** $X_0 = \begin{pmatrix} 2.3000000 \\ -1.1000000 \end{pmatrix}$

### 2. Chi tiết bước lặp 1 (k = 0)
* **Giá trị hàm số tại $X_0$:**
  $$F(X_0) = \begin{pmatrix} -0.2900000 \\ -0.4800000 \end{pmatrix}$$
* **Ma trận Jacobi tại $X_0$:**
  $$J(X_0) = \begin{pmatrix} 4.6000000 & -4.4000000 \\ -1.1000000 & -8.7000000 \end{pmatrix}$$
* **Ma trận nghịch đảo Jacobi (cố định cho mọi bước):**
  $$J(X_0)^{-1} = \begin{pmatrix} 0.1939367 & -0.0980829 \\ -0.0245207 & -0.1025412 \end{pmatrix}$$
* **Tính toán nghiệm xấp xỉ $X_1$:**
  $$X_1 = X_0 - J(X_0)^{-1} \cdot F(X_0)$$
  $$X_1 = \begin{pmatrix} 2.3000000 \\ -1.1000000 \end{pmatrix} - \begin{pmatrix} 0.1939367 & -0.0980829 \\ -0.0245207 & -0.1025412 \end{pmatrix} \cdot \begin{pmatrix} -0.2900000 \\ -0.4800000 \end{pmatrix}$$
  $$X_1 = \begin{pmatrix} 2.3091618 \\ -1.1563308 \end{pmatrix}$$

### 3. Các bước lặp tiếp theo
Sử dụng công thức lặp Newton cải tiến: $X_{k+1} = X_k - J(X_0)^{-1} F(X_k)$

* **Tại bước $k=1$ (Tìm $X_2$):**
  - $F(X_1) = \begin{pmatrix} 0.0064303 \\ 0.0153497 \end{pmatrix}^T$
  - $X_2 = \begin{pmatrix} 2.3094203 \\ -1.1545992 \end{pmatrix}^T$

### 4. Bảng tổng hợp các giá trị xấp xỉ
| k | $x_1$ | $x_2$ |
| :--- | :--- | :--- |
| 0 | 2.3000000 | -1.1000000 |
| 1 | 2.3091618 | -1.1563308 |
| 2 | 2.3094203 | -1.1545992 |
| 3 | 2.3093999 | -1.1547068 |
| 4 | 2.3094012 | -1.1547002 |
| 5 | 2.3094011 | -1.1547006 |

=> **KẾT LUẬN:** Nghiệm xấp xỉ của hệ sau 5 bước lặp là:
  $$X \approx \begin{pmatrix} 2.3094011 \\ -1.1547006 \end{pmatrix}$$
# SVD Solver

Tool giai SVD theo dung huong bai giang: quy ve bai toan tri rieng cua `A^T A` hoac `A A^T`, sau do dung phuong phap luy thua va xuong thang.

## Chay giao dien

```bash
python svd_gui.py
```

Trong giao dien:

- Nhap ma tran va tham so truc tiep.
- Hoac bam `Mo input bang Notepad`, sua file `svd_input.txt`, roi bam `Doc tu file input`.
- Bam `Xem truoc` de xem Markdown.
- Bam `Giai va xuat .md` de tao loi giai.

## Chay bang dong lenh

```bash
python svd_gui.py --input svd_input.txt --output svd_solution.md
```

## Dinh dang input

```text
matrix:
3 0
4 0
0 2

tol: 1e-7
rank_tol: 1e-9
max_iter: 100
components: 2
mode: exam
pinv: true
condition_number: true
```

## Dang dang xu ly

- Ma tran dung `m >= n`: lap `B = A^T A`, tim `v_i`, roi tinh `u_i = A v_i / sigma_i`.
- Ma tran ngang `m < n`: lap `C = A A^T`, tim `u_i`, roi tinh `v_i = A^T u_i / sigma_i`.
- SVD rut gon: `A = U_r Sigma_r V_r^T`.
- Nghich dao suy rong: `A^dagger = V_r Sigma_r^{-1} U_r^T`.
- So dieu kien: `cond(A) = sigma_max / sigma_min`.

## Kiem tra trong output

- Vector rieng dung de xuong thang duoc chuan hoa theo chuan 2: `||q_i||_2 = 1`.
- Sai so du cua cap tri rieng la `r_i = ||M q_i - lambda_i q_i||_2`.
- Sai so tai tao la `E = ||A - U_r Sigma_r V_r^T||_F`.
- Neu `mode: full`, output hien them `U_r^T U_r` va `V_r^T V_r`.
- Neu lay du so thanh phan can thiet, output dung dau `=`; neu chi lay mot phan, output dung dau `approx`.
- Neu `components < min(m,n)`, chuong trinh khong ket luan so dieu kien vi chua biet gia tri ky di nho nhat.

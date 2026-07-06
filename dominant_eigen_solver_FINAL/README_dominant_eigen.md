# Dominant Eigenvalue Solver

Tool giai bai toan tim gia tri rieng troi bang phuong phap luy thua va xuong thang.

## Chay giao dien

```bash
python dominant_eigen_gui.py
```

Trong giao dien:

- Nhap ma tran va tham so truc tiep.
- Hoac bam `Mo input bang Notepad`, sua file `dominant_input.txt`, roi bam `Doc tu file input`.
- Bam `Xem truoc` de xem Markdown.
- Bam `Giai va xuat .md` de tao file loi giai.

## Chay bang dong lenh

```bash
python dominant_eigen_gui.py --input dominant_input.txt --output dominant_solution.md
```
## Đơn giản hơn thì chỉ cần vào file có chữ gui rồi ấn run là được ##
## Dinh dang input

```text
matrix:
2 1 0
1 2 1
0 1 2

x0:
1 1 1

tol: 1e-7
max_iter: 30
deflation_count: 2
mode: exam
show_special_checks: auto
```

Trong do:

- `matrix`: ma tran vuong A.
- `x0`: vector ban dau. Neu bo qua, chuong trinh dung vector toan 1.
- `tol`: sai so dung.
- `max_iter`: so lan lap toi da.
- `deflation_count`: so gia tri rieng muon xuong thang. Dat 1 neu chi can gia tri rieng troi dau tien.
- `mode`: `exam` de xuat loi giai gon kieu bai thi, `full` de hien them kiem tra phu.
- `show_special_checks`: `auto`, `always` hoac `never`.

## Quy tac lam tron

- Phep tinh ben trong dung gia tri goc.
- Ma tran trong output lam tron 4 chu so thap phan.
- Cac so lieu khac lam tron 7 chu so thap phan.

## Cac truong hop dang xu ly

Chuong trinh kiem tra theo dung thu tu:

```text
TH1 -> TH2 -> TH3
```

- `TH1`: day luy thua hoi tu truc tiep. Output ket luan gia tri rieng troi va vector rieng tu bang lap.
- `TH2`: day khong hoi tu truc tiep nhung co chu ky 2, va uoc luong `lambda^2` la so thuc duong. Output dung luy thua chan de suy ra cap `lambda` va `-lambda`.
- `TH3`: TH1 va TH2 khong thoa. Output lap he cu the `z_{n+2}-p z_{n+1}+q z_n=0`, giai `p,q`, roi giai phuong trinh `t^2 - pt + q = 0`.

Khi kiem tra TH2, chuong trinh chi so sanh cac thanh phan troi cua vector de tranh viec cac thanh phan rat nho lam sai ket luan.

## Vector rieng, sai so va xuong thang

- Cot `Sai so` trong bang lap la chuan vector du `||Ax - lambda x||_2`.
- Output co ghi cong thuc `E1` cho TH1 va `E2` cho TH2 truoc khi thay so.
- Voi TH2 va TH3, chuong trinh tinh vector rieng bang cach giai he `(A - lambda I)v = 0`.
- Neu `deflation_count > 1`, chuong trinh thuc hien xuong thang bang cong thuc:

```text
A1 = A - lambda/(w^T v) v w^T
```

trong do `v` la vector rieng phai va `w` la vector rieng trai.

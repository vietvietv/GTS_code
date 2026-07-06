# Danielevsky Solver

Tool nho de giai bai toan Danielevsky va xuat loi giai dang Markdown.

## Cach chay giao dien

```bash
python danielevsky_gui.py
```

Trong giao dien:

- Nhap ma tran truc tiep vao o ben trai.
- Hoac bam `Mo input bang Notepad`, sua file `matrix_input.txt`, sau do bam `Doc tu file input`.
- Bam `Xem truoc` de xem Markdown.
- Bam `Giai va xuat .md` de tao file loi giai.

## Cach chay bang dong lenh

```bash
python danielevsky_gui.py --input matrix_input.txt --output danielevsky_solution.md
```

## Dinh dang input

Nhap tung hang cua ma tran tren mot dong:

```text
4 -1 -1
0 1 0
2 2 1
```

Co the dung dau cach, dau phay hoac dau cham phay de tach so.

## Quy tac lam tron

- Phep tinh ben trong dung gia tri goc, khong dung gia tri da lam tron.
- Ma tran trong output lam tron 4 chu so thap phan.
- Cac so lieu khac lam tron 7 chu so thap phan.

## Luu y

Ban hien tai xu ly cac truong hop:

- TH1: phan tu tru khac 0.
- TH2: phan tu tru bang 0 nhung ben trai con phan tu khac 0, can hoan vi.
- TH3: hang dang xet co phan ben trai bang 0, can tach khoi Frobenius.

Co the thu nhanh TH3 bang file:

```bash
python danielevsky_gui.py --input matrix_input_th3.txt --output danielevsky_solution_th3.md
```

import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import io

# ==========================================
# PHẦN 1: THUẬT TOÁN GAUSS-SEIDEL (DỪNG BẰNG HẬU NGHIỆM)
# ==========================================

def matrix_to_latex_core(np_mat, decimals=None):
    """Chuyển ma trận numpy thành chuỗi LaTeX pmatrix"""
    mat = np.copy(np_mat).astype(float)
    mat[np.abs(mat) < 1e-10] = 0.0
    
    latex_str = [r"\begin{pmatrix}"]
    if mat.ndim == 1:
        if decimals is not None:
            latex_str.append(" \\\\ ".join([f"{val:.{decimals}f}" for val in mat]) + r" \\")
        else:
            latex_str.append(" \\\\ ".join([f"{val:g}" for val in mat]) + r" \\")
    else:
        for row in mat:
            if decimals is not None:
                row_str = " & ".join([f"{val:.{decimals}f}" for val in row])
            else:
                row_str = " & ".join([f"{val:g}" for val in row])
            latex_str.append("  " + row_str + r" \\")
    latex_str.append(r"\end{pmatrix}")
    return "\n".join(latex_str)

def format_vec(vec, decimals):
    """Format vector nghiệm với số chữ số làm tròn cố định"""
    return r"\begin{pmatrix}" + " \\\\ ".join([f"{val:.{decimals}f}" for val in vec]) + r" \\ \end{pmatrix}"

def inline_vec(vec, decimals=None):
    """Format vector thành dạng (x1, x2, ...)^T cho text inline"""
    if decimals is not None:
        elements = ", ".join([f"{val:.{decimals}f}" for val in vec])
    else:
        elements = ", ".join([f"{val:g}" for val in vec])
    return f"({elements})^T"

def raw_fraction_latex(num, den):
    """Ghép số thực thành phân số LaTeX tường minh CHƯA RÚT GỌN (Giữ nguyên dấu vết)"""
    if abs(num) < 1e-10:
        return "0"
    n_val = float(num)
    d_val = float(den)
    
    # Xử lý dấu trừ cho đẹp
    if d_val < 0:
        n_val = -n_val
        d_val = -d_val
        
    if n_val < 0:
        return f"-\\frac{{{-n_val:g}}}{{{d_val:g}}}"
    else:
        return f"\\frac{{{n_val:g}}}{{{d_val:g}}}"

def str_matrix_to_latex(str_mat):
    """Chuyển ma trận chứa string phân số thành chuỗi LaTeX"""
    latex_str = [r"\begin{pmatrix}"]
    if isinstance(str_mat[0], str): 
        latex_str.append(" \\\\ ".join(str_mat) + r" \\")
    else:
        for row in str_mat:
            latex_str.append("  " + " & ".join(row) + r" \\")
    latex_str.append(r"\end{pmatrix}")
    return "\n".join(latex_str)

def run_gauss_seidel(aug_mat, eps, x0, decimals, n_cols_b):
    m, total_cols = aug_mat.shape
    n = total_cols - n_cols_b
    
    if m != n:
        return f"# LỖI: Ma trận hệ số A phải là ma trận vuông. Kích thước hiện tại của A là {m}x{n}."
        
    A = aug_mat[:, :n]
    B = aug_mat[:, n:]
    
    md = []
    md.append("# CHƯƠNG 2. Phương pháp lặp Gauss-Seidel\n")
    md.append(f"Lưu ý: Các số liệu trong bài làm lấy đến {decimals} chữ số thập phân, riêng các ma trận chỉ cần ghi đến 4 chữ số thập phân.\n")
    
    # --- KIỂM TRA ĐIỀU KIỆN HỘI TỤ ---
    is_row_dominant = True
    is_col_dominant = True
    
    row_check_str = []
    for i in range(n):
        diag = abs(A[i, i])
        off_diag = sum(abs(A[i, j]) for j in range(n) if i != j)
        status = "(Tm)" if diag > off_diag else "(KTM)"
        if diag <= off_diag: is_row_dominant = False
        sum_str = " + ".join([f"|{A[i,j]:g}|" for j in range(n) if i != j])
        row_check_str.append(f"* Hàng {i+1}: $|{A[i,i]:g}| > {sum_str} = {off_diag:g}$ {status}")
        
    col_check_str = []
    for j in range(n):
        diag = abs(A[j, j])
        off_diag = sum(abs(A[i, j]) for i in range(n) if i != j)
        status = "(Tm)" if diag > off_diag else "(KTM)"
        if diag <= off_diag: is_col_dominant = False
        sum_str = " + ".join([f"|{A[i,j]:g}|" for i in range(n) if i != j])
        col_check_str.append(f"* Cột {j+1}: $|{A[j,j]:g}| > {sum_str} = {off_diag:g}$ {status}")

    mode = ""
    title_suffix = ""
    if is_row_dominant:
        mode = "row"
        title_suffix = "Chéo trội hàng"
    elif is_col_dominant:
        mode = "col"
        title_suffix = "Chéo trội cột"
    else:
        return "# LỖI: Ma trận không chéo trội hàng cũng không chéo trội cột ngặt.\nThuật toán không đảm bảo hội tụ theo chuẩn cơ bản."

    for k in range(n_cols_b):
        b = B[:, k]
        
        md.append(f"\n## 2.1 Ví dụ {k+1}: {title_suffix}\n")
        md.append("### 2.1.1 Điều kiện hội tụ\n")
        md.append(f"Xét $AX = b$:\n$$ A = {matrix_to_latex_core(A)}, \\quad b = {matrix_to_latex_core(b)} $$\n")
        
        if mode == "row":
            md.append("Kiểm tra đk chéo trội hàng $A$:\n\n" + "\n".join(row_check_str) + "\n\n")
            md.append(f"Ma trận $A$ chéo trội hàng, PP Gauss-Seidel hội tụ với mọi xấp xỉ ban đầu $X^{(0)}$.\n")
        else:
            md.append("KT chéo trội hàng (KTM). Ta kiểm tra đk chéo trội cột ngặt:\n\n" + "\n".join(col_check_str) + "\n\n")
            md.append(f"Ma trận $A$ chéo trội cột ngặt, PP Gauss-Seidel hội tụ với mọi xấp xỉ ban đầu $X^{(0)}$.\n")

        # --- TẠO MA TRẬN B1, B2, d TƯỜNG MINH CHƯA RÚT GỌN ---
        B1_str_mat = []
        B2_str_mat = []
        d_str_vec = []
        for i in range(n):
            d_str_vec.append(raw_fraction_latex(b[i], A[i, i]))
            row_b1 = []
            row_b2 = []
            for j in range(n):
                if j < i:
                    row_b1.append(raw_fraction_latex(-A[i, j], A[i, i]))
                    row_b2.append("0")
                elif j > i:
                    row_b1.append("0")
                    row_b2.append(raw_fraction_latex(-A[i, j], A[i, i]))
                else:
                    row_b1.append("0")
                    row_b2.append("0")
            B1_str_mat.append(row_b1)
            B2_str_mat.append(row_b2)

        md.append("\n### 2.1.2 Ma trận lặp và số vòng lặp (Tiên nghiệm)\n")
        md.append("CT lặp $X^{(k+1)} = B_1 X^{(k+1)} + B_2 X^{(k)} + d$:\n")
        md.append(f"$$ B_1 = {str_matrix_to_latex(B1_str_mat)}, \\quad B_2 = {str_matrix_to_latex(B2_str_mat)}, \\quad d = {str_matrix_to_latex(d_str_vec)} $$\n")
        
        # --- TÍNH MU VÀ S ---
        if mode == "row":
            mu_vals = []
            for i in range(n):
                sum_lt = sum(abs(A[i, j]) for j in range(i))       
                sum_gt = sum(abs(A[i, j]) for j in range(i+1, n))  
                mu_vals.append(sum_lt / (abs(A[i, i]) - sum_gt))
            mu = max(mu_vals)
            s = 0 
            norm_sym = r"\infty"
            md.append("Dùng công thức hệ số co (chuẩn vô cùng):\n")
            md.append(f"$$ \\mu = \\max_i \\left\\{{ \\frac{{\\sum_{{j<i}} |a_{{ij}}|}}{{|a_{{ii}}| - \\sum_{{j>i}} |a_{{ij}}|}} \\right\\}} = {mu:g} < 1 $$\n")
        else:
            mu_vals = []
            s_vals = []
            for j in range(n):
                sum_lt = sum(abs(A[i, j]) for i in range(j))       
                sum_gt = sum(abs(A[i, j]) for i in range(j+1, n))  
                mu_vals.append(sum_lt / (abs(A[j, j]) - sum_gt))
                s_vals.append(sum_gt / abs(A[j, j]))
            mu = max(mu_vals)
            s = max(s_vals)
            norm_sym = "1"
            md.append("Dùng công thức hệ số co (chuẩn 1):\n")
            md.append(f"$$ \\mu = \\max_j \\left\\{{ \\frac{{\\sum_{{i<j}} |a_{{ij}}|}}{{|a_{{jj}}| - \\sum_{{i>j}} |a_{{ij}}|}} \\right\\}} = {mu:g} < 1 $$\n")
            md.append(f"$$ s = \\max_j \\left\\{{ \\frac{{\\sum_{{i>j}} |a_{{ij}}|}}{{|a_{{jj}}|}} \\right\\}} = {s:g} $$\n")

        # --- TÍNH X(1) TIÊN NGHIỆM BẰNG JACOBI ---
        x1_tien_nghiem = np.zeros(n)
        for i in range(n):
            sum_val = b[i]
            for j in range(n):
                if i != j: sum_val -= A[i, j] * x0[j]
            x1_tien_nghiem[i] = sum_val / A[i, i]
            
        norm_type = np.inf if mode == "row" else 1
        delta_x0 = np.linalg.norm(x1_tien_nghiem - x0, norm_type)
        
        md.append(f"Chọn $X^{(0)} = {inline_vec(x0)} \\implies X^{(1)} = {inline_vec(x1_tien_nghiem, decimals)}$.\n")
        md.append(f"Độ lệch bước đầu: $\\|X^{(1)} - X^{(0)}\\|_{norm_sym} = {delta_x0:.{decimals}f}$.\n")
        
        # Chỉ giữ lại tiên nghiệm để tham khảo
        md.append(f"\n*Dự đoán* số bước lặp với $\\epsilon = 10^{{{int(np.log10(eps))}}}$:\n")
        if mode == "row":
            n_steps = np.log(eps * (1 - mu) / delta_x0) / np.log(mu)
            md.append(f"$$ \\frac{{\\mu^n}}{{1 - \\mu}} \\|X^{(1)} - X^{(0)}\\|_{norm_sym} \\le \\epsilon \\implies n \\ge {n_steps:.4f} $$\n")
        else:
            n_steps = np.log(eps * (1 - s) * (1 - mu) / delta_x0) / np.log(mu)
            md.append(f"$$ \\frac{{\\mu^n}}{{(1-s)(1-\\mu)}} \\|X^{(1)} - X^{(0)}\\|_{norm_sym} \\le \\epsilon \\implies n \\ge {n_steps:.4f} $$\n")

        # --- QUÁ TRÌNH LẶP (Gauss-Seidel) ---
        md.append("\n### 2.1.3 Quá trình lặp\n")
        md.append(f"Ta thực hiện lặp và sẽ dừng thuật toán khi **sai số hậu nghiệm** $\\le \\epsilon = {eps}$:\n")
        
        # Tính x1 của Gauss-Seidel để dùng cho công thức mẫu bước 2
        x1_gs = np.zeros(n)
        for i in range(n):
            sum_val = b[i]
            for j in range(n):
                if j < i: sum_val -= A[i, j] * x1_gs[j]
                elif j > i: sum_val -= A[i, j] * x0[j]
            x1_gs[i] = sum_val / A[i, i]

        # Khai báo công thức bước 2 tường minh
        x2_gs = np.zeros(n)
        md.append("Lặp thứ 2 ($X^{(2)}$):\n$$\n\\begin{aligned}\n")
        for i in range(n):
            sum_val = b[i]
            eq_str = f"x^{{(2)}}_{{{i+1}}} &= \\frac{{1}}{{{A[i,i]:g}}} \\left( {b[i]:g} "
            for j in range(n):
                if j < i: 
                    sum_val -= A[i, j] * x2_gs[j]
                    eq_str += f" - ({A[i,j]:g})({x2_gs[j]:.{decimals}f})"
                elif j > i: 
                    sum_val -= A[i, j] * x1_gs[j]
                    eq_str += f" - ({A[i,j]:g})({x1_gs[j]:.{decimals}f})"
            x2_gs[i] = sum_val / A[i, i]
            eq_str += f" \\right) = {x2_gs[i]:.{decimals}f} \\\\\n"
            md.append(eq_str)
        md.append("\\end{aligned}\n$$\n")

        md.append(f"\nBảng lặp (làm tròn đến {decimals} chữ số thập phân):\n")
        md.append("| $n$ | " + " | ".join([f"$x_{i+1}^{(n)}$" for i in range(n)]) + f" | $\\|\\Delta X^{(n)}\\|_{norm_sym}$ |\n")
        md.append("|:---:|" + ":---:|" * (n + 1) + "\n")
        
        curr_x = x0.copy()
        md.append(f"| 0 | " + " | ".join([f"{v:.{decimals}f}" for v in curr_x]) + " | - |\n")
        
        # VÒNG LẶP ĐỘNG: KIỂM TRA SAI SỐ HẬU NGHIỆM Ở MỖI BƯỚC
        step = 1
        final_err = 0
        final_diff = 0
        MAX_ITER = 1000 # Chặn an toàn
        
        while step <= MAX_ITER:
            prev_x = curr_x.copy()
            for i in range(n):
                sum_val = b[i]
                for j in range(n):
                    if j < i: sum_val -= A[i, j] * curr_x[j]
                    elif j > i: sum_val -= A[i, j] * prev_x[j]
                curr_x[i] = sum_val / A[i, i]
                
            diff = np.linalg.norm(curr_x - prev_x, norm_type)
            md.append(f"| {step} | " + " | ".join([f"{v:.{decimals}f}" for v in curr_x]) + f" | {diff:.{decimals}f} |\n")
            
            # Tính sai số hậu nghiệm hiện tại
            if mode == "row":
                err_post = (mu / (1 - mu)) * diff
            else:
                err_post = (mu / ((1 - s) * (1 - mu))) * diff
                
            # ĐIỀU KIỆN DỪNG THEO SAI SỐ HẬU NGHIỆM
            if err_post <= eps:
                final_err = err_post
                final_diff = diff
                break
            step += 1

        # --- ĐÁNH GIÁ HẬU NGHIỆM ---
        md.append("\n### 2.1.4 Đánh giá sai số hậu nghiệm\n")
        md.append(f"Tại bước lặp $n = {step}$, độ lệch giữa hai bước lặp cuối cùng là:\n")
        
        diff_vec = curr_x - prev_x
        md.append(f"$$ X^{({step})} - X^{({step-1})} = {format_vec(diff_vec, decimals)} \\implies \\|X^{({step})} - X^{({step-1})}\|_{norm_sym} = {final_diff:.{decimals}f} $$\n")
        
        if mode == "row":
            md.append(f"Kiểm tra sai số hậu nghiệm:\n$$ \\|X^{({step})} - X^*\\|_{norm_sym} \\le \\frac{{\\mu}}{{1 - \\mu}} \\|X^{({step})} - X^{({step-1})}\|_{norm_sym} = \\frac{{{mu:g}}}{{{1-mu:g}}} \\times {final_diff:.{decimals}f} \\approx {final_err:.8f} $$\n")
        else:
            md.append(f"Kiểm tra sai số hậu nghiệm:\n$$ \\|X^{({step})} - X^*\\|_{norm_sym} \\le \\frac{{\\mu}}{{(1-s)(1-\\mu)}} \\|X^{({step})} - X^{({step-1})}\|_{norm_sym} = \\frac{{{mu:g}}}{{({1-s:g})({1-mu:g})}} \\times {final_diff:.{decimals}f} \\approx {final_err:.8f} $$\n")

        md.append(f"Vì ${final_err:.8f} \\le {eps}$, thỏa mãn điều kiện dừng. \nNghiệm xấp xỉ của hệ là $X \\approx {inline_vec(curr_x, decimals)}$.\n")

    return "".join(md)

# ==========================================
# PHẦN 2: GIAO DIỆN TKINTER
# ==========================================

class MatrixApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Gauss-Seidel Dừng Bằng Hậu Nghiệm - Toán Tin BKHN")
        self.root.geometry("650x650")
        
        config_frame = ttk.LabelFrame(root, text="1. Cấu hình thông số", padding=10)
        config_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Label(config_frame, text="Số cột ma trận b:").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.cols_b_var = tk.IntVar(value=1)
        ttk.Entry(config_frame, textvariable=self.cols_b_var, width=10).grid(row=0, column=1, sticky="w")
        
        ttk.Label(config_frame, text="Sai số (eps):").grid(row=0, column=2, padx=5, pady=5, sticky="e")
        self.eps_var = tk.DoubleVar(value=0.000001)
        ttk.Entry(config_frame, textvariable=self.eps_var, width=15).grid(row=0, column=3, sticky="w")
        
        ttk.Label(config_frame, text="Số chữ số làm tròn:").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.decimals_var = tk.IntVar(value=6)
        ttk.Entry(config_frame, textvariable=self.decimals_var, width=10).grid(row=1, column=1, sticky="w")
        
        ttk.Label(config_frame, text="Vector X(0):").grid(row=1, column=2, padx=5, pady=5, sticky="e")
        self.x0_var = tk.StringVar(value="1 1 1 1")
        ttk.Entry(config_frame, textvariable=self.x0_var, width=15).grid(row=1, column=3, sticky="w")

        input_frame = ttk.LabelFrame(root, text="2. Nhập ma trận mở rộng [A | b]", padding=10)
        input_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        toolbar = ttk.Frame(input_frame)
        toolbar.pack(fill="x", pady=(0, 5))
        ttk.Button(toolbar, text="Mở file .txt", command=self.load_file).pack(side="left")
        
        self.text_area = tk.Text(input_frame, wrap="none", font=("Consolas", 11))
        self.text_area.pack(fill="both", expand=True)
        # Ma trận mẫu chéo trội hàng với số thập phân
        sample_data = "18 -2 3 -5 24.5\n3 24 -4 4 40.3\n2 7 -15 0 -6.2\n-1 2 -12 32 -47.2"
        self.text_area.insert("1.0", sample_data) 
        
        btn_frame = ttk.Frame(root, padding=10)
        btn_frame.pack(fill="x")
        ttk.Button(btn_frame, text="GIẢI GAUSS-SEIDEL & XUẤT MD", command=self.solve).pack(pady=10)

    def load_file(self):
        filepath = filedialog.askopenfilename(filetypes=[("Text Files", "*.txt")])
        if filepath:
            with open(filepath, 'r') as f:
                self.text_area.delete("1.0", tk.END)
                self.text_area.insert(tk.END, f.read())

    def solve(self):
        try:
            n_cols_b = self.cols_b_var.get()
            if n_cols_b <= 0:
                messagebox.showerror("Lỗi", "Số cột của b phải lớn hơn 0!")
                return
                
            raw_text = self.text_area.get("1.0", tk.END).strip()
            aug_mat = np.loadtxt(io.StringIO(raw_text))
            if aug_mat.ndim == 1: aug_mat = aug_mat.reshape(1, -1)
            
            x0_str = self.x0_var.get().strip()
            x0 = np.array([float(x) for x in x0_str.split()])
            
            total_cols = aug_mat.shape[1]
            n = total_cols - n_cols_b
            
            if len(x0) != n:
                messagebox.showerror("Lỗi", f"Số lượng phần tử của X0 ({len(x0)}) phải bằng số ẩn n ({n})!\nVui lòng nhập lại X0.")
                return
            
            save_path = filedialog.asksaveasfilename(defaultextension=".md", filetypes=[("Markdown", "*.md")])
            if save_path:
                md_content = run_gauss_seidel(
                    aug_mat=aug_mat, 
                    eps=self.eps_var.get(), 
                    x0=x0, 
                    decimals=self.decimals_var.get(),
                    n_cols_b=n_cols_b
                )
                with open(save_path, 'w', encoding='utf-8') as f:
                    f.write(md_content)
                messagebox.showinfo("Thành công", "Bảng lặp tự động dừng khi thỏa mãn sai số hậu nghiệm!")
        except Exception as e:
            messagebox.showerror("Lỗi", f"Đã xảy ra lỗi: {str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = MatrixApp(root)
    root.mainloop()
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import io

# ==========================================
# THUẬT TOÁN: JACOBI GIẢI HỆ AX = B
# ==========================================

def format_num(val, decimals=4):
    """Làm tròn số thập phân, tự động bỏ các số 0 thừa ở đuôi"""
    if abs(val) < 1e-15: return "0"
    if abs(val - round(val)) < 1e-15: return f"{int(round(val))}"
    return f"{val:.{decimals}f}".rstrip('0').rstrip('.')

def matrix_to_latex(np_mat, decimals=4):
    """Chuyển ma trận 1D/2D sang LaTeX matrix"""
    mat = np.array(np_mat)
    if mat.ndim == 1:
        mat = mat.reshape(-1, 1)
        
    latex_str = [r"\begin{pmatrix}"]
    for row in mat:
        row_str = " & ".join([format_num(val, decimals) for val in row])
        latex_str.append("  " + row_str + r" \\")
    latex_str.append(r"\end{pmatrix}")
    return "\n".join(latex_str)

def run_jacobi_solver(A, B, X0, epsilon, decimals):
    n = A.shape[0]
    if A.shape[0] != A.shape[1]:
        return "# LỖI: Ma trận A phải là ma trận vuông."
    
    # Đưa B và X0 về dạng vector cột
    B_mat = np.array(B)
    if B_mat.ndim == 1: B_mat = B_mat.reshape(-1, 1)
        
    X0_mat = np.array(X0)
    if X0_mat.ndim == 1: X0_mat = X0_mat.reshape(-1, 1)
        
    if B_mat.shape[0] != n or X0_mat.shape[0] != n:
        return f"# LỖI: Kích thước của B và X(0) phải khớp với cấp của A ({n})."
        
    md = []
    md.append("## Thuật toán: Jacobi giải hệ phương trình $AX = B$\n\n")
    md.append("**Đầu vào (Input):**\n")
    md.append(f"- Ma trận vuông $A$ cấp $n={n}$:\n$$ A = {matrix_to_latex(A, decimals)} $$\n")
    md.append(f"- Ma trận $B$:\n$$ B = {matrix_to_latex(B_mat, decimals)} $$\n")
    md.append(f"- Sai số cho trước: $\\epsilon = {epsilon}$\n")
    md.append(f"- Véc-tơ xấp xỉ ban đầu $X^{(0)}$:\n$$ X^{(0)} = {matrix_to_latex(X0_mat, decimals)} $$\n\n")
    md.append("**Các bước thực hiện:**\n\n")
    
    # Kiểm tra phần tử trên đường chéo bằng 0
    diag_A = np.diag(A)
    if np.any(diag_A == 0):
        md.append("Phát hiện phần tử trên đường chéo chính bằng 0. Không thể áp dụng thuật toán Jacobi.\n")
        return "".join(md)

    # 1. Kiểm tra tính chéo trội
    md.append("### 1. Kiểm tra tính chéo trội của ma trận $A$ và thực hiện gán.\n")
    
    p, lam, q = None, None, None
    norm_symbol = ""
    
    # Tính q cho chéo trội hàng
    sum_row_off_diag = np.sum(np.abs(A), axis=1) - np.abs(diag_A)
    q_row_arr = sum_row_off_diag / np.abs(diag_A)
    q_row = np.max(q_row_arr)
    
    # Tính q cho chéo trội cột
    sum_col_off_diag = np.sum(np.abs(A), axis=0) - np.abs(diag_A)
    q_col_arr = sum_col_off_diag / np.abs(diag_A)
    q_col = np.max(q_col_arr)

    md.append(f"Kiểm tra chéo trội hàng: $\\max_{{1 \\le i \\le n}} \\frac{{1}}{{|a_{{ii}}|}} \\sum_{{j \\neq i}} |a_{{ij}}| = {format_num(q_row, decimals)}$.\n\n")
    
    if q_row < 1:
        p = np.inf
        lam = 1.0
        q = q_row
        norm_symbol = "\\infty"
        md.append("Vì ma trận $A$ chéo trội hàng, ta đặt:\n")
        md.append(f"$$ p = \\infty, \\quad \\lambda = 1, \\quad q = {format_num(q, decimals)} $$\n\n")
    else:
        md.append(f"Ma trận $A$ không chéo trội hàng ($q \\ge 1$). Ta kiểm tra chéo trội cột:\n")
        md.append(f"$$ \\max_{{1 \\le j \\le n}} \\frac{{1}}{{|a_{{jj}}|}} \\sum_{{i \\neq j}} |a_{{ij}}| = {format_num(q_col, decimals)} $$\n\n")
        
        if q_col < 1:
            p = 1
            min_diag = np.min(np.abs(diag_A))
            max_diag = np.max(np.abs(diag_A))
            lam = max_diag / min_diag  # Sửa lại thành max/min theo đúng chuẩn
            q = q_col
            norm_symbol = "1"
            md.append("Vì ma trận $A$ chéo trội cột, ta đặt:\n")
            md.append(f"$$ p = 1, \\quad \\lambda = \\frac{{\\max |a_{{ii}}|}}{{\\min |a_{{ii}}|}} = \\frac{{{format_num(max_diag, decimals)}}}{{{format_num(min_diag, decimals)}}} = {format_num(lam, decimals)}, \\quad q = {format_num(q, decimals)} $$\n\n")
        else:
            md.append(f"**DỪNG CHƯƠNG TRÌNH:** Ma trận $A$ không chéo trội hàng cũng không chéo trội cột ($q \\ge 1$).\n")
            return "".join(md)

    # 2. Khởi tạo
    md.append("### 2. Khởi tạo.\n")
    eps_0 = ((1 - q) / (lam * q)) * epsilon
    md.append(f"Đặt $k = 1$ và tính sai số dừng $\\epsilon_0$:\n")
    md.append(f"$$ \\epsilon_0 = \\frac{{1 - q}}{{\\lambda q}} \\epsilon = \\frac{{1 - {format_num(q, decimals)}}}{{({format_num(lam, decimals)})({format_num(q, decimals)})}} \\times {epsilon} = {format_num(eps_0, decimals)} $$\n\n")
    
    T_diag = 1.0 / diag_A
    T = np.diag(T_diag)
    md.append("Khởi tạo $I$ là ma trận đơn vị cấp $n$ và ma trận $T$:\n")
    md.append(f"$$ T = \\text{{diag}}\\left(\\frac{{1}}{{a_{{11}}}}, \\dots, \\frac{{1}}{{a_{{nn}}}}\\right) = {matrix_to_latex(T, decimals)} $$\n\n")

    # 3. Tính ma trận lặp
    md.append("### 3. Tính ma trận lặp.\n")
    I = np.eye(n)
    C = I - T @ A
    D = T @ B_mat
    md.append(f"$$ C = I - TA = {matrix_to_latex(C, decimals)} $$\n")
    md.append(f"$$ D = TB = {matrix_to_latex(D, decimals)} $$\n\n")

    # 4. Thực hiện lặp
    md.append("### 4. Thực hiện lặp.\n")
    md.append("Ở mỗi bước, tính $X^{(k)} = CX^{(k-1)} + D$. Ta lập bảng lặp như sau:\n\n")
    
    # Bảng lặp
    header = ["$k$"]
    for i in range(n):
        header.append(f"$x_{{{i+1}}}$")
    header.append(f"$\\|X^{{(k)}} - X^{{(k-1)}}\\|_{norm_symbol}$")
    
    md.append("| " + " | ".join(header) + " |\n")
    md.append("|" + "|".join(["---"] * len(header)) + "|\n")
    
    # Dòng k = 0
    row_0 = ["0"]
    for val in X0_mat.flatten():
        row_0.append(format_num(val, decimals))
    row_0.append("-")
    md.append("| " + " | ".join(row_0) + " |\n")
    
    k = 1
    X_prev = X0_mat.copy()
    
    while True:
        X_k = C @ X_prev + D
        diff = X_k - X_prev
        norm_diff = np.linalg.norm(diff, ord=p)
        
        row_k = [str(k)]
        for val in X_k.flatten():
            row_k.append(format_num(val, decimals))
        row_k.append(format_num(norm_diff, decimals))
        md.append("| " + " | ".join(row_k) + " |\n")
        
        if norm_diff <= eps_0:
            md.append("\n")
            md.append(f"Tại bước lặp $k = {k}$, ta có $\\|X^{{({k})}} - X^{{({k-1})}}\\|_{norm_symbol} = {format_num(norm_diff, decimals)} \\le \\epsilon_0 = {format_num(eps_0, decimals)}$.\n\n")
            md.append(f"**Kết luận:** Dừng chương trình và in ra nghiệm xấp xỉ:\n")
            md.append(f"$$ X \\approx X^{{({k})}} = {matrix_to_latex(X_k, decimals)} $$\n")
            break
            
        X_prev = X_k
        k += 1
        
        if k > 500:
            md.append("\n**CẢNH BÁO:** Quá 500 vòng lặp chưa hội tụ. Thuật toán tự dừng.\n")
            break

    return "".join(md)

# ==========================================
# GIAO DIỆN TKINTER
# ==========================================

class JacobiApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Thuật toán Jacobi giải hệ AX = B")
        self.root.geometry("750x700")
        
        # 1. Cấu hình
        config_frame = ttk.LabelFrame(root, text="1. Cấu hình thông số", padding=10)
        config_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Label(config_frame, text="Sai số (ε):").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.epsilon_var = tk.DoubleVar(value=0.001)
        ttk.Entry(config_frame, textvariable=self.epsilon_var, width=15).grid(row=0, column=1, sticky="w")
        
        ttk.Label(config_frame, text="Số chữ số làm tròn:").grid(row=0, column=2, padx=15, pady=5, sticky="e")
        self.decimals_var = tk.IntVar(value=4)
        ttk.Entry(config_frame, textvariable=self.decimals_var, width=10).grid(row=0, column=3, sticky="w")

        # Layout ma trận
        matrices_frame = ttk.Frame(root)
        matrices_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        # 2. Input A
        frame_a = ttk.LabelFrame(matrices_frame, text="2. Ma trận vuông A", padding=10)
        frame_a.pack(side="left", fill="both", expand=True, padx=(0, 5))
        self.text_a = tk.Text(frame_a, wrap="none", font=("Consolas", 11), width=20)
        self.text_a.pack(fill="both", expand=True)
        
        # 3. Input B
        frame_b = ttk.LabelFrame(matrices_frame, text="3. Ma trận/Vector B", padding=10)
        frame_b.pack(side="left", fill="both", expand=True, padx=5)
        self.text_b = tk.Text(frame_b, wrap="none", font=("Consolas", 11), width=15)
        self.text_b.pack(fill="both", expand=True)

        # 4. Input X0
        frame_x0 = ttk.LabelFrame(matrices_frame, text="4. Xấp xỉ ban đầu X(0)", padding=10)
        frame_x0.pack(side="left", fill="both", expand=True, padx=(5, 0))
        self.text_x0 = tk.Text(frame_x0, wrap="none", font=("Consolas", 11), width=15)
        self.text_x0.pack(fill="both", expand=True)

        # Dữ liệu mẫu (A chéo trội hàng)
        sample_a = (
            "10 -1 2 0\n"
            "-1 11 -1 3\n"
            "2 -1 10 -1\n"
            "0 3 -1 8"
        )
        sample_b = (
            "6\n"
            "25\n"
            "-11\n"
            "15"
        )
        sample_x0 = (
            "0\n"
            "0\n"
            "0\n"
            "0"
        )
        self.text_a.insert("1.0", sample_a) 
        self.text_b.insert("1.0", sample_b) 
        self.text_x0.insert("1.0", sample_x0) 
        
        btn_frame = ttk.Frame(root, padding=10)
        btn_frame.pack(fill="x")
        ttk.Button(btn_frame, text="GIẢI HỆ & XUẤT MARKDOWN", command=self.solve).pack(pady=10)

    def solve(self):
        try:
            raw_a = self.text_a.get("1.0", tk.END).strip()
            raw_b = self.text_b.get("1.0", tk.END).strip()
            raw_x0 = self.text_x0.get("1.0", tk.END).strip()
            
            if not raw_a or not raw_b or not raw_x0:
                messagebox.showwarning("Cảnh báo", "Vui lòng nhập đầy đủ ma trận A, B và X(0)!")
                return
                
            A = np.loadtxt(io.StringIO(raw_a))
            B = np.loadtxt(io.StringIO(raw_b))
            X0 = np.loadtxt(io.StringIO(raw_x0))
            
            if A.ndim == 1:
                A = A.reshape(1, -1)
                
            epsilon = self.epsilon_var.get()
            decimals = self.decimals_var.get()
            
            save_path = filedialog.asksaveasfilename(
                defaultextension=".md", 
                filetypes=[("Markdown", "*.md")], 
                initialfile="LoiGiai_Jacobi.md"
            )
            
            if save_path:
                md_content = run_jacobi_solver(A, B, X0, epsilon, decimals)
                with open(save_path, 'w', encoding='utf-8') as f:
                    f.write(md_content)
                messagebox.showinfo("Thành công", "Đã xuất file trình bày lời giải Markdown thành công!")
                
        except Exception as e:
            messagebox.showerror("Lỗi", f"Đã xảy ra lỗi:\n{str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = JacobiApp(root)
    root.mainloop()
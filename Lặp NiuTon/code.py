import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import io

# ==========================================
# HÀM XUẤT MA TRẬN RA LATEX
# ==========================================
def matrix_to_latex_core(np_mat, decimals=None):
    mat = np.copy(np_mat).astype(float)
    mat[np.abs(mat) < 1e-10] = 0.0
    
    latex_str = [r"\begin{bmatrix}"]
    for row in mat:
        if decimals is not None:
            row_str = " & ".join([f"{val:.{decimals}f}".rstrip('0').rstrip('.') if '.' in f"{val:.{decimals}f}" else f"{val:.{decimals}f}" for val in row])
        else:
            row_str = " & ".join([f"{val:g}" for val in row])
        latex_str.append("  " + row_str + r" \\")
    latex_str.append(r"\end{bmatrix}")
    return "\n".join(latex_str)

# ==========================================
# THUẬT TOÁN NEWTON & XUẤT FORM ĐỀ THI
# ==========================================
def run_newton_exam_style(mat, decimals, epsilon, tolerance):
    A = np.array(mat, dtype=float)
    n = A.shape[0]
    if n != A.shape[1]:
        return "# LỖI: Ma trận A phải là ma trận vuông."
        
    md = []
    
    # 1. Phần Lý thuyết & Khởi tạo 
    norm_A = np.linalg.norm(A, np.inf)
    norm_AT = np.linalg.norm(A.T, np.inf)
    
    md.append("a) Khởi tạo tham số và ma trận ban đầu:\n\n")
    md.append(f"Xét ma trận $A = {matrix_to_latex_core(A, decimals)}$. Ta có chuẩn vô cùng:\n")
    md.append(f"$$ \\|A\\|_\\infty = {norm_A:g}, \\quad \\|A^T\\|_\\infty = {norm_AT:g} $$\n\n")
    
    X0 = A.T / (norm_A * norm_AT)
    md.append("Khởi tạo ma trận xấp xỉ ban đầu để đảm bảo thuật toán hội tụ:\n")
    md.append(f"$$ X_0 = \\frac{{A^T}}{{\\|A\\|_\\infty \\cdot \\|A^T\\|_\\infty}} = {matrix_to_latex_core(X0, decimals)} $$\n\n")
    
    md.append("Từ lý thuyết phương pháp lặp Newton, ta có công thức lặp và sai số:\n")
    md.append("$$ X_{k+1} = X_k(2E - AX_k) $$\n")
    md.append("$$ \\Delta = \\|X_{k+1} - X_k\\|_\\infty \\le \\epsilon $$\n\n")
    
    # 2. Vòng lặp tính toán
    md.append(f"Ta sử dụng phương pháp lặp Newton, với điều kiện dừng $\\Delta \\le {epsilon:g}$:\n\n")
    
    X_k = X0.copy()
    E = np.eye(n)
    steps = []
    
    for k in range(1, 101):
        X_next = X_k @ (2 * E - A @ X_k)
        delta = np.linalg.norm(X_next - X_k, np.inf)
        
        steps.append({
            'k': k,
            'X': X_next,
            'delta': delta
        })
        
        if delta <= epsilon:
            break
        X_k = X_next
        
    total_steps = len(steps)
    
    for i, step in enumerate(steps):
        k = step['k']
        is_first_two = (i < 2)
        is_last_two = (i >= total_steps - 2)
        
        if is_first_two or is_last_two:
            md.append(f"---\n**🔹 Lần lặp $k = {k}$:**\n")
            md.append(f"$$ X_{{{k}}} = {matrix_to_latex_core(step['X'], decimals)} $$\n")
            md.append(f"$$ \\Delta_{{{k}}} = \\|X_{{{k}}} - X_{{{k-1}}}\\|_\\infty = {step['delta']:.{decimals}e} $$\n\n")
        elif i == 2 and total_steps > 4:
            md.append(f"---\n**🔹 Các bước lặp từ $k = 3$ đến $k = {total_steps - 2}$ tính toán tương tự...**\n\n")

    # ==========================================
    # 3. GIÁM ĐỊNH KẾT QUẢ & PHÁN QUYẾT
    # ==========================================
    last_step = steps[-1]
    residual = np.linalg.norm(E - A @ X_k, np.inf)
        
    if residual > tolerance:
        md.append(f"---\n\n**CẢNH BÁO:** Sau {last_step['k']} lần lặp, thuật toán hội tụ về ma trận giả nghịch đảo (Moore-Penrose).\n")
        md.append(f"Tuy nhiên, kiểm tra phần dư $\\|E - AX\\|_\\infty = {residual:.4f} > {tolerance:g}$.\n")
        md.append("**Kết luận:** Ma trận $A$ suy biến (không khả nghịch). Thuật toán không thể tìm được ma trận nghịch đảo thực sự.\n")
    else:
        md.append(f"---\n\nSau {last_step['k']} lần lặp, thuật toán hội tụ. Kiểm tra phần dư $\\|E - AX\\|_\\infty = {residual:.2e} \\le {tolerance:g}$, thỏa mãn điều kiện khả nghịch.\n")
        md.append(f"Nghiệm ma trận nghịch đảo $A^{{-1}} \\approx X_{{{last_step['k']}}}$ với sai số lặp $\\Delta: {last_step['delta']:.{decimals}e} \\le {epsilon:g}$.\n")
    
    return "\n".join(md)

# ==========================================
# GIAO DIỆN TKINTER
# ==========================================
class NewtonExamApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Newton - Báo Cáo Chép Tay")
        self.root.geometry("650x580")
        
        # Frame cấu hình
        config_frame = ttk.LabelFrame(root, text="1. Cấu hình tham số", padding=10)
        config_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Label(config_frame, text="Số chữ số thập phân:").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.decimals_var = tk.IntVar(value=6)
        ttk.Entry(config_frame, textvariable=self.decimals_var, width=15).grid(row=0, column=1, sticky="w")
        
        ttk.Label(config_frame, text="Sai số dừng (ε):").grid(row=0, column=2, padx=5, pady=5, sticky="e")
        self.epsilon_var = tk.StringVar(value="1e-5")
        ttk.Entry(config_frame, textvariable=self.epsilon_var, width=15).grid(row=0, column=3, sticky="w")
        
        ttk.Label(config_frame, text="Ngưỡng phần dư (tolerance):").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        self.tolerance_var = tk.StringVar(value="1e-2")
        ttk.Entry(config_frame, textvariable=self.tolerance_var, width=15).grid(row=1, column=1, sticky="w")
        
        # Frame nhập ma trận
        input_frame = ttk.LabelFrame(root, text="2. Nhập ma trận vuông A", padding=10)
        input_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        toolbar = ttk.Frame(input_frame)
        toolbar.pack(fill="x", pady=(0, 5))
        ttk.Button(toolbar, text="Mở file .txt", command=self.load_file).pack(side="left")
        
        self.text_area = tk.Text(input_frame, wrap="none", font=("Consolas", 12))
        self.text_area.pack(fill="both", expand=True)
        
        self.text_area.insert("1.0", "3 5 15 12 11 11 4\n14 9 8 15 12 1 2\n2 15 13 10 12 5 14\n14 15 3 1 6 1 16\n4 6 7 12 10 9 1\n11 1 5 1 2 13 3\n5 9 15 15 3 15 10") 
        
        # Nút tạo báo cáo
        ttk.Button(root, text="TẠO FILE MARKDOWN", command=self.solve).pack(pady=15)

    def load_file(self):
        filepath = filedialog.askopenfilename(filetypes=[("Text Files", "*.txt")])
        if filepath:
            with open(filepath, 'r') as f:
                self.text_area.delete("1.0", tk.END)
                self.text_area.insert(tk.END, f.read())

    def solve(self):
        try:
            A = np.loadtxt(io.StringIO(self.text_area.get("1.0", tk.END).strip()))
            if A.ndim == 1: 
                messagebox.showerror("Lỗi", "Hãy nhập ma trận vuông!")
                return
                
            epsilon_val = float(self.epsilon_var.get())
            tolerance_val = float(self.tolerance_var.get())
            decimals_val = self.decimals_var.get()
            
            save_path = filedialog.asksaveasfilename(defaultextension=".md", filetypes=[("Markdown", "*.md")])
            if save_path:
                md_content = run_newton_exam_style(A, decimals_val, epsilon_val, tolerance_val)
                with open(save_path, 'w', encoding='utf-8') as f:
                    f.write(md_content)
                messagebox.showinfo("Thành công", f"Đã xuất file markdown")
        except ValueError:
            messagebox.showerror("Lỗi", "Vui lòng nhập đúng định dạng số cho ma trận và sai số (VD: 1e-5 hoặc 0.00001).")
        except Exception as e:
            messagebox.showerror("Lỗi", str(e))

if __name__ == "__main__":
    root = tk.Tk()
    app = NewtonExamApp(root)
    root.mainloop()
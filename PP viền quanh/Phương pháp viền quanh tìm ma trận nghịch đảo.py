import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import io

def matrix_to_latex_core(np_mat, decimals=None):
    mat = np.copy(np_mat).astype(float)
    mat[np.abs(mat) < 1e-10] = 0.0
    
    latex_str = [r"\begin{bmatrix}"]
    if mat.ndim == 1:
        if decimals is not None:
            latex_str.append(" \\\\ ".join([f"{val:.{decimals}f}".rstrip('0').rstrip('.') if '.' in f"{val:.{decimals}f}" else f"{val:.{decimals}f}" for val in mat]) + r" \\")
        else:
            latex_str.append(" \\\\ ".join([f"{val:g}" for val in mat]) + r" \\")
    else:
        for row in mat:
            if decimals is not None:
                row_str = " & ".join([f"{val:.{decimals}f}".rstrip('0').rstrip('.') if '.' in f"{val:.{decimals}f}" else f"{val:.{decimals}f}" for val in row])
            else:
                row_str = " & ".join([f"{val:g}" for val in row])
            latex_str.append("  " + row_str + r" \\")
    latex_str.append(r"\end{bmatrix}")
    return "\n".join(latex_str)

def run_bordering_method_final(mat, decimals):
    A = np.array(mat, dtype=float)
    n = A.shape[0]
    if n != A.shape[1]:
        return "# LỖI: Ma trận A phải là ma trận vuông."
        
    md = []
    md.append("## Giải bài toán tìm ma trận nghịch đảo bằng phương pháp Viền quanh\n\n")
    md.append(f"**Ma trận đầu vào $A$:**\n$$ A = {matrix_to_latex_core(A, decimals)} $$\n\n")
    
    # BƯỚC 1: Khởi tạo M = A^T * A
    M = A.T @ A
    md.append("### Bước 1:\n")
    md.append(f"Khởi tạo ma trận $M = A^T A$:\n$$ M = {matrix_to_latex_core(M, decimals)} $$\n")
    
    m11 = M[0, 0]
    md.append(f"Kiểm tra $m_{{11}} = {m11:g}$. ")
    if abs(m11) < 1e-10:
        md.append("Do $m_{11} = 0$, thuật toán dừng (Ma trận A suy biến).\n")
        return "".join(md)
    else:
        md.append("Do $m_{11} \\neq 0$, tiếp tục sang Bước 2.\n\n")
        
    # BƯỚC 2: Khởi tạo M_1^-1
    M_inv = np.array([[1.0 / m11]])
    md.append(f"### Bước 2:\n")
    md.append(f"Khởi tạo $M_1^{{-1}} = \\left[ \\frac{{1}}{{{m11:g}}} \\right] = {matrix_to_latex_core(M_inv, decimals)}$. Gán $k=2$.\n\n")
    
    md.append(f"### Bước 3: Vòng lặp viền quanh cho ma trận $M$\n")
    
    # BƯỚC 3: Lặp k
    for k in range(2, n + 1):
        u = M[:k-1, k-1:k]  
        v = M[k-1:k, :k-1]  
        c = M[k-1, k-1]     
        
        # Biến trung gian tính ngầm để tối ưu code
        Minv_u = M_inv @ u
        v_Minv = v @ M_inv
        theta = c - (v @ Minv_u)[0, 0] 
        
        # BẮT LỖI TẠI ĐÂY: Ép hiển thị chi tiết nếu phát hiện lỗi theta = 0
        is_error = abs(theta) < 1e-10
        show_details = (k <= 3) or (k >= n - 1) or is_error
        
        if show_details:
            md.append(f"---\n**🔹 Lặp với $k = {k}$:**\n\n")
            md.append(f"*3.1. Lấy các khối từ ma trận $M$:*\n")
            md.append(f"$$ u_{k} = {matrix_to_latex_core(u, decimals)}, \\quad v_{k} = {matrix_to_latex_core(v, decimals)}, \\quad m_{{{k}{k}}} = {c:g} $$\n\n")
            
            md.append(f"*3.2. Tính đại lượng trung gian:*\n")
            md.append(f"$$ \\theta = m_{{{k}{k}}} - v_{k} M_{{{k-1}}}^{{-1}} u_{k} = {c:g} - {(v @ Minv_u)[0,0]:.{decimals}f} = {theta:.{decimals}f} $$\n\n")
            
            md.append(f"*3.3. Kiểm tra $\\theta$:* $\\theta = {theta:.{decimals}f} \\neq 0$, chuyển sang 3.4.\n\n")

        if abs(theta) < 1e-10:
            md.append(f"**LỖI:** $\\theta = 0$. Thuật toán dừng.\n")
            return "".join(md)
            
        # 3.4
        m_prime = 1.0 / theta
        u_prime = -m_prime * Minv_u
        v_prime = -m_prime * v_Minv
        B_mat = M_inv - (u_prime @ v_Minv)
        
        if show_details:
            md.append(f"*3.4. Xây dựng khối nghịch đảo mới $M_{k}^{{-1}}$:*\n")
            md.append(f"$$ m'_{{{k}{k}}} = \\frac{{1}}{{\\theta}} = {m_prime:.{decimals}f} $$\n")
            md.append(f"$$ u'_{k} = -\\frac{{1}}{{\\theta}} M_{{{k-1}}}^{{-1}} u_{k} = {matrix_to_latex_core(u_prime, decimals)} $$\n")
            md.append(f"$$ v'_{k} = -\\frac{{1}}{{\\theta}} v_{k} M_{{{k-1}}}^{{-1}} = {matrix_to_latex_core(v_prime, decimals)} $$\n")
            md.append(f"$$ B_{{{k-1}}}^{{-1}} = M_{{{k-1}}}^{{-1}} - u'_{k} v_{k} M_{{{k-1}}}^{{-1}} = {matrix_to_latex_core(B_mat, decimals)} $$\n\n")
            
        top = np.hstack((B_mat, u_prime))
        bottom = np.hstack((v_prime, [[m_prime]]))
        M_inv = np.vstack((top, bottom))
        
        if show_details:
            md.append(f"Ghép lại ta được:\n")
            md.append(f"$$ M_{k}^{{-1}} = \\begin{{bmatrix}} B_{{{k-1}}}^{{-1}} & u'_{k} \\\\ v'_{k} & m'_{{{k}{k}}} \\end{{bmatrix}} = {matrix_to_latex_core(M_inv, decimals)} $$\n\n")
        elif k == 4 and n > 4:
            md.append(f"---\n**🔹 Các bước lặp từ $k = 4$ đến $k = {n-2}$ tính toán tương tự...**\n\n")
            
    # BƯỚC 4
    Final_Inv = M_inv @ A.T
    md.append(f"### Bước 4: Xuất kết quả ma trận nghịch đảo\n")
    md.append(f"$$ A^{{-1}} = M_{n}^{{-1}} A^T = {matrix_to_latex_core(Final_Inv, decimals)} $$\n")
    
    return "".join(md)

# ==========================================
# GIAO DIỆN TKINTER
# ==========================================

class FinalBorderingApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Viền Quanh - Bản Chốt M=A^TA")
        self.root.geometry("600x500")
        
        config_frame = ttk.LabelFrame(root, text="1. Cấu hình", padding=10)
        config_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Label(config_frame, text="Số chữ số làm tròn:").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.decimals_var = tk.IntVar(value=4)
        ttk.Entry(config_frame, textvariable=self.decimals_var, width=10).grid(row=0, column=1, sticky="w")
        
        input_frame = ttk.LabelFrame(root, text="2. Nhập ma trận vuông A", padding=10)
        input_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        toolbar = ttk.Frame(input_frame)
        toolbar.pack(fill="x", pady=(0, 5))
        ttk.Button(toolbar, text="Mở file .txt", command=self.load_file).pack(side="left")
        
        self.text_area = tk.Text(input_frame, wrap="none", font=("Consolas", 12))
        self.text_area.pack(fill="both", expand=True)
        
        self.text_area.insert("1.0", "3 5 15 12 11 11 4\n14 9 8 15 12 1 2\n2 15 13 10 12 5 14\n14 15 3 1 6 1 16\n4 6 7 12 10 9 1\n11 1 5 1 2 13 3\n5 9 15 15 3 15 10") 
        
        ttk.Button(root, text="Tạo file markdown", command=self.solve).pack(pady=15)

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
            
            save_path = filedialog.asksaveasfilename(defaultextension=".md", filetypes=[("Markdown", "*.md")])
            if save_path:
                md_content = run_bordering_method_final(A, self.decimals_var.get())
                with open(save_path, 'w', encoding='utf-8') as f:
                    f.write(md_content)
                messagebox.showinfo("Thành công", "Đã xuất file markdown!")
        except Exception as e:
            messagebox.showerror("Lỗi", str(e))

if __name__ == "__main__":
    root = tk.Tk()
    app = FinalBorderingApp(root)
    root.mainloop()
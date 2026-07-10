import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import io

# ==========================================
# THUẬT TOÁN: LẶP ĐƠN GIẢI HỆ X = BX + D
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

def run_lap_don(B, D, X0, epsilon, decimals):
    n = B.shape[0]
    if B.shape[0] != B.shape[1]:
        return "# LỖI: Ma trận B phải là ma trận vuông."
    
    # Đưa D và X0 về dạng ma trận cột nếu là mảng 1D
    D_mat = np.array(D)
    if D_mat.ndim == 1: D_mat = D_mat.reshape(-1, 1)
        
    X0_mat = np.array(X0)
    if X0_mat.ndim == 1: X0_mat = X0_mat.reshape(-1, 1)
        
    if D_mat.shape[0] != n or X0_mat.shape[0] != n:
        return f"# LỖI: Số hàng của D và X(0) phải bằng kích thước cấp của B ({n})."
        
    md = []
    md.append("## Thuật toán lặp đơn dạng $X = BX + D$\n\n")
    md.append("**Đầu vào (Input):**\n")
    md.append(f"- Ma trận vuông $B$ cấp $n={n}$:\n$$ B = {matrix_to_latex(B, decimals)} $$\n")
    md.append(f"- Ma trận $D$:\n$$ D = {matrix_to_latex(D_mat, decimals)} $$\n")
    md.append(f"- Sai số cho phép: $\\epsilon = {epsilon}$\n")
    md.append(f"- Xấp xỉ ban đầu $X^{(0)}$:\n$$ X^{(0)} = {matrix_to_latex(X0_mat, decimals)} $$\n\n")
    md.append("**Các bước thực hiện:**\n\n")
    
    # 1. Kiểm tra điều kiện hội tụ
    md.append("### 1. Kiểm tra điều kiện hội tụ.\n")
    md.append("Tính $q = \\|B\\|$. Ta xét ưu tiên theo thứ tự: Chuẩn hàng $\\to$ Chuẩn cột $\\to$ Chuẩn Frobenius.\n\n")
    
    q = None
    norm_type = None
    norm_symbol = ""
    
    norm_inf = np.linalg.norm(B, np.inf)
    md.append(f"- Kiểm tra chuẩn hàng: $\\|B\\|_\\infty = {format_num(norm_inf, decimals)}$.\n")
    if norm_inf < 1:
        q = norm_inf
        norm_type = np.inf
        norm_symbol = "\\infty"
        md.append(f"Vì $\\|B\\|_\\infty < 1$ nên ta chọn chuẩn hàng, suy ra $q = {format_num(q, decimals)}$.\n\n")
    else:
        norm_1 = np.linalg.norm(B, 1)
        md.append(f"Vì $\\|B\\|_\\infty \\ge 1$, ta xét tiếp chuẩn cột: $\\|B\\|_1 = {format_num(norm_1, decimals)}$.\n")
        if norm_1 < 1:
            q = norm_1
            norm_type = 1
            norm_symbol = "1"
            md.append(f"Vì $\\|B\\|_1 < 1$ nên ta chọn chuẩn cột, suy ra $q = {format_num(q, decimals)}$.\n\n")
        else:
            norm_fro = np.linalg.norm(B, 'fro')
            md.append(f"Vì $\\|B\\|_1 \\ge 1$, ta xét tiếp chuẩn Frobenius: $\\|B\\|_F = {format_num(norm_fro, decimals)}$.\n")
            if norm_fro < 1:
                q = norm_fro
                norm_type = 'fro'
                norm_symbol = "F"
                md.append(f"Vì $\\|B\\|_F < 1$ nên ta chọn chuẩn Frobenius, suy ra $q = {format_num(q, decimals)}$.\n\n")
            else:
                md.append("**DỪNG CHƯƠNG TRÌNH:** Không có chuẩn nào của ma trận $B$ nhỏ hơn 1 ($q \\ge 1$). Phương pháp lặp đơn không đảm bảo hội tụ.\n")
                return "".join(md)
                
    # 2. Khởi tạo sai số dừng
    eps_0 = ((1 - q) / q) * epsilon
    md.append("### 2. Khởi tạo sai số dừng.\n")
    md.append("Tính sai số dừng $\\epsilon_0$:\n")
    md.append(f"$$ \\epsilon_0 = \\frac{{1 - q}}{{q}} \\epsilon = \\frac{{1 - {format_num(q, decimals)}}}{{{format_num(q, decimals)}}} \\times {epsilon} = {format_num(eps_0, decimals)} $$\n")
    md.append("Đặt $k = 1$.\n\n")
    
    # 3. Lặp (Sử dụng Bảng)
    md.append("### 3. Lặp.\n")
    md.append(f"Sử dụng công thức lặp $X^{{(k)}} = BX^{{(k-1)}} + D$. Ta lập bảng các bước lặp như sau:\n\n")
    
    # Tạo tiêu đề bảng
    header = ["$k$"]
    for i in range(n):
        for j in range(X0_mat.shape[1]):
            if X0_mat.shape[1] == 1:
                header.append(f"$x_{{{i+1}}}$")
            else:
                header.append(f"$x_{{{i+1}{j+1}}}$")
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
    X_k_minus_1 = X0_mat.copy()
    
    while True:
        # Tính X^(k)
        X_k = B @ X_k_minus_1 + D_mat
        
        # Tính sai số
        diff = X_k - X_k_minus_1
        norm_diff = np.linalg.norm(diff, ord=norm_type)
        
        # Thêm dòng vào bảng
        row_k = [str(k)]
        for val in X_k.flatten():
            row_k.append(format_num(val, decimals))
        row_k.append(format_num(norm_diff, decimals))
        md.append("| " + " | ".join(row_k) + " |\n")
        
        # Kiểm tra điều kiện dừng
        if norm_diff <= eps_0:
            md.append("\n")
            md.append(f"Tại bước lặp $k = {k}$, ta thấy sai số $\\|X^{{({k})}} - X^{{({k-1})}}\\|_{norm_symbol} = {format_num(norm_diff, decimals)} \\le \\epsilon_0 = {format_num(eps_0, decimals)}$.\n\n")
            md.append(f"**Kết luận:** Thỏa mãn điều kiện dừng. Dừng chương trình và in ra nghiệm gần đúng $X = X^{{({k})}}$:\n")
            md.append(f"$$ X \\approx {matrix_to_latex(X_k, decimals)} $$\n")
            break
            
        X_k_minus_1 = X_k
        k += 1
        
        # Điều kiện chống lặp vô hạn
        if k > 500:
            md.append("\n**CẢNH BÁO:** Đã vượt quá 500 vòng lặp nhưng chưa hội tụ. Thuật toán tự động dừng.\n")
            break

    return "".join(md)

# ==========================================
# GIAO DIỆN TKINTER
# ==========================================

class LapDonApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Thuật toán Lặp Đơn: X = BX + D")
        self.root.geometry("750x700")
        
        # 1. Cấu hình
        config_frame = ttk.LabelFrame(root, text="1. Cấu hình thông số", padding=10)
        config_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Label(config_frame, text="Sai số (ε):").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.epsilon_var = tk.DoubleVar(value=0.01)
        ttk.Entry(config_frame, textvariable=self.epsilon_var, width=15).grid(row=0, column=1, sticky="w")
        
        ttk.Label(config_frame, text="Số chữ số làm tròn:").grid(row=0, column=2, padx=15, pady=5, sticky="e")
        self.decimals_var = tk.IntVar(value=4)
        ttk.Entry(config_frame, textvariable=self.decimals_var, width=10).grid(row=0, column=3, sticky="w")

        # Layout cột cho ma trận
        matrices_frame = ttk.Frame(root)
        matrices_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        # 2. Input B
        frame_b = ttk.LabelFrame(matrices_frame, text="2. Ma trận vuông B", padding=10)
        frame_b.pack(side="left", fill="both", expand=True, padx=(0, 5))
        self.text_b = tk.Text(frame_b, wrap="none", font=("Consolas", 11), width=20)
        self.text_b.pack(fill="both", expand=True)
        
        # 3. Input D
        frame_d = ttk.LabelFrame(matrices_frame, text="3. Ma trận/Vector D", padding=10)
        frame_d.pack(side="left", fill="both", expand=True, padx=5)
        self.text_d = tk.Text(frame_d, wrap="none", font=("Consolas", 11), width=15)
        self.text_d.pack(fill="both", expand=True)

        # 4. Input X0
        frame_x0 = ttk.LabelFrame(matrices_frame, text="4. Xấp xỉ ban đầu X(0)", padding=10)
        frame_x0.pack(side="left", fill="both", expand=True, padx=(5, 0))
        self.text_x0 = tk.Text(frame_x0, wrap="none", font=("Consolas", 11), width=15)
        self.text_x0.pack(fill="both", expand=True)

        # Dữ liệu mẫu
        sample_b = (
            "0.1 0.2 0.3\n"
            "0.2 0.1 -0.2\n"
            "0.1 0.1 0.1"
        )
        sample_d = (
            "1\n"
            "2\n"
            "3"
        )
        sample_x0 = (
            "1\n"
            "2\n"
            "3"
        )
        self.text_b.insert("1.0", sample_b) 
        self.text_d.insert("1.0", sample_d) 
        self.text_x0.insert("1.0", sample_x0) 
        
        btn_frame = ttk.Frame(root, padding=10)
        btn_frame.pack(fill="x")
        ttk.Button(btn_frame, text="GIẢI HỆ & XUẤT MARKDOWN", command=self.solve).pack(pady=10)

    def solve(self):
        try:
            raw_b = self.text_b.get("1.0", tk.END).strip()
            raw_d = self.text_d.get("1.0", tk.END).strip()
            raw_x0 = self.text_x0.get("1.0", tk.END).strip()
            
            if not raw_b or not raw_d or not raw_x0:
                messagebox.showwarning("Cảnh báo", "Vui lòng nhập đầy đủ ma trận B, D và X(0)!")
                return
                
            B = np.loadtxt(io.StringIO(raw_b))
            D = np.loadtxt(io.StringIO(raw_d))
            X0 = np.loadtxt(io.StringIO(raw_x0))
            
            if B.ndim == 1:
                B = B.reshape(1, -1)
                
            epsilon = self.epsilon_var.get()
            decimals = self.decimals_var.get()
            
            save_path = filedialog.asksaveasfilename(
                defaultextension=".md", 
                filetypes=[("Markdown", "*.md")], 
                initialfile="LoiGiai_LapDon.md"
            )
            
            if save_path:
                md_content = run_lap_don(B, D, X0, epsilon, decimals)
                with open(save_path, 'w', encoding='utf-8') as f:
                    f.write(md_content)
                messagebox.showinfo("Thành công", "Đã xuất file trình bày lời giải Markdown thành công!")
                
        except Exception as e:
            messagebox.showerror("Lỗi", f"Đã xảy ra lỗi:\n{str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = LapDonApp(root)
    root.mainloop()
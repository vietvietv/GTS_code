from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from dominant_eigen_solver import save_markdown_from_text, solve_text_to_markdown


DEFAULT_INPUT = """matrix:
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
"""


def open_external_editor(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(DEFAULT_INPUT, encoding="utf-8")

    if os.name == "nt":
        subprocess.Popen(["notepad.exe", str(path)])
    elif sys.platform == "darwin":
        subprocess.Popen(["open", str(path)])
    else:
        subprocess.Popen(["xdg-open", str(path)])


class DominantEigenApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Dominant Eigenvalue Solver")
        self.geometry("1040x740")
        self.minsize(860, 580)

        self.input_path = Path("dominant_input.txt")
        self.output_path = Path("dominant_solution.md")

        self._build_widgets()
        self.input_text.insert("1.0", DEFAULT_INPUT)

    def _build_widgets(self) -> None:
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        top = ttk.Frame(self, padding=10)
        top.grid(row=0, column=0, sticky="ew")
        top.columnconfigure(1, weight=1)

        ttk.Label(top, text="File output (.md):").grid(row=0, column=0, sticky="w")
        self.output_entry = ttk.Entry(top)
        self.output_entry.insert(0, str(self.output_path))
        self.output_entry.grid(row=0, column=1, sticky="ew", padx=8)
        ttk.Button(top, text="Chọn file...", command=self.choose_output).grid(row=0, column=2)

        main = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        main.grid(row=1, column=0, sticky="nsew", padx=10, pady=(0, 10))

        left = ttk.Frame(main)
        right = ttk.Frame(main)
        main.add(left, weight=1)
        main.add(right, weight=1)

        left.rowconfigure(1, weight=1)
        left.columnconfigure(0, weight=1)
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)

        ttk.Label(left, text="Nhập ma trận và tham số").grid(row=0, column=0, sticky="w")
        self.input_text = tk.Text(left, wrap="none", undo=True, font=("Consolas", 12))
        self.input_text.grid(row=1, column=0, sticky="nsew", pady=(6, 0))
        left_scroll = ttk.Scrollbar(left, orient="vertical", command=self.input_text.yview)
        left_scroll.grid(row=1, column=1, sticky="ns", pady=(6, 0))
        self.input_text.configure(yscrollcommand=left_scroll.set)

        ttk.Label(right, text="Preview Markdown").grid(row=0, column=0, sticky="w")
        self.preview_text = tk.Text(right, wrap="word", undo=False, font=("Consolas", 10))
        self.preview_text.grid(row=1, column=0, sticky="nsew", pady=(6, 0))
        right_scroll = ttk.Scrollbar(right, orient="vertical", command=self.preview_text.yview)
        right_scroll.grid(row=1, column=1, sticky="ns", pady=(6, 0))
        self.preview_text.configure(yscrollcommand=right_scroll.set)

        bottom = ttk.Frame(self, padding=(10, 0, 10, 10))
        bottom.grid(row=2, column=0, sticky="ew")

        ttk.Button(bottom, text="Mở input bằng Notepad", command=self.open_notepad_input).pack(side=tk.LEFT)
        ttk.Button(bottom, text="Đọc từ file input", command=self.load_input_file).pack(side=tk.LEFT, padx=6)
        ttk.Button(bottom, text="Lưu input", command=self.save_input_file).pack(side=tk.LEFT)
        ttk.Button(bottom, text="Xem trước", command=self.preview).pack(side=tk.RIGHT, padx=6)
        ttk.Button(bottom, text="Giải và xuất .md", command=self.solve_and_export).pack(side=tk.RIGHT)

    def choose_output(self) -> None:
        path = filedialog.asksaveasfilename(
            defaultextension=".md",
            filetypes=[("Markdown", "*.md"), ("All files", "*.*")],
            initialfile=self.output_entry.get() or "dominant_solution.md",
        )
        if path:
            self.output_entry.delete(0, tk.END)
            self.output_entry.insert(0, path)

    def open_notepad_input(self) -> None:
        self.save_input_file(show_message=False)
        try:
            open_external_editor(self.input_path)
        except Exception as exc:
            messagebox.showerror("Lỗi", f"Không mở được editor ngoài:\n{exc}")

    def load_input_file(self) -> None:
        path = filedialog.askopenfilename(
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            initialfile=str(self.input_path),
        )
        if not path:
            return
        self.input_path = Path(path)
        self.input_text.delete("1.0", tk.END)
        self.input_text.insert("1.0", self.input_path.read_text(encoding="utf-8"))

    def save_input_file(self, show_message: bool = True) -> None:
        self.input_path.write_text(self.input_text.get("1.0", tk.END).strip() + "\n", encoding="utf-8")
        if show_message:
            messagebox.showinfo("Đã lưu", f"Đã lưu input vào {self.input_path}")

    def preview(self) -> None:
        try:
            markdown = solve_text_to_markdown(self.input_text.get("1.0", tk.END))
        except Exception as exc:
            messagebox.showerror("Lỗi", str(exc))
            return
        self.preview_text.delete("1.0", tk.END)
        self.preview_text.insert("1.0", markdown)

    def solve_and_export(self) -> None:
        output = Path(self.output_entry.get().strip() or "dominant_solution.md")
        try:
            save_markdown_from_text(self.input_text.get("1.0", tk.END), str(output))
        except Exception as exc:
            messagebox.showerror("Lỗi", str(exc))
            return
        messagebox.showinfo("Xong", f"Đã xuất lời giải ra:\n{output}")


def run_cli() -> int:
    parser = argparse.ArgumentParser(description="Tìm giá trị riêng trội và xuất Markdown.")
    parser.add_argument("--input", "-i", help="File text chứa ma trận và tham số.")
    parser.add_argument("--output", "-o", default="dominant_solution.md", help="File Markdown output.")
    parser.add_argument("--gui", action="store_true", help="Mở giao diện tkinter.")
    args = parser.parse_args()

    if args.gui or not args.input:
        app = DominantEigenApp()
        app.mainloop()
        return 0

    matrix_text = Path(args.input).read_text(encoding="utf-8")
    save_markdown_from_text(matrix_text, args.output)
    print(f"Đã xuất lời giải ra {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run_cli())

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from danielevsky_solver import save_markdown_from_text, solve_text_to_markdown


DEFAULT_INPUT = """4 -1 -1
0 1 0
2 2 1
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


class DanielevskyApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Danielevsky Solver")
        self.geometry("980x720")
        self.minsize(820, 560)

        self.input_path = Path("matrix_input.txt")
        self.output_path = Path("danielevsky_solution.md")

        self._build_widgets()
        self.matrix_text.insert("1.0", DEFAULT_INPUT)

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

        ttk.Button(top, text="Chon file...", command=self.choose_output).grid(row=0, column=2)

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

        ttk.Label(left, text="Nhap ma tran").grid(row=0, column=0, sticky="w")
        self.matrix_text = tk.Text(left, wrap="none", undo=True, font=("Consolas", 12))
        self.matrix_text.grid(row=1, column=0, sticky="nsew", pady=(6, 0))

        y_scroll_left = ttk.Scrollbar(left, orient="vertical", command=self.matrix_text.yview)
        y_scroll_left.grid(row=1, column=1, sticky="ns", pady=(6, 0))
        self.matrix_text.configure(yscrollcommand=y_scroll_left.set)

        ttk.Label(right, text="Preview Markdown").grid(row=0, column=0, sticky="w")
        self.preview_text = tk.Text(right, wrap="word", undo=False, font=("Consolas", 10))
        self.preview_text.grid(row=1, column=0, sticky="nsew", pady=(6, 0))

        y_scroll_right = ttk.Scrollbar(right, orient="vertical", command=self.preview_text.yview)
        y_scroll_right.grid(row=1, column=1, sticky="ns", pady=(6, 0))
        self.preview_text.configure(yscrollcommand=y_scroll_right.set)

        bottom = ttk.Frame(self, padding=(10, 0, 10, 10))
        bottom.grid(row=2, column=0, sticky="ew")

        ttk.Button(bottom, text="Mo input bang Notepad", command=self.open_notepad_input).pack(side=tk.LEFT)
        ttk.Button(bottom, text="Doc tu file input", command=self.load_input_file).pack(side=tk.LEFT, padx=6)
        ttk.Button(bottom, text="Luu input", command=self.save_input_file).pack(side=tk.LEFT)
        ttk.Button(bottom, text="Xem truoc", command=self.preview).pack(side=tk.RIGHT, padx=6)
        ttk.Button(bottom, text="Giai va xuat .md", command=self.solve_and_export).pack(side=tk.RIGHT)

    def choose_output(self) -> None:
        path = filedialog.asksaveasfilename(
            defaultextension=".md",
            filetypes=[("Markdown", "*.md"), ("All files", "*.*")],
            initialfile=self.output_entry.get() or "danielevsky_solution.md",
        )
        if path:
            self.output_entry.delete(0, tk.END)
            self.output_entry.insert(0, path)

    def open_notepad_input(self) -> None:
        self.save_input_file(show_message=False)
        try:
            open_external_editor(self.input_path)
        except Exception as exc:
            messagebox.showerror("Loi", f"Khong mo duoc editor ngoai:\n{exc}")

    def load_input_file(self) -> None:
        path = filedialog.askopenfilename(
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            initialfile=str(self.input_path),
        )
        if not path:
            return
        self.input_path = Path(path)
        self.matrix_text.delete("1.0", tk.END)
        self.matrix_text.insert("1.0", self.input_path.read_text(encoding="utf-8"))

    def save_input_file(self, show_message: bool = True) -> None:
        self.input_path.write_text(self.matrix_text.get("1.0", tk.END).strip() + "\n", encoding="utf-8")
        if show_message:
            messagebox.showinfo("Da luu", f"Da luu input vao {self.input_path}")

    def preview(self) -> None:
        try:
            markdown = solve_text_to_markdown(self.matrix_text.get("1.0", tk.END))
        except Exception as exc:
            messagebox.showerror("Loi", str(exc))
            return

        self.preview_text.delete("1.0", tk.END)
        self.preview_text.insert("1.0", markdown)

    def solve_and_export(self) -> None:
        output = Path(self.output_entry.get().strip() or "danielevsky_solution.md")
        try:
            save_markdown_from_text(self.matrix_text.get("1.0", tk.END), str(output))
        except Exception as exc:
            messagebox.showerror("Loi", str(exc))
            return
        messagebox.showinfo("Xong", f"Da xuat loi giai ra:\n{output}")


def run_cli() -> int:
    parser = argparse.ArgumentParser(description="Giai bai toan Danielevsky va xuat Markdown.")
    parser.add_argument("--input", "-i", help="File text chua ma tran.")
    parser.add_argument("--output", "-o", default="danielevsky_solution.md", help="File Markdown output.")
    parser.add_argument("--gui", action="store_true", help="Mo giao dien tkinter.")
    args = parser.parse_args()

    if args.gui or not args.input:
        app = DanielevskyApp()
        app.mainloop()
        return 0

    matrix_text = Path(args.input).read_text(encoding="utf-8")
    save_markdown_from_text(matrix_text, args.output)
    print(f"Da xuat loi giai ra {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run_cli())

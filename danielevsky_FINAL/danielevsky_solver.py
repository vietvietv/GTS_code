from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import List, Sequence


EPS = 1e-10


Number = float
Matrix = List[List[Number]]
Vector = List[complex]


@dataclass
class Step:
    title: str
    description: List[str]
    matrices: List[tuple[str, Matrix]]
    scalars: List[tuple[str, float]]


@dataclass
class FrobeniusBlock:
    start: int
    size: int
    coefficients: List[float]
    matrix: Matrix


@dataclass
class Solution:
    original: Matrix
    frobenius: Matrix
    transform: Matrix
    inverse_transform: Matrix
    coefficients: List[float]
    blocks: List[FrobeniusBlock]
    eigenvalues: List[complex]
    frobenius_vectors: List[Vector]
    original_vectors: List[Vector]
    steps: List[Step]
    warnings: List[str]


def parse_matrix(text: str) -> Matrix:
    rows: Matrix = []
    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        line = (
            line.replace("[", " ")
            .replace("]", " ")
            .replace("(", " ")
            .replace(")", " ")
            .replace(";", " ")
            .replace(",", " ")
        )
        values = [float(token) for token in line.split()]
        rows.append(values)

    if not rows:
        raise ValueError("Chua nhap ma tran.")

    n = len(rows)
    if any(len(row) != n for row in rows):
        raise ValueError("Ma tran phai vuong: so hang va so cot phai bang nhau.")

    return rows


def identity(n: int) -> Matrix:
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def copy_matrix(a: Matrix) -> Matrix:
    return [row[:] for row in a]


def submatrix(a: Matrix, start: int, size: int) -> Matrix:
    return [row[start : start + size] for row in a[start : start + size]]


def mat_mul(a: Matrix, b: Matrix) -> Matrix:
    rows = len(a)
    cols = len(b[0])
    inner = len(b)
    return [[sum(a[i][k] * b[k][j] for k in range(inner)) for j in range(cols)] for i in range(rows)]


def mat_vec_mul(a: Matrix, v: Sequence[complex]) -> Vector:
    return [sum(a[i][j] * v[j] for j in range(len(v))) for i in range(len(a))]


def swap_similarity(a: Matrix, i: int, j: int) -> Matrix:
    p = identity(len(a))
    p[i], p[j] = p[j], p[i]
    return p


def inverse(a: Matrix) -> Matrix:
    n = len(a)
    aug = [a[i][:] + identity(n)[i] for i in range(n)]

    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot][col]) < EPS:
            raise ValueError("Ma tran bien doi bi suy bien, khong the nghich dao.")
        aug[col], aug[pivot] = aug[pivot], aug[col]

        div = aug[col][col]
        aug[col] = [value / div for value in aug[col]]

        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            if abs(factor) > EPS:
                aug[row] = [aug[row][c] - factor * aug[col][c] for c in range(2 * n)]

    return [row[n:] for row in aug]


def clean_zero(value: complex, eps: float = 1e-9) -> complex:
    real = 0.0 if abs(value.real) < eps else value.real
    imag = 0.0 if abs(value.imag) < eps else value.imag
    return complex(real, imag)


def poly_multiply(a: Sequence[float], b: Sequence[float]) -> List[float]:
    result = [0.0 for _ in range(len(a) + len(b) - 1)]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] += ai * bj
    return result


def eval_poly(coefficients: Sequence[float], x: complex) -> complex:
    value = 0j
    for coefficient in coefficients:
        value = value * x + coefficient
    return value


def polynomial_roots(coefficients: Sequence[float]) -> List[complex]:
    if len(coefficients) < 2:
        return []

    lead = coefficients[0]
    if abs(lead) < EPS:
        raise ValueError("He so bac cao nhat cua da thuc bang 0.")
    coefficients = [c / lead for c in coefficients]
    degree = len(coefficients) - 1

    if degree == 1:
        return [complex(-coefficients[1])]

    if degree == 2:
        _, b, c = coefficients
        delta = complex(b * b - 4 * c)
        return [(-b + cmath.sqrt(delta)) / 2, (-b - cmath.sqrt(delta)) / 2]

    radius = 1.0 + max(abs(c) for c in coefficients[1:])
    roots = [
        radius * cmath.exp(2j * math.pi * (i + 0.25) / degree)
        for i in range(degree)
    ]

    for _ in range(3000):
        max_change = 0.0
        next_roots = roots[:]
        for i, root in enumerate(roots):
            denominator = 1 + 0j
            for j, other in enumerate(roots):
                if i != j:
                    denominator *= root - other
            if abs(denominator) < EPS:
                denominator = complex(EPS, EPS)
            delta = eval_poly(coefficients, root) / denominator
            next_roots[i] = root - delta
            max_change = max(max_change, abs(delta))
        roots = next_roots
        if max_change < 1e-12:
            break

    roots = [clean_zero(root) for root in roots]
    roots.sort(key=lambda z: (round(z.real, 10), round(z.imag, 10)))
    return roots


def nullspace_vector(a: List[List[complex]]) -> Vector:
    rows = len(a)
    cols = len(a[0])
    rref = [[complex(value) for value in row] for row in a]
    pivot_cols: List[int] = []
    row = 0

    for col in range(cols):
        pivot = None
        best = EPS
        for candidate in range(row, rows):
            size = abs(rref[candidate][col])
            if size > best:
                best = size
                pivot = candidate
        if pivot is None:
            continue

        rref[row], rref[pivot] = rref[pivot], rref[row]
        div = rref[row][col]
        rref[row] = [value / div for value in rref[row]]

        for other in range(rows):
            if other == row:
                continue
            factor = rref[other][col]
            if abs(factor) > EPS:
                rref[other] = [
                    rref[other][j] - factor * rref[row][j] for j in range(cols)
                ]

        pivot_cols.append(col)
        row += 1
        if row == rows:
            break

    free_cols = [col for col in range(cols) if col not in pivot_cols]
    if not free_cols:
        free_cols = [cols - 1]

    free_col = free_cols[-1]
    vector = [0j for _ in range(cols)]
    vector[free_col] = 1 + 0j

    for pivot_row in range(len(pivot_cols) - 1, -1, -1):
        pivot_col = pivot_cols[pivot_row]
        vector[pivot_col] = -sum(
            rref[pivot_row][j] * vector[j] for j in range(pivot_col + 1, cols)
        )

    max_abs = max(abs(value) for value in vector)
    if max_abs > EPS:
        vector = [value / max_abs for value in vector]
    return [clean_zero(value) for value in vector]


def eigenvector_for_matrix(a: Matrix, eigenvalue: complex) -> Vector:
    n = len(a)
    matrix = [
        [
            complex(a[i][j]) - (eigenvalue if i == j else 0j)
            for j in range(n)
        ]
        for i in range(n)
    ]
    return nullspace_vector(matrix)


def solve_danielevsky(a: Matrix) -> Solution:
    n = len(a)
    current = copy_matrix(a)
    transform = identity(n)
    steps: List[Step] = []
    warnings: List[str] = []
    blocks: List[FrobeniusBlock] = []

    def add_block(start: int, size: int) -> None:
        block_matrix = submatrix(current, start, size)
        coefficients = [-block_matrix[0][j] for j in range(size)]
        blocks.append(
            FrobeniusBlock(
                start=start,
                size=size,
                coefficients=coefficients,
                matrix=block_matrix,
            )
        )

    def solve_prefix(size: int) -> None:
        nonlocal current, transform
        if size <= 0:
            return
        if size == 1:
            add_block(0, 1)
            return

        k = size - 1
        while k > 0:
            pivot = current[k][k - 1]

            if abs(pivot) < EPS:
                swap_col = None
                for s in range(k - 1):
                    if abs(current[k][s]) >= EPS:
                        swap_col = s
                        break

                if swap_col is None:
                    bottom_size = size - k
                    add_block(k, bottom_size)
                    before = copy_matrix(current)
                    steps.append(
                        Step(
                            title=f"Buoc tach khoi tai hang {k + 1}",
                            description=[
                                f"Xét hàng {k + 1}, ta có "
                                f"\\(a_{{{k + 1},1}}=\\cdots=a_{{{k + 1},{k}}}=0\\).",
                                "Do đó thuộc trường hợp 3 của phương pháp Danielevsky.",
                                "Ma trận khi đó có dạng khối tam giác trên "
                                "\\(\\begin{pmatrix}A_1&B\\\\0&F_1\\end{pmatrix}\\).",
                                "Suy ra đa thức đặc trưng bằng tích đa thức đặc trưng của các khối đường chéo.",
                                f"Ta tách khối Frobenius cấp {bottom_size} ở góc dưới bên phải "
                                f"và tiếp tục xử lý ma trận con cấp {k} ở phía trên.",
                            ],
                            matrices=[
                                ("A_{\\mathrm{hiện\\ tại}}", before),
                                ("F_1", submatrix(before, k, bottom_size)),
                            ],
                            scalars=[],
                        )
                    )
                    solve_prefix(k)
                    return

                permutation = swap_similarity(current, swap_col, k - 1)
                before = copy_matrix(current)
                current = mat_mul(mat_mul(permutation, current), permutation)
                transform = mat_mul(permutation, transform)
                steps.append(
                    Step(
                        title=f"Buoc hoan vi tai hang {k + 1}",
                        description=[
                            f"Vì \\(a_{{{k + 1},{k}}}=0\\) nhưng "
                            f"\\(a_{{{k + 1},{swap_col + 1}}}\\ne0\\), "
                            f"ta hoán vị đồng thời hàng/cột {swap_col + 1} và {k}.",
                            "Phép biến đổi này là biến đổi đồng dạng nên không làm thay đổi giá trị riêng.",
                        ],
                        matrices=[
                            ("A_{\\mathrm{trước}}", before),
                            ("C", permutation),
                            ("A_{\\mathrm{sau}}", current),
                        ],
                        scalars=[],
                    )
                )
                pivot = current[k][k - 1]

            if abs(pivot) < EPS:
                raise ValueError(f"Không tìm được phần tử trụ tại hàng {k + 1}.")

            m = identity(n)
            for col in range(n):
                m[k - 1][col] = 0.0
            for col in range(size):
                m[k - 1][col] = current[k][col]
            m_inv = inverse(m)
            before = copy_matrix(current)
            current = mat_mul(mat_mul(m, current), m_inv)
            transform = mat_mul(m, transform)

            steps.append(
                Step(
                    title=f"Buoc dua hang {k + 1} ve dang Frobenius",
                    description=[
                        f"Xét hàng {k + 1}, ta có "
                        f"\\(a_{{{k + 1},{k}}}={pivot:.7f}\\ne0\\).",
                        "Do đó thuộc trường hợp 1 của phương pháp Danielevsky.",
                        f"Chọn \\(M\\) sao cho hàng {k} của \\(M\\) bằng phần tương ứng của hàng {k + 1}.",
                        "Khi đó \\(A_{\\mathrm{mới}}=MA_{\\mathrm{cũ}}M^{-1}\\).",
                    ],
                    matrices=[
                        ("A_{\\mathrm{cũ}}", before),
                        ("M", m),
                        ("M^{-1}", m_inv),
                        ("A_{\\mathrm{mới}}", current),
                    ],
                    scalars=[(f"a_{{{k + 1},{k}}}", pivot)],
                )
            )
            k -= 1

        add_block(0, size)

    solve_prefix(n)
    blocks.sort(key=lambda block: block.start)

    frobenius = current
    polynomial_coefficients = [1.0]
    for block in blocks:
        polynomial_coefficients = poly_multiply(polynomial_coefficients, [1.0] + block.coefficients)
    coefficients = polynomial_coefficients[1:]
    roots = polynomial_roots(polynomial_coefficients)
    transform_inverse = inverse(transform)

    frobenius_vectors: List[Vector] = []
    original_vectors: List[Vector] = []
    for value in roots:
        if len(blocks) == 1 and blocks[0].start == 0:
            u = [value ** power for power in range(n - 1, -1, -1)]
        else:
            u = eigenvector_for_matrix(frobenius, value)
        x = mat_vec_mul(transform_inverse, u)
        frobenius_vectors.append([clean_zero(v) for v in u])
        original_vectors.append([clean_zero(v) for v in x])

    return Solution(
        original=copy_matrix(a),
        frobenius=frobenius,
        transform=transform,
        inverse_transform=transform_inverse,
        coefficients=coefficients,
        blocks=blocks,
        eigenvalues=roots,
        frobenius_vectors=frobenius_vectors,
        original_vectors=original_vectors,
        steps=steps,
        warnings=warnings,
    )


def format_real(value: float, digits: int) -> str:
    if abs(value) < 0.5 * 10 ** (-digits):
        value = 0.0
    text = f"{value:.{digits}f}"
    return "0" if text == f"-{0:.{digits}f}" else text


def format_complex(value: complex, digits: int = 7) -> str:
    value = clean_zero(value)
    if abs(value.imag) < EPS:
        return format_real(value.real, digits)
    if abs(value.real) < EPS:
        return f"{format_real(value.imag, digits)}i"
    sign = "+" if value.imag >= 0 else "-"
    return f"{format_real(value.real, digits)} {sign} {format_real(abs(value.imag), digits)}i"


def matrix_latex(a: Matrix, digits: int = 4) -> str:
    rows = []
    for row in a:
        rows.append(" & ".join(format_real(value, digits) for value in row))
    return "\\begin{pmatrix}\n" + "\\\\\n".join(rows) + "\n\\end{pmatrix}"


def vector_latex(v: Sequence[complex], digits: int = 7) -> str:
    rows = "\\\\\n".join(format_complex(value, digits) for value in v)
    return "\\begin{pmatrix}\n" + rows + "\n\\end{pmatrix}"


def polynomial_latex(coefficients: Sequence[float]) -> str:
    n = len(coefficients)
    pieces = ["\\lambda^{" + str(n) + "}"]
    for index, coefficient in enumerate(coefficients):
        power = n - index - 1
        if abs(coefficient) < EPS:
            continue
        sign = "+" if coefficient > 0 else "-"
        abs_coeff = abs(coefficient)
        coeff_text = "" if abs(abs_coeff - 1.0) < EPS and power > 0 else format_real(abs_coeff, 7)
        if power == 0:
            term = coeff_text or "1"
        elif power == 1:
            term = f"{coeff_text}\\lambda" if coeff_text else "\\lambda"
        else:
            term = f"{coeff_text}\\lambda^{{{power}}}" if coeff_text else f"\\lambda^{{{power}}}"
        pieces.append(f"{sign} {term}")
    return " ".join(pieces)


def render_markdown(solution: Solution) -> str:
    lines: List[str] = []
    n = len(solution.original)

    lines.append("# Lời giải bằng phương pháp Danielevsky")
    lines.append("")
    lines.append("> Các phép tính bên trong sử dụng giá trị gốc, không dùng giá trị đã làm tròn để tính tiếp.")
    lines.append("> Khi trình bày, ma trận được làm tròn 4 chữ số thập phân; các số liệu khác làm tròn 7 chữ số thập phân.")
    lines.append("")
    lines.append("## Bài toán")
    lines.append("")
    lines.append("Cho ma trận")
    lines.append("")
    lines.append("$$")
    lines.append("A=" + matrix_latex(solution.original, 4))
    lines.append("$$")
    lines.append("")
    lines.append("Dùng phương pháp Danielevsky tìm đa thức đặc trưng, giá trị riêng và véc tơ riêng của ma trận \\(A\\).")
    lines.append("")

    if solution.warnings:
        lines.append("## Lưu ý")
        lines.append("")
        for warning in solution.warnings:
            lines.append(f"- {warning}")
        lines.append("")

    lines.append("## Quá trình biến đổi")
    lines.append("")
    lines.append("Đặt")
    lines.append("")
    lines.append("$$")
    lines.append("A^{(1)}=A.")
    lines.append("$$")
    lines.append("")

    for index, step in enumerate(solution.steps, start=1):
        title = (
            step.title.replace("Buoc", "Bước")
            .replace("dua", "đưa")
            .replace("hang", "hàng")
            .replace("ve dang", "về dạng")
            .replace("hoan vi", "hoán vị")
            .replace("tach khoi", "tách khối")
            .replace("tai", "tại")
        )
        lines.append(f"### {index}. {title}")
        lines.append("")
        for item in step.description:
            lines.append(item)
            lines.append("")
        for name, matrix in step.matrices:
            lines.append("$$")
            lines.append(f"{name}=" + matrix_latex(matrix, 4))
            lines.append("$$")
            lines.append("")

    if len(solution.blocks) == 1:
        lines.append("Sau các phép biến đổi đồng dạng, thu được ma trận Frobenius")
    else:
        lines.append("Sau các phép biến đổi đồng dạng, thu được ma trận Frobenius dạng khối")
    lines.append("")
    lines.append("$$")
    lines.append("F=" + matrix_latex(solution.frobenius, 4))
    lines.append("$$")
    lines.append("")

    lines.append("Ma trận biến đổi tổng hợp là")
    lines.append("")
    lines.append("$$")
    lines.append("P=" + matrix_latex(solution.transform, 4))
    lines.append("$$")
    lines.append("")
    lines.append("và")
    lines.append("")
    lines.append("$$")
    lines.append("P^{-1}=" + matrix_latex(solution.inverse_transform, 4))
    lines.append("$$")
    lines.append("")
    lines.append("Khi đó")
    lines.append("")
    lines.append("$$")
    lines.append("F=PAP^{-1}.")
    lines.append("$$")
    lines.append("")

    lines.append("## Đa thức đặc trưng")
    lines.append("")
    if len(solution.blocks) == 1:
        lines.append("Với ma trận Frobenius dạng")
        lines.append("")
        lines.append("$$")
        lines.append("F=\\begin{pmatrix}-p_1&-p_2&\\cdots&-p_n\\\\1&0&\\cdots&0\\\\0&1&\\cdots&0\\\\\\vdots&\\vdots&\\ddots&\\vdots\\end{pmatrix},")
        lines.append("$$")
        lines.append("")
        lines.append("ta có")
        lines.append("")
        lines.append("$$")
        lines.append("p(\\lambda)=\\lambda^n+p_1\\lambda^{n-1}+\\cdots+p_n.")
        lines.append("$$")
        lines.append("")
        for i, coefficient in enumerate(solution.coefficients, start=1):
            lines.append(f"\\(p_{i}={format_real(coefficient, 7)}\\).")
        lines.append("")
    else:
        lines.append("Vì ma trận cuối có dạng khối tam giác trên, nên đa thức đặc trưng bằng tích đa thức đặc trưng của các khối đường chéo.")
        lines.append("")
        for index, block in enumerate(solution.blocks, start=1):
            block_poly = polynomial_latex(block.coefficients)
            lines.append(f"Khối \\(F_{index}\\) cấp {block.size}:")
            lines.append("")
            lines.append("$$")
            lines.append(f"F_{index}=" + matrix_latex(block.matrix, 4))
            lines.append("$$")
            lines.append("")
            lines.append("$$")
            lines.append(f"p_{index}(\\lambda)={block_poly}.")
            lines.append("$$")
            lines.append("")
    polynomial = polynomial_latex(solution.coefficients)
    lines.append("Do đó")
    lines.append("")
    lines.append("$$")
    lines.append(f"p(\\lambda)={polynomial}.")
    lines.append("$$")
    lines.append("")

    lines.append("## Giá trị riêng và véc tơ riêng")
    lines.append("")
    lines.append("Giải phương trình")
    lines.append("")
    lines.append("$$")
    lines.append(f"{polynomial}=0")
    lines.append("$$")
    lines.append("")
    lines.append("thu được các giá trị riêng:")
    lines.append("")
    for i, value in enumerate(solution.eigenvalues, start=1):
        lines.append(f"- \\(\\lambda_{i}={format_complex(value, 7)}\\)")
    lines.append("")
    if len(solution.blocks) == 1:
        lines.append("Với mỗi giá trị riêng \\(\\lambda\\), véc tơ riêng của khối Frobenius có dạng")
        lines.append("")
        lines.append("$$")
        lines.append("u_\\lambda=\\begin{pmatrix}\\lambda^{n-1}\\\\\\lambda^{n-2}\\\\\\vdots\\\\1\\end{pmatrix}.")
        lines.append("$$")
        lines.append("")
    else:
        lines.append("Với mỗi giá trị riêng \\(\\lambda\\), ta tìm véc tơ riêng \\(u_\\lambda\\) của ma trận Frobenius khối bằng cách giải")
        lines.append("")
        lines.append("$$")
        lines.append("(F-\\lambda I)u_\\lambda=0.")
        lines.append("$$")
        lines.append("")
    lines.append("Vì \\(F=PAP^{-1}\\), véc tơ riêng của ma trận ban đầu là")
    lines.append("")
    lines.append("$$")
    lines.append("x_\\lambda=P^{-1}u_\\lambda.")
    lines.append("$$")
    lines.append("")

    for i, (value, u, x) in enumerate(
        zip(solution.eigenvalues, solution.frobenius_vectors, solution.original_vectors),
        start=1,
    ):
        lines.append(f"Với \\(\\lambda_{i}={format_complex(value, 7)}\\):")
        lines.append("")
        lines.append("$$")
        lines.append(f"u_{i}=" + vector_latex(u, 7))
        lines.append("$$")
        lines.append("")
        lines.append("$$")
        lines.append(f"x_{i}=P^{{-1}}u_{i}=" + vector_latex(x, 7))
        lines.append("$$")
        lines.append("")

    lines.append("## Kết luận")
    lines.append("")
    lines.append("Đa thức đặc trưng của ma trận là")
    lines.append("")
    lines.append("$$")
    lines.append(f"p(\\lambda)={polynomial}.")
    lines.append("$$")
    lines.append("")
    lines.append("Các cặp giá trị riêng và véc tơ riêng tương ứng là:")
    lines.append("")
    for i, (value, x) in enumerate(zip(solution.eigenvalues, solution.original_vectors), start=1):
        lines.append("$$")
        lines.append(f"\\lambda_{i}={format_complex(value, 7)},\\qquad x_{i}=" + vector_latex(x, 7) + ".")
        lines.append("$$")
        lines.append("")

    return "\n".join(lines)


def solve_text_to_markdown(text: str) -> str:
    matrix = parse_matrix(text)
    solution = solve_danielevsky(matrix)
    return render_markdown(solution)


def save_markdown_from_text(matrix_text: str, output_path: str) -> None:
    markdown = solve_text_to_markdown(matrix_text)
    with open(output_path, "w", encoding="utf-8") as file:
        file.write(markdown)

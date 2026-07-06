from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Sequence


EPS = 1e-12

Matrix = List[List[float]]
Vector = List[float]


@dataclass
class SVDOptions:
    tol: float = 1e-7
    rank_tol: float = 1e-9
    max_iter: int = 100
    components: int | None = None
    mode: str = "exam"
    full_svd: bool = False
    pinv: bool = True
    condition_number: bool = True


@dataclass
class EigenStep:
    index: int
    eigenvalue: float
    eigenvector: Vector
    iterations: int
    residual: float


@dataclass
class SVDSolution:
    matrix: Matrix
    options: SVDOptions
    branch: str
    gram_matrix: Matrix
    eigen_steps: List[EigenStep]
    singular_values: Vector
    left_vectors: List[Vector]
    right_vectors: List[Vector]
    rank: int
    reduced_sigma: Matrix
    pseudo_inverse: Matrix | None
    condition_number: float | None
    reconstruction: Matrix
    reconstruction_error: float


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y", "co", "có"}


def parse_input(text: str) -> tuple[Matrix, SVDOptions]:
    rows: Matrix = []
    options = SVDOptions()
    section = "matrix"

    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        lower = line.lower().replace(" ", "")

        if lower in {"matrix:", "a:"}:
            section = "matrix"
            continue
        if lower.startswith("tol:") or lower.startswith("eps:"):
            options.tol = float(line.split(":", 1)[1].strip())
            continue
        if lower.startswith("rank_tol:") or lower.startswith("ranktol:"):
            options.rank_tol = float(line.split(":", 1)[1].strip())
            continue
        if lower.startswith("max_iter:") or lower.startswith("maxiter:"):
            options.max_iter = int(float(line.split(":", 1)[1].strip()))
            continue
        if lower.startswith("components:") or lower.startswith("k:"):
            options.components = int(float(line.split(":", 1)[1].strip()))
            continue
        if lower.startswith("mode:"):
            options.mode = line.split(":", 1)[1].strip().lower()
            continue
        if lower.startswith("full_svd:") or lower.startswith("full:"):
            options.full_svd = parse_bool(line.split(":", 1)[1])
            continue
        if lower.startswith("pinv:") or lower.startswith("pseudo_inverse:"):
            options.pinv = parse_bool(line.split(":", 1)[1])
            continue
        if lower.startswith("condition_number:") or lower.startswith("cond:"):
            options.condition_number = parse_bool(line.split(":", 1)[1])
            continue

        if section == "matrix":
            values = [
                float(token)
                for token in line.replace("[", " ")
                .replace("]", " ")
                .replace("(", " ")
                .replace(")", " ")
                .replace(",", " ")
                .replace(";", " ")
                .split()
            ]
            rows.append(values)

    if not rows:
        raise ValueError("Chưa nhập ma trận.")
    width = len(rows[0])
    if width == 0 or any(len(row) != width for row in rows):
        raise ValueError("Các hàng của ma trận phải có cùng số phần tử.")
    if options.mode not in {"exam", "full"}:
        raise ValueError("mode chỉ nhận 'exam' hoặc 'full'.")
    return rows, options


def transpose(a: Matrix) -> Matrix:
    return [[a[i][j] for i in range(len(a))] for j in range(len(a[0]))]


def mat_mul(a: Matrix, b: Matrix) -> Matrix:
    rows, inner, cols = len(a), len(b), len(b[0])
    return [[sum(a[i][k] * b[k][j] for k in range(inner)) for j in range(cols)] for i in range(rows)]


def mat_sub(a: Matrix, b: Matrix) -> Matrix:
    return [[a[i][j] - b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def mat_vec_mul(a: Matrix, x: Sequence[float]) -> Vector:
    return [sum(a[i][j] * x[j] for j in range(len(x))) for i in range(len(a))]


def dot(a: Sequence[float], b: Sequence[float]) -> float:
    return sum(a[i] * b[i] for i in range(len(a)))


def norm(x: Sequence[float]) -> float:
    return math.sqrt(sum(value * value for value in x))


def frobenius_norm(a: Matrix) -> float:
    return math.sqrt(sum(value * value for row in a for value in row))


def normalize(x: Sequence[float]) -> Vector:
    size = norm(x)
    if size < EPS:
        raise ValueError("Gặp vector gần bằng 0.")
    result = [value / size for value in x]
    return fix_vector_sign(result)


def fix_vector_sign(x: Sequence[float]) -> Vector:
    result = list(x)
    for value in result:
        if abs(value) > 1e-10:
            if value < 0:
                result = [-item for item in result]
            break
    return [clean_zero(value) for value in result]


def clean_zero(value: float, eps: float = 1e-10) -> float:
    return 0.0 if abs(value) < eps else value


def residual_norm(a: Matrix, eigenvalue: float, eigenvector: Sequence[float]) -> float:
    ax = mat_vec_mul(a, eigenvector)
    return norm([ax[i] - eigenvalue * eigenvector[i] for i in range(len(ax))])


def default_start_vector(n: int, shift: int = 0) -> Vector:
    return [float(((i + shift) % n) + 1) for i in range(n)]


def power_method_symmetric(a: Matrix, tol: float, max_iter: int, shift: int = 0) -> EigenStep:
    n = len(a)
    x = normalize(default_start_vector(n, shift))
    old_lambda = 0.0
    eigenvalue = 0.0
    residual = float("inf")

    for iteration in range(1, max_iter + 1):
        y = mat_vec_mul(a, x)
        if norm(y) < EPS:
            return EigenStep(0, 0.0, x, iteration, 0.0)
        x = normalize(y)
        ax = mat_vec_mul(a, x)
        eigenvalue = dot(x, ax)
        residual = residual_norm(a, eigenvalue, x)
        rel_change = abs(eigenvalue - old_lambda) / max(1.0, abs(eigenvalue))
        if iteration > 1 and rel_change <= tol and residual <= max(tol, tol * abs(eigenvalue)):
            break
        old_lambda = eigenvalue

    return EigenStep(0, clean_zero(eigenvalue), x, iteration, residual)


def rank_one_deflate(a: Matrix, eigenvalue: float, vector: Sequence[float]) -> Matrix:
    n = len(a)
    return [
        [clean_zero(a[i][j] - eigenvalue * vector[i] * vector[j]) for j in range(n)]
        for i in range(n)
    ]


def symmetric_eigen_by_power(a: Matrix, count: int, options: SVDOptions) -> List[EigenStep]:
    working = [row[:] for row in a]
    steps: List[EigenStep] = []

    for index in range(1, count + 1):
        step = power_method_symmetric(working, options.tol, options.max_iter, shift=index - 1)
        value = 0.0 if abs(step.eigenvalue) < options.rank_tol else step.eigenvalue
        if value < 0 and abs(value) < 100 * options.rank_tol:
            value = 0.0
        step = EigenStep(index, clean_zero(value), step.eigenvector, step.iterations, step.residual)
        steps.append(step)
        if abs(value) <= options.rank_tol:
            break
        working = rank_one_deflate(working, value, step.eigenvector)

    steps.sort(key=lambda item: item.eigenvalue, reverse=True)
    for index, step in enumerate(steps, start=1):
        step.index = index
    return steps


def diagonal(values: Sequence[float]) -> Matrix:
    n = len(values)
    return [[values[i] if i == j else 0.0 for j in range(n)] for i in range(n)]


def reduced_sigma(values: Sequence[float]) -> Matrix:
    return diagonal(values)


def pseudo_inverse_from_svd(singular_values: Sequence[float], left: Sequence[Vector], right: Sequence[Vector], n: int, m: int) -> Matrix:
    result = [[0.0 for _ in range(m)] for _ in range(n)]
    for sigma, u, v in zip(singular_values, left, right):
        if sigma < EPS:
            continue
        factor = 1.0 / sigma
        for i in range(n):
            for j in range(m):
                result[i][j] += factor * v[i] * u[j]
    return [[clean_zero(value) for value in row] for row in result]


def reconstruct_from_svd(singular_values: Sequence[float], left: Sequence[Vector], right: Sequence[Vector], m: int, n: int) -> Matrix:
    result = [[0.0 for _ in range(n)] for _ in range(m)]
    for sigma, u, v in zip(singular_values, left, right):
        for i in range(m):
            for j in range(n):
                result[i][j] += sigma * u[i] * v[j]
    return [[clean_zero(value) for value in row] for row in result]


def column_gram(vectors: Sequence[Vector]) -> Matrix:
    r = len(vectors)
    return [[clean_zero(dot(vectors[i], vectors[j])) for j in range(r)] for i in range(r)]


def condition_number(values: Sequence[float], rank_tol: float, total_count: int) -> float | None:
    positive = [value for value in values if value > rank_tol]
    if not positive:
        return None
    if len(positive) < total_count:
        return math.inf
    return max(positive) / min(positive)


def solve_svd(text: str) -> SVDSolution:
    a, options = parse_input(text)
    m, n = len(a), len(a[0])
    max_components = min(m, n)
    count = options.components if options.components is not None else max_components
    count = max(1, min(count, max_components))

    if m >= n:
        branch = "right"
        gram = mat_mul(transpose(a), a)
        eigen_steps = symmetric_eigen_by_power(gram, count, options)
        singular_values: Vector = []
        right_vectors: List[Vector] = []
        left_vectors: List[Vector] = []
        for step in eigen_steps:
            if step.eigenvalue <= options.rank_tol:
                continue
            sigma = math.sqrt(max(0.0, step.eigenvalue))
            v = normalize(step.eigenvector)
            av = mat_vec_mul(a, v)
            u = normalize([value / sigma for value in av])
            singular_values.append(clean_zero(sigma))
            right_vectors.append(v)
            left_vectors.append(u)
    else:
        branch = "left"
        gram = mat_mul(a, transpose(a))
        eigen_steps = symmetric_eigen_by_power(gram, count, options)
        singular_values = []
        right_vectors = []
        left_vectors = []
        at = transpose(a)
        for step in eigen_steps:
            if step.eigenvalue <= options.rank_tol:
                continue
            sigma = math.sqrt(max(0.0, step.eigenvalue))
            u = normalize(step.eigenvector)
            atu = mat_vec_mul(at, u)
            v = normalize([value / sigma for value in atu])
            singular_values.append(clean_zero(sigma))
            left_vectors.append(u)
            right_vectors.append(v)

    rank = len(singular_values)
    pinv = pseudo_inverse_from_svd(singular_values, left_vectors, right_vectors, n, m) if options.pinv else None
    cond = None
    if options.condition_number and count == max_components:
        cond = condition_number(singular_values, options.rank_tol, max_components)
    reconstruction = reconstruct_from_svd(singular_values, left_vectors, right_vectors, m, n)
    reconstruction_error = frobenius_norm(mat_sub(a, reconstruction))

    return SVDSolution(
        matrix=a,
        options=options,
        branch=branch,
        gram_matrix=gram,
        eigen_steps=eigen_steps,
        singular_values=singular_values,
        left_vectors=left_vectors,
        right_vectors=right_vectors,
        rank=rank,
        reduced_sigma=reduced_sigma(singular_values),
        pseudo_inverse=pinv,
        condition_number=cond,
        reconstruction=reconstruction,
        reconstruction_error=reconstruction_error,
    )


def format_real(value: float, digits: int = 7) -> str:
    if math.isinf(value):
        return "+\\infty"
    if abs(value) < 0.5 * 10 ** (-digits):
        value = 0.0
    return f"{value:.{digits}f}"


def matrix_latex(a: Matrix, digits: int = 4) -> str:
    rows = [" & ".join(format_real(value, digits) for value in row) for row in a]
    return "\\begin{pmatrix}\n" + "\\\\\n".join(rows) + "\n\\end{pmatrix}"


def vector_latex(v: Sequence[float], digits: int = 7) -> str:
    return "\\begin{pmatrix}\n" + "\\\\\n".join(format_real(value, digits) for value in v) + "\n\\end{pmatrix}"


def vectors_as_matrix(vectors: Sequence[Vector], rows: int) -> Matrix:
    if not vectors:
        return [[] for _ in range(rows)]
    return [[vectors[j][i] for j in range(len(vectors))] for i in range(rows)]


def render_markdown(solution: SVDSolution) -> str:
    a = solution.matrix
    m, n = len(a), len(a[0])
    max_components = min(m, n)
    requested_components = solution.options.components if solution.options.components is not None else max_components
    requested_components = max(1, min(requested_components, max_components))
    relation = "=" if requested_components == max_components else "\\approx"
    relation_text = "bằng" if relation == "=" else "xấp xỉ"
    lines: List[str] = []
    lines.append("# Lời giải khai triển kỳ dị SVD")
    lines.append("")
    lines.append("> Các phép tính bên trong dùng giá trị gốc. Ma trận trình bày 4 chữ số thập phân; các số liệu khác lấy 7 chữ số thập phân.")
    lines.append("")
    lines.append("## Bài toán")
    lines.append("")
    lines.append("Cho ma trận")
    lines.append("")
    lines.append("$$")
    lines.append("A=" + matrix_latex(a, 4))
    lines.append("$$")
    lines.append("")
    lines.append(f"Ở đây \\(A\\in\\mathbb{{R}}^{{{m}\\times {n}}}\\).")
    lines.append("")

    lines.append("## Bước 1. Quy SVD về bài toán trị riêng")
    lines.append("")
    if solution.branch == "right":
        lines.append("Vì \\(m\\ge n\\), ta lập")
        lines.append("")
        lines.append("$$")
        lines.append("B=A^TA.")
        lines.append("$$")
        lines.append("")
        lines.append("Khi đó các véc tơ kỳ dị phải \\(v_i\\) là véc tơ riêng của \\(B\\):")
        lines.append("")
        lines.append("$$")
        lines.append("Bv_i=\\lambda_i v_i,\\qquad \\sigma_i=\\sqrt{\\lambda_i},\\qquad u_i=\\frac{Av_i}{\\sigma_i}.")
        lines.append("$$")
    else:
        lines.append("Vì \\(m<n\\), ta lập")
        lines.append("")
        lines.append("$$")
        lines.append("C=AA^T.")
        lines.append("$$")
        lines.append("")
        lines.append("Khi đó các véc tơ kỳ dị trái \\(u_i\\) là véc tơ riêng của \\(C\\):")
        lines.append("")
        lines.append("$$")
        lines.append("Cu_i=\\lambda_i u_i,\\qquad \\sigma_i=\\sqrt{\\lambda_i},\\qquad v_i=\\frac{A^Tu_i}{\\sigma_i}.")
        lines.append("$$")
    lines.append("")
    lines.append("Ta có")
    lines.append("")
    lines.append("$$")
    matrix_name = "B" if solution.branch == "right" else "C"
    lines.append(matrix_name + "=" + matrix_latex(solution.gram_matrix, 4))
    lines.append("$$")
    lines.append("")

    lines.append("## Bước 2. Tìm trị riêng bằng phương pháp lũy thừa và xuống thang")
    lines.append("")
    lines.append("Áp dụng phương pháp lũy thừa cho ma trận đối xứng ở trên, sau mỗi trị riêng thì xuống thang")
    lines.append("")
    lines.append("Trước khi xuống thang, véc tơ riêng được chuẩn hóa theo chuẩn 2:")
    lines.append("")
    lines.append("$$")
    lines.append("q_k\\leftarrow \\frac{q_k}{\\|q_k\\|_2},\\qquad \\|q_k\\|_2=1.")
    lines.append("$$")
    lines.append("")
    lines.append("Do đó công thức xuống thang dùng trong bài là")
    lines.append("")
    lines.append("$$")
    lines.append("M_{k+1}=M_k-\\lambda_k q_kq_k^T")
    lines.append("$$")
    lines.append("")
    lines.append("để tìm trị riêng tiếp theo.")
    lines.append("")
    lines.append("Sai số dư trong bảng được tính bởi")
    lines.append("")
    lines.append("$$")
    lines.append("r_i=\\|Mq_i-\\lambda_iq_i\\|_2.")
    lines.append("$$")
    lines.append("")
    lines.append("Nếu \\(r_i\\) càng nhỏ thì cặp \\((\\lambda_i,q_i)\\) càng gần với trị riêng và véc tơ riêng đúng.")
    lines.append("")
    lines.append("| \\(i\\) | \\(\\lambda_i\\) | \\(\\sigma_i=\\sqrt{\\lambda_i}\\) | Số lặp | Sai số dư |")
    lines.append("|---:|---:|---:|---:|---:|")
    for index, step in enumerate(solution.eigen_steps, start=1):
        sigma = math.sqrt(max(0.0, step.eigenvalue)) if step.eigenvalue > solution.options.rank_tol else 0.0
        lines.append(
            f"| {index} | {format_real(step.eigenvalue, 7)} | {format_real(sigma, 7)} | {step.iterations} | {format_real(step.residual, 7)} |"
        )
    lines.append("")

    lines.append("## Bước 3. Giá trị kỳ dị và hạng")
    lines.append("")
    lines.append("Với mỗi trị riêng của ma trận đối xứng vừa xét:")
    lines.append("")
    lines.append("$$")
    lines.append("\\lambda_i>0\\Rightarrow \\sigma_i=\\sqrt{\\lambda_i}>0,\\qquad \\lambda_i=0\\Rightarrow \\sigma_i=0.")
    lines.append("$$")
    lines.append("")
    lines.append("Số giá trị kỳ dị khác 0 chính là hạng của ma trận.")
    lines.append("")
    if solution.singular_values:
        lines.append("Các giá trị kỳ dị khác 0 là")
        lines.append("")
        lines.append("$$")
        lines.append(
            ",\\quad ".join(
                f"\\sigma_{index}={format_real(value, 7)}"
                for index, value in enumerate(solution.singular_values, start=1)
            )
            + "."
        )
        lines.append("$$")
    else:
        lines.append("Không tìm được giá trị kỳ dị khác 0 theo ngưỡng đã chọn.")
    lines.append("")
    lines.append(f"Suy ra \\(\\operatorname{{rank}}(A)={solution.rank}\\).")
    lines.append("")

    lines.append("## Bước 4. Vector kỳ dị")
    lines.append("")
    for index, sigma in enumerate(solution.singular_values, start=1):
        u = solution.left_vectors[index - 1]
        v = solution.right_vectors[index - 1]
        lines.append(f"Với \\(\\sigma_{index}={format_real(sigma, 7)}\\), ta có")
        lines.append("")
        lines.append("$$")
        lines.append("v_" + str(index) + "\\approx " + vector_latex(v, 7) + ",\\qquad u_" + str(index) + "\\approx " + vector_latex(u, 7) + ".")
        lines.append("$$")
        lines.append("")
        lines.append("Các véc tơ này đã được chuẩn hóa:")
        lines.append("")
        lines.append("$$")
        lines.append("\\|v_" + str(index) + "\\|_2=" + format_real(norm(v), 7) + ",\\qquad \\|u_" + str(index) + "\\|_2=" + format_real(norm(u), 7) + ".")
        lines.append("$$")
        lines.append("")

    lines.append("## Bước 5. Khai triển SVD rút gọn")
    lines.append("")
    u_matrix = vectors_as_matrix(solution.left_vectors, m)
    v_matrix = vectors_as_matrix(solution.right_vectors, n)
    lines.append("Đặt")
    lines.append("")
    lines.append("$$")
    lines.append("U_r=" + matrix_latex(u_matrix, 4))
    lines.append("$$")
    lines.append("")
    lines.append("$$")
    lines.append("\\Sigma_r=" + matrix_latex(solution.reduced_sigma, 4))
    lines.append("$$")
    lines.append("")
    lines.append("$$")
    lines.append("V_r=" + matrix_latex(v_matrix, 4))
    lines.append("$$")
    lines.append("")
    lines.append("Khi đó")
    lines.append("")
    lines.append("$$")
    lines.append("A" + relation + "U_r\\Sigma_rV_r^T.")
    lines.append("$$")
    lines.append("")
    if solution.rank:
        lines.append("Tương đương")
        lines.append("")
        lines.append("$$")
        lines.append(
            "A"
            + relation
            + "+".join(
                f"\\sigma_{index}u_{index}v_{index}^T"
                for index in range(1, solution.rank + 1)
            )
            + "."
        )
        lines.append("$$")
        lines.append("")

    lines.append("## Kiểm tra tái tạo")
    lines.append("")
    lines.append("Từ các thành phần SVD đã tìm được, ta tính lại")
    lines.append("")
    lines.append("$$")
    lines.append("\\widehat A=U_r\\Sigma_rV_r^T=" + matrix_latex(solution.reconstruction, 4))
    lines.append("$$")
    lines.append("")
    lines.append("Sai số tái tạo theo chuẩn Frobenius là")
    lines.append("")
    lines.append("$$")
    lines.append("E=\\|A-\\widehat A\\|_F=" + format_real(solution.reconstruction_error, 7) + ".")
    lines.append("$$")
    lines.append("")
    if relation == "=":
        lines.append("Vì đã lấy đủ số thành phần cần thiết, về lý thuyết ta có khai triển đúng của ma trận.")
    else:
        lines.append("Vì chỉ lấy một số thành phần đầu, đây là khai triển xấp xỉ của ma trận.")
    lines.append("")

    if solution.options.mode == "full":
        lines.append("## Kiểm tra trực chuẩn")
        lines.append("")
        lines.append("Ta kiểm tra các cột của \\(U_r\\) và \\(V_r\\):")
        lines.append("")
        lines.append("$$")
        lines.append("U_r^TU_r=" + matrix_latex(column_gram(solution.left_vectors), 4))
        lines.append("$$")
        lines.append("")
        lines.append("$$")
        lines.append("V_r^TV_r=" + matrix_latex(column_gram(solution.right_vectors), 4))
        lines.append("$$")
        lines.append("")

    if solution.pseudo_inverse is not None:
        lines.append("## Nghịch đảo suy rộng")
        lines.append("")
        lines.append("Theo công thức Moore-Penrose,")
        lines.append("")
        lines.append("$$")
        lines.append("A^\\dagger=V_r\\Sigma_r^{-1}U_r^T.")
        lines.append("$$")
        lines.append("")
        lines.append("Suy ra")
        lines.append("")
        lines.append("$$")
        lines.append("A^\\dagger=" + matrix_latex(solution.pseudo_inverse, 4))
        lines.append("$$")
        lines.append("")

    if solution.condition_number is not None:
        lines.append("## Số điều kiện")
        lines.append("")
        lines.append("Nếu các giá trị kỳ dị khác 0 phủ đủ số chiều cần xét thì")
        lines.append("")
        lines.append("$$")
        lines.append("\\operatorname{cond}(A)=\\frac{\\sigma_{\\max}}{\\sigma_{\\min}}.")
        lines.append("$$")
        lines.append("")
        lines.append("Nếu tồn tại \\(\\sigma_i=0\\), ma trận bị suy biến theo nghĩa SVD và số điều kiện bằng \\(+\\infty\\).")
        lines.append("")
        lines.append("Trong bài này ta có")
        lines.append("")
        lines.append("$$")
        lines.append("\\operatorname{cond}(A)=\\frac{\\sigma_{\\max}}{\\sigma_{\\min}}=" + format_real(solution.condition_number, 7) + ".")
        lines.append("$$")
        lines.append("")

    lines.append("## Kết luận")
    lines.append("")
    lines.append("Khai triển kỳ dị rút gọn của ma trận là")
    lines.append("")
    lines.append("$$")
    lines.append("A" + relation + "U_r\\Sigma_rV_r^T")
    lines.append("$$")
    lines.append("")
    lines.append("với các ma trận \\(U_r,\\Sigma_r,V_r\\) đã tính ở trên; sai số tái tạo là \\(E=" + format_real(solution.reconstruction_error, 7) + "\\).")
    lines.append("")
    return "\n".join(lines)


def solve_text_to_markdown(text: str) -> str:
    return render_markdown(solve_svd(text))


def save_markdown_from_text(text: str, output_path: str) -> None:
    markdown = solve_text_to_markdown(text)
    with open(output_path, "w", encoding="utf-8") as file:
        file.write(markdown)

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import List, Sequence


EPS = 1e-12

Matrix = List[List[float]]
Vector = List[complex]


@dataclass
class IterationStep:
    k: int
    x: Vector
    y: Vector
    pivot_index: int
    scale: complex
    next_x: Vector
    lambda_estimate: complex
    residual_norm: float


@dataclass
class SolverOptions:
    mode: str = "exam"
    show_special_checks: str = "auto"


@dataclass
class EvenSquareAnalysis:
    from_power: int
    to_power: int
    component_index: int
    z_from: Vector
    z_to: Vector
    value: complex


@dataclass
class QuadraticAnalysis:
    n_index: int
    row_indices: tuple[int, int]
    z_n: Vector
    z_n1: Vector
    z_n2: Vector
    p: complex
    q: complex
    roots: tuple[complex, complex]


@dataclass
class ClassificationReport:
    normal_error: float | None
    period2_error: float | None
    normal_threshold: float
    period2_threshold: float
    even_square_value: complex | None
    case: str
    reason: str


@dataclass
class PowerResult:
    matrix: Matrix
    initial_vector: Vector
    tol: float
    max_iter: int
    steps: List[IterationStep]
    eigenvalue: complex
    eigenvector: Vector
    converged: bool
    ratio_estimates: List[complex]
    even_square_estimate: complex | None
    even_square_analysis: EvenSquareAnalysis | None
    quadratic_coefficients: tuple[complex, complex] | None
    quadratic_roots: tuple[complex, complex] | None
    quadratic_analysis: QuadraticAnalysis | None
    case: str
    classification: ClassificationReport


@dataclass
class DeflationResult:
    index: int
    matrix: Matrix
    eigenvalue: complex
    right_vector: Vector
    left_vector: Vector
    denominator: complex
    deflated_matrix: Matrix
    next_result: PowerResult | None
    source_case: str


@dataclass
class DominantEigenSolution:
    matrix: Matrix
    initial_vector: Vector
    tol: float
    max_iter: int
    deflation_count: int
    options: SolverOptions
    main_result: PowerResult
    deflations: List[DeflationResult]


def parse_input(text: str) -> tuple[Matrix, Vector, float, int, int, SolverOptions]:
    matrix_rows: Matrix = []
    x0: Vector = []
    tol = 1e-7
    max_iter = 100
    deflation_count = 1
    mode = "exam"
    show_special_checks = "auto"
    section = "matrix"

    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        lower = line.lower().replace(" ", "")

        if lower in {"matrix:", "a:"}:
            section = "matrix"
            continue
        if lower in {"x0:", "vector:", "vecto:", "vectơ:"}:
            section = "x0"
            continue
        if lower.startswith("tol:") or lower.startswith("eps:"):
            tol = float(line.split(":", 1)[1].strip())
            continue
        if lower.startswith("max_iter:") or lower.startswith("maxiter:"):
            max_iter = int(float(line.split(":", 1)[1].strip()))
            continue
        if lower.startswith("deflation_count:") or lower.startswith("count:"):
            deflation_count = int(float(line.split(":", 1)[1].strip()))
            continue
        if lower.startswith("mode:"):
            mode = line.split(":", 1)[1].strip().lower()
            continue
        if lower.startswith("show_special_checks:") or lower.startswith("special:"):
            show_special_checks = line.split(":", 1)[1].strip().lower()
            continue

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
        if section == "matrix":
            matrix_rows.append(values)
        else:
            x0.extend(complex(value) for value in values)

    if not matrix_rows:
        raise ValueError("Chưa nhập ma trận.")

    n = len(matrix_rows)
    if any(len(row) != n for row in matrix_rows):
        raise ValueError("Ma trận phải vuông.")

    if not x0:
        x0 = [1 + 0j for _ in range(n)]
    if len(x0) != n:
        raise ValueError("Số phần tử của x0 phải bằng cấp của ma trận.")
    if vector_norm(x0) < EPS:
        raise ValueError("Vector ban đầu x0 phải khác 0.")

    if mode not in {"exam", "full"}:
        raise ValueError("mode chỉ nhận 'exam' hoặc 'full'.")
    if show_special_checks not in {"auto", "always", "never"}:
        raise ValueError("show_special_checks chỉ nhận 'auto', 'always' hoặc 'never'.")

    return matrix_rows, x0, tol, max_iter, deflation_count, SolverOptions(mode, show_special_checks)


def default_start_vector(n: int) -> Vector:
    return [complex(i + 1) for i in range(n)]


def copy_matrix(a: Matrix) -> Matrix:
    return [row[:] for row in a]


def transpose(a: Matrix) -> Matrix:
    return [[a[i][j] for i in range(len(a))] for j in range(len(a))]


def mat_vec_mul(a: Matrix, x: Sequence[complex]) -> Vector:
    return [sum(a[i][j] * x[j] for j in range(len(x))) for i in range(len(a))]


def mat_sub_rank_one(a: Matrix, scalar: complex, v: Sequence[complex], w: Sequence[complex]) -> Matrix:
    n = len(a)
    result: Matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            value = complex(a[i][j]) - scalar * v[i] * w[j]
            row.append(value.real if abs(value.imag) < 1e-10 else value)
        result.append(row)  # type: ignore[arg-type]
    return result


def dot_plain(a: Sequence[complex], b: Sequence[complex]) -> complex:
    return sum(a[i] * b[i] for i in range(len(a)))


def vector_norm(x: Sequence[complex]) -> float:
    return math.sqrt(sum(abs(value) ** 2 for value in x))


def normalize_by_l2(x: Sequence[complex]) -> Vector:
    size = vector_norm(x)
    if size < EPS:
        raise ValueError("Không thể chuẩn hóa vector gần bằng 0.")
    return [clean_zero(value / size) for value in x]




def vector_distance(a: Sequence[complex], b: Sequence[complex]) -> float:
    diff = [a[i] - b[i] for i in range(len(a))]
    return vector_norm(diff) / max(1.0, vector_norm(a), vector_norm(b))


def dominant_vector_distance(a: Sequence[complex], b: Sequence[complex]) -> float:
    max_abs = max(max(abs(value) for value in a), max(abs(value) for value in b))
    if max_abs < EPS:
        return 0.0

    indices = [
        i
        for i in range(len(a))
        if max(abs(a[i]), abs(b[i])) >= 0.1 * max_abs
    ]
    if not indices:
        return vector_distance(a, b)

    diff_norm = math.sqrt(sum(abs(a[i] - b[i]) ** 2 for i in indices))
    base_norm = max(
        1.0,
        math.sqrt(sum(abs(a[i]) ** 2 for i in indices)),
        math.sqrt(sum(abs(b[i]) ** 2 for i in indices)),
    )
    return diff_norm / base_norm


def normalize_by_max_component(x: Sequence[complex]) -> tuple[Vector, complex, int]:
    pivot = max(range(len(x)), key=lambda i: abs(x[i]))
    scale = x[pivot]
    if abs(scale) < EPS:
        raise ValueError("Gặp vector gần bằng 0 trong quá trình lặp.")
    return [value / scale for value in x], scale, pivot


def residual_norm(a: Matrix, eigenvalue: complex, eigenvector: Sequence[complex]) -> float:
    ax = mat_vec_mul(a, eigenvector)
    residual = [ax[i] - eigenvalue * eigenvector[i] for i in range(len(a))]
    return vector_norm(residual)


def rayleigh_estimate(a: Matrix, x: Sequence[complex]) -> complex:
    ax = mat_vec_mul(a, x)
    numerator = dot_plain(x, ax)
    denominator = dot_plain(x, x)
    if abs(denominator) < EPS:
        return 0j
    return numerator / denominator


def raw_powers(a: Matrix, x0: Sequence[complex], count: int) -> List[Vector]:
    powers = [[complex(value) for value in x0]]
    for _ in range(count):
        powers.append(mat_vec_mul(a, powers[-1]))
    return powers


def component_ratios(previous: Sequence[complex], current: Sequence[complex]) -> List[complex]:
    ratios = []
    for old, new in zip(previous, current):
        if abs(old) > EPS:
            ratios.append(new / old)
    return ratios


def dominant_component_ratios(previous: Sequence[complex], current: Sequence[complex]) -> List[complex]:
    max_previous = max(abs(value) for value in previous)
    if max_previous < EPS:
        return []
    ratios = []
    for old, new in zip(previous, current):
        if abs(old) >= 0.1 * max_previous:
            ratios.append(new / old)
    return ratios


def average(values: Sequence[complex]) -> complex:
    if not values:
        return 0j
    return sum(values) / len(values)


def estimate_even_square(powers: Sequence[Vector]) -> complex | None:
    analysis = estimate_even_square_analysis(powers)
    return analysis.value if analysis is not None else None


def estimate_even_square_analysis(powers: Sequence[Vector]) -> EvenSquareAnalysis | None:
    if len(powers) < 5:
        return None
    even_last = len(powers) - 1
    if even_last % 2 == 1:
        even_last -= 1
    if even_last < 2:
        return None
    z_from = powers[even_last - 2]
    z_to = powers[even_last]
    max_previous = max(abs(value) for value in z_from)
    if max_previous < EPS:
        return None
    component = max(range(len(z_from)), key=lambda i: abs(z_from[i]))
    value = z_to[component] / z_from[component]
    return EvenSquareAnalysis(
        from_power=even_last - 2,
        to_power=even_last,
        component_index=component,
        z_from=[clean_zero(value) for value in z_from],
        z_to=[clean_zero(value) for value in z_to],
        value=clean_zero(value),
    )


def solve_quadratic_from_powers(powers: Sequence[Vector]) -> QuadraticAnalysis | None:
    if len(powers) < 3:
        return None
    z0, z1, z2 = powers[-3], powers[-2], powers[-1]
    n_index = len(powers) - 3
    n = len(z0)
    best_pair = None
    best_det = 0.0

    for r in range(n):
        for s in range(r + 1, n):
            det = z1[r] * (-z0[s]) - z1[s] * (-z0[r])
            if abs(det) > best_det:
                best_det = abs(det)
                best_pair = (r, s, det)

    if best_pair is None or best_det < EPS:
        return None

    r, s, det = best_pair
    p = (z2[r] * (-z0[s]) - z2[s] * (-z0[r])) / det
    q = (z1[r] * z2[s] - z1[s] * z2[r]) / det
    delta = p * p - 4 * q
    roots = ((p + cmath.sqrt(delta)) / 2, (p - cmath.sqrt(delta)) / 2)
    return QuadraticAnalysis(
        n_index=n_index,
        row_indices=(r, s),
        z_n=[clean_zero(value) for value in z0],
        z_n1=[clean_zero(value) for value in z1],
        z_n2=[clean_zero(value) for value in z2],
        p=clean_zero(p),
        q=clean_zero(q),
        roots=(clean_zero(roots[0]), clean_zero(roots[1])),
    )


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

    if vector_norm(vector) > EPS:
        vector = normalize_by_l2(vector)
    return [clean_zero(value) for value in vector]


def solve_square_system(a: List[List[complex]], b: List[complex]) -> Vector:
    n = len(a)
    aug = [a[i][:] + [b[i]] for i in range(n)]

    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot][col]) < 1e-10:
            raise ValueError("Singular system")
        aug[col], aug[pivot] = aug[pivot], aug[col]

        div = aug[col][col]
        aug[col] = [value / div for value in aug[col]]

        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            if abs(factor) > 1e-12:
                aug[row] = [aug[row][j] - factor * aug[col][j] for j in range(n + 1)]

    return [clean_zero(row[-1]) for row in aug]


def approximate_null_vector(a: List[List[complex]]) -> Vector:
    n = len(a)
    best_vector: Vector | None = None
    best_residual = float("inf")

    for free_col in range(n):
        columns = [col for col in range(n) if col != free_col]
        for skipped_row in range(n):
            rows = [row for row in range(n) if row != skipped_row]
            square = [[a[row][col] for col in columns] for row in rows]
            rhs = [-a[row][free_col] for row in rows]
            try:
                solution = solve_square_system(square, rhs)
            except ValueError:
                continue

            candidate = [0j for _ in range(n)]
            candidate[free_col] = 1 + 0j
            for col, value in zip(columns, solution):
                candidate[col] = value

            residual = vector_norm(
                [sum(a[i][j] * candidate[j] for j in range(n)) for i in range(n)]
            )
            if residual < best_residual:
                best_residual = residual
                best_vector = candidate

    if best_vector is None:
        return [0j for _ in range(n)]

    if vector_norm(best_vector) > EPS:
        best_vector = normalize_by_l2(best_vector)
    return [clean_zero(value) for value in best_vector]


def eigenvector_for_matrix(a: Matrix, eigenvalue: complex) -> Vector:
    n = len(a)
    matrix = [
        [
            complex(a[i][j]) - (eigenvalue if i == j else 0j)
            for j in range(n)
        ]
        for i in range(n)
    ]
    vector = nullspace_vector(matrix)
    if vector_norm(vector) < EPS:
        vector = approximate_null_vector(matrix)
    return normalize_by_l2(vector)


def dominant_roots_from_result(result: PowerResult) -> tuple[complex, complex] | None:
    if result.quadratic_roots:
        return result.quadratic_roots
    if result.case == "opposite" and result.even_square_estimate is not None:
        root = cmath.sqrt(result.even_square_estimate)
        return clean_zero(root), clean_zero(-root)
    return None


def eigenpair_for_deflation(result: PowerResult, matrix: Matrix) -> tuple[complex, Vector, str] | None:
    if result.converged:
        return result.eigenvalue, result.eigenvector, "TH1"

    roots = dominant_roots_from_result(result)
    if roots:
        eigenvalue = max(roots, key=lambda value: abs(value))
        eigenvector = eigenvector_for_matrix(matrix, eigenvalue)
        if result.case == "opposite":
            return eigenvalue, eigenvector, "TH2"
        if result.case == "quadratic":
            return eigenvalue, eigenvector, "TH3"

    return None


def classify_power_case(
    steps: Sequence[IterationStep],
    converged: bool,
    roots: tuple[complex, complex] | None,
    even_square: complex | None,
    tol: float,
) -> ClassificationReport:
    normal_threshold = max(1e-7, 10 * tol)
    period2_threshold = max(1e-5, 1000 * tol)

    normal_error: float | None = None
    if len(steps) >= 2:
        last = steps[-1]
        prev = steps[-2]
        lambda_error = abs(last.lambda_estimate - prev.lambda_estimate) / max(
            1.0, abs(last.lambda_estimate)
        )
        residual_error = last.residual_norm / max(1.0, abs(last.lambda_estimate))
        normal_error = max(lambda_error, residual_error)

    period2_error: float | None = None
    if len(steps) >= 3:
        period2_error = dominant_vector_distance(steps[-1].next_x, steps[-3].next_x)

    if converged:
        return ClassificationReport(
            normal_error=normal_error,
            period2_error=period2_error,
            normal_threshold=normal_threshold,
            period2_threshold=period2_threshold,
            even_square_value=even_square,
            case="normal",
            reason="TH1: dãy lặp hội tụ trực tiếp.",
        )

    is_period2 = period2_error is not None and period2_error <= period2_threshold
    even_square_is_positive_real = (
        even_square is not None
        and abs(even_square.imag) <= 1e-7 * max(1.0, abs(even_square.real))
        and even_square.real > 0
    )
    if is_period2 and even_square_is_positive_real:
        return ClassificationReport(
            normal_error=normal_error,
            period2_error=period2_error,
            normal_threshold=normal_threshold,
            period2_threshold=period2_threshold,
            even_square_value=even_square,
            case="opposite",
            reason="TH2: dãy không hội tụ trực tiếp nhưng hội tụ theo chu kỳ 2 và \\(\\lambda^2>0\\).",
        )

    if roots is not None:
        return ClassificationReport(
            normal_error=normal_error,
            period2_error=period2_error,
            normal_threshold=normal_threshold,
            period2_threshold=period2_threshold,
            even_square_value=even_square,
            case="quadratic",
            reason="TH3: TH1 và TH2 không thỏa, lập hệ bậc hai để tìm cặp trị riêng trội.",
        )

    return ClassificationReport(
        normal_error=normal_error,
        period2_error=period2_error,
        normal_threshold=normal_threshold,
        period2_threshold=period2_threshold,
        even_square_value=even_square,
        case="undetermined",
        reason="Chưa đủ điều kiện kết luận sau các kiểm tra.",
    )


def power_method(a: Matrix, x0: Sequence[complex], tol: float = 1e-7, max_iter: int = 30) -> PowerResult:
    x, _, _ = normalize_by_max_component(x0)
    steps: List[IterationStep] = []
    converged = False
    eigenvalue = 0j

    for k in range(max_iter):
        y = mat_vec_mul(a, x)
        next_x, scale, pivot_index = normalize_by_max_component(y)
        eigenvalue = scale
        residual = residual_norm(a, eigenvalue, next_x)
        steps.append(
            IterationStep(
                k=k,
                x=x,
                y=y,
                pivot_index=pivot_index,
                scale=scale,
                next_x=next_x,
                lambda_estimate=eigenvalue,
                residual_norm=residual,
            )
        )
        if k > 0:
            previous = steps[-2].lambda_estimate
            if abs(eigenvalue - previous) <= tol and residual <= max(tol, 10 * tol * abs(eigenvalue)):
                converged = True
                break
        x = next_x

    eigenvector = steps[-1].next_x if steps else x
    eigenvector = normalize_by_l2(eigenvector)
    if converged:
        eigenvalue = rayleigh_estimate(a, eigenvector)
    elif steps:
        eigenvalue = steps[-1].lambda_estimate
    if steps and converged:
        steps[-1].lambda_estimate = eigenvalue

    powers = raw_powers(a, x0, max(6, min(max_iter, 20)))
    ratios = component_ratios(powers[-2], powers[-1])
    even_square_analysis = estimate_even_square_analysis(powers)
    even_square = even_square_analysis.value if even_square_analysis is not None else None
    quadratic_analysis = solve_quadratic_from_powers(powers)
    if quadratic_analysis is None:
        coefficients = None
        roots = None
    else:
        coefficients = (quadratic_analysis.p, quadratic_analysis.q)
        roots = quadratic_analysis.roots
    classification = classify_power_case(steps, converged, roots, even_square, tol)
    case = classification.case

    return PowerResult(
        matrix=copy_matrix(a),
        initial_vector=[complex(value) for value in x0],
        tol=tol,
        max_iter=max_iter,
        steps=steps,
        eigenvalue=clean_zero(eigenvalue),
        eigenvector=[clean_zero(value) for value in eigenvector],
        converged=converged,
        ratio_estimates=[clean_zero(value) for value in ratios],
        even_square_estimate=clean_zero(even_square) if even_square is not None else None,
        even_square_analysis=even_square_analysis,
        quadratic_coefficients=tuple(clean_zero(v) for v in coefficients) if coefficients else None,
        quadratic_roots=tuple(clean_zero(v) for v in roots) if roots else None,
        quadratic_analysis=quadratic_analysis,
        case=case,
        classification=classification,
    )


def solve_dominant_eigen(text: str) -> DominantEigenSolution:
    matrix, x0, tol, max_iter, deflation_count, options = parse_input(text)
    main_result = power_method(matrix, x0, tol, max_iter)
    deflations: List[DeflationResult] = []
    current = copy_matrix(matrix)
    current_x0 = x0[:]
    current_result = main_result

    for index in range(1, max(1, deflation_count)):
        pair = eigenpair_for_deflation(current_result, current)
        if pair is None:
            break
        eigenvalue, right_vector, source_case = pair
        left_vector = eigenvector_for_matrix(transpose(current), eigenvalue)
        denominator = dot_plain(left_vector, right_vector)
        if abs(denominator) < EPS:
            break
        scalar = eigenvalue / denominator
        next_matrix = mat_sub_rank_one(current, scalar, right_vector, left_vector)
        next_result = power_method(next_matrix, default_start_vector(len(matrix)), tol, max_iter)
        deflations.append(
            DeflationResult(
                index=index,
                matrix=copy_matrix(current),
                eigenvalue=clean_zero(eigenvalue),
                right_vector=[clean_zero(value) for value in right_vector],
                left_vector=[clean_zero(value) for value in left_vector],
                denominator=clean_zero(denominator),
                deflated_matrix=next_matrix,
                next_result=next_result,
                source_case=source_case,
            )
        )
        current = next_matrix
        current_x0 = default_start_vector(len(matrix))
        current_result = next_result

    return DominantEigenSolution(
        matrix=matrix,
        initial_vector=x0,
        tol=tol,
        max_iter=max_iter,
        deflation_count=deflation_count,
        options=options,
        main_result=main_result,
        deflations=deflations,
    )


def clean_zero(value: complex | None, eps: float = 1e-9) -> complex:
    if value is None:
        return 0j
    real = 0.0 if abs(value.real) < eps else value.real
    imag = 0.0 if abs(value.imag) < eps else value.imag
    return complex(real, imag)


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
        rows.append(" & ".join(format_complex(complex(value), digits) for value in row))
    return "\\begin{pmatrix}\n" + "\\\\\n".join(rows) + "\n\\end{pmatrix}"


def vector_latex(v: Sequence[complex], digits: int = 7) -> str:
    rows = "\\\\\n".join(format_complex(value, digits) for value in v)
    return "\\begin{pmatrix}\n" + rows + "\n\\end{pmatrix}"


def render_iteration_table(steps: Sequence[IterationStep]) -> List[str]:
    lines = []
    lines.append("| \\(k\\) | \\(y^{(k)}=Ax^{(k)}\\) | Chuẩn hóa | \\(x^{(k+1)}\\) | \\(\\lambda^{(k)}\\) | Sai số |")
    lines.append("|---:|---|---|---|---:|---:|")
    for step in steps:
        pivot_label = f"\\(s={step.pivot_index + 1},\\ y_s={format_complex(step.scale, 7)}\\)"
        lines.append(
            "| "
            + str(step.k + 1)
            + " | "
            + "\\("
            + vector_inline(step.y, 7)
            + "\\)"
            + " | "
            + pivot_label
            + " | "
            + "\\("
            + vector_inline(step.next_x, 7)
            + "\\)"
            + " | "
            + format_complex(step.lambda_estimate, 7)
            + " | "
            + format_real(step.residual_norm, 7)
            + " |"
        )
    return lines


def vector_inline(v: Sequence[complex], digits: int = 7) -> str:
    return "(" + ", ".join(format_complex(value, digits) for value in v) + ")^T"


def render_even_square_analysis(analysis: EvenSquareAnalysis) -> List[str]:
    i = analysis.component_index + 1
    lines: List[str] = []
    lines.append("Dãy lặp không hội tụ trực tiếp mà có dạng dao động. Ta xét lũy thừa chẵn.")
    lines.append("")
    lines.append("Ta có")
    lines.append("")
    lines.append("$$")
    lines.append(f"z_{{{analysis.from_power}}}=A^{{{analysis.from_power}}}x^{{(0)}}=" + vector_latex(analysis.z_from, 7))
    lines.append("$$")
    lines.append("")
    lines.append("$$")
    lines.append(f"z_{{{analysis.to_power}}}=A^{{{analysis.to_power}}}x^{{(0)}}=" + vector_latex(analysis.z_to, 7))
    lines.append("$$")
    lines.append("")
    lines.append(f"Chọn thành phần thứ {i}, suy ra")
    lines.append("")
    lines.append("$$")
    lines.append(
        "\\lambda^2\\approx "
        + f"\\frac{{(z_{{{analysis.to_power}}})_{i}}}{{(z_{{{analysis.from_power}}})_{i}}}"
        + "="
        + f"\\frac{{{format_complex(analysis.z_to[analysis.component_index], 7)}}}{{{format_complex(analysis.z_from[analysis.component_index], 7)}}}"
        + "="
        + format_complex(analysis.value, 7)
        + "."
    )
    lines.append("$$")
    lines.append("")
    lines.append("Do đó")
    lines.append("")
    lines.append("$$")
    lines.append(
        "\\lambda_{1,2}\\approx \\pm\\sqrt{"
        + format_complex(analysis.value, 7)
        + "}."
    )
    lines.append("$$")
    lines.append("")
    return lines


def render_quadratic_analysis(analysis: QuadraticAnalysis, matrix: Matrix) -> List[str]:
    r, s = analysis.row_indices
    rows = [r, s]
    lines: List[str] = []
    lines.append("Đặt")
    lines.append("")
    lines.append("$$")
    lines.append("z_k=A^kx^{(0)}.")
    lines.append("$$")
    lines.append("")
    lines.append("Ta lấy ba véc tơ liên tiếp:")
    lines.append("")
    lines.append("$$")
    lines.append(f"z_{{{analysis.n_index}}}=" + vector_latex(analysis.z_n, 7))
    lines.append("$$")
    lines.append("")
    lines.append("$$")
    lines.append(f"z_{{{analysis.n_index + 1}}}=" + vector_latex(analysis.z_n1, 7))
    lines.append("$$")
    lines.append("")
    lines.append("$$")
    lines.append(f"z_{{{analysis.n_index + 2}}}=" + vector_latex(analysis.z_n2, 7))
    lines.append("$$")
    lines.append("")
    lines.append("Theo bài giảng, với cặp trị riêng trội \\(\\lambda_1,\\lambda_2\\), đặt")
    lines.append("")
    lines.append("$$")
    lines.append("p=\\lambda_1+\\lambda_2,\\qquad q=\\lambda_1\\lambda_2.")
    lines.append("$$")
    lines.append("")
    lines.append("Khi đó")
    lines.append("")
    lines.append("$$")
    lines.append(f"z_{{{analysis.n_index + 2}}}-pz_{{{analysis.n_index + 1}}}+qz_{{{analysis.n_index}}}=0.")
    lines.append("$$")
    lines.append("")
    lines.append(f"Chọn hai thành phần độc lập là thành phần {r + 1} và {s + 1}, ta được hệ:")
    lines.append("")
    lines.append("$$")
    lines.append("\\begin{cases}")
    for idx, row in enumerate(rows):
        suffix = "\\\\" if idx < len(rows) - 1 else ""
        lines.append(
            f"{format_complex(analysis.z_n2[row], 7)}"
            f"-p\\left({format_complex(analysis.z_n1[row], 7)}\\right)"
            f"+q\\left({format_complex(analysis.z_n[row], 7)}\\right)=0"
            + suffix
        )
    lines.append("\\end{cases}")
    lines.append("$$")
    lines.append("")
    lines.append("Giải hệ trên được")
    lines.append("")
    lines.append("$$")
    lines.append("p=" + format_complex(analysis.p, 7) + ",\\qquad q=" + format_complex(analysis.q, 7) + ".")
    lines.append("$$")
    lines.append("")
    lines.append("Do đó các trị riêng trội là nghiệm của phương trình")
    lines.append("")
    lines.append("$$")
    lines.append("t^2-pt+q=0.")
    lines.append("$$")
    lines.append("")
    lines.append("Thay \\(p,q\\) vào, ta có")
    lines.append("")
    lines.append("$$")
    lines.append(
        "t^2-("
        + format_complex(analysis.p, 7)
        + ")t+("
        + format_complex(analysis.q, 7)
        + ")=0."
    )
    lines.append("$$")
    lines.append("")
    r1, r2 = analysis.roots
    lines.append("Suy ra")
    lines.append("")
    lines.append("$$")
    lines.append("t_1=" + format_complex(r1, 7) + ",\\qquad t_2=" + format_complex(r2, 7) + ".")
    lines.append("$$")
    lines.append("")
    lines.append("Véc tơ riêng tương ứng được tìm từ \\((A-tI)v=0\\):")
    lines.append("")
    for idx, root in enumerate((r1, r2), start=1):
        special_vector = eigenvector_for_matrix(matrix, root)
        lines.append("$$")
        lines.append(
            "t_"
            + str(idx)
            + "="
            + format_complex(root, 7)
            + ",\\qquad v_"
            + str(idx)
            + "\\approx "
            + vector_latex(special_vector, 7)
            + "."
        )
        lines.append("$$")
        lines.append("")
    return lines


def render_eigenvectors_for_roots(matrix: Matrix, roots: Sequence[complex]) -> List[str]:
    lines: List[str] = []
    lines.append("Với mỗi giá trị riêng, véc tơ riêng được tìm từ hệ")
    lines.append("")
    lines.append("$$")
    lines.append("(A-tI)v=0.")
    lines.append("$$")
    lines.append("")
    for idx, root in enumerate(roots, start=1):
        vector = eigenvector_for_matrix(matrix, root)
        lines.append("$$")
        lines.append(
            "t_"
            + str(idx)
            + "="
            + format_complex(root, 7)
            + ",\\qquad v_"
            + str(idx)
            + "\\approx "
            + vector_latex(vector, 7)
            + "."
        )
        lines.append("$$")
        lines.append("")
    return lines


def render_power_result(result: PowerResult, title: str = "Phương pháp lũy thừa") -> List[str]:
    lines: List[str] = []
    lines.append(f"## {title}")
    lines.append("")
    lines.append("Ta chọn véc tơ ban đầu")
    lines.append("")
    lines.append("$$")
    lines.append("x^{(0)}=" + vector_latex(result.initial_vector, 7))
    lines.append("$$")
    lines.append("")
    lines.append("Tại mỗi bước, tính")
    lines.append("")
    lines.append("$$")
    lines.append("y^{(k)}=Ax^{(k)}.")
    lines.append("$$")
    lines.append("")
    lines.append("Chọn chỉ số \\(s\\) sao cho")
    lines.append("")
    lines.append("$$")
    lines.append("|y_s^{(k)}|=\\max_i |y_i^{(k)}|,")
    lines.append("$$")
    lines.append("")
    lines.append("rồi chuẩn hóa")
    lines.append("")
    lines.append("$$")
    lines.append("x^{(k+1)}=\\frac{y^{(k)}}{y_s^{(k)}}.")
    lines.append("$$")
    lines.append("")
    lines.append("Cột sai số trong bảng được tính bằng chuẩn của vector dư:")
    lines.append("")
    lines.append("$$")
    lines.append("r^{(k)}=\\left\\|Ax^{(k+1)}-\\lambda^{(k)}x^{(k+1)}\\right\\|_2.")
    lines.append("$$")
    lines.append("")
    lines.extend(render_iteration_table(result.steps))
    lines.append("")
    if result.converged:
        lines.append("Dãy lặp hội tụ theo sai số đã chọn.")
        lines.append("")
        lines.append("Suy ra giá trị riêng trội xấp xỉ")
        lines.append("")
        lines.append("$$")
        lines.append("\\lambda_1\\approx " + format_complex(result.eigenvalue, 7) + ".")
        lines.append("$$")
        lines.append("")
        lines.append("Véc tơ riêng tương ứng có thể lấy là")
        lines.append("")
        lines.append("$$")
        lines.append("v_1\\approx " + vector_latex(result.eigenvector, 7) + ".")
        lines.append("$$")
        lines.append("")
    else:
        lines.append("Dãy lặp chưa thỏa điều kiện dừng trong số bước đã chọn.")
        lines.append("")
        lines.append("Khi đó chưa kết luận trực tiếp từ bảng lặp; cần xét trường hợp đặc biệt.")
        lines.append("")
    return lines


def render_classification_report(result: PowerResult) -> List[str]:
    report = result.classification
    lines: List[str] = []
    lines.append("## Kiểm tra trường hợp")
    lines.append("")
    lines.append("Ta kiểm tra lần lượt theo thứ tự:")
    lines.append("")
    lines.append("$$")
    lines.append("\\text{TH1} \\rightarrow \\text{TH2} \\rightarrow \\text{TH3}.")
    lines.append("$$")
    lines.append("")

    lines.append("### Kiểm tra TH1")
    lines.append("")
    lines.append("TH1 xảy ra khi dãy lặp hội tụ trực tiếp, tức là \\(x^{(k+1)}\\approx x^{(k)}\\) và \\(\\lambda^{(k)}\\) ổn định.")
    lines.append("")
    lines.append("Ta dùng đại lượng kiểm tra")
    lines.append("")
    lines.append("$$")
    lines.append(
        "E_1=\\max\\left\\{"
        "\\frac{|\\lambda^{(k)}-\\lambda^{(k-1)}|}{\\max(1,|\\lambda^{(k)}|)},"
        "\\frac{\\left\\|Ax^{(k+1)}-\\lambda^{(k)}x^{(k+1)}\\right\\|_2}{\\max(1,|\\lambda^{(k)}|)}"
        "\\right\\}."
    )
    lines.append("$$")
    lines.append("")
    if report.normal_error is not None:
        lines.append("Ở bước cuối, sai số kiểm tra xấp xỉ:")
        lines.append("")
        lines.append("$$")
        lines.append(
            "E_1="
            + format_real(report.normal_error, 7)
            + ",\\qquad \\varepsilon_1="
            + format_real(report.normal_threshold, 7)
            + "."
        )
        lines.append("$$")
        lines.append("")
    if result.case == "normal":
        lines.append("Vì \\(E_1\\le \\varepsilon_1\\), bài toán thuộc TH1.")
        lines.append("")
        return lines
    lines.append("TH1 không thỏa, nên tiếp tục kiểm tra TH2.")
    lines.append("")

    lines.append("### Kiểm tra TH2")
    lines.append("")
    lines.append("TH2 xảy ra khi dãy không hội tụ trực tiếp nhưng có chu kỳ 2, tức là \\(x^{(k+2)}\\approx x^{(k)}\\), đồng thời \\(\\lambda^2\\) ước lượng được là số thực dương.")
    lines.append("")
    lines.append("Gọi \\(D\\) là tập các thành phần trội được dùng để so sánh chu kỳ 2. Ta dùng")
    lines.append("")
    lines.append("$$")
    lines.append(
        "E_2=\\frac{\\left\\|x_D^{(k)}-x_D^{(k-2)}\\right\\|_2}"
        "{\\max\\left(1,\\left\\|x_D^{(k)}\\right\\|_2,\\left\\|x_D^{(k-2)}\\right\\|_2\\right)}."
    )
    lines.append("$$")
    lines.append("")
    if report.period2_error is not None:
        lines.append("Ở bước cuối, sai số chu kỳ 2 là:")
        lines.append("")
        lines.append("$$")
        lines.append(
            "E_2="
            + format_real(report.period2_error, 7)
            + ",\\qquad \\varepsilon_2="
            + format_real(report.period2_threshold, 7)
            + "."
        )
        lines.append("$$")
        lines.append("")
    if report.even_square_value is not None:
        lines.append("Ước lượng \\(\\lambda^2\\) từ lũy thừa chẵn:")
        lines.append("")
        lines.append("$$")
        lines.append("\\lambda^2\\approx " + format_complex(report.even_square_value, 7) + ".")
        lines.append("$$")
        lines.append("")
    if result.case == "opposite":
        lines.append("Vì \\(E_2\\le \\varepsilon_2\\) và \\(\\lambda^2>0\\), bài toán thuộc TH2.")
        lines.append("")
        return lines
    lines.append("TH2 không thỏa, nên chuyển sang TH3.")
    lines.append("")

    lines.append("### Kiểm tra TH3")
    lines.append("")
    if result.case == "quadratic":
        lines.append("Vì TH1 và TH2 đều không thỏa, ta dùng hệ bậc hai")
        lines.append("")
        lines.append("$$")
        lines.append("z_{n+2}-pz_{n+1}+qz_n=0")
        lines.append("$$")
        lines.append("")
        lines.append("để tìm cặp trị riêng trội.")
    else:
        lines.append("Chưa đủ dữ liệu để kết luận chắc chắn TH3.")
    lines.append("")
    return lines


def should_show_special_checks(solution: DominantEigenSolution) -> bool:
    option = solution.options.show_special_checks
    if option == "always":
        return True
    if option == "never":
        return False
    if solution.options.mode == "full":
        return True
    return solution.main_result.case != "normal"


def render_markdown(solution: DominantEigenSolution) -> str:
    lines: List[str] = []
    result = solution.main_result

    lines.append("# Lời giải tìm giá trị riêng trội")
    lines.append("")
    lines.append("> Các phép tính bên trong dùng giá trị gốc. Ma trận trình bày 4 chữ số thập phân; các số liệu khác 7 chữ số thập phân.")
    lines.append("")
    lines.append("## Bài toán")
    lines.append("")
    lines.append("Cho ma trận")
    lines.append("")
    lines.append("$$")
    lines.append("A=" + matrix_latex(solution.matrix, 4))
    lines.append("$$")
    lines.append("")
    lines.append("Dùng phương pháp lũy thừa để tìm giá trị riêng trội và véc tơ riêng tương ứng.")
    lines.append("")
    lines.append("Các tham số sử dụng:")
    lines.append("")
    lines.append(f"- Sai số: \\(\\varepsilon={format_real(solution.tol, 7)}\\)")
    lines.append(f"- Số bước lặp tối đa: \\(N={solution.max_iter}\\)")
    lines.append(f"- Chế độ xuất lời giải: `{solution.options.mode}`")
    lines.append("")
    lines.extend(render_power_result(result))
    lines.extend(render_classification_report(result))

    if should_show_special_checks(solution):
        if result.case == "normal":
            lines.append("## Kiểm tra phụ")
            lines.append("")
        elif result.case == "opposite":
            lines.append("## Trường hợp cặp đối dấu")
            lines.append("")
        elif result.case == "quadratic":
            lines.append("## Trường hợp cặp nghiệm bậc hai")
            lines.append("")
        else:
            lines.append("## Kiểm tra trường hợp đặc biệt")
            lines.append("")

        if result.ratio_estimates and solution.options.mode == "full":
            lines.append("Với lũy thừa thô \\(z_k=A^kx^{(0)}\\), theo slide có thể ước lượng")
            lines.append("")
            lines.append("$$")
            lines.append("\\lambda_1\\approx \\frac{(z_{k+1})_i}{(z_k)_i}.")
            lines.append("$$")
            lines.append("")
            lines.append("Các tỉ số thành phần ở bước cuối là:")
            lines.append("")
            for index, ratio in enumerate(result.ratio_estimates, start=1):
                lines.append(f"- Thành phần {index}: \\({format_complex(ratio, 7)}\\)")
            lines.append("")

        if result.case == "opposite":
            if result.even_square_analysis is not None:
                lines.extend(render_even_square_analysis(result.even_square_analysis))
            roots = dominant_roots_from_result(result)
            if roots is not None:
                lines.extend(render_eigenvectors_for_roots(result.matrix, roots))
            if solution.options.mode == "full" and result.quadratic_analysis is not None:
                lines.append("Có thể kiểm tra lại bằng hệ bậc hai:")
                lines.append("")
                lines.extend(render_quadratic_analysis(result.quadratic_analysis, result.matrix))

        if result.case == "quadratic" and result.quadratic_analysis is not None:
            lines.extend(render_quadratic_analysis(result.quadratic_analysis, result.matrix))

    if solution.deflations:
        lines.append("## Phương pháp xuống thang")
        lines.append("")
        lines.append("Sau khi tìm được \\((\\lambda_1,v_1)\\), tìm véc tơ riêng trái \\(w_1\\) từ \\(A^T w_1=\\lambda_1w_1\\), rồi đặt")
        lines.append("")
        lines.append("$$")
        lines.append("A_1=A-\\frac{\\lambda_1}{w_1^Tv_1}v_1w_1^T.")
        lines.append("$$")
        lines.append("")
        for item in solution.deflations:
            lines.append(f"### Xuống thang lần {item.index}")
            lines.append("")
            lines.append(f"Từ {item.source_case}, ta lấy trị riêng và véc tơ riêng phải:")
            lines.append("")
            lines.append("$$")
            lines.append(
                "\\lambda_"
                + str(item.index)
                + "\\approx "
                + format_complex(item.eigenvalue, 7)
                + ",\\qquad v_"
                + str(item.index)
                + "\\approx "
                + vector_latex(item.right_vector, 7)
            )
            lines.append("$$")
            lines.append("")
            lines.append("Tìm véc tơ riêng trái từ hệ \\((A^T-\\lambda I)w=0\\):")
            lines.append("")
            lines.append("$$")
            lines.append("w_" + str(item.index) + "\\approx " + vector_latex(item.left_vector, 7))
            lines.append("$$")
            lines.append("")
            lines.append("Mẫu số")
            lines.append("")
            lines.append("$$")
            lines.append("w_" + str(item.index) + "^Tv_" + str(item.index) + "=" + format_complex(item.denominator, 7) + ".")
            lines.append("$$")
            lines.append("")
            lines.append("Ma trận sau khi xuống thang:")
            lines.append("")
            lines.append("$$")
            lines.append("A_" + str(item.index) + "=" + matrix_latex(item.deflated_matrix, 4))
            lines.append("$$")
            lines.append("")
            if item.next_result is not None:
                lines.append(f"Tiếp tục áp dụng phương pháp lũy thừa cho \\(A_{item.index}\\), ta thu được giá trị riêng trội tiếp theo")
                lines.append("")
                lines.append("$$")
                lines.append(
                    "\\lambda_"
                    + str(item.index + 1)
                    + "\\approx "
                    + format_complex(item.next_result.eigenvalue, 7)
                    + "."
                )
                lines.append("$$")
                lines.append("")
                lines.append("Véc tơ riêng tương ứng trong bài toán sau xuống thang là")
                lines.append("")
                lines.append("$$")
                lines.append(
                    "\\tilde v_"
                    + str(item.index + 1)
                    + "\\approx "
                    + vector_latex(item.next_result.eigenvector, 7)
                    + "."
                )
                lines.append("$$")
                lines.append("")

    lines.append("## Kết luận")
    lines.append("")
    if result.converged:
        lines.append("Giá trị riêng trội tìm được là")
        lines.append("")
        lines.append("$$")
        lines.append("\\lambda_1\\approx " + format_complex(result.eigenvalue, 7) + ".")
        lines.append("$$")
        lines.append("")
        lines.append("Véc tơ riêng tương ứng là")
        lines.append("")
        lines.append("$$")
        lines.append("v_1\\approx " + vector_latex(result.eigenvector, 7) + ".")
        lines.append("$$")
        lines.append("")
    elif result.case == "opposite" and dominant_roots_from_result(result):
        r1, r2 = dominant_roots_from_result(result) or (0j, 0j)
        lines.append("Do dãy lũy thừa không hội tụ trực tiếp mà dao động theo chu kỳ 2, ta dùng lũy thừa chẵn để lấy cặp giá trị riêng trội đối dấu:")
        lines.append("")
        lines.append("$$")
        lines.append("t_1=" + format_complex(r1, 7) + ",\\qquad t_2=" + format_complex(r2, 7) + ".")
        lines.append("$$")
        lines.append("")
    elif result.case == "quadratic" and dominant_roots_from_result(result):
        r1, r2 = dominant_roots_from_result(result) or (0j, 0j)
        lines.append("Do dãy lũy thừa không hội tụ trực tiếp, dùng hệ bậc hai ở trên để lấy cặp giá trị riêng trội:")
        lines.append("")
        lines.append("$$")
        lines.append("t_1=" + format_complex(r1, 7) + ",\\qquad t_2=" + format_complex(r2, 7) + ".")
        lines.append("$$")
        lines.append("")
    else:
        lines.append("Dãy lặp chưa hội tụ trong số bước đã chọn, nên chưa thể kết luận chắc chắn giá trị riêng trội.")
        lines.append("")
    if solution.deflations:
        lines.append("Các giá trị riêng tiếp theo tìm được bằng xuống thang:")
        lines.append("")
        for item in solution.deflations:
            if item.next_result is None:
                continue
            lines.append(
                f"- \\(\\lambda_{item.index + 1}\\approx "
                + format_complex(item.next_result.eigenvalue, 7)
                + "\\)"
            )
        lines.append("")
    return "\n".join(lines)


def solve_text_to_markdown(text: str) -> str:
    return render_markdown(solve_dominant_eigen(text))


def save_markdown_from_text(text: str, output_path: str) -> None:
    markdown = solve_text_to_markdown(text)
    with open(output_path, "w", encoding="utf-8") as file:
        file.write(markdown)

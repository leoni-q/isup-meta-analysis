"""Render REPORT.md as an A4 PDF using Python Markdown and local Chrome."""

from __future__ import annotations

import html
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import markdown  # type: ignore[import-not-found]


ANALYSIS_DIR = Path(__file__).resolve().parent
REPORT_PATH = ANALYSIS_DIR / "REPORT.md"
PDF_PATH = ANALYSIS_DIR / "REPORT.pdf"

RR_VARIANTS = (
    ("invitation", "screening invitation vs no organised screening"),
    ("diagnostic", "diagnostic/biopsy strategy"),
    ("timing", "PROBASE timing-of-screening contrast"),
)
PROPORTION_VARIANTS = (
    ("no_screening_usual_care", "no screening / usual care"),
    ("psa_dre", "PSA + DRE"),
    ("psa", "PSA"),
    ("psa_mri", "PSA → MRI"),
    ("biomarker_mri", "biomarker → MRI"),
    ("mri_first", "MRI-first"),
)
# Biopsy yield has no biomarker → MRI unit: STHLM3-MRI is excluded from those
# outcomes and ProScreen from every analysis.
BIOPSY_YIELD_VARIANTS = tuple(
    variant
    for variant in PROPORTION_VARIANTS[2:]
    if variant[0] != "biomarker_mri"
)
FIGURE_GROUPS = {
    "fig1_rr_gg1_detection": (
        "GG1 detection per randomized/invited participant",
        RR_VARIANTS,
    ),
    "fig2_rr_gg1_share_of_pca": (
        "GG1 among all detected prostate cancers",
        RR_VARIANTS,
    ),
    "fig3_rr_gg2plus_detection": (
        "GG≥2 detection per randomized/invited participant",
        RR_VARIANTS,
    ),
    "fig4_rr_gg2plus_share_of_pca": (
        "GG≥2 among all detected prostate cancers",
        RR_VARIANTS,
    ),
    "fig5_rr_gg1_per_biopsy": (
        "GG1 among participants biopsied",
        RR_VARIANTS[1:],
    ),
    "fig6_rr_gg2plus_per_biopsy": (
        "GG≥2 among participants biopsied",
        RR_VARIANTS[1:],
    ),
    "fig7_rr_biopsies_per_participant": (
        "Biopsies performed per randomized/invited participant",
        RR_VARIANTS[1:],
    ),
    "fig8_prop_gg1_per_screened": (
        "GG1 per screened participant",
        PROPORTION_VARIANTS,
    ),
    "fig9_prop_gg1_per_randomized": (
        "GG1 per randomized/invited participant",
        PROPORTION_VARIANTS,
    ),
    "fig10_prop_gg1_share_of_pca": (
        "GG1 among all detected prostate cancers",
        PROPORTION_VARIANTS,
    ),
    "fig11_prop_gg2plus_per_screened": (
        "GG≥2 per screened participant",
        PROPORTION_VARIANTS,
    ),
    "fig12_prop_gg2plus_per_randomized": (
        "GG≥2 per randomized/invited participant",
        PROPORTION_VARIANTS,
    ),
    "fig13_prop_gg2plus_share_of_pca": (
        "GG≥2 among all detected prostate cancers",
        PROPORTION_VARIANTS,
    ),
    "fig14_prop_gg1_per_biopsy": (
        "GG1 among participants biopsied",
        BIOPSY_YIELD_VARIANTS,
    ),
    "fig15_prop_gg2plus_per_biopsy": (
        "GG≥2 among participants biopsied",
        BIOPSY_YIELD_VARIANTS,
    ),
    "fig16_prop_biopsies_per_screened": (
        "Biopsies per screened/assessed participant",
        PROPORTION_VARIANTS[2:],
    ),
}
FIGURE_MARKER_PATTERN = re.compile(
    r"<!--\s*figures:\s*([a-z0-9_]+)\s*-->"
)


def render_figure_group(prefix: str) -> str:
    title, variants = FIGURE_GROUPS[prefix]
    figure_number_match = re.match(r"fig(\d+)_", prefix)
    if figure_number_match is None:
        raise RuntimeError(f"Invalid figure prefix: {prefix}")
    figure_number = figure_number_match.group(1)

    rendered_figures = []
    for panel_index, (variant, variant_label) in enumerate(variants):
        image_path = ANALYSIS_DIR / "output" / f"{prefix}_{variant}.png"
        if not image_path.is_file():
            raise FileNotFoundError(
                f"Missing generated figure: {image_path}. "
                "Run meta_analysis.R before rendering the report."
            )

        panel_letter = chr(ord("a") + panel_index)
        caption = (
            f"Figure {figure_number}{panel_letter}. Forest plot for {title} — "
            f"{variant_label}."
        )
        rendered_figures.append(
            "\n".join(
                [
                    '<figure class="report-figure">',
                    (
                        f'<img src="output/{html.escape(image_path.name)}" '
                        f'alt="{html.escape(caption)}">'
                    ),
                    f"<figcaption>{html.escape(caption)}</figcaption>",
                    "</figure>",
                ]
            )
        )
    return "\n\n".join(rendered_figures)


def inject_figure_groups(markdown_text: str) -> str:
    markers = FIGURE_MARKER_PATTERN.findall(markdown_text)
    expected_markers = set(FIGURE_GROUPS)
    if set(markers) != expected_markers or len(markers) != len(expected_markers):
        missing = sorted(expected_markers - set(markers))
        unexpected = sorted(set(markers) - expected_markers)
        raise RuntimeError(
            "REPORT.md figure markers do not match the renderer configuration. "
            f"Missing: {missing or 'none'}; unexpected: {unexpected or 'none'}."
        )

    return FIGURE_MARKER_PATTERN.sub(
        lambda match: render_figure_group(match.group(1)),
        markdown_text,
    )


def find_chrome() -> str:
    candidates = [
        "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
        "/Applications/Chromium.app/Contents/MacOS/Chromium",
        shutil.which("google-chrome"),
        shutil.which("chromium"),
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return str(candidate)
    raise RuntimeError("Google Chrome or Chromium is required to render REPORT.pdf")


def build_html(markdown_text: str) -> str:
    body = markdown.markdown(
        inject_figure_groups(markdown_text),
        extensions=["tables", "fenced_code", "sane_lists"],
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<base href="{ANALYSIS_DIR.as_uri()}/">
<style>
  @page {{
    size: A4;
    margin: 18mm 16mm 18mm;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0 auto;
    max-width: 180mm;
    color: #1f2933;
    font: 10.2pt/1.45 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  }}
  h1, h2, h3 {{
    color: #102a43;
    line-height: 1.2;
    page-break-after: avoid;
  }}
  h1 {{ font-size: 23pt; margin: 0 0 10mm; }}
  h2 {{
    font-size: 16pt;
    margin: 9mm 0 4mm;
    border-bottom: 1px solid #bcccdc;
    padding-bottom: 2mm;
  }}
  h3 {{ font-size: 12.5pt; margin: 6mm 0 3mm; }}
  p, li {{ orphans: 3; widows: 3; }}
  a {{ color: #0b69a3; text-decoration: none; }}
  hr {{ border: 0; border-top: 1px solid #d9e2ec; margin: 8mm 0; }}
  table {{
    width: 100%;
    border-collapse: collapse;
    margin: 3mm 0 6mm;
    font-size: 8.6pt;
    page-break-inside: avoid;
  }}
  th, td {{
    border: 1px solid #bcccdc;
    padding: 1.8mm 2mm;
    text-align: left;
    vertical-align: top;
  }}
  th {{ background: #eaf1f8; color: #102a43; }}
  tr:nth-child(even) td {{ background: #f7f9fb; }}
  figure.report-figure {{
    margin: 5mm 0 8mm;
    break-inside: avoid;
    page-break-inside: avoid;
  }}
  figure.report-figure img {{
    display: block;
    width: 100%;
    max-height: 235mm;
    object-fit: contain;
  }}
  figure.report-figure figcaption {{
    margin-top: 2mm;
    color: #334e68;
    font-size: 8.8pt;
    line-height: 1.35;
    text-align: center;
  }}
  code {{
    background: #f0f4f8;
    border-radius: 2px;
    padding: 0.2mm 0.8mm;
    font-size: 0.92em;
  }}
  ul, ol {{ padding-left: 6mm; }}
</style>
</head>
<body>{body}</body>
</html>
"""


def main() -> None:
    html = build_html(REPORT_PATH.read_text(encoding="utf-8"))
    chrome = find_chrome()
    temporary_path: Path | None = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".html",
            dir=ANALYSIS_DIR,
            encoding="utf-8",
            delete=False,
        ) as temporary_file:
            temporary_file.write(html)
            temporary_path = Path(temporary_file.name)

        subprocess.run(
            [
                chrome,
                "--headless",
                "--disable-gpu",
                "--no-sandbox",
                "--allow-file-access-from-files",
                "--no-pdf-header-footer",
                f"--print-to-pdf={PDF_PATH}",
                temporary_path.as_uri(),
            ],
            check=True,
        )
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    print(f"Wrote {PDF_PATH}")


if __name__ == "__main__":
    main()

import pytask

from src.config import BLD
from src.config import SRC

# Run application.R (creates figures)

dependencies_task_run_r_script = [
    "application.R",
    SRC / "application" / "data" / "lee_2008.dta",
]

application_figures = [
    SRC / "application" / "figures" / "application_figure_01.pdf",
    SRC / "application" / "figures" / "application_figure_02.pdf",
    SRC / "application" / "figures" / "application_figure_03a.pdf",
    SRC / "application" / "figures" / "application_figure_03b.pdf",
    SRC / "application" / "figures" / "application_figure_04a.pdf",
    SRC / "application" / "figures" / "application_figure_04b.pdf",
    SRC / "application" / "figures" / "application_figure_05.pdf",
]


@pytask.mark.r
@pytask.mark.depends_on(dependencies_task_run_r_script)
@pytask.mark.produces(application_figures)
def task_run_r_script():
    pass


# Compile application.tex (creates PDF)

dependencies_task_compile_tex_file = [
    "application.tex",
    "application.bib",
    SRC / "application" / "figures" / "application_figure_01.pdf",
    SRC / "application" / "figures" / "application_figure_02.pdf",
    SRC / "application" / "figures" / "application_figure_03a.pdf",
    SRC / "application" / "figures" / "application_figure_03b.pdf",
    SRC / "application" / "figures" / "application_figure_04a.pdf",
    SRC / "application" / "figures" / "application_figure_04b.pdf",
    SRC / "application" / "figures" / "application_figure_05.pdf",
]


@pytask.mark.latex
@pytask.mark.depends_on(dependencies_task_compile_tex_file)
@pytask.mark.produces(BLD / "application.pdf")
def task_compile_tex_file():
    pass

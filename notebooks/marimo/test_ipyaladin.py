import marimo

__generated_with = "0.10.0"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    from ipyaladin import Aladin
    return Aladin, mo


@app.cell
def _(Aladin, mo):
    aladin = Aladin(target="M1", fov=2)
    widget = mo.ui.anywidget(aladin)
    return aladin, widget


@app.cell
def _(widget):
    widget
    return


if __name__ == "__main__":
    app.run()

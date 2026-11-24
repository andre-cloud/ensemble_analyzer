from pathlib import Path


def append_to_input(calculator, text: str, writer: str, input_name) -> None:
    """Append ``text`` to the input file once ``calculator.<writer>`` has written it.

    ``input_name`` is a callable returning the input file name inside
    ``calculator.directory``. Earlier patches are chained, not replaced.
    """
    original = getattr(calculator, writer)

    def patched(*args, **kwargs):
        original(*args, **kwargs)
        with open(Path(calculator.directory) / input_name(), "a") as f:
            f.write("\n" + text + "\n")

    setattr(calculator, writer, patched)

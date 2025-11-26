@ -1,40 +0,0 @@
import logging
import os

LOG_FORMAT = "%(message)s"

DEBUG = bool(os.getenv("DEBUG"))


def ordinal(n):
    return "%d-%s" % (n, "tsnrhtdd"[(n // 10 % 10 != 1) * (n % 10 < 4) * n % 10 :: 4])


def create_log(output):
    """
    Creating an logger instance.

    :param output: output filename
    :type output: str

    :return: logger instance
    :rtype: logging
    """

    print(output)
    logging.basicConfig(
        filename=output,
        level=logging.DEBUG if DEBUG else logging.INFO,
        format=LOG_FORMAT,
        filemode="w",
    )

    log = logging.getLogger()

    log.warning(f"DEBUG mode: {DEBUG}")

    return logging.getLogger()


if __name__ == "__main__":  # pragma: no cover:
    log = create_log("output_test.out")
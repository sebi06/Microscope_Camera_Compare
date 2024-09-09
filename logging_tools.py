from loguru import logger as loguru_logger
import sys


def set_logging():

    loguru_logger.remove()
    loguru_logger.add(
        sys.stdout,
        colorize=True,
        format="<green>{time}s</green> - <level>{level}</level> - <level>{message}</level>",
    )

    return loguru_logger

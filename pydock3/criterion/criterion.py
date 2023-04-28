from typing import NoReturn


class Criterion(object):
    def __init__(self):
        pass

    @property
    def name(self) -> NoReturn:
        raise NotImplementedError

    def calculate(self, *args, **kwargs) -> NoReturn:
        raise NotImplementedError

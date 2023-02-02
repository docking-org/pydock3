from typing import NoReturn, Iterable, Union

from pydock3.dockopt.roc import ROC


class Criterion(object):
    def __init__(self):
        pass

    @property
    def name(self) -> NoReturn:
        raise NotImplementedError

    def calculate(self, *args, **kwargs) -> NoReturn:
        raise NotImplementedError


class EnrichmentScore(Criterion):
    def __init__(self):
        super().__init__()

    @property
    def name(self) -> str:
        return "enrichment_score"

    def calculate(
        self,
        booleans: Iterable[bool],
        image_save_path: Union[None, str] = None
    ) -> float:
        roc = ROC(booleans)

        # save plot
        if image_save_path is not None:
            roc.plot(save_path=image_save_path)

        return roc.enrichment_score

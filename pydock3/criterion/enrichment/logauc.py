from typing import Iterable, Union

from pydock3.criterion.criterion import Criterion
from pydock3.criterion.enrichment.roc import ROC


class NormalizedLogAUC(Criterion):
    def __init__(self):
        super().__init__()

    @property
    def name(self) -> str:
        return "normalized_log_auc"

    def calculate(
        self,
        booleans: Iterable[bool],
        image_save_path: Union[None, str] = None
    ) -> float:
        roc = ROC(booleans)

        # save plot
        if image_save_path is not None:
            roc.plot(save_path=image_save_path)

        return roc.normalized_log_auc

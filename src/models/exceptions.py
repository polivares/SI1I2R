__author__ = "Patricio Andrés Olivares Roncagliolo"
__email__ = "patricio.olivaresr@usm.cl"

class NotEvaluatedError(Exception):
    """Exception raised when trying to use some functions on models that are not yet
    evaluated
    """
    
    def __init__(self, *args: object) -> None:
        super().__init__(*args)
    